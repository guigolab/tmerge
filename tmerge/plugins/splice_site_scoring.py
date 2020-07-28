from __future__ import division, print_function
import os, struct, math, sys
from pyfaidx import Fasta
import logging


"""
Format of coefficient files:

The donor coefficient file consists of 4^7 numbers stored in little-endian double format.
The index into this file for a particular splice candidate is computed from the string
of 3 bases 5' of the GT concatenated to the 4 bases 3' of the GT. The index is a base-4
number created from the indices of the bases of this base string in ACGT with the first
being most significant. For example, the index for a potential splice site context of:
ACG-GT-TACG is 0*4^6 + 1*4^5 + 2*4^4 + 3*4^3 + 0*4^2 + 1*4^1 + 2*4^0
Once the coefficient is extracted, the score is
	log(16.302010666666664 * coeff, 2)

The acceptor coefficient file consists of 9 concatenated sequences containing
4^7, 4^7, 4^7, 4^7, 4^7, 4^3, 4^4, 4^3, and 4^4 numbers. The relevant context
is 18 bases 5' of the AG concatenated to 3 bases 3'. The 9 coefficients are extracted
from these arrays using indices computed from the substrings starting at the following
offsets with the following lengths:
0 7, 7 7, 14 7 (including the 3 bases 3' of the AG), 4 7, 11 7, 4 3, 7 4, 11 3, 14 4
Once these 9 coefficients are found, the score is:
    log(16.3482025 * (prod 1st 5 coeffs) / (prod other 4 coeffs), 2)
"""

class DonorPredictor() :
    def __init__(self, donor_path) :
        self.file = open(os.path.abspath(os.path.expanduser(donor_path)), "rb")
        self.file.seek(0, os.SEEK_END)
        assert self.file.tell() / RecordSize == 16384, \
            '%s is not a valid donor coefficients file.' % donor_path

    def __call__(self, prev3bases, next4bases) :
        # Score potential GT splice donor site given 3 prev bases and 4 next ones.
        assert len(prev3bases) == 3, len(prev3bases)
        assert len(next4bases) == 4, len(next4bases)
        index = _bases_to_number(prev3bases + next4bases)
        self.file.seek(index * RecordSize)
        coeff = struct.unpack(RecordFormat, self.file.read(RecordSize))[0]
        return math.log(16.302010666666664 * coeff, 2)

class AcceptorPredictor() :
    def __init__(self, acceptor_path) :
        self.file = open(os.path.abspath(os.path.expanduser(acceptor_path)), "rb")
        self.file.seek(0, os.SEEK_END)
        assert self.file.tell() / RecordSize == 82560, \
            '%s is not a valid acceptor coefficients file.' % acceptor_path

    def __call__(self, prev18bases, next3bases) :
        # Score potential AG splice acceptor site given 18 prev bases and 3 next ones.
        assert len(prev18bases) == 18, len(prev18bases)
        assert len(next3bases) == 3, len(next3bases)
        bases = prev18bases + next3bases
        coeffsCombination = 1
        for ii, (start, end) in enumerate(AcceptorStartEnds) :
            index = _bases_to_number(bases[start : end + 1])
            self.file.seek((AcceptorArrayLengthSums[ii] + index) * RecordSize)
            coeff = struct.unpack(RecordFormat, self.file.read(RecordSize))[0]
            if ii < 5 :
                coeffsCombination *= coeff
            else :
                coeffsCombination /= coeff
        return math.log(16.3482025 * coeffsCombination, 2)

def _bases_to_number(bases) :
    "Convert a string of DNA bases to a base-4 number."
    BaseMap = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}
    result = 0
    for b in bases :
        result = 4 * result + BaseMap[b]
    return result

RecordFormat = '<d' # little-endian double
RecordSize = struct.calcsize(RecordFormat) # Typically 8

# Relevant intervals for acceptor site prediction
AcceptorStartEnds = [(0, 6), (7, 13), (14, 20), (4, 10), (11, 17),
                     (4, 6), (7, 10), (11, 13), (14, 17)]
# Lengths of the coefficient arrays for acceptor sites.
AcceptorArrayLengths = [4 ** (end - start + 1) for start, end in AcceptorStartEnds]
AcceptorArrayLengthSums = [sum(AcceptorArrayLengths[:ii])
                           for ii in range(len(AcceptorArrayLengths))]


class SpliceSiteScoring():
    """
    Given the context of a GT or AG in a string of DNA bases and precomputed coefficients (see above comment for format of coefficient files)
    files, add a score to the transcript model that indicates how likely it is to be a splice donor or acceptor,
    respectively, using the maximum entropy method [1]. Filters out reads that have a donor or accetor site below the valid_donor/acceptor threshold.

    [1] Yeo, G., & Burge, C. B. (2004). Maximum entropy modeling of short sequence motifs with
    applications to RNA splicing signals. Journal of Computational Biology : a Journal of
    Computational Molecular Cell Biology, 11(2-3), 377-394. doi:10.1089/1066527041410418
    """
    def __init__(self, hooks, splice_scoring, donor_path, acceptor_path, fasta_path, valid_donor, valid_acceptor, **kwargs):
        if not splice_scoring:
            return
        
        self.donorPredictor = DonorPredictor(donor_path)
        self.acceptorPredictor = AcceptorPredictor(acceptor_path)
        self.chromosomes = Fasta(fasta_path)

        self.valid_donor = valid_donor
        self.valid_acceptor = valid_acceptor

        self.min_donor_score = -19
        self.min_acceptor_score = -24

        if not self.valid_donor:
            self.valid_donor = self.min_donor_score
        if not self.valid_acceptor:
            self.valid_acceptor = self.min_acceptor_score

        self.skipped_chroms = set()

        hooks["contig_merged"].tap(self.add_splice_site_score)
        hooks["complete"].tap(self.print_skipped_chromosomes)

    def add_splice_site_score(self, transcripts):
        # Transcripts in "transcripts" are from same contig, so same chromosome, so can stop here if the chromosome is not in the genome
        if len(transcripts) == 0:
            return
        
        if transcripts[0].chromosome not in self.chromosomes:
                self.skipped_chroms.add(transcripts[0].chromosome)
                return

        for transcript in transcripts:
            if transcript.monoexonic:
                continue

            chrom = transcript.chromosome
            strand = transcript.strand

            min_donor_score = None
            min_acceptor_score = None
            all_canonical = True
            donor_scores = []
            acceptor_scores = []

            transcript.meta["DEBUG"] = []

            for junction in transcript.sorted_junctions:
                intron_start = junction[0] + 1
                intron_end = junction[1] - 1
                
                left_start = intron_start - 1
                right_start = intron_end - 2
                left_end = left_start + 2
                right_end = right_start + 2

                if strand == '+':
                    donor = self.chromosomes[chrom][left_start-3:left_end+4].seq.upper()
                    acceptor = self.chromosomes[chrom][right_start-18:right_end+3].seq.upper()
                elif strand == '-':
                    donor = self.chromosomes[chrom][right_start-4:right_end+3].reverse.complement.seq.upper()
                    acceptor = self.chromosomes[chrom][left_start-3:left_end+18].reverse.complement.seq.upper()

                donorDiNt= donor[3:5]
                acceptorDiNt=acceptor[18:20]

                junction_canonical = donorDiNt == 'GT' and acceptorDiNt == 'AG'
                correct_lengths = len(donor) == 9 and len(acceptor) == 23
                no_unmatched_bases = 'N' not in donor and 'N' not in acceptor
                
                if not junction_canonical:
                    all_canonical = False

                if junction_canonical and correct_lengths and no_unmatched_bases:
                    donorFlank1=donor[:3]
                    donorFlank2=donor[5:]
                    acceptorFlank1=acceptor[:18]
                    acceptorFlank2=acceptor[20:]
                    scoreDonor=self.donorPredictor(donorFlank1, donorFlank2)
                    scoreAcceptor=self.acceptorPredictor(acceptorFlank1, acceptorFlank2)
                else:
                    scoreDonor=self.min_donor_score
                    scoreAcceptor=self.min_acceptor_score

                if not correct_lengths:
                    transcript.meta["DEBUG"].append("not correct lengths")
                if not no_unmatched_bases:
                    transcript.meta["DEBUG"].append("unmatched_bases")

                donor_scores.append(scoreDonor)
                acceptor_scores.append(scoreAcceptor)
                
                if not min_donor_score or scoreDonor < min_donor_score:
                    min_donor_score = scoreDonor
                if not min_acceptor_score or scoreAcceptor < min_acceptor_score:
                    min_acceptor_score = scoreAcceptor

                # If any one of the junctions does not meet score requirements, remove it
                if scoreDonor < self.valid_donor or scoreAcceptor < self.valid_acceptor:
                    transcript.remove()
                    break

            transcript.meta["min_acceptor_score"] = min_acceptor_score
            transcript.meta["min_donor_score"] = min_donor_score
            transcript.meta["all_canonical"] = all_canonical
            transcript.meta["acceptor_scores"] = acceptor_scores
            transcript.meta["donor_scores"] = donor_scores
    
    def print_skipped_chromosomes(self):
        if len(self.skipped_chroms) > 0:
            logging.warn("Splice site scoring: some input chromosomes were not in the given genomes. Skipped splice site scoring on the following chromsomes: ")
            logging.warn(", ".join(self.skipped_chroms))
