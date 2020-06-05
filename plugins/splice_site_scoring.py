from __future__ import division, print_function
import os, struct, math, sys
from pyfaidx import Fasta

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
    def __init__(self, hooks, splice_scoring, donor_path, acceptor_path, fasta_path, valid_donor = 4, valid_acceptor = 4, **kwargs):
        if not splice_scoring:
            return
        
        self.donorPredictor = DonorPredictor(donor_path)
        self.acceptorPredictor = AcceptorPredictor(acceptor_path)
        self.chromosomes = Fasta(fasta_path)

        self.valid_donor = valid_donor
        self.valid_acceptor = valid_acceptor

        hooks["contig_merged"].tap(self.add_splice_site_score)

    def add_splice_site_score(self, transcripts):
        for transcript in transcripts:
            chrom = transcript.chromosome
            strand = transcript.strand

            cum_donor_score = 0
            cum_acceptor_score = 0
            max_donor_score = None
            max_acceptor_score = None
            min_donor_score = None
            min_acceptor_score = None

            for junction in transcript.junctions:
                intron_start = junction[0] + 1
                intron_end = junction[1] - 1
                
                left_start = intron_start - 1
                right_start = intron_end - 2
                left_end = left_start + 2
                right_end = right_start + 2

                if transcript.strand == "+":
                    donor_start = left_start - 3
                    donor_end = left_end + 4
                    acceptor_start = right_start - 18 # WHy the second - 18??? 
                    acceptor_end = right_end + 3
                else:
                    donor_start = right_start - 4
                    donor_end = right_end + 3
                    acceptor_start = left_start - 3
                    acceptor_end = left_end + 18

                donor = self.chromosomes[chrom][donor_start:donor_end].seq.upper()
                acceptor = self.chromosomes[chrom][acceptor_start:acceptor_end].seq.upper()
                
                donorDiNt= donor[3:5]
                acceptorDiNt=acceptor[18:20]
                
                if (donorDiNt == 'GT' and acceptorDiNt == 'AG' and len(donor) == 9 and len(acceptor) == 23 and 'N' not in donor and 'N' not in acceptor):
                    donorFlank1=donor[:3]
                    donorFlank2=donor[5:]
                    acceptorFlank1=acceptor[:18]
                    acceptorFlank2=acceptor[20:]
                    scoreDonor=self.donorPredictor(donorFlank1, donorFlank2)
                    scoreAcceptor=self.acceptorPredictor(acceptorFlank1, acceptorFlank2)

                    if scoreDonor < self.valid_donor:
                        transcript.remove()
                        break
                        max_donor_score = scoreDonor
                    if scoreAcceptor < self.valid_acceptor:
                        transcript.remove()
                        break
                    
                    if not max_donor_score or scoreDonor > max_donor_score:
                        max_donor_score = scoreDonor
                    if not min_donor_score or scoreDonor < min_donor_score:
                        min_donor_score = scoreDonor
                    
                    if not max_acceptor_score or scoreAcceptor > max_acceptor_score:
                        max_acceptor_score = scoreAcceptor
                    if not min_acceptor_score or scoreAcceptor < min_acceptor_score:
                        min_acceptor_score = scoreAcceptor

                    cum_donor_score += scoreDonor
                    cum_acceptor_score += scoreAcceptor

                else:
                    # Not canonical, remove anyway
                    transcript.remove()
                    break

            transcript.meta["max_acceptor_score"] = max_acceptor_score
            transcript.meta["max_donor_score"] = max_donor_score
            transcript.meta["min_acceptor_score"] = min_acceptor_score
            transcript.meta["min_donor_score"] = min_donor_score
            transcript.meta["cum_acceptor_score"] = cum_acceptor_score
            transcript.meta["cum_donor_score"] = cum_donor_score
