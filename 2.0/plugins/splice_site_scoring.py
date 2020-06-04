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
    def __init__(self, merger, donor_path, acceptor_path, fa_path):
        self.donorPredictor = DonorPredictor(donor_path)
        self.acceptorPredictor = AcceptorPredictor(acceptor_path)
        self.chromosomes = Fasta(fa_path)

        self.minDonorScore=-19
        self.minAcceptorScore=-24

        merger.hooks["contig_merged"].tap(self.add_splice_site_score)

    def add_splice_site_score(self, transcripts):
        for transcript in transcripts:
            chrom = transcript.chromosome
            strand = transcript.strand

            max_scor_donor = self.minDonorScore
            max_scor_acceptor = self.minAcceptorScore

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

                    if scoreDonor > max_scor_donor:
                        max_scor_donor = scoreDonor
                    if scoreAcceptor > max_scor_acceptor:
                        max_scor_acceptor = scoreAcceptor

            if max_scor_donor < 4 or max_scor_acceptor < 4:
                transcript.remove()

            transcript.meta["acceptor_score"] = max_scor_acceptor
            transcript.meta["donor_score"] = max_scor_donor