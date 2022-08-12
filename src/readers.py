import sys
import gzip

""" All Classes assume that queries are made in sorted order (first
alphabetically by chromosome and then numerically by pos) and that input files
are sorted """

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class BedReader():
    def _open_file(self):
        if self.fname.endswith('.bed.gz'):
            self.f = gzip.open(self.fname, 'rt')
        elif self.fname.endswith('.bed'):
            self.f = open(self.fname)
        else:
            raise Exception('bed file should end with ".bed" or ".bed.gz"')

    def __init__(self,fname):
        self.fname = fname
        self._open_file()

    def reset(self):
        self.f.close()
        self._open_file()

    def __iter__(self):
        for line in self.f:
            L = line.split()
            chrom, start, end = L[:3]
            for pos in range(int(start), int(end)):
                yield chrom, pos, 1

class AllAutoReader():
    def __init__(self, tb):
        self.tb = tb
        if any(x.startswith('chr') for x in tb.chroms()):
            prefix = 'chr'
        else:
            prefix = ''
        self.autosomes = [prefix + str(x) for x in range(1,23)]
    def __iter__(self):
        for chrom in self.autosomes:
            start = 0
            end = self.tb.chroms()[chrom]
            for pos in range(int(start), int(end)):
                yield chrom, pos, 1


# class SortedBed:
#     def __init__(self, bedname, toint=False):
#         self.name = bedname
#         self.f = open(bedname)
#         self.bed_chrom = 'chr0'
#         self.bed_start = -1
#         self.bed_end = -1
#         self.toint = toint

#     def reset(self):
#         self.f.close()
#         self.f = open(self.name)
#         self.bed_chrom = 'chr0'
#         self.bed_start = -1
#         self.bed_end = -1
        
#     def query(self, chrom, pos):
#         # TODO: check if the input is sorted
#         while self.bed_chrom < chrom or (self.bed_chrom == chrom and self.bed_end <= pos):
#             bedL = self.f.readline().split()
#             if len(bedL) == 0:
#                 break
#             if (bedL[0] == self.bed_chrom and
#                 ((int(bedL[1])+1) < self.bed_start or
#                  (int(bedL[2])+1) < self.bed_end)):
#                 eprint("bed positions not sorted:", self.name)
#                 eprint(' '.join(bedL))
#                 eprint( self.bed_chrom, self.bed_start, self.bed_end)
#                 sys.exit()
#             if bedL[0] < self.bed_chrom:
#                 eprint("bed chromosomes not sorted:", self.name)
#                 eprint(' '.join(bedL))
#                 eprint(self.bed_chrom, self.bed_start, self.bed_end)
#                 sys.exit()
#             self.bed_chrom = bedL[0]
#             self.bed_start = int(bedL[1])+1
#             self.bed_end = int(bedL[2])+1
#         if (chrom == self.bed_chrom) and ( self.bed_start <= pos) and (pos < self.bed_end):
#             if self.toint:
#                 return 1
#             else:
#                 return True
#         else:
#             if self.toint:
#                 return 0
#             else:
#                 return False



# class SortedBedGraph:
#     def __init__(self, bedname, discretizer):
#         self.name = bedname
#         self.f = open(bedname)
#         self.discretizer = discretizer
#         self.bed_chrom = 'chr0'
#         self.bed_start = -1
#         self.bed_end = -1
#         self.bed_value = "NA"
#         self.prev_chrom = '0'
#         self.prev_pos = -1

#     def query(self, chrom, pos):
#         if chrom < self.prev_chrom or (chrom==self.prev_chrom and pos < self.prev_pos):
#             eprint("BedGraph-queries should be in sorted order!")
#             eprint("chrom",chrom, "prev_chrom", self.prev_chrom, "pos", pos, "prev_pos", self.prev_pos)
#             sys.exit()
#         self.prev_pos = pos
#         self.prev_chrom = chrom
#         while self.bed_chrom < chrom or (self.bed_chrom == chrom and self.bed_end <= pos):
#             bedL = self.f.readline().split()
#             if len(bedL) == 0:
#                 break
#             if (bedL[0] == self.bed_chrom and
#                 ((int(bedL[1])+1) < self.bed_start or
#                  (int(bedL[2])+1) < self.bed_end)):
#                 eprint("bed positions not sorted:", self.name)
#                 eprint(' '.join(bedL))
#                 eprint( self.bed_chrom, self.bed_start, self.bed_end)
#                 sys.exit()
#             if bedL[0] < self.bed_chrom:
#                 eprint("bed chromosomes not sorted:", self.name)
#                 eprint(' '.join(bedL))
#                 eprint(self.bed_chrom, self.bed_start, self.bed_end)
#                 sys.exit()
#             self.bed_chrom = bedL[0]
#             self.bed_start = int(bedL[1])+1
#             self.bed_end = int(bedL[2])+1
#             self.bed_value = bedL[3]
#         if (chrom == self.bed_chrom) and ( self.bed_start <= pos) and (pos < self.bed_end):
#             return self.discretizer[self.bed_value]
#         else:
#             return 'NA'
