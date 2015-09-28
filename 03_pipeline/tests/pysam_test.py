__author__ = 'sturm'
import numpy as np
import pysam
import timeit
from tempfile import NamedTemporaryFile
SIGNIFICANCE = 0.05
from multiprocessing import Pool
import pysam
import pandas

def readstats(args, ref_seq=" " * 30000, total_reads=10):
        """statistics for a single read.

        Args:
            read: pysam.AlignedRead
        """
        filename, read_id = args
        samfile = pysam.AlignmentFile(filename)
        generator = samfile.fetch()
        read = None
        for i in range(read_id+1):
            read = generator.next()


        s = {
            "read_len": read.get_tag("ZQ"),
            "mapping_quality": read.mapping_quality,
            "aln_len": read.alen,
            "aln_score": read.get_tag("AS"),
            "aln_evalue": read.get_tag("ZE"),
            "aln_editdistance": read.get_tag("NM"),
            "mapped_nts": sum([l for op, l in read.cigartuples if op == 0]),
            "ins": sum(l for op, l in read.cigartuples if op == 1),
            "del": sum(l for op, l in read.cigartuples if op == 2),
            "subst": sum(1 for query, ref in read.get_aligned_pairs(matches_only=True)
                         if read.seq[query] != ref_seq[ref]),
            "is_significant": e2p(read.get_tag("ZE")) < SIGNIFICANCE / total_reads
        }
        assert s["read_len"] >= s["mapped_nts"]
        assert s["aln_len"] >= s["mapped_nts"]
        assert len(read.query_sequence) == s["read_len"]
        return s

class samstats:
    """generate statistics for a samfile"""
    def __init__(self, samfile, ref_seq):
        """
        Args:
            samfile: pysam.AlignmentFile
            ref_seq: (str) reference sequence
        """
        self.samfile=samfile
        self.ref_seq = ref_seq
        total_reads = [500]
        self.file_stats = {
            "total_reads": int(total_reads[0]),
            "mapped_reads": samfile.mapped,
            "total_nts": sum(r.get_tag("ZQ") for r in samfile.fetch())
        }
        start_time = timeit.default_timer()
        p = Pool(1)
        filename = "../david_calling.called.sorted.bam"
        self.read_stats = p.map(readstats, [(filename, i) for i, r in enumerate(samfile.fetch()) if not r.is_unmapped])
        print(timeit.default_timer() - start_time)

    def sumstat(self, stat):
        """sum of the stat <stat> for all reads."""
        return sum(r[stat] for r in self.read_stats)

    def print_summary(self):
        lines = [
            ["mapped_reads/total_reads", self.file_stats["mapped_reads"], self.file_stats["total_reads"]],
            ["significant_reads/total_reads", self.sumstat("is_significant"), self.file_stats["total_reads"]],
            ["mapped_nts/total_nts", self.sumstat("mapped_nts"), self.file_stats["total_nts"]],
            ["editdistance/alignment_length", self.sumstat("aln_editdistance"), self.sumstat("aln_len")],
            ["alignment_score/alignment_length", self.sumstat("aln_score"), self.sumstat("aln_len")],
            ["SNPs/mapped_nts", self.sumstat("subst"), self.sumstat("mapped_nts")],
            ["ins/mapped_nts", self.sumstat("ins"), self.sumstat("mapped_nts")],
            ["del/mapped_nts", self.sumstat("del"), self.sumstat("mapped_nts")],
        ]
        def process_line(line):
            try:
                return [str(x) for x in (line + ["{0:%}".format(float(line[1])/line[2])])]
            except ZeroDivisionError:
                return ["*"] * 3
        return [process_line(line) for line in lines]


def e2p(e):
    """convert Evalue to Pvalue (of Alignment)"""
    return 1-np.exp(-e)

ref_file = "../../../../david_eccles_bc_ideas/mouse_ref.fa"
with open(ref_file) as f:
    lines = [l for l in f.readlines()]
    ref = lines[1]

print(ref[:100])

samfile = pysam.AlignmentFile("../david_calling.called.sorted.bam")
sst = samstats(samfile, ref)
print(pandas.DataFrame(sst.print_summary()))