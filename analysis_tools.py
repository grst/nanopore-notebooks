#!/usr/bin/env python
# -*- coding: utf-8 -*-

from skbio.alignment import StripedSmithWaterman, local_pairwise_align_ssw
from itertools import repeat

"""
Contains helper classes for the analysis in the ipython notebooks.
"""

class KmerAligner():
    """
    Align a read to a know reference sequence.
    Performs sequence statistics (s.a. sequence identity) and
    returns the set of identical kmers.
    """

    def __init__(self, ref_seq, verbose=False, nmer=5):
        """
        Args:
            ref_seq: reference nucleotide sequence
            nmers: length of nmers to analyse (max. 5)
        """
        assert nmer <= 5
        self.nmer = nmer
        self.ref_seq = ref_seq
        self.verbose = verbose
        self._mk_aligner();
        self.counters = dict(zip(['nt_identity', "target_skipped", "read_length", "aligned_read_length",  "true_events"], repeat(0)));

        print ("Query Length {0}".format(len(ref_seq)));

    def _vprint(self, msg):
        if (self.verbose):
            print(msg)

    def _mk_aligner(self):
        self.align_ref = StripedSmithWaterman(self.ref_seq,
                 match_score = 6,
                 mismatch_score = 0,
                 gap_open_penalty=6,
                 gap_extend_penalty=6) #query

    def find_seq_start(self, events, start):
        """
        find the kmer that corresponds to the first kmer of the aligned sequences
        """
        moved = 0
        for i, ev in enumerate(events):
            if moved + ev["move"] >= start:
                offset = moved - start
                return i, offset
            else:
                moved += ev["move"]

    def gapmove(self, i_seq, ts, move):
        """
        moves by *move* steps, respecting gaps.
        """
        moved = 0
        if i_seq < 0:
            return move
        for i, c in enumerate(ts[i_seq:]):
            if moved == move:
                return i
            if c != "-": moved +=1

    def sequence_identity(self, alignment):
        """
        get the sequence identity of the alignment.
        Returns:
            (correctly mapped, target-length)
        """
        nt_identity = 0
        aqs = str(alignment.aligned_query_sequence)
        ats = str(alignment.aligned_target_sequence)
        for i in range(len(alignment.aligned_query_sequence)):
            if aqs[i] == ats[i]:
                nt_identity += 1
        return nt_identity

    def alignment_statistics(self, alignment):
        """
        calculate stats such as seq identity and add
        them to the global counter.
        """
        counters = dict(zip(self.counters.keys(), repeat(0)));
        """ check if part of the read was omitted due to poor alignment """
        t_beg = alignment.target_begin
        t_end = alignment.target_end_optimal
        len_read = len(alignment.target_sequence)
        t_diff = len_read - (t_end - t_beg)

        counters["read_length"] = len_read
        counters["aligned_read_length"] = len_read - t_diff
        counters["target_skipped"] = t_diff
        counters["nt_identity"] = self.sequence_identity(alignment)

        for key in counters.keys():
            self.counters[key] += counters[key]

        self.print_event_statistics(alignment, counters)

        return counters

    def print_event_statistics(self, alignment, counters):
        msg = """
        Start/End Query {0}/{1} Target {2}/{3} Skipped {4}
        Length {5}/{6} Identity {7:.5f}
        """.format(
            alignment.query_begin,
            alignment.query_end,
            alignment.target_begin,
            alignment.target_end_optimal,
            counters["target_skipped"],
            counters["aligned_read_length"],
            counters["read_length"],
            counters["nt_identity"]/counters["aligned_read_length"]
        )
        self._vprint(msg)

    def print_statistics(self):
        msg = """
        length of all reads {0}
        length of aligned parts of the reads {1}
        ratio nucleotides not aligned {2:.5f}
        sequence identity of aligned reads {3:.5f}
        fully correct kmers {4}
        ratio fully correct kmers {5:.5f}
        """.format(
            self.counters["read_length"],
            self.counters["aligned_read_length"],
            self.counters["target_skipped"]/self.counters["read_length"],
            self.counters["nt_identity"]/self.counters["aligned_read_length"],
            self.counters["true_events"],
            self.counters["true_events"]/self.counters["read_length"]
        )
        self._vprint(msg)


    def align_file(self, file_obj):
        """
        Align a (processed) fast5 file to the reference.

        Args:
            file_obj: dict ["called_seq", "events" => [...], "channel", "file_id"]

        Returns:
            list of correctly mapped (identity 100%) events.

        """

        self._vprint("processing ch {0} file {1}".format(file_obj["channel"], file_obj["file_id"]))
        alignment = self.align_ref(file_obj["called_seq"]) #target
        self.alignment_statistics(alignment)
        true_events = []

        events = file_obj["events"]
        pos_kmer, offset = self.find_seq_start(events, alignment.target_begin)
        i_seq = offset #in case the last kmer was shifted by 2, i = sequence index of the target sequence with gaps
        qs = alignment.aligned_query_sequence #the reference
        ts = alignment.aligned_target_sequence #the called seq
        for ev in events[pos_kmer:]:
            i_seq += self.gapmove(i_seq, ts, ev["move"])
            ev["channel"] = file_obj["channel"]

            q_kmer = qs[i_seq:i_seq+self.nmer]
            t_kmer = ts[i_seq:i_seq+self.nmer]
            t_kmer_gapless = ts[i_seq:].replace("-","")[:self.nmer]

            if len(t_kmer_gapless) < self.nmer:
                """end of sequence"""
                break

            assert ev["kmer"][:self.nmer] == t_kmer_gapless, "kmer from event file does not match aligned target sequence. "
            # print(i_seq)
            # print(len(ts))
            # print(" " * i_seq + t_kmer)
            # print(" " * i_seq + q_kmer)
            # print(" " * i_seq + t_kmer_gapless)

            if q_kmer == t_kmer:
                true_events.append(ev)


        self.counters["true_events"] += len(true_events);

        return true_events;

