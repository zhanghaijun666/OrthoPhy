#!/usr/bin/python3

import sys
import csv
from collections import defaultdict
import fileinput
import subprocess
import re
from pathlib import Path
import pickle
import multiprocessing
# import networkx as nx
import dask.dataframe as dd
import pandas as pd


def get_hit_bitscore(blast_result):
    # 0: rbh 1: bh 2: others
    last_query     = ""
    all_hit_tuples = {}
    with open(str(blast_result), "r") as f:
        res = blast_result.stem
        query_file, subject_file = res.split("_x_")
        tsv_reader = csv.reader(f, delimiter='\t')
        for line in tsv_reader:
            query    = line[0]
            subject  = line[1]
            bitscore = float(line[-1])
            if query == subject:
                continue
            if query != last_query and query_file != subject_file:
                last_query = query
                all_hit_tuples[(query, subject)] = (bitscore, 0)
            else:
                all_hit_tuples[(query, subject)] = (bitscore, 1)
    return all_hit_tuples


def extract_hit_bitscore(blast_results, cpu=multiprocessing.cpu_count()):
    from contextlib import closing
    ori2num        = {}
    count          = 0
    all_hit_tuples = {}
    with closing(multiprocessing.get_context('spawn').Pool(cpu)) as pool:
        for hits in pool.imap_unordered(get_hit_bitscore, blast_results):
            hits_num = {}
            for key, valeus in hits.items():
                query, subject = key
                try:
                    query = ori2num[query]
                except KeyError:
                    ori2num[query] = str(count)
                    query          = str(count)
                    count          += 1
                try:
                    subject = ori2num[subject]
                except KeyError:
                    ori2num[subject]  = str(count)
                    subject           = str(count)
                    count            += 1
                hits_num[(query, subject)] = valeus
            all_hit_tuples.update(hits_num)
    return all_hit_tuples, ori2num


class NormalizedClustering():
    orthologs           = []
    non_ortholog_groups = []
    called              = 0

    def __init__(self, fasdir, indir, outdir, last_result, inflation=1.5, otu=4, cpu=multiprocessing.cpu_count()):
        NormalizedClustering.called              += 1
        NormalizedClustering.orthologs            = []
        NormalizedClustering.non_ortholog_groups  = []
        self.fasdir                     = Path(fasdir)
        self.indir                      = Path(indir)
        self.base                       = Path(outdir)
        self.last                       = Path(last_result)
        self.otu                        = otu
        self.cpu                        = cpu
        self.logdir                     = self.base / "log"
        self.blast_result_dir           = self.logdir / str(inflation)
        self.all_orthologs_pickle       = self.logdir / "ortho.bco"
        self.non_ortholog_groups_pickle = self.logdir / "leftovers.bi"
        self.outdir                     = self.base / "result_{0}".format(inflation)
        self.orthologs                  = []
        self.non_ortholog_groups        = []
        self.inflation                  = inflation
        self.logdir.mkdir(parents=True, exist_ok=True)
        self.blast_result_dir.mkdir(parents=True, exist_ok=True)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.init_pickle(self.all_orthologs_pickle, [])
        self.init_pickle(self.non_ortholog_groups_pickle, [])
        self.get_last_result()
        self.generate_clustering_input_with_dask()

    @staticmethod
    def exclude_orthologs_with_pandas(argments):
        blast_result = argments[0]
        outdir = argments[1]
        blast_result = Path(blast_result)
        df = pd.read_csv(str(blast_result), names=('query', 'subject', 'qs', 'qe', 'ss', 'se', 'nbit'), sep="\t",
                         header=None)
        outdir = Path(outdir)
        values = NormalizedClustering.orthologs
        df = df.query('query not in @values & subject not in @values')
        outfile = outdir / blast_result.name
        if outfile.exists():
            return
        df.to_csv(str(outfile), index=False, sep="\t", header=None)

    @staticmethod
    def exclude_orthologs(argments):
        blast_result = argments[0]
        outdir = argments[1]
        blast_result = Path(blast_result)
        outdir = Path(outdir)
        outfile = outdir / blast_result.name
        excluded_blast = ""
        if outfile.exists():
            return
        with open(str(blast_result), "r") as f, open(str(outfile), "w", encoding="UTF-8") as o:
            all_balst_result_lines = f.read().splitlines()
            for line in all_balst_result_lines:
                line = line.rstrip()
                line_eles = line.split("\t")
                query = line_eles[0]
                subject = line_eles[1]
                if (query in NormalizedClustering.orthologs) or (subject in NormalizedClustering.orthologs):
                    continue
                line = line + "\n"
                excluded_blast += line
            o.write(excluded_blast)

    @staticmethod
    def generate_pandas_frame(normalized_result):
        try:
            df = dd.read_csv(str(normalized_result), names=('query', 'subject', 'qs', 'qe', 'ss', 'se', 'nbit'),
                             sep="\t", blocksize=None)
        except:
            df = dd.read_csv(str(normalized_result), names=('query', 'subject', 'nbit'), sep='\t', blocksize=None)
        return df, normalized_result

    def generate_clustering_input_with_dask(self):
        blast_results = self.indir.glob("*.norm")
        args_list = [blast for blast in blast_results]
        if NormalizedClustering.called == 1:
            used_columns = ['query', 'subject', 'nbit']
        else:
            used_columns = ['query', 'subject', 'qs']
        with multiprocessing.get_context('spawn').Pool(self.cpu) as pool:
            for df, blast_result in pool.imap_unordered(NormalizedClustering.generate_pandas_frame, args_list):
                df           = df.compute()
                blast_result = Path(blast_result)
                outdir       = Path(self.blast_result_dir)
                values       = NormalizedClustering.orthologs
                df           = df.query('query not in @values & subject not in @values')
                outfile      = outdir / blast_result.name
                if outfile.exists():
                    return
                df.to_csv(str(outfile), index=False, sep="\t", header=None, columns=used_columns)


    def init_pickle(self, object_file, init_value=None):
        object_file = Path(object_file)
        if object_file.exists() and NormalizedClustering.called == 1:
            object_file.unlink()
        if object_file.exists():
            return
        if init_value is None:
            init_value = []
        with open(str(object_file), "wb") as bo:
            pickle.dump(init_value, bo)

    def get_last_result(self):
        last_orthologs = self.last / "log" / "ortho.bco"
        if last_orthologs.exists() and self.last != Path(""):
            with open(str(last_orthologs), "rb") as bf:
                orthologs = pickle.load(bf)
                orthologs = list(set(orthologs))
                NormalizedClustering.orthologs.extend(orthologs)
        with open(str(self.all_orthologs_pickle), "rb") as bf:
            all_orthologs = pickle.load(bf)
            NormalizedClustering.orthologs.extend(all_orthologs)

        last_non_ortholog_groups = self.last / "log" / "leftovers.bi"
        if last_non_ortholog_groups.exists() and self.last != Path(""):
            with open(str(last_non_ortholog_groups), "rb") as bf:
                non_orthologs = pickle.load(bf)
                non_orthologs = [non_ortholog for non_ortholog in non_orthologs
                                 if non_ortholog not in NormalizedClustering.non_ortholog_groups]
                NormalizedClustering.non_ortholog_groups.extend(non_orthologs)
        with open(str(self.non_ortholog_groups_pickle), "rb") as bf:
            all_orthologs = pickle.load(bf)
            NormalizedClustering.non_ortholog_groups.extend(all_orthologs)

        with open(str(self.all_orthologs_pickle), "wb") as orthologs, open(str(self.non_ortholog_groups_pickle),
                                                                           "wb") as non_orthologs:
            pickle.dump(NormalizedClustering.orthologs, orthologs)
            pickle.dump(NormalizedClustering.non_ortholog_groups, non_orthologs)

    def rbh_clustering(self):
        rbh = RBH(fasta_dir=self.fasdir, blast_result_dir=self.blast_result_dir, outdir=self.outdir, otu=self.otu)
        rbh.extract_hit_bitscore()
        rbh.extract_rbh()
        rbh.calc_min_rbh_bitscore()
        rbh.append_bbh()
        rbh.generate_mcl_input()
        rbh.exec_mcl()
        rbh.mcl_groups()
        rbh.make_mcl_group_fasta()


class RBH():
    Seq2Dict = defaultdict(str)

    def __init__(self, fasta_dir, blast_result_dir, outdir, inflation=1.5, cpu=multiprocessing.cpu_count(),
                 otu=4):
        self.fasta_dir        = Path(fasta_dir)
        self.blast_dir        = Path(blast_result_dir)
        self.outdir           = Path(outdir)
        self.otu              = otu
        self.besthit_bitscore = {}
        self.hit_bitscore     = {}
        self.rbhs             = []
        self.rbh_groups       = []
        self.cpu              = cpu
        self.inflation        = inflation
        self.ori2num          = {}
        self.group_id2seq_ids = defaultdict(list)
        self.ortholog_dir     = self.outdir / "groups"
        self.mcl_wd           = self.outdir / "mcl"
        self.mcl_file         = self.mcl_wd / "mcl.mcl"
        self.mcl_out          = self.mcl_wd / "mcl_result.mcl"
        self.min_rbh_bitscore = defaultdict(lambda: float("inf"))
        self.query2besthits   = defaultdict(list)
        self.rbh2score        = {}
        self.all_hits         = {}
        self.outdir.mkdir(exist_ok=True, parents=True)
        self.mcl_wd.mkdir(exist_ok=True, parents=True)
        self.ortholog_dir.mkdir(exist_ok=True, parents=True)
        self.make_seq2dict()

    @staticmethod
    def get_hit_bitscore(blast_result):
        last_query = ""
        besthit_bitscore = {}
        hit_bitscore = {}
        all_hits = {}
        with open(str(blast_result), "r") as f:
            res = blast_result.stem
            query_file, subject_file = res.split("_x_")
            for line in f:
                line = line.rstrip()
                line_ele = line.split("\t")
                query = line_ele[0]
                subject = line_ele[1]
                bitscore = float(line_ele[-1])
                if query == subject:
                    continue
                all_hits[(query, subject)] = bitscore
                if query != last_query and query_file != subject_file:
                    last_query = query
                    besthit_bitscore[(query, subject)] = bitscore
                else:
                    hit_bitscore[(query, subject)] = bitscore
        return besthit_bitscore, hit_bitscore, all_hits

    def extract_hit_bitscore(self):
        self.ori2num       = {}
        blast_results      = self.blast_dir.glob("*.*")
        cpu                = self.cpu
        hits, self.ori2num = extract_hit_bitscore(blast_results, cpu=cpu)
        for key in sorted(hits.keys()):
            query, subject  = key
            bit, identifier = hits[key]
            if identifier == 0:
                self.besthit_bitscore[(query, subject)] = bit
                self.hit_bitscore[(query, subject)] = bit
            else:
                self.hit_bitscore[(query, subject)] = bit
            del hits[key]

    def extract_rbh(self):
        best_hits = self.besthit_bitscore.keys()
        for hit in best_hits:
            query, subject = hit
            reverse_hit    = (subject, query)
            if reverse_hit in best_hits:
                self.rbhs.append(hit)
                del self.hit_bitscore[hit]
        del best_hits

    def calc_min_rbh_bitscore(self):
        for rbh in self.rbhs:
            query, subject = rbh
            bitscore       = self.besthit_bitscore[(query, subject)]
            if bitscore < self.min_rbh_bitscore[query]:
                self.min_rbh_bitscore[query] = bitscore

    def append_bbh(self):
        for hit in list(self.hit_bitscore.keys()):
            bitscore       = self.hit_bitscore[hit]
            query, subject = hit
            min_bitscore   = self.min_rbh_bitscore[query]
            if bitscore >= min_bitscore:
                self.rbhs.append(hit)
                self.besthit_bitscore[(query, subject)] = bitscore
                del self.hit_bitscore[(query, subject)]


    @staticmethod
    def output_fasta_files(argments):
        group, outfile = argments
        fasta = ""
        for seq_id in group:
            seq = RBH.Seq2Dict[seq_id]
            ori_seq_id = seq_id
            fas = ">{0}\n{1}\n".format(ori_seq_id, seq)
            fasta += fas
        with open(str(outfile), "w", encoding="UTF-8") as out:
            out.write(fasta)

    @staticmethod
    def output_mcl_fasta_files(argments):
        fasta = ""
        group, outfile = argments
        for seq_id in group:
            seq         = RBH.Seq2Dict[seq_id]
            ori_seq_id  = seq_id
            fas         = ">{0}\n{1}\n".format(ori_seq_id, seq)
            fasta      += fas
        with open(str(outfile), "w", encoding="UTF-8") as out:
            out.write(fasta)

    def make_seq2dict(self):
        fastas = self.fasta_dir.glob("**/*.*")
        fastas = [str(fasta) for fasta in fastas]
        last_id = ""
        with fileinput.input(files=fastas) as f:
            for line in f:
                line = line.rstrip()
                if line.startswith(">"):
                    last_id = line[1:]
                else:
                    RBH.Seq2Dict[last_id] += line

    def make_group_fasta(self):
        ortholog_dir = self.ortholog_dir
        ortholog_dir.mkdir(exist_ok=True, parents=True)
        groups = len(self.rbh_groups)
        digit = len(str(groups))

        args = []
        cpu = self.cpu
        for index, group in enumerate(self.rbh_groups):
            group_id = str(index).zfill(digit)
            gruop_file = "OG{0}.fa".format(group_id)
            outfile = ortholog_dir / gruop_file
            args.append((group, outfile))
        with multiprocessing.Pool(cpu) as pool:
            for _ in pool.imap_unordered(RBH.output_fasta_files, args):
                pass

    def generate_mcl_input(self):
        mcl      = self.mcl_file
        besthits = set()
        for query, subject in self.rbhs:
            self.query2besthits[query].append(subject)
            besthits.add(query)
            besthits.add(subject)
        besthits = sorted(besthits)
        if not self.ori2num:
            self.ori2num = {ori_name: str(new_name) for new_name, ori_name in enumerate(besthits)}
        besthits_len = len(self.ori2num.keys())
        mcl_header = "(mclheader\nmcltype matrix\ndimensions {0}x{0}\n)\n\n(mclmatrix\nbegin\n\n".format(
            besthits_len)
        seq2probs = defaultdict(list)

        for query, num in sorted(self.ori2num.items()):
            query = num
            hits  = self.query2besthits[query]
            for hit in hits:
                subject = hit
                try:
                    score = self.hit_bitscore[(query, subject)]
                except KeyError:
                    score = self.besthit_bitscore[(query, subject)]
                query_num   = query
                subject_num = subject
                prob        = "{0}:{1:.3f}".format(subject_num, score)
                seq2probs[query_num].append(prob)
        del self.hit_bitscore
        del self.besthit_bitscore
        if self.mcl_file.exists():
            return
        with open(str(mcl), "w", encoding="UTF-8") as o:
            o.write(mcl_header)
            for query, hits in sorted(seq2probs.items(),
                                      key=lambda x: (int(x[0]), x[1])):
                line = "{0}\t{1}$".format(query, " ".join(hits))
                print(line, file=o)
            print(")", file=o)

    def exec_mcl(self):
        cmd = "mcl {0} -I {1} -t {2} -o {3}".format(str(self.mcl_file), str(self.inflation), self.cpu,
                                                    str(self.mcl_out))
        if not self.mcl_out.exists():
            subprocess.run(cmd, shell=True)

    def mcl_groups(self):
        mcl_out = self.mcl_out
        with open(str(mcl_out), "r", encoding="UTF-8") as f:
            while True:
                line = next(f)
                if line.startswith("begin"):
                    break
            mcl_group_line = ""
            for line in f:
                if line.startswith(")"):
                    break
                line = line.rstrip()
                mcl_group_line += line
                if "$" in line:
                    ids            = re.split(r"\s+", mcl_group_line)
                    mcl_group_line = ""
                    seq_ids        = ids[1:-1]
                    group_id       = ids[0]
                    self.group_id2seq_ids[group_id] = seq_ids

    def make_mcl_group_fasta(self):
        ortholog_dir = self.ortholog_dir
        ortholog_dir.mkdir(exist_ok=True, parents=True)
        groups = len(self.group_id2seq_ids.keys())
        digit  = len(str(groups))

        args    = []
        cpu     = self.cpu
        num2ori = {num : ori for ori, num in self.ori2num.items()}
        for group_id, seq_ids in self.group_id2seq_ids.items():
            group = [num2ori[seq_id] for seq_id in seq_ids]
            if (len(group) < self.otu) or (set(group) in NormalizedClustering.non_ortholog_groups):
                continue
            group_id   = str(group_id).zfill(digit)
            gruop_file = "OG{0}.fa".format(group_id)
            outfile    = ortholog_dir / gruop_file
            if not outfile.exists():
                args.append((group, outfile))
        with multiprocessing.Pool(cpu) as pool:
            for _ in pool.imap_unordered(RBH.output_mcl_fasta_files, args):
                pass