#!/usr/bin/python3

from collections import defaultdict
import multiprocessing
import subprocess
import sys
import shutil
from pathlib import Path
from tqdm import tqdm
from Bio import SeqIO

seq2spe = {}
NO_FILTERING = False
LOG_FILE = ""
duplicates = 0


def exec_reconstruct(instance):
    instance.auto_exec()


class TreeReconstructor():
    def __init__(self, infile, align_out, trim_out, tree_out, dup_out):
        self.seq          = Path(infile)
        self.mafft_dir    = Path(align_out)
        self.trimal_dir   = Path(trim_out)
        self.fasttree_dir = Path(tree_out)
        self.prefix       = self.seq.stem
        self.mafft_out    = self.mafft_dir / "{0}.faa".format(self.prefix)
        self.trimal_out   = self.trimal_dir / "{0}.trim".format(self.prefix)
        self.fasttree_out = self.fasttree_dir / "{0}.nwk".format(self.prefix)
        self.dup_out      = dup_out / "{0}.tsv".format(self.prefix)

    def auto_exec(self):
        self.no_duplicate_mafft()
        if NO_FILTERING:
            self.trimal_out = self.mafft_out
        else:
            self.trimal()
        self.fasttree()

    def no_duplicate_mafft(self):
        if self.mafft_out.exists():
            return
        infile       = str(self.seq)
        outfile      = str(self.mafft_out)
        dup_seqs_tsv = str(self.dup_out)
        cmd          = "seqkit rmdup -s --quiet {0} -D {1} | mafft --anysymbol --quiet --auto - > {2}".format(infile, dup_seqs_tsv, outfile)
        subprocess.call(cmd, shell=True)

    def mafft(self):
        if self.mafft_out.exists():
            return
        infile = str(self.seq)
        outfile = str(self.mafft_out)
        cmd = "mafft --anysymbol --quiet --auto {0} > {1}".format(infile, outfile)
        subprocess.call(cmd, shell=True)

    def trimal(self):
        if self.trimal_out.exists():
            return
        infile = str(self.mafft_out)
        outfile = str(self.trimal_out)
        cmd = "trimal -automated1 -in {0} -out {1} > /dev/null 2>&1".format(infile, outfile)
        subprocess.call(cmd, shell=True)
        cmd = 'grep -v -c ">" {0}'.format(outfile)
        non_header_lines = subprocess.getoutput(cmd)
        if non_header_lines == "0":
            shutil.copy(infile, outfile)

    def fasttree(self):
        if self.fasttree_out.exists():
            return
        infile = str(self.trimal_out)
        outfile = str(self.fasttree_out)
        cmd = "fasttree -lg -quiet -nosupport {0} > {1} 2> /dev/null".format(infile, outfile)
        sc = subprocess.call(cmd, shell=True)
        if sc != 0:
            cmd = "fasttree -quiet -nosupport {0} > {1} 2> /dev/null".format(infile, outfile)
            sc = subprocess.call(cmd, shell=True)
            if sc != 0:
                Path(self.fasttree_out).unlink()


class GeneTreeReconstructor():
    LOG_FILE = ""

    def __init__(self, indir, align_out, trim_out, tree_out, cpu, log, seq_to_spe, no_filtering=False, dup_out=None):
        self.cpu          = int(cpu)
        self.indir        = Path(indir)
        self.mafft_out    = Path(align_out)
        self.trimal_out   = Path(trim_out)
        self.fasttree_out = Path(tree_out)
        self.dup_out      = Path(dup_out)
        self.outdirs      = [self.mafft_out, self.trimal_out, self.fasttree_out, self.dup_out]
        self.make_outdirs()
        self.log = str(log)
        GeneTreeReconstructor.LOG_FILE = self.log

        global seq2spe, NO_FILTERING
        seq2spe      = seq_to_spe
        NO_FILTERING = no_filtering

    @staticmethod
    def logging(msg):
        print(msg)
        with open(GeneTreeReconstructor.LOG_FILE, "a", encoding="UTF-8") as out:
            print(msg, file=out)

    def make_outdirs(self):
        for sub_dir in self.outdirs:
            sub_dir.mkdir(parents=True, exist_ok=True)

    def gene_trees(self):
        cpu = self.cpu
        fastas = sorted(self.indir.glob("*.*"))
        self.logging("mcl clustering groups: {0}".format(len(fastas)))

        reconstructors = [TreeReconstructor(infile, self.mafft_out, self.trimal_out, self.fasttree_out, self.dup_out)
                          for infile in fastas]
        all_trees = len(reconstructors)
        self.logging("reconstruct {0} trees".format(str(all_trees)))
        cpu = 9
        with multiprocessing.Pool(cpu) as pool:
            for _ in tqdm(pool.imap_unordered(exec_reconstruct, reconstructors), total=all_trees):
                pass

    def count_otu(self, fasta):
        cmd = 'grep -c ">" {0}'.format(fasta)
        count = subprocess.getoutput(cmd)
        count = int(count)
        return count

    def multicore_reconstruction(self):
        cpu    = self.cpu
        fastas = sorted(self.indir.glob("*.*"))
        self.logging("mcl clustering groups: {0}".format(len(fastas)))
        reconstructors = [TreeReconstructor(infile, self.mafft_out, self.trimal_out, self.fasttree_out, self.dup_out) for infile in
                          fastas]
        all_trees = len(reconstructors)

        self.logging("reconstruct {0} trees".format(str(all_trees)))
        #with multiprocessing.Pool(cpu) as pool:
        with multiprocessing.Pool(9) as pool:
            for _ in tqdm(pool.imap_unordered(exec_reconstruct, reconstructors), total=all_trees):
                pass
