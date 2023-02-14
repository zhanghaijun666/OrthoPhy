#!/usr/bin/python3

import multiprocessing
import os
from pathlib import Path
import re
import shutil
import subprocess
from Bio import SeqIO


def exec_hgt_prediction(instance):
    instance.exec_dimob()
    no_hgt = instance.extract_hgt_protein()
    return no_hgt


class Acceptor():
    def __init__(self, indir, outdir, hgt_dir="", input_type="", cpu=1):
        self.indir    = Path(indir)
        self.base_dir = self.indir.resolve().name
        self.infiles  = sorted(Path(self.indir).glob("**/*.*"))
        self.out      = Path(outdir)
        self.suffix   = ".faa"
        self.hgt_dir  = Path(hgt_dir)
        self.dirtype  = input_type
        self.cpu      = cpu
        self.spe2seq  = {}

    def generate_orthologfinder_input(self):
        self.make_input()
        self.taxonomy()
        return self.spe2seq

    def make_input(self):
        if self.dirtype == "GH":
            self.predict_hgt()
        elif self.dirtype == "GN":
            self.extract_cds()
        elif self.dirtype == "FA":
            self.copy_fasta()
        return self.spe2seq

    def taxonomy(self):
        fasta_dir = Path(self.out).resolve()
        fastas    = sorted(fasta_dir.glob("*.*"))
        sources   = self.infiles
        base      = self.indir.resolve()
        for fasta in fastas:
            for source in sources:
                is_same = self.is_same_prefix(source, fasta)
                if is_same:
                    abs_path   = source.resolve()
                    relativity = abs_path.relative_to(base)
                    dest       = fasta_dir / relativity
                    dest_dir   = Path(dest.parent)
                    dest_path  = dest_dir / str(fasta.name)
                    if dest_path.exists():
                        fasta.unlink()
                    else:
                        dest_dir.mkdir(parents=True, exist_ok=True)
                        shutil.move(str(fasta), str(dest_path))
                    break

    def is_same_prefix(self, p1, p2):
        name1 = p1.stem
        name2 = p2.stem
        name1 = re.sub(r'\W', "_", name1)
        name2 = re.sub(r'\W', "_", name2)
        if name1 == name2:
            return True
        else:
            return False

    def copy_fasta(self):
        indir = Path(self.indir)
        outdir = Path(self.out)
        for fasta in indir.glob("**/*.*"):
            species = fasta.stem
            species = re.sub(r'\W', "_", species)
            dest = outdir / "{0}.faa".format(species)
            with open(str(dest), "w", encoding="UTF-8") as o:
                for record in SeqIO.parse(str(fasta), "fasta"):
                    header = record.description
                    # replace invalid chars to _
                    header = re.sub(r'\W', "_", header)
                    seq = str(record.seq)
                    self.spe2seq[header] = species
                    new_rec = ">{0}\n{1}".format(header, seq)
                    print(new_rec, file=o)
        return self.spe2seq

    def predict_hgt(self):
        outdir        = self.hgt_dir
        predictors    = (LGTFinder(gbk, outdir) for gbk in self.infiles)
        cpu           = self.cpu
        no_hgt_fastas = set()
        with multiprocessing.Pool(cpu) as pool:
            for no_hgt in pool.imap_unordered(exec_hgt_prediction, predictors):
                no_hgt_fastas.add(no_hgt)
        dest = self.out
        for nohgt in no_hgt_fastas:
            species = nohgt.parent.name
            suffix = self.suffix
            outfile = dest / "{0}{1}".format(species, suffix)
            shutil.copy(str(nohgt), str(outfile))

        fasta_dir = Path(dest)
        for fasta in fasta_dir.glob("*.*"):
            species = fasta.stem
            # The name of input files must contain only alphabet and numbers.
            species = re.sub(r'\W', "_", species)
            with open(str(fasta), "r", encoding="UTF-8") as o:
                for record in SeqIO.parse(str(fasta), "fasta"):
                    header = record.description
                    # replace invalid chars to _
                    header = re.sub(r'\W', "_", header)
                    self.spe2seq[header] = species
        return self.spe2seq

    def extract_cds(self):
        for gbk in self.infiles:
            species = gbk.stem
            species = re.sub(r'\W', "_", species)
            suffix = self.suffix
            fasta = species + suffix
            outpath = str(self.out / fasta)
            all_cds = 0
            with open(outpath, "w", encoding="UTF-8") as outfile:
                for record in SeqIO.parse(str(gbk), "genbank"):
                    for feature in record.features:
                        if not feature.type == "CDS":
                            continue
                        try:
                            seq = feature.qualifiers["translation"][0]
                        except KeyError:
                            continue
                        all_cds += 1
                        prot_id = feature.qualifiers["protein_id"][-1]
                        anno = prot_id
                        anno = re.sub(r'\W', "_", anno)
                        self.spe2seq[anno] = species
                        print(">{0}\n{1}".format(anno, seq), file=outfile)
                if all_cds == 0:
                    raise ValueError(
                        "The following file has no CDS information: {0}".format(
                            gbk))
        return self.spe2seq


class LGTFinder():
    def __init__(self, gbk, out):
        self.gbk       = Path(gbk)
        self.species   = self.gbk.stem
        self.species   = re.sub(r'\W', "_", self.species)
        self.outdir    = Path(out)      / self.species
        self.dimob_dir = self.outdir    / "IslandPath-DIMOB"
        self.no_hgt    = self.outdir    / "no_hgt.faa"
        self.hgt       = self.outdir    / "hgt.faa"
        self.summary   = self.outdir    / "summary.txt"
        self.dimob_out = self.dimob_dir / "{0}_GIs.txt".format(self.species)
        self.GIs       = set()
        self.dimob_dir.mkdir(parents=True, exist_ok=True)
        self.init_output()

    def init_output(self):
        outfiles = [self.no_hgt, self.hgt, self.summary]
        for outfile in outfiles:
            outfile.touch()

    def get_prefix(self, path):
        base   = os.path.basename(path)
        prefix = os.path.splitext(base)[0]
        return prefix

    def exec_dimob(self):
        infile  = str(self.gbk)
        outfile = str(self.dimob_out)
        log     = str(self.dimob_dir / "Dimob.log")
        cmd     = "Dimob.pl {0} {1} > {2} 2>&1".format(infile, outfile, log)
        if not self.dimob_out.exists():
            subprocess.call(cmd, shell=True)
        tmpfile = Path("./Dimob.log")
        if tmpfile.exists():
            os.remove(str(tmpfile))

    def GI_location(self):
        dimob_out = str(self.dimob_out)
        with open(dimob_out, "r", encoding="UTF-8") as f:
            for line in f:
                if not line:
                    continue
                line = line.rstrip()
                _, first, last = line.split("\t")
                self.GIs.add(range(int(first), int(last) + 1))

    def extract_hgt_protein(self):
        self.GI_location()
        if self.no_hgt.stat().st_size != 0:
            return self.no_hgt
        gbk           = str(self.gbk)
        hgt_out       = str(self.hgt)
        nohgt_out     = str(self.no_hgt)
        summary_out   = str(self.summary)
        record        = SeqIO.read(gbk, "genbank")
        features      = record.features
        all_cds       = 0
        hgt_prots     = 0
        non_hgt_prots = 0
        with open(hgt_out, "w", encoding="UTF-8") as hgt, open(nohgt_out, "w", encoding="UTF-8") as no_hgt:
            for feature in features:
                if feature.type != "CDS":
                    continue
                try:
                    trans = feature.qualifiers["translation"][0]
                except KeyError:
                    continue

                all_cds += 1
                prot_id  = feature.qualifiers["protein_id"][0]
                anno     = prot_id
                # replace invalid chars to _
                anno     = re.sub(r'\W', "_", anno)
                location = feature.location
                try:
                    start = location.parts[0].start.position + 1
                    end   = location.parts[-1].end.position
                except AttributeError:
                    non_hgt_prots += 1
                    print(">{0}\n{1}".format(anno, trans), file=no_hgt)
                else:
                    if self.is_in_genomic_island(start, end):
                        hgt_prots += 1
                        print(">{}\n{}".format(anno, trans), file=hgt)
                    else:
                        non_hgt_prots += 1
                        print(">{0}\n{1}".format(anno, trans), file=no_hgt)

        with open(summary_out, "w", encoding="UTF-8") as f:
            hgt_per     = hgt_prots / all_cds * 100
            non_hgt_per = non_hgt_prots / all_cds * 100
            print("ALL_sequence\t{0}".format(all_cds), file=f)
            print("HGT_sequence\t{0}({1:.2f}%)".format(hgt_prots, hgt_per), file=f)
            print("NON-HGT_sequence\t{0}({1:.2f}%)".format(non_hgt_prots, non_hgt_per), file=f)
        return self.no_hgt

    def is_in_genomic_island(self, start, end):
        prot_reg = range(start, end)
        for gi_reg in self.GIs:
            is_overlap = set(gi_reg) & set(prot_reg)
            if is_overlap:
                return True
        return False
