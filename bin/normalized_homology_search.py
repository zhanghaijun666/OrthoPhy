#!/usr/bin/python3

from collections import defaultdict
import itertools
import math
import multiprocessing
import os
from pathlib import Path
from scipy.stats import linregress
import subprocess


class NormalizedHomologySearch():
    def __init__(self, indir, outdir, coverage=50,
                 cpu=multiprocessing.cpu_count()):
        self.indir              = Path(indir)
        self.outdir             = Path(outdir)
        self.coverage           = coverage
        self.cpu                = cpu
        self.converted_fastadir = self.outdir / "fasta"
        self.normalized_dir     = self.outdir / "normalized"
        self.blast_result       = self.outdir / "diamond_results"

    def normalized_homology_search(self):
        diamond = SpeciesToSpeciesHGS(self.indir, self.outdir, self.coverage)
        diamond.make_diamond_db()
        diamond.all_vs_all_diamond()
        norm = Normalized(fasdir=self.indir, indir=self.blast_result,
                          outdir=self.normalized_dir, cpu=self.cpu)
        norm.exec_multiprocess_norm()


class SpeciesToSpeciesHGS():
    def __init__(self, fasta_dir, out_dir, coverage=50,
                 cpu=multiprocessing.cpu_count()):
        self.fasta_dir        = Path(fasta_dir)
        self.fastas           = sorted(self.fasta_dir.glob("**/*.*"))
        self.hgs_result_dir   = Path(out_dir) / "diamond_results"
        self.db_dir           = Path(out_dir) / "diamond_database"
        self.coverage         = coverage
        self.cpu              = cpu
        self.hgs_result_dir.mkdir(exist_ok=True, parents=True)
        self.fasta_dir.mkdir(exist_ok=True, parents=True)
        self.db_dir.mkdir(exist_ok=True, parents=True)

    def make_diamond_db(self):
        fastas = self.fasta_dir.glob("**/*.*")
        for fasta in fastas:
            species = fasta.stem
            db_file = self.db_dir / species
            cmd = 'diamond makedb --quiet --in "{0}" --db "{1}"'.format(str(fasta),
                                                                str(db_file))
            db_path = db_file.with_suffix(".dmnd")
            if not db_path.exists():
                subprocess.run(cmd, shell=True)

    @staticmethod
    def execute_diamond(cmds):
        subprocess.run(cmds, shell=True)

    def all_vs_all_diamond(self):
        queries = self.fasta_dir.glob("**/*.*")
        databases = sorted(self.db_dir.glob("*.*"))
        cmds = []
        for query, database in itertools.product(queries, databases):
            query_species = query.stem
            database_species = database.stem
            outfile = self.hgs_result_dir / "{0}_x_{1}.diamond".format(
                query_species, database_species)
            if outfile.exists():
                continue
            cmd = 'diamond blastp --quiet --threads 1 --more-sensitive --max-hsps 1 -q "{0}" -d "{1}" -o "{2}" --query-cover  {3} --subject-cover {3} ' \
                  '--outfmt 6 qseqid sseqid qstart qend sstart send bitscore'.format(
                str(query), str(database), str(outfile), str(self.coverage))
            cmd = 'diamond blastp --quiet --threads 1 --more-sensitive --max-hsps 1 -q "{0}" -d "{1}" -o "{2}" --query-cover  {3} --subject-cover {3} ' \
                  '--outfmt 6 qseqid sseqid qstart qend sstart send bitscore'.format(
                str(query), str(database), str(outfile), str(25))
            cmds.append(cmd)
        with multiprocessing.Pool(self.cpu) as pool:
            for _ in pool.imap_unordered(SpeciesToSpeciesHGS.execute_diamond, cmds):
                pass


class NormalizedBlast():

    def __init__(self, blast, outfile, cpu=multiprocessing.cpu_count()):
        self.blast       = Path(blast)
        self.outfile     = Path(outfile)
        self.slope       = 0.0
        self.intercept   = 0.0
        self.log_len2bit = []
        self.cpu         = cpu

    def generate_len2bit_with_log(self):
        self.log_len2bit = []
        with open(str(self.blast), "r") as f:
            for line in f:
                line          = line.rstrip()
                line          = line.split("\t")
                bit_score     = float(line[-1])
                log_bit_score = math.log10(bit_score)
                query_start   = int(line[2])
                query_end     = int(line[3])
                query_len     = abs(query_start - query_end) + 1
                subject_start = int(line[4])
                subject_end   = int(line[5])
                subject_len   = abs(subject_start - subject_end) + 1
                hit_len       = query_len * subject_len
                log_hit_len   = math.log10(hit_len)
                self.log_len2bit.append((log_hit_len, log_bit_score))
        return self.log_len2bit

    def normalized_parameter(self):
        bins = defaultdict(list)
        self.log_len2bit.sort()
        all_hits = len(self.log_len2bit)
        if all_hits >= 5000:
            bin_size = 1000
        else:
            bin_size = 200

        for index, xy in enumerate(self.log_len2bit):
            counter = index // bin_size
            bins[counter].append(list(xy))

        top_5percents = []
        for count, xy in bins.items():
            sort_by_bits = sorted(xy, key=lambda x: (x[1], x[0]), reverse=True)
            tops_num     = int(math.floor(len(sort_by_bits) * 0.05))
            if tops_num < 10:
                tops_num = 10
            tops = sort_by_bits[:tops_num]
            top_5percents.extend(tops)

        slope, intercept, rvalue, pvalue, stderr = linregress(top_5percents)
        self.slope     = slope
        self.intercept = intercept

    def normalize_bitscore(self, bit_score, hit_len):
        normalized_bit = bit_score / (
                (10 ** self.intercept) * (hit_len ** self.slope))
        normalized_bit = round(normalized_bit, 6)
        return normalized_bit

    def write_nomarlized_blast(self):
        outfile          = self.outfile
        query_lines      = []
        normalized_blast = ""
        with open(str(self.blast), "r") as f, open(str(outfile), "w", encoding="UTF-8") as out:
            all_blast_lines = f.read().splitlines()
            line            = all_blast_lines.pop(0)
            line            = line.rstrip()
            line            = line.split("\t")
            query           = line[0]
            last_query      = query
            query_start     = int(line[2])
            query_end       = int(line[3])
            query_len       = abs(query_start - query_end) + 1
            subject_start   = int(line[4])
            subject_end     = int(line[5])
            subject_len     = abs(subject_start - subject_end) + 1
            hit_len         = query_len * subject_len
            bit_score       = float(line[-1])
            normalized_bitscore = self.normalize_bitscore(bit_score, hit_len)
            new_line = line[:-1] + [normalized_bitscore]
            query_lines.append(new_line)

            for line in all_blast_lines:
                line = line.rstrip()
                line = line.split("\t")
                query = line[0]
                if query != last_query:
                    query_lines = sorted(query_lines, key=lambda x: x[-1], reverse=True)
                    for ql in query_lines:
                        ql = list(map(str, ql))
                        ql = "\t".join(ql[:])
                        ql = ql + "\n"
                        normalized_blast += ql
                        # print(ql, file=out)
                    query_lines = []
                    last_query = query

                query_start         = int(line[2])
                query_end           = int(line[3])
                query_len           = abs(query_start - query_end) + 1
                subject_start       = int(line[4])
                subject_end         = int(line[5])
                subject_len         = abs(subject_start - subject_end) + 1
                hit_len             = query_len * subject_len
                bit_score           = float(line[-1])
                normalized_bitscore = self.normalize_bitscore(bit_score   , hit_len)
                new_line            = line[:-1] + [normalized_bitscore]
                query_lines.append(new_line)
            query_lines = sorted(query_lines, key=lambda x: x[-1], reverse=True)
            for ql in query_lines:
                ql = list(map(str, ql))
                ql = "\t".join(ql[:])
                ql = ql + "\n"
                normalized_blast += ql
            out.write(normalized_blast)


class Normalized():

    def __init__(self, fasdir, indir, outdir, cpu=multiprocessing.cpu_count()):
        self.indir                 = Path(indir)
        self.fasdir                = Path(fasdir)
        self.diamonds              = self.indir.glob("**/*.diamond")
        self.outdir                = Path(outdir)
        self.normalizers           = []
        self.cpu                   = cpu
        self.outdir.mkdir(exist_ok=True, parents=True)
        self.init_normalizers()

    def exec_multiprocess_norm(self):
        with multiprocessing.Pool(self.cpu) as pool:
            for _ in pool.imap_unordered(Normalized.multiprocess_norm,
                                         self.normalizers):
                pass

    @staticmethod
    def multiprocess_norm(instance: NormalizedBlast):
        instance.generate_len2bit_with_log()
        instance.normalized_parameter()
        instance.write_nomarlized_blast()

    def init_normalizers(self):
        for diamond in self.diamonds:
            infile = diamond
            file_size = os.path.getsize(str(infile))
            if file_size == 0:
                continue
            name    = infile.stem
            outfile = self.outdir / "{0}.norm".format(name)
            if outfile.exists():
                continue
            normalizer = NormalizedBlast(infile, outfile)
            self.normalizers.append(normalizer)
