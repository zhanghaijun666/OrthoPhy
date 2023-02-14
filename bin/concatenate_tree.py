#!/usr/bin/python3

from collections import defaultdict
import math
from pathlib import Path
import shutil
import subprocess
import multiprocessing
from Bio import SeqIO
import os


def multiple_alignment(instance):
    instance.auto_exec()

def multiple_gene_tree_reconstruction(instance):
    instance.auto_exec()
    instance.fasttree()



class MultipleAlignment():
    def __init__(self, infile, align_out, trim_out, tree_out=""):
        self.seq          = Path(infile)
        self.prefix       = self.seq.stem
        self.parent       = self.seq.parent.name
        self.mafft_dir    = Path(align_out)
        self.mafft_out    = self.mafft_dir / "{0}-{1}.aln".format(self.prefix         , self.parent)
        self.trimal_dir   = Path(trim_out)
        self.trimal_out   = self.trimal_dir / "{0}.trim".format(self.mafft_out.stem)
        self.tree_out_dir = Path(tree_out)
        self.tree_out     = self.tree_out_dir / "{0}.nwk".format(self.trimal_out.stem)

    def fasttree(self):
        if self.tree_out.exists():
            return
        infile  = str(self.trimal_out)
        outfile = str(self.tree_out)
        cmd     = "fasttree -lg -quiet -nosupport {0} > {1} 2> /dev/null".format(infile, outfile)
        subprocess.call(cmd, shell=True)

    def mafft(self):
        if self.mafft_out.exists():
            return
        infile  = str(self.seq)
        outfile = str(self.mafft_out)
        cmd     = "mafft --anysymbol --quiet --maxiterate 1000 --localpair  {0} > {1}".format(infile, outfile)
        subprocess.call(cmd, shell=True)

    def trimal(self):
        if self.trimal_out.exists():
            return
        infile  = str(self.mafft_out)
        outfile = str(self.trimal_out)
        cmd     = "trimal -automated1 -in {0} -out {1}".format(infile, outfile)
        subprocess.call(cmd, shell=True)

        cmd = 'grep -c -v ">" {0}'.format(outfile)
        non_header_lines = int(subprocess.getoutput(cmd))
        if non_header_lines == 0:
            shutil.copy(infile, outfile)

    def auto_exec(self):
        self.mafft()
        self.trimal()


class ConcatenateTree():
    def __init__(self, indir, outdir, species_dict, otu=4, log=None, share=0.5, cpu=1):
        self.indir          = Path(indir)
        self.outdir         = Path(outdir)
        self.align_out      = self.outdir / "mafft"
        self.trim_out       = self.outdir / "trimal"
        self.tree_out       = self.outdir / "tree"
        self.gene_trees_out = self.outdir / "ortholog_trees"
        self.conc_seq       = self.tree_out / "concatenate.trim"
        self.conc_tree      = self.tree_out / "concatenate.nwk"
        self.sub_dirs       = [self.align_out, self.trim_out, self.tree_out]
        self.orthologues    = sorted(self.indir.glob("**/*.ortholog"))
        self.otu            = otu
        self.share          = share
        self.log            = log
        self.seq2spe        = species_dict
        self.species        = self.get_all_species()
        self.cpu            = cpu
        self.make_outdirs()

    def astral(self):
        self.gene_trees_out.mkdir(parents=True, exist_ok=True)
        rate        = self.share
        species_num = len(self.species)
        threshold   = math.ceil(species_num * rate)
        if threshold < 4:
            threshold = 4
        orthologs = [ortholog for ortholog in self.orthologues
                     if self.count_otu(ortholog) >= threshold]
        self.ortholog_dir = self.outdir / "orthologs"
        for ortholog in orthologs:
            spe       = set()
            inflation = ortholog.parent.name
            outdir    = self.ortholog_dir / inflation
            outdir.mkdir(exist_ok=True, parents=True)
            outfile = outdir / (ortholog.stem + ".fas")
            with open(str(outfile), "w", encoding="UTF-8") as f:
                for record in SeqIO.parse(str(ortholog), "fasta"):
                    header  = record.description
                    species = self.seq2spe[header]
                    if species in spe:
                        continue
                    else:
                        spe.add(species)
                    seq = str(record.seq)
                    fas = ">{0}\n{1}".format(species, seq)
                    print(fas, file=f)
        fastas    = self.ortholog_dir.glob("**/*.fas")
        instances = [MultipleAlignment(seq, self.align_out, self.trim_out, self.gene_trees_out) for seq in fastas]
        cpu       = self.cpu
        with multiprocessing.Pool(cpu) as pool:
            for _ in pool.imap_unordered(multiple_gene_tree_reconstruction, instances):
                pass
        all_gene_trees     = self.tree_out.resolve() / "gene_trees.nwk"
        cmd                = f"cat {self.gene_trees_out.resolve()}/*.nwk >> {all_gene_trees}"
        subprocess.run(cmd, shell=True)
        cwd                = Path.cwd()
        home_dir           = Path(os.environ['HOME'])
        astral_dir         = home_dir / "bin" / "Astral"
        conc_tree_out      = str(self.conc_tree.resolve())
        os.chdir(str(astral_dir))
        cmd = f'java -D"java.library.path=lib/" -jar astral.5.15.4.jar -i {all_gene_trees} -o {conc_tree_out}'
        subprocess.run(cmd, shell=True)
        os.chdir(str(cwd))

    def iqtree2(self):
        rate        = self.share
        species_num = len(self.species)
        threshold   = math.ceil(species_num * rate)
        orthologs   = [ortholog for ortholog in self.orthologues if self.count_otu(ortholog) >= threshold]
        self.ortholog_dir = self.outdir / "orthologs"
        for ortholog in orthologs:
            spe = set()
            inflation = ortholog.parent.name
            outdir    = self.ortholog_dir / inflation
            outdir.mkdir(exist_ok=True, parents=True)
            outfile = outdir / (ortholog.stem + ".fas")
            with open(str(outfile), "w", encoding="UTF-8") as f:
                for record in SeqIO.parse(str(ortholog), "fasta"):
                    header  = record.description
                    species = self.seq2spe[header]
                    if species in spe:
                        continue
                    else:
                        spe.add(species)
                    seq = str(record.seq)
                    fas = ">{0}\n{1}".format(species, seq)
                    print(fas, file=f)
        fastas    = self.ortholog_dir.glob("**/*.fas")
        instances = [MultipleAlignment(seq, self.align_out, self.trim_out) for seq in fastas]
        cpu       = self.cpu
        with multiprocessing.Pool(cpu) as pool:
            for _ in pool.imap_unordered(multiple_alignment, instances):
                pass
        prefix = "tree"
        cmd    = ["iqtree2", "-Q", str(self.trim_out), "-bb", "1000",
                    "-pre", prefix, "-T", str(self.cpu)]
        subprocess.run(cmd, cwd=str(self.tree_out))
        species_tree_file                 = self.tree_out / "tree.treefile"
        normalized_species_tree_file_name = self.tree_out / "concatenate.nwk"
        shutil.copy(str(species_tree_file), str(normalized_species_tree_file_name))


    def make_outdirs(self):
        for sub in self.sub_dirs:
            sub.mkdir(parents=True, exist_ok=True)

    def get_all_species(self):
        all_species = set(self.seq2spe.values())
        return all_species

    def concatenate_iqtree(self):
        self.alignment_all_orthologs()
        self.concatenate()
        self.c_iqtree()

    def c_iqtree(self):
        seq = str(self.conc_seq)
        cpu = self.cpu
        cmd = "iqtree -s {0} -bb 1000 -nt {1}".format(seq, cpu)
        subprocess.call(cmd, shell=True)

    def reconstruct_concatenate_orthologs_tree_by_fasttree(self):
        self.alignment_all_orthologs()
        self.concatenate()
        self.fasttree()

    def count_otu(self, seq):
        cmd = 'grep -c ">" {0}'.format(seq)
        otus = int(subprocess.getoutput(cmd))
        return otus

    def alignment_all_orthologs(self):
        rate = self.share
        species_num = len(self.species)
        threshold = math.ceil(species_num * rate)
        if threshold < 4:
            threshold = 4
        instances = [MultipleAlignment(seq, self.align_out, self.trim_out) for seq in self.orthologues
                     if self.count_otu(seq) >= threshold]
        cpu = self.cpu
        with multiprocessing.Pool(cpu) as pool:
            for _ in pool.imap_unordered(multiple_alignment, instances):
                pass

    def concatenate(self):
        if self.conc_seq.exists():
            return
        trim_dir    = self.trim_out
        fastas      = sorted(trim_dir.glob("*.trim"))
        all_species = self.species
        species2seq = defaultdict(str)
        for fasta in fastas:
            spe2seq = {}
            seq_len = 0
            for record in SeqIO.parse(str(fasta), "fasta"):
                header  = record.id
                species = self.seq2spe[header]
                if species in spe2seq.keys():
                    continue
                seq              = str(record.seq)
                seq_len          = len(seq)
                spe2seq[species] = seq

            appear = set(spe2seq.keys())
            for sp, sq in spe2seq.items():
                species2seq[sp] += sq

            not_exist = all_species - appear
            gap       = "-" * seq_len
            for sp in not_exist:
                species2seq[sp] += gap

        outfile = str(self.conc_seq)
        with open(outfile, "w", encoding="UTF-8") as out:
            for species, seq in sorted(species2seq.items()):
                fasta = ">{0}\n{1}".format(species, seq)
                print(fasta, file=out)

    def fasttree(self):
        if self.conc_tree.exists():
            return
        seq     = str(self.conc_seq)
        outfile = str(self.conc_tree)
        log     = str(self.log)
        cmd     = "FastTreeMP -bionj -lg {0} > {1} 2>> {2}".format(seq, outfile, log)
        subprocess.call(cmd, shell=True)
