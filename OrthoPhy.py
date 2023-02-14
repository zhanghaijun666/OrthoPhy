#!/usr/bin/python3

"""\
OrthoPhy: A program to construct ortholog datasets using taxonomic information

USAGE: OrthoPhy.py --indir <dir> --out <dir> [options]

[OPTIONS]:
  --indir        <dir>    Input fasta or genbank format files directory
  --out          <dir>    Output directory
  --Inflation|-I <float>  the initial mcl inflation value [default:1.2]
  --delta        <float>  the increment value of the inflation [default:0.2]
  --end          <float>  the max inflation value that is to be accept (may not be last value) [default:2.0]
  --List|-l      <list>   User specific inflation range (not allow duplicate)
  --group|-g     <int>    the minimum taxonomic group count that is contained in each ortholog group [default:2]
  --otu|-O       <int>    the minimum otu number that is contained in each ortholog group [default:4]
  --species|-spe <int>    Minimum sharing rate of orthologs used for concatenated phylogenetic tree reconstruction [default: 0.5]
  --threads|-T   <int>    the max analysis threads [default:'max']
  --hgt                   execute hgt prediction if this option is specified [default:off]
  --help|-h               show this help message and exit
  --version|-V            show program's version number and exit
  --force|-f              Delete the exist output directory even if which is not empty
  --accurate              Use IQ-TREE in concatenate orthologs tree reconstruction
  --seperate              Estimate the concatenated tree by estimating the substitution matrix for each gene in IQ-TREE (i.e. use partition model).
  --no_filtering          Do not filter unreliable alignment region [default: off]
  --no_concatenate        do not infer concatenate tree [default:off]
"""

import argparse
import datetime
import multiprocessing
import os
import platform
import subprocess
import shutil
import sys
import textwrap
import traceback
from decimal import Decimal
from distutils.version import StrictVersion
from pathlib import Path
from glob import glob
from Bio import SeqIO
from bin import accept_input
from bin import reconstruct_tree
from bin import tree_cut
from bin import concatenate_tree
from bin import normalized_homology_search
from bin import normalized_mcl


class ArgumentParser(argparse.ArgumentParser):
    def format_help(self):
        return __doc__


def generate_opt_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--indir", required=True,
        help="Input fasta or genbank format files directory"
    )
    parser.add_argument(
        "--out", required=True,
        help="Output directory"
    )
    parser.add_argument(
        "--hgt", default=False, action="store_true",
        help="execute hgt prediction if this option is specified [default:off]"
    )
    parser.add_argument(
        "--no_concatenate", default=False, action="store_true",
        help="do not reconstruct concatenate tree [default:off]"
    )
    parser.add_argument(
        "--coverage", "-c", type=float, default=0.5,
        help="obsolete announced [default:0.5]"
    )
    parser.add_argument(
        "--Inflation", "-I", type=str, default="1.2",
        help="the initial inflation value [default:1.2]"
    )
    parser.add_argument(
        "--delta", type=str, default="0.2",
        help="the increment value of the inflation [default:0.2]"
    )
    parser.add_argument(
        "--end", type=str, default="2.0",
        help="the max inflation value that is to be accept (may not be last value) [defalut:2.0]"
    )
    parser.add_argument(
        "--group", "-g", type=int, default=2,
        help="the minimum taxonomic group count that is contained in each ortholog group [default:2]"
    )
    parser.add_argument(
        "--otu", "-O", type=int, default=4,
        help="the minimum otu count that is contained in each ortholog group [default:4]"
    )
    parser.add_argument(
        "--threads", "-T", type=int, default=multiprocessing.cpu_count(),
        help="the max"
             " analysis threads [default:{0}]".format(
            str(multiprocessing.cpu_count()))
    )
    parser.add_argument(
        "--species", "-sp", type=float, default=0.5,
        help="Minimum sharing rate of orthologous genes used for concatenated phylogenetic tree reconstruction"
    )
    parser.add_argument(
        "--force", "-f", action="store_true",
        help="Delete the exist output directory even if which is not empty"
    )
    parser.add_argument(
        "--retry", "-r", action="store_true",
        help="Reuse previous result and skip some process"
    )
    parser.add_argument(
        "--accurate", action="store_true",
        help="Use IQ-TREE in concatenate orthologs tree reconstruction"
    )
    parser.add_argument(
        "--seperate", action="store_true",
        help="Estimate the concatenated tree by estimating the substitution matrix for each gene in IQ-TREE (i.e. use partition model)"
    )
    parser.add_argument(
        "--nofiltering", action="store_true",
        help="Do not filtering unreliable alignment region [default:off]"
    )
    parser.add_argument(
        "--version", "-V",
        action="version", version="Ortholog-Finder2 version 1.0 (beta)",
        default="1.0"
    )
    parser.add_argument(
        "--List", "-l",
        default=None, nargs="*",
        help="User specific inflation range (not allow duplicate)"
    )
    args = parser.parse_args()
    return args


class OrthologFinder():
    def __init__(self, indir, outdir, cpu,
                 hgt, otu, group, inflation, delta, end,
                 coverage, no_concatenate, share, force, accurate,
                 seperate, version, no_filtering, inflation_list, retry):
        self.start          = datetime.datetime.today().replace(microsecond = 0)
        self.version        = version
        self.indir          = Path(indir)
        self.out            = Path(outdir)
        self.predict        = hgt
        self.otu            = otu
        self.group          = group
        self.share          = share
        self.cpu            = cpu
        self.inflation      = Decimal(inflation)
        self.delta          = Decimal(delta)
        self.end            = Decimal(end)
        self.no_concatenate = no_concatenate
        self.accurate       = accurate
        self.seperate       = seperate
        self.retry          = retry
        self.inflations     = []
        if inflation_list is None:
            inf = self.inflation
            while inf <= self.end:
                self.inflations.append(inf)
                inf += self.delta
        else:
            self.inflations = [Decimal(inf) for inf in inflation_list]
        self.coverage                = coverage * 100
        self.no_filtering            = no_filtering
        self.dirtype                 = self.validate_input()
        self.last                    = ""
        self.last_inflation          = ""
        self.force                   = force
        self.spe2seq                 = {}
        self.hgt_out                 = self.out / "00_hgt_prediction"
        self.fasta_out               = self.out / "01_fasta"
        self.homology_search_out     = self.out / "02_homology_search"
        self.clustering_out          = self.out / "03_clustering"
        self.gene_trees_out          = self.out / "04_gene_trees"
        self.ortholog_out            = self.out / "05_orthologs"
        self.species_tree_out        = self.out / "06_species_tree"
        self.log                     = self.out / "Ortholog-Finder2.log"
        self.ortholog_log            = self.ortholog_out.glob("**/log")
        self.sub_dirs                = [self.fasta_out, self.homology_search_out, self.clustering_out,
                                        self.gene_trees_out, self.ortholog_out]
        if self.predict:
            self.sub_dirs.append(self.hgt_out)
        if not self.no_concatenate:
            self.sub_dirs.append(self.species_tree_out)

    def validate_param(self):
        if not self.otu >= 2:
            raise ValueError("The minimum otu count must be 2 or larger")
        if not (0 < self.share <= 1.0):
            raise ValueError(
                "The minimum share rate must be between 0 and 1.0 or 1.0")
        if not self.group > 0:
            raise ValueError("The minimum group count must be natural number")
        if self.end < self.inflation:
            self.inflation = self.end
            self.end       = self.inflation
        if not self.inflation >= 1.0:
            raise ValueError(
                "The initial inflation value must be over 1.0 or 1.0")
        if not self.delta > 0:
            raise ValueError("The inflation increment must be natural number")
        if not (0 < self.coverage <= 100):
            raise ValueError(
                "The minimum coverage rate must be between 0 and 1.0 or 1.0")

    def inferring_ortholog(self):
        self.install_check()
        self.make_outdirs()
        self.analysis_log()

        if self.dirtype == "GH":
            self.logging("***** HGT prediction *****")
        elif self.dirtype == "GN":
            self.logging("***** Extract gene sequences *****")
        elif self.dirtype == "FA":
            self.logging("***** Accepted FASTA files *****")
        else:
            raise TypeError("Input file type is wrong!")
        self.make_input_fasta_dir()
        if self.dirtype == "GH":
            ortho_start     = datetime.datetime.today()
            ortho_start     = ortho_start.replace(microsecond = 0)
            hgt_time        = ortho_start - self.start
            ortho_start_str = ortho_start.isoformat()
            self.logging("{0} (Consume {1})".format(ortho_start_str, hgt_time))
        else:
            ortho_start = self.start
            ortho_start = ortho_start.replace(microsecond=0)

        self.logging("***** Homology search *****")
        self.exec_normalized_homology_search()

        for inf in self.inflations:
            self.logging("***** Inflation: {0} *****".format(str(inf)))
            str_inf      = str(inf)
            mcl_out      = self.clustering_out / "result_{0}".format(str_inf) / "groups"
            mafft_out    = self.gene_trees_out / str_inf / "mafft"
            trimal_out   = self.gene_trees_out / str_inf / "trimal"
            fasttree_out = self.gene_trees_out / str_inf / "fasttree"
            dup_out      = self.gene_trees_out / str_inf / "dup"
            ortholog_out = self.ortholog_out / str_inf

            self.logging("***** Clustering *****")
            self.normalized_clustering(inf)

            self.logging("***** Reconstruct gene trees *****")
            self.reconstruct_gene_tree(mcl_out, mafft_out, trimal_out,
                                       fasttree_out, dup_out)

            self.logging("***** Reconstruct orthologue groups *****")
            self.extract_ortholog(fasttree_out, ortholog_out, seq2spe=self.spe2seq)

            self.last_inflation = str(inf)

        ortho_end     = datetime.datetime.today()
        ortho_end     = ortho_end.replace(microsecond = 0)
        ortho_time    = ortho_end - ortho_start
        ortho_end_str = ortho_end.isoformat()
        self.logging(
            "!!!!! {0}: Ortholog inferrence was finished (Consume {1})!!!!!".format(
                ortho_end_str, ortho_time))
        if not self.no_concatenate:
            self.logging("***** Reconcstruct species tree *****")
            self.reconstruct_concatenate_tree()
            phylo_end     = datetime.datetime.today()
            phylo_end     = phylo_end.replace(microsecond = 0)
            phylo_time    = phylo_end - ortho_end
            phylo_end_str = phylo_end.isoformat()
            self.logging("Species tree inference has been finished {0} (Consume {1})".format(phylo_end_str, phylo_time))
        end = datetime.datetime.today()
        end = end.replace(microsecond=0)
        tol_time = end - self.start
        end_str  = end.isoformat()
        self.logging("End: {0} (Consume {1})".format(end_str, tol_time))

    def logging(self, text):
        print(text)
        with open(str(self.log), "a", encoding="UTF-8") as out:
            print(text, file=out)

    def validate_input(self):
        if not self.indir.exists():
            raise ValueError("Not found: {0}".format(str(self.indir)))
        direct_files = [direct_file for direct_file in self.indir.glob("*.*")
                        if direct_file.is_file()]
        if direct_files:
            raise ValueError(
                "Species seq file must be contained in any subdirectory")

        files = list(self.indir.glob("**/*.*"))
        if self.predict:
            args = [(gbk, "genbank", True) for gbk in files]
            with multiprocessing.Pool(self.cpu) as pool:
                for is_completed, gbk in pool.imap_unordered(OrthologFinder.validate_file, args):
                    print(is_completed)
                    if not is_completed:
                        raise ValueError(
                            "{0} is not a complete genome genbank file".format(
                                gbk))
            return "GH"

        args = [(infile, "genbank", False) for infile in files]
        guesses = []
        with multiprocessing.Pool(self.cpu) as pool:
            for file_type, gbk in pool.imap_unordered(OrthologFinder.validate_file, args):
                guesses.append(file_type)
        if all(guesses):
            return "GN"
        args = [(infile, "fasta", False) for infile in files]
        is_fastas = []
        with multiprocessing.Pool(self.cpu) as pool:
            for file_type, fasta in pool.imap_unordered(OrthologFinder.validate_file, args):
                is_fastas.append(file_type)
        if all(is_fastas):
            return "FA"
        else:
            raise ValueError("The all file type must be the same")

    @staticmethod
    def validate_file(args):
        infile, filetype, single = args
        infile = str(infile)
        if single:
            func = SeqIO.read
        else:
            func = SeqIO.parse
        try:
            record = func(infile, filetype)
        except:
            return False, infile
        else:
            return any(record), infile

    def install_check(self):
        version = StrictVersion("3.4")
        python  = StrictVersion(platform.python_version())
        if python < version:
            sys.exit("Required python 3.4 or later")

        dependencies = ["fasttree", "diamond", "mcl", "mafft", "trimal",
                        "seqkit"]
        need = []
        if self.predict:
            cmd = "Dimob.pl -h > /dev/null 2>&1"
            sc  = subprocess.call(cmd           , shell = True)
            if sc != 0:
                need.append("Dimob.pl (IslandPath-DIMOB)")
        if self.accurate or self.seperate:
            dependencies.append("iqtree2")
        for program in dependencies:
            is_executable = shutil.which(program)
            if not is_executable:
                need.append(program)
        if need:
            need = " ".join(need)
            msg  = "Following programs are not executable: {0}".format(need)
            sys.exit(msg)

    def make_outdirs(self):
        if str(self.indir) == str(self.out):
            raise ValueError("Output directory must not be input directory")
        if self.force and self.out.exists():
            shutil.rmtree(str(self.out))
        self.out.mkdir(parents=True, exist_ok=True)

        if (not self.is_empty()) and (not self.retry):
            raise ValueError("The output directory already exists and not empty.")
        if self.retry:
            for log in self.ortholog_log:
                shutil.rmtree(str(log))

        for sub in self.sub_dirs:
            sub.mkdir(exist_ok=True)

    def is_empty(self):
        children = list(self.out.glob("*"))
        if not children:
            return True

    def analysis_log(self):
        ver_info     = "Ortholog-Finder2 analysis log (version {0})".format(self.version)
        program_path = str(Path(sys.argv[0]).resolve())
        today        = self.start.isoformat()
        indir        = str(self.indir.resolve())
        outdir       = str(self.out.resolve())
        command_line = []
        command_line.append(program_path)
        for arg in sys.argv[1:]:
            if " " in arg:
                arg = "'{0}'".format(arg)
            command_line.append(arg)
        command_line = " ".join(command_line)
        inf_range    = [str(inf) for inf in self.inflations]
        inflation    = " ".join(inf_range)
        coverage     = self.coverage / 100
        threads      = self.cpu
        otu          = self.otu
        group        = self.group

        if self.predict:
            prediction = "ON"
        else:
            prediction = "OFF"

        if self.no_concatenate:
            concatenated = "OFF"
        else:
            concatenated = "ON (minimum species ortholog share rate among species: {0:.1%})".format(
                self.share)

        msg = textwrap.dedent('''\
                                {0}
                                {1}
                                Start             : {2}
                                Input directory   : {3}
                                Output directory  : {4}
                                HGT Prediction    : {5}
                                Inflation         : {6}
                                Coverage          : {7:.1%}
                                Minimum otu       : {8}
                                Minimum group     : {9}
                                Concatenated tree : {10}
                                Number of threads : {11}
                                '''.format(ver_info, command_line, today, indir, outdir,
                                           prediction, inflation,
                                           coverage,
                                           otu, group, concatenated, threads, ))
        self.logging(msg)

    def make_input_fasta_dir(self):
        acceptor = accept_input.Acceptor(self.indir,
                                         self.fasta_out, self.hgt_out,
                                         self.dirtype, self.cpu)
        retval = acceptor.generate_orthologfinder_input()
        self.spe2seq.update(retval)

    def exec_normalized_homology_search(self):
        fasdir         = self.fasta_out
        outdir         = self.homology_search_out
        normalized_hgs = normalized_homology_search.NormalizedHomologySearch(
            fasdir, outdir, coverage=self.coverage, cpu=self.cpu)
        normalized_hgs.normalized_homology_search()
        self.all_sequences_fasta = Path(outdir) / "db.faa"
        subprocess.call(
            'find {0} -name "*.faa" -print0 | xargs -0 cat > {1}'.format(
                fasdir, self.all_sequences_fasta), shell=True)

    def ancestor_clustering(self, inflation):
        anc_clustering_dir = self.ancestor_clustering_out / str(inflation)
        anc_clustering_dir.mkdir(parents=True, exist_ok=True)
        anc_out = anc_clustering_dir / "connected"
        mcl_groups = self.clustering_out / "result_{0}".format(inflation) / "groups"
        cmd = "clustering_ancestor_sequence.py {0} {1}".format(mcl_groups, anc_clustering_dir)
        subprocess.run(cmd, shell=True)
        return anc_out

    def normalized_clustering(self, inflation):
        fasdir     = self.fasta_out
        if self.last_inflation == "":
            results = self.homology_search_out / "normalized"
        else:
            results = self.clustering_out / "log" / self.last_inflation
        self.logging(str(results))
        clustering = normalized_mcl.NormalizedClustering(fasdir=fasdir,
                                                         indir=results,
                                                         last_result=self.last,
                                                         outdir=self.clustering_out,
                                                         otu=self.otu,
                                                         inflation=inflation)
        clustering.rbh_clustering()
        self.last = self.ortholog_out / str(inflation)

    def reconstruct_gene_tree(self, fas_dir, align_out, trim_out, tree_out, dup_out=None):
        indir        = Path(fas_dir)
        mafft_out    = Path(align_out)
        trimal_out   = Path(trim_out)
        fasttree_out = Path(tree_out)
        dup_out      = Path(dup_out)
        log          = str(self.log)
        cpu          = int(self.cpu)

        gene_tree_reconstructor = reconstruct_tree.GeneTreeReconstructor(indir,
                                                                         mafft_out,
                                                                         trimal_out,
                                                                         fasttree_out,
                                                                         cpu,
                                                                         log,
                                                                         self.spe2seq,
                                                                         self.no_filtering,
                                                                         dup_out=dup_out)
        gene_tree_reconstructor.multicore_reconstruction()

    def extract_ortholog(self, nwk_dir, outdir, seq2spe):
        nwk_dir = str(nwk_dir)
        taxon   = self.fasta_out
        outdir  = str(outdir)
        all_fas = self.all_sequences_fasta
        otu     = self.otu
        group   = self.group
        cpu     = self.cpu
        log     = self.log
        ortholog_groups = tree_cut.TreeCutter(nwkdir=nwk_dir, taxondir=taxon,
                                              cutTree_output=outdir,
                                              allfasta=all_fas, otu_limit=otu,
                                              group_limit=group,
                                              cpu=cpu, log=log, s2s=seq2spe)
        ortholog_groups.generate_ortholog_groups()
        orthologues = os.path.join(outdir   , "*.ortholog")
        count       = len(glob(orthologues))
        self.logging("Ortholog groups: {0}".format(str(count)))

    def reconstruct_concatenate_tree(self):
        reconstructor = concatenate_tree.ConcatenateTree(self.ortholog_out,
                                                         self.species_tree_out,
                                                         self.spe2seq,
                                                         self.otu, self.log,
                                                         self.share, self.cpu)
        if self.accurate:
            reconstructor.concatenate_iqtree()
        elif self.seperate:
            reconstructor.iqtree2()
        else:
            #reconstructor.reconstruct_concatenate_orthologs_tree_by_fasttree()
            reconstructor.astral()


if __name__ == '__main__':
    st_t           = datetime.datetime.today()
    args           = generate_opt_parser()
    indir          = args.indir
    outdir         = args.out
    prediction     = args.hgt
    coverage       = args.coverage
    inflation      = args.Inflation
    delta          = args.delta
    end            = args.end
    otu            = args.otu
    group          = args.group
    cpu            = args.threads
    concatenate    = args.no_concatenate
    share          = args.species
    force          = args.force
    retry          = args.retry
    version        = args.version
    accurate       = args.accurate
    seperate       = args.seperate
    no_filtering   = args.nofiltering
    inflation_list = args.List
    if args.threads is None:
        cpu = multiprocessing.cpu_count()
    try:
        ortholog_finder = OrthologFinder(indir=indir, outdir=outdir, cpu=cpu,
                                         hgt=prediction, otu=otu, group=group,
                                         inflation=inflation, delta=delta,
                                         end=end, coverage=coverage,
                                         no_concatenate=concatenate,
                                         share=share, force=force,
                                         accurate=accurate, seperate=seperate,
                                         version=version,
                                         no_filtering=no_filtering, inflation_list=inflation_list,
                                         retry=retry)
        ortholog_finder.inferring_ortholog()
    except (Exception, KeyboardInterrupt) as e:
        err = traceback.format_exc()
        log = list(Path(outdir).glob("Ortholog-Finder2.log"))
        if not log:
            print(err, end="")
            sys.exit(1)
        log   = Path(log[0])
        mtime = log.stat().st_mtime
        mdate = datetime.datetime.fromtimestamp(mtime)
        if mdate >= st_t:
            end_time = datetime.datetime.today()
            end      = end_time.isoformat()
            err_log  = str(log)
            print(err, end="")
            with open(err_log, "a", encoding="UTF-8") as o:
                o.write(err)
                print(
                    "The analysis was aborted due to unexpected error ({0})".format(
                        end), file=o)
        else:
            print(err, end="")
        sys.exit(1)
