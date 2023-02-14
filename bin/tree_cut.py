#!/usr/bin/env python3

from itertools import product
from collections import defaultdict
import re
import sys
import pickle
import itertools
import subprocess

from Bio import SeqIO
from pathlib import Path
from decimal import Decimal, ROUND_DOWN
from multiprocessing import Pool, cpu_count

from bin.convert_newick2node import Newick2Node
from bin.distanceotu import DistanceOTU

newick2node = Newick2Node(branch_seq="&").newick2node


def logging(msg):
    with open(str(LOG_FILE), "a", encoding="UTF-8") as log:
        print(msg)
        print(msg, file=log)


# create a dictionary of grouping taxon every hierarchy
def generate_index(tax_file):
    tax    = Path(tax_file)
    gindex = {record.id: (tax.stem + "-" + str(i).zfill(6)) for i, record in enumerate(SeqIO.parse(str(tax), 'fasta'))}
    return gindex


# read number of taxon at selected hierarchy
def read_hierarchy_info(hierarchy):
    global hierarchies
    # create a dictionary of Pseudo-fractionation of classification information and reverse dictionary
    decimal_tax = {}  # Pseudo-fractionation of classification information
    only_species_num = max([i[-1] for i in hierarchies.values()]) + 1  # number of unclassified single line species

    for sp_num, (sp_file, tax_num_list) in enumerate(hierarchies.items()):
        sp = sp_file.stem

        if hierarchy < len(tax_num_list):
            tax_num = tax_num_list[hierarchy]
        else:
            # decimal_tax[sp] = "{0}.{1}".format(only_species_num, str(sp_num).rjust(4, '0'))
            # only_species_num += 1
            continue

        # No taxon data spicies labeled minus
        if "__Ignore__" in str(sp_file):
            decimal_tax[sp] = "-" + str(tax_num) + "." + str(sp_num).rjust(4, '0')
        elif "__Unclassification__" in str(sp_file):
            decimal_tax[sp] = "{0}.{1}".format(only_species_num, str(sp_num).rjust(4, '0'))
            only_species_num += 1
        else:
            decimal_tax[sp] = str(tax_num) + "." + str(sp_num).rjust(4, '0')

    r_decimal_tax = {str(v): k for k, v in decimal_tax.items()}  # Reverse dictionary
    return decimal_tax, r_decimal_tax


# judge lowest threshold which of OTU or group
def judge_over(nl, tree_file):
    node_list = nl[:]
    mono_list = list_to_flat([re.split('[,&]', x) for x in node_list if x])
    group_list = set(map(lambda x: int(Decimal(x)), mono_list))  # set object with int(group)
    species_list = set(map(lambda x: str(Decimal(x).quantize(Decimal('0.0001'), rounding=ROUND_DOWN)),
                           mono_list))  # set object with 4th decimal place(tax)
    if len(group_list) < limit_group:
        msg = "-*" * 40 + \
              "\nThis newick:{0} doesn't meet the condition of the number of groups({1}).\n".format(tree_file.stem,
                                                                                                    limit_group) + \
              "-*" * 40
        logging(msg)
        write_log_file(tree_file, "condition")
        return True

    elif len(species_list) < limit_otu:
        msg = "-*" * 40 + \
              "\nThis newick:{0} doesn't meet the condition of the number of OTU(1).\n".format(tree_file.stem,
                                                                                               limit_otu) + \
              "-*" * 40
        logging(msg)
        write_log_file(tree_file, "condition")
        return True
    else:
        return False


def newick_2_decimal(tree_file, decimal_tax):
    with open(str(tree_file), "r") as TF:
        newick = TF.read().replace("\n", "").replace(";", "")

    # OTU chenges to decimal
    d_new = newick
    for otu in re.sub(r"\):\d+\.\d+|\(|\)", "", d_new).split(","):
        bl = re.search(r":\d+\.\d+", otu).group()
        TODO = d_new
        d_new = d_new.replace(otu, numbering_seq[re.sub(":.+", "", otu)] + bl)

    while re.search(r"([,(])(\w+?-\d+)", d_new):
        m = re.search(r"([,\(])(\w+?-\d+)", d_new)
        replaced = m.group()
        replaced = replaced.replace("(", "\(")
        point = m.group(2)
        symbol = m.group(1)
        tax, seq_num = point.split("-")
        if tax in decimal_tax:
            decimal_otu = decimal_tax[tax] + seq_num.zfill(6)
        else:
            logging("no data of spicies or tax")
            decimal_otu = "-9999.0000" + seq_num.zfill(6)
        chenged = r"{0}{1}".format(symbol, decimal_otu)
        d_new = re.sub(replaced, chenged, d_new)
        # d_new = d_new.replace(point, decimal_otu)

    length_new = d_new + ";"  # decimal newick with branch length
    while re.search(r"\):\d+\.\d+|\)\d+\.\d+", d_new):
        d_new = re.sub(r"\):\d+\.\d+|\)\d+\.\d+", ")", d_new)
    d_new = re.sub(r":\d+\.\d+", "", d_new)

    # newick format chenges to node format
    node_list = Newick2Node(branch_seq="&").newick2node(d_new)
    return node_list, length_new


# write log files which saving ortholog candidate
def write_log_file(tree_file, log_type, otu_list=None):
    log_dic = {"condition": "not_meet_conditions.log",
               "otu": "under_{0}OTU.log".format(limit_otu),
               "group": "under_{0}group.log".format(limit_otu),
               "cut": "none_cut.log",
               "ortho": "ortho.log"}
    file_name = log_dic[log_type]
    if otu_list == None:
        with open(str(tree_file), "r") as TF:
            newick = TF.read().replace("\n", "").replace(";", "")
        otu_list = [x for x in re.sub(r"\(|\)|:\d\.\d+|;", "", newick).split(",")]
        try:
            otu_species = {seq2spe[ort] for ort in otu_list}
        except KeyError:
            print(tree_file.stem)
            print(seq2spe)
        # 除去した重複配列の復元
        dup_orthologs = []
        for otu in otu_list:
            dups = seq2dups.get(otu, [])
            for d in dups:
                sp = seq2spe[d]
                if sp not in otu_species:
                    dup_orthologs.append(d)
        otu_list += dup_orthologs

    with open(str(log_dir / file_name), "a") as o:
        o.write("{0} :".format(tree_file.stem))
        for x in otu_list:
            o.write("\t" + x)
        o.write("\n")


def list_to_flat(list_list, del_duplication=0):
    """
    list which has list -> list which has no list

    parameters
    ----------
    list_list : list which has list

    returns
    -------
    flat_list : list which has no list
        :param del_duplication:
        remove duplication (use:1, not use:0, [default:0])
    """
    flat_list = []  # list wiich has no list
    for flat in list_list:

        if type(flat) is set:
            flat = list(flat)
        if type(flat) is list:
            flat_list.extend(list_to_flat(flat))
        else:
            flat_list.append(flat)
    if del_duplication:
        flat_list = list(set(flat_list))
    return flat_list


# remove paralogs from tree by node format
# def cut_clustar(node_list, paralogs):
def cut_clustar(nl, paralogs):
    node_list = nl[:]
    # select cutting paralog
    for paralog in paralogs:
        for nn, node in enumerate(node_list):
            if node is None:
                continue

            node = [x.split("&") for x in node.split(",")]
            # if cluster in node is same paralog and tree has over 2 OTU, the cluster is removed.
            if ([paralog] in node) and len(node_list) != 1:
                node_list[nn] = None
                continue

            # select cluster in every node
            for bn, cluster in enumerate(node):
                # if cluster has paralog and tree has over 2 OTU, the paralog removed cluster.
                if paralog in cluster:
                    cluster.remove(paralog)
                    cluster.sort()
                    if len(cluster) > 0:
                        node = ["&".join(x) for x in node]
                        node_list[nn] = ",".join(node)
                    else:
                        if len(node_list) != 1:
                            node_list[nn] = None
                        else:
                            node = ["&".join(x) for x in node if x]
                            node_list[nn] = ",".join(node)
                    break
        node_list = [x for x in node_list if x]
    node_list = [x for x in node_list if x]

    # return [x for x in node_list if x]
    # if None in node_list:
    #     cmd = "echo '{0}'".format(node_list)
    #     subprocess.run(cmd)
    return node_list


# find inparalog tax
def find_inparalog(nl, distance=None):
    node_list = nl[:]
    inparalog = []
    for node in node_list:
        #cluster_list = [x.split("&") for x in node.split(",") if
                        #"-" not in x]  # list with 3 cluster at a node removing not ignore species
        cluster_list = [x.split("&") for x in node.split(",")]  # list with 3 cluster at a node removing not ignore species
        for cluster in cluster_list:
            remaining_cluster   = list_to_flat([x for x in cluster_list if x != cluster])  # aren't selected clusters
            sp_branch           = list(map(lambda x: str(Decimal(x).quantize(Decimal('0.0001'), rounding=ROUND_DOWN)),
                                    cluster))  # list with 4th decimal place(tax)
            remaining_sp_branch = list(map(lambda x: str(Decimal(x).quantize(Decimal('0.0001'), rounding=ROUND_DOWN)),
                                      remaining_cluster))  # list with 4th decimal place(tax)

            if len(set(sp_branch)) == 1 and len(sp_branch) != 1:
                other_leaf          = remaining_cluster[0]
                min_len             = float("inf")
                min_inparalog_index = ""
                for index, candidate_inparalog in enumerate(cluster):
                    branch_length = distance.otu_distance(other_leaf, candidate_inparalog)
                    if branch_length < min_len:
                        min_inparalog_index = index
                inps = cluster[:]
                inps.pop(min_inparalog_index)
                inparalog += inps
            if len(set(remaining_sp_branch)) == 1 and len(remaining_sp_branch) != 1:
                other_leaf          = cluster[0]
                min_len             = float("inf")
                min_inparalog_index = ""
                for index, candidate_inparalog in enumerate(remaining_cluster):
                    branch_length = distance.otu_distance(other_leaf, candidate_inparalog)
                    if branch_length < min_len:
                        min_inparalog_index = index
                inps = remaining_cluster[:]
                inps.pop(min_inparalog_index)
                inparalog += inps
    return set(inparalog)


# choose a cluster to remove?by the number of otu or the shortness of branches
def get_farther_leaf(cluster0, cluster1, distance):
    if type(distance.node_list[0]) is list:
        in_tree = list_to_flat(distance.node_list[0])
    elif type(distance.node_list[0]) is str:
        in_tree = re.split("[,&]", distance.node_list[0])
    else:
        print(distance.node_list)

    not_in = [otu for otu in (cluster0 + cluster1) if otu not in in_tree]

    if len(not_in) > 0:
        cluster0 = [otu for otu in cluster0 if otu not in not_in]
        cluster1 = [otu for otu in cluster1 if otu not in not_in]
        if len(cluster0) * len(cluster1) == 0:
            return None

    cluster_distance = distance.cluster_distance(cluster0, cluster1, midway_expression=1)
    if cluster_distance[0] >= cluster_distance[1]:
        return cluster0
    else:
        return cluster1


# find visible out paralog
def find_out_paralog(cluster_list, distance=None, removed_all=False):
    out_paralogs = set()
    leaves       = [len(cluster) > 1 for cluster in cluster_list]
    leaves       = sum(leaves)
    if leaves < 2:
        return False, None
    for x, y in itertools.combinations(cluster_list, 2):
        # list with 4th decimal place(tax)
        x_sp = list(map(lambda xx: str(Decimal(xx).quantize(Decimal('0.0001'), rounding=ROUND_DOWN)), x))
        y_sp = list(map(lambda xx: str(Decimal(xx).quantize(Decimal('0.0001'), rounding=ROUND_DOWN)), y))
        # set object with integer group
        x_group, y_group = set(map(lambda xx: int(Decimal(xx)), x)), set(
            map(lambda yy: int(Decimal(yy)), y))

        # If 2 cluster has a same taxon which has overlap tax,
        # the cluster judged causing gene duplication (visible out paralog)
        # If group has taxon which has minus number (is meanning No taxon data), don't judge
        if (x_group == y_group and len(x_group)) == 1 and (not removed_all):
            # If group has taxon which has minus number (is meanning No taxon data), don't judge
            if "-" in x_sp[0]:
                return False, None
            x_sp, y_sp = set(x_sp), set(y_sp)
            # judged cluster has overlapping tax
            if len(x_sp & y_sp) >= 1:
                try:
                    visible_out_paralog = get_farther_leaf(x, y, distance)
                except Exception:
                    if len(x_sp) >= len(y_sp):
                        visible_out_paralog = y
                    else:
                        visible_out_paralog = x
                # choose a cluster to remove
                if not visible_out_paralog:
                    continue
                out_paralogs = out_paralogs | set(visible_out_paralog)
            # judged cluster hasn't overlapping tax
            elif len(x_sp & y_sp) == 0:
                continue
        elif removed_all:
            # If group has taxon which has minus number (is meanning No taxon data), don't judge
            if "-" in x_sp[0]:
                return False, None
            x_sp, y_sp = set(x_sp), set(y_sp)
            if len(x_sp & y_sp) >= 1:  # judged cluster has overlapping tax
                # choose a cluster to remove
                try:
                    visible_out_paralog = get_farther_leaf(x, y, distance)
                except Exception:
                    if len(x_sp) >= len(y_sp):
                        visible_out_paralog = y
                    else:
                        visible_out_paralog = x
                if not visible_out_paralog:
                    continue
                out_paralogs = out_paralogs | set(visible_out_paralog)
            elif len(x_sp & y_sp) == 0:  # judged cluster hasn't overlapping tax
                continue
    if out_paralogs:
        return True, out_paralogs
    else:
        return False, None


# judge both clusters contain all species, when a phylogenetic tree is cut at an internal branch

# if clusters contain all species, cut phylogenetic tree at the internal branch
def find_invisible_outparalog(nl, distance, root=False):
    node_list = nl[:]
    # list with invisible out paralogs
    invisible_out_paralog = list()
    for node in node_list:
        tax_list   = [x.split("&") for x in node.split(",")]  # list with clusters at a node
        group_list = [set(map(lambda xx: int(Decimal(xx)), x)) for x in tax_list]  # list has set object with int(group)
        # choose a branch
        for i, cluster_x in enumerate(group_list):
            remaining_taxon    = [otu for n, otu in enumerate(tax_list) if n != i]
            remaining_clusters = [cluster for n, cluster in enumerate(group_list) if n != i]
            cluster_y          = set(list_to_flat(remaining_clusters))  # not chosen clusters

            if len(remaining_clusters) == 1:
                continue

            # If group has taxon which has minus number (is meanning No taxon data), don't judge
            cluster_x = set([x for x in cluster_x if x > 0])
            cluster_y = set([y for y in cluster_y if y > 0])

            # if either cluster has only no data taxon or both cluster has same a tax, continue
            if len(cluster_x) * len(cluster_y) in (0, 1):
                continue

            # if taxon in selected cluster are same the taxon in one of the other cluster,
            # the cluster is judged invisible out paralogs
            if cluster_x == cluster_y:
                if (cluster_x != max(remaining_clusters, key=len)) or (all(cluster_x == cl for cl in remaining_clusters)):
                    try:
                        fc = get_farther_leaf(tax_list[i], list_to_flat(remaining_taxon), distance)
                        invisible_out_paralog.append(fc)
                    except Exception:
                        invisible_out_paralog.append(tax_list[i])

            # if all the taxon in the selected cluster are contained in the other cluster,
            # the cluster may be out paralogs
            elif ((cluster_y & cluster_x) == cluster_x) and not root:
                # if selected cluster contains exactly the same taxon as one of the other,
                # the cluster is judged invisible out paralogs
                if cluster_x in remaining_clusters and len(cluster_x) > 1:
                    # if taxon in one of the other cluster is contained the other,
                    # the cluster is judged invisible out paralogs
                    if remaining_clusters[0] & remaining_clusters[1] == min(remaining_clusters, key=len):
                        distance.node_list = node_list
                        try:
                            farther_cluster = get_farther_leaf(min(remaining_taxon, key=len), tax_list[i], distance)
                        except Exception:
                            farther_cluster = tax_list[i]
                        invisible_out_paralog.append(farther_cluster)
                    else:
                        invisible_out_paralog.append(tax_list[i])

            # if all the taxon in one of the other cluster are contained in the other,
            # the cluster is judged invisible out paralogs when user wish it
            elif (len(cluster_x & cluster_y) == 0 and detect_invisible == 1) and not root:
                if (remaining_clusters[0] & remaining_clusters[1]) in remaining_clusters and len(
                        max(remaining_clusters, key=len)) > 1:
                    try:
                        fc = get_farther_leaf(remaining_taxon[0], remaining_taxon[-1], distance)
                    except Exception:
                        fc = min(remaining_taxon, key=len)
                    invisible_out_paralog.append(fc)
        if invisible_out_paralog:
            return set(list_to_flat(invisible_out_paralog))


    invisible_out_paralog = list_to_flat(invisible_out_paralog)
    all_invisible_out_paralogs = set(invisible_out_paralog)
    return all_invisible_out_paralogs


# detect out paralogs when tree has root
def find_detect_outparalog(node_list, distance):
    invisible_out_paralog = list()  # list with invisible out paralogs
    for node in node_list:
        # seeing other taxon as roots
        node       = re.sub(r"-\d+\.\d+", "-888", node)
        node       = re.sub(r"(-888&)+", "-888&", node)
        tax_list   = [x.split("&") for x in node.split(",")]  # list with clusters at a node
        group_list = [set(map(lambda xx: int(Decimal(xx)), x)) for x in tax_list]  # list has set object with int(group)
        if len(group_list) > 3:
            continue
        # choose a branch
        for i, root_cluster in enumerate(group_list):
            # check if "cluster_x" is root
            if -888 not in root_cluster:
                continue

            # list with otu (decimal)
            taxon_x, taxon_y     = [otu for n, otu in enumerate(tax_list) if n != i]
            # set object with group number (int)
            cluster_x, cluster_y = [cluster for n, cluster in enumerate(group_list) if n != i]
            cluster_x, cluster_y = set([x for x in cluster_x if x > 0]), set([y for y in cluster_y if y > 0])

            # if either cluster has only no data taxon or same a tax, continue
            if len(cluster_x) * len(cluster_y) in (0, 1):
                continue

            elif len(cluster_x & cluster_y) > 0:

                # select cluster which contain more species
                if len(taxon_x) > len(taxon_y):
                    invisible_out_paralog.append(taxon_y)
                elif len(taxon_y) > len(taxon_x):
                    invisible_out_paralog.append(taxon_x)

                # both cluster contain the same number
                else:
                    # select cluster which contain more taxon
                    if len(cluster_x) > len(cluster_y):
                        invisible_out_paralog.append(taxon_y)
                    elif len(cluster_y) > len(cluster_x):
                        invisible_out_paralog.append(taxon_x)
                    else:
                        distance.node_list = node_list
                        distance.key_stone()
                        farther_cluster = get_farther_leaf(taxon_x, taxon_y, distance)
                        invisible_out_paralog.append(farther_cluster)
    return list_to_flat(invisible_out_paralog, del_duplication=1)


# OTU revert from decimal or OTU convert to decimal in node list
def revert_conevert_otu(node_list, convert_dict, mode=0):
    convert_node_list = list()
    for node in node_list:
        # OTU revert from decimal
        if type(node) == str:
            convert_cluster = list()
            for cluster in node.split(","):
                convert_otu = list()
                for otu in cluster.split("&"):
                    sp = str(Decimal(otu).quantize(Decimal('0.0001'), rounding=ROUND_DOWN))
                    seq = otu.replace(sp, "")
                    if sp in convert_dict:
                        convert_otu.append(convert_dict[sp] + "-" + seq)
                    elif mode:
                        convert_otu.append(otu)
                convert_cluster.append(convert_otu)
            convert_node_list.append(convert_cluster)

        # OTU convert to decimal
        else:
            convert_cluster = list()
            for cluster in node:
                convert_otu = list()
                for otu in cluster:
                    if re.search(r".+-\d", otu):
                        sp, seq = otu.split("-")
                        convert_otu.append(convert_dict[sp] + seq)
                    elif mode:
                        convert_otu.append(otu)
                convert_cluster.append("&".join(convert_otu))
            convert_node_list.append(",".join(convert_cluster))

    return convert_node_list


def find_ortholog(tree_file):
    # numbering the top layer classification
    decimal_tax, r_decimal_tax = read_hierarchy_info(0)
    # newick format change to node format and OTU change to decimal
    try:
        node_list, length_new = newick_2_decimal(str(tree_file), decimal_tax)
    except Exception:
        return
    try:
        distance = DistanceOTU(length_new)
    except IndexError:
        write_log_file(tree_file, "cut")
        return
    if judge_over(node_list, tree_file):
        return None


    # organize into small groups
    # cut inparalog
    distance.node_list = node_list
    distance.key_stone()
    if inparalog_cut:
        inparalog = find_inparalog(node_list, distance)  # inparalog taxon list
        if len(inparalog) > 0:
            # cut inparalog from tree by node format
            node_list = cut_clustar(node_list, inparalog)
            if judge_over(node_list, tree_file):
                return None

    # cut visible out paralog(taxon are duplicated)
    if invisible_outparalog_cut:
        distance.node_list   = node_list
        distance.key_stone()
        invisible_outparalog = find_invisible_outparalog(node_list, distance, root=True)
        node_list            = cut_clustar(node_list, invisible_outparalog)
        while invisible_outparalog:
            distance.node_list = node_list
            distance.key_stone()
            invisible_outparalog = find_invisible_outparalog(node_list, distance, root=True)
            node_list            = cut_clustar(node_list, invisible_outparalog)
        if judge_over(node_list, tree_file):
            return None

    # cut visible out paralog(species are duplicated)
    distance.node_list = node_list
    distance.key_stone()
    CUTTED_NODE_LIST   = node_list[:]
    OUTPARALOG_FLAG    = True
    while OUTPARALOG_FLAG:
        OUTPARALOG_FLAG = False
        if visible_outparalog_cut:
            distance.node_list = CUTTED_NODE_LIST
            distance.key_stone()
            for cluster in CUTTED_NODE_LIST:
                if cluster is None:
                    continue
                # list with 3 cluster at a node
                cluster_list = [x.split("&") for x in cluster.split(",")]
                try:
                    judge_cut, visible_out_paralog = find_out_paralog(cluster_list, distance)
                except UnboundLocalError:
                    OUTPARALOG_FLAG = False
                    break
                if visible_out_paralog:
                    CUTTED_NODE_LIST = cut_clustar(CUTTED_NODE_LIST, visible_out_paralog)
                    OUTPARALOG_FLAG = True
                    if judge_over(CUTTED_NODE_LIST, tree_file):
                        return None
                    break

    node_list = CUTTED_NODE_LIST[:]
    # tree by node format revert to numbering the top layer classification
    node_list                  = revert_conevert_otu(node_list, r_decimal_tax)
    decimal_tax, r_decimal_tax = read_hierarchy_info(0)
    node_list                  = revert_conevert_otu(node_list, decimal_tax)

    # if tree contain unclassified spicies more than 2, cut tree when user want to do.
    if visible_outparalog_cut and len(node_list) > 0:
        # list with custer presumed to be gene duplication
        duple_cluster  = list()
        duaple_unknown = list(map(lambda x: str(Decimal(x).quantize(Decimal('0.0001'), rounding=ROUND_DOWN)),
                                  re.split("[,&]", node_list[0])))
        # list with unclassified OTU causing gene duplication
        duaple_unknown = [x for x in duaple_unknown if duaple_unknown.count(x) > 1]

        # search cluster which contain unclassified OTU causing gene duplication
        for unknown in set(duaple_unknown):
            for node in node_list:
                # cluster which contain unclassified OTU causing gene duplication
                unknown_cluster = list_to_flat([x.split("&") for x in node.split(",") if unknown in x])
                # whether there is a dropout of duplicate genes
                if len([x for x in unknown_cluster if unknown in x]) == len(duaple_unknown):
                    duple_cluster.append(unknown_cluster)
            if len(duple_cluster) > 0:
                # sallest of cluster which contain unclassified OTU causing gene duplication
                duple_cluster = min(duple_cluster, key=len)
                node_list = cut_clustar(node_list, duple_cluster)

    if invisible_outparalog_cut:
        distance.node_list = node_list
        distance.key_stone()
        invisible_outparalog = find_invisible_outparalog(node_list, distance, root=False)
        node_list            = cut_clustar(node_list, invisible_outparalog)
        while invisible_outparalog:
            distance.node_list   = node_list
            distance.key_stone()
            invisible_outparalog = find_invisible_outparalog(node_list, distance, root=False)
            node_list            = cut_clustar(node_list, invisible_outparalog)
        if judge_over(node_list, tree_file):
            return None

        ori_decimal_tax   = decimal_tax
        ori_r_decimal_tax = r_decimal_tax
        for hierarchy in range(max_hierarchy):
            if hierarchy == 0:
                continue

            next_decimal = read_hierarchy_info(hierarchy)
            decimal_tax.update(next_decimal[0])  # decimal dict update to next hierarchy tax number
            # decimal tree by node format update to next hierarchy tax number
            next_node_list = revert_conevert_otu(node_list, r_decimal_tax)

            next_node_list = revert_conevert_otu(next_node_list, decimal_tax, mode=1)

            # set with low-level taxon number to verify
            target_group_numbers = set([int(Decimal(x)) for x in re.split("[,&]", next_node_list[0])])

            # set with taxon number in tree
            tree_contein_numbers = set(map(lambda x: int(Decimal(x)), next_decimal[0].values()))

            if len(target_group_numbers & tree_contein_numbers) == 0:
                continue

            length_new = newick_2_decimal(tree_file, decimal_tax)[1]
            distance   = DistanceOTU(length_new)

            for positive_num in target_group_numbers:
                # not targetting taxon number chenge to minus
                minus_numbers = [otu for otu in re.split("[,&]", node_list[0]) if
                                 (int(Decimal(otu)) != positive_num) and "-" not in otu]

                root_node_list = [re.sub("({0})".format("|".join(minus_numbers)), r"-\1", cluster) for cluster in node_list]
                # OTU revert from decimal in target taxon
                root_node_list = revert_conevert_otu(root_node_list, r_decimal_tax, mode=1)
                # OTU chenge to next hierarchy decimal in target taxon
                root_node_list = revert_conevert_otu(root_node_list, decimal_tax, mode=1)

                # if tree has minus number
                distance.node_list = node_list
                distance.key_stone()
                if len(minus_numbers) > 0:
                    try:
                        invisible_outparalog = find_detect_outparalog(root_node_list, distance)
                    except ValueError:
                        continue
                else:
                    invisible_outparalog = find_invisible_outparalog(root_node_list, distance, root=False)
                    next_node_list       = cut_clustar(root_node_list, invisible_outparalog)
                    while invisible_outparalog:
                        try:
                            invisible_outparalog = find_invisible_outparalog(root_node_list, distance, root=False)
                        except Exception:
                            break
                        node_list = cut_clustar(root_node_list, invisible_outparalog)
                    if judge_over(node_list, tree_file):
                        return None

                next_node_list = cut_clustar(next_node_list, invisible_outparalog)

            node_list = next_node_list
            r_decimal_tax.update(next_decimal[1])

    decimal_tax, r_decimal_tax = ori_decimal_tax, ori_r_decimal_tax
    CUTTED_NODE_LIST   = node_list[:]
    distance.node_list = CUTTED_NODE_LIST
    distance.key_stone()
    OUTPARALOG_FLAG = True
    while OUTPARALOG_FLAG:
        OUTPARALOG_FLAG = False
        if visible_outparalog_cut:
            distance.node_list = CUTTED_NODE_LIST
            distance.key_stone()
            for cluster in CUTTED_NODE_LIST:
                if cluster is None:
                    OUTPARALOG_FLAG = False
                    continue
                cluster_list = [x.split("&") for x in cluster.split(",")]  # list with 3 cluster at a node
                try:
                    judge_cut, visible_out_paralog = find_out_paralog(cluster_list, distance, removed_all=True)
                except UnboundLocalError:
                    OUTPARALOG_FLAG = False
                    if judge_over(CUTTED_NODE_LIST, tree_file):
                        return None
                    break
                if visible_out_paralog:
                    CUTTED_NODE_LIST = cut_clustar(CUTTED_NODE_LIST, visible_out_paralog)
                    OUTPARALOG_FLAG  = True
                    if judge_over(CUTTED_NODE_LIST, tree_file):
                        return None
                    break

    node_list = CUTTED_NODE_LIST[:]

    # check monophyly
    TAXA = set()
    MONOPHYLY_TAXA = set()
    IS_MONOPHYLY = False
    for node in node_list:
        tmp = [list(set(list(map(lambda x: int(Decimal(x)), (x.split("&")))))) for x in node.split(",")]
        cl = []
        for n in tmp:
            if n in cl and len(n) == 1:
                continue
            else:
                cl.append(n)
        clusters = list_to_flat(cl)
        # If group has taxon which has minus number (is meanning No taxon data), don't judge
        clusters = [x for x in clusters if x > 0]

        TAXA = TAXA | set(clusters)
        # if list's length dosen't same set object' length, the tree isn't monophyly
        # if len(clusters) != len(set(clusters)) and len(set(clusters)) > 2:
        #     write_log_file(tree_file, "cut")
        #     break
        if len(clusters) != len(set(clusters)):
            # write_log_file(tree_file, "cut")
            continue
        else:
            for c in cl:
                if len(c) == 1:
                    MONOPHYLY_TAXA.add(c[0])

    if TAXA == MONOPHYLY_TAXA:
        #decimal_tax, r_decimal_tax = read_hierarchy_info(0)
        decimal_tax, r_decimal_tax = ori_decimal_tax, ori_r_decimal_tax
        binary_orthologs = re.split("[,&]", node)  # list with?digitized otu
        group_list = set(list(map(lambda otu: int(Decimal(otu)), binary_orthologs)))  # list with taxon

        orthologs = list_to_flat(revert_conevert_otu(node_list, r_decimal_tax), del_duplication=1)  # list with otu
        orthologs = [seq_id[x] for x in orthologs]
        orthologs_species = {seq2spe[ort] for ort in orthologs}
        #
        dup_orthologs = []
        for ortholog in orthologs:
            dups = seq2dups.get(ortholog, [])
            for d in dups:
                sp = seq2spe[d]
                if sp not in orthologs_species:
                    orthologs_species.add(sp)
                    dup_orthologs.append(d)

        orthologs += dup_orthologs

        # check meeting conditions
        if limit_group > len(group_list):
            write_log_file(tree_file, "group", orthologs)
        elif limit_otu > len(orthologs):
            write_log_file(tree_file, "otu", orthologs)
        else:
            write_log_file(tree_file, "ortho", orthologs)

            # create multi fasta file of ortholog
            with open(str(out_dir / (tree_file.stem + ".ortholog")), "w") as o:
                for seq in orthologs:
                    seq = "^{0}$".format(str(seq))
                    cmd = ["seqkit", "grep", "-nrp", seq, str(all_fasta)]
                    o.write(subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0].decode('utf-8'))
    else:
        write_log_file(tree_file, "cut")


# *******************************************************************************

class TreeCutter():
    def __init__(self, nwkdir, taxondir, cutTree_output, allfasta,
                 otu_limit, group_limit=2, cpu=1, log="", s2s=None):
        self.nwk                = nwkdir
        self.taxon              = taxondir
        self.out                = cutTree_output
        self.all_sequence_fasta = allfasta
        self.otu                = int(otu_limit)
        self.group              = group_limit
        self.cpu                = cpu

        global nwk_dir, taxon_dir, out_dir, all_fasta
        nwk_dir       = Path(self.nwk)
        self.dup_dir  = Path(nwk_dir.parent) / "dup"
        self.seq2dups = {}
        taxon_dir     = Path(self.taxon)
        out_dir       = Path(self.out)
        all_fasta     = Path(self.all_sequence_fasta)

        global process_num, limit_group, limit_otu, LOG_FILE, log_dir
        limit_otu   = self.otu
        limit_group = self.group
        process_num = self.cpu
        LOG_FILE    = log
        log_dir     = Path(self.out) / "log"

        global inparalog_cut, visible_outparalog_cut, invisible_outparalog_cut, detect_invisible
        inparalog_cut            = 1
        visible_outparalog_cut   = 1
        invisible_outparalog_cut = 1
        detect_invisible         = 1
        global seq2spe
        if s2s is None:
            seq2spe = {}
        else:
            seq2spe = s2s
            self.extract_all_duplicates()
        global seq2dups
        seq2dups = self.seq2dups

    def extract_all_duplicates(self):
        dup_files = "{0}/*.tsv".format(self.dup_dir.resolve())
        all_dups  = self.dup_dir / "all.dup"
        cmd       = "cat {0} > {1}".format(dup_files, all_dups)
        subprocess.run(cmd, shell=True)
        with open(str(all_dups), "r") as f:
            for line in f:
                line = line.rstrip()
                line = line.replace(" ", "")
                dup_num, dup_seqs = line.split("\t")
                # In seqkit, the first array of FASTA is left behind. Also, the sequence left in the file listing the removed sequence names is written at the left end
                dup_seqs = dup_seqs.split(",")
                remain   = dup_seqs[0]
                self.seq2dups[remain] = dup_seqs[1:]

    def generate_ortholog_groups(self):
        global numbering_seq, process_num, seq_id
        log_dir.mkdir(parents=True, exist_ok=True)
        numbering_seq = dict()
        with Pool(process_num) as p:
            for d in p.imap_unordered(generate_index, taxon_dir.glob("**/*.faa")):
                numbering_seq.update(d)
        seq_id = {v: k for k, v in numbering_seq.items()}  # reverse dictionary

        global taxons, hierarchies, max_hierarchy
        taxons      = [tax_dir for tax_dir in taxon_dir.glob("**/*") if tax_dir.is_dir()]
        taxons      = sorted(taxons, key=lambda x: (str(x).count("/"), str(x)))
        hierarchies = dict()
        for path in Path(taxon_dir).glob("**/*.faa"):
            tax_numbers       = set(path.parents) & set(taxons)
            tax_numbers       = sorted([taxons.index(path) + 1 for path in tax_numbers])
            hierarchies[path] = tax_numbers
        max_hierarchy = len(max(hierarchies.values(), key=len))

        with Pool(process_num) as pool:
            for _ in pool.imap_unordered(find_ortholog, sorted(nwk_dir.glob("*.nwk"))):
                pass

        otu_in_ortholog = list()  # list with otu which is judged ortholog
        leftovers       = list()  # list with otu which isn't judged ortholog

        if (log_dir / "ortho.log").exists():
            with open(str(log_dir / "ortho.log"), "r") as OT:
                for line in OT.readlines():
                    line = line.replace("\n", "").split("\t")[1:]
                    otu_in_ortholog += line

        #  "ortho.bco" is binary file to use mcl process
        with open(str(log_dir / "ortho.bco"), "wb") as bco:
            pickle.dump(otu_in_ortholog, bco)

        for none_ortholog_file in ["under_{0}OTU.log".format(limit_otu),
                                   "under_{0}group.log".format(limit_otu), "none_cut.log", "not_meet_conditions.log"]:
            if (log_dir / none_ortholog_file).exists():
                with open(str(log_dir / none_ortholog_file), "r") as NO:
                    for line in NO.readlines():
                        leftovers.append(set(line.replace("\n", "").split("\t")[1:]))

        #  "leftovers.bi" is binary file to use mcl process
        with open(str(log_dir / "leftovers.bi"), "wb") as lbi:
            pickle.dump(leftovers, lbi)
