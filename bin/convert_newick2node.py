#!/usr/bin/env python

import subprocess

import re
import argparse

class Newick2Node:
    def __init__(self, newick_file=None, branch_seq="&"):
        self.newick_file = str(newick_file)
        self.branch_seq = branch_seq

    def translate(self, translate_dic, strings):
        ORI = strings
        for k, v in translate_dic.items():
            pat = re.compile(r"\b{0}\b".format(k))
            if pat.search(strings):
                strings = pat.sub(v, strings)
            else:
                continue
        if strings is None:
            print(translate_dic, ORI)
        return strings

    def numbering_duaple(self, newick):
        new_split = newick.split(",")
        otu_numbering = dict()

        is_not_duplicated_name = True
        for n, otu_tuple in enumerate(new_split):
            otu = re.sub("\(|\)", "", otu_tuple)
            if otu in otu_numbering:
                otu_numbering[otu] += 1
            else:
                otu_numbering[otu] = 0
            if otu_numbering[otu] > 0:
                is_not_duplicated_name = False
                new_split[n] = re.sub(otu, "{0}#{1}".format(
                    otu, otu_numbering[otu]), otu_tuple)
        if is_not_duplicated_name:
            pass
        else:
            newick = ",".join(new_split)
        return newick

    def sort_node(self, node_list):
        for ni, node in enumerate(node_list):
            cluster_list = node.split(",")
            for ci, cluster in enumerate(cluster_list):
                cluster = self.branch_seq.join(sorted(cluster.split(self.branch_seq)))
                cluster_list[ci] = cluster

            node = ",".join(sorted(cluster_list, key=lambda x: (len(x), x)))
            if node is None:
                cmd = "echo {0}".format(node_list)
                subprocess.run(cmd, shell=True)

            node_list[ni] = node

        node_list = sorted(node_list, key=lambda x: (len(x), x))
        return node_list
    
    def numbering_otu(self, newick):
        self.numbering_otu_dic = dict()
        numbering_new = newick
        new_split = re.sub("\(|\)", "", newick).split(",")
        
        for num, otu in enumerate (new_split):
            num = str(num).zfill(6)
            self.numbering_otu_dic[num] = otu
            numbering_new = re.sub(otu, num, numbering_new)
        return numbering_new
            

    def newick2node(self, newick=None):
        node_list = list()  # list with tree by node format

        # read tree by newick format
        if newick is None:
            newick = open(self.newick_file, "r").read()

        newick = re.sub(";|\n", "", newick)
        if re.search(":\d\.\d+", newick):
            newick = re.sub(":\d\.\d+", "", newick)

        newick = self.numbering_duaple(newick)
        newick = self.numbering_otu(newick)
        ORI_NEWICK = newick
        # print(newick)
        otu_list = [re.sub(r"[()]", "", otu) for otu in newick.split(",")]
        # print(otu_list)

        # newick to node
        cluster_expression = r"\(([^\(\)]+)\)"
        while re.search(cluster_expression, newick):
            cluster_tuple = re.search(cluster_expression, newick).group()  #
            # if only one "(" in newick, break this roop.
            if newick.count("(") == 1:
                node = ",".join(re.sub("\(|\)", "", cluster_tuple).split(","))
                node_list.append(node)
                break

            in_cluster = re.sub(r"[()]", "", cluster_tuple).split(",")
            otu_in_cluster = [otu for group in in_cluster for otu in group.split(self.branch_seq)]
            out_cluster = [x for x in otu_list if x not in otu_in_cluster]

            # rewritting to node
            if len(out_cluster) > 0:
                node = "{0},{1}".format(
                    ",".join(in_cluster), self.branch_seq.join(out_cluster))
            else:
                node = ",".join(in_cluster)
            if node is None:
                print(ORI_NEWICK)
                raise Exception
            node_list.append(node)

            #  chenge cluster to one group
            cluster_tuple = cluster_tuple.replace(
                "(", "").replace(",", self.branch_seq).replace(")", "").strip()
            newick = re.sub(cluster_expression, cluster_tuple, newick, 1)
        # print(node_list)
        node_list = [self.translate(self.numbering_otu_dic, node) for node in node_list]
        # print(node_list)
        return self.sort_node(node_list)

    def node2newick(self, node_list):
        node_tuple_list = list()  # list with node in tuple
        cluster_dict = dict()     # dict with key:cluster, value:node number
        reverse_dict = dict()     # dict with key:node number, value:cluster

        newick = list()
        rooted_newick = [node for node in node_list if node.count(",") == 1]

        if len(node_list) == 2:
            newick = [node for node in node_list[0].split(",") if self.branch_seq not in node]
            in_cluster = "({0})".format([node for node in node_list[0].split(",") if node not in newick][0].replace(self.branch_seq, ","))
            newick.append(in_cluster)
            newick = "({0})".format(",".join(newick))
            return newick


        while len(newick) != 1:
            for node in node_list:
                min_cluster = [otu for otu in node.split(",") if self.branch_seq not in otu]
                if len(min_cluster) >= 2:
                    num = str(len(cluster_dict))
                    cluster_dict[",".join(min_cluster)] = "node@" + num
                    cluster_dict[self.branch_seq.join(min_cluster)] = "node@" + num
                    reverse_dict["node@" + num] = "({0})".format(",".join(min_cluster))

            node_list = list(map(lambda x: self.translate(cluster_dict, x), node_list))
            rooted_newick = list(map(lambda x: self.translate(x), rooted_newick))
            newick = [node for node in node_list if self.branch_seq not in node]
            
            end_root = [node for node in node_list if node.count(self.branch_seq) == 1]
            if len(end_root) > 0:
                last_cluster = [cluster for cluster in end_root[0].split(",") if self.branch_seq in cluster][0]
                newick = [re.sub(last_cluster, "({0})".format(last_cluster.replace(self.branch_seq, ",")), end_root[0])]

        if len(rooted_newick) > 0:
            newick = list()
            for node in rooted_newick[0].split(","):
                if self.branch_seq in node:
                    node = "({0})".format(node.replace(self.branch_seq, ","))
                newick.append(node)
            newick = "({0})".format(",".join(newick))


        else:
            newick = "({0})".format(newick[0])

        while re.search("node@", newick):
            newick = self.translate(reverse_dict, newick)
        return newick



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script is . . . .")
    parser.add_argument("newick_file", help="A tree file written by newick format")
    args = parser.parse_args()
    newick_file = args.newick_file
    nn = Newick2Node(newick_file)
    print(nn.newick2node())
