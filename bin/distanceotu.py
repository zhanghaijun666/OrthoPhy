import re
import argparse
import numpy as np
import pandas as pd
from decimal import Decimal
from itertools import combinations

from bin.convert_newick2node import Newick2Node

newick2node = Newick2Node(branch_seq="&").newick2node

class DistanceOTU:
    def __init__(self, newick_tree, node_list=None):
        self.nwk_tree = newick_tree
        self.node_list = node_list
        self.distance_table()
        self.key_stone()
            
    
    def distance_table(self, nwk_tree=None):
        if nwk_tree is None:
            nwk_tree = self.nwk_tree
        tree = nwk_tree
        otu_list = sorted([otu for otu in re.sub("\(|\)|:\d+\.\d+|;", "", tree).split(",")], 
                         key=lambda x: (len(x), x))
        otu_dic = {otu: [None]*len(otu_list) for otu in otu_list}
        length_dic = {}
        df = pd.DataFrame(data=otu_dic, index=otu_list)

        match =  "\(([^\(\)]+)\)"
        while re.search(match, tree):
            if tree.count("(") == 1:
                break
            node = re.search(match, tree).group(1)
            length = Decimal(re.search(r"\({0}\):(\d+\.\d+)".format(node), tree).group(1))
            nodes = tuple(sorted(re.split(",|/", (re.sub(":\d+\.\d+|<|>", "", node)))))
            length_dic[nodes] = length
            new_node = "<{0}>".format(node).replace(",", "/")
            tree = re.sub(r"\({0}\)".format(node), new_node, tree)
        for otu in re.sub("\(|;|\):\d+\.\d+|\)", "", nwk_tree).split(","):
            otu, length = otu.split(":")
            df.loc[otu, otu] = Decimal(length)
        
        if self.node_list is None:
            self.node_list = newick2node(nwk_tree)
            for node in self.node_list:
                if node is None:
                    print(nwk_tree)
                    raise Exception
        for node in self.node_list:
            cluster_list = node.split(",")
            internal_nodes = tuple(sorted(cluster_list[0].split("&") + cluster_list[1].split("&")))
            if len(internal_nodes) == 2:
                f_otu, t_otu = internal_nodes
                df.loc[f_otu, t_otu] = df.loc[f_otu, f_otu] + df.loc[t_otu, t_otu]
                df.loc[t_otu, f_otu] = df.loc[f_otu, f_otu] + df.loc[t_otu, t_otu]

            if internal_nodes in length_dic:
                outgroup = tuple(sorted(cluster_list[2].split("&")))
                for f_otu in internal_nodes:
                    for t_otu in outgroup:
                        
                        outnodes = []
                        for nodes, length in length_dic.items():
                            if (f_otu in nodes) and (t_otu not in nodes):
                                outnodes.append(length)
                            elif (f_otu not in nodes) and (t_otu in nodes):
                                outnodes.append(length)

                        if len(outnodes) > 0:
                            add_length = sum(outnodes)
                        else:
                            add_length = 0
                        
                        add_length += df.loc[f_otu, f_otu]
                        add_length += df.loc[t_otu, t_otu]
                        df.loc[t_otu, f_otu] = add_length
                        df.loc[f_otu, t_otu] = add_length

        for f_otu in otu_list:
            for t_otu in otu_list:
                if df.loc[f_otu, t_otu] is None:
                    df.loc[f_otu, t_otu] = df.loc[f_otu, f_otu] + df.loc[t_otu, t_otu]
        self.distance = df
        
    def list_to_flat(self, list_list, del_duplication=0, tuple_ok=0):
        """
        list which has list -> list which has no list

        parameters
        ----------
        list_list : list which has list

        returns
        -------list
        flat_list : list which has no list
            :param del_duplication:
            remove duplication (use:1, not use:0, [default:0])
        """
        flat_list = []  # list wiich has no list
        for flat in list_list:

            if type(flat) is set:
                flat = list(flat)
            if type(flat) is list:
                flat_list.extend(self.list_to_flat(flat))
            elif tuple_ok and type(flat) is tuple:
                flat_list.extend(self.list_to_flat(flat))
            else:
                flat_list.append(flat)
        if del_duplication:
            flat_list = list(set(flat_list))
        return flat_list
        
    def key_stone(self):
        key_cluster = [node.split(",") for node in self.node_list]
        key_cluster = [[(cluster.split("&")) for cluster in node] for node in key_cluster]
        key_cluster = [list(combinations(cluster, 2)) for cluster in key_cluster]
        key_cluster = self.list_to_flat(key_cluster)
        self.key_cluster = key_cluster
        
    def otu_distance(self, leafa, leafk, distance=None):
        if distance is None:
            distance = self.distance
        elif type(distance) is str:
            distance = self.distance_table(tree)

        length = distance.loc[leafa, leafk] - (distance.loc[leafa, leafa] + distance.loc[leafk, leafk])
        return length
    
    def cluster_distance(self, cluster0, cluster1, another_key_OTU=None, midway_expression=False, distance=None):
        cluster0 = cluster0[:]
        cluster1 = cluster1[:]
        if distance is None:
            distance = self.distance
        elif type(distance) is str:
            distance = self.distance_table(tree)
        
        if type(cluster0) in (list, tuple):
            cluster0 = sorted(list(cluster0))
            cluster1 = sorted(list(cluster1))
        elif type(cluster0) == str:
            cluster0 = sorted(cluster0.split("&"))
            cluster1 = sorted(cluster1.split("&"))
        else:
            print("Error")
            
        if another_key_OTU is None:
            if type(self.node_list[0]) is str:
                another_list = [node.split(",") for node in self.node_list]
            elif type(self.node_list[0]) is list:
                another_list = self.node_list

            another_list = [[sorted(cluster.split("&")) for cluster in node] for node in another_list]
            ak_OTU = list(map((lambda x: x[2]  if cluster0 in x else None), another_list))
            ak_OTU = [x for x in ak_OTU if x][0][0]

        else:
            ak_OTU = another_key_OTU.split("&")[0]
            
        if len(cluster0) == 1:
            a0, k0 = cluster0[0], cluster0[0]
        if len(cluster1) == 1:
            a1, k1 = cluster1[0], cluster1[0]
        
        flag = 0
        for key in self.key_cluster:
            if cluster0 == sorted(key[0] + key[1]):
                a0, k0 = map(lambda x: x[0], key)
                flag += 1
            if cluster1 == sorted(key[0] + key[1]):
                a1, k1 = map(lambda x: x[0], key)
                flag += 1
            if flag == 2:
                break
        length0 = (self.otu_distance(ak_OTU, a0) + self.otu_distance(ak_OTU, k0) - self.otu_distance(a0, k0)) / 2
        length1 = (self.otu_distance(ak_OTU, a1) + self.otu_distance(ak_OTU, k1) - self.otu_distance(a1, k1)) / 2
        length01 = ((self.otu_distance(ak_OTU, a0) + self.otu_distance(ak_OTU, k1)) - self.otu_distance(a0, k1)) / 2
        length0 -= length01
        length1 -= length01
        
        if len(cluster0) == 1:
            length0 *= 2
        if len(cluster1) == 1:
            length1 *= 2
        
        if midway_expression:
            return length0, length1
        else:
            return length0 + length1

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script is . . . .")
    parser.add_argument("newick_file", help="A tree file written by newick format with branch length.")
    parser.add_argument("f_otu", help="One of OTU or cluster which you want to know the length.(if it is cluster, you join OTU by ',')")
    parser.add_argument("t_otu", help="The other OTU or cluster which you want to know the length.")
    parser.add_argument("--midway_expression", "-me", help="If you want to know the length which each branch, you use this.", action='store_true')
    args = parser.parse_args()
    newick_file = args.newick_file
    f_otu        = args.f_otu
    t_otu        = args.t_otu
    midway_expression = args.midway_expression
    print(midway_expression)
    
    if "," in f_otu:
        f_otu = f_otu.split(",")
    if "," in t_otu:
        t_otu = t_otu.split(",")
    tree = open(newick_file, "r").read()
	
    distance = DistanceOTU(tree)
    
    if type(f_otu) is list or type(t_otu) is list:
        cluster_distance = distance.cluster_distance(f_otu, t_otu, midway_expression=midway_expression)
    else:
        cluster_distance = distance.otu_distance(f_otu, t_otu)
        
    print(cluster_distance)
