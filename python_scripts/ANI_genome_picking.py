#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 12:36:23 2021

@author: earl
"""
"""
Usage:
    ANI_genome_picking.py fast_ANI.out genome_list number_of_shortlist_clusters
        or
    ANI_genome_picking.py TEST
"""
from collections import Counter
import copy
import sys

tester = sys.argv[1]
if tester == "TEST":
    ANI_out="/Users/earlm/gits/ortho_phylo/TESTER/ANI_out"
    genome_names="/Users/earlm/gits/ortho_phylo/TESTER/genome_names"
    n_clusters=10
elif tester == "local":
    ANI_out="/Users/earlm/Projects/orthophylo/ecoli/ANI/ANI_out"
    genome_names="/Users/earlm/Projects/orthophylo/ecoli/ANI/genome_names"
    n_clusters=10
else:
    ANI_out=sys.argv[1]
    genome_names=sys.argv[2]
    n_clusters=int(sys.argv[3])

genome_list = set()
f = open(genome_names)
for line in f:
    l = line.replace('"',"")
    l = l.replace("'","")
    l1 = l.replace(")", "")
    l1 = l1.replace("(", "")
    l2 = l1.replace("\t\n", "")
    l2 = l2.replace("\n", "")
    genome_list.add(l2)

full_dict = {}
for I in genome_list:
    full_dict[I] = {}
    for J in genome_list:
        # set the minimum ANI value for comparisons that are too divergent
        full_dict[I][J] = int(50)

f = open(ANI_out)
for line in f:
    #print (line)
    l = line.replace("../genomes/","")
    l = line.replace('"',"")
    l = l.replace("'","")
    l1 = l.replace(")", "")
    l1 = l1.replace("(", "")
    l2 = l1.replace("\t\n", "")
    l2 = l2.replace("\n", "")
    l3 = l2.replace("\t\t", "\t")
    l4 = l3.split("\t")
    #store the names of selected genomes in a dict with their avg ANI to other genomes
    if l4[0] != l4[1]:
        full_dict[l4[0]][l4[1]] = float(l4[2])
f.close()


#print all pairwise ANIs that were not calculated (too divergent)
num_too_divergent = 0
total_comparison = 0
for I in full_dict:
    for J in full_dict:
        total_comparison += 1
        if full_dict[I][J] == 50:
            print ("I = " + I + " and J = " + J)
            num_too_divergent += 1
print ("##################################################")
print (str(num_too_divergent) + " comparisons were too divergent to calculate ANI")
print ("Of " + str(total_comparison) + " total comparisons")
print ("##################################################")

#rm self entries
#for key in full_dict:
#    del full_dict[key][key]

full_dict_copy = copy.deepcopy(full_dict)

#dont remember what this is for...
def average_dict_values(in_dict):
    res = 0
    for val in in_dict.values():
        res += val
    # using len() to get total keys for mean computation
    mean = res / len(in_dict)
    return mean

def cluster(in_dict,first_iteration):
    ANI_max={}
    for key1,val in full_dict.items():
        #print (key1,val)
        max_key = max(val, key=val.get)
        ANI_max[max_key,key1]=full_dict[key1][max_key]

    max_ANI_pair = max(ANI_max, key=ANI_max.get)
    max_val = ANI_max[max_ANI_pair]
    #sum ANI from both in pair and create a dictionary
    sums = dict(Counter(full_dict[max_ANI_pair[0]]) + Counter(full_dict[max_ANI_pair[1]]))
    # create a new "full_dict entry for the clustered clade (max_ANI_pair)
    #  Average distances, add 1/2 the max_ANI_pair's distance to get new distances
    full_dict[max_ANI_pair] = {k: (sums[k] / float((k in full_dict[max_ANI_pair[0]]) + (k in full_dict[max_ANI_pair[1]])))+((100-max_val)/2) for k in sums}
    print (max_ANI_pair)
    #add other values for other side of matrix for merged genomes
    for I in full_dict[max_ANI_pair]:
        full_dict[I][max_ANI_pair] = full_dict[max_ANI_pair][I]

    #print (full_dict.keys())
    for I in max_ANI_pair:
        #print (full_dict[I])
        del full_dict[I]
        #del full_dict[max_ANI_pair][I]
    for I in max_ANI_pair:
        for key in full_dict:
            if I in full_dict[key]:
                #print (key)
                #print ("This is I   " + str(I))
                del full_dict[key][I]
    return full_dict
    return ANI_max
    return max_ANI_pair


first_iteration=0
while len(full_dict) > n_clusters:
    cluster(full_dict,first_iteration)
    first_iteration+=1

NN_distances={}
for I in full_dict:
    NN_distances[I]=str(I).replace(")","").replace("(","").replace(",,",",").replace("'","").replace(" ","").split(",")

dict_dist={}
for I in NN_distances:
    dict_dist[I]={}
    for J in NN_distances[I]:
        dict_dist[I][J]=average_dict_values(full_dict_copy[J])

species_list=[]
for I in dict_dist:
     species_list.append(max(dict_dist[I], key=dict_dist[I].get))

outfile=open("Species_shortlist", 'w')
for k in species_list:
    outline=str(k)
    outfile.write(outline + "\n")
outfile.close()
