#!~/miniconda3/bin/python
# -*- coding: utf-8 -*-

"""
Created on Fri Jul 30 17:50:53 2021

@author: earlyevol
"""
import sys

OG_file = sys.argv[1]
homologues = int(sys.argv[2])


OG_single=[]

#just a list to use later to exclude all OGs with paralogs (not 0 or 1)
one_zero = [1,0]

# open ortholog count file
f = open(OG_file)
for line in f:
    # ignore header line
    if "Orthogroup" in line:
        None
    else:
        #print (line)
        l = line.replace('"',"")
        l = l.replace("'","")
        l1 = l.replace(",\n", "")
        l2 = l1.replace("\n", "")
        OGline = l2.split("\t")

        i = 1
        one_zero_OG=True
        while i < len(OGline)-1 and one_zero_OG == True:
            #print (OGline[i])
            if int(OGline[i]) in one_zero:
                #print (OGline[i])
                one_zero_OG=True
            else:
                one_zero_OG=False
            i+=1
        if one_zero_OG==True:
            #print (OGline)
            if int(OGline[-1]) >= homologues:
                #print (OGline)
                OG_single.append(OGline[0])
outfile = open("SCO_%s" %(str(homologues)), "w")
for I in OG_single:
    outfile.write(I + "\n")


outfile.close()

