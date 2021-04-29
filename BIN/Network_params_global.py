#!/usr/bin/env python
import math
import sys
import os
from collections import defaultdict
# sys.path.append('BIN/')
bin_path = '/lustre/scratch117/cellgen/team297/kt16/BCRSeq/BIN/'
lib_path = '/lustre/scratch117/cellgen/team297/kt16/BCRSeq/LIBRARY/'
if not os.path.exists(bin_path):
  bin_path = os.getcwd() + '/BIN/'
if not os.path.exists(bin_path):
  raise OSError('Cannot locate path to BIN folder')
if not os.path.exists(lib_path):
  lib_path = os.getcwd() + '/LIBRARY/'
if not os.path.exists(lib_path):
  raise OSError('Cannot locate path to LIBRARY folder')
sys.path.append(bin_path)
import Functions
from Functions import *
import re
import time
import commands
import numpy as np
from numpy import outer
from operator import itemgetter, attrgetter
import copy

def Renyi_entropy(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads):
  vrenyi = 0
  crenyi = 0
  tv=totalreads*totalreads*1.0
  tc=totalv*totalv*1.0
  for i in range(0,len(vpoints)):
    vrenyi = vrenyi + vvdf[i]*(vpoints[i]*vpoints[i]/tv)
  for i in range(0,len(cpoints)):
    crenyi = crenyi + cvdf[i]*(cpoints[i]*cpoints[i]/tc)
  return(vrenyi,crenyi)

def Uniq(v):
  C=set(v)
  return list(C)

def VDF (n):
  points=sorted(Uniq(n))
  vdf=[]
  for i in range(0,len(points)):
    vdf.append(n.count(points[i]))
  return (points,vdf)

def Gini_index(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads): 
  (vgini)=Get_Gini(vpoints,vvdf)
  (cgini)=Get_Gini(cpoints,cvdf)
  return(vgini, cgini)

def Get_Gini(n,v):
  values=[]
  for i in range(0,len(n)):
    for j in range(0,v[i]):
      values.append(n[i])
  n = len(values)
  assert(n > 0), 'Empty list of values'
  sortedValues = sorted(values) #Sort smallest to largest
  cumm = [0]
  for i in range(n):
    cumm.append(sum(sortedValues[0:(i + 1)]))
  LorenzPoints = [[], []]
  sumYs = 0           #Some of all y values
  robinHoodIdx = -1   #Robin Hood index max(x_i, y_i)
  for i in range(1, n + 2):
    x = 100.0 * (i - 1)/n
    y = 100.0 * (cumm[i - 1]/float(cumm[n]))
    LorenzPoints[0].append(x)
    LorenzPoints[1].append(y)
    sumYs += y
    maxX_Y = x - y
    if maxX_Y > robinHoodIdx: robinHoodIdx = maxX_Y   
  giniIdx = 100 + (100 - 2 * sumYs)/n #Gini index 
  return(giniIdx/100)

def Get_cluster_vertex_distributions(cluster_file):
  fh=open(cluster_file,"r")
  cluster,vertices=Functions.Tree(), Functions.Tree()
  index, totalc, totalv, totalreads, sizesv, c_sizes,vertices_in_max_cluster = 0,0,0,0,[],{},0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cluster[l[1]][l[2]].value=1
      size = int(l[3])
      sizesv.append(size)
      totalv=totalv+1
      totalreads=totalreads+size
      if(int(l[1])==1):vertices_in_max_cluster=vertices_in_max_cluster+1
      if(l[1] in c_sizes):c_sizes[l[1]]=c_sizes[l[1]]+size
      else:c_sizes[l[1]]=size
  fh.close()
  sizes=[] 
  totalc=len(cluster)
  for c in cluster:
    sizes.append(len(cluster[c]))
  (cpoints,cvdf)=VDF(sizes)
  (vpoints,vvdf)=VDF(sizesv)
  return(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads,c_sizes, vertices_in_max_cluster)
  
def Proportional_measures(c_sizes, totalreads):
  sizes=[]
  for c in c_sizes:
    sizes.append((c,c_sizes[c]))
  s=sorted(sizes, key=itemgetter(1), reverse=True)
  (max_pop, max_1_pop)=(s[0][1]*100.0/totalreads, s[1][1]*100.0/totalreads)
  return(max_pop, max_1_pop)

def Get_differential_stats(cluster_file, id, dir,network_statistics,gene,species,gene_freq_file):
  fh=open(network_statistics,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      totalreads, vgini, cgini, max_pop, max_1_pop, vertices_in_max_cluster, vrenyi, crenyi = l[2],l[4],l[5],l[6],l[7],l[8],l[9],l[10]
  fh.close()
  # largest v (%)
  fh = open(cluster_file, "r")
  max_cluster,clusters, max_id, max_f,max_c ={}, Functions.Tree(),'',0,0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id1,f,c = l[2], int(l[3]),l[1]
      if(f>max_f):max_id, max_f,max_c =id1,f,c
      clusters[c][id1].value=1
  fh.close()
  max_vertex = max_f*100.0/int(totalreads)
  for id1 in clusters[max_c]:
    max_cluster[id1]=max_c
  # number of edges/total poss number and edges in maximum cluster
  edge_file = dir.replace("ANNOTATIONS/","NETWORKS/")+"Edges_"+id+".txt" 
  fh=open(edge_file, "r")
  total_e, max_c_e = 0,0
  for l in fh:
    l=l.strip().split()
    if(l[0] in max_cluster or l[1] in max_cluster):
      max_c_e=max_c_e+1
    total_e=total_e+1
  fh.close()
  # max gene combination
  fh=open(gene_freq_file,"r")
  gene_f = []
  for l in fh:
    l=l.strip().split()
    gene_f.append([l[1], int(l[2])])
  fh.close()
  gene_f=sorted(gene_f, key=itemgetter(1), reverse=True)
  gene_fs = [gene_f[0][1], gene_f[1][1]]
  out="#Id\tAnalysis\tN reads\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\t% Vertices in largest cluster\tVertex Renyi\tCluster Renyi\tMax. vertex (%)\tTotal # edges\t# Edges in max. cluster\tMax. VJ Gene freq (%)\t2nd max. VJ Gene freq\n"
  out=out+id+"\tALL_REPERTOIRE_PARAMS\t"+str(totalreads)+"\t"+str(vgini)+"\t"+str(cgini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\t"+str(vertices_in_max_cluster)+"\t"+str(vrenyi)+"\t"+str(crenyi)+"\t"+str(max_vertex)+"\t"+str(total_e)+"\t"+str(max_c_e)+"\t"+str(gene_fs[0])+"\t"+str(gene_fs[1])+"\n"
  print out
  return()

def Get_network_statistics(cluster_file, id, dir,network_statistics):
  (cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads,c_sizes,vertices_in_max_cluster)=Get_cluster_vertex_distributions(cluster_file)
  (vrenyi,crenyi)=Renyi_entropy(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads)
  (vgini, cgini)=Gini_index(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads)
  (max_pop, max_1_pop)=Proportional_measures(c_sizes, totalreads)
  out="#Id\tAnalysis\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\t% Vertices in largest cluster\tVertex Renyi\tCluster Renyi\n"
  out=out+str(id)+"\tMAX_CLUSTER_REMOVED\t"+str(totalreads)+"\t"+str(totalv)+"\t"+str(vgini)+"\t"+str(cgini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\t"+str(vertices_in_max_cluster*100.0/totalv)+"\t"+str(vrenyi)+"\t"+str(crenyi)+"\n"
  #out="#Id\tAnalysis\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\t% Vertices in largest cluster\n"
  #out=out+str(id)+"\tMAX_CLUSTER_REMOVED\t"+str(totalreads)+"\t"+str(totalv)+"\t"+str(vgini)+"\t"+str(cgini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\t"+str(vertices_in_max_cluster*100.0/totalv)+"\n"
  fh=open(network_statistics, "a")
  fh.write(out)
  fh.close()
  return()

def Get_large_clusters(cluster_file,threshold):
  fh=open(cluster_file,"r")
  total,clusters,clust_size=0,Functions.Tree(),{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      c,id,f = l[1], l[2], int(l[3])
      id = id.split("__")[0]+"__"+l[3]
      total = total+f
      if(c in clust_size):clust_size[c]=clust_size[c]+f
      else:clust_size[c]=f
      clusters[c][id][f].value=1
  fh.close()
  sig_clust_info,ids = [],{}
  for c in clust_size:
    prop= clust_size[c]*100.0/total
    if(prop>=threshold):
      perc_network, number_vertices, number_reads, max_seq, reads_in_max_seq = prop, len(clusters[c]),0,'',0
      for id in clusters[c]:
        ids[id]=c
        for f in clusters[c][id]:
          number_reads = number_reads+f
          if(f>reads_in_max_seq):
            max_seq, reads_in_max_seq = id,f
      sig_clust_info.append([c,perc_network, number_vertices, number_reads, max_seq, reads_in_max_seq])
  del clusters
  sig_clust_info =sorted(sig_clust_info,key=itemgetter(1),reverse=True)
  return(sig_clust_info,ids)

def Get_max_sequences(seq_file,ids):
  fh=open(seq_file,"r")
  seqs={}
  for header,sequence in fasta_iterator(fh):
    if(header in ids):
      seqs[header.split("__")[0]]=sequence
  fh.close()
  return(seqs)

def Initalise_file(file):
  fh=open(file,"w")
  fh.close()
  return()

def Write_out(out, file):
  fh = open (file,"a")
  fh.write(out)
  fh.close()
  return()

###########################
cluster_file = sys.argv[1]              ### clustering file
network_statistics = sys.argv[2]        ### output file
id = sys.argv[3]
dir = sys.argv[4]
######################### Commands
#Initalise_file(network_statistics)
Get_network_statistics(cluster_file, id, dir,network_statistics)




