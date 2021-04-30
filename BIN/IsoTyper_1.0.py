#!/usr/bin/python
import math
import sys
from collections import defaultdict
sys.path.append('./')
import Functions
from Functions import *
import os
import subprocess
import sys
from operator import itemgetter

from operator import add
import numpy as np
import copy
import random
from collections import Counter
from tqdm import tqdm

class Tree(defaultdict):
  def __init__(self, value=None):
    super(Tree, self).__init__(Tree)
    self.value = value

def fasta_iterator(fh):
  while True:
    line = fh.readline()
    if line.startswith('>'): break	
    elif(len(line)==0):break
  while True:
    header = line[1:-1].rstrip()
    sequence = fh.readline().rstrip()
    while True:
      line = fh.readline()
      if not line: break
      if line.startswith('>'): break
      sequence += line.rstrip()
    yield(header, sequence)
    if not line: return

def Trim_sequences(s1,s2,l1,l2,max_mm):
  p, segment=0, len(s1)/max_mm
  for i in range(segment, len(s1)):
    sample1=s1[i-segment:i]
    fin = s2.count(sample1)
    if(fin!=-1):
      inde=s2.find(sample1)
      if (inde!=-1):
        if(inde>i):
          if(inde-i <=max_mm):s2=s2[inde-i:l2]
        else:
          if(i-inde <=max_mm):s1=s1[i-inde:l1]
        min_len=min([len(s1),len(s2)])
        if ((max([len(s1),len(s2)]) - min_len) <max_mm):
          s1=s1[0:min_len]
          s2=s2[0:min_len]
    #      print "\t",s1
    #      print "\t",s2
          p=1
  return (s1, s2, p)

def Get_annotations(annot_file):
  CDR3s,mutations = Tree(),{}
  fh=open(annot_file,"r")
  for l in fh:
    l=l.strip().split("\t")
    if(l[0]!='Sequence number'):
      if(len(l)>=21 and l[2].count('productive')!=0):
        id,v_muts, j_muts,cdr3 = l[1].split("__")[0],l[6], l[12], l[20].replace("#","")
        v,j = l[3], l[9]
        if(v_muts.count('nt')!=0 and j_muts.count('nt')!=0 and len(cdr3)>1):
          v = v.split()[1].split("*")[0]
          j = j.split()[1].split("*")[0]
          cdr3 = cdr3.split()[0]+":"+v+":"+j
          v_muts, j_muts = list(map(int,v_muts.split()[0].split("/"))), list(map(int,j_muts.split()[0].split("/")))
          v_muts, j_muts = v_muts[1]-v_muts[0], j_muts[1]-j_muts[0]
          mutations[id] = [v_muts, j_muts,cdr3, v+"\t"+j]
          CDR3s[cdr3][id].value = 1
  fh.close()
  return(CDR3s,mutations)

def Deconvolute_same_array (tree):
  decon = Tree()
  inv = {}
  inde = 0
  for i in tree:
    if (i not in inv):
      inde = inde+1
      decon[i][i].value = 1
      inv[i]=i
      array = []
      for j in tree[i]:
        decon[i][j].value = 1
        inv[j]=i
        array.append(j)
      array_new=[]
      found = 0
      if (len(array)>0):
        found = 1
      while (found == 1):
        array_new =[]
        for k in array:
          for l in tree[k]:
            if (l not in inv):
              array_new.append(l)
              inv[l]=i
              decon[i][l].value = 1
        array = array_new
        if (len(array_new)==0):
          found = 0
  return (decon, inv)

def Merge_clusters(CDR3s,mutations, cluster_file, merged_cluster_file,seqs,freqs): 
  clusters,CDR3_clusters,cluster_ids, max_id = Tree(),Tree(),{},{}
  fh=open(cluster_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cluster, id = l[1],l[2].split("__")[0]
      if(id in mutations):
        clusters[cluster][id].value = 1
        cluster_ids[id] = cluster
        CDR3_clusters[mutations[id][2]][cluster].value = 1
        f = sum(freqs[id])
        if(cluster in max_id):
          if(f>max_id[cluster][1]):max_id[cluster] = [seqs[id],f]
        else:max_id[cluster] = [seqs[id],f]
  fh.close()
  same,max_mm = Tree(),8
  for cdr3 in CDR3_clusters:
    #if(len(CDR3_clusters[cdr3])>1):
    #  print ''
    #  for c in CDR3_clusters[cdr3]:
    #    print cdr3, c
    cs = []
    for c in CDR3_clusters[cdr3]:
      cs.append(c)
    for i in range(0,len(cs)):
      max1 = max_id[cs[i]][0]
      for j in range(i,len(cs)):
        if(i<j):
          max2 = max_id[cs[j]][0]
          aln1, aln2, p = Trim_sequences(max1, max2, len(max1), len(max2), max_mm)
          if(p==1):
            same[cs[i]][cs[j]].value = 1
            same[cs[j]][cs[i]].value = 1
  decon, inv = Deconvolute_same_array(same)
  fh=open(merged_cluster_file, "w")
  fh.close()
  out, ind = '#index\tmerged cluster id\tseq id\tv mutations\tj mutations\tCDR3 region\tv gene\tj gene\told cluster id\n',0
  done_clusters = {}
  inde = 0
  for c in decon:
    for c1 in decon[c]:
      done_clusters[c1] = 1
      for id in clusters[c1]:
        mut = mutations[id]
        out=out+str(inde)+"\tM"+str(c)+"\t"+id+"\t"+str(mut[0])+"\t"+str(mut[1])+"\t"+str(mut[2])+"\t"+str(mut[3])+"\t"+str(c1)+"\n"
        ind = ind+1
        if(ind>100):
          Write_out(out, merged_cluster_file)
          out, ind = '',0
  Write_out(out, merged_cluster_file)
  out, ind = '',0
  for c in clusters:
    inde = inde+1
    if(c not in done_clusters):
      for id in clusters[c]:
        mut = mutations[id]
        out=out+str(inde)+"\tM"+str(c)+"\t"+id+"\t"+str(mut[0])+"\t"+str(mut[1])+"\t"+str(mut[2])+"\t"+str(mut[3])+"\t"+str(c)+"\n"
        ind = ind+1
        if(ind>100):
          Write_out(out, merged_cluster_file)
          out, ind = '',0
  Write_out(out, merged_cluster_file)
  out, ind = '',0
  return()

def Compare_mutations_IMGT_internal(annot_file_internal,annot_file_IMGT,mutation_correlation_file,sample,seq_file):
  CDR3s,mutations = Get_annotations(annot_file_IMGT)
  annots = Get_annotation(annot_file_internal)
  seqs,freqs,alias = Get_sequences(seq_file)
  mut_correlation = {}
  for id in mutations:
    if(id in annots):
      name = str(annots[id][2])+"\t"+str(mutations[id][0])
      if(mutations[id][0]-annots[id][2]>1):
        #print ">", id, mutations[id][0],annots[id][2]
        #print seqs[id]
        print(id+"\t"+str(mutations[id][0])+"\t"+str(annots[id][2]))
      if(name in mut_correlation):mut_correlation[name]= mut_correlation[name]+1
      else:mut_correlation[name]=1
  out="#sample\tinternal mutations\tIMGT mutations\tfrequecy\n"
  for n in mut_correlation:
    out=out+sample+"\t"+n+"\t"+str(mut_correlation[n])+"\n"
  fh=open(mutation_correlation_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_developmental_classification(developmental_classification_file):
  classification = {}
  fh = open(developmental_classification_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      classification[l[1]] = l[2].split("|")
  fh.close()
  return(classification)

def Get_sequences(seq_file):
  fh=open(seq_file,"r")
  seqs,freqs,alias = {},{},{}
  for header,sequence in fasta_iterator(fh):
    seqs[header.split("__")[0]]=sequence
    f = list(map(int, header.split("__")[1].split("|")[0].split("_")))
    freqs[header.split("__")[0]] = f
    alias[header.split("__")[0]] = header
  fh.close()
  return(seqs, freqs, alias)

def Get_closest_links(clusters, c,seqs,freqs,merged_cluster_nearest_neighbours,alias):
  local_clusters,max_freq, max_seqs = [],{},[]
  for c1 in clusters[c]:
    local_clusters.append(c1)
  max_mm = 8
  out,ind = '',0
  for i1 in range(0,len(local_clusters)):
    min_distance = [100000,'','']
    for j in range(0,len(local_clusters)):
      if(i1!=j):
        for id1 in  clusters[c][local_clusters[i1]]:
          s1 = seqs[id1]
          for id2 in  clusters[c][local_clusters[j]]:
            s2 = seqs[id2]
            aln1, aln2, p = Trim_sequences(s1, s2, len(s1), len(s2), max_mm)
            if(p==1):
              d = [i for i,(a1i,a2i)  in enumerate(zip(aln1,aln2)) if a1i!=a2i]
              if(len(d)<min_distance[0]):
                min_distance = [len(d), id1,id2]
    if(min_distance[0]<8):
      out=out+alias[min_distance[1]]+"\t"+alias[min_distance[2]]+"\t"+str(min_distance[0])+"\n"
      ind = ind+1
      if(ind>100):
        Write_out(out, merged_cluster_nearest_neighbours)
        out, ind = '',0
  Write_out(out, merged_cluster_nearest_neighbours)
  out, ind = '',0
  return()

def something():
  aligned_seqs = {}
  kmers,max_seqs = {},seqs[max_freq[c][0]]
  indent = len(max_seq)/max_mm
  excluded_clusters={}
  for i in range(indent,len(max_seq)+1):
    kmers[max_seq[i-indent:i]] = i-indent
  fin = 0
  for c1 in clusters[c]:
    for id in clusters[c][c1]:
      seq,seq_aligned = seqs[id],''
      passed = 1
      for kmer in kmers:
        if(seq.count(kmer)!=0):
          index=seq.find(kmer)
          d = index-kmers[kmer]
          if(d==0):
            seq_aligned = seq
            break
          elif(max([d,-d])<max_mm):
            if(d<0):
              seq_aligned = "-"*(-1*d)+seq
            else:seq_aligned = seq[d:len(seq)]
            break
          else:passed = 0
      len_diff = len(seq_aligned)-len(max_seq)
      if(passed==1):
        if(max([len_diff,-len_diff])<max_mm):
          aligned_seqs[id] = seq_aligned
        else:passed = 0
      if(passed==0):
        excluded_clusters[id] = 1 
  print(len(aligned_seqs), len(excluded_clusters))
  return()

def Calculate_merge_cluster_edges(seqs, merged_cluster_file,merged_cluster_nearest_neighbours, edge_file,freqs,alias):
  fh=open(merged_cluster_file,"r")
  clusters,cluster_ids = Tree(),{}
  mutations = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cluster_merged, id, cluster_old = l[1],l[2],l[8]
      v_mut, j_mut = int(l[3]), int(l[4])
      clusters[cluster_merged][cluster_old][id].value = 1
      mutations[id] = [v_mut, j_mut]
  fh.close()
  ## also get sequence with lowest distnace to GL
  closest_GL_sequences={}
  fh=open(merged_cluster_nearest_neighbours, "w")
  fh.close()
  edges = {}
  out,ind = "#Id1\tId2\tEdge_length\n",0
  fh=open(merged_cluster_nearest_neighbours,"w")
  fh.write(out)
  fh.close()
  ind1 = 0
  for c in clusters:
    ind1 = ind1+1
    if (len(clusters[c])>1):
      print("\r",c, len(clusters[c]), ind1*100.0/len(clusters), end=' ')
      Get_closest_links(clusters, c,seqs,freqs,merged_cluster_nearest_neighbours,alias)
  os.system("cat "+edge_file+" >> "+merged_cluster_nearest_neighbours)
  return()

def Get_nearest_neighbours(edge_file, merged_cluster_file, seq_file,merged_cluster_nearest_neighbours):
  seqs,freqs,alias = Get_sequences(seq_file)
  Calculate_merge_cluster_edges(seqs, merged_cluster_file,merged_cluster_nearest_neighbours, edge_file,freqs,alias)
  return()

def Merge_clusters_on_CDR3_region(seq_file, annot_file, tmp_file1, merged_cluster_file, merged_cluster_nearest_neighbours,cluster_file,edge_file):
  CDR3s,mutations = Get_annotations(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  Merge_clusters(CDR3s,mutations, cluster_file, merged_cluster_file,seqs,freqs)
  Get_nearest_neighbours(edge_file, merged_cluster_file, seq_file,merged_cluster_nearest_neighbours)
  return()

def Get_probability_connected_nodes(cluster_file, connectivity_probability_file,seq_file,sample):
  seqs,freqs,alias = Get_sequences(seq_file)
  reclassify, classes,classification = Get_reclassification(alias)
  fh=open(cluster_file, "r")
  clusters = Tree()
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cluster, id = l[1],l[2].split("__")[0]
      f = freqs[id]
      non_zero = [i for i in range(len(f)) if f[i]!=0]
      for i in non_zero:
        clusters[classification[i].split("*")[0]][cluster][id][f[i]].value = 1
  fh.close()
  all_chains, linked_chains, unlinked_chains = {},{},{}
  for chain in clusters:
    all_chains[chain], linked_chains[chain], unlinked_chains[chain] = [],[],[]
  for chain in clusters:
    for c in clusters[chain]:
      f,n = 0,0
      for id in clusters[chain][c]:
        for f1 in clusters[chain][c][id]:
          f = f+f1
      if(len(clusters[chain][c])==1):unlinked_chains[chain]= unlinked_chains[chain]+[f]
      else:linked_chains[chain] = linked_chains[chain] + [f]
      all_chains[chain] = all_chains[chain]+[f]
  out = "#ID\tChain\tall BCRs\tlinked BCRs\tunlinked BCRS\n"
  for chain in clusters:
    out=out+sample+"\t"+chain+"\t"+str(sum(all_chains[chain]))+"\t"+str(sum(linked_chains[chain]))+"\t"+str(sum(unlinked_chains[chain]))+"\n"
  fh=open(connectivity_probability_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_aa_CDR(annot_file4):
  cdrs = {}
  fh=open(annot_file4,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split("\t")
      if(len(l)>=18):
        cdr3 = l[13].upper()
        cdr2 = l[11].upper()
        cdrs[l[1].split("__")[0]] = [cdr2,cdr3]
  fh.close()
  return(cdrs)

def Get_CDR3s(annot_file):
  fh=open(annot_file,"r")
  cdr3s = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split("\t")
      if(len(l)>=20):
        cdr3 = l[20].upper()
        if(len(cdr3)>3):
          cdr3s[l[1].split("__")[0]] = cdr3
  fh.close()
  return(cdr3s)

def Get_annotation(annot_file):
  fh=open(annot_file,"r")
  annots = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(len(l)>=19):
        v,j,v_mm,j_mm = l[1],l[13],int(l[len(l)-4]), int(l[len(l)-3])
        if(v_mm+j_mm>60):print(l)
        annots[l[0].split("__")[0]] = [v,j,v_mm,j_mm]
  fh.close()
  return(annots)

def Get_cluster_occupancy(merged_cluster_file):
  clusters, inv_clusters= Tree(),{}
  fh=open(merged_cluster_file, "r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cluster, id = l[1],l[2].split("__")[0]
      clusters[cluster][id].value = 1
      inv_clusters[id]=cluster
  fh.close()
  return(clusters, inv_clusters)

def Assess_transitional_frequencies(merged_cluster_nearest_neighbours,clusters, inv_clusters,gl_mutations,transition_classification_file_directed):
  seqs,freqs,alias = Get_sequences(seq_file)
  reclassify, classes,classification = Get_reclassification(alias)
  fh=open(merged_cluster_nearest_neighbours,"r")
  out,ind = '',0
  mutations_unweighted, mutations_weighted, mutations_most_likely = {},{},{}
  all_transitions = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id1, id2, dist = l[0], l[1], int(l[2])
      if(id1.split("__")[0] in gl_mutations and id2.split("__")[0] in gl_mutations):
        glmut1, glmut2 = sum(gl_mutations[id1.split("__")[0]]), sum(gl_mutations[id2.split("__")[0]])
        if(glmut1<glmut2):ida, idb = id1, id2
        else:ida, idb = id2, id1
        fa, fb = freqs[ida.split("__")[0]], freqs[idb.split("__")[0]]
        f1, f2 = Reclassify_sequences(fa, reclassify, classification,classes), Reclassify_sequences(fb, reclassify, classification,classes)
        non_zero1, non_zero2 = [i for i in range(len(f1)) if f1[i]!=0], [i for i in range(len(f2)) if f2[i]!=0]
   #     if(non_zero1==non_zero2 and len(non_zero1)==1):
   #       for i in non_zero1:
   #         name = classes[i]+"\t"+classes[i]
   #         all_transitions[name] = 1
   #         for X in [mutations_unweighted, mutations_weighted, mutations_most_likely]:
   #           if(name in X):X[name] = X[name] = 1
   #           else: X[name] = 1
   #     else:
        if(1==1):
          for i in non_zero1:
            for j in non_zero2:
              weight1 = (f1[i]*1.0/sum(f1))*(f2[j]*1.0/sum(f2))
              weight2 = (1.0/len(non_zero1)) *(1.0/len(non_zero2))
              name = classes[i]+"\t"+classes[j]
              print(classes[i]+"|"+classes[j], "\t",np.array(classes)[non_zero1],np.array(classes)[non_zero2])
              if(name in mutations_weighted):mutations_weighted[name] = mutations_weighted[name]+ weight1
              else:mutations_weighted[name]= weight1
              if(name in mutations_unweighted):mutations_unweighted[name] = mutations_unweighted[name] + weight2
              else:mutations_unweighted[name] = weight2
              all_transitions[name] = 1
          ## most likely
          for j in non_zero2:
            if(j in non_zero1):
              name = classes[j]+"\t"+classes[j]
              if(name in mutations_most_likely): mutations_most_likely[name] = mutations_most_likely[name] + 1
              else:mutations_most_likely[name] = 1
              all_transitions[name] = 1
            else:
              for i in non_zero1:
                name = classes[i]+"\t"+classes[j]
                weight = 1.0/len(non_zero1)
                if(name in mutations_most_likely): mutations_most_likely[name] = mutations_most_likely[name] + weight
                else:mutations_most_likely[name] = weight
                all_transitions[name] = 1
  fh.close()
  out = "#From\tTo\tmutations_unweighted\tmutations_weighted\tmutations_most_likely\n"
  for name in all_transitions:
    out = out+name
    for X in [mutations_unweighted, mutations_weighted, mutations_most_likely]:
      if(name in X):
        out=out+"\t"+str(X[name])
      else:out=out+"\t0"
    out=out+"\n"
  fh=open(transition_classification_file_directed, "w")
  fh.write(out)
  fh.close()
  return()

def Get_reclassification(alias):
  for id in alias:
    break
  classification = alias[id].split("|")[1].split("_")
  reclassify,classes = {},[]
  for c in classification:
    c1 = c.split("*")[0]
    if(c1 not in classes):classes.append(c1)
    reclassify[c] = classes.index(c1)
  return(reclassify, classes,classification)

def Assess_parental_frequencies1(annots,parental_classification_file_directed,edge_file,freqs,alias,sample,cluster_file):
  clusters, inv_clusters,gl_mutations = Get_cluster_occupancy(cluster_file)
  for id in alias:
    classes = alias[id].split("|")[1].split("_")
    break
  print(classes)
  for c in clusters:
    if(len(clusters[c])>1):
      f = [0]*len(classes)
      for id in clusters[c]:
        f = list(map(add, f, freqs[id]))
      class_usage = Get_class_usage(f,classes)
      print(class_usage)


def Assess_parental_frequencies(annots,parental_classification_file_directed,edge_file,freqs,alias,sample,cluster_file):
  total = ''
  for id in freqs:
    if(len(total)==0):total=freqs[id]
    else:total = list(map(add, total, freqs[id]))
  print(total, sum(total))
  totals = {}
  for id in alias:
    classes = alias[id].split("|")[1].split("_")
    break
  for i in range(0,len(classes)):
    clas = classes[i].split("*")[0]
    if(clas in totals):totals[clas] = totals[clas]+total[i]
    else:totals[clas]=total[i]
  fh = open(edge_file,"r")
  parent_counts,child_counts,connected_counts={},{},{}
  done_par, done_chil,done_conn = {},{},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id1,id2 = l[0].split("__")[0],l[1].split("__")[0]
      if(id1 in annots and id2 in annots):
        mm1,mm2 = annots[id1][2],annots[id2][2]
        if(mm1!=mm2):
          ids,mms = [id1,id2],[mm1,mm2]
          if(mm1>mm2):
            ids,mms = [id2,id1],[mm2,mm1]
          f1 = freqs[ids[0]]
          #class_usage = Get_class_usage(freqs[ids[0]],classes)
          #print class_usage
          classes = alias[ids[0]].split("|")[1].split("_")
          nz = [i for i in range(len(f1)) if f1[i]!=0]
          for i1 in range(0,len(nz)):
            i = nz[i1]
            clas = classes[i].split("*")[0]
            if(ids[0] not in done_par):
              if(clas in parent_counts):parent_counts[clas] = parent_counts[clas]+f1[i]
              else:parent_counts[clas]=f1[i]
            if(ids[0] not in done_conn):
              if(clas in connected_counts):connected_counts[clas] = connected_counts[clas]+f1[i]
              else:connected_counts[clas]=f1[i]
          done_par[ids[0]],done_conn[ids[0]] = 1,1
          f1 = freqs[ids[1]]
          nz = [i for i in range(len(f1)) if f1[i]!=0]
          for i1 in range(0,len(nz)):
            i = nz[i1]
            clas = classes[i].split("*")[0]
            if(ids[1] not in done_chil):
              if(clas in child_counts):child_counts[clas]=child_counts[clas]+f1[i]
              else:child_counts[clas]=f1[i]
            if(ids[1] not in done_conn):
              if(clas in connected_counts):connected_counts[clas] = connected_counts[clas]+f1[i]
              else:connected_counts[clas]=f1[i]
          done_chil[ids[0]], done_conn[ids[0]] = 1,1
  fh.close()
  rep_total = str(sum(total))
  out= "#ID\tClass\tnumber BCRs in class\tnumber of parents\tnumber of children\tnumber connected\ttotal in repertoire\n"
  for clas in totals:
    if(clas in parent_counts):par = parent_counts[clas]
    else:par = 0
    if(clas in child_counts):chil = child_counts[clas]
    else:chil = 0
    if(clas in connected_counts):conn = connected_counts[clas]
    else:conn = 0
    out=out+sample+"\t"+clas+"\t"+str(totals[clas])+"\t"+str(par)+"\t"+str(chil)+"\t"+str(conn)+"\t"+rep_total+"\n"
  fh=open(parental_classification_file_directed,"w")
  fh.write(out)
  fh.close()
  return()

def Get_class_usage(freq,classes):
  nz = [i for i in range(len(freq)) if freq[i]!=0]
  class_usage={}
  for i in range(0,len(nz)):
    clas = classes[i].split("*")[0]
    if(clas in class_usage):class_usage[clas] = class_usage[clas]+freq[i]
    else:class_usage[clas]=freq[i]
  return(class_usage)

def Mutational_region_isotypes(cluster_file, annot_file2, seq_file, Mutational_region_isotype_file):
  seqs,freqs,alias = Get_sequences(seq_file)
  for id in alias:
    chains= alias[id].split("|")[1].split("_")
    break
  chains_all = []
  for i in range(0,len(chains)):
    if(chains[i].split("*")[0] not in chains_all):chains_all = chains_all+[chains[i].split("*")[0]]
  chains_all.sort()
  indexing = [0]*len(chains)
  for i in range(0,len(chains)):
    indexing[i] =chains_all.index(chains[i].split("*")[0])
  heading_file = "Headings_mutations_2.txt"
  heading_use,columns,heading_out = [],[],[]
  fh=open(Mutational_region_isotype_file,"w")
  fh.close()
  fh=open(heading_file,"r")
  for l in fh:
    l=l.strip()
    heading_use.append(l)
    name = l.replace(" Nb of mutations","_N_mutations").replace("-IMGT","")
    #if(l.count("nonsilent")!=0):name = l.split('-')[0]+"-Nonsilent"
    #else:name = l.split('-')[0]+"-Silent"
    heading_out.append(name)
  fh.close()
  fh=open(annot_file2, "r")
  lens = len(heading_use)
  out,ind='#ID\t'+'\t'.join(heading_out)+"\t"+"\t".join(chains_all)+"\n",0
  for l in fh:
    l=l.strip().split("\t")
    if(l[0]=='Sequence number'):
      for i in heading_use:
        columns = columns +[l.index(i)]
    else:
      if(l[2]=="productive" and len(l)>max(columns)):
        muts = [0]*lens
        for i in range(0,lens):
          m = int(l[columns[i]].split()[0])
          muts[i] = m
        f = freqs[l[1].split("__")[0]]
        f1 = [0]*len(chains_all)
        for i in range(0,len(chains)):
          f1[indexing[i]] = f1[indexing[i]]+f[i]
        out=out+l[1].split("__")[0]+"\t"+"\t".join(map(str, muts))+"\t"+"\t".join(map(str, f1))+"\n"
        ind =ind+1
        if(ind>100):
          Write_out(out, Mutational_region_isotype_file)
          out, ind = '',0
  fh.close()
  Write_out(out, Mutational_region_isotype_file)
  out, ind = '',0
  return()

def Get_VJ_gene_usage_subtype_specific(annot_file_IMGT,annot_file, seq_file,subtybe_VJ_gene_usage_file,sample):
  CDR3s,annots = Get_annotations(annot_file_IMGT)
  seqs,freqs,alias = Get_sequences(seq_file)
  #annots = Get_annotation(annot_file)
  for id in alias:
    chains= alias[id].split("|")[1].split("_")
    break
  for i in range(0,len(chains)):
    chains[i] = chains[i].split("*")[0]
  VJ_usage,VJ_usage_reads = {},{}
  total1 = 0
  for id in annots:
    if(id in seqs):
      f = freqs[id]
      nz = [i for i in range(len(f)) if f[i]!=0]
      total1 = total1+sum(f)
      for i in nz:
        c = chains[i]
        if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
        elif(c in ["IGHG1","IGHG2"]):c= "IGHG1/2"
        v,j =  annots[id][3].split("\t")[0], annots[id][3].split("\t")[1]
        typ = c+"\t"+v+"\t"+j
        if(typ in VJ_usage):
          VJ_usage[typ] = VJ_usage[typ]+1
          VJ_usage_reads[typ] = VJ_usage_reads[typ]+f[i]
        else:
          VJ_usage[typ]=1
          VJ_usage_reads[typ]=f[i]
        if(c in ["IGHD","IGHM"]):
          c= "IGHD/M"
          muts = annots[id][0]
          if(muts<=3):c = c+"_Unmutated"
          else:c = c+"_Mutated"
          typ = c+"\t"+v+"\t"+j
          if(typ in VJ_usage):
            VJ_usage[typ] = VJ_usage[typ]+1
            VJ_usage_reads[typ] = VJ_usage_reads[typ]+f[i]
          else:
            VJ_usage[typ]=1
            VJ_usage_reads[typ]=f[i]
  out="#Sample\tIsotype\tV\tJ\tfrequency_uniq\tfrequency_reads\n"
  total2 = 0
  for typ in VJ_usage:
    out=out+sample+"\t"+typ+"\t"+str(VJ_usage[typ])+"\t"+str(VJ_usage_reads[typ])+"\n"
    total2 = total2+VJ_usage[typ]
  fh=open(subtybe_VJ_gene_usage_file,"w")
  fh.write(out)
  fh.close()
  print(total1, total2)
  return()

def Get_maximum_sequence_clusters(clusters, seqs, maximum_cluster_file,freqs,alias,maximum_cluster_file_filtered):
  threshold = 2 ######## filtering threshold number of reads
  for f in [maximum_cluster_file, maximum_cluster_file_filtered]:
    fh=open(f, "w")
    fh.close()
  for id in freqs:
    break
  l = len(freqs[id])
  totals = []
  for c in clusters:
    if(len(clusters[c])>10):
      tot = [0]*l
      for id in clusters[c]:
        tot = list(map(add,tot, freqs[id]))
      totals.append([c,sum(tot)])
  totals.sort(key=lambda x: x[1],reverse = True)
  out,ind,out1 = '',0,''
  c_max = totals[0][0]
  tot_to_align = 0
  for id in clusters[c_max]:
    out=out+">"+alias[id]+"\n"+seqs[id]+"\n"
    if(sum(freqs[id])>=threshold):
      out1=out1+">"+alias[id]+"\n"+seqs[id]+"\n"
      tot_to_align = tot_to_align+1
    ind = ind+1
    if(ind>500):
      Write_out(out, maximum_cluster_file)
      Write_out(out1, maximum_cluster_file_filtered)
      out,ind,out1 = '',0,''
  print(totals)
  Write_out(out, maximum_cluster_file)
  Write_out(out1, maximum_cluster_file_filtered)
  insert = ''
  if(tot_to_align > 5000):insert = '--parttree'
  command0 = "/software/pubseq/bin/mafft-6.857/bin/mafft --retree 2 "+insert+" "+maximum_cluster_file_filtered+" > "+maximum_cluster_file_filtered+".aligned"
  os.system(command0)
  return()

def Get_max_cluster_mutational_distances(maximum_cluster_file, freqs, alias, distance_from_central_file, sample,annots,seqs):
  fh=open(distance_from_central_file,"w")
  fh.close()
  ids,max_node,max_freq = {},'',[0]
  fh = open(maximum_cluster_file,"r")
  for header,seq in fasta_iterator(fh):
    id = header.split("__")[0]
    ids[id] = freqs[id]
    if(sum(freqs[id])>sum(max_freq)):max_node,max_freq = id, freqs[id]
  fh.close()
  out="#ID\tsequence_ID\tmutations from central\tchain\tfrequency\n"
  ind = 0
  max_seq = seqs[max_node]
  classes = alias[max_node].split("|")[1].split("_")
  for i in range(0,len(classes)):
    classes[i] = classes[i].split("*")[0]
  for id in ids:
    s1 = seqs[id]
    a1,a2 = Do_align(max_seq, s1)
    diffs = [i for i in range(len(a1)) if a1[i]!=a2[i]]
    mm = 0
    for i in diffs:
      if(a1[i] !=None and a2[i] !=None):
        mm = mm+1
    freq = freqs[id]
    nz = [i for i in range(len(freq)) if freq[i]!=0]
    class_freq = {}
    for i in nz:
      if(classes[i] in class_freq):class_freq[classes[i]] = class_freq[classes[i]]+freq[i]
      else:class_freq[classes[i]]=freq[i]
    for c in class_freq:
      out=out+sample+"\t"+id+"\t"+str(mm)+"\t"+c+"\t"+str(class_freq[c])+"\n"
      ind = ind+1
      if(ind>100):
        Write_out(out, distance_from_central_file)
        out, ind = '',0
  Write_out(out, distance_from_central_file)
  return()

def Get_maximum_cluster_information(cluster_file, seq_file, sample, edge_file, annot_file, maximum_cluster_file ,distance_from_central_file,maximum_cluster_file_filtered):
    clusters, inv_clusters,gl_mutations = Get_cluster_occupancy(cluster_file)
    annots = Get_annotation(annot_file)
    seqs,freqs,alias = Get_sequences(seq_file)
    Get_maximum_sequence_clusters(clusters, seqs, maximum_cluster_file,freqs,alias,maximum_cluster_file_filtered)
    Get_max_cluster_mutational_distances(maximum_cluster_file, freqs, alias, distance_from_central_file, sample,annots,seqs)
    return()

def Clonality_per_sample(seq_file,sample,clonality_file,cluster_file,reverse_primer_group):
  seqs,freqs,alias = Get_sequences(seq_file)
  cluster_total = {}
  fh=open(cluster_file, "r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cluster, id = l[1],l[2].split("__")[0]
      if(id in freqs):
        f = sum(freqs[id])
        if(cluster in cluster_total):cluster_total[cluster] = cluster_total[cluster]+f
        else:cluster_total[cluster] = f
  fh.close()
  sizes = {}
  for c in cluster_total:
    siz = cluster_total[c]
    if(siz in sizes):sizes[siz] = sizes[siz]+1
    else:sizes[siz]=1
  out="#sample\tsize\tfrequency\n"
  for siz in sizes:
    out=out+sample+"\t"+str(siz)+"\t"+str(sizes[siz])+"\n"
  fh=open(clonality_file,"w")
  fh.write(out)
  fh.close()
  return()

def Classify_sequences_into_developmental_stages_per_cluster(sample, annot_file_IMGT, cluster_file,per_cluster_developmental_classification_file,reverse_primer_group):
  clusters, inv_clusters = Get_cluster_occupancy(cluster_file)
  CDR3s,mutations = Get_annotations(annot_file_IMGT)
  #annots = Get_annotation(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  total = 0 
  for id in freqs:
    if(id in inv_clusters):
      l=len(freqs[id])
      chains= alias[id].split("|")[1].split("_")
      total = total+sum(freqs[id])
  print(total)
  for i in range(0,len(chains)):
    chains[i] = chains[i].split("*")[0]
    if(reverse_primer_group !="ISOTYPER"):
      if(chains[i] in ["IGHA1","IGHA2"]):chains[i] = "IGHA1/2"
      elif(chains[i] in ["IGHG1","IGHG2"]):chains[i] = "IGHG1/2"
  cluster_freqs,chains = {}, np.array(chains)
  fh=open(per_cluster_developmental_classification_file,"w")
  fh.close()
  out,ind = '#ID\tsequence\tclassifiation\tall_classes\tV\tJ\tmutation\tCDR3\tnode size\tnode %\tcluster ID\tcluster size\n',0
  print(len(alias))
  for c in clusters: 
    if(len(clusters[c])>0):
      freq,mm,vj = [0]*l,[],{}
      for id in clusters[c]:
        freq = list(map(add, freq, freqs[id]))
        if(id in mutations):
          mm = mm+[mutations[id][0]+mutations[id][1]]
          vj_class = mutations[id][3]
          if(vj_class in vj):vj[vj_class] = vj[vj_class]+sum(freqs[id])
          else:vj[vj_class]=sum(freqs[id])
      max_vj_f,max_vj = 0,'-\t-'
      for i in vj:
        if(max_vj_f<vj[i]):max_vj_f,max_vj = vj[i], i
      if(len(mm)==0):mm=-1
      else:mm=mean(mm)
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      chain = ",".join(np.sort(np.unique(chains[nz])))
      chains_used = np.sort(np.unique(chains[nz]))
      classifiations=[chain]
      if(chain in ["IGHD,IGHM","IGHD","IGHM"] and mm <=4):classifiations.append("IGHD,IGHM_unmutated")
      elif(chain in ["IGHD,IGHM","IGHD","IGHM"] and mm > 4):classifiations.append("IGHD,IGHM_mutated")
      #if(mm > 0):classifiations.append("mutated")
      MD = 0
      if("IGHD" in chains_used):MD = MD+1
      if("IGHM" in chains_used):MD = MD+1
      if(len(chains_used)>MD):classifiations.append("Class_switched")
      classifiations.append(str(len(chains_used))+"_isotype(s)")
      if(len(chains_used)>1):
        for i in range(0,len(chains_used)):
          classifiations.append(chains_used[i])
      for id in clusters[c]:
        out=out+sample+"\t"+id+"\t"+"|".join(classifiations)+"\t"+chain+"\t"+max_vj+"\t"+str(mm)+"\t"+"CDR3"+"\t"+str(sum(freq))+"\t"+str(sum(freqs[id])*100.0/total)+"\t"+c+"\t"+str(len(clusters[c]))+"\n"
        ind = ind+1
        if(ind>500):
          Write_out(out,per_cluster_developmental_classification_file)
          out, ind = '',0
  Write_out(out,per_cluster_developmental_classification_file)
  return()

def Get_nearest_neighbour_distributions(sample, annot_file, edge_file, nearest_neighbour_file):
  #CDR3s,annots = Get_annotations(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  for id in freqs:
    chains= alias[id].split("|")[1].split("_")
    break
  fh = open(edge_file, "r")
  nneighbours = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(int(l[2])==1):
        v1,v2 = l[0].split("__")[0], l[1].split("__")[0]
        vs = [v1,v2]
        classes,len_classes =[],[]
        for i in v1,v2:
          x,freq = [],freqs[i]
          nz = [i for i in range(len(freq)) if freq[i]!=0]
          for i in nz:
            c = chains[i].split("*")[0]
            if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
            elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
            elif(c in ["IGHD","IGHM"]):c = "IGHD/M"
            if(c not in x):x.append(c)
          classes.append(x)
          len_classes .append(len(x))
        if(max(len_classes)==1):
          n = classes[0][0]+"\t"+classes[1][0]
          if(n in nneighbours):nneighbours[n] = nneighbours[n]+1
          else:nneighbours[n] = 1
          n = classes[1][0]+"\t"+classes[0][0]
          if(n in nneighbours):nneighbours[n] = nneighbours[n]+1
          else:nneighbours[n] = 1
  fh.close()
  out = "#sample\tisotype_from\tisotype_to\tfrequency_of_nearest_neighbours\n"
  for n in nneighbours:
    out = out+"\t".join(map(str, [sample,n,nneighbours[n]]))+"\n"
  fh=open(nearest_neighbour_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_nn_CDR3_regions(annot_file3):
  nn_CDR3s = {}
  fh=open(annot_file3,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split("\t")
      if(l[2]!="unproductive" and l[2]!='No results' and len(l)>14):
        id, CDR3 = l[1].split("__")[0], l[14].upper()
        nn_CDR3s[id] = CDR3
  fh.close()
  return(nn_CDR3s)

def Get_seq_type(per_cluster_developmental_classification_file):
  seq_types = {}
  fh=open(per_cluster_developmental_classification_file,"r")
  clusters, clusters_all = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split("\t")
      clust = l[10]
      id, clas = l[1],l[2]#.split("|")
      c = "ALL"
      if(clas.count("Class_switched")!=0):c = "Class_switched"
      elif(clas.count("IGHD,IGHM_unmutated")!=0): c = "IGHD/M_unmutated"
      elif(clas.count("IGHD,IGHM_mutated")!=0): c = "IGHD/M_mutated"
      if(clust not in clusters):clusters[clust] = [id, c] 
      clusters_all[id] = c
      #if(c in clusters_all):clusters_all[c] = clusters_all[c]+[id]  
      #else:clusters_all[c]=[id] 
  fh.close()
  for c in clusters:
    seq_types[clusters[c][0]] = clusters[c][1]
  del clusters
  return(seq_types, clusters_all)

def Get_public_CDR3_percentages(subsample_file,annot_file, per_cluster_developmental_classification_file, unique_CDR3_regions_per_isotype_group_file,annot_file3,sample,public_CDR3_file):
  public_CDR3_file_all = "/lustre/scratch118/infgen/team146/rbr1/MISEQ/AUTOIMMUNITY_AAV_SLE_HEALTHY1/ORIENTATED_SEQUENCES/Public_BCRs1.txt"
  #Pre_filter_public_BCRs(subsample_file,annot_file, per_cluster_developmental_classification_file, unique_CDR3_regions_per_isotype_group_file,annot_file3,sample,public_CDR3_file,public_CDR3_file_all)
  Public_private_CDR3_lengths(public_CDR3_file_all, subsample_file, annot_file,public_CDR3_file)
  return()

def Public_private_CDR3_lengths(public_CDR3_file_all, subsample_file, annot_file,public_CDR3_file):
  fh=open(public_CDR3_file_all,"r")
  publics = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(int(l[1])>=5):
        publics[l[0]]=int(l[1])
  fh.close()
  CDR3s,annots = Get_annotations(annot_file)
  public_CDR3s, private_CDR3s = {},{}
  for id in annots:
    CDR3_aa = annots[id][2].split(":")[0]
    if(len(CDR3_aa)>4):
      CDR3_aa = CDR3_aa[1:len(CDR3_aa)]
      length = len(CDR3_aa)
      if(CDR3_aa in publics):
        if(length in public_CDR3s):public_CDR3s[length] = public_CDR3s[length]+[CDR3_aa]
        else:public_CDR3s[length]=[CDR3_aa]
      else:
        if(length in private_CDR3s):private_CDR3s[length] = private_CDR3s[length]+[CDR3_aa]
        else:private_CDR3s[length]=[CDR3_aa]
  out = "#sample\tCDR3_length\tpublic CDR3s\tprivate_CDR3s\tuniq_public CDR3s\tuniq_private_CDR3s\n"
  for l in range(4,60):
    f_public,f_private,uf_public, uf_private = 0,0,0,0
    if(l in public_CDR3s):f_public,uf_public = len(public_CDR3s[l]), len(set(public_CDR3s[l]))
    if(l in private_CDR3s):f_private,uf_private = len(private_CDR3s[l]), len(set(private_CDR3s[l]))
    if(f_public+ f_private ):
      out=out+"\t".join(map(str, [sample,l,f_public,f_private,uf_public, uf_private ]))+"\n"
  fh=open(public_CDR3_file.replace("CDR3_file_","CDR3_PRIVATE_proportion_length"),"w")
  fh.write(out)
  fh.close()
  return()

def Pre_filter_public_BCRs(subsample_file,annot_file, per_cluster_developmental_classification_file, unique_CDR3_regions_per_isotype_group_file,annot_file3,sample,public_CDR3_file,public_CDR3_file_all):
  fh=open(public_CDR3_file_all,"r")
  publics = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      publics[l[0]]=int(l[1])
  fh.close()
  CDR3s,annots = Get_annotations(annot_file)
  seq_types,clusters_all = Get_seq_type(per_cluster_developmental_classification_file)
  public_counts = {}
  public_sequence = {}
  for id in annots:
    CDR3_aa = annots[id][2].split(":")[0]
    if(len(CDR3_aa)>4):
      CDR3_aa = CDR3_aa[1:len(CDR3_aa)]
      #print CDR3_aa
      if(CDR3_aa in publics):
        t = "PUBLIC\t"+str(publics[CDR3_aa])
        t1 = CDR3_aa+"\t"+str(publics[CDR3_aa])
        if(t1 not in public_sequence):public_sequence[t1] = [id]
        else: public_sequence[t1] =public_sequence[t1] +[id]
      else:t = "PRIVATE\t0"
      seq_type = clusters_all[id]
      name = seq_type+"\t"+t
      if(name in public_counts):public_counts[name] = public_counts[name]+1
      else:public_counts[name] = 1
      if(t.count("PUBLIC")==1):
        t = "PUBLIC\t>1"
        name = seq_type+"\t"+t
        if(name in public_counts):public_counts[name] = public_counts[name]+1
        else:public_counts[name] = 1
  out="sample\tCDR3\tcount\tids\n"
  for t1 in public_sequence:
    out = out+"\t".join(map(str, [sample, t1]+public_sequence[t1]))+"\n"
  fh=open(public_CDR3_file.replace("CDR3_file_","CDR3_IDS_"),"w")
  fh.write(out)
  fh.close()
  out="#sample\tisotype\tstatus\tshared_between_n_samples\tnumber of sequences\n"
  for n in public_counts:
    out=out+"\t".join(map(str, [sample, n,public_counts[n]]))+"\n"
  fh=open(public_CDR3_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_unique_CDR3_regions_per_isotype_group(subsample_file,annot_file, per_cluster_developmental_classification_file, unique_CDR3_regions_per_isotype_group_file,annot_file3,sample,CDR3_length_v_genes_file):
  nn_CDR3s=Get_nn_CDR3_regions(annot_file3)
  CDR3s,annots = Get_annotations(annot_file)
  seq_types,clusters_all = Get_seq_type(per_cluster_developmental_classification_file)
  types,CDR3s_refined = {},{}
  for id in seq_types:
    if(id in annots):
      CDR3_aa = annots[id][2].split(":")[0]
      if(id in nn_CDR3s):
        CDR3_nn = nn_CDR3s[id]
        CDR3_aa = CDR3_aa[1:len(CDR3_aa)-1]
        if(len(CDR3_nn)>9 and len(CDR3_aa)>3):
          CDR3s_refined[id] = [CDR3_nn,CDR3_aa]
          if(seq_types[id] in types):types[seq_types[id]] = types[seq_types[id]]+[id]#, CDR3_nn, CDR3_aa]]
          else:types[seq_types[id]]=[id]#, CDR3_nn, CDR3_aa]]
          t = "ALL"
          if(t in  types):types[t] = types[t]+[id]
          else:types[t]=[id]
  #fh=open(subsample_file,"r")
  n_repeats = 1000
  depths = [10,20,50,100,200,500,1000,2000,5000]
  samples = Tree()
  out = "sample\tsequence_group\tn_clusters_total\tsample_depth\tmean_unique_CDR3aa\tmean_unique_CDR3nn\tmean_ratio_CDR3aa/CDR3nn\n"
  for t in types:
    for depth in depths:
      if(len(types[t])>=depth):
        meanaa, meannn,meanratio = [],[],[]
        clust_sub =types[t]
        for r in range(0,n_repeats):
          draw = np.random.choice(clust_sub, depth)#, p=probability_distribution)
          Caa, Cnn = {},{}
          for id in draw:
            Caa[CDR3s_refined[id][1]], Cnn[CDR3s_refined[id][0]] = 1,1
          meanaa, meannn,meanratio = meanaa+[len(Caa)], meannn+[len(Cnn)],meanratio+[len(Caa)*1.0/len(Cnn)]
        out=out+ "\t".join(map(str,  [sample, t,len(types[t]), depth, mean(meanaa), mean(meannn), mean(meanratio)]))+"\n"
  print(unique_CDR3_regions_per_isotype_group_file)
  fh=open(unique_CDR3_regions_per_isotype_group_file,"w")
  fh.write(out)
  fh.close()
  return()

def sdfg():
  #for l in fh:
  #  if(l[0]!="#"):
  #    l=l.strip().split("\t")
  #    subs, id = l[0], l[2].split("__")[0]
  #    samples[id][subs].value = 1
  fh.close()
  #seqs,freqs,alias = Get_sequences(seq_file)
  CDR3s_all = Tree()
  mean_CDR3_lengths = {}
  for id in annots:
    CDR3_aa = annots[id][2].split(":")[0]
    if(id in nn_CDR3s and id in seq_types):
      v = annots[id][2].split(":")[1]
      CDR3_nn = nn_CDR3s[id]
      CDR3_aa = CDR3_aa[1:len(CDR3_aa)-1]
      if(len(CDR3_nn)>9 and len(CDR3_aa)>3):
        t = seq_types[id]+["ALL"]
        for t1 in t:
          pas = 0
          if(t1 in ["3_isotype(s)","4_isotype(s)","5_isotype(s)","6_isotype(s)"]):
            t1 = ">=3_isotype(s)"
            pas = 1
          elif(t1.count("_")!=0):pas = 1
          elif(t1.count(",")==0):pas = 1
          if(pas==1):
            t2 = t1+"\t"+v
            if(t2 in mean_CDR3_lengths):mean_CDR3_lengths[t2] = mean_CDR3_lengths[t2]+[len(CDR3_aa)]
            else:mean_CDR3_lengths[t2] = [len(CDR3_aa)]
            for subs in samples[id]:
              CDR3s_all[subs][t1]["AA"][CDR3_aa].value = 1
              CDR3s_all[subs][t1]["NN"][CDR3_nn].value = 1
  out ="#sample\tsequence_group\tV_gene\tn_seqs\tmean_CDR3_length\n"
  for t in mean_CDR3_lengths:
    out=out+ "\t".join(map(str, [sample,t, len(mean_CDR3_lengths[t]), mean(mean_CDR3_lengths[t])]))+"\n"
  fh=open(CDR3_length_v_genes_file,"w")
  fh.write(out)
  fh.close()
  out="#sample\tsubsample\tsequence_group\tn_aa_CDR3s\tn_nn_CDR3s\tratio_NN/AA\n"
  for subs in CDR3s_all:
    for t in CDR3s_all[subs]:
      out=out+ "\t".join(map(str, [sample,subs, t, len(CDR3s_all[subs][t]["AA"]), len(CDR3s_all[subs][t]["NN"]), len(CDR3s_all[subs][t]["NN"])*1.0/len(CDR3s_all[subs][t]["AA"])]))+"\n"
  fh=open(unique_CDR3_regions_per_isotype_group_file,"w")
  fh.write(out)
  fh.close()
  return()

def Classify_sequences_into_developmental_stages(sample, annot_file, cluster_file, developmental_classification_file,reverse_primer_group):
  CDR3s,annots = Get_annotations(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  total = 0
  for id in freqs:
    chains= alias[id].split("|")[1].split("_")
    total = total+sum(freqs[id])
  for i in range(0,len(chains)):
    chains[i] = chains[i].split("*")[0]
  chains = np.array(chains)
  out,ind = '#ID\tsequence\tclassifiation\tall_classes\tV\tJ\tmutation\tCDR3\tnode size\tnode %\n',0
  fh=open(developmental_classification_file, "w")
  fh.close()
  for id in freqs:
    mm,info = -1,"-\t-\t-1"
    if(id in annots):
      #info = "\t".join(map(str,annots[id]))
      info = annots[id][3]+"\t"+str(annots[id][0])+"\t"+annots[id][2].split(":")[0]
      mm = annots[id][0]
    freq = freqs[id]
    nz = [i for i in range(len(freq)) if freq[i]!=0]
    chain = ",".join(np.sort(np.unique(chains[nz])))
    chains_used = np.sort(np.unique(chains[nz]))
    classifiations=[chain]
    if(chain in ["IGHD,IGHM","IGHD","IGHM"] and mm == 0):classifiations.append("IGHD,IGHM_unmutated")
    elif(chain in ["IGHD,IGHM","IGHD","IGHM"] and mm > 0):classifiations.append("IGHD,IGHM_mutated")
    #if(mm > 0):classifiations.append("mutated")
    MD = 0
    if("IGHD" in chains_used):MD = MD+1
    if("IGHM" in chains_used):MD = MD+1
    if(len(chains_used)>MD):classifiations.append("Class_switched")
    if(len(chains_used)>1):
      for i in range(0,len(chains_used)):
        classifiations.append(chains_used[i])
    out=out+sample+"\t"+id+"\t"+"|".join(classifiations)+"\t"+chain+"\t"+info+"\t"+str(sum(freq))+"\t"+str(sum(freq)*100.0/total)+"\n"
    ind = ind+1
    if(ind>500):
      Write_out(out, developmental_classification_file)
      out, ind = '',0
  Write_out(out, developmental_classification_file)
  return()


def Get_pre_post_class_switching(cluster_file, sample, annot_file, pre_post_class_switching_file,seq_file):
  clusters, inv_clusters,gl_mutations = Get_cluster_occupancy(cluster_file)
  annots = Get_annotation(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  for id in alias:
    chains= alias[id].split("|")[1].split("_")
    break
  for i in range(0,len(chains)):
    chains[i] = chains[i].split("*")[0]
  chains = np.array(chains)
  fh=open(pre_post_class_switching_file,"w")
  fh.close()
  out, ind = '#ID\tcluster ID\tcluster_unique_seqs\ttotal_reads\tunique sequences of chain type\tchain(s) used\tmin mm\tmax mm\tmm distribution\n',0
  for c in clusters: 
    if(len(clusters[c])>0):
      chain_mutations = {}
      n_sequences,n_sequences_per_chain,n_reads = len(clusters[c]),{},[]
      for id in clusters[c]:
        if(id in annots):
          freq = freqs[id]
          if(len(n_reads)==0):n_reads = freq
          else:n_reads = list(map(add,n_reads, freq))
          nz = [i for i in range(len(freq)) if freq[i]!=0]
          chain = ",".join(np.sort(np.unique(chains[nz])))
          if(chain in chain_mutations):
            chain_mutations[chain] = chain_mutations[chain]+[annots[id][2]]
            n_sequences_per_chain[chain] = n_sequences_per_chain[chain]+sum(freq)
          else:
            chain_mutations[chain]=[annots[id][2]]
            n_sequences_per_chain[chain]=sum(freq)
      for chain in chain_mutations:
        out=out+sample+"\t"+c+"\t"+str(n_sequences)+"\t"+str(sum(n_reads))+"\t"+str(n_sequences_per_chain[chain])+"\t"+chain+"\t"+str(min(chain_mutations[chain]))+"\t"+str(max(chain_mutations[chain]))+"\t"+",".join(map(str, sort(chain_mutations[chain])))+"\n"
        ind = ind+1
        if(ind>300):
          Write_out(out, pre_post_class_switching_file)
          out, ind = '',0
  Write_out(out, pre_post_class_switching_file)
  return()

def Get_cluster_mutation_sharing_stats(cluster_file, annot_file, seq_file, cluster_mutation_sharing_probability,sample):
  #clusters, inv_clusters,gl_mutations = Get_cluster_occupancy(merged_cluster_file)
  clusters, inv_clusters = Get_cluster_occupancy(cluster_file)
  annots = Get_annotation(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  for id in alias:
    chains= alias[id].split("|")[1].split("_")
    break
  chains_all = []
  for i in range(0,len(chains)):
    #chains[i] = chains[i].split("*")[0]
    if(chains[i].split("*")[0] not in chains_all):chains_all = chains_all+[chains[i].split("*")[0]]
  chains_all.sort()
  print(chains_all)
  fh=open(cluster_mutation_sharing_probability,"w")
  fh.close()
  indexing = [0]*len(chains)
  for i in range(0,len(chains)):
    indexing[i] =chains_all.index(chains[i].split("*")[0])
  out, ind = '#ID\tmean mutations\tnumber of chains shared\tchain combination\tnumber sequences\tnumber reads\n',0
  mutation_isotype = {}
  number_reads = {}
  for c in clusters: 
    f = [0]*len(chains)
    muts = []
    for id in clusters[c]:
      f = list(map(add, f, freqs[id]))
      if(id in annots):
        muts = muts+[annots[id][2]]*sum(freqs[id])
    f1 = [0]*len(chains_all)
    for i in range(0,len(chains)):
      f1[indexing[i]] = f1[indexing[i]]+f[i]
    chains_used = []
    for i in range(0,len(chains_all)):
      if(f1[i]!=0):
        chains_used.append(chains_all[i])
    mean_mut = "NA"
    shared = len(f1)-f1.count(0)
    if(len(muts)>0):mean_mut = str(mean(muts))
    name = mean_mut+"\t"+str(shared)+"\t"+",".join(map(str, chains_used))
    if(name in mutation_isotype):
      mutation_isotype[name] = mutation_isotype[name]+1
      number_reads[name] =number_reads[name]+[sum(f1)]
    else:
      mutation_isotype[name]=1
      number_reads[name]=[sum(f1)]
  for name in mutation_isotype:
    out=out+sample+"\t"+name+"\t"+str(mutation_isotype[name])+"\t"+",".join(map(str,number_reads[name]))+"\n"
    #out=out+sample+"\t"+c+"\t"+str(len(clusters[c]))+"\t"+mean_mut+"\t"+str(shared)+"\t"+"\t".join(map(str, chains_used))+"\n"
    ind = ind+1
    if(ind>100):
      Write_out(out, cluster_mutation_sharing_probability)
      out, ind = '',0
  Write_out(out, cluster_mutation_sharing_probability)
  return()

def Classify_transitions(seq_file,edge_file, transition_classification_file_directed,cluster_file, children_classification_file_directed,parental_classification_file_directed,annot_file,sample):
  #print transition_classification_file_directed
  #print children_classification_file_directed
  #print parental_classification_file_directed
  annots = Get_annotation(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  Assess_transitional_frequencies(merged_cluster_nearest_neighbours,clusters, inv_clusters,gl_mutations,transition_classification_file_directed)
  Assess_parental_frequencies(annots,parental_classification_file_directed,edge_file,freqs,alias,sample,cluster_file)
  return()

def Reclassify_sequences(f, reclassify, classification,classes):
  f1 = [0]*len(classes)
  non_zero = [i for i in range(len(f)) if f[i]!=0]
  for i in non_zero:
    f1[reclassify[classification[i]]] = f1[reclassify[classification[i]]] + f[i]
  return(f1)

def Get_classes(file):
  fh=open(file,"r")
  classes,classes_inv,clusters = Tree(),{},Tree()
  definitions_use = ['Class_switched','IGHD,IGHM_mutated','IGHD,IGHM_unmutated']
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      classes_all = l[2].split("|")
      classes_use = []
      for c in classes_all:
        if(c in definitions_use):
          #classes[c][l[1]].value = 1
          classes_use.append(c)
          #clusters[l[10]][c].value = 1
      classes_inv[l[1]] = classes_use
  fh.close()
  return(classes,classes_inv,definitions_use)

def Get_class_overlap_subclass(seq_file,class_overlap_subtype_file,sample,developmental_classification_file,cluster_file,reverse_primer_group):
  classes,classes_inv,definitions_use = Get_classes(developmental_classification_file)
  #clusters, inv_clusters,gl_mutations = Get_cluster_occupancy(cluster_file)
  clusters, inv_clusters= Get_cluster_occupancy(cluster_file)
  cluster_types, cluster_sizes = Tree(),{}
  count_classes = {}
  fh=open(cluster_file, "r")
  total = 0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      clust, id, freq  = l[1],l[2].split("__")[0], int(l[3])
      classes = classes_inv[id]
      for clas in classes:
        cluster_types[clust][clas].value = 1
        if(clust in cluster_sizes):cluster_sizes[clust]=cluster_sizes[clust]+freq
        else:cluster_sizes[clust]=freq
        total = total + freq
  fh.close()
  overlap,cluster_counts={},{}
  for c in cluster_types:
    if(len(cluster_types[c])>1):
      class_use = []
      for clas in cluster_types[c]:
        class_use.append(clas)
      class_use.sort()
      for i in range(0,len(class_use)):
        for j in range(i,len(class_use)):
          if(i<j):
            name = class_use[i]+"\t"+class_use[j]
            if(name in overlap):
              overlap[name] = overlap[name] + cluster_sizes[c]
              cluster_counts[name] = cluster_counts[name]+1
            else:
              overlap[name] = cluster_sizes[c]
              cluster_counts[name]=1
  out="#ID\tclass1\tclass2\tread_overlap\tread_total\tcluster_overlap\tcluster_total\n"
  for n in overlap:
    out=out+sample+"\t"+n+"\t"+str(overlap[n])+"\t"+str(total)+"\t"+str(cluster_counts[n])+"\t"+str(len(cluster_sizes))+"\n"
  fh=open(class_overlap_subtype_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_class_overlap(seq_file, class_overlap_file,sample):
  seqs,freqs,alias = Get_sequences(seq_file)
  totals,overlaps = {},{}
  for id in seqs:
    classification = alias[id].split("|")[1].split("_")
    break
  reclassify,classes = {},[]
  for c in classification:
    c1 = c.split("*")[0]
    if(c1 not in classes):classes.append(c1)
    reclassify[c] = classes.index(c1)
  for id in seqs:
    f = freqs[id]
    classification = alias[id].split("|")[1].split("_")
    f1 = Reclassify_sequences(f, reclassify, classification,classes)
    non_zero = [i for i in range(len(f1)) if f1[i]!=0]
    for i in range(0,len(non_zero)):
      chain = classes[non_zero[i]]
      if(chain in totals):totals[chain] = totals[chain]+1
      else:totals[chain] = 1
      if(len(non_zero)>1):
        for j in range(i,len(non_zero)):
          if(i<j):
            over = classes[non_zero[i]]+"\t"+classes[non_zero[j]]
            if(over in overlaps):overlaps[over] = overlaps[over]+1
            else:overlaps[over] =1
  out = "#Id\tfrom\tto\tfrequency\n"
  for i in totals:
    out=out+sample+"\t"+i+"\tALL\t"+str(totals[i])+"\n"
  for i in overlaps:
    out=out+sample+"\t"+i+"\t"+str(overlaps[i])+"\n"
  out = out+"ALL\tALL\t"+str(len(seqs))+"\n"
  fh=open(class_overlap_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_GL_distance_class_distributions(seq_file, annot_file,GL_distance_class_distribution_file):
  clusters, inv_clusters,gl_mutations = Get_cluster_occupancy(merged_cluster_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  fh=open(GL_distance_class_distribution_file,"w")
  fh.close()
  total_distribution,cluster_sizes,cluster_distribution = {},{},{}
  for c in clusters:
    for id in clusters[c]:
      f = freqs[id]
      classification = alias[id].split("|")[1].split("_")
      non_zero = [i for i in range(len(f)) if f[i]!=0]
      mut = gl_mutations[id]
      if(c in cluster_sizes):cluster_sizes[c] = cluster_sizes[c]+len(non_zero)
      else:cluster_sizes[c] = len(non_zero)
      for i in non_zero:
        typ = classification[i].split("*")[0]+"\t"+str(sum(mut))
        if(c+"\t"+typ in cluster_distribution):cluster_distribution[c+"\t"+typ] = cluster_distribution[c+"\t"+typ]+1
        else:cluster_distribution[c+"\t"+typ]=1
        if(typ in total_distribution):total_distribution[typ] = total_distribution[typ]+1
        else:total_distribution[typ] = 1
  out,ind = "#Cluster\tClass\tVJ_Mutations\tFrequency\n", 0
  for t in total_distribution:
    out = out+"ALL\t"+t+"\t"+str(total_distribution[t])+"\n"
    ind = ind+1
    if(ind>100):
      Write_out(out, GL_distance_class_distribution_file)
      out, ind = '',0
  Write_out(out, GL_distance_class_distribution_file)
  out, ind = '',0
  for t in cluster_distribution:
    c = t.split("\t")[0]
    if(cluster_sizes[c]>=10):
      out=out+t+"\t"+str(cluster_distribution[t])+"\n"
      print(t+"\t"+str(cluster_distribution[t]))
      ind = ind+1
      if(ind>100):
        Write_out(out, GL_distance_class_distribution_file)
        out, ind = '',0
  Write_out(out, GL_distance_class_distribution_file)
  out, ind = '',0
  return()

def IgD_IgM_ratio(annot_file, seq_file,IgD_IgM_ratio_file,cluster_file,sample):
  seqs,freqs,alias = Get_sequences(seq_file)
  annots = Get_annotation(annot_file)
  clusters, inv_clusters,gl_mutations = Get_cluster_occupancy(cluster_file)
  for id in seqs:
    classification = alias[id].split("|")[1].split("_")
    break
  col_types = []
  chain_types = ["IGHM","IGHD"]
  for i1 in range(0,len(chain_types)):
    ind = [i for i in range(len(classification)) if classification[i].count(chain_types[i1])!=0]
    col_types.append(ind)
  type_freq = {}
  for id in seqs:
    if(id in annots):
      freq = freqs[id]
      freq1 = [0]*len(chain_types)
      for i in range(0,len(col_types)):
        for j in col_types[i]:
          freq1[i] = freq1[i]+freq[j]
      if(sum(freq)==sum(freq1)):
      #if(sum(freq1)>=2):
        mut = annots[id][2]
        name = str(sum(freq1))+"\t"+str((freq1[0]-freq1[1])/(freq1[0]+freq1[1]))+"\t"+str(mut)
        if(name in type_freq):type_freq[name] =type_freq[name]+sum(freq1)
        else:type_freq[name] = sum(freq1)
  out="#ID\tTotal_reads\tIgD-IgM ratio\tN mutations\tnumber of sequences\n"
  for name in type_freq:
    out=out+sample+"\t"+name+"\t"+str(type_freq[name])+"\n"
  print(IgD_IgM_ratio_file)
  fh=open(IgD_IgM_ratio_file,"w")
  fh.write(out)
  fh.close()
  return()


def Cluster_size_connection_sharing_probability(annot_file, seq_file, merged_cluster_nearest_neighbours, cluster_size_connection_sharing_probability_file):
  seqs,freqs,alias = Get_sequences(seq_file)
  clusters, inv_clusters,gl_mutations = Get_cluster_occupancy1(merged_cluster_file)
  out, ind = '#Cluster_ID\tNumber_chains_shared\tmean_number_of_mutations\ttotal_frequency\n',0
  fh=open(cluster_size_connection_sharing_probability_file,"w")
  fh.close()
  for c in clusters:
    muts_cluster,freq = [],[]
    for id in clusters[c]:
      if(id in gl_mutations):
        if(len(freq)==0):freq = freqs[id]
        else:freq = list(map(add, freq, freqs[id]))
        muts_cluster = muts_cluster+[sum(gl_mutations[id])]
    if(len(muts_cluster)>0):
      number_shared = len([i for i in range(len(freq)) if freq[i]!=0])
      out=out+c+"\t"+str(number_shared)+"\t"+str(mean(muts_cluster))+"\t"+str(sum(freq))+"\n"
      ind = ind+1
      if(ind>300):
        Write_out(out, cluster_size_connection_sharing_probability_file)
        out,ind = '',0
  Write_out(out, cluster_size_connection_sharing_probability_file)
  return()

def Clone_size_of_mutated_and_unmutated(annot_file, seq_file,mean_node_size_mutated_unmutated,sample):
  annots = Get_annotation(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  sizes,all_sizes = {},[]
  for id in annots:
    mut = annots[id][2]
    f = freqs[id]
    classification = alias[id].split("|")[1].split("_")
    non_zero = [i for i in range(len(f)) if f[i]!=0]
    mut_stat = "Unmutated"
    if(mut!=0):mut_stat = "Mutated"
    for i in non_zero:
      name = classification[i]+"\t"+mut_stat
      if(name in sizes):sizes[name] = sizes[name]+[f[i]]
      else:sizes[name]=[f[i]]
      if(f[i] not in all_sizes):all_sizes = all_sizes+[f[i]]
  all_sizes = sort(all_sizes)
  out="#ID\tchain\tmutational_status\tfreqs\n"
  for name in sizes:
    out= out+sample+"\t"+name+"\t"+",".join(map(str,sizes[name]))+"\n"
  fh=open(mean_node_size_mutated_unmutated,"w")
  fh.write(out)
  fh.close()
  return()

def Get_residue_positivity():
  file = "/lustre/scratch118/infgen/team146/rbr1/REFERENCE_FILES/Amino_acid_side_chain_charges.txt"
  fh=open(file, "r")
  pka = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[1]=="Pos"):pka[l[0]]=1
      else:pka[l[0]]=-1
  fh.close()
  return(pka)

def Get_positivity(seq,pka):
  #m = 0
  #aas = list(set(seq))
  #for aa in pka:
  #  print aa, pka[aa]
  #  m = m + pka[aa]* seq.count(aa)
  m_acidic = seq.count("D")+seq.count("E")
  m_basic = seq.count("K")+seq.count("R")+seq.count("H")
  m = m_basic-m_acidic
  return(m)

def Theil_indices_fast(x):
  xi = np.array(x)
  n = len(xi)
  mu = np.mean(xi)
  T0 = (1.0/n)* sum((xi/mu)*np.log(xi/mu))
  T1 = (1.0/n)* sum(np.log(mu/xi))
  return(T0,T1)

def Theil_indices(points,vdf):
  x = np.array(vdf)
  p = np.array(points)
  xi = np.repeat(x, p)
  n = len(xi)
  mu = np.mean(xi)
  T0 = (1.0/n)* sum((xi/mu)*np.log(xi/mu))
  T1 = (1.0/n)* sum(np.log(mu/xi))
  return(T0,T1)

def Gini_index(cpoints,cvdf, vpoints,vvdf): 
  (vgini)=Get_Gini(vpoints,vvdf)
  (cgini)=Get_Gini(cpoints,cvdf)
  return(vgini, cgini)

def Get_Renyi(x):
  x1 = np.array(x)
  if(sum(x1)==0): print(x1)
  x1 = x1*1.0/sum(x1)
  r = -1*sum(x1*np.log(x1))
  return(r)

def Get_Gini_fast(x, w=None):# w = weight if required
    # The rest of the code requires numpy arrays.
    x = np.asarray(x)
    if w is not None:
        w = np.asarray(w)
        sorted_indices = np.argsort(x)
        sorted_x = x[sorted_indices]
        sorted_w = w[sorted_indices]
        # Force float dtype to avoid overflows
        cumw = np.cumsum(sorted_w, dtype=float)
        cumxw = np.cumsum(sorted_x * sorted_w, dtype=float)
        return (np.sum(cumxw[1:] * cumw[:-1] - cumxw[:-1] * cumw[1:]) / 
                (cumxw[-1] * cumw[-1]))
    else:
        sorted_x = np.sort(x)
        n = len(x)
        cumx = np.cumsum(sorted_x, dtype=float)
        # The above formula, with all weights equal to 1 simplifies to:
        return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n

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

def Identify_example_clusters(sample, annot_file, cluster_file,per_cluster_developmental_classification_file,example_cluster_annot_file,seq_file):
  #Get_clusters_specific(sample, annot_file, cluster_file,per_cluster_developmental_classification_file,example_cluster_annot_file,seq_file)
  Get_largest_cluster(sample, annot_file, cluster_file,per_cluster_developmental_classification_file,example_cluster_annot_file,seq_file)
  return()

def Get_largest_cluster(sample, annot_file, cluster_file,per_cluster_developmental_classification_file,example_cluster_annot_file,seq_file):
  clusters,cluster_inv,cluster_information = Tree(),{},{}
  fh=open(cluster_file,"r")
  n_isotypes = []
  for i in range(4,12):
    n_isotypes.append(str(i)+"_isotype(s)")
  n_vertices = 0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,cluster = l[2],l[1]
      clusters[cluster][id].value = 1
      n_vertices = n_vertices+1
  fh.close()
  classes = l[2].split("|")[1].split("_")
  clusters_count = []
  max_c, max_size = '',0
  for c in clusters:
    if(len(clusters[c])>200):
      clusters_count.append([c, len(clusters[c])])
  clusters_count.sort(key=lambda x: x[1], reverse = True)
  fh=open(example_cluster_annot_file,"w")
  fh.close()
  out = "#sample\tcluster_ID\tcluster_n_unique_reads\ttotal_reads\tisotype\tfrequency_uniq\tfrequency_total\n"
  ind = 0
  n_c = 0
  for c in clusters_count: 
    n_c = n_c+1
    if(n_c>10):break
    f_uniq,f = [0]*len(classes), [0]*len(classes)
    for id in clusters[c[0]]:
      count = list(map(int,id.split("__")[1].split("|")[0].split("_")))
      nz = [i for i in range(len(count)) if count[i]!=0]
      for i in nz:
        f_uniq[i],f[i] =f_uniq[i]+1,f[i]+count[i]
    nz = [i for i in range(len(f)) if f[i]!=0]
    for i in nz:
      out=out+"\t".join(map(str, [sample,c[0],len(clusters[c[0]]), n_vertices,classes[i],f_uniq[i],f[i]]))+"\n"
      ind= ind+1
      if(ind>100):
        Write_out(out, example_cluster_annot_file)
        out, ind = '',0
  Write_out(out, example_cluster_annot_file)
  out, ind = '',0
  print(example_cluster_annot_file)
  return()

def dfgs():
  fh=open(seq_file,"r")
  out=''
  for header,sequence in fasta_iterator(fh):
    if(header.split("__")[0] in clusters[max_c]):
      print(header)
      out=out+">"+header+"\n"+sequence+"\n"
  fh.close()
  out_file = example_cluster_annot_file+"_Clust_"+max_c+".txt"
  fh=open(out_file,"w")
  fh.write(out)
  fh.close()
  print(out_file)
  return()

def Get_clusters_specific(sample, annot_file, cluster_file,per_cluster_developmental_classification_file,example_cluster_annot_file,seq_file):
  Initialise_clusters (sample, annot_file, cluster_file,per_cluster_developmental_classification_file,example_cluster_annot_file,seq_file)
  Generate_trees(sample, annot_file, cluster_file,per_cluster_developmental_classification_file,example_cluster_annot_file,seq_file)
  return()

def Generate_trees(sample, annot_file, cluster_file,per_cluster_developmental_classification_file,example_cluster_annot_file,seq_file):
  clusters_of_isotype_of_interest = example_cluster_annot_file+"_ALL.txt"
  fh=open(clusters_of_isotype_of_interest,"r")
  clusters =  []
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      clust= l[1]
      clusters.append(clust)
  fh.close()
  annots = Get_annotation(annot_file)
  out = "#sample\tcluster\tmin_mut_id\tmin_mut\n"
  clusters_of_isotype_of_interest_mutations = example_cluster_annot_file+"_MUTS.txt"
  fh = open(clusters_of_isotype_of_interest_mutations,"w")
  fh.write(out)
  fh.close()
  for c in clusters:
    infile = example_cluster_annot_file+"_"+c+".txt"
    aligned = example_cluster_annot_file+"_"+c+"_ALIGNED.fasta"
    tree_file = example_cluster_annot_file+"_"+c+"_TREE.txt"
    command="ALIGN"
    if(command=="ALIGN"):
      insert = ''
      #if( len(clusters[c]) > 5000):insert = '--parttree'
      command = "/software/pubseq/bin/mafft-6.857/bin/mafft --retree 2 --op 100 "+insert+" "+infile+" > "+aligned
      subprocess.getoutput(command)
    command="TREE"
    if(command=="TREE"):
      command = "FastTree -nt "+aligned+" > "+tree_file
      subprocess.getoutput(command)
    command = "MIN_MUTS"
    if(command=="MIN_MUTS"):
      ids = []
      fh=open(infile, "r")
      for header,sequence in fasta_iterator(fh):
        id = header.split("__")[0]
        if(id in annots):
          mm = annots[id][3]+annots[id][2]
          ids.append([id, mm])
      fh.close()
      ids.sort(key=lambda x: x[1])
      min_mut = ids[0][1]
      out = ''
      for i in range(0,len(ids)):
        if(ids[i][1]==min_mut):
          out = out+sample+"\t"+c+"\t"+ids[i][0]+"\t"+str(ids[i][1])+"\n"
      Write_out(out, clusters_of_isotype_of_interest_mutations)
    command="iMRCA"
    if(command=="MRCA"):
      from Bio import Phylo
      import networkx, pylab
      tree = Phylo.read(tree_file, 'newick')
      print(tree)
      #net = Phylo.to_networkx(tree)
      #networkx.draw(net)
      #pylab.show()
      break
    print(c)
  return()

def Initialise_clusters (sample, annot_file, cluster_file,per_cluster_developmental_classification_file,example_cluster_annot_file,seq_file):
  seqs,freqs,alias = Get_sequences(seq_file)
  clusters,cluster_inv,cluster_information = Tree(),{},{}
  fh=open(per_cluster_developmental_classification_file,"r")
  n_isotypes = []
  for i in range(4,12):
    n_isotypes.append(str(i)+"_isotype(s)")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cluster = l[10]
      id = l[1]
      classes = l[2].split("|")
      p1,p2 = 0,0
      for t in classes:
        if(t.count(",")==0 and t.count("IGH")==1):
          clusters[cluster][t][id].value = 1
          #if(t.count("IGHE")!=0):p1 = 1
          #if(t.count("IGHG3")!=0):p2 = 1
        #if(t.count("IGHD")!=0 and t.count("IGHM")!=0):
        #  p2 = 1
        #if(t in n_isotypes):
        #  p1 = 1
  fh.close()
  out1 = "#sample\tcluster\n"
  # All clusters with multiple IgE and IgG3
  for c in clusters:
    if(len(clusters[c])>=2):
      if("IGHE" in clusters[c]):# and "IGHG3" in clusters[c]):
        count= 0
        ids = {}
        for t in  clusters[c]:
          for id in  clusters[c][t]:
            count = count+1#len(clusters[c][t])
            ids[id] = 1
        if(len( clusters[c]["IGHE"])>=5 and len( clusters[c]["IGHA"])>=5):
          if(len(ids)<500):
            out1=out1+sample+"\t"+c+"\n"
            out_file1 = example_cluster_annot_file+"_"+c+".txt"
            out = ''
            for id in ids:
              out=out+">"+alias[id]+"\n"+seqs[id]+"\n"
            fh=open(out_file1,"w")
            fh.write(out)
            fh.close()
            print(c, count)#, clusters[c]
  out_file1 = example_cluster_annot_file+"_ALL.txt"
  fh=open(out_file1,"w")
  fh.write(out1)
  fh.close()
  return()

def Get_network_parameters_per_sequence_classification(sample, cluster_file, developmental_classification_file, per_sequence_classification_network_parameters,reverse_primer_group):
  classification_ids = Get_classification(developmental_classification_file)
  fh = open(cluster_file,"r")
  chains = []
  chain_counts, total = {},0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,freq, id_short,clust = l[2],list(map(int, l[2].split("__")[1].split("|")[0].split("_"))),l[2].split("__")[0],l[1]
      if(len(chains)==0):chains = l[2].split("|")[1].split("_")
      id_short = id.split("__")[0]
      if(id_short in classification_ids):
        total = total+sum(freq)
        nz = [i for i in range(len(freq)) if freq[i]!=0]
        freqs = {}
        for i in nz:
          c,f = chains[i],freq[i]
          c = c.split("*")[0]
          if(reverse_primer_group!="ISOTYPER"):
            if(c in ["IGHA1","IGHA2"]):c = "IGHA"
            elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
          if(c in freqs):freqs[c] = freqs[c]+freq[i]
          else:freqs[c]=freq[i]
          if(c not in ["IGHM","IGHD"]):
            c1 = "Class_switched"
            if(c1 in freqs):freqs[c1] = freqs[c1]+freq[i]
            else:freqs[c1]=freq[i]
          else:
            if("IGHD,IGHM_unmutated" in classification_ids[id_short]):c1 = "IGHD,IGHM_unmutated"
            elif("IGHD,IGHM_mutated" in classification_ids[id_short]):c1 = "IGHD,IGHM_mutated"
            else:print(c, id, classification_ids[id_short])
            if(c1 in freqs):freqs[c1] = freqs[c1]+freq[i]
            else:freqs[c1]=freq[i]
        for c in freqs:
          if(c in chain_counts):chain_counts[c]=chain_counts[c]+[freqs[c]]
          else:chain_counts[c]=[freqs[c]]
          #if(len(chain_counts[c])>1000):break ################### remove
  fh.close()
  out="#Id\tIsotype\tN reads\tN vertices\tVertex Gini Index\ttotal_reads\n"
  for c in chain_counts:
    (vpoints,vvdf)=VDF(chain_counts[c])
    (vgini)=Get_Gini(vpoints,vvdf)
    out=out+sample+"\t"+c+"\t"+str(sum(chain_counts[c]))+"\t"+str(len(chain_counts[c]))+"\t"+str(vgini)+"\t"+str(total)+"\n"
  fh=open(per_sequence_classification_network_parameters,"w")
  fh.write(out)
  fh.close()
  return()

def Get_classification(developmental_classification_file):
  fh = open(developmental_classification_file,"r")
  classification_ids = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,classes = l[1].split("__")[0],l[2].split("|")#,int(l[10])
      for c in classes:
        if(c in ["IGHA1","IGHA2"]):c = "IGHA"
        elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        if(id in classification_ids):
          if(c not in classification_ids[id]):
            classification_ids[id] = classification_ids[id]+[c]
        else:classification_ids[id]=[c]
      if("IGHM" in classes or "IGHD" in classes):
        if("IGHD,IGHM_unmutated" not in classes and "IGHD,IGHM_mutated" not in classes):
          mut = int(l[6])
          if(mut==0):c = "IGHD,IGHM_unmutated"
          else:c = "IGHD,IGHM_mutated"
          classification_ids[id] = classification_ids[id]+[c]
  fh.close()
  return(classification_ids)

def Subsample_repertoire(subsample_file, sample_depth, seq_file,samplex,cluster_file):
  n_repeats = 40
  clusters = {}
  fh=open(cluster_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      clusters[l[2].split("__")[0]] = l[1]
  fh.close()
  fh=open(seq_file,"r")
  sequences,id_seqs = [],{}
  for header,sequence in fasta_iterator(fh):
    id = header.split("__")[0]
    if(id in clusters):
      freq,classes = list(map(int, header.split("__")[1].split("|")[0].split("_"))), header.split("|")[1].split("_")
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      for i in nz: 
        for j in range(0,freq[i]):
          sequences.append([sequence,classes[i],id])
  fh.close()
  fh=open(subsample_file,"w")
  fh.write("#repeat\tcluster\tid\ttotal frequency\n")
  fh.close()
  if(sample_depth<len(sequences)):
    classes_dict = {}
    for i in range(0,len(classes)):
      classes_dict[classes[i]] = i
    classes_lab = "_".join(classes)
    for r in range(0,n_repeats):
      rep_id = "REP"+str(r)+"_"+samplex
      print("\t",sample_depth)
      sample = random.sample(sequences, sample_depth)
      seqs = Functions.Tree()
      out,ind='',0
      for s in sample:
        seqs[s[0]][s[1]][s[2]][ind].value = 1
        ind = ind+1
      for s in seqs:
        f=[0]*len(classes)
        for clas in seqs[s]:
          for header in seqs[s][clas]:
            f[classes_dict[clas]] = f[classes_dict[clas]]+len(seqs[s][clas][header])
            break
        cluster = clusters[header]
        header = header+"__"+"_".join(map(str, f))+"|"+classes_lab
        out=out+rep_id+"\t"+cluster+"\t"+header+"\t"+str(sum(f))+"\n"
      Write_out(out, subsample_file)
      out=''
  return()

def Get_network_parameters_per_classificiation_subsample_clones(sample,per_cluster_developmental_classification_file, per_cluster_developmental_classification_network_parameters_subsampled_CLONE,reverse_primer_group,cluster_file,subsample_depth_file,developmental_classification_file):
  Raw_values_Get_network_parameters_per_classificiation_subsample_clones(sample,per_cluster_developmental_classification_file, per_cluster_developmental_classification_network_parameters_subsampled_CLONE,reverse_primer_group,cluster_file,subsample_depth_file,developmental_classification_file)
  Summarise_Get_network_parameters_per_classificiation_subsample_clones(sample,per_cluster_developmental_classification_file, per_cluster_developmental_classification_network_parameters_subsampled_CLONE,reverse_primer_group,cluster_file,subsample_depth_file,developmental_classification_file)
  return()

def Get_random_sample(arr, n_iter=None, sample_size=10, fast=True):
    if fast:
        # find the index we last sampled from
        start_idx = (n_iter * sample_size) % n
        if start_idx + sample_size >= n:
            # shuffle array if we have reached the end and repeat again
            np.random.shuffle(arr)
        return arr[start_idx:start_idx+sample_size] 
    else:
        return np.random.choice(arr, sample_size, replace=False)
    
def Collect_samples(arr,sample_size,n_samples,fast=False):
    samples = np.zeros((n_samples + 1, sample_size), np.int32)
    for sample_n in range(0, n_samples):
        sample = Get_random_sample(arr, n_iter=sample_n,sample_size=sample_size,fast=fast)
        samples[sample_n] = sample
    return samples

def Summarise_Get_network_parameters_per_classificiation_subsample_clones(sample,per_cluster_developmental_classification_file, per_cluster_developmental_classification_network_parameters_subsampled_CLONE,reverse_primer_group,cluster_file,subsample_depth_file,developmental_classification_file):
  fh=open(per_cluster_developmental_classification_network_parameters_subsampled_CLONE, "r")
  X,ids = [],{}
  ind = 0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id = l[0]+"\t"+l[2]
      values = list(map(float, l[3:len(l)]))
      X.append(values)
      if(id in ids):ids[id] = ids[id]+[ind]
      else:ids[id]=[ind]
      ind = ind+1
  fh.close()
  X = np.array(X)
  out="#Id\tIsotype\tN repeats\tN reads\tN vertices\tN_clusters\tVertex Gini Index\tCluster Gini Index\tCluster Reads Gini Index\tmean_total_cluster_size\tmean_vertex_size\tmax_clust_pop\tmax_vertex_pop\tcluster_Thiel0\tCluster_Thiel1\tsd_Vertex Gini Index\tsd_Cluster Gini Index\n"
  out="#Id\tIsotype\tN repeats\tN reads\tN vertices\tN_clusters\tVertex Gini Index\tCluster Gini Index\tmean_total_cluster_size\tmean_vertex_size\tmax_clust_pop\tmax_vertex_pop\tVertex Reyni\tCluster_Renyi\tThielC1\tThielC2\n"
  for id in ids:
    w = ids[id]
    w = X[w,:]
    values = []
    for i in range(0,len(w[0,:])):
      values.append(np.mean(w[:,i]))
    k = 2
    inde = len(w[0,:])-1
    #weighted_mean_c=  sum(w[:,5]*(w[:,1])**k)/sum((w[:,1])**k)
    #weighted_mean_c=  sum(w[:, inde]*(w[:,2])**k)/sum((w[:,2])**k)
    #weighted_mean_v=  sum(w[:,4]*w[:,0])/sum(w[:,0])
    #values = values+[np.std(w[:,4]), np.std(w[:,5])]
    #values = values+[weighted_mean_v, weighted_mean_c]
    out=out+"\t".join(map(str, [id,len(w[:,0])]+values))+"\n"
  fh=open(per_cluster_developmental_classification_network_parameters_subsampled_CLONE.replace("parameters_SUBSAMPLED", "parameters_SUMMARY_SUBSAMPLED"), "w")
  fh.write(out)
  fh.close()
  print(out)
  return()

def Raw_values_Get_network_parameters_per_classificiation_subsample_clones(sample,per_cluster_developmental_classification_file, per_cluster_developmental_classification_network_parameters_subsampled_CLONE,reverse_primer_group,cluster_file,subsample_depth_file,developmental_classification_file):
  sample1 = sample
  fh = open(developmental_classification_file,"r")
  muts1,muts = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id, mut = l[1], int(l[6])
      muts1[id] = mut
  fh.close()
  fh = open(per_cluster_developmental_classification_file,"r")
  classification_ids,classification_clusters = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id, cluster, mut, types = l[1],l[10], float(l[6]),l[2].replace("IGHA","IGHA1/2").replace("IGHA1/21/2","IGHA1/2").split("|")
  #    classification_clusters[cluster] = l[2].split("|")+["ALL"]
      if(id in muts1):
        muts[id] = [muts1[id], cluster, types]
  #    else:print id,"FAIL"
  fh.close()
  del muts1
  print("READ DEVELOP FILE")
  fh=open(subsample_depth_file,"r")
  subsample_depth,subsample_depth_clone = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      subsample_depth[l[0]] = int(l[2])
      if(len(l)>3):subsample_depth_clone[l[0]] = int(l[3])
  fh.close()
  subsample_depth["ALL"] = 1250
  subsample_depth_clone["ALL"] = 2000
  if(sample.count("MUR")!=0):
    subsample_depth["ALL"] = 1000
    subsample_depth_clone["ALL"] = 1000
  fh=open(cluster_file,"r")
  v_sizes_all={}
  clusters = Tree()
  clusters_array, ids_clusters1,clas_clusters = [],[],[]
  index = 0
  for l in fh:
    index=index+1
    if (index>1):
      l=l.strip().split()
      id,freq, id_short,clust = l[2],list(map(int, l[2].split("__")[1].split("|")[0].split("_"))),l[2].split("__")[0], l[1]
      class_raw = l[2].split("|")[1].split("_")
      #classification = classification_clusters[clust]
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      freq_clas = {}
      for i in nz:
        c=class_raw[i].split("*")[0]
        #if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
        #if(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        if(c in freq_clas):freq_clas[c] = freq_clas[c]+freq[i]
        else:freq_clas[c]=freq[i]
        if(c not in ["IGHD","IGHM"]):
          c = "Class_switched"
          if(c in freq_clas):freq_clas[c] = freq_clas[c]+freq[i]
          else:freq_clas[c]=freq[i]
        else:
          m = muts[id_short][0]
          if(m<=3): c = "IGHD_IGHM_unmutated"
          else:c = "IGHD_IGHM_mutated"
          if(c in freq_clas):freq_clas[c] = freq_clas[c]+freq[i]
          else:freq_clas[c]=freq[i]
      freq_clas["ALL"]=sum(freq)
      for clas in freq_clas:
        if(clas in subsample_depth):
          clusters_array.append(clust)
          ids_clusters1.append(id_short)
          clas_clusters.append(clas)
          v_sizes_all[clas+":"+id_short] = freq_clas[clas]
      #if(index>5000):break
  fh.close()
  clusters_array, ids_clusters1,clas_clusters = np.array(clusters_array), np.array(ids_clusters1),np.array(clas_clusters)
  clas_clusters_unique = np.unique(clas_clusters)
  clusters, ids_clusters = {},{}
  for clas in clas_clusters_unique:
    ii = np.where(clas_clusters == clas)[0]
    clusters[clas] = clusters_array[ii]
    ids_clusters[clas] = ids_clusters1[ii]
  del clusters_array, ids_clusters1,clas_clusters
  fh=open(per_cluster_developmental_classification_network_parameters_subsampled_CLONE, "w")
  fh.close()
  out="#Id\trepeat\tIsotype\tN reads\tN vertices\tN_clusters\tVertex Gini Index\tCluster Gini Index\tmean_total_cluster_size\tmean_vertex_size\tmax_clust_pop\tmax_vertex_pop\tVertex Reyni\tCluster_Renyi\tThielC1\tThielC2\n"
  ind = 0
  repeats =100
  for clas in clusters:
    clas1 = clas
    if(clas=="IGHA"):clas1="IGHA1/2"
    if(clas1 in subsample_depth):
      print(clas,"\t",len(clusters[clas]), len(np.unique(clusters[clas])))
      if(subsample_depth[clas1] <= len(np.unique(clusters[clas]))):
        clusters_sub = np.array(clusters[clas])
        ids_sub = np.array(ids_clusters[clas])
        depth = subsample_depth[clas1]
        sample_indices = np.arange(len(ids_sub)) 
        unique_elements, c_sizes = np.unique(np.array(clusters_sub), return_counts=True)
        print(len(np.unique(clusters_sub)))
        for r in range(0,repeats):
          #print r
          n = len(sample_indices)
          sample = np.random.choice(sample_indices,n,replace = False)
          #sample = Collect_samples(sample_indices, n, 1, fast=True)
          start,end_vertex,v_finished = 0,0,0
          clusters_sampled = np.array([])
          length = len(sample)
          for i in range(depth,len(sample)):
            end = i
            #new = np.array(clusters_sub)[start:end]
            clusters_sampled = np.unique(clusters_sub[sample[0:end]])
            #if(v_finished==0):
            #  vertex_sampled = np.unique(ids_sub[sample[0:end]])
            #  if(len(vertex_sampled)>=depth):
            #    end_vertex = end
            #    v_finished = 1
            #start = end
            #print len(clusters_sampled),depth#, len(vertex_sampled)
            if(len(clusters_sampled)>=depth):break
          sample = sample[0:end]
          #print clas, len(clusters_sampled), len(sample), depth
          clusters_samples = clusters_sub[sample[0:end]]
          vertices_samples = ids_sub[sample[0:depth]]
          #print len(np.unique(vertices_samples))
          c_sizes,v_sizes = [],[]
          unique, c_sizes = np.unique(clusters_samples, return_counts=True)
          for v in range(0,len(vertices_samples)):
            v_sizes.append(v_sizes_all[clas1+":"+vertices_samples[v]])
          #unique_elements, c_sizes = np.unique(np.array(clusters_samples), return_counts=True)
          #clusts_reads = {}
          #for i in range(0,len(clusters_samples)):
          #  v_sizes.append(v_sizes_all[clas1+":"+vertices_samples[i]])
          #  c = clusters_samples[i]
          #  if(c in clusts_reads):clusts_reads[c] = clusts_reads[c]+v_sizes_all[clas1+":"+vertices_samples[i]]
          #  else:clusts_reads[c] = v_sizes_all[clas1+":"+vertices_samples[i]]
          #for c in clusts_reads:
          #  cluster_size_reads.append(clusts_reads[c])
          #print clas,sort( c_sizes).tolist()
          vgini = Get_Gini_fast(v_sizes)
          cgini = Get_Gini_fast(c_sizes)
          #c_readsgini = Get_Gini_fast(cluster_size_reads)
          #print clas, "\t",cgini, "\t",cgini_clone_sampling
          renyi_v = Get_Renyi(v_sizes)
          renyi_c = Get_Renyi(c_sizes)
          T0c, T1c = Theil_indices_fast(c_sizes)
          mean_total_cluster_size = mean(c_sizes)
          mean_vertex_size = mean(v_sizes)
          max_cpop, max_vpop =max(v_sizes)*100.0/sum(v_sizes), max(c_sizes)*100.0/sum(c_sizes)
          #prop_clusters_with_vertices_greater_than_5 = len([1 for i in max_c_size if i > 5])*100.0/len(max_c_size)
          #prop_clusters_with_vertices_greater_than_10 = len([1 for i in max_c_size if i > 10])*100.0/len(max_c_size)
          #prop_clusters_with_vertices_greater_than_20 = len([1 for i in max_c_size if i > 20])*100.0/len(max_c_size)
          out=out+"\t".join(map(str,[sample1,"REP"+str(r), clas1, sum(v_sizes),len(v_sizes), len(c_sizes), vgini,cgini, mean_total_cluster_size, mean_vertex_size, max_cpop, max_vpop, renyi_v, renyi_c,T0c, T1c ]))+"\n"
          #out=out+"\t".join(map(str,[sample1,"REP"+str(r), clas1, sum(v_sizes),len(v_sizes), len(c_sizes), vgini,cgini,c_readsgini, mean_total_cluster_size, mean_vertex_size, max_cpop, max_vpop, prop_clusters_with_vertices_greater_than_5, prop_clusters_with_vertices_greater_than_10, prop_clusters_with_vertices_greater_than_20 , T0c, T1c ]))+"\n"
          ind = ind+1
          if(ind>500):
            Write_out(out, per_cluster_developmental_classification_network_parameters_subsampled_CLONE)
            out, ind = '',0
          #del c_sizes,v_sizes,max_c_size,cluster_size_reads
  Write_out(out, per_cluster_developmental_classification_network_parameters_subsampled_CLONE)
  out, ind = '',0
  return()

def Get_network_parameters_per_classificiation_subsample(samplex,per_cluster_developmental_classification_file,per_cluster_developmental_classification_network_parameters,reverse_primer_group,cluster_file,subsample_depth_file,developmental_classification_file):
  fh = open(developmental_classification_file,"r")
  muts1,muts = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id, mut = l[1], int(l[6])
      muts1[id] = mut
  fh.close()
  fh = open(per_cluster_developmental_classification_file,"r")
  classification_ids,classification_clusters = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id, cluster, mut, types = l[1],l[10], float(l[6]),l[2].replace("IGHA","IGHA1/2").replace("IGHA1/21/2","IGHA1/2").split("|")
      if(id in muts1):
        muts[id] = [muts1[id], cluster, types]
      else:
        print(id,"FAIL")
  fh.close()
  del muts1
  print("READ DEVELOP FILE")
  fh=open(cluster_file,"r")
  v_sizes,c_sizes = {},Tree()
  index = 0
  sample_tree,reads_per_sample = Tree(),{}
  cluster_size_cumul = {}
  inv_cluster ={}
  cluster_per_chain ={}
  for l in fh:
    index=index+1
    if (index>1):
      l=l.strip().split()
      id,freq, id_short,clust = l[2],list(map(int, l[2].split("__")[1].split("|")[0].split("_"))),l[2].split("__")[0], l[1]
      clust = l[1]
      classification = l[2].split("|")[1].split("_")
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      classes = {}
      inv_cluster[id_short]=clust
      for i in nz:
        c = classification[i].split("*")[0]
        if(reverse_primer_group!="ISOTYPER"):
          if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
          elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
          elif(c.count("TRBC")==1):c = "TRBC"
        if(c in classes):classes[c] = classes[c]+freq[i]
        else:classes[c]=freq[i]
        if(c in ["IGHD", "IGHM"]):
          m = muts[id_short][0]
          if(m<4):c1 = "IGHD,IGHM_unmutated"
          else:c1 ="IGHD,IGHM_mutated"
        else:c1 = "Class_switched"
        if(c1 in classes):classes[c1] = classes[c1]+freq[i]
        else:classes[c1]=freq[i]
      for c1 in classes:
        c = c1
        if(c in v_sizes):
          v_sizes[c] = v_sizes[c]+[classes[c1]]
          reads_per_sample[c] = reads_per_sample[c] +classes[c1]
        else:
          reads_per_sample[c]=classes[c1]
          v_sizes[c]=[classes[c1]]
        if(c in cluster_per_chain):cluster_per_chain[c]=cluster_per_chain[c]+[id_short]*classes[c1]
        else:cluster_per_chain[c]=[id_short]*classes[c1]
      classes1 = muts[id_short][2]
      for c1 in classes1:
        c = c1
  fh.close()
  fh=open(subsample_depth_file,"r")
  subsample_depth = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      subsample_depth[l[0]] = int(l[1])
  fh.close()
  out="#Id\tIsotype\tN reads\tN vertices\tVertex Gini Index\tCluster Unique Gini Index\tCluster Total Gini Index\tLargest Cluster (total %)\t2nd Largest Cluster (total %)\tLargest Cluster (unique %)\t2nd Largest Cluster (unique %)\tnumber of unique clusters\ttotal reads per sample\tsample\tmean_unique_cluster_size\tmean_total_cluster_size\n"
  repeats =20
  for c in v_sizes:
    if(c in subsample_depth):
      if(subsample_depth[c]>0):
        depth = subsample_depth[c]
        v = v_sizes[c]
        print(sum(v), depth)
        if(sum(v)>depth):
          prob = np.array(v)/(sum(v)*1.0)
          xlen = len(v)
          for r in range(0,repeats):
            sample = np.random.choice(xlen, depth,p=prob)
            unique, counts = np.unique(sample, return_counts=True)
            vpoints,vvdf = np.unique(counts, return_counts=True)
            vgini=Get_Gini(vpoints,vvdf)

            if(c in cluster_per_chain):
              cluster_sizes_sub = cluster_per_chain[c]
              sample_cs = np.random.choice(cluster_sizes_sub, depth)
              unique, counts = np.unique(sample_cs, return_counts=True)
              unique_BCRs_per_cluster, total_BCRs_per_cluster = {},{}
              unique_cluster_sizes, total_cluster_sizes = [],[]
              for i in range(0,len(unique)):
                id = unique[i]
                c1 = inv_cluster[id]
                if(c1 in unique_BCRs_per_cluster):
                  unique_BCRs_per_cluster[c1] = unique_BCRs_per_cluster[c1]+1
                  total_BCRs_per_cluster[c1] = total_BCRs_per_cluster[c1]+counts[i]
                else:
                  unique_BCRs_per_cluster[c1]=1
                  total_BCRs_per_cluster[c1]=counts[i]
              for c1 in unique_BCRs_per_cluster:
                unique_cluster_sizes = unique_cluster_sizes+[unique_BCRs_per_cluster[c1]]
                total_cluster_sizes=total_cluster_sizes+[total_BCRs_per_cluster[c1]]
              mean_unique_cluster_size,mean_total_cluster_size = mean(unique_cluster_sizes), mean(total_cluster_sizes)
              unique_uniq_cs, counts_uniq_cs = np.unique(unique_cluster_sizes, return_counts=True)
              unique_total_cs, counts_total_cs = np.unique(total_cluster_sizes, return_counts=True)
              uniq_clust_gini=Get_Gini(unique_uniq_cs,counts_uniq_cs)
              total_clust_gini=Get_Gini(unique_total_cs, counts_total_cs)
              max_pop, max_1_pop = unique_total_cs[len(unique_total_cs)-1]*100.0/sum(total_cluster_sizes), unique_total_cs[len(unique_total_cs)-2]*100.0/sum(total_cluster_sizes)
              max_pop_uniq, max_1_pop_uniq = unique_uniq_cs[len(unique_uniq_cs)-1]*100.0/sum(unique_cluster_sizes), unique_uniq_cs[len(unique_uniq_cs)-2]*100.0/sum(unique_cluster_sizes)
            else:
              max_pop, max_1_pop,max_pop_uniq, max_1_pop_uniq  = -1,-1,-1,-1
              uniq_clust_gini,total_clust_gini = -1,-1
            out = out+"REP"+str(r)+"\t"+c+"\t"+str(len(sample))+"\t"+str(len(unique))+"\t"+str(vgini)+"\t"+str(uniq_clust_gini)+"\t"+str(total_clust_gini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\t"+str(max_pop_uniq)+"\t"+str(max_1_pop_uniq)+"\t"+str(len(unique))+"\t"+str(sum(v))+"\t"+samplex+"\t"+str(mean_unique_cluster_size)+"\t"+str(mean_total_cluster_size)+"\n"
  fh=open(per_cluster_developmental_classification_network_parameters, "w")
  fh.write(out)
  fh.close()
  return()

def Get_expanded_cluster_summary(cluster_SUBSAMPLE_information_file,cluster_isotype_expansion_SUBSAMPLE_file,sample,developmental_classification_file,per_cluster_developmental_classification_file):
  fh = open(developmental_classification_file,"r")
  muts,muts_ids = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,mut = l[1],float(l[6])
      muts_ids[id] = mut
  fh.close()
  fh = open(per_cluster_developmental_classification_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,clust = l[1],l[10]
      if(clust in muts):muts[clust] =muts[clust] +[muts_ids[id]]
      else:muts[clust] = [muts_ids[id]]
  fh.close()
  del muts_ids
  fh=open(cluster_SUBSAMPLE_information_file,"r")
  total_BCRs_per_isotype,unique_BCRs_per_isotype, singletons_per_isotype = {},{},{}
  expanded_clonal_range, expanded_clonal_SD={},{}
  perc_exp_3, perc_exp_5, perc_exp_10, perc_exp_20 = {},{},{},{}
  unmutated_singletons, unmutated_perc_exp_3 = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      rep, isotype =l[1],l[3]
      reads_per_isotype = int(l[7])
      range_muts, n_unique_BCRs,sd_muts, total_cluster_size = int(l[5]), int(l[8]),float(l[7]),int(l[9])
      typ = rep+"\t"+isotype
      mut = -1
      clust = l[2]
      if(clust in muts):mut = mean(muts[clust])
      if(typ in total_BCRs_per_isotype):
        total_BCRs_per_isotype[typ] = total_BCRs_per_isotype[typ] + reads_per_isotype
        unique_BCRs_per_isotype[typ] = unique_BCRs_per_isotype[typ] + n_unique_BCRs
      else:
        total_BCRs_per_isotype[typ] = reads_per_isotype
        unique_BCRs_per_isotype[typ] = n_unique_BCRs
      if(total_cluster_size==reads_per_isotype and n_unique_BCRs==1):
        if(typ in singletons_per_isotype):singletons_per_isotype[typ] = singletons_per_isotype[typ]+1
        else:singletons_per_isotype[typ] = 1
        if(mut!=-1 and mut<=4):
          if(typ in unmutated_singletons):unmutated_singletons[typ] = unmutated_singletons[typ]+1
          else:unmutated_singletons[typ] = 1
      if(n_unique_BCRs>=3):
        if(mut!=-1 and mut<=4):
          if(typ in unmutated_perc_exp_3):unmutated_perc_exp_3[typ] = unmutated_perc_exp_3[typ]+reads_per_isotype
          else:unmutated_perc_exp_3[typ]=reads_per_isotype
        if(typ in expanded_clonal_range):
          expanded_clonal_range[typ] = expanded_clonal_range[typ]+[range_muts*1.0/n_unique_BCRs]
          expanded_clonal_SD[typ] = expanded_clonal_SD[typ]+[sd_muts]
          perc_exp_3[typ] = perc_exp_3[typ]+reads_per_isotype
        else:
          expanded_clonal_range[typ]=[range_muts*1.0/n_unique_BCRs]
          expanded_clonal_SD[typ]=[sd_muts]
          perc_exp_3[typ] = reads_per_isotype
        if(n_unique_BCRs>=5):
          if(typ in perc_exp_5):perc_exp_5[typ] = perc_exp_5[typ]+reads_per_isotype
          else:perc_exp_5[typ]=reads_per_isotype
        if(n_unique_BCRs>=10):
          if(typ in perc_exp_10):perc_exp_10[typ] = perc_exp_10[typ]+reads_per_isotype
          else:perc_exp_10[typ]=reads_per_isotype
        if(n_unique_BCRs>=20):
          if(typ in perc_exp_20):perc_exp_20[typ] = perc_exp_20[typ]+reads_per_isotype
          else:perc_exp_20[typ]=reads_per_isotype
  fh.close()
  out = "#sample\tsubsample\tisotype\t% singletons (total BCRs)\t% singletons (unique)\tnumber of expanded clusters\tmean_cluster mutational range (bp/unique sequences)\tSD cluster mutational range\ttotal_BCRs_per_isotype\tunique_BCRs_per_isotype\tn_in_expanded_3_clusters\tn_in_expanded_5_clusters\tn_in_expanded_10_clusters\tn_in_expanded_20_clusters\tunmutated_singletons\tunmutated_perc_expanded_3_clusters\n"
  for typ in total_BCRs_per_isotype:
    expanded_clonal_ranges, expanded_clonal_SDs,length = -1,-1,0
    singletons1,singletons2 = 0,0
    if(typ in expanded_clonal_range):
      length = len(expanded_clonal_range[typ])
      expanded_clonal_ranges = mean(expanded_clonal_range[typ])
      expanded_clonal_SDs = mean(expanded_clonal_SD[typ])
    if(typ in singletons_per_isotype):
      singletons1 = singletons_per_isotype[typ]*100.0/total_BCRs_per_isotype[typ]
      singletons2 = singletons_per_isotype[typ]*100.0/unique_BCRs_per_isotype[typ]
    array = []
    for x in [perc_exp_3, perc_exp_5, perc_exp_10, perc_exp_20,unmutated_singletons, unmutated_perc_exp_3 ]:
      if(typ in x):array.append(x[typ])
      else:array.append(0)
    out=out+sample+"\t"+typ+"\t"+"\t".join(map(str, [singletons1, singletons2, length, expanded_clonal_ranges, expanded_clonal_SDs,total_BCRs_per_isotype[typ],unique_BCRs_per_isotype[typ]]+array))+"\n"
  fh=open(cluster_isotype_expansion_SUBSAMPLE_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_V_gene_isotype_frequency_grouped(sample, annot_file, V_gene_isotype_frequency_file_grouped,per_cluster_developmental_classification_file,CDR3_length_file_grouped,per_V_gene_cluster_file,per_V_gene_cluster_summary_file):
  print(V_gene_isotype_frequency_file_grouped)
  seqs,freqs,alias = Get_sequences(seq_file)
  CDR3s,mutations = Get_annotations(annot_file)
  counts,CDR3_lengths,counts_uniq = {},{},{}
  d_counts, j_counts = {},{}
  cluster_v = {}
  fh=open(per_cluster_developmental_classification_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id, alias1,cluster = l[1],l[2].split("|"),l[10]
      v1,j1 = l[4],l[5]
      if(id in mutations):
        freq = sum(map(int, alias[id].split("__")[1].split("|")[0].split("_")))
        v = mutations[id][3].split("\t")
        cdr3= mutations[id][2].split(":")[0]
        cdr3 = len(cdr3)-2
        v,j = v[0], v[1]
        #if(v!=v1 or j!=j1):print id, v,j, v1,j1
        freq1 = 1
        c1 = "ALL"+"\t"+v
        if(c1 in counts):
          counts[c1] = counts[c1]+freq
          counts_uniq[c1] = counts_uniq[c1]+1
        else:
          counts[c1]=freq 
          counts_uniq[c1] = 1
        c1 = "ALL"
        if(c1 in CDR3_lengths):CDR3_lengths[c1] = CDR3_lengths[c1] + [cdr3]*freq1
        else:CDR3_lengths[c1] = [cdr3]*freq1
        for a in alias1: 
          p = 1
          if(a.count("isotype")!=0):p = 0
          if(a.count(",")!=0 and a.count("mut")==0):p=0
          if(p==1):
            c1 = a+"\t"+v
            if(c1 in counts):
              counts[c1] = counts[c1]+freq
              counts_uniq[c1] = counts_uniq[c1]+1
            else:
              counts[c1]=freq
              counts_uniq[c1] = 1
            if(c1 in cluster_v):cluster_v[c1] = cluster_v[c1]+[cluster]
            else:cluster_v[c1]=[cluster]
            c1 = a+"\t"+j
            if(c1 in j_counts):j_counts[c1] = j_counts[c1]+1
            else:j_counts[c1] = 1
            c1 = a
            if(c1 in CDR3_lengths):CDR3_lengths[c1] = CDR3_lengths[c1] + [cdr3]*freq1
            else:CDR3_lengths[c1] = [cdr3]*freq1
        if("IGHD,IGHM_unmutated" in alias1):
          if(int(l[8])<=5 and int(l[8])==freq):
            c1 = "IGHD,IGHM_unmutated_singleton"+"\t"+v
            if(c1 in counts):
              counts[c1] = counts[c1]+freq
              counts_uniq[c1] = counts_uniq[c1]+1
            else:
              counts[c1]=freq
              counts_uniq[c1] = 1
            if(c1 in cluster_v):cluster_v[c1] = cluster_v[c1]+[cluster]
            else:cluster_v[c1]=[cluster]
  fh.close()
  out= "#sample\tclass\tV gene\tn_seq\tn_clusters\tmax_cluster_size\tmean_cluster_size\n"
  for c in cluster_v:
    unique, counts1 = np.unique(cluster_v[c], return_counts=True)
    n_seq,n_clust = len(cluster_v[c]), len(unique)
    max_cluster = max(counts1)*100.0/n_seq
    mean_cluster = mean(counts1)*100.0/n_seq
    out = out+"\t".join(map(str,[sample,c,n_seq,n_clust,max_cluster,mean_cluster]))+"\n"
  fh=open(per_V_gene_cluster_summary_file,"w")
  fh.write(out)
  fh.close()
  out = "#sample\tclass\tV gene\tfrequency\tuniq_read_freq\n"
  for c in counts:
    out=out+sample+"\t"+c+"\t"+str(counts[c])+"\t"+str(counts_uniq[c])+"\n"
  fh=open(V_gene_isotype_frequency_file_grouped,"w")
  fh.write(out)
  fh.close()
  out = "#sample\tclass\tmean_CDR3 length\tcount\n"
  for c in CDR3_lengths:
    out=out+sample+"\t"+c+"\t"+str(mean(CDR3_lengths[c]))+"\t"+str(len(CDR3_lengths[c]))+"\n"
  fh=open(CDR3_length_file_grouped,"w")
  fh.write(out)
  fh.close()
  out = "#sample\tclass\tJ gene\tfrequency\tuniq_read_freq\n"
  for c in j_counts:
    out=out+sample+"\t"+c+"\t"+str("NA")+"\t"+str(j_counts[c])+"\n"
  fh=open(V_gene_isotype_frequency_file_grouped.replace("V_gene","J_gene"),"w")
  fh.write(out)
  fh.close()
  return()

def Get_V_gene_isotype_frequency(sample, annot_file, V_gene_isotype_frequency_file,per_cluster_developmental_classification_file,V_gene_isotype_frequency_file_expanded):
  #Get_total_V_gene_usage(sample, annot_file, V_gene_isotype_frequency_file,per_cluster_developmental_classification_file)
  Get_expanded_V_gene_usage(sample, per_cluster_developmental_classification_file,V_gene_isotype_frequency_file_expanded)

def Get_expanded_V_gene_usage(sample, per_cluster_developmental_classification_file,V_gene_isotype_frequency_file_expanded):
  counts_expanded = Tree()
  threshold = 3 # number of BCRs per cluster to be considered expanded
  fh=open(per_cluster_developmental_classification_file,"r")
  out,ind = '#sample\tclass\tV gene\tJ gene\tfrequency\tmean_mutations\n',0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cluster_size=int(l[11])
      if(cluster_size>=threshold):
        id,clas,v = l[1],l[2].split("|"),l[4]
        for c in clas: 
          if(c.count(",")==0):
            counts_expanded[c][v][id].value = 1
  fh.close()
  out = '#sample\tclass\tV gene\tuniq_BCRs\n'
  for c in counts_expanded:
    for v in counts_expanded[c]:
      out = out+"\t".join(map(str,[sample, c, v, len(counts_expanded[c][v])]))+"\n"
  fh=open(V_gene_isotype_frequency_file_expanded,"w")
  fh.write(out)
  fh.close()
  return()

def Get_total_V_gene_usage(sample, annot_file, V_gene_isotype_frequency_file,per_cluster_developmental_classification_file):
  seqs,freqs,alias = Get_sequences(seq_file)
  CDR3s,mutations = Get_annotations(annot_file)
  counts,counts_expnaded = {},{}
  for id in alias:
    classification = alias[id].split("|")[1].split("_")
    break
  for id in mutations:
    vj,mut = mutations[id][3], mutations[id][0]
    freq = freqs[id]
    nz = [i for i in range(len(freq)) if freq[i]!=0]
    for i in nz:
      c = classification[i].split("*")[0]
      if(reverse_primer_group!="ISOTYPER"):
        if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
        elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
      c1 = c+"\t"+vj
      if(c1 in counts):counts[c1] = counts[c1]+[mut]*freq[i]
      else:counts[c1]=[mut]*freq[i]
      c1 = "ALL\t"+vj
      if(c1 in counts):counts[c1] = counts[c1]+[mut]*freq[i]
      else:counts[c1]=[mut]*freq[i]
  out = "#sample\tclass\tV gene\tJ gene\tfrequency\tmean_mutations\n"
  for c in counts:
    out=out+sample+"\t"+c+"\t"+str(len(counts[c]))+"\t"+str(mean(counts[c]))+"\n"
  fh=open(V_gene_isotype_frequency_file,"w")
  fh.write(out)
  fh.close()
  return()

def Isotype_sharing_subsampled(sample,reverse_primer_group, subsample_file,isotype_sharing_vertices_SUBSAMPLED,cluster_file,annot_file,isotype_sharing_SUBSAMPLED_summary):
  CDR3s,annots = Get_annotations(annot_file)
  Isotypes_shared_with_one_other_subsampled(sample,reverse_primer_group, subsample_file,isotype_sharing_vertices_SUBSAMPLED,cluster_file,annots,isotype_sharing_SUBSAMPLED_summary)
  #Per_cluster_analysis_isotype_sharing_subsampled(sample,reverse_primer_group, subsample_file,isotype_sharing_vertices_SUBSAMPLED,cluster_file,annot_file)

def Isotypes_shared_with_one_other_subsampled(sample,reverse_primer_group, subsample_file,isotype_sharing_vertices_SUBSAMPLED,cluster_file,annots,isotype_sharing_SUBSAMPLED_summary):
  fh=open(subsample_file,"r")
  freqs,clusters = {}, Tree()
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,freq, id_short,clust = l[2].split("__")[0],list(map(int, l[2].split("__")[1].split("|")[0].split("_"))),l[2].split("__")[0], l[1]
      clust,f = l[1],sum(freq)
      classification = l[2].split("|")[1].split("_")
      subsample = l[0]
      #if(subsample in total):total[subsample] = total[subsample]+f
      #else:total[subsample]=f
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      classes = {}
      for i in nz:
        c = classification[i].split("*")[0]
        if(reverse_primer_group!="ISOTYPER"):
          if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
          elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        if(c in ["IGHD", "IGHM"]):c = "IGHD/M"
        if(c in classes):classes[c] = classes[c]+freq[i]
        else:classes[c]=freq[i]
        clusters[subsample][id][c].value = 1
      for c in classes:
        n = id+"\t"+c
        freqs[n] = classes[c]
      #if(len(freqs)>10000):break ##################3
  fh.close()
  out="#Sample\tsubsample\tn_isotypes\tisotypes\tcounts\tmean_mutations(nn)\tmean_CDR3_region(aa)\tn_BCRs_counted_for_muts\n"
  counts,mutations,CDR3_lengths = {},{},{}
  for subsample in clusters:
    for id in clusters[subsample]:
      if(len(clusters[subsample][id])>0):
        if(id in annots):
          muts, CDR3 = annots[id][0],len(annots[id][2].split(":")[0])
        iso = []
        for c in clusters[subsample][id]:
          iso.append(c)
        iso.sort()
        name = subsample+"\t"+str(len(iso))+"\t"+",".join(iso)
        if( name in counts): counts[name] = counts[name]+1
        else:counts[name]= 1
        if(id in annots):
          muts, CDR3 = annots[id][0],len(annots[id][2].split(":")[0])
          if( name in mutations):mutations[name],CDR3_lengths[name] = mutations[name] +[muts], CDR3_lengths[name] +[CDR3]
          else:mutations[name],CDR3_lengths[name] =[muts],[CDR3]
  for n in counts:
    num, muts, CDR3s = -1,[-1],[-1]
    if(n in mutations):
      num, muts, CDR3s = len(mutations[n]), mean(mutations[n]), mean(CDR3_lengths[n])
      out=out+"\t".join(map(str,[sample, n, counts[n], mean(muts), mean(CDR3s),num ] ))+"\n"
  fh=open(isotype_sharing_SUBSAMPLED_summary,"w")
  fh.write(out)
  fh.close()
  return()

def Per_cluster_analysis_isotype_sharing_subsampled(sample,reverse_primer_group, subsample_file,isotype_sharing_vertices_SUBSAMPLED,cluster_file,annot_file):
  fh=open(subsample_file,"r")
  clusters,vertices,uniq_vertices,total_vertices = Tree(),{},{},{}
  total = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,freq, id_short,clust = l[2].split("__")[0],list(map(int, l[2].split("__")[1].split("|")[0].split("_"))),l[2].split("__")[0], l[1]
      clust,f = l[1],sum(freq)
      classification = l[2].split("|")[1].split("_")
      subsample = l[0]
      if(subsample in total):total[subsample] = total[subsample]+f
      else:total[subsample]=f
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      classes = []
      for i in nz:
        c = classification[i].split("*")[0]
        if(reverse_primer_group!="ISOTYPER"):
          if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
          elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        if(c in ["IGHD", "IGHM"]):c = "IGHD/M"
        if(c not in classes):classes = classes+[c]
        clusters[subsample][clust][c].value = 1
      c = subsample
      if(c in total_vertices):total_vertices[c] = total_vertices[c]+1
      else:total_vertices[c] = 1
      if(len(classes)>1):
        classes.sort()
        for c1 in range(0,len(classes)):
          for c2 in range(c1,len(classes)):
            if(c1!=c2):
              c = subsample+"\t"+classes[c1]+"-"+classes[c2]
              if(c in vertices): 
                vertices[c] = vertices[c]+f
                uniq_vertices[c] = uniq_vertices[c]+1
              else:
                vertices[c] = f
                uniq_vertices[c] =1 
  fh.close()
  fh=open(cluster_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,freq, id_short,clust = l[2].split("__")[0],list(map(int, l[2].split("__")[1].split("|")[0].split("_"))),l[2].split("__")[0], l[1]
      clust,f = l[1],sum(freq)
      classification = l[2].split("|")[1].split("_")
      subsample = "ALL"
      if(subsample in total):total[subsample] = total[subsample]+f
      else:total[subsample]=f
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      classes = []
      for i in nz:
        c = classification[i].split("*")[0]
        if(reverse_primer_group!="ISOTYPER"):
          if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
          elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        if(c in ["IGHD", "IGHM"]):c = "IGHD/M"
        if(c not in classes):classes = classes+[c]
        clusters[subsample][clust][c].value = 1
      c = subsample
      if(c in total_vertices):total_vertices[c] = total_vertices[c]+1
      else:total_vertices[c] = 1
      if(len(classes)>1):
        classes.sort()
        for c1 in range(0,len(classes)):
          for c2 in range(c1,len(classes)):
            if(c1!=c2):
              c = subsample+"\t"+classes[c1]+"-"+classes[c2]
              if(c in vertices): 
                vertices[c] = vertices[c]+f
                uniq_vertices[c] = uniq_vertices[c]+1
              else:
                vertices[c] = f
                uniq_vertices[c] =1
  cluster_sharing = {}
  for subsample in clusters:
    for clust in clusters[subsample]:
      if(len(clusters[subsample][clust])>1):
        classes = []
        for c in clusters[subsample][clust]:
          classes.append(c)
        classes.sort()
        for c1 in range(0,len(classes)):
          for c2 in range(c1,len(classes)):
            if(c1!=c2):
              c = subsample+"\t"+classes[c1]+"-"+classes[c2]
              if(c in cluster_sharing):cluster_sharing[c] = cluster_sharing[c]+1
              else:cluster_sharing[c]=1
  out = "#ID\tisotype\tvertex sharing frequency\tcluster_sharing_frequency\tuniq_vertex_sharing_frequency\ttotal reads\ttotal_vertices\tsample\n"
  for c in cluster_sharing:
    subsample = c.split("\t")[0]
    v,vu = 0,0
    if(c in vertices):v,vu = vertices[c],uniq_vertices[c]
    total_vertice = total_vertices[subsample]
    out = out+c+"\t"+str(v)+"\t"+str(cluster_sharing[c])+"\t"+str(vu)+"\t"+str(total[subsample])+"\t"+str(total_vertice)+"\t"+sample+"\n"
  fh=open(isotype_sharing_vertices_SUBSAMPLED,"w")
  fh.write(out)
  fh.close()
  return()

def Isotype_sharing_subsampled1(sample,reverse_primer_group, subsample_depth_file,isotype_sharing_vertices_SUBSAMPLED,cluster_file,isotype_sharing_SUBSAMPLED_summary):
  Get_repeated_raw_values(sample,reverse_primer_group, subsample_depth_file,isotype_sharing_vertices_SUBSAMPLED,cluster_file)
  #Get_summary_isotype_sharing_subsampled(isotype_sharing_SUBSAMPLED_summary,isotype_sharing_vertices_SUBSAMPLED)

def Get_summary_isotype_sharing_subsampled(isotype_sharing_SUBSAMPLED_summary,isotype_sharing_vertices_SUBSAMPLED):
  fh=open(isotype_sharing_vertices_SUBSAMPLED, "r")
  counts1, counts2, muts_overlap1, muts_overlap2,muts_non_overlap1,muts_non_overlap2 = {},{},{},{},{},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(int(l[5])>0 and int(l[6])>0):
        name = "\t".join(l[0:7])
        ovlap1, ovlap2,mut_in1, mut_in2, mutout1, mutout2 = int(l[7]),int(l[8]),float(l[9]),float(l[10]),float(l[11]),float(l[12])
        if(name in counts1):
          counts1[name] = counts1[name]+[ovlap1]
          counts2[name] = counts2[name]+[ovlap2]
        else:
          counts1[name]=[ovlap1]
          counts2[name]=[ovlap2]
        if(mutout1!=-1 and mutout2!=-1):
          if(name in muts_non_overlap1):
            muts_non_overlap1[name],muts_non_overlap2[name] = muts_non_overlap1[name]+[mutout1],muts_non_overlap2[name]+[mutout2]
          else:
            muts_non_overlap1[name],muts_non_overlap2[name] = [mutout1],[mutout2]
        if(mut_in1!=-1 and mut_in2!=-1):
          if(name in muts_overlap1):
            muts_overlap1[name], muts_overlap2[name] = muts_overlap1[name]+[mut_in1], muts_overlap2[name]+[mut_in2]
          else:
            muts_overlap1[name], muts_overlap2[name] = [mut_in1],[mut_in2]
  fh.close()
  out="#sample\tisotype1\tisotype2\ttotal1\ttotal2\tdepth1\tdepth2\toverlap1\toverlap2\tmean_muts_overlap1\tmean_muts_overlap2\tmean_muts_no_overlap1\tmean_muts_no_overlap2\n"
  for name in counts1:
    mut_in1, mut_in2, mutout1, mutout2 = -1,-1,-1,-1
    if(name in muts_overlap1): 
      if(len(muts_overlap1[name])>=5):
        mut_in1, mut_in2 = mean(muts_overlap1[name]), mean(muts_overlap2[name])
    if(name in muts_non_overlap1):
      if(len(muts_non_overlap1[name])>=5):
        mutout1, mutout2 = mean(muts_non_overlap1[name]), mean(muts_non_overlap2[name])
    values = [mean(counts1[name]), mean(counts2[name]), mut_in1, mut_in2, mutout1, mutout2]
    out=out+name+"\t"+"\t".join(map(str, values))+"\n"
  fh=open(isotype_sharing_SUBSAMPLED_summary,"w")
  fh.write(out)
  fh.close()
  return()

def Get_repeated_raw_values(sample,reverse_primer_group, subsample_depth_file,isotype_sharing_vertices_SUBSAMPLED,cluster_file):
  CDR3s,mutations = Get_annotations(annot_file)
  fh=open(subsample_depth_file,"r")
  depth = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      depth[l[0]] = int(l[1])
  fh.close()
  fh=open(cluster_file,"r")
  class_IDs, clusters = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      freq, id_short,clust =list(map(int, l[2].split("__")[1].split("|")[0].split("_"))),l[2].split("__")[0], int(l[1])
      classification = l[2].split("|")[1].split("_")
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      for i in nz:
        c = classification[i].split("*")[0]
        if(reverse_primer_group!="ISOTYPER"):
          if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
          elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        if(c in ["IGHD", "IGHM"]):c = "IGHD/M"
        if(c.count("IGHD")!=0):c = "IGHD/M"
        if(c.count("IGHM")!=0):c = "IGHD/M"
        if(c in class_IDs):class_IDs[c] = class_IDs[c]+[id_short]*freq[i]
        else:class_IDs[c]=[id_short]*freq[i]
        # if(len(class_IDs[c])>1000):break ########################
      clusters[id_short] = clust
  fh.close()
  classes = []
  isotype_sharing_vertices_UNSAMPLED = isotype_sharing_vertices_SUBSAMPLED.replace("SUBSAMPLED","TOTAL")
  out = "Chain\tnumber of reads\tthreshold reads for sampling\tpass/fail\n"
  for c in class_IDs:
    if(c in depth):
      if(len(class_IDs[c])>=depth[c]):
        classes.append(c)
        out=out+c+"\t"+str(len(class_IDs[c]))+"\t"+str(depth[c])+"\tPASS\n"
      else:
        out=out+c+"\t"+str(len(class_IDs[c]))+"\t"+str(depth[c])+"\tFAIL: NOT ENOUGH READS\n"
    else:
      out=out+c+"\t"+str(len(class_IDs[c]))+"\tNA\tFAIL: NOT DEFINED IN REFERENCE FILE\n"
  classes.sort()
  print(out)
  print(classes)
  out="#sample\tisotype1\tisotype2\ttotal1\ttotal2\tdepth1\tdepth2\toverlap1\toverlap2\tmean_muts_overlap1\tmean_muts_overlap2\tmean_muts_no_overlap1\tmean_muts_no_overlap2\n"
  out1, ind ='', 0
  for f in [isotype_sharing_vertices_SUBSAMPLED,isotype_sharing_vertices_UNSAMPLED]:
    fh=open(f,"w")
    fh.write(out)
    fh.close()
  if(len(classes)>1):
    repeats = 1000
    for c1 in range(0,len(classes)):
      for c2 in range(c1,len(classes)):
        if(c1!=c2):
          ######## unsampled
          depth1, depth2 = len(class_IDs[classes[c1]]), len(class_IDs[classes[c2]])
          clusts = [[],[]]
          muts = [[],[]]
          sample1,sample2 = class_IDs[classes[c1]], class_IDs[classes[c2]]
          for i in [0,1]:
            X = [sample1,sample2][i]
            clustx,mutx = [],[]
            for id in X:
              clustx.append(clusters[id])
              if(id in mutations):
                mutx.append(mutations[id][0])
              else:mutx.append(-1)
            clusts[i] = clustx
            muts[i] = mutx
          intersection = set(clusts[0]) & set(clusts[1])
          muts_overlap1,muts_overlap2,muts_not_overlap1,muts_not_overlap2 = [],[],[],[]
          overlap1, overlap2 = 0,0
          for i in range(0,depth1):
            if(clusts[0][i] in intersection):
              overlap1 = overlap1+1
              if(muts[0][i] !=-1):muts_overlap1.append(muts[0][i])
            else:
              if(muts[0][i] !=-1):muts_not_overlap1.append(muts[0][i])
          for i in range(0,depth2):
            if(clusts[1][i] in intersection):
              overlap2 = overlap2+1
              if(muts[1][i] !=-1):muts_overlap2.append(muts[1][i])
            else:
              if(muts[1][i] !=-1):muts_not_overlap2.append(muts[1][i])
          muts_array = [muts_overlap1,muts_overlap2,muts_not_overlap1,muts_not_overlap2]
          means = []
          for i in muts_array: 
            if(len(i)>0):mean1 = mean(i)
            else:mean1=-1
            means.append(mean1)
          out1=out1+sample+"\t"+classes[c1]+"\t"+classes[c2]+"\t"+str(len(class_IDs[classes[c1]]))+"\t"+str(len(class_IDs[classes[c2]]))+"\t"+str(depth1)+"\t"+str(depth2)+"\t"+str(overlap1)+"\t"+str(overlap2)+"\t"+"\t".join(map(str, means))+"\n"
        if(c1!=c2 and 1==2):
          for r in range(0,repeats):
            depth1, depth2 = depth[classes[c1]], depth[classes[c2]]
            sample1 = random.sample(class_IDs[classes[c1]], depth1)
            sample2 = random.sample(class_IDs[classes[c2]], depth2)
            clusts = [[],[]]
            muts = [[],[]]
            for i in [0,1]:
              X = [sample1,sample2][i]
              clustx,mutx = [],[]
              for id in X:
                clustx.append(clusters[id])
                if(id in mutations):
                  mutx.append(mutations[id][0])
                else:mutx.append(-1)
              clusts[i] = clustx
              muts[i] = mutx
            intersection = set(clusts[0]) & set(clusts[1])
            muts_overlap1,muts_overlap2,muts_not_overlap1,muts_not_overlap2 = [],[],[],[]
            overlap1, overlap2 = 0,0
            for i in range(0,depth1):
              if(clusts[0][i] in intersection):
                overlap1 = overlap1+1
                if(muts[0][i] !=-1):muts_overlap1.append(muts[0][i])
              else:
                if(muts[0][i] !=-1):muts_not_overlap1.append(muts[0][i])
            for i in range(0,depth2):
              if(clusts[1][i] in intersection):
                overlap2 = overlap2+1
                if(muts[1][i] !=-1):muts_overlap2.append(muts[1][i])
              else:
                if(muts[1][i] !=-1):muts_not_overlap2.append(muts[1][i])
            muts_array = [muts_overlap1,muts_overlap2,muts_not_overlap1,muts_not_overlap2]
            means = []
            for i in muts_array:
              if(len(i)>0):mean1 = mean(i)
              else:mean1=-1
              means.append(mean1)
            out=out+sample+"\t"+classes[c1]+"\t"+classes[c2]+"\t"+str(len(class_IDs[classes[c1]]))+"\t"+str(len(class_IDs[classes[c2]]))+"\t"+str(depth1)+"\t"+str(depth2)+"\t"+str(overlap1)+"\t"+str(overlap2)+"\t"+"\t".join(map(str, means))+"\n"
            ind = ind+1
            if(ind>500):
              Write_out(out,isotype_sharing_vertices_SUBSAMPLED)
              out, ind = '',0
    Write_out(out,isotype_sharing_vertices_SUBSAMPLED)
  Write_out(out1,isotype_sharing_vertices_UNSAMPLED)
  return()

def Get_network_parameters_per_classificiation(sample, cluster_file, per_cluster_developmental_classification_file, per_cluster_developmental_classification_network_parameters,reverse_primer_group):
  fh = open(per_cluster_developmental_classification_file,"r")
  classification_ids,classification_clusters = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,classes,clust = l[1],l[2].split("|"),int(l[10])
      classification_ids[id] = classes
      for c in classes:
        if(c in classification_clusters):
          if(clust not in classification_clusters[c]):classification_clusters[c] = classification_clusters[c]+[clust]
        else:classification_clusters[c]=[clust]
  fh.close()
  fh=open(cluster_file,"r")
  v_sizes,c_sizes = {},{}
  index = 0
  for l in fh:
    index=index+1
    if (index>1):
      l=l.strip().split()
      id,freq, id_short,clust = l[2],sum(map(int, l[2].split("__")[1].split("|")[0].split("_"))),l[2].split("__")[0], int(l[1])
      if(id_short in classification_ids):
        for c in classification_ids[id_short]:
          if(c in v_sizes):
            v_sizes[c] = v_sizes[c]+[freq]
          else:v_sizes[c]=[freq]
      if(clust in c_sizes):c_sizes[clust] = c_sizes[clust]+freq
      else:c_sizes[clust]=freq
  fh.close()
  out="#Id\tIsotype\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\n"
  for c in v_sizes:
    c_sizes_classification = []
    for clust in classification_clusters[c]:
      c_sizes_classification = c_sizes_classification+[c_sizes[clust]]
    (vpoints,vvdf)=VDF(v_sizes[c])
    (cpoints,cvdf)=VDF(c_sizes_classification)
    vgini, cgini=Gini_index(cpoints,cvdf, vpoints,vvdf)
    max_pop, max_1_pop = cpoints[len(cpoints)-1]*100.0/sum(v_sizes[c]), cpoints[len(cpoints)-2]*100.0/sum(v_sizes[c])
    out = out+str(sample)+"\t"+c+"\t"+str(sum(v_sizes[c]))+"\t"+str(len(v_sizes[c]))+"\t"+str(vgini)+"\t"+str(cgini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\n"
  fh=open(per_cluster_developmental_classification_network_parameters, "w")
  fh.write(out)
  fh.close()
  return()

def Uniq(v):
  C=set(v)
  return list(C)

def VDF (n):
  points=sorted(Uniq(n))
  vdf=[]
  for i in range(0,len(points)):
    vdf.append(n.count(points[i]))
  return (points,vdf)

def Get_distribution_mutational_differences(tmp_file1,len_cluster):
  aligned = tmp_file1+".Aln"
  insert = ''
  if(len_cluster > 1000):insert = '--parttree'
  command0 = "/software/pubseq/bin/mafft-6.857/bin/mafft --retree 2 "+insert+" "+tmp_file1+" > "+aligned
  subprocess.getoutput(command0)
  fh=open(aligned,"r")
  sub_seqs =[]
  for header,sequence in fasta_iterator(fh):
    sub_seqs.append(sequence)
  fh.close()
  dist = {}
  t = 0
  for i1 in range(0,len(sub_seqs)):
    for j in range(i1,len(sub_seqs)):
      if(i1<j):
        t = t+1
        s1,s2 = sub_seqs[i1],sub_seqs[j]
        diffs = len([i for i in range(len(s1)) if s1[i]!=s2[i] and s1[i]!="-" and s2[i]!="-"])
        if(diffs in dist):dist[diffs] = dist[diffs]+1
        else:dist[diffs] = 1
  return(dist)

def Get_mutational_distance_distribution_per_cluster(seq_file, per_cluster_developmental_classification_file, mutational_distance_distribution_per_cluster,tmp_file1,sample):
  print(per_cluster_developmental_classification_file)
  fh = open(per_cluster_developmental_classification_file,"r")
  clusters,cluster_info = Tree(),{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,cluster = l[1],l[10]
      clusters[cluster][id].value = 1
      cluster_info[cluster] = l[2]
  fh.close()
  ids,seqs = {},{}
  fh=open(seq_file,"r")
  for header,sequence in fasta_iterator(fh):
    seqs[header.split("__")[0]] = sequence
  fh.close()
  print("N IDs:",len(ids),"N seqs:", len(seqs))
  fh=open(mutational_distance_distribution_per_cluster,"w")
  fh.close()
  out,ind= "#sample\tcluster_ID\tcluster_size\tcluster_definition\tn_differences\tfrequency\n",0
  for c in clusters:
    if(len(clusters[c])>20 and len(clusters[c])<1000):
      out1 = ''
      for id in clusters[c]:
        out1 = out1+">"+id+"\n"+seqs[id]+"\n"
      fh=open(tmp_file1,"w")
      fh.write(out1)
      fh.close()
      dist = Get_distribution_mutational_differences(tmp_file1,len(clusters[c]))
      clust_info = str(len(clusters[c]))+"\t"+cluster_info[c]
      for d in dist:
        out=out+sample+"\t"+c+"\t"+clust_info+"\t"+str(d)+"\t"+str(dist[d])+"\n"
        ind = ind+1
        if(ind>200):
          Write_out(out, mutational_distance_distribution_per_cluster)
          out, ind= '',0
  Write_out(out, mutational_distance_distribution_per_cluster)
  out, ind= '',0
  print(mutational_distance_distribution_per_cluster)
  return()

def Isotype_per_cluster_sharing_probabilities(seq_file,per_cluster_developmental_classification_file, per_cluster_inter_isotype_sharing,sample,reverse_primer_group):
  seqs,freqs,alias = Get_sequences(seq_file)
  fh = open(per_cluster_developmental_classification_file,"r")
  total = 0
  class_count = {}
  cluster_sizes = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id = l[1]
      f = sum(freqs[id])
      classes,total = l[3].split(","),total+f
      cluster_id = l[10]
      #if(cluster_id in cluster_count):
      #  cluster_count[cluster_id] = cluster_count[cluster_id]+f
      #else:cluster_count[cluster_id]=f
      classes.sort()
      cluster_sizes[cluster_id] = classes
  fh.close()
  isotype_sharing = {}
  for c in cluster_sizes:
    classes = cluster_sizes[c]
    if(len(classes)>1):
      for i in range(0,len(classes)):
        for j in range(i,len(classes)):
          if(i<j):
            name = classes[i]+"-"+classes[j]
            if(name.count("P")==0):
              if(name in isotype_sharing):isotype_sharing[name] = isotype_sharing[name]+1
              else: isotype_sharing[name]=1
    else:
      name = classes[0]
      if(name.count("P")==0):
        if(name in isotype_sharing):isotype_sharing[name] = isotype_sharing[name]+1
        else: isotype_sharing[name]=1
  out = "#sample\tcluster_class_overlap\tfrequency\ttotal_clusters\n"
  n_clusters = len(cluster_sizes)
  for c in isotype_sharing:
    if(reverse_primer_group !="ISOTYPER"):
      out=out+sample+"\t"+c.replace("IGHA","IGHA1/2")+"\t"+str(isotype_sharing[c])+"\t"+str(n_clusters)+"\n"
    else:
      out=out+sample+"\t"+c+"\t"+str(isotype_sharing[c])+"\t"+str(n_clusters)+"\n"
  fh=open(per_cluster_inter_isotype_sharing,"w")
  fh.write(out)
  fh.close()
  return()

def Get_cluster_properties_per_isotype(cluster_properties_per_isotype_SUBSAMPLED, sample, per_cluster_developmental_classification_file,subsample_file,cluster_size_file):
  fh = open(subsample_file,"r")
  clusters, id_types = Tree(),{}
  overall_cluster_sizes = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      subsample, cluster,freq,classification = l[0],l[1],list(map(int, l[2].split("__")[1].split("|")[0].split("_"))), l[2].split("|")[1].split("_")
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      classes = {}
      for i in nz:
        c = classification[i].split("*")[0]
        if(c in ["IGHA1", "IGHA2"]):c = "IGHA1/2"
        if(c in ["IGHG1", "IGHG2"]):c = "IGHG1/2"
        if(c.count("P")==0):
          if(c in classes):classes[c] = classes[c]+freq[i]
          else:classes[c]=freq[i]
      c1 = subsample+"\t"+cluster
      if(c1 in overall_cluster_sizes):overall_cluster_sizes[c1] = overall_cluster_sizes[c1]+1
      else:overall_cluster_sizes[c1] = 1
      for c in classes:
        clusters[subsample][cluster][c].value = 1
        c1 = subsample+"\t"+cluster+"\t"+c
        if(c1 in id_types):id_types[c1] = id_types[c1]+classes[c]
        else:id_types[c1]=classes[c]
  fh.close()
  cluster_size_dist = {}
  for c in overall_cluster_sizes:
    subsample = c.split("\t")[0]
    if(subsample in cluster_size_dist):cluster_size_dist[subsample] = cluster_size_dist[subsample]+[overall_cluster_sizes[c]]
    else:cluster_size_dist[subsample]=[overall_cluster_sizes[c]]
  out = "#sample\tsubsample_repeats\tcluster_size\tfrequency\n"
  mean_freq = {}
  for subsample in cluster_size_dist:
    l = cluster_size_dist[subsample]
    list_freq = [[x,l.count(x)] for x in set(l)]
    for i in list_freq:
      if(i[0] in mean_freq):mean_freq[i[0]] = mean_freq[i[0]]+[i[1]]
      else:mean_freq[i[0]]=[i[1]]
  repeats = len(cluster_size_dist)
  for size in mean_freq:
    freq = sum(mean_freq[size])*1.0/repeats
    out = out+sample+"\t"+str(repeats)+"\t"+str(size)+"\t"+str(freq)+"\n"
  fh=open(cluster_size_file,"w")
  fh.write(out)
  fh.close()


def dfgdfg():
  ### For all clusters with isotype X, proportion of Y: Y_given_X
  Y_given_X, All_X= {},{}
  for subsample in clusters:
    for cluster in clusters[subsample]:
      array,classes = [],[]
      for c in clusters[subsample][cluster]:
        array,classes = array+[id_types[subsample+"\t"+cluster+"\t"+c]],classes+[c]
      array = np.array(array)
      cluster_size = sum(array)
      array = array*1.0/sum(array)
      for c1 in range(0,len(array)):
        c = subsample+"\t"+classes[c1]
        if(c in All_X):All_X[c] = All_X[c]+1
        else:All_X[c]=1
        for c2 in range(0,len(array)):
            c = subsample+"\t"+classes[c1]+"\t"+classes[c2] # X given Y
            if(c in Y_given_X):Y_given_X[c] = Y_given_X[c]+[array[c1]]
            else:Y_given_X[c]=[array[c1]]
  del id_types, clusters
  out="#sample\tsubsample\tX\tX:Y\tmean(Y given X)\tnumber of clusters\tnumber of clusters with isotype Y\tnumber of clusters with isotype X\n"
  for c in Y_given_X:
    c1 = c.split("\t")
    c2 = c1[0]+"\t"+c1[1]
    c3 = c1[0]+"\t"+c1[2]
    out=out+sample+"\t"+c1[0]+"\t"+c1[2]+"\t"+c1[2]+":"+c1[1]+"\t"+"\t".join(map(str, [mean(Y_given_X[c]), len(Y_given_X[c]), All_X[c2],All_X[c3]]))+"\n"
  fh=open(cluster_properties_per_isotype_SUBSAMPLED,"w")
  fh.write(out)
  fh.close()
  return()


def Proportion_clusters_isotype_overlapping(sample, reverse_primer_group, subsample_file, clusters_isotype_overlapping_SUBSAMPLED):
  fh = open(subsample_file,"r")
  clusters,isotypes = {},Tree()
  total = 0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      subsample, cluster, freq, classification,id = l[0],l[1],list(map(int, l[2].split("__")[1].split("|")[0].split("_"))), l[2].split("|")[1].split("_"), l[2].split("__")[0]
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      total = total+sum(freq)
      for i in nz:
        c = classification[i].split("*")[0]
        if(c in ["IGHD","IGHM"]):c = "IGHD/M"
        if(reverse_primer_group!="ISOTYPER"):
          if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
          if(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        isotypes[subsample][cluster][c].value = 1
      cluster = subsample+":"+cluster
      if(cluster in clusters):clusters[cluster] = clusters[cluster]+sum(freq)
      else:clusters[cluster]=sum(freq)
  fh.close()
  n_isotypes = {}
  for subsample in isotypes:
    for c in isotypes[subsample]:
      n = len(isotypes[subsample][c])
      freq = clusters[subsample+":"+c]
      n = subsample+"\t"+str(n)
      if(n in n_isotypes):n_isotypes[n] = n_isotypes[n]+freq
      else:n_isotypes[n]=freq
  out,ind="#sample\tsubsample\tnumber of isotypes\tnumber of reads\ttotal\n",0
  for n in n_isotypes:
    out=out+sample+"\t"+n+"\t"+str(n_isotypes[n])+"\t"+str(total)+"\n"
  fh=open(clusters_isotype_overlapping_SUBSAMPLED,"w")
  fh.write(out)
  fh.close()
  return()

def Get_proportion_read_mutiple_subtypes(seq_file,per_cluster_developmental_classification_file, shared_isotype_counts,sample,reverse_primer_group):
  seqs,freqs,alias = Get_sequences(seq_file)
  fh = open(per_cluster_developmental_classification_file,"r")
  total = 0
  class_count = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id = l[1]
      f = sum(freqs[id])
      classes,total = l[2].split("|"),total+f
      for c in classes:
        if(c.count("P")==0):
          if(c in class_count):class_count[c] = class_count[c]+f
          else:class_count[c]=f
  fh.close()
  out = "#sample\tisotype_definition\treads\ttotal reads\tpercentage\n"
  for c in class_count:
    out=out+sample+"\t"+c+"\t"+str(class_count[c])+"\t"+str(total)+"\t"+str(class_count[c]*100.0/total)+"\n"
  fh=open(shared_isotype_counts,"w")
  fh.write(out)
  fh.close()
  return()

def Get_CDR3_lengths_charges_per_classification(sample,seq_unmutated_file, annot_file, developmental_classification_file,developmental_classification_CDR3_lengths_charges,annot_file_internal,reverse_primer_group):
  Get_charges_IMGT_annotation(sample, annot_file, developmental_classification_file,developmental_classification_CDR3_lengths_charges,seq_unmutated_file)
  #Get_charges_internal_annotation(sample,annot_file_internal, developmental_classification_file,developmental_classification_CDR3_lengths_charges)

def Get_charges_internal_annotation(sample,annot_file_internal, developmental_classification_file,developmental_classification_CDR3_lengths_charges):
  fh=open(annot_file_internal,"r")
  cdr3s= {}
  CDR3, CDR3_VJ = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(len(l)>17):
        id = l[0].split("__")[0]
        cdr3,v,j= l[17],l[1].split("*")[0],l[13].split("*")[0]
        if(cdr3 != 'CDR3_UNDETERMINED' and len(cdr3)>5):
          CDR3[id] = len(cdr3)
          v_fam = v[0:5]
          v = v_fam
          CDR3_VJ[id] = str(len(cdr3))+"\t"+v+"\t"+j+"\t"+v+":"+j+"\t"+v_fam+"\t"+v_fam+":"+j
  fh.close()
  print(len(CDR3_VJ) , len(CDR3))
  fh = open(developmental_classification_file,"r") 
  CDR3_lengths,CDR3_charge,total,clusters_done = {},{},0,{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id = l[1]
      if(id in CDR3_VJ and l[10] not in clusters_done):
        clusters_done[l[10]] = 1 
        classes,total = l[2].split("|"),total+1
        length =CDR3_VJ[id]
        for c in classes:
          if(c.count("isotype")==0):
            name_charge = c+"\t"+str(length)
            if(name_charge in CDR3_charge):CDR3_charge[name_charge] = CDR3_charge[name_charge]+1
            else:CDR3_charge[name_charge] = 1
  fh.close()
  out = "#ID\tclassification\tCDR3 length\tV\tJ\tVJ\tVfam\tVfamJ\tnumber of clusters\ttotal number of clusters\n"
  for c in CDR3_charge:
    out=out+sample+"\t"+c+"\t"+str(CDR3_charge[c])+"\t"+str(total)+"\n"
  fh=open(developmental_classification_CDR3_lengths_charges,"w")
  fh.write(out)
  fh.close()
  return()

def Get_charges_IMGT_annotation(sample, annot_file, developmental_classification_file,developmental_classification_CDR3_lengths_charges,seq_file):
  #c = commands.getoutput("ls -l "+annot_file).split()
  #if(c[0]!='ls:' and len(c[4])>1):
  #  print "PASS", c
  if(1==1):
    pka = Get_residue_positivity()
    seqs,freqs,alias = Get_sequences(seq_file)
    cdr3s = Get_CDR3s(annot_file)
    fh = open(developmental_classification_file,"r") 
    CDR3_lengths,CDR3_charge,total,clusters_done = {},{},0,{}
    for l in fh:
      if(l[0]!="#"):
        l=l.strip().split()
        if(l[1] in cdr3s and l[10] not in clusters_done and l[1] in seqs):
          clusters_done[l[10]] = 1
          classes,total = l[2].split("|"),total+1
          length =len(cdr3s[l[1]])
          pos = Get_positivity(cdr3s[l[1]],pka)
          f = sum(freqs[l[1]])
          for c in classes:
            name_charge = c+"\t"+str(length)+"\t"+str(pos)
            if(name_charge in CDR3_charge):CDR3_charge[name_charge] = CDR3_charge[name_charge]+f
            else:CDR3_charge[name_charge] = f
            #name_length = c+"\t"+str(length)
            #if(name_length in CDR3_lengths):CDR3_lengths[name_length] = CDR3_lengths[name_length] +1
            #else:CDR3_lengths[name_length] = 1
    out = "#ID\tclassification\tCDR3 length\tCDR3 charge\tnumber of clusters\ttotal number of clusters\n"
    for c in CDR3_charge:
      out=out+sample+"\t"+c+"\t"+str(CDR3_charge[c])+"\t"+str(total)+"\n"
    fh=open(developmental_classification_CDR3_lengths_charges,"w")
    fh.write(out)
    fh.close()
  else:
    print("\tFAIL",c)
  return()

def Get_secondary_rearrangement_potential(sample, annot_file,annot_file3, secondary_rearrangement_file,seq_file,raw_secondary_rearrangement_file, secondary_rearrangement_cdr3_lengths,mapped_secondary_rearrangements,tmp_file1,v_replacement_transition_counts,secondary_rearrangement_file_SAMPLED,cluster_file,reverse_primer_group):
  Get_secondary_rearrangement_sequences(sample, annot_file,annot_file3, secondary_rearrangement_file,seq_file,raw_secondary_rearrangement_file,reverse_primer_group)
  sample_depth = 250 ## for each FCGR2B sample
  sample_depth = 250 ## for each MURab sample
  Subsample_secondary_rearrangments(raw_secondary_rearrangement_file, sample_depth, secondary_rearrangement_file_SAMPLED,cluster_file,seq_file,sample,annot_file)
  #Subsample_per_Vgene(raw_secondary_rearrangement_file, sample_depth, secondary_rearrangement_file_SAMPLED,cluster_file,seq_file,sample,annot_file)
  #Secondary_rearrangement_CDR3_lengths(sample, annot_file,raw_secondary_rearrangement_file,secondary_rearrangement_cdr3_lengths,seq_file,cluster_file)
  #Map_secondary_rearrangement_joins(seq_file,sample, annot_file3, raw_secondary_rearrangement_file, mapped_secondary_rearrangements,tmp_file1,v_replacement_transition_counts)

def Map_secondary_rearrangement_joins(seq_file,sample, annot_file3, raw_secondary_rearrangement_file, mapped_secondary_rearrangements,tmp_file1,v_replacement_transition_counts):
  V_families = Get_reference_sequence_families()
  for f in [v_replacement_transition_counts,mapped_secondary_rearrangements]:
    fh=open(f,"w")
    fh.close()
  seqs,freqs,alias = Get_sequences(seq_file)
  fh=open(annot_file3,"r")
  regions = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split("\t")
      if(l[2]=="productive"):
        if(l[1].split("__")[0] in freqs):
          v_region, nd_region, j_region = l[8],l[18],l[40]
          v= l[3].split()[1].split("*")[0]
          f = sum(freqs[l[1].split("__")[0]])
          n1,n2 = l[21],l[28]
          regions[l[1].split("__")[0]] = [v, v_region, nd_region, j_region,f,n1,n2]
  fh.close()
  clustered,cluster_number = Tree(),1
  fh=open(raw_secondary_rearrangement_file, "r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[0] in regions):
        v = regions[l[0]][0]
        f = regions[l[0]][4]
        v1 = V_families[v]
        clustered[l[1]][v1][v][f][l[0]].value = 1
  fh.close()
  cluster_number = 1
  out_all,ind = '#sample\tcluster_id\tv_gene\tsequence\n',0
  transition = {}
  for c in clustered: 
    if(len(clustered[c])>1):
      ids,vs = [],[]
      vs_all = []
      for v1 in clustered[c]:
        ids_sub,max_freq,v_max = '',0,''
        for v in clustered[c][v1]:
          vs_all.append(v)
          count = 0
          for f in clustered[c][v1][v]:
            if(count>5):break
            if(max_freq<f):
              for id in clustered[c][v1][v][f]:
                ids_sub,max_freq,v_max = id,f,v
                #ids.append(id)
                #vs.append(v)
                count =count+1
        ids.append(ids_sub)
        vs.append(v_max) 
      #full_sequences,seq_lengths = [],[]
      #v_len, ndn_len,j_len = [],[],[]
      #for id in ids:
      #  full_sequences.append(seqs[id])
        #seq_lengths.append(len(seqs[id]))
        #v_len.append(len(regions[id][1]))
        #ndn_len.append(len(regions[id][2]))
        #j_len.append(len(regions[id][3]))
      #print v_len, ndn_len,j_len
      if(len(vs_all)>1):
        trans = ",".join(vs_all)
        if(trans in transition):transition[trans] = transition[trans]+1
        else:transition[trans] = 1
      if(len(ids)>1 and len(clustered[c])>2):
        out = ''
        for i in range(0,len(ids)):
          out=out+">"+ ids[i]+"\n"+seqs[ids[i]]+"\n"
        aln_seqs,region_starts = Get_aligned_sequences(out, tmp_file1,regions)
        if(len(region_starts)==3):
          annot_out= sample+"\t"+str(cluster_number)+"\tV_gene\t"
          #FULL ANNOTATION
          annot_out = annot_out+"V"*region_starts["N1"][0]+"N"*(region_starts["N1"][1]-region_starts["N1"][0])
          annot_out = annot_out+"D"*(region_starts["N2"][0]-region_starts["N1"][1])
          annot_out = annot_out+"M"*(region_starts["N2"][1]-region_starts["N2"][0])
          annot_out = annot_out+"J"*(region_starts["J"][1]-region_starts["J"][0])+"\n"
          for i in range(0,len(ids)): 
            annot_out= annot_out+sample+"\t"+str(cluster_number)+"\t"+vs[i]+"\t"+aln_seqs[ids[i]]+"\n"
          ind = ind+1
          out_all = out_all+annot_out
          print(annot_out)
          print(cluster_number)
          if(ind>50):
            Write_out(out_all, mapped_secondary_rearrangements)
            out_all,ind = '',0
          cluster_number = cluster_number+1
        #else:
          #BRIEF ANNOTATION
  Write_out(out_all, mapped_secondary_rearrangements)
  out = "#sample\ttransition\tfrequency\n"
  for trans in transition:
    out=out+sample+"\t"+trans+"\t"+str(transition[trans])+"\n"
  print(out)
  fh=open(v_replacement_transition_counts,"w")
  fh.write(out)
  fh.close()
  return()

def Get_aligned_sequences(out, tmp_file1,regions):
  fh=open(tmp_file1,"w")
  fh.write(out)
  fh.close()
  insert = ''
  command1 = "/software/pubseq/bin/mafft-6.857/bin/mafft --retree 2 "+insert+" "+tmp_file1+" > "+tmp_file1+".aligned"
  subprocess.getoutput(command1)
  fh=open(tmp_file1+".aligned","r")
  aln_seqs = {}
  region_starts = {}
  for header,sequence in fasta_iterator(fh):
    aln_seqs[header] = sequence.upper()
    n1,n2,j_region = regions[header][5],regions[header][6],regions[header][3]
    names, seq_reg = ["N1","N2","J"],[n1,n2,j_region]
    for i in range(0,len(names)):
      if(names[i] not in region_starts):
        if(sequence.count(seq_reg[i])==1):
          pos = sequence.index(seq_reg[i])
          region_starts[names[i]] = [pos,pos+len(seq_reg[i])]
  fh.close()
  return(aln_seqs,region_starts)

def Get_chromosomal_locations():
  file = "LIBRARY/18_06_25_Chromosomal_locations_IGHV_Ensembl.txt"
  fh=open(file,"r")
  #gene_loc, pos_loc = [],[]
  pos_loc = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      fail = 0
      if(len(l)>14):
        if(l[14]=="No"):
          fail = 1
      if(fail==0):
        start, v = int(l[3]), l[13]
        pos_loc[v] = start
        #gene_loc, pos_loc =gene_loc+[v], pos_loc +[start]
  fh.close()
  #gene_loc, pos_loc = np.array(gene_loc), np.array(pos_loc)
  return(pos_loc)

def Secondary_rearrangement_CDR3_lengths(sample, annot_file,raw_secondary_rearrangement_file,secondary_rearrangement_cdr3_lengths,seq_file,cluster_file):
  pos_loc = Get_chromosomal_locations()
  CDR3s,mutations1 = Get_annotations(annot_file)
  charged, polar = Get_charged_polar_positions()
  seqs,freqs,alias = Get_sequences(seq_file)
  fh=open(cluster_file, "r")
  cluster_sizes, cluster_ids,n_uniq_BCRs = {},{},0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      clust, id = l[1], l[2].split("__")[0]
      cluster_ids[id] =clust
      if(clust in cluster_sizes):cluster_sizes[clust] = cluster_sizes[clust]+1
      else:cluster_sizes[clust] = 1
      n_uniq_BCRs = n_uniq_BCRs+1
  fh.close()
  cdr3s = Get_CDR3s(annot_file)
  ids_rearranged,stems,longest_stems = {},{},{}
  fh=open(raw_secondary_rearrangement_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[0] in mutations1):
        ids_rearranged[l[0].split("__")[0]] = l[1]
        v_gene = mutations1[l[0]][3].split("\t")[0]
        if(v_gene in pos_loc):
          if( l[1] in stems):
            if(v_gene not in  stems[l[1]]):
              stems[l[1]] = stems[l[1]]+[v_gene]
          else:stems[l[1]]=[v_gene]
  fh.close()
  clostest_to_DJ_per_group = {}
  for stem in stems:
    if(len(stems[stem])>1):
      a = np.array(sort(stems[stem]))
      #higher chromosomal location = furtherst from IGHD/J
      locations = [pos_loc[a[i]] for i in range(len(a))]
      closest_to_DJ = a[np.where(np.array(locations)==min(locations))[0]]
      clostest_to_DJ_per_group[stem] = closest_to_DJ
  V_replaced_charged_polar_length = {}
  SHM,CDR3_len,CDR3_polar,CDR3_charge,R_residues,V_residues,cluster_sizes_per_group = {},{},{},{},{},{},{}
  overall_CDR3_length = {}
  for id in seqs:
    if(id in mutations1 and id  in cdr3s):
      freq, chains = freqs[id],alias[id].split("|")[1].split("_")
      mut = mutations1[id][0]+mutations1[id][1]
      v_gene = mutations1[id][3].split("\t")[0]
      cdr3 = cdr3s[id]
      length = len(cdr3)
      cdr3 = cdr3[0:8]
      n_polar,n_charged = 0,0
      for aa in charged:
        n_charged = n_charged+cdr3.count(aa)
      for aa in polar:
        n_polar = n_polar+cdr3.count(aa)
      if(id.split("__")[0] in ids_rearranged):
        stem= ids_rearranged[id.split("__")[0]]
        subtype = "NA"
        if(stem in clostest_to_DJ_per_group):
          if(v_gene==clostest_to_DJ_per_group[stem]):
            subtype = "Voriginal"
          else:subtype = "Vreplaced"
        typ = "REPLACED\t"+subtype
      else:typ = "NOT_REPLACED\tNone"
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      classes = {}
      for i in nz:
        c = chains[i].split("*")[0]
        if(c in ["IGHA1", "IGHA2"]):c = "IGHA1/2"
        if(c in ["IGHG1", "IGHG2"]):c = "IGHG1/2"
        if(c in ["IGHM", "IGHD"]):c = "IGHD/M"
        if(c in ["IGHD/M"]):
          if(mut<=4):c = c+"_unmutated"
          else:c = c+"_mutated"
        else:c = "Class_switched"
        if(c in classes):classes[c] = classes[c]+1#freq[i]
        else:classes[c] = 1#freq[i]
      for c in classes:
        name = c+"\t"+typ
        clust_s = cluster_sizes[cluster_ids[id.split("__")[0]]]
        clust_s = clust_s*100.0/n_uniq_BCRs
        if(name in SHM):
          SHM[name]=SHM[name]+[mut]*classes[c]
          CDR3_len[name] = CDR3_len[name]+[length]*classes[c]
          CDR3_polar[name] =CDR3_polar[name]+[n_polar]*classes[c]
          CDR3_charge[name] =CDR3_charge[name]+[n_charged]*classes[c]
          R_residues[name] = R_residues[name]+[cdr3.count("R")]*classes[c]
          V_residues[name] = V_residues[name]+[cdr3.count("A")]*classes[c]
          cluster_sizes_per_group[name] = cluster_sizes_per_group[name]+[clust_s]*classes[c]
        else:
          SHM[name] = [mut]*classes[c]
          CDR3_len[name]=[length]*classes[c]
          CDR3_polar[name]=[n_polar]*classes[c]
          CDR3_charge[name]=[n_charged]*classes[c]
          R_residues[name]=[cdr3.count("R")]*classes[c]
          V_residues[name]=[cdr3.count("A")]*classes[c]
          cluster_sizes_per_group[name]=[clust_s]*classes[c]
        if(c in overall_CDR3_length):overall_CDR3_length[c] = overall_CDR3_length[c]+[length]*classes[c]
        else:overall_CDR3_length[c]=[length]*classes[c]
  fh.close()
  out="#sample\tIsotype\tReplacement_status\tIGHV status\tmean_SHM\tnumber of sequences\tmean CDR3 length\tmean polar\tmean charged\toverall_mean_CDR3_lengths\tmean(R)residues\tmean(A)residues\tmean_cluster_sizes\n"
  for n in SHM:
    out=out+sample+"\t"+n+"\t"+str(mean(SHM[n]))+"\t"+str(len(SHM[n]))+"\t"+str(mean(CDR3_len[n]))+"\t"+str(mean(CDR3_polar[n]))+"\t"+str(mean(CDR3_charge[n]))+"\t"+str(mean(overall_CDR3_length[n.split("\t")[0]]))+"\t"+str(mean(R_residues[n]))+"\t"+str(mean(V_residues[n]))+"\t"+str(mean(cluster_sizes_per_group[n]))+"\n"
  fh=open(secondary_rearrangement_cdr3_lengths,"w")
  fh.write(out)
  fh.close()
  return()

def sdgfD():
  for id in seqs:
    if(id in cdr3s):
      freq, chains = freqs[id],alias[id].split("|")[1].split("_")
      cdr3 = cdr3s[id]
      length = len(cdr3)
      n_polar,n_charged = 0,0
      for aa in charged:
        n_charged = n_charged+cdr3.count(aa)
      for aa in polar:
        n_polar = n_polar+cdr3.count(aa)
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      for i in nz:
        c = chains[i].split("*")[0]
        if(c in ["IGHA1", "IGHA2"]):c = "IGHA1/2"
        if(c in ["IGHG1", "IGHG2"]):c = "IGHG1/2"
        if(id in ids_rearranged):typ = "REPLACED"
        else:typ = "NOT_REPLACED"
        name = c+"\t"+typ+"\t"+str(length)+"\t"+str(n_charged)+"\t"+str(n_polar)
        if(name in V_replaced_charged_polar_length):V_replaced_charged_polar_length[name] = V_replaced_charged_polar_length[name]+freq[i]
        else:V_replaced_charged_polar_length[name]=freq[i]
  out="#sample\tisotype\tReplacement_status\tCDR3_length\tN_aa_charged\tN_aa_polar\tfrequency\n"
  for name in V_replaced_charged_polar_length:
    out=out+sample+"\t"+name+"\t"+str(V_replaced_charged_polar_length[name])+"\n"
  fh=open(secondary_rearrangement_cdr3_lengths,"w")
  fh.write(out)
  fh.close()
  return()

def Get_reference_sequence_families():
  V_families = {}
  fh=open("LIBRARY/CD_hit_grouping_HOMO_SAPIENS_IGHV.txt","r")
  diff_cutoff_choice = 0.95
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(float(l[0])==diff_cutoff_choice): 
        V_families[l[2].split("*")[0]] = int(l[1])
  fh.close()
  return(V_families)

def Get_secondary_rearrangement_sequences1(sample, annot_file,annot_file3, secondary_rearrangement_file,seq_file,raw_secondary_rearrangement_file):
  V_families = Get_reference_sequence_families()
  fh=open(raw_secondary_rearrangement_file,"w")
  fh.close()
  fh=open(annot_file3,"r")
  seqs,freqs,alias = Get_sequences(seq_file)
  J_regions,stems = Tree(),{}
  for id in freqs:
    break
  total= [0]*len(freqs[id])
  for l in fh:
    l=l.strip().split("\t")
    if(l[0]!='Sequence number'):
      if(len(l)>=18):
        id, ndj_region =  l[1].split("__")[0], l[17]
        n1,n2 = l[21],l[28]
        if(id in freqs):
          total = list(map(add, total, freqs[id]))
          if(min([len(n1),len(n2)]) > 3):
            classification = alias[l[1].split("__")[0]].split("|")[1].split("_")
            if(len(ndj_region)>8):
              ndj_region = ndj_region[3:len(ndj_region)]
              v= l[3].split()[1].split("*")[0]
              v1 = v
              j = l[4].split()[1].split("*")[0]
              if(v in V_families):v = V_families[v]
              J_regions[j][id].value = 1
              v_region = l[8]
              #stem[j][ndj_region][v][id].value = 1
              stems[id] = [v,j,ndj_region, len(ndj_region), l[17],v1,id,v_region]
  fh.close()
  stem,done = Tree(),{}
  for j in J_regions:
    print("J",j)
    stem_regions = []
    for id in J_regions[j]:
      stem_regions.append(stems[id])
    stem_regions.sort(key=lambda x: x[3])
    for i1 in range(0,len(stem_regions)):
      if(stem_regions[i1][6] not in done):
        #done[stem_regions[i1][6]]=1
        v_end = stem_regions[i1][7]
        v_end = v_end[len(v_end)-8:len(v_end)].upper()
        stem_id = stem_regions[i1][4]+":"+j
        stem_seq = stem_regions[i1][4].upper()
        for i2 in range(0,len(stem_regions)):
          #if(i1<i2 and stem_regions[i2][6] not in done):
          if(stem_regions[i2][6] not in done):
            if(stem_regions[i1][0]!=stem_regions[i2][0]):### v region different
            #if(stem_regions[i1][3]<stem_regions[i2][3]): ### different ndn region length
              if(seqs[stem_regions[i2][6]].count(stem_seq)!=0): ### short stem preset in sequence
                #if(seqs[stem_regions[i2][6]].count(v_end)==0): ## different end of V region
                  #if(stem_regions[i2][7].count(stem_regions[i1][7])==0):
                    #print stem_regions[i1][4],stem_regions[i1][7]
                    #print stem_regions[i2][4],stem_regions[i2][7]
                    #print '\t',v_end
                stem[stem_id][stem_regions[i1][0]][stem_regions[i1][6]].value = 1
                stem[stem_id][stem_regions[i2][0]][stem_regions[i2][6]].value = 1
                done[stem_regions[i2][6]] = 1
                done[stem_regions[i1][6]] = 1
                #print "\t",stem_regions[i1][4],stem_regions[i1][5]
                #print "\t",stem_regions[i2][4],stem_regions[i2][5]
    #print j,len(stem_regions)
  #stem_length,lengths = Tree(),[]
  #for stem_reg in stem:
  #  stem_length[len(stem_reg)][stem_reg].value = 1
  #  if(len(stem_reg) not in lengths):lengths.append(len(stem_reg))
  #lengths.sort()
  #stem_combined = Tree()
  done = {}
  #for l in range(0,len(lengths)):
  #  for cdr31 in stem_length[l]:
  #    for l2 in range(l,len(lengths)):
  #      if(l2>l):
  #        for cdr32 in stem_length[l2]:
  ##          if(cdr32 not in done):
  #            if(cdr32.count(cdr31)!=0):
  #              stem_combined[cdr31][cdr32].value = 1
  #              done[cdr32]=1
  #for cdr31 in stem_combined:
    #print "\t",cdr31,len(stem_combined[cdr31])
  #  for cdr32 in stem_combined[cdr31]:
  #    for v in stem[cdr32]:
  #      for id in stem[cdr32][v]:
  #        stem[cdr31][v][id].value = 1
  stem_count = {}
  out,ind="#id_secondary_rearranged\tstem\n",0
  for s in stem:
    if(s not in done):
      if(len(stem[s])>1):
        for v in stem[s]:
          for id in stem[s][v]:
            out = out+id+"\t"+s+"\n"
            ind = ind+1
            if(ind>500):
              Write_out(out, raw_secondary_rearrangement_file)
              out, ind = '',0
            freq = freqs[id]
            nz = [i for i in range(len(freq)) if freq[i]!=0]
            for i in nz:
              c = classification[i].split("*")[0]
              if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
              elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
              if(c in stem_count):stem_count[c] = stem_count[c]+freq[i]
              else:stem_count[c]=freq[i]
  totals = {}
  Write_out(out, raw_secondary_rearrangement_file)
  out, ind = '',0
  for i in range(0,len(total)):
    c = classification[i].split("*")[0]
    if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
    elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
    if(c in totals): totals[c] = totals[c]+total[i]
    else:totals[c] = total[i]
  out="#sample\tchain\tsecondary rearrangement count\ttotal\tpercentage\n"
  for c in totals:
    if(totals[c] !=0):
      if(c not in stem_count):stem_count[c] = 0
      out=out+sample+"\t"+c+"\t"+str(stem_count[c])+"\t"+str(totals[c])+"\t"+str(stem_count[c]*100.0/totals[c])+"\n"
  fh=open(secondary_rearrangement_file,"w")
  fh.write(out)
  fh.close()
  print(out)
  return()

def Get_secondary_rearrangement_sequences(sample, annot_file,annot_file3, secondary_rearrangement_file,seq_file,raw_secondary_rearrangement_file,reverse_primer_group):
  V_families = Get_reference_sequence_families()
  fh=open(raw_secondary_rearrangement_file,"w")
  fh.close()
  fh=open(secondary_rearrangement_file,"w")
  fh.close()
  fh=open(annot_file3,"r")
  seqs,freqs,alias = Get_sequences(seq_file)
  J_regions,stems = Tree(),{}
  for id in freqs:
    break
  total= [0]*len(freqs[id])
  total_all, total_IGHDM = 0,0
  for l in fh:
    l=l.strip().split("\t")
    if(l[0]!='Sequence number'):
      if(len(l)>=18):
        id, ndj_region =  l[1].split("__")[0], l[39]
        junction = l[14]
        n1,n2 = l[21],l[28]
        if(id in freqs):
          if(min([len(n2)]) > 4 ):
            nz = [i for i in range(len(freqs[id])) if freqs[id][i]!=0]
            clas = alias[id].split("|")[1].split("_")
            total_all = total_all+1
            IGHDM = 0
            for i in nz :
              total[i] = total[i]+1
              if(clas[i].count("IGHD")==1):IGHDM =1
              if(clas[i].count("IGHM")==1):IGHDM =1
            if(IGHDM==1):total_IGHDM = total_IGHDM+1
            #total = map(add, total, freqs[id])
            classification = alias[l[1].split("__")[0]].split("|")[1].split("_")
            if(len(ndj_region)>8):
              ndj_region = ndj_region[3:len(ndj_region)]
              v= l[3].split()[1].split("*")[0]
              v1 = v
              j = l[4].split()[1].split("*")[0]
              if(v in V_families):v = V_families[v]
              J_regions[j][id].value = 1
              v_end = l[8]
              #stem[j][ndj_region][v][id].value = 1
              v_region = v_end[len(v_end)-2:len(v_end)].upper()
              stems[id] = [v,j,ndj_region, len(ndj_region), ndj_region.upper() ,v1,id,v_region,junction.upper()+v_region]
  fh.close()
  stem,done = Tree(),{}
  print(len(stems))
  for j in J_regions:
    print("J",j)
    stem_regions = []
    for id in J_regions[j]:
      stem_regions.append(stems[id])
    stem_regions.sort(key=lambda x: x[3])
    for i1 in range(0,len(stem_regions)):
      if(stem_regions[i1][6] not in done):
        #done[stem_regions[i1][6]]=1
        v_end = stem_regions[i1][7]
        ### v_end = v_end[len(v_end)-8:len(v_end)].upper()
        stem_id = stem_regions[i1][4]+":"+j
        stem_seq = stem_regions[i1][4]
        for i2 in range(0,len(stem_regions)):
          #if(i1<i2 and stem_regions[i2][6] not in done):
          if(stem_regions[i2][6] not in done):
            if(stem_regions[i1][0]!=stem_regions[i2][0]):### v region different
              if(stem_regions[i1][8]!=stem_regions[i2][8]): #### junction is different
                if(seqs[stem_regions[i2][6]].count(stem_seq)!=0): ### short stem preset in sequence
                  stem[stem_id][stem_regions[i1][0]][stem_regions[i1][6]].value = 1
                  stem[stem_id][stem_regions[i2][0]][stem_regions[i2][6]].value = 1
                  done[stem_regions[i2][6]] = 1
                  done[stem_regions[i1][6]] = 1
  done = {}
  stem_count = {}
  out,ind="#id_secondary_rearranged\tstem\n",0
  for s in stem:
    if(s not in done):
      if(len(stem[s])>1):
        for v in stem[s]:
          for id in stem[s][v]:
            out = out+id+"\t"+s+"\n"
            ind = ind+1
            if(ind>500):
              Write_out(out, raw_secondary_rearrangement_file)
              out, ind = '',0
            freq = freqs[id]
            nz = [i for i in range(len(freq)) if freq[i]!=0]
            c = "ALL"
            if(c in stem_count):stem_count[c] = stem_count[c]+1#freq[i]
            else:stem_count[c]=1#freq[i]
            IGHDM = 0
            for i in nz:
              c = classification[i].split("*")[0]
              if(reverse_primer_group!="ISOTYPER"):
                if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
                elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
              if(c in stem_count):stem_count[c] = stem_count[c]+1#freq[i]
              else:stem_count[c]=1#freq[i]
              if(c in ["IGHD","IGHM"]):IGHDM = 1
            if(IGHDM==1):
              c = "IGHD_M"
              if(c in stem_count):stem_count[c] = stem_count[c]+1#freq[i]
              else:stem_count[c]=1#freq[i]
  totals = {}
  Write_out(out, raw_secondary_rearrangement_file)
  out, ind = '',0
  for i in range(0,len(total)):
    c = classification[i].split("*")[0]
    if(reverse_primer_group!="ISOTYPER"):
      if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
      elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
    if(c in totals): totals[c] = totals[c]+total[i]
    else:totals[c] = total[i]
  totals["ALL"] =total_all
  totals["IGHD_M"] =total_IGHDM
  out="#sample\tchain\tsecondary rearrangement count\ttotal\tpercentage\n"
  for c in totals:
    if(totals[c] !=0):
      if(c not in stem_count):stem_count[c] = 0
      out=out+sample+"\t"+c+"\t"+str(stem_count[c])+"\t"+str(totals[c])+"\t"+str(stem_count[c]*100.0/totals[c])+"\n"
  fh=open(secondary_rearrangement_file,"w")
  fh.write(out)
  fh.close()
  print(out)
  return()

def Subsample_per_Vgene(raw_secondary_rearrangement_file, sample_depth, secondary_rearrangement_file_SAMPLED,cluster_file,seq_file,sample,annot_file):
  CDR3s,annots = Get_annotations(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  for i in alias:
    classification = alias[i].split("|")[1].split("_")
    break
  fh=open(raw_secondary_rearrangement_file,"r")
  inv_stem ={}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      stem, id = l[1],l[0]
      #stems[stem][id].value = 1
      inv_stem[id] =stem
  fh.close()
  Vgene_array,genes_inv = {},{}
  subsample_groups = Tree()
  for id in annots:
    V = annots[id][3].split("\t")[0]
    freq = freqs[id]
    nz = [i for i in range(len(freq)) if freq[i]!=0]
    mut =  annots[id][0]
    classes = {}
    for i in nz:
      c = classification[i].split("*")[0]
      if(c in ["IGHD","IGHM"]):
        if(mut<=3):c = "IGHD/M_unmutated"
        else:c = "IGHD/M_mutated"
      else:c = "Class_switched"
      classes[c] = 1
    if("Class_switched" in classes):c = "Class_switched"
    else:
      for c in classes:
        break
    c1 = c+"\t"+V
    subsample_groups[c][c1].value = 1
    if(c1 in Vgene_array):Vgene_array[c1] = Vgene_array[c1]+[id]
    else:Vgene_array[c1]=[id]
    genes_inv[id] = c1
  sample_depth = 20 ## depth of sampling per IGHV gene
  repeats = 10000
  out = "#sample\tisotype\tV_gene\tsample depth per V gene\tn repeats\tn v genes assessed\tuniq_id_replacement_freq per v gene\n"
  for c in subsample_groups:
    v_genes_sample = []
    for v in subsample_groups[c]:
      if(len(Vgene_array[v])>=sample_depth):
        v_genes_sample.append(v)
    print(len(v_genes_sample), len(subsample_groups[c]))
    v_gene_editing_frequencies_overall = {}
    for v in v_genes_sample:
      v_gene_editing_frequencies_overall[v] = []
    for r in range(0,repeats):
      ids_sampled = []
      for v in v_genes_sample:
        rand_ids = np.random.choice(Vgene_array[v], sample_depth,replace=False)
        ids_sampled = ids_sampled+rand_ids.tolist()
      sub_stems,V_editing_frequencies = Tree(),{}
      for id in ids_sampled:
        if(id in inv_stem):
          sub_stems[inv_stem[id]][genes_inv[id]][id].value = 1
      for v in v_genes_sample:
        V_editing_frequencies[v] = 0
      for stem in sub_stems:
        if(len(sub_stems[stem])>1):
          for v in sub_stems[stem]:
            V_editing_frequencies[v]= V_editing_frequencies[v]+1
      for v in V_editing_frequencies:
        v_gene_editing_frequencies_overall[v] = v_gene_editing_frequencies_overall[v]+[V_editing_frequencies[v]]
    for v in v_gene_editing_frequencies_overall:
      out=out+ "\t".join(map(str, [sample, v,sample_depth,repeats,len(v_genes_sample), mean(v_gene_editing_frequencies_overall[v])]))+"\n"
  print(out)
  fh=open(secondary_rearrangement_file_SAMPLED.replace("Secondary_rearrangement_file_","Secondary_NORMALISED_PER_GENE_rearrangement_file_"),"w")
  fh.write(out)
  fh.close()
  return()

def Subsample_secondary_rearrangments(raw_secondary_rearrangement_file, sample_depth, secondary_rearrangement_file_SAMPLED,cluster_file,seq_file,sample,annot_file):
  CDR3s,annots = Get_annotations(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  for i in alias:
    classification = alias[i].split("|")[1].split("_")
    break
  fh=open(raw_secondary_rearrangement_file,"r")
  stems,inv_stem = Tree(),{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      stem, id = l[1],l[0]
      stems[stem][id].value = 1
      inv_stem[id] =stem
  fh.close()
  ids_array = []
  id_annotation = {}
  total = 0
  classes_all = Tree()
  for id in seqs:
    if(id in annots):
      classes = {}
      freq = freqs[id]
      ids_array = ids_array+[id]*sum(freq)
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      mut =  annots[id][0]
      for i in nz:
        c = classification[i].split("*")[0]
        if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
        elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        elif(c in ["IGHD","IGHM"]):c = "IGHD/M"
        if(c in classes):classes[c] = classes[c]+freq[i]
        else:classes[c]=freq[i]
        if(c in ["IGHD/M"]):
          if(mut<=3):c = "IGHD/M_unmutated"
          else:c = "IGHD/M_mutated"
          if(c in classes):classes[c] = classes[c]+freq[i]
          else:classes[c]=freq[i]
          if(c == "IGHD/M_unmutated" and sum(freq)<=5):
            c = "IGHD/M_unmutated_unexpanded"
            if(c in classes):classes[c] = classes[c]+freq[i]
            else:classes[c]=freq[i]
      c = "All"
      if(c in classes):classes[c] = classes[c]+sum(freq)
      else:classes[c]=sum(freq)
      array = []
      for c in classes:
        array.append([c,classes[c]])
        classes_all[c][id].value = 1
      array.append(["All",sum(freq)])
      id_annotation[id] = array
      total = total+sum(freq)
  repeats = 1000
  fh=open(secondary_rearrangement_file_SAMPLED,"w")
  fh.close()
  out = '#sample\tsample_ID\tisotype\tsample depth\tisotype_total\tid_replacement_freq\tuniq_id_replacement_freq\treplacement_v_gene_incidence\treplacement_incidence\n'
  print(len(ids_array),total)
  for clas in classes_all:
    if(len(classes_all[clas])>sample_depth):
      ids_array1 = []
      for id in classes_all[clas]:
        ids_array1.append(id)
      for r in range(0,repeats):
        rep_id = "REP"+str(r)+"_"+sample
        sample1 = random.sample(ids_array1, sample_depth)
        stems_sample = Tree()
        classes_total={}
        for id in sample1:
          if(id in annots):
            for c1 in id_annotation[id]:
              c = c1[0]
              if(c in classes_total):classes_total[c] = classes_total[c]+c1[1]
              else:classes_total[c]=c1[1]
            if(id in inv_stem):
              v_gene = annots[id][3].split("\t")[0]
              stems_sample[inv_stem[id]][v_gene][id].value = 1
        stem_count, stem_count_uniq, stem_count_incidence = {},{},{}
        stem_count_incidence["All"],stem_count_incidence["v"] =0,0
        for stem in stems_sample:
          if(len(stems_sample[stem])>1):
            for v in stems_sample[stem]:
              for id in stems_sample[stem][v]:
                c = clas
                if(c in stem_count):stem_count[c] = stem_count[c]+freq[i]
                else:stem_count[c]=freq[i]
                if(c in stem_count_uniq):stem_count_uniq[c] = stem_count_uniq[c]+1
                else:stem_count_uniq[c] = 1
              stem_count_incidence["v"] = stem_count_incidence["v"]+1
              stem_count_incidence["All"] = stem_count_incidence["All"]+1
        for c in stem_count:
          array, X = [],[stem_count, stem_count_uniq]
          for i in range(0,len(X)):
            x = X[i]
            if(c in x):array.append(x[c])
            else:array.append(0)
          out=out+"\t".join(map(str, [sample, rep_id, c, sample_depth, len(classes_all[clas])]+array+[stem_count_incidence["v"],stem_count_incidence["All"]]))+"\n"
    fh=open(secondary_rearrangement_file_SAMPLED,"w")
    fh.write(out)
    fh.close()
  return()

def Screen_auto_Ab_sequences(annot_file4, sample, seq_file, captured_auto_antibodies_file):
  file_reference_autoAb = "/nfs/users/nfs_r/rbr1/IMMUNOLOGY/NETWORK_METHODS/STANDARD_BCR_PROCESSING/LIBRARY/Auto_reactive_IgGpos_BCR_reference.txt"
  fh=open(file_reference_autoAb,"r")
  sequences_Ab, seq_type = {},Tree()
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,v,j,CDR3 = l[0],"IGHV"+l[1],"IGHJ"+l[4], l[6].replace("_","")
      sequences_Ab[id] = [v,j,CDR3]
      seq_type[v+":"+j][CDR3][id].value = 1
  fh.close()
  seqs,freqs,alias = Get_sequences(seq_file)
  fh = open(annot_file4 ,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split("\t")
      if(l[1].split("__")[0] in freqs):
        if(l[2]=="productive" and len(l)>=18):
          v= l[3].split()[1].split("*")[0]
          j= l[4].split()[1].split("*")[0]
          if(v+":"+j in seq_type):
            cdr3 = l[14]
            print(cdr3)
            if(cdr3 in seq_type[v+":"+j]):
              print(l)
  fh.close()



def Get_percentage_public_repertoire(sample,seq_unmutated_file,annot_file,public_repertoire):
  fh = open("/lustre/scratch108/viruses/rbr1/MISEQ/AUTOIMMUNITY_AAV_SLE_HEALTHY1/ORIENTATED_SEQUENCES/ANNOTATIONS/Public_CDR3_PUBLIC_BCR_healthy.txt","r")
  public_CDR3 = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(int(l[len(l)-1])>3):
        public_CDR3[l[0]] = int(l[len(l)-1])
  fh.close()
  cdr3s = Get_CDR3s(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  total,public = {},{}
  fh = open(per_cluster_developmental_classification_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[1] in cdr3s and l[1] in seqs):
        CDR3 = cdr3s[l[1]]
        classification = alias[l[1].split("__")[0]].split("|")[1].split("_")
        freq = freqs[l[1].split("__")[0]]
        nz = [i for i in range(len(freq)) if freq[i]!=0]
        for i in nz:
          c = classification[i].split("*")[0]
          if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
          elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
          if(c in total):total[c] = total[c]+freq[i]
          else:total[c]=freq[i]
          if(CDR3 in public_CDR3):
            if(c in public):public[c] = public[c]+freq[i]
            else:public[c]=freq[i]
  fh.close()
  out = "#sample\tclassification\ttotal\tpublic\tpercent public\n"
  for c in total:
    if(c in public):pub = public[c]
    else:pub = 0
    out=out+sample+"\t"+c+"\t"+str(total[c])+"\t"+str(pub)+"\t"+str(pub*100.0/total[c])+"\n"
  fh=open(public_repertoire,"w")
  fh.write(out)
  fh.close()
  return()

def Get_CDR3_region_composition(sample,seq_file,annot_file, per_cluster_developmental_classification_file,CDR3_region_composition,CDR3_5prime_charge,CDR3_length_file,annot_file4,CDR3_charge_file,CDR23_charge_file,mean_IgD_IgM_ratio,CDR3_5prime_charge_long_short,CDR3_length_vjgenes):
  #Amino_acid_compistion(sample,seq_file,annot_file, per_cluster_developmental_classification_file,CDR3_region_composition)
  #Get_CDR3_5prime_charge(CDR3_5prime_charge, sample,seq_file,annot_file,CDR3_5prime_charge_long_short)
  #Get_CDR3_lengths(seq_file,annot_file,CDR3_length_file,sample,annot_file4,CDR3_charge_file,CDR23_charge_file)
  Get_CDR3_length_VJ_genes(seq_file,annot_file,CDR3_length_vjgenes,sample)
  #Get_mean_IgD_IgM_ratio(mean_IgD_IgM_ratio, seq_file,sample)
  return()

def Get_mean_IgD_IgM_ratio(mean_IgD_IgM_ratio, seq_file,sample):
  seqs,freqs,alias = Get_sequences(seq_file)
  ratio = {}
  for id in seqs:
    f, chains = freqs[id],alias[id].split("|")[1].split("_")
    nz = [i for i in range(len(f)) if f[i]!=0]
    IGHD,IGHM = 0,0
    for i in nz:
      c = chains[i].split("*")[0]
      if(c =="IGHD"):IGHD = f[i]
      if(c =="IGHM"):IGHM = f[i]
    if(IGHD+IGHM >0):
      tot = IGHD+IGHM
      if(tot in ratio):ratio[tot] = [ratio[tot][0]+IGHD, ratio[tot][1]+tot]
      else:ratio[tot] = [IGHD, tot]
  cumul_ratios,totals = {},[]
  for i in ratio:
    totals.append(i)
  totals.sort()
  for tot in ratio:
    cumul_ratios[tot] = [0,0]
    for t1 in ratio:
      if(t1>=tot):
         cumul_ratios[tot] = [cumul_ratios[tot][0]+ratio[t1][0], cumul_ratios[tot][1]+ratio[t1][1]]
  out="#sample\tmin_freq\tmean_IGHD_proportion\n"
  for tot in cumul_ratios:
    out=out+sample+"\t"+str(tot)+"\t"+str(cumul_ratios[tot][0]*100.0/cumul_ratios[tot][1])+"\n"
  fh=open(mean_IgD_IgM_ratio,"w")
  fh.write(out)
  fh.close()
  return()

def Get_CDR3_length_VJ_genes(seq_file,annot_file,CDR3_length_vjgenes,sample):
  CDR3_length_vjgenes1 = CDR3_length_vjgenes.replace("CDR3_VJ_genes_","CDR3_length_array_")
  seqs,freqs,alias = Get_sequences(seq_file)
  for id in alias:
    break
  chains = []
  classifications = alias[id].split("|")[1].split("_")
  for c in classifications:
    c = c.split("*")[0]
    if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
    if(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
    chains.append(c)
  fh=open(annot_file,"r")
  cdr3s = {}
  cdr3_array_template = [0]*60
  cdr3_lengths_array = {}
  mutations = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split("\t")
      if(l[2].count('productive')):
        if(len(l)>=20):
          cdr3 = l[20].upper()
          v,j = l[3].split()[1],l[9].split()[1]
          if(len(cdr3)>3 and l[1].split("__")[0] in freqs):
            v,j, = v.split("*")[0], j.split("*")[0]
            name = v+"\t"+j+"|"+str(len(cdr3))
            freq = freqs[l[1].split("__")[0]]
            nz = [i for i in range(len(freq)) if freq[i]!=0]
            chains_sub = {}
            v_muts, j_muts,mut_n = l[6], l[12] ,0
            if(v_muts.count('nt')!=0 and j_muts.count('nt')!=0):
              v_muts, j_muts = list(map(int,v_muts.split()[0].split("/"))), list(map(int,j_muts.split()[0].split("/")))
              v_muts, j_muts = v_muts[1]-v_muts[0], j_muts[1]-j_muts[0]
              mut_n = 1
            for i in nz:
              chains_sub[chains[i]] = 1
            #j1 = "NA"
            for c in chains_sub:
              name1 = c+"\t"+name
              if(name1 in cdr3s):cdr3s[name1] = cdr3s[name1]+1
              else:cdr3s[name1] = 1
              name1,value = c,len(cdr3)
              if(name1 in cdr3_lengths_array):
                array = copy.deepcopy(cdr3_lengths_array[name1])
                if(value<60):
                  array[value] = array[value]+1
                cdr3_lengths_array[name1] = array
              else:
                array = copy.deepcopy(cdr3_array_template)
                if(value<40):array[value] = array[value]+1
                cdr3_lengths_array[name1] = array
              name1 = c+"\t"+v+"\t"+j
              if(name1 in mutations):mutations[name1]= mutations[name1]+[v_muts]
              else:mutations[name1]=[v_muts]
            if("IGHD" in chains_sub or 'IGHM' in chains_sub):
              if(mut_n==1):
                c = "IGHD/M_Unmutated"
                if(v_muts+j_muts > 1):c = "IGHD/M_Mutated"
                name1 = c+"\t"+name
                if(name1 in cdr3s):cdr3s[name1] = cdr3s[name1]+1
                else:cdr3s[name1] = 1
                name1 = c+"\t"+v+"\t"+j
                if(name1 in mutations):mutations[name1]= mutations[name1]+[v_muts]
                else:mutations[name1]=[v_muts]
            else:
              c = "Class_switched"
              name1 = c+"\t"+v+"\t"+j
              if(name1 in mutations):mutations[name1]= mutations[name1]+[v_muts]
              else:mutations[name1]=[v_muts]
              name1 = c+"\t"+name
              if(name1 in cdr3s):cdr3s[name1] = cdr3s[name1]+1
              else:cdr3s[name1] = 1
            #if(len(mutations[name1])>10):break
  fh.close()
  out = ''
  mutations_vjgenes= CDR3_length_vjgenes.replace("CDR3_","Mutations_")
  out, ind=  '#sample\tisotype\tv\tj\tmean_mutations\tn_unique_BCRs\n',0
  for n in mutations:
    out=out+sample+"\t"+n+"\t"+str(mean(mutations[n]))+"\t"+str(len(mutations[n]))+"\n"
  fh=open(mutations_vjgenes,"w")
  fh.write(out)
  fh.close()
  for a in cdr3_lengths_array:
    out=out+sample+"\t"+"\t".join(map(str, cdr3_lengths_array[a]))+"\n"
  fh=open(CDR3_length_vjgenes1, "w")
  fh.write(out)
  fh.close()
  mean_CDR3_lengths = {}
  for name in cdr3s:
    length = int(name.split("|")[1])
    vj=name.split("|")[0]
    if(vj in mean_CDR3_lengths):mean_CDR3_lengths[vj] = mean_CDR3_lengths[vj]+[length]*cdr3s[name]
    else:mean_CDR3_lengths[vj]=[length]*cdr3s[name]
  del cdr3s
  fh=open(CDR3_length_vjgenes,"w")
  fh.close()
  out, ind=  '#sample\tisotype\tv\tj\tmean_CDR3_length\tmax_length\tmin_length\t95th_percentile\t5percentile\tnumber_f_unique_sequences\n',0
  for name in mean_CDR3_lengths:
    perc1,perc2 = np.percentile(np.array(mean_CDR3_lengths[name]), 5), np.percentile(np.array(mean_CDR3_lengths[name]), 95)
    out=out+sample+"\t"+name+"\t"+str(mean(mean_CDR3_lengths[name]))+"\t"+str(max(mean_CDR3_lengths[name]))+"\t"+str(min(mean_CDR3_lengths[name]))+"\t"+str(perc2)+"\t"+str(perc1)+"\t"+str(len(mean_CDR3_lengths[name]))+"\n"
    ind = ind+1
    if(ind>300):
      Write_out(out, CDR3_length_vjgenes)
      out, ind = '',0
  Write_out(out, CDR3_length_vjgenes)
  return()

def Get_autoreactive_V4_34_motif_frequency(annot_file4, V4_34_quantification_file, sample,seq_file,reverse_primer_group,annot_file):
  seqs,freqs,alias = Get_sequences(seq_file)
  CDR3s,mutations = Get_annotations(annot_file)
  for id in alias:
    classes = alias[id].split("|")[1].split("_")
    break
  fh=open(annot_file4, "r")
  print(annot_file4)
  avy, nhs, avy_nhs,total,v4_34 = {},{},{},{},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split("\t")
      id = l[1].split("__")[0]
      if(len(l)>=18 and id in mutations):
        v = l[3].split()[1].split("*")[0]
        FWR1,CDR2 = l[9],l[12]
        freq = freqs[id]
        nz = [i for i in range(len(freq)) if freq[i]!=0]
        ##### between isotypes
        classes1 = {}
        for i in nz:
          c = classes[i].split("*")[0]
          if(reverse_primer_group!="ISOTYPER"):
            if(c in ["IGHA1","IGHA2"]):c="IGHA1/2"
            if(c in ["IGHG1","IGHG2"]):c="IGHG1/2"
          if(c  in ["IGHD","IGHM"]):
            mut = mutations[id]
            if(mut[0]+mut[1]<=3):c = "IGHD/M_unmutated"
            else:c = "IGHD/M_mutated"
          else:c = "class_switched"
          classes1[c] = 1
        if("class_switched" in classes1):c = "class_switched"
        else:
          for c in classes1:
            break
        #print l[1], c, classes
        classes1 = [c,"all"]
        for c in classes1:
          f1 = sum(freq)
          if(c in total):
            a = [total[c][0]+f1, total[c][1]+1]
            total[c] = a
          else:total[c]=[f1,1]
          for x in [avy, nhs, avy_nhs]:
            if(c not in x):
              x[c] = [0,0,0,0]
          for x in [v4_34]:
            if(c not in x):
              x[c] = [0,0]
          if(v=="IGHV4-34"):
            if(c in v4_34):f = v4_34[c]
            else:f = [0,0]
            f = [f[0]+freq[i], f[1]+1]
            v4_34[c] = f
            if(len(FWR1)>=3):
              if(c in avy):f = avy[c]
              else:f = [0,0,0,0]
              if(FWR1.count("AVY")==1):f = [f[0]+freq[i], f[1]+1,f[2]+freq[i], f[3]+1 ] ########### total with annotation, total with unmutated annotation
              else:f = [f[0]+freq[i], f[1]+1,f[2],f[3]]
              avy[c] = f
            if(len(CDR2)>=3):
              if(c in nhs):f = nhs[c]
              else:f = [0,0,0,0]
              if(CDR2.count("NHS")==1):f = [f[0]+freq[i], f[1]+1,f[2]+freq[i], f[3]+1 ] 
              else:f = [f[0]+freq[i], f[1]+1,f[2],f[3]]
              nhs[c] = f
              if(len(FWR1)>=3):
                if(c in avy_nhs):f = avy_nhs[c]
                else:f = [0,0,0,0]
                if(FWR1.count("AVY")==1 and CDR2.count("NHS")==1):f = [f[0]+freq[i], f[1]+1,f[2]+freq[i], f[3]+1 ]
                else:f = [f[0]+freq[i], f[1]+1,f[2],f[3]]
                avy_nhs[c] = f
  fh.close()
  out = "#sample\tisotype\ttotal_all\ttotal_uniq\tV4_34_total\tV4_34_unique\tV4_34_AVY_total\tV4_34_AVY_uniq\tV4_34_AVY_total_unmut\tV4_34_AVY_uniq_unmut\tV4_34_NHS_total\tV4_34_NHS_uniq\tV4_34_NHS_total_unmut\tV4_34_NHS_uniq_unmut\tV4_34_AVY_NHS_total\tV4_34_NHS_AVY_uniq\tV4_34_AVY_NHS_total_unmut\tV4_34_NHS_AVY_uniq_unmut\n"
  for c in total:
    array = total[c]+v4_34[c]+avy[c]+nhs[c]+avy_nhs[c]
    out=out+sample+"\t"+c+"\t"+"\t".join(map(str, array))+"\n"
  fh=open(V4_34_quantification_file, "w")
  fh.write(out)
  fh.close()
  return()

def Get_total_network_summary(sample, annot_file, isotype_usages_SUBSAMPLED_clone,reverse_primer_group):
  seqs,freqs,alias = Get_sequences(seq_file)
  CDR3s,annots = Get_annotations(annot_file)
  mutations, counts = {},{}
  for id in seqs:
    if(id in annots):
      f, chains = freqs[id],alias[id].split("|")[1].split("_")
      nz = [i for i in range(len(f)) if f[i]!=0]
      chains_sub = {}
      class_switched = 0
      for i in nz:
        c = chains[i].split("*")[0]
        if(reverse_primer_group!="ISOTYPER"):
          if(c in ["IGHA1", "IGHA2"]):c,class_switched = "IGHA1/2",1
          if(c in ["IGHG1", "IGHG2"]):c,class_switched = "IGHG1/2",1
        if(c not in ["IGHD","IGHM"]):class_switched = 1
        if(c in chains_sub):chains_sub[c] = chains_sub[c]+f[i]
        else: chains_sub[c]=f[i]
      mut =  annots[id][0]#+annots[id][1]
      if(class_switched==0):
        c1 = ["IGHD/M"]
        if(mut<=5):
          c1.append("IGHD/M_unmutated")
          if(sum(f)<5):c1.append("IGHD/M_unmutated_unexpanded")
        else:c1.append("IGHD/M_mutated")
      else:
        c1 = ["Class_switched"]
      for c in c1: 
        chains_sub[c]=1
      for c in chains_sub:
        if(c in mutations):
          mutations[c] = mutations[c]+[mut]
          counts[c] = counts[c]+1
        else:
          mutations[c]=[mut]
          counts[c]=1
  out = "#Id\tsubsample\tIsotype\tuniq_BCR_frequency\tmean_mutations\n"
  for c in counts:
    out=out+"\t".join(map(str, [sample, c, counts[c], mean(mutations[c])]))+"\n"
  print(out)
  fh=open(isotype_usages_SUBSAMPLED_clone,"w")
  fh.write(out)
  fh.close()
  return()

def Get_CDR3_lengths(seq_file,annot_file,CDR3_length_file,sample,annot_file4,CDR3_charge_file,CDR23_charge_file,reverse_primer_group,CDR3_distribution):
  cdr3s = Get_CDR3s(annot_file)
  cdrs = Get_aa_CDR(annot_file4)
  seqs,freqs,alias = Get_sequences(seq_file)
  lengths,R_K_residues_CDR3,R_K_residues_CDR2_3 = {},{},{}
  for id in seqs:
    if(id in cdr3s): 
      f, chains = freqs[id],alias[id].split("|")[1].split("_")
      cdr3 = cdr3s[id]
      length = len(cdr3)
      RK_CDR3 = cdrs[id][1].count("R")+cdrs[id][1].count("K")
      RK_CDR2_3 = RK_CDR3 +cdrs[id][0].count("R")+cdrs[id][0].count("K")
      nz = [i for i in range(len(f)) if f[i]!=0]
      for i in nz:
        c = chains[i].split("*")[0]
        if(reverse_primer_group!="ISOTYPER"):
          if(c in ["IGHA1", "IGHA2"]):c = "IGHA1/2"
          if(c in ["IGHG1", "IGHG2"]):c = "IGHG1/2"
        name = c
        if(name in lengths):
          lengths[name] = lengths[name]+[length]*1
          R_K_residues_CDR3[name] = R_K_residues_CDR3[name]+[RK_CDR3]*f[i]
          R_K_residues_CDR2_3[name] = R_K_residues_CDR2_3[name]+[RK_CDR2_3]*f[i]
        else:
          lengths[name]=[length]*1
          R_K_residues_CDR3[name]=[RK_CDR3]*f[i]
          R_K_residues_CDR2_3[name]=[RK_CDR2_3]*f[i]
  out="#sample\tisotype\tmean_CDR3_length\tNumber_of_BCRs\n"
  out1 = "#sample\tisotype\tCDR3_length\tNumber_of_BCRs\n"
  for c in lengths:
    out=out+sample+"\t"+c+"\t"+str(mean(lengths[c]))+"\t"+str(len(lengths[c]))+"\n"
    length = lengths[c]
    counts= [ (i,length.count(i)) for i in set(length) ]
    for i in counts:
      out1=out1+sample+"\t"+c+"\t"+str(i[0])+"\t"+str(i[1])+"\n"
  fh=open(CDR3_length_file,"w")
  fh.write(out)
  fh.close()
  fh=open(CDR3_distribution,"w")
  fh.write(out1)
  fh.close()
  out="#sample\tisotype\tmean_CDR3_R_K_residues\tNumber_of_BCRs\n"
  for c in R_K_residues_CDR3:
    out=out+sample+"\t"+c+"\t"+str(mean(R_K_residues_CDR3[c]))+"\t"+str(len(R_K_residues_CDR3[c]))+"\n"
  fh=open(CDR3_charge_file,"w")
  fh.write(out)
  fh.close()
  out="#sample\tisotype\tmean_CDR2/3_R_K_residues\tNumber_of_BCRs\n"
  for c in R_K_residues_CDR2_3:
    out=out+sample+"\t"+c+"\t"+str(mean(R_K_residues_CDR2_3[c]))+"\t"+str(len(R_K_residues_CDR2_3[c]))+"\n"
  fh=open(CDR23_charge_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_charged_polar_positions():
  file = "/lustre/scratch118/infgen/team146/rbr1/REFERENCE_FILES/Amino_acid_properties4.txt"
  fh=open(file, "r")
  charged, polar = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      charge, pol = l[5],l[6]
      if(charge !="No"):charged[l[0]] = charge
      if(pol !="No"):polar[l[0]] = pol
  fh.close()
  return(charged, polar)

def Get_CDR3_5prime_charge(CDR3_5prime_charge, sample,seq_file,annot_file,CDR3_5prime_charge_long_short):
  charged, polar = Get_charged_polar_positions()
  cdr3s = Get_CDR3s(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  charged_polar = {}
  lengths,charge_per_length = {},{}
  for id in seqs:
    if(id in cdr3s):
      freq, chains = freqs[id],alias[id].split("|")[1].split("_")
      cdr3 = cdr3s[id]
      length = len(cdr3)
      if(len(cdr3)>8):cdr3 = cdr3[0:8]
      n_polar,n_charged = 0,0
      for aa in charged:
        n_charged = n_charged+cdr3.count(aa)
      for aa in polar: 
        n_polar = n_polar+cdr3.count(aa)
      lengths[length]=1
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      for i in nz:
        c = chains[i].split("*")[0]
        if(c in ["IGHA1", "IGHA2"]):c = "IGHA1/2"
        if(c in ["IGHG1", "IGHG2"]):c = "IGHG1/2"
        name = c+"\t"+str(n_charged)+"\t"+str(n_polar)
        if(name in charged_polar):charged_polar[name] = charged_polar[name]+freq[i]
        else:charged_polar[name]=freq[i]
        name = c+"\t"+str(length)+"\t"+str(n_charged)+"\t"+str(n_polar)
        if(name in charge_per_length):charge_per_length[name] = charge_per_length[name] + freq[i]
        else:charge_per_length[name] = freq[i]
  length = []
  for l in lengths:
    length.append(l)
  length = np.array(length)
  lower,upper = np.percentile(length, 25), np.percentile(length, 75)
  lower,upper = 15, 22
  polar_length = {}
  charge_length = {}
  for name in charge_per_length:
    info = name.split("\t")
    clas, length, n_charged,n_polar = info[0],int(info[1]),float(info[2]),float(info[3])
    number = charge_per_length[name]
    group = "MID"
    if(length<lower):group = "LOWER"
    elif(length>=upper):group = "UPPER"
    if(group in polar_length):
      polar_length[group] = polar_length[group]+[n_polar]*number
      charge_length[group] = charge_length[group]+[n_charged]*number
    else:
      polar_length[group]=[n_polar]*number
      charge_length[group]=[n_charged]*number
  out="#sample\tCDR3_length_group\tlength_param\tmean_polar_residues\tmean_charged_residues\tnumber_compared\n"
  params,groups = ["<"+str(lower)+"aa", str(lower)+"-"+str(upper)+"aa", ">"+str(upper)+"aa"],["LOWER", "MID", "UPPER"]
  for group in groups:
    param = params[groups.index(group)]
    out=out+sample+"\t"+group+"\t"+param+"\t"+str(mean(polar_length[group]))+"\t"+str(mean(charge_length[group]))+"\t"+str(len(polar_length[group]))+"\n"
  fh=open(CDR3_5prime_charge_long_short,"w")
  fh.write(out)
  fh.close()
  out="#sample\tisotype\tN_aa_charged\tN_aa_polar\tfrequency\n"
  for name in charged_polar:
    out=out+sample+"\t"+name+"\t"+str(charged_polar[name])+"\n"
  fh=open(CDR3_5prime_charge,"w")
  fh.write(out)
  fh.close()
  return()

def Amino_acid_compistion(sample,seq_file,annot_file, per_cluster_developmental_classification_file,CDR3_region_composition):
  aas = []
  fh = open("/lustre/scratch108/viruses/rbr1/REFERENCE_FILES/Amino_acid_properties2.txt","r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      aas.append(l[0])
  fh.close()
  cdr3s = Get_CDR3s(annot_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  fh = open(per_cluster_developmental_classification_file,"r")
  composition,total = {},0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[1] in cdr3s and l[1] in seqs):
        classes,total = l[2].split("|"),total+1
        CDR3 = cdr3s[l[1]]
        freq_aa = [0]*len(aas)
        for i in range(0,len(aas)):
          freq_aa[i] = CDR3.count(aas[i])
        freq_aa = "\t".join(map(str,freq_aa))
        f = sum(freqs[l[1]])
        for c in classes:
          name = c+":"+freq_aa
        if(name in composition):composition[name] = composition[name]+f
        else:composition[name]=f
  fh.close()
  out = "#ID\tclassification\tfreq_code("+",".join(aas)+")\t"+"\t".join(aas)+"\tfreq\n"
  for n in composition:
    out=out+sample+"\t"+n.split(":")[0]+"\t"+n.split(":")[1].replace("\t",",")+"\t"+n.split(":")[1]+"\t"+str(composition[n])+"\n"
  fh=open(CDR3_region_composition,"w")
  fh.write(out)
  fh.close()
  return()

def Get_mutational_frequencies_per_classification_per_cluster(sample, annot_file,per_cluster_developmental_classification_file,developmental_classification_mutational_file, developmental_classification_cluster_size_file,developmental_classification_vj_usage_file,reverse_primer_group):
  seqs,freqs,alias = Get_sequences(seq_file)
  fh = open(per_cluster_developmental_classification_file,"r") 
  mutations,gene_combinations,cluster_size,total = {},{},{},0
  mutations_all,cluster_size_all,clusters_done = {},{},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[10] not in clusters_done):
        clusters_done[l[10]] = 1
        classes,total = l[2].split("|"),total+1
        for c in classes:
          if(l[6] !="-1"):
            f = 1#sum(freqs[l[1]])
            #name_mut = c+"\t"+l[6]
            name_mut = c
            if(name_mut in mutations):mutations[name_mut] = mutations[name_mut]+[float(l[6])]*f
            else:mutations[name_mut]=[float(l[6])]*f
            #if(name_mut in mutations):mutations[name_mut] = mutations[name_mut]+f
            #else:mutations[name_mut]=f
            if(c in mutations_all):
              mutations_all[c] = mutations_all[c]+[float(l[6])]*f
              cluster_size_all[c] = cluster_size_all[c]+[int(l[8])]*f
            else:
              mutations_all[c]=[float(l[6])]
              cluster_size_all[c]=[int(l[8])]
            name_gene = c+"\t"+l[4].split("*")[0]+"\t"+l[5].split("*")[0]+"\t"+l[4].split("*")[0]+"|"+l[5].split("*")[0]
            if(name_gene in gene_combinations):gene_combinations[name_gene] = gene_combinations[name_gene]+f
            else:gene_combinations[name_gene] = f
            name_cluster_size = c+"\t"+l[8]
            if(name_cluster_size in cluster_size):cluster_size[name_cluster_size] = cluster_size[name_cluster_size] +f
            else:cluster_size[name_cluster_size]=f
  fh.close()
  Print_developmental_counts(mutations, developmental_classification_mutational_file, total,gene_combinations, developmental_classification_vj_usage_file, cluster_size, developmental_classification_cluster_size_file,mutations_all,cluster_size_all)
  return()

def Get_mutational_frequencies_per_classification(sample, annot_file_internal, developmental_classification_file,developmental_classification_mutational_file, developmental_classification_cluster_size_file,developmental_classification_vj_usage_file,reverse_primer_group):
  seqs,freqs,alias = Get_sequences(seq_file)
  #classification = Get_developmental_classification(developmental_classification_file)
  fh = open(developmental_classification_file,"r")
  mutations,gene_combinations,cluster_size,total = {},{},{},0
  mutations_all,cluster_size_all = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      classes = l[2].split("|")
      total = total+int(l[9])
      for c in classes:
        if(l[6] !="-1"):
          f = sum(freqs[l[1]])
          #name_mut = c+"\t"+l[6]
          name_mut = c
          if(name_mut in mutations):mutations[name_mut] = mutations[name_mut]+[int(l[6])]*f
          else:mutations[name_mut]=[int(l[6])]*f
          if(c in mutations_all):
            mutations_all[c] = mutations_all[c]+[int(l[6])]*f
            cluster_size_all[c] = cluster_size_all[c]+[int(l[9])]*f
          else:
            mutations_all[c]=[int(l[6])]
            cluster_size_all[c]=[int(l[9])]
        name_gene = c+"\t"+l[4].split("*")[0]+"\t"+l[5].split("*")[0]+"\t"+l[4].split("*")[0]+"|"+l[5].split("*")[0]
        if(name_gene in gene_combinations):gene_combinations[name_gene] = gene_combinations[name_gene]+f
        else:gene_combinations[name_gene] = f
        name_cluster_size = c+"\t"+l[9]
        if(name_cluster_size in cluster_size):cluster_size[name_cluster_size] = cluster_size[name_cluster_size] +f
        else:cluster_size[name_cluster_size]=f 
  fh.close()
  Print_developmental_counts(mutations, developmental_classification_mutational_file, total,gene_combinations, developmental_classification_vj_usage_file, cluster_size, developmental_classification_cluster_size_file,mutations_all,cluster_size_all)
  return()

def Print_developmental_counts(mutations, developmental_classification_mutational_file, total,gene_combinations, developmental_classification_vj_usage_file, cluster_size, developmental_classification_cluster_size_file,mutations_all,cluster_size_all):
  out = "#ID\tclassification\tmean number of mutations\tnumber of sequences\ttotal number of reads in repertoire\n"
  for c in mutations:
    out=out+sample+"\t"+c+"\t"+str(mean(mutations[c]))+"\t"+str(len(mutations[c]))+"\t"+str(total)+"\n"
  fh=open(developmental_classification_mutational_file,"w")
  fh.write(out)
  fh.close()
  out = "#ID\tclassification\tV\tJ\tVJ\tnumber of sequences\t% reads/clusters\ttotal reads\clusters\n"
  for c in gene_combinations:
    out=out+sample+"\t"+c+"\t"+str(gene_combinations[c])+"\t"+str(gene_combinations[c]*100.0/total)+"\t"+str(total)+"\n"
  fh=open(developmental_classification_vj_usage_file,"w")
  fh.write(out)
  fh.close()
  out = "#ID\tclassification\tcluster size(%)\tnumber of sequences\n"
  for c in cluster_size:
    s = int(c.split("\t")[1])
    out=out+sample+"\t"+c.split("\t")[0]+"\t"+str(s*100.0/total)+"\t"+str(cluster_size[c])+"\n"
  fh=open(developmental_classification_cluster_size_file,"w")
  fh.write(out)
  fh.close()
  ###########################
  out = "#ID\tclassification\tnumber of unique sequences\tmean mutations\tmean cluster size (%)\t% of repertoire\n"
  for c in mutations_all:
    out=out+sample+"\t"+c+"\t"+str(len(mutations_all[c]))+"\t"+ str(mean(mutations_all[c]))+"\t"+str(mean(cluster_size_all[c])*100.0/total)+"\t"+str(sum(cluster_size_all[c])*100.0/total)+"\n"
  #fh=open(developmental_classification_mutational_file,"w")
  #fh.write(out)
  #f.close()
  return()

def Get_gene_mutation_frequencies(seq_file, annot_file,mutational_gene_freq_file,cluster_file,sample,annot_file_IMGT):
  seqs,freqs,alias = Get_sequences(seq_file)
  #annots = Get_annotation(annot_file)
  CDR3s,mutations = Get_annotations(annot_file_IMGT)
  print(len(mutations))
  gene_mutations = {}
  for id in mutations:
    mut = mutations[id][0]
    f = freqs[id]
    classification = alias[id].split("|")[1].split("_")
    v,j = mutations[id][3].split("\t")[0],mutations[id][3].split("\t")[1]
    #v = v[0:5]#+":"+j[0:5]
    nz = [i for i in range(len(f)) if f[i]!=0]
    name1 = v
    for i in nz:
      c = classification[i].split("*")[0]
      if(c in ["IGHA1","IGHA2"]):c="IGHA1/2"
      if(c in ["IGHG1","IGHG2"]):c="IGHG1/2"
      name =name1+"\t"+c
      if(name in gene_mutations):gene_mutations[name] = gene_mutations[name]+[mut]*f[i]
      else:gene_mutations[name] = [mut]*f[i]
      if(c in ["IGHD","IGHM"]):
        c = "IGHD/M"
        name =name1+"\t"+c
        if(name in gene_mutations):gene_mutations[name] = gene_mutations[name]+[mut]*f[i]
        else:gene_mutations[name] = [mut]*f[i]
  out ='#sample\tv\tclass\tmean number of mutations\tnumber_BCRs\n'
  for name in gene_mutations:
    out=out+sample+"\t"+name+"\t"+str(mean(gene_mutations[name]))+"\t"+str(len(gene_mutations[name]))+"\n"
  fh=open(mutational_gene_freq_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_unexpanded_cluster_proportions(seq_file, sample,unexpanded_cluster_file,cluster_file):
  seqs,freqs,alias = Get_sequences(seq_file)
  fh=open(cluster_file, "r")
  clusters = Tree()
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,cluster = l[2].split("__")[0],l[1]
      clusters[cluster][id].value = 1
  fh.close()
  unexpanded_cluster_props,all_unexpanded_cluster_props,expanded_cluster_props = {},{},{}
  classification = alias[id.split("__")[0]].split("|")[1].split("_")
  classes = []
  threshold_number_of_sequences = 2
  for i in range(0,len(classification)):
    c = classification[i].split("*")[0]
    if(c in ["IGHA1","IGHA2"]):c="IGHA1/2"
    if(c in ["IGHG1","IGHG2"]):c="IGHG1/2"
    classes.append(c)
  for c in clusters:
    total_freqs,per_cluster_freqs = 0,{}
    for id in clusters[c]:
      f = freqs[id]
      total_freqs = total_freqs+sum(f)
      #if(total_freqs>threshold_number_of_sequences):
      #  break
      nz = [i for i in range(len(f)) if f[i]!=0]
      for i in nz:
        if(classes[i].count("P")==0):
          if(classes[i] in per_cluster_freqs):per_cluster_freqs[classes[i]] = per_cluster_freqs[classes[i]]+f[i]
          else:per_cluster_freqs[classes[i]]=f[i]
          if(classes[i] in all_unexpanded_cluster_props):all_unexpanded_cluster_props[classes[i]]=all_unexpanded_cluster_props[classes[i]]+f[i]
          else:all_unexpanded_cluster_props[classes[i]]=f[i]
    if(total_freqs<=threshold_number_of_sequences):
      for c in per_cluster_freqs:
        if(c in unexpanded_cluster_props):unexpanded_cluster_props[c] = unexpanded_cluster_props[c]+per_cluster_freqs[c]
        else:unexpanded_cluster_props[c]=per_cluster_freqs[c]
    else:
      for c in per_cluster_freqs:
        if(c in expanded_cluster_props):expanded_cluster_props[c] = expanded_cluster_props[c]+per_cluster_freqs[c]
        else:expanded_cluster_props[c]=per_cluster_freqs[c]
  total,total_all = 0,0
  for c in unexpanded_cluster_props:
    total = total+unexpanded_cluster_props[c]
  for c in all_unexpanded_cluster_props:
    total_all = total_all+all_unexpanded_cluster_props[c]
  out="#Sample\tisotype\tnumber of reads (unexpanded)\t% reads (unexpanded)\tnumber of reads (all)\t% reads (all)\tnumber of reads (expanded)\t% reads (expanded)\n"
  for c in all_unexpanded_cluster_props:
    if(c not in unexpanded_cluster_props):unexp = 0
    else:unexp = unexpanded_cluster_props[c]
    if(c not in expanded_cluster_props):exp = 0
    else:exp = expanded_cluster_props[c]
    out=out+sample+"\t"+c+"\t"+str(unexp)+"\t"+str(unexp*100.0/total)+"\t"+str(all_unexpanded_cluster_props[c])+"\t"+str(all_unexpanded_cluster_props[c]*100.0/total_all)+"\t"+str(exp)+"\t"+str(exp*100.0/total)+"\n"
  fh=open(unexpanded_cluster_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_cluster_ids_expanded(cluster_file,freqs,alias,threshold):
  fh=open(cluster_file, "r")
  clusters = Tree()
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,cluster = l[2],l[1]
      clusters[cluster][id].value = 1
  fh.close()
  expanded_cluster_ids = {}
  expanded_cluster_isotypes = Tree()
  expanded_clusters = Tree()
  expanded_cluster_all_sequences = Tree()
  classification = alias[id.split("__")[0]].split("|")[1].split("_")
  for c in clusters:
    if(len(clusters[c])>=threshold):
      total_f = 0
      for id in clusters[c]:
        id = id.split("__")[0]
        f = freqs[id]
        total_f=total_f+sum(f)
        nz = [i for i in range(len(f)) if f[i]!=0]
        for i in nz:
          c1 = classification[i].split("*")[0]
          if(c1 in ["IGHA1","IGHA2"]):c1="IGHA1/2"
          if(c1 in ["IGHG1","IGHG2"]):c1="IGHG1/2"
          name = c+"\t"+c1
          expanded_cluster_isotypes[name][f[i]][id+":"+str(i)].value = 1
        expanded_cluster_all_sequences[c][id].value = 1
        expanded_cluster_ids[id] = c
      expanded_clusters[c] = total_f
  return(expanded_cluster_ids, clusters, expanded_clusters, expanded_cluster_isotypes,expanded_cluster_all_sequences,)

def Number_of_unique_clusters(seq_file,annot_file_IMGT, unique_cluster_isotypes_file, sample):
  seqs,freqs,alias = Get_sequences(seq_file)
  CDR3s,mutations1 = Get_annotations(annot_file_IMGT)
  for id in alias:
    break
  classification = alias[id.split("__")[0]].split("|")[1].split("_")
  clusters = {}
  for CDR in CDR3s:
    classes_used = {}
    for id in CDR3s[CDR]:
      f = freqs[id]
      nz = [i for i in range(len(f)) if f[i]!=0]
      for i in nz:
        c1 = classification[i].split("*")[0]
        if(c1 in ["IGHA1","IGHA2"]):c1="IGHA1/2"
        if(c1 in ["IGHG1","IGHG2"]):c1="IGHG1/2"
        classes_used[c1] = 1
    for c in classes_used:
      if(c1 in clusters):clusters[c1] = clusters[c1]+1
      else:clusters[c1]=1
  total_clusters = len(CDR3s)
  out="#Sample\tisotype\tNumber_of_clusters\tPercentage_clusters\tTotal_clusters\n"
  for c in clusters:
    if(c.count("P")==0):
      out=out+sample+"\t"+c+"\t"+str(clusters[c])+"\t"+str(clusters[c]*100.0/total_clusters)+"\t"+str(total_clusters)+"\n"
  fh=open(unique_cluster_isotypes_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_mutational_positions_per_gene(seq_file, developmental_classification_file, mutational_positions_per_gene_file,annot_file7,sample,mutational_positions_per_gene_file_summary):
  fh=open(developmental_classification_file,"r")
  classification = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      classification[l[1]] = l[2].split("|")
  fh.close()
  print(len(classification))
  mutational_variation = {}
  mutational_variation_summary = Tree()
  fh=open(annot_file7,"r")
  for l in fh:
    l=l.strip().split("\t")
    if(l[1].split("__")[0] in classification):
      if(len(l)>=5):
        classes = classification[l[1].split("__")[0]]
        classes = classes+["ALL"]
        v = l[3].split()[1].split("*")[0]
        v = v[0:5]
        muts = l[4].split("|")
        for c in classes:
          if(c.count(",")==0):
            if(c in ["IGHA1","IGHA2","IGHA"]):c = "IGHA1/2"
            if(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
            for m in muts: 
              if(m.count(",")!=0):
                m = m.split(",")[1].split("(")[0].replace(">","\t")
                m=m[0]+"\t"+m[1:len(m)]
                m1 = m.split("\t")
                m1=[m1[0]+"\t"+m1[1],m1[2]]
                name = c+"\t"+v+"\t"+m
                if(name in mutational_variation):mutational_variation[name] = mutational_variation[name]+1
                else:mutational_variation[name] = 1
                name = c+"\t"+v+"\t"+m1[0]
                mutational_variation_summary[name][m1[1]].value = 1
  fh.close()
  for f in [mutational_positions_per_gene_file_summary,mutational_positions_per_gene_file]:
    fh=open(f,"w")
    fh.close()
  out = "#sample\tIsotype\tV_gene\tAA\tPosition\tAA_variant_number\n"
  for c in mutational_variation_summary:
    out=out+sample+"\t"+c+"\t"+str(len(mutational_variation_summary[c]))+"\n"
  Write_out(out, mutational_positions_per_gene_file_summary)
  out = "#sample\tIsotype\tV_gene\tAA\tPosition\tAA_variant\tfrequency\n"
  ind = 0
  for c in mutational_variation:
    out=out+sample+"\t"+c+"\t"+str(mutational_variation[c])+"\n"
    ind = ind+1
    if(ind>300):
      Write_out(out, mutational_positions_per_gene_file)
      out,ind = '',0
  Write_out(out, mutational_positions_per_gene_file)
  return()

def Get_mutations_and_switching_of_expanded_clusters_specific_isotypes(seq_file, expanded_cluster_annotation, cluster_file, annot_file2,reverse_primer_group):
  fh=open(expanded_cluster_annotation,"w")
  fh.close()
  seqs,freqs,alias = Get_sequences(seq_file)
  threshold = 3
  CDR3s,mutations1 = Get_annotations(annot_file)
  expanded_cluster_ids, clusters,expanded_clusters, expanded_cluster_isotypes,expanded_cluster_all_sequences = Get_cluster_ids_expanded(cluster_file,freqs,alias, threshold)
  isotypes_of_interest = ["IGHE","IGHA1/2", "IGHA1","IGHA2"]
  for id in alias:
    break
  classification = alias[id.split("__")[0]].split("|")[1].split("_")
  out = '#sample\tcluster\tisotype(s)\tmutations\tfrequency\tnumber_of_isotypes\tcluster_unique_sequences\n'
  ind=  0
  for c in expanded_cluster_all_sequences:
    done = 0
    for iso in isotypes_of_interest:
      name = c+"\t"+iso
      if(done==1):break
      if(name in expanded_cluster_isotypes):
        done = 1
        sequence_types = {}
        total_isotypes = {}
        for id in expanded_cluster_all_sequences[c]:
          if(id in mutations1):
            muts = mutations1[id]
            muts = muts[0]+muts[1]
            classes = []
            f = freqs[id]
            nz = [i for i in range(len(f)) if f[i]!=0]
            for i in nz:
              c1 = classification[i].split("*")[0]
              if(reverse_primer_group !="ISOTYPER"):
                if(c1 in ["IGHA1","IGHA2"]):c1="IGHA1/2"
                if(c1 in ["IGHG1","IGHG2"]):c1="IGHG1/2"
              if(c1 not in classes):classes.append(c1)
              total_isotypes[c1] = 1
            classes.sort()
            for c1 in classes:
              name = c1+"\t"+str(muts)
              if(name in sequence_types):sequence_types[name] = sequence_types[name]+1
              else:sequence_types[name]=1
            classes = ",".join(classes)
            name = classes+"\t"+str(muts)
            if(name in sequence_types):sequence_types[name] = sequence_types[name]+1
            else:sequence_types[name]=1
          else:
            print("ERROR",id)
        #print c, name
        if(len(sequence_types)>=2):
          number_isotypes = len(total_isotypes)
          for name in sequence_types:
            out=out+sample+"\t"+c+"\t"+name+"\t"+str(sequence_types[name])+"\t"+str(number_isotypes)+"\t"+str(len(expanded_cluster_all_sequences[c]))+"\n"
            ind = ind+1
            if(ind>500):
              Write_out(out,expanded_cluster_annotation) 
              out, ind = '',0
  Write_out(out,expanded_cluster_annotation)
  out, ind = '',0
  return()

def Get_mutations_per_expanded_cluster(seq_file, annot_file2, sample,mutations_per_expanded_cluster_file,cluster_file,reverse_primer_group):
  seqs,freqs,alias = Get_sequences(seq_file)
  threshold = 3
  expanded_cluster_ids, clusters,expanded_clusters, expanded_cluster_isotypes,expanded_cluster_all_sequences = Get_cluster_ids_expanded(cluster_file,freqs,alias,threshold)
  fh=open(annot_file2,"r")
  regions = ['FR1','CDR1','FR2','CDR2','FR3','CDR3']
  n_nucleotides = [23,41,59,77,95,113]
  col_total =     [25,43,61,79,97,115]
  col_silent =    [26,44,62,80,98,116]
  col_nonsilent = [27,45,63,81,99,117]
  cols_all = [n_nucleotides,col_total,col_silent,col_nonsilent]
  m_CDR_mm,m_FWR_mm,m_CDR_FWR_ratio,m_silent,m_nonsilent,m_silent_nonsilent_ratio = {},{},{},{},{},{}
  counts_all = [m_CDR_mm,m_FWR_mm,m_CDR_FWR_ratio,m_silent,m_nonsilent,m_silent_nonsilent_ratio]
  count_type = ['mean CDR_mm per BCR','mean FWR_mm per BCR','mean CDR_FWR_ratio','mean silent per BCR','mean nonsilent per BCR','mean silent_nonsilent_ratio']
  mutation_information = {}
  for l in fh:
    l=l.strip().split("\t")
    id = l[1].split("__")[0]
    if(l[2]=="productive" and len(l) >= 130 and id in expanded_cluster_ids):
      array = []
      for i in cols_all:
        array1 = []
        for j in i:
          array1 = array1 + [int(l[j].split()[0])]
        array = array+[array1]
      array = np.array(array)
      if(sum(array[1,])>=0):
        FWR_mm = sum(array[1,[0,2,4]])
        CDR_mm = sum(array[1,[1,3,5]])
        FWR_len, CDR3_len = sum(array[0,[0,2,4]]), sum(array[0,[1,3,5]])
        silent_mm, non_silent_mm = sum(array[2,]), sum(array[3,])
        total_mm = int(l[7].split()[0])
        info = [FWR_mm, CDR_mm, silent_mm, non_silent_mm,total_mm]
        mutation_information[id] = info
  fh.close()
  count_names = ['FWR_mm', 'CDR_mm', 'silent_mm', 'non_silent_mm','total_mm']
  out="#sample\tisotype\tnumber_of_clusters\tnumber_of_reads_isotype\tnumber_reads_total_cluster\tmean_mutations\ttotal_silent\ttotal_non_silent\ttotal_CDR\ttotal_FWR\n"
  isotypes_information = {}
  for c1 in expanded_cluster_isotypes:
    c, clas = c1.split("\t")[0], c1.split("\t")[1]
    total_cluster_size = expanded_clusters[c]
    isotype_cluster_size = 0
    ids = {}
    for f1 in expanded_cluster_isotypes[c1]:
      isotype_cluster_size = isotype_cluster_size+(f1*len(expanded_cluster_isotypes[c1][f1]))
      for id in expanded_cluster_isotypes[c1][f1]:
        if(id.split(":")[0] in mutation_information):
          ids[id.split(":")[0]]=f1
    if(len(ids)>=3):
      values,totals = {},0
      for id in ids:
        for i in range(0,len(count_names)):
          if(count_type[i] in values):values[count_names[i]] = values[count_names[i]] + mutation_information[id][i]*ids[id]
          else:values[count_names[i]]=mutation_information[id][i]*ids[id]
        totals = totals+ids[id]
      values_print = []
      for i in range(0,len(count_names)):
        values_print.append(values[count_names[i]]*1.0/totals)
      info = [isotype_cluster_size, total_cluster_size,values['total_mm'],values['silent_mm'],values['non_silent_mm'], values['CDR_mm'], values['FWR_mm']]
      if(clas in isotypes_information):
        isotypes_information[clas] = isotypes_information[clas]+[info]
      else:isotypes_information[clas]=[info]
  names_isotypes_information = ['isotype_cluster_size','total_cluster_size','total_mm','silent_mm','non_silent_mm','CDR_mm','FWR_mm']
  for clas in isotypes_information:
    values_print = []
    for i in range(0,len(names_isotypes_information)):
      values = []
      for info in isotypes_information[clas]:
        values.append(info[i])
      if(names_isotypes_information[i] in ['total_mm','silent_mm','non_silent_mm','CDR_mm','FWR_mm']):
        values = mean(values)
      else:values = sum(values)
      values_print.append(values)
    out=out+sample+"\t"+clas+"\t"+str(len(isotypes_information[clas]))+"\t"+"\t".join(map(str, values_print))+"\n"
  fh=open(mutations_per_expanded_cluster_file,"w")
  fh.write(out)
  fh.close()
  return()

def Indel_analysis_per_chain(seq_file, sample, indel_file,annot_file,reverse_primer_group):
  seqs,freqs,alias = Get_sequences(seq_file)
  print(annot_file)
  fh=open(annot_file,"r")
  indels_position,totals = {},{}
  for l in fh:
    l1 = l
    l=l.strip().split("\t")
    if(l[2].count("productive")!=0):
      id = l[1].split("__")[0]
      f = freqs[id]
      classification = alias[id].split("|")[1].split("_")
      nz = [i for i in range(len(f)) if f[i]!=0]
      classes = {}
      for i in nz:
        c = classification[i].split("*")[0]
        if(reverse_primer_group !="ISOTYPER"):
          if(c in ["IGHA1","IGHA2"]):c="IGHA1/2"
          if(c in ["IGHG1","IGHG2"]):c="IGHG1/2"
        classes[c] = 1
      for c in classes:
        if(c in totals):totals[c] = totals[c]+1
        else:totals[c] = 1
      indel_events, V_insertions, V_deletions = l[8], l[26],l[27]
      if(len(indel_events)!=0 or l1.count("frameshift")!=0):
        if(len(V_insertions)!=0):
          region = V_insertions.split()[1].replace("-IMGT",'')
          codon = V_insertions.split("codon ")[1].split()[0]
          nucleotides = V_insertions.split("(")[1].split()[0]
          frameshift = V_insertions.count('(do not cause frameshift)')
          if(frameshift==1):frameshift = "NO_FRAMESHIFT"
          else:frameshift = "FRAMESHIFT"
          #types = [region, region+",CODON:"+codon, region+",LENGTH:"+nucleotides, region+","+frameshift]
          types = [region, region+","+frameshift, region+",CODON:"+codon]
          for t in types:
            for c1 in classes:
              c = c1+"\t"+"INSERTION:"+t
              if(c in indels_position):indels_position[c] = indels_position[c]+1
              else:indels_position[c] = 1
        if(len(V_deletions)!=0):
          region = V_deletions.split()[1].replace("-IMGT",'').replace(",","")
          codon = V_deletions.split("codon ")[1].split()[0]
          nucleotides = V_deletions.split("V-REGION: ")[1].split()[0]
          frameshift = V_deletions.count('(do not cause frameshift)')
          if(frameshift==1):frameshift = "NO_FRAMESHIFT"
          else:frameshift = "FRAMESHIFT"
          #types = [region, region+",CODON:"+codon, region+",LENGTH:"+nucleotides, region+","+frameshift]
          types = [region, region+","+frameshift, region+",CODON:"+codon]
          for t in types:
            for c1 in classes:
              c = c1+"\t"+"DELETION:"+t
              if(c in indels_position):indels_position[c] = indels_position[c]+1
              else:indels_position[c] = 1
  fh.close()
  out= '#sample\tisotype\tindel_type\ttotal_BCRs_of_isotype\tfrequency of indel\tpercentage indel\n'
  for t in indels_position:
    c = t.split()[0]
    out = out+"\t".join(map(str, [sample,t,totals[c], indels_position[t], indels_position[t]*100.0/totals[c]]))+"\n"
  fh=open(indel_file,"w")
  fh.write(out)
  fh.close()
  return()

def Mutation_selection_per_chain(seq_file, annot_file2, sample, syn_non_syn_mutations,syn_non_syn_mutations_summary,reverse_primer_group):
  seqs,freqs,alias = Get_sequences(seq_file)
  fh=open(annot_file2,"r")
  print(annot_file2)
  regions = ['FR1','CDR1','FR2','CDR2','FR3','CDR3']
  n_nucleotides = [23,41,59,77,95,113]
  col_total =     [25,43,61,79,97,115]
  col_silent =    [26,44,62,80,98,116]
  col_nonsilent = [27,45,63,81,99,117]
  cols_all = [n_nucleotides,col_total,col_silent,col_nonsilent]
  mutations,ids_done = {},{}
  total,total_sub = 0,0
  m_CDR_mm,m_FWR_mm,m_CDR_FWR_ratio,m_silent,m_nonsilent,m_silent_nonsilent_ratio,FWR3 = {},{},{},{},{},{},{}
  counts_all = [m_CDR_mm,m_FWR_mm,m_CDR_FWR_ratio,m_silent,m_nonsilent,m_silent_nonsilent_ratio, FWR3]
  count_type = ['mean CDR_mm per BCR','mean FWR_mm per BCR','mean CDR_FWR_ratio','mean silent per BCR','mean nonsilent per BCR','mean silent_nonsilent_ratio',"FWR3_mm"]
  print(len(freqs))
  for l in fh:
    l=l.strip().split("\t")
    if(l[2].count("productive")!=0 and len(l) >= 130):
      array = []
      for i in cols_all:
        array1 = []
        for j in i:
          array1 = array1 + [int(l[j].split()[0])]
        array = array+[array1]
      array = np.array(array)
      if(sum(array[1,])>0):
        id = l[1].split("__")[0]
        f = freqs[id]
        #total = total+sum(f)
        classification = alias[id].split("|")[1].split("_")
        nz = [i for i in range(len(f)) if f[i]!=0]
        FWR_mm = sum(array[1,[0,2,4]])
        CDR_mm = sum(array[1,[1,3,5]])
        FWR_len, CDR3_len = sum(array[0,[0,2,4]]), sum(array[0,[1,3,5]])
        fwr3 = sum(array[1,[4]])
        silent_mm, non_silent_mm = sum(array[2,]), sum(array[3,])
        info = [FWR_mm, CDR_mm, FWR_len, CDR3_len, silent_mm, non_silent_mm]
        ratio1,ratio2 = -1,-1
        if((CDR_mm+FWR_mm)>0):ratio1 = (CDR_mm)*1.0/(CDR_mm+FWR_mm)
        if((silent_mm+non_silent_mm)>0):ratio2 =(silent_mm)*1.0/(silent_mm+non_silent_mm)
        info1 = [CDR_mm, FWR_mm, ratio1, silent_mm, non_silent_mm,ratio2,fwr3]
        #info = "\t".join(map(str, info))
        for i in nz:
          c = classification[i].split("*")[0]
          if(reverse_primer_group !="ISOTYPER"):
            if(c in ["IGHA1","IGHA2"]):c="IGHA1/2"
            if(c in ["IGHG1","IGHG2"]):c="IGHG1/2"
          #name = c+"\t"+info
          #if(name in mutations):mutations[name] = mutations[name]+f[i]
          #else:mutations[name]=f[i]
          #total_sub = total_sub+f[i]
          c1 = [c]
          if(CDR_mm+FWR_mm>=2  and c in ["IGHD","IGHM"]):c1 = c1+["IGHD_IGHM_mutated"]
          for c in c1:
            for i1 in range(0,len(info1)):
              if(info1[i1]!=-1):
                name = c
                if(name in counts_all[i1]):
                  counts_all[i1][name] = [counts_all[i1][name][0]+info1[i1]*1, counts_all[i1][name][1]+1]
                else:counts_all[i1][name]=[info1[i1]*1, 1]
          #print len(counts_all[i1][name])
          #if(len(counts_all[i1][name])>100):break
  fh.close()
  out = "#sample\tcount_type\tchain\tmean value\ttotal_BCRs_counted\n"
  for i1 in range(0,len(counts_all)):
    X = counts_all[i1]
    for c in X:
      mean1, length = X[c][0]*1.0/X[c][1], X[c][1]
      out=out+sample+"\t"+count_type[i1]+"\t"+c+"\t"+str(mean1)+"\t"+str(length)+"\n"
  fh=open(syn_non_syn_mutations_summary,"w")
  fh.write(out)
  fh.close()
  return()

def Get_proportion_non_functional_reads(seq_file, non_functional_rearrangments_count_file, non_ORF_filtered_reads,sample):
  total1, total2 = {},{}
  fh = open(non_ORF_filtered_reads, "r")
  for header,sequence in fasta_iterator(fh):
    freq, classification = list(map(int, header.split("__")[1].split("|")[0].split("_"))), header.split("|")[1].split("_")
    nz = [i for i in range(len(freq)) if freq[i]!=0]
    for i in nz:
      c = classification[i].split("*")[0]
      if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
      elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
      if(c in total1):total1[c] = total1[c]+freq[i]
      else:total1[c]=freq[i]
  fh.close()
  fh = open(seq_file,"r")
  for header,sequence in fasta_iterator(fh):
    freq, classification = list(map(int, header.split("__")[1].split("|")[0].split("_"))), header.split("|")[1].split("_")
    nz = [i for i in range(len(freq)) if freq[i]!=0]
    for i in nz:
      c = classification[i].split("*")[0]
      if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
      elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
      if(c in total2):total2[c] = total2[c]+freq[i]
      else:total2[c]=freq[i]
  out="#sample\tchain\treads_pre_ORF_filtering\treads_post_ORF_filtering\tpercentage_non_functional\n"
  for c in total1:
    f = 0
    if(c in total2):f = total2[c]
    out=out+sample+"\t"+c+"\t"+str(total1[c])+"\t"+str(f)+"\t"+str(f*100.0/total1[c])+"\n"
  fh=open(non_functional_rearrangments_count_file,"w")
  fh.write(out)
  fh.close()
  return()

def Get_clonality_above_n_mutations(seq_file, annot_file_IMGT, clonality_mutation_correlation,cluster_file,sample):
  for f in [clonality_mutation_correlation]:
    fh=open(f,"w")
    fh.write("#Id\tMin_n_mutations\tIsotype\tN reads\tN vertices\tVertex Gini Index\tunique sequences in subset\ttotal_reads\n")
    fh.close()
  seqs,freqs,alias = Get_sequences(seq_file)
  CDR3s,mutations1 = Get_annotations(annot_file_IMGT)
  for mm in range(0,35):
    chain_counts, total = {},0
    for id in mutations1:
      if(mutations1[id][0]>=mm):
        freq,chains = freqs[id],alias[id].split("|")[1].split("_")
        nz = [i for i in range(len(freq)) if freq[i]!=0]
        total = total+sum(freq)
        f = {}
        for i in nz:
          c= chains[i]
          c = c.split("*")[0]
          if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
          elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
          if(c in f):f[c] = f[c]+freq[i]
          else:f[c]=freq[i]
        for c in f:
          if(c in chain_counts):chain_counts[c]=chain_counts[c]+[f[c]]
          else:chain_counts[c]=[f[c]]
        c = "ALL"
        if(c in chain_counts):chain_counts[c]=chain_counts[c]+[sum(freq)]
        else:chain_counts[c]=[sum(freq)]
    out=''
    for c in chain_counts:
      (vpoints,vvdf)=VDF(chain_counts[c])
      (vgini)=Get_Gini(vpoints,vvdf)
      out=out+sample+"\t"+str(mm)+"\t"+c+"\t"+str(sum(chain_counts[c]))+"\t"+str(len(chain_counts[c]))+"\t"+str(vgini)+"\t"+str(len(chain_counts[c]))+"\t"+str(total)+"\n"
    Write_out(out, clonality_mutation_correlation)
    print("\t",sample,mm)
  return()

def Get_per_cluster_per_isotype_degree(seq_file, edge_file, cluster_file, per_cluster_per_isotype_degree_file,sample,annot_file_IMGT):
  fh=open(per_cluster_per_isotype_degree_file,"w")
  fh.close()
  clusters, inv_clusters,gl_mutations = Get_cluster_occupancy(cluster_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  CDR3s,mutations1 = Get_annotations(annot_file_IMGT)
  del seqs, CDR3s
  connections = Tree()
  fh = open(edge_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id1,id2 = l[0].split("__")[0],l[1].split("__")[0]
      connections[id1][id2].value = 1
      connections[id2][id1].value = 1
  fh.close()
  cluster_sizes = {}
  for id in alias:
    length = len(alias[id])
    chains =alias[id].split("|")[1].split("_")
    break
  for c in clusters:
    f = 0#[0]*length
    for id in clusters[c]:
      f = f+sum(freqs[id])#map(add, f, freqs[id])
    cluster_sizes[c] = f
  cluster_degree_freq = {}
  for id in freqs:
    if(id in mutations1):
      freq = freqs[id]
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      degree = len(connections[id])
      mut = mutations1[id][0]+mutations1[id][1]
      size = cluster_sizes[inv_clusters[id]]
      if(size>50):size = ">50"
      elif(size>20):size = "20-50"
      elif(size>10):size = "10-20"
      elif(size>5):size = "5-10"
      elif(size>1):size = "1-5"
      else:size = "1"
      n = sample+"\t"+str(degree)+"\t"+str(mut)+"\t"+str(size)
      for i in nz:
        c= chains[i].split("*")[0]
        if(c in ["IGHA1","IGHA2"]):c = "IGHA1/2"
        elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        name = n+"\t"+c
        if(name in cluster_degree_freq):cluster_degree_freq[name] = cluster_degree_freq[name]+1
        else:cluster_degree_freq[name]=1
  del connections
  out,ind = "#sample\tdegree\tmutations\tcluster_size\tchain\tnumber_of_sequences\n",0
  for name in cluster_degree_freq:
    out=out+name+"\t"+str(cluster_degree_freq[name])+"\n"
    ind = ind+1
    if(ind>300):
      Write_out(out, per_cluster_per_isotype_degree_file)
      out, ind = '',0
  Write_out(out, per_cluster_per_isotype_degree_file)
  return()

def Get_unmutated_repertoire(seq_file, annot_file_IMGT, seq_unmutated_file,clusters_unmutated_file,cluster_file): 
  CDR3s,mutations1 = Get_annotations(annot_file_IMGT)
  seqs,freqs,alias = Get_sequences(seq_file)
  out, ind,out1 = '',0,''
  for f in [seq_unmutated_file,clusters_unmutated_file]:
    fh=open(f,"w")
    fh.close()
  ids_use = {}
  total,filtered = 0,0
  for id in mutations1:
    muts= [mutations1[id][0],mutations1[id][1]]
    freq,chains = freqs[id],alias[id].split("|")[1].split("_")
    total = total+sum(freq)
    if(muts[0]<=3):# and muts[1]<=2):
      filtered = filtered+sum(freq)
      id_full = alias[id]
      nz = [i for i in range(len(freq)) if freq[i]!=0]
      passes = 0
      for i in nz:
        c = chains[i].split("*")[0]
        if(c in ["IGHM","IGHD"]):
          passes = 1
          break
      if(passes ==1):
        ids_use[id]=1
        out=out+">"+id_full+"\n"+seqs[id]+"\n"
        ind = ind+1
        if(ind>100):
          Write_out(out, seq_unmutated_file)
          out, ind = '',0
    else:
      print(freq, muts)
  Write_out(out, seq_unmutated_file)
  fh=open(cluster_file,"r")
  ind = 0
  for l in fh:
    if(l[0]!="#"):
      l1=l.strip().split()
      if(l1[2].split("__")[0] in ids_use):
        out1 = out1+l
        ind = ind+1
        if(ind>100):
          Write_out(out1,clusters_unmutated_file)
          out1,ind = '',0
  fh.close()
  Write_out(out1,clusters_unmutated_file)
  print("FILTERING\t"+str(total)+"\t"+str(filtered))
  return()

def Get_example_sequence(seq_file, annot_file):
  seqs,freqs,alias = Get_sequences(seq_file)
  CDR3s,mutations1 = Get_annotations(annot_file)
  for id in mutations1:
    if( mutations1[id][0]+ mutations1[id][1]==0):
      if(sum(freqs[id])>=5):
        if(len(seqs[id])>370):
          print(mutations1[id], id, sum(freqs[id]), freqs[id])
          print(len(seqs[id]))
          print(seqs[id])
  return()

def Get_mutational_frequencies(seq_file, annot_file,mutational_freq_file,cluster_file,sample,annot_file_IMGT,reverse_primer_group,n_unmutated_file):
  #clusters, inv_clusters,gl_mutations = Get_cluster_occupancy(cluster_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  #annots = Get_annotation(annot_file)
  CDR3s,mutations1 = Get_annotations(annot_file_IMGT)
  fh=open(mutational_freq_file,"w")
  fh.close()
  mutations,max_mut = {},0
  number_unmutated_BCRs,number_BCRs,mutated_20={},{},{}
  mean_mutations = {}
  #for id in annots:
  for id in mutations1:
    if(id in freqs):
      f = freqs[id]
      classification = alias[id].split("|")[1].split("_")
      #mut = sum(gl_mutations[id])
      #mut = annots[id][2]
      mut= mutations1[id][0]
      f = f+[sum(f)]
      classification = classification+["All"]
      if(mut>max_mut):max_mut = mut
      for i in range(0,len(f)):
        if(f[i]!=0):
          clas = classification[i].split("*")[0]
          if(reverse_primer_group !="ISOTYPER"):
            if(clas in ["IGHA1","IGHA2"]):clas = "IGHA1/2"
            if(clas in ["IGHG1","IGHG2"]):clas = "IGHG1/2"
          if(clas.count("P")==0):
            if(clas in number_BCRs):
              number_BCRs[clas] = number_BCRs[clas]+f[i]
              mean_mutations[clas] = mean_mutations[clas] +[mut]
            else:
              number_BCRs[clas]=f[i]
              mean_mutations[clas] =[mut]
            if(mut<=2):
              if(clas in number_unmutated_BCRs):number_unmutated_BCRs[clas] = number_unmutated_BCRs[clas]+f[i]
              else:number_unmutated_BCRs[clas]=f[i]
            if(mut>10):
              if(clas in mutated_20):mutated_20[clas] = mutated_20[clas] +f[i]
              else:mutated_20[clas] =f[i]
            #clas = clas+"\t"+str(mut)
            #if(clas in mutations):mutations[clas] = mutations[clas]+f[i]
            #else:mutations[clas]=f[i]
  out = "#Sample\tisotype\ttotal_sequences\tperc_unumtated\tperc_mutated>10bp\tmean mutations\tnumber of unique sequences\n"
  for c in number_BCRs:
    if(c in number_unmutated_BCRs):n_unmut = number_unmutated_BCRs[c]
    else:n_unmut = 0
    if(c in mutated_20):mut20 = mutated_20[c]
    else:mut20 = 0
    out=out+sample+"\t"+c+"\t"+str(number_BCRs[c])+"\t"+str(n_unmut*100.0/number_BCRs[c])+"\t"+str(mut20*100.0/number_BCRs[c])+"\t"+str(mean(mean_mutations[c]))+"\t"+str(len(mean_mutations[c]))+"\n"
  fh=open(n_unmutated_file,"w")
  fh.write(out)
  fh.close()
  #classes = []
  #for i in range(0,len(classification)):
  #  if(classification[i].split("*")[0] not in classes):
  #    classes.append(classification[i].split("*")[0])
  #out = "#ID\tMutations\t"+"\t".join(classes)+"\n"
  #out = "#Sample\tisotype\tmutations\tfrequency\n"
  #for c in mutations:
  #  out=out+sample+"\t"+c+"\t"+str(mutations[c])+"\n"
#  for i in range(0,max_mut+1):
#    out=out+sample+"\t"+str(i)
#    for clas in classes:
#      if(clas in mutations):
#        out=out+"\t"+str(mutations[clas].count(i))
#      else:out=out+"\t0"
#    out=out+"\n"
  #fh=open(mutational_freq_file,"w")
  #fh.write(out)
  #fh.close()
  return()

def Make_trees_from_clusters(cluster_annotation_file,merged_cluster_file,seq_file,overall_cluster_large_file,id_sample, output_dir):
  fh=open(cluster_annotation_file,"r")
  cluster_annot,done = {},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      cluster, freqs = l[0], list(map(int, l[1].split(",")))
      if(cluster not in done and cluster!="ALL"):
        done[cluster] = 1
        non_zero = [i for i in range(len(freqs)) if freqs[i]!=0]
        if(len(non_zero)>=3 and sum(freqs)>30):
          cluster_annot[cluster] = l[1]
  fh.close()
  clusters, inv_clusters,gl_mutations = Get_cluster_occupancy(merged_cluster_file)
  seqs,freqs,alias = Get_sequences(seq_file)
  out_all='#Cluster ID\tSequence frequencies\n'
  for c in cluster_annot:
    out, ind = '',0
    file = output_dir+"Unaligned_cluster_sequences_"+id_sample+"_"+c+".txt"
    for id in clusters[c]:
      out=out+">"+alias[id]+"\n"+seqs[id]+"\n"
    fh=open(file, "w")
    fh.write(out)
    fh.close
    out_all = out_all+c+"\t"+cluster_annot[c]+"\n"
  fh=open(overall_cluster_large_file,"w")
  fh.write(out_all)
  fh.close()
  return()

def Get_all_chain_sequences_separately(seq_file, separated_sequence_files):
  seqs,freqs,alias = Get_sequences(seq_file)
  for id in alias:
    chains= alias[id].split("|")[1].split("_")
    break
  chains_all = []
  for i in range(0,len(chains)):
    chains[i] = chains[i].split("*")[0]
    if(chains[i].split("*")[0] not in chains_all):chains_all = chains_all+[chains[i].split("*")[0]]
  chains_all.sort()
  files,outs = [],{}
  for i in range(0,len(chains)):
    outs[chains[i]] = ''
    files.append(separated_sequence_files+chains[i]+".fasta")
    fh=open(separated_sequence_files+chains[i]+".fasta","w")
    fh.close()
  for id in seqs:
    f1 = freqs[id]
    nz = [i for i in range(len(f1)) if f1[i]!=0]
    for i in nz:
      outs[chains[i]] = outs[chains[i]]+">"+id+"\n"+seqs[id]+"\n"
      if(outs[chains[i]].count(">")>100):
        Write_out(outs[chains[i]], files[i])
        outs[chains[i]]=''
  for i in range(0,len(chains)):
    Write_out(outs[chains[i]], files[i])
    outs[chains[i]]=''
  return()

def Get_gene_distribution_isotype_specific(seq_file, annot_file,isotype_specific_gene_frequency_file,sample,annot_file_IMGT):
  gene_combintation = Tree()
  seqs,freqs,alias = Get_sequences(seq_file)
  CDR3s,mutations = Get_annotations(annot_file_IMGT)
  #annots = Get_annotation(annot_file)
  for id in alias:
    chains= alias[id].split("|")[1].split("_")
    break
  chains_all = []
  for i in range(0,len(chains)):
    chains[i] = chains[i].split("*")[0]
    if(chains[i].split("*")[0] not in chains_all):chains_all = chains_all+[chains[i].split("*")[0]]
  chains_all.sort()
  for id in mutations:
    if(id in freqs):
      f1 = freqs[id]
      nz = [i for i in range(len(f1)) if f1[i]!=0]
      vj =  mutations[id][3]
      #vj = annots[id][0].split("*")[0]+"\t"+annots[id][1].split("*")[0]
      for i in nz:
        c = chains[i].split("*")[0]
        if(c in ["IGHA1", "IGHA2"]):c = "IGHA1/2"
        if(c in ["IGHG1", "IGHG2"]):c = "IGHG1/2"
        gene_combintation[vj][c][id][f1[i]].value = 1
        if(c in ["IGHD","IGHM"]):
          c = "IGHD/M"
          muts = mutations[id][0]+mutations[id][1]
          if(muts>4):c = c+"_MUTATED"
          else:c = c+"_UNMUTATED"
          gene_combintation[vj][c][id][f1[i]].value = 1
  out="#Sample\tV\tJ\tchain\tN_reads\tVJ\n"
  for vj in gene_combintation: 
    for chain in gene_combintation[vj]:
      f = 0
      for id in gene_combintation[vj][chain]:
        for f1 in gene_combintation[vj][chain][id]:
          f = f+f1
      out=out+sample+"\t"+vj+"\t"+chain+"\t"+str(f)+"\t"+vj.replace("\t",":")+"\n"
    #for chain in chains_all:
    #  if(chain in gene_combintation[vj]):
    #    f = 0
    #    for id in gene_combintation[vj][chain]:
    #      for f1 in gene_combintation[vj][chain][id]:
    #        f = f+f1
    #    out=out+sample+"\t"+vj+"\t"+chain+"\t"+str(f)+"\t"+vj.replace("\t",":")+"\n"
    #out=out+"\t"+vj.replace("\t",":")+"\n"
  fh=open(isotype_specific_gene_frequency_file, "w")
  fh.write(out)
  fh.close()
  return()

def Annotating_clusters_for_analysis(merged_cluster_file, cluster_annotation_file, annot_file,seq_file):
  fh=open(cluster_annotation_file,"w")
  fh.close()
  seqs,freqs,alias = Get_sequences(seq_file)
  clusters, inv_clusters = Get_cluster_occupancy(merged_cluster_file)
  CDR3s,mutations = Get_annotations(annot_file)
  chains,total = '',''
  out,ind = '',0
  mutations_all = list(range(100))
  out = "#Cluster ID\tcluster_totals\tChain\tV gene\tJ gene\t"+"\t".join(map(str, mutations_all))+"\n"
  for c in clusters: 
    f,muts = {},{}
    vj = ''
    for id in clusters[c]:
      if(id in mutations):
        if(len(chains)==0):chains = alias[id].split("|")[1].split("_")
        if(len(f)==0):f = freqs[id]
        else:f = list(map(add, f, freqs[id]))
        mut = mutations[id][0]+mutations[id][1]
        vj = mutations[id][3]
        f1 = freqs[id]
        nz = [i for i in range(len(f1)) if f1[i]!=0]
        for i in nz:
          names = chains[i]
          if(names in muts):muts[names][mut] = muts[names][mut]+f1[i]
          else:
            muts[names] = [0]*len(mutations_all)
            muts[names][mut] = muts[names][mut]+f1[i]
    for chain in muts:
      out=out+c+"\t"+",".join(map(str,f))+"\t"+chain+"\t"+vj+"\t"+"\t".join(map(str, muts[chain]))+"\n"
      ind = ind+1
      if(ind>100):
        Write_out(out, cluster_annotation_file)
        out ,ind = '',0
    if(len(total)==0):total = f
    else:total =list(map(add,f, total))
  for chain in chains:
    out=out+"ALL\t"+",".join(map(str,total))+"\t"+chain+"\tNA\t"+"\t".join(map(str,[0]*len(mutations_all)))+"\n"
  Write_out(out, cluster_annotation_file)
  out, ind = '',0
  return()

def Allele_specific_plots (annot_file, seq_file, cluster_file, allele_specific_file1, allele_specific_file2):
  seqs,freqs,alias = Get_sequences(seq_file)
  CDR3s,mutations = Get_annotations(annot_file)
  allele_group = Tree()
  reclassify, classes,classification = Get_reclassification(alias)
  for id in freqs:
    if(id in mutations):
      f = freqs[id]
      non_zero = [i for i in range(len(f)) if f[i]!=0]
      for i in non_zero:
        allele_group[classification[i].split("*")[0]][classification[i]][id][f[i]].value = 1
  out="#isotype\tallele\tnumber_seqs\tv\tj\n"
  out1 = "#isotype\tallele\tnumber_seqs\tmutations\n"
  for isotype in allele_group:
    if(len(allele_group[isotype])>1):
      for allele in allele_group[isotype]:
        gene_allele, mutation_allele = {},{}
        total = 0 
        for id in allele_group[isotype][allele]:
          for f in allele_group[isotype][allele][id]:
            total = total+f
          vj, mut = mutations[id][3], mutations[id][0]+mutations[id][1]
          if(vj in gene_allele):gene_allele[vj] = gene_allele[vj]+f
          else:gene_allele[vj]=f
          if(mut in mutation_allele):mutation_allele[mut] = mutation_allele[mut]+f
          else:mutation_allele[mut]=f
        for vj in gene_allele:
          out=out+isotype+"\t"+allele+"\t"+str(gene_allele[vj])+"\t"+vj+"\n"
        for mut in mutation_allele:
          out1 = out1+isotype+"\t"+allele+"\t"+str(mutation_allele[mut])+"\t"+str(mut)+"\n"
  fh=open(allele_specific_file1, "w")
  fh.write(out)
  fh.close()
  fh=open(allele_specific_file2, "w")
  fh.write(out1)
  fh.close()
  return()

def Write_out(out, file):
  fh = open (file,"a")
  fh.write(out)
  fh.close()
  return()

###########################
id = sys.argv[1]
group = sys.argv[2]
input_dir = sys.argv[3]+"ORIENTATED_SEQUENCES/"
output_dir = input_dir
species = sys.argv[4]
reverse_primer_group = sys.argv[5]
print("reverse_primer_group ",reverse_primer_group)
command_source = ["1"]
###########################
id1 = id
if(id.count("unmutated_")!=0):
  id= id.replace("unmutated_","")
###########################
seq_file = input_dir+"NETWORKS/Fully_reduced_"+id+".fasta"
annot_file = input_dir+"ANNOTATIONS/IMGT_ISOTYPER/IMGT_"+id+"_1_Summary.txt"
annot_file2 = input_dir+"ANNOTATIONS/IMGT_ISOTYPER/IMGT_"+id+"_8_V-REGION-nt-mutation-statistics.txt"
annot_file3 = input_dir+"ANNOTATIONS/IMGT_ISOTYPER/IMGT_"+id+"_3_Nt-sequences.txt"
annot_file4 = input_dir+"ANNOTATIONS/IMGT_ISOTYPER/IMGT_"+id+"_5_AA-sequences.txt"
annot_file7 = input_dir+"ANNOTATIONS/IMGT_ISOTYPER/IMGT_"+id+"_7_V-REGION-mutation-and-AA-change-table.txt"
annot_file_internal = input_dir+"ANNOTATIONS/TMP/Annotation_"+id+".txt"
#cluster_file = input_dir+"NETWORKS/Cluster_CDR3_defined_identities_"+id+".txt"
cluster_file = input_dir+"NETWORKS/Cluster_identities_"+id+".txt"
if(id.count("alpha")+id.count("beta")>0):
  cluster_file = input_dir+"NETWORKS/Cluster_identities_"+id+".txt"
  seq_file = input_dir+"NETWORKS/Fully_reduced_"+id.replace("PUBLIC_REMOVED_","").replace("PUBLIC_ONLY_","")+".fasta"
  #cluster_file = input_dir+"NETWORKS/Cluster_identities_PUBLIC_REMOVED_"+id+".txt"
edge_file = input_dir+"NETWORKS/Edges_"+id+".txt"
tmp_file1 = output_dir+"ISOTYPER/Tmp_cluster_"+id
non_ORF_filtered_reads = input_dir+"TMP/Filtered_J_"+id+".fasta"
merged_cluster_file  = output_dir+"ISOTYPER/Merged_cluster_identities_"+id+".txt"
merged_cluster_nearest_neighbours = output_dir+"ISOTYPER/Merged_nearest_neighbours_"+id+".txt"
transition_classification_file_directed = output_dir+"ISOTYPER/Transition_classification_directed_"+id+".txt"
mutational_freq_file = output_dir+"ISOTYPER/Raw_mutational_frequencies_undirected_"+id+".txt"
mutational_gene_freq_file = output_dir+"ISOTYPER/Class_mutation_gene_frequencies_"+id+".txt"
GL_distance_class_distribution_file = output_dir+"ISOTYPER/GL_distance_class_distribution_"+id+".txt"
class_overlap_file= output_dir+"ISOTYPER/Isotype_overlap_probabilities_"+id+".txt"
class_overlap_subtype_file= output_dir+"ISOTYPER/Isotype_overlap_subtype_probabilities_"+id+".txt"
children_classification_file_directed = output_dir+"ISOTYPER/End_of_mutational_chain_class_distribution_"+id+".txt"
cluster_annotation_file = output_dir+"ISOTYPER/Cluster_annotation_file_"+id+".txt"
overall_cluster_large_file = output_dir+"ISOTYPER/Isotyper_cluster_large_file_"+id+".txt"
allele_specific_file1 = output_dir+"ISOTYPER/Allele_specific_analysis_genes_"+id+".txt"
allele_specific_file2 = output_dir+"ISOTYPER/Allele_specific_analysis_mutations_"+id+".txt"
connectivity_probability_file = output_dir+"ISOTYPER/Connectivity_probabilities_"+id+".txt"
parental_classification_file_directed = output_dir+"ISOTYPER/Parental_probabilities_"+id+".txt"
mean_node_size_mutated_unmutated = output_dir+"ISOTYPER/Mean_node_size_mutated_unmutated_"+id+".txt"
cluster_size_connection_sharing_probability_file = output_dir+"ISOTYPER/Cluster_size_connection_sharing_probability_"+id+".txt"
IgD_IgM_ratio_file = output_dir+"ISOTYPER/IgD_IgM_ratio_"+id+".txt"
isotype_specific_gene_frequency_file = output_dir+"ISOTYPER/Isotype_specific_gene_frequencies_"+id+".txt"
separated_sequence_files = output_dir+"ISOTYPER/Isotype_fully_reduced_"+id+"_"
cluster_mutation_sharing_probability = output_dir+"ISOTYPER/Cluster_mutation_sharing_probability_"+id+".txt"
Mutational_region_isotype_file = output_dir+"ISOTYPER/Mutational_region_isotypes_"+id+".txt"
subtybe_VJ_gene_usage_file = output_dir+"ISOTYPER/Subtybe_VJ_gene_usage_freqs_"+id+".txt"
maximum_cluster_file = output_dir+"ISOTYPER/MAXIMUM_CLUSTERS/Maximum_cluster_"+id+".fasta"
maximum_cluster_file_filtered = output_dir+"ISOTYPER/MAXIMUM_CLUSTERS/Maximum_filtered_cluster_"+id+".fasta"
distance_from_central_file = output_dir+"ISOTYPER/MAXIMUM_CLUSTERS/Distance_from_central_mutations_"+id+".txt"
pre_post_class_switching_file = output_dir+"ISOTYPER/Pre_post_class_switching_"+id+".txt"
developmental_classification_file = output_dir+"ISOTYPER/Developmental_classifiation_"+id+".txt"
per_cluster_developmental_classification_file = output_dir+"ISOTYPER/Developmental_per_cluster_classifiation_"+id+".txt"
developmental_classification_mutational_file =  output_dir+"ISOTYPER/Mutational_development_classifiation_"+id+".txt"
developmental_classification_cluster_size_file =  output_dir+"ISOTYPER/Cluster_size_development_classifiation_"+id+".txt"
developmental_classification_vj_usage_file = output_dir+"ISOTYPER/VJ_gene_usage_freqs_development_classifiation_"+id+".txt"
developmental_classification_CDR3_lengths_charges = output_dir+"ISOTYPER/CDR3_length_charge_development_classifiation_"+id+".txt"
per_cluster_developmental_classification_mutational_file =  output_dir+"ISOTYPER/per_cluster_Mutational_development_classifiation_"+id+".txt"
per_cluster_developmental_classification_cluster_size_file =  output_dir+"ISOTYPER/per_cluster_Cluster_size_development_classifiation_"+id+".txt"
per_cluster_developmental_classification_vj_usage_file = output_dir+"ISOTYPER/per_cluster_VJ_gene_usage_freqs_development_classifiation_"+id+".txt"
per_cluster_developmental_classification_CDR3_lengths_charges = output_dir+"ISOTYPER/per_cluster_CDR3_length_charge_development_classifiation_"+id+".txt"
per_cluster_developmental_classification_network_parameters = output_dir+"ISOTYPER/per_cluster_network_parameters_"+id+".txt"
per_sequence_classification_network_parameters = output_dir+"ISOTYPER/per_sequence_network_parameters_"+id+".txt"
example_cluster_annot_file = output_dir+"ISOTYPER/EXAMPLE_CLUSTERS/Example_clusters_"+id+".txt"
mutation_correlation_file = output_dir+"ISOTYPER/Mutation_correlation_file_"+id+".txt"
syn_non_syn_mutations = output_dir+"ISOTYPER/Mutation_selection_"+id+".txt"
syn_non_syn_mutations_summary = output_dir+"ISOTYPER/Mutation_summmary_selection_"+id+".txt"
indel_file = output_dir+"ISOTYPER/Indel_summary_"+id+".txt"
clonality_mutation_correlation = output_dir+"ISOTYPER/Clonality_above_n_mutations_"+id+".txt"
non_functional_rearrangments_count_file = output_dir+"ISOTYPER/Non_functional_counts_"+id+".txt"
per_cluster_per_isotype_degree_file = output_dir+"ISOTYPER/Per_cluster_per_isotype_degree_"+id+".txt"
CDR3_region_composition = output_dir+"ISOTYPER/CDR3_composition_"+id+".txt"
CDR3_5prime_charge = output_dir+"ISOTYPER/CDR3_5prime_charge_"+id+".txt"
CDR3_5prime_charge_long_short = output_dir+"ISOTYPER/CDR3_5prime_long_short_charge_"+id+".txt"
public_repertoire = output_dir+"ISOTYPER/Public_repertoire_"+id+".txt"
seq_unmutated_file = input_dir+"NETWORKS/Fully_reduced_unmutated_"+id+".fasta"
clusters_unmutated_file = input_dir+"NETWORKS/Cluster_identities_unmutated_"+id+".txt"
secondary_rearrangement_file = output_dir+"ISOTYPER/Secondary_rearrangements_"+id+".txt"
raw_secondary_rearrangement_file = output_dir+"ISOTYPER/Raw_secondary_rearrangements_"+id+".txt"
secondary_rearrangement_cdr3_lengths = output_dir+"ISOTYPER/CDR3_length_secondary_rearrangements_"+id+".txt"
clonality_file = output_dir+"ISOTYPER/Clonality_per_sample_"+id+".txt"
mapped_secondary_rearrangements = output_dir+"ISOTYPER/Mapped_secondary_rearrangements_"+id+".txt"
v_replacement_transition_counts =  output_dir+"ISOTYPER/V_replacement_transition_counts_"+id+".txt"
CDR3_length_file = output_dir+"ISOTYPER/Overall_CDR3_length_"+id+".txt"
CDR3_length_file_grouped = output_dir+"ISOTYPER/Overall_CDR3_grouped_length_"+id+".txt"
CDR3_charge_file= output_dir+"ISOTYPER/CDR3_charge_"+id+".txt"
CDR23_charge_file = output_dir+"ISOTYPER/CDR23_charge_"+id+".txt"
CDR3_length_vjgenes = output_dir+"ISOTYPER/CDR3_VJ_genes_"+id+".txt"
captured_auto_antibodies_file = output_dir+"ISOTYPER/Captured_auto_antibodies_file_"+id+".txt"
mean_IgD_IgM_ratio = output_dir+"ISOTYPER/IGD_M_ratios_"+id+".txt"
shared_isotype_counts = output_dir+"ISOTYPER/Shared_isotype_counts_"+id+".txt"
per_cluster_inter_isotype_sharing = output_dir+"ISOTYPER/Per_cluster_inter_isotype_sharing_"+id+".txt"
mutations_per_expanded_cluster_file = output_dir+"ISOTYPER/Mutations_per_expanded_cluster_"+id+".txt"
unique_cluster_isotypes_file = output_dir+"ISOTYPER/Unique_cluster_isotypes_"+id+".txt"
unexpanded_cluster_file = output_dir+"ISOTYPER/Unexpanded_cluster_isotypes_"+id+".txt"
expanded_cluster_annotation = output_dir+"ISOTYPER/Expanded_cluster_isotypes_mutations_"+id+".txt"
mutational_positions_per_gene_file = output_dir+"ISOTYPER/Mutational_positions_per_gene_"+id+".txt"
mutational_positions_per_gene_file_summary = output_dir+"ISOTYPER/Mutational_positions_summmary_per_gene_"+id+".txt"
mutational_distance_distribution_per_cluster = output_dir+"ISOTYPER/Distribution_of_mutational_distance_within_clusters_"+id+".txt"
V_gene_isotype_frequency_file = output_dir+"ISOTYPER/V_gene_isotype_frequency_"+id+".txt"
V_gene_isotype_frequency_file_expanded = output_dir+"ISOTYPER/V_expanded_clusters_gene_isotype_frequency_"+id+".txt"
V_gene_isotype_frequency_file_grouped = output_dir+"ISOTYPER/V_gene_grouped_isotype_frequency_"+id+".txt"
V4_34_quantification_file = output_dir+"ISOTYPER/V4_34_quantification_"+id+".txt"
n_unmutated_file = output_dir+"ISOTYPER/Unmutated_sequences_"+id+".txt"
CDR3_distribution = output_dir+"ISOTYPER/CDR3_length_distribution_"+id+".txt"
nearest_neighbour_file = output_dir+"ISOTYPER/Nearest_neighbour_"+id+".txt"
CDR3_length_v_genes_file = output_dir+"ISOTYPER/CDR3_length_v_genes_"+id+".txt"
per_V_gene_cluster_file = output_dir+"ISOTYPER/per_V_gene_cluster_file_"+id+".txt"
per_V_gene_cluster_summary_file= output_dir+"ISOTYPER/per_V_gene_summary_cluster_file_"+id+".txt"
public_CDR3_file =  output_dir+"ISOTYPER/Public_CDR3_file_"+id+".txt"
########################## subsampled
subsample_file = output_dir+"ISOTYPER/Subsample_"+id+".txt"
clonality_file_SUBSAMPLED = output_dir+"ISOTYPER/Clonality_per_sample_SUBSAMPLED_"+id+".txt"
isotype_usages_SUBSAMPLED = output_dir+"ISOTYPER/Isotype_usages_SUBSAMPLED_"+id+".txt"
per_cluster_developmental_classification_network_parameters_subsampled = output_dir+"ISOTYPER/per_cluster_network_parameters_SUBSAMPLED"+id+".txt"
isotype_sharing_SUBSAMPLED =  output_dir+"ISOTYPER/Sharing_isotype_SUBSAMPLED_"+id+".txt"
isotype_sharing_SUBSAMPLED_summary=  output_dir+"ISOTYPER/Summary_Sharing_isotype_SUBSAMPLED_"+id+".txt"
clusters_isotype_overlapping_SUBSAMPLED = output_dir+"ISOTYPER/Overlapping_isotype_SUBSAMPLED_"+id+".txt"
cluster_SUBSAMPLE_information_file = output_dir+"ISOTYPER/Cluster_SUBSAMPLE_information_"+id+".txt"
cluster_isotype_expansion_SUBSAMPLE_file = output_dir+"ISOTYPER/Expansion_cluster_isotype_SUBSAMPLE_"+id+".txt"
secondary_rearrangement_file_SAMPLED = output_dir+"ISOTYPER/Secondary_rearrangement_file_SAMPLED_"+id+".txt"
cluster_properties_per_isotype_SUBSAMPLED = output_dir+"ISOTYPER/Cluster_properties_per_isotype_SUBSAMPLED_"+id+".txt"
cluster_size_file = output_dir+"ISOTYPER/Cluster_size_SUBSAMPLED_"+id+".txt"
unique_CDR3_regions_per_isotype_group_file = output_dir+"ISOTYPER/Unique_CDR3_regions_per_isotype_group_"+id+".txt"
########################## subsampled clones
subsample_file_clones =output_dir+"ISOTYPER/CLONES_Subsample_"+id+".txt"
clonality_file_SUBSAMPLED_clone = output_dir+"ISOTYPER/CLONES_Clonality_per_sample_SUBSAMPLED_"+id+".txt"
isotype_usages_SUBSAMPLED_clone=output_dir+"ISOTYPER/CLONES_Isotype_usages_SUBSAMPLED_"+id+".txt"
cluster_SUBSAMPLE_information_file_clone=output_dir+"ISOTYPER/CLONES_Cluster_SUBSAMPLE_information_"+id+".txt"
cluster_isotype_expansion_SUBSAMPLE_file_clone=output_dir+"ISOTYPER/CLONES_Expansion_cluster_isotype_SUBSAMPLE_"+id+".txt"
per_cluster_developmental_classification_network_parameters_subsampled_CLONE = output_dir+"ISOTYPER/CLONES_per_cluster_network_parameters_SUBSAMPLED"+id+".txt"
###########################
sample = id
########################### TO RUN ON ANY NEW SAMPLE
type ="PBMC"
#type = "ISOTYPE_SWITCHING_PAPER"
#type = "other"
#type = "SORT_FCR2G"
#type = "SORT_SW.MEM"
if(id.count("CD19.SW.MEM")!=0):
  type = "SORT_SW.MEM"
elif(id.count("SORT")!=0):
  type = "SORT_FCR2G"
elif(id.count("MUR")!=0):
  type = "MURINE"
#type = "ISOTYPE_SWITCHING_PAPER"
if(id.count("alpha")+id.count("beta")>0):
  type = "TCR_CARLOS_CALDAS"


if(type =="PBMC"):
  subsample_depth_file = "/lustre/scratch118/infgen/team146/rbr1/MISEQ/AUTOIMMUNITY_AAV_SLE_HEALTHY1/ORIENTATED_SEQUENCES/ISOTYPER/Sampling_depth_per_isotype_PBMCs.txt"
  subsample_depth = 2000
elif(type =="SORT_FCR2G"):
  subsample_depth_file = "/lustre/scratch118/infgen/team146/rbr1/MISEQ/AUTOIMMUNITY_AAV_SLE_HEALTHY1/ORIENTATED_SEQUENCES/ISOTYPER/Sampling_depth_per_isotype_PBMCs.txt"
  subsample_depth = 2000
elif(type=="other"): 
  subsample_depth_file = "/lustre/scratch118/infgen/team146/rbr1/MISEQ/AUTOIMMUNITY_AAV_SLE_HEALTHY1/ORIENTATED_SEQUENCES/ISOTYPER/Sampling_depth_per_isotype_PBMCs.txt"
  subsample_depth = 2000
elif(type =="SORT_SW.MEM"):
  subsample_depth_file = "/lustre/scratch118/infgen/team146/rbr1/MISEQ/AUTOIMMUNITY_AAV_SLE_HEALTHY1/ORIENTATED_SEQUENCES/ISOTYPER/Sampling_depth_per_isotype_CD19.SW.MEM.txt"
  subsample_depth = 100
elif(type=="ISOTYPE_SWITCHING_PAPER"):
  subsample_depth_file = "/lustre/scratch118/infgen/team146/rbr1/MISEQ/AUTOIMMUNITY_AAV_SLE_HEALTHY1/ORIENTATED_SEQUENCES/ISOTYPER/Sampling_depth_per_isotype_PBMCs_class_switching.txt"
  subsample_depth = 2000
elif(type=="TCR_CARLOS_CALDAS"):
  subsample_depth_file = "/lustre/scratch118/infgen/team146/rbr1/MISEQ/TCR_CARLOS_CALDOS_BREAST_CANCER/ORIENTATED_SEQUENCES/ISOTYPER/Sampling_depth_per_isotype_TCRs.txt"
  subsample_depth = 500
  if(id.count("alpha")!=0):# and id.count("PUBLIC_ONLY")!=0):
    subsample_depth = 500
elif(type =="MURINE"):
   subsample_depth_file= "/lustre/scratch118/infgen/team146/rbr1/MISEQ/AUTOIMMUNITY_AAV_SLE_HEALTHY1/ORIENTATED_SEQUENCES/ISOTYPER/Sampling_depth_per_isotype_MURINE1.txt"
   subsample_depth = 500


command="STANDARD"
#command="PRINT FILES"
#command = "NONE"

if(command == "PRINT FILES"):
  a = [developmental_classification_mutational_file,
      developmental_classification_cluster_size_file,
      per_cluster_developmental_classification_mutational_file,
      per_cluster_developmental_classification_cluster_size_file, 
      per_cluster_developmental_classification_network_parameters,
      class_overlap_subtype_file,
      syn_non_syn_mutations_summary,
      clonality_file,
      shared_isotype_counts,
      per_cluster_inter_isotype_sharing,
      expanded_cluster_annotation,
      mutational_freq_file,
      developmental_classification_CDR3_lengths_charges]
  a = [developmental_classification_cluster_size_file,
      per_cluster_developmental_classification_mutational_file,
      per_cluster_developmental_classification_network_parameters,
      syn_non_syn_mutations_summary,
      mutational_freq_file,
      CDR3_length_file,
      clonality_file_SUBSAMPLED,
      isotype_usages_SUBSAMPLED,
      per_cluster_developmental_classification_network_parameters_subsampled,
      isotype_sharing_SUBSAMPLED,clusters_isotype_overlapping_SUBSAMPLED]
  for f in a:
  #  print "cat "+f.replace(id,"REP*")+"> "+f.replace(id,"REP").replace("ISOTYPER/","ISOTYPER/All_")
    print("cat "+f.replace(id,"*")+" > "+f.replace(id,"").replace("ISOTYPER/","ISOTYPER/All_").replace("_.txt",".txt"))

#command = 11
if(command == "STANDARD"):
  print("subsample depth: ",subsample_depth)
  #### TOP
  #Merge_clusters_on_CDR3_region(seq_file, annot_file, tmp_file1, merged_cluster_file, merged_cluster_nearest_neighbours,cluster_file,edge_file)
  #Classify_sequences_into_developmental_stages(sample, annot_file, cluster_file, developmental_classification_file,reverse_primer_group)
  #Classify_sequences_into_developmental_stages_per_cluster(sample, annot_file, cluster_file,per_cluster_developmental_classification_file,reverse_primer_group)

  #### Run on whole data
  #Get_mutational_frequencies(seq_file, annot_file_internal,mutational_freq_file,cluster_file,sample,annot_file,reverse_primer_group,n_unmutated_file)
  #Get_CDR3_lengths(seq_file,annot_file,CDR3_length_file,sample,annot_file4,CDR3_charge_file,CDR23_charge_file,reverse_primer_group,CDR3_distribution)
  ###Get_CDR3_region_composition(sample,seq_file,annot_file, per_cluster_developmental_classification_file,CDR3_region_composition,CDR3_5prime_charge,CDR3_length_file,annot_file4,CDR3_charge_file,CDR23_charge_file,mean_IgD_IgM_ratio,CDR3_5prime_charge_long_short, CDR3_length_vjgenes)
  #Get_secondary_rearrangement_potential(sample, annot_file,annot_file3, secondary_rearrangement_file,seq_file,raw_secondary_rearrangement_file,secondary_rearrangement_cdr3_lengths,mapped_secondary_rearrangements,tmp_file1,v_replacement_transition_counts,secondary_rearrangement_file_SAMPLED,cluster_file,reverse_primer_group)
  #Get_total_network_summary(sample, annot_file, isotype_usages_SUBSAMPLED_clone,reverse_primer_group)

  #Get_V_gene_isotype_frequency_grouped(sample, annot_file, V_gene_isotype_frequency_file_grouped,per_cluster_developmental_classification_file,CDR3_length_file_grouped,per_V_gene_cluster_file,per_V_gene_cluster_summary_file)
  #Get_network_parameters_per_classificiation_subsample_clones(sample,per_cluster_developmental_classification_file, per_cluster_developmental_classification_network_parameters_subsampled_CLONE,reverse_primer_group,cluster_file,subsample_depth_file,developmental_classification_file)
  #Subsample_repertoire(subsample_file, subsample_depth, seq_file,sample,cluster_file)
  #Get_V_gene_isotype_frequency(sample, annot_file, V_gene_isotype_frequency_file,per_cluster_developmental_classification_file,V_gene_isotype_frequency_file_expanded)









  #### optional
  #Get_mutational_frequencies_per_classification_per_cluster(sample, annot_file_internal,per_cluster_developmental_classification_file,per_cluster_developmental_classification_mutational_file, per_cluster_developmental_classification_cluster_size_file,per_cluster_developmental_classification_vj_usage_file,reverse_primer_group)
  Isotype_per_cluster_sharing_probabilities(seq_file,per_cluster_developmental_classification_file, per_cluster_inter_isotype_sharing,sample,reverse_primer_group)
  #Get_mutations_per_expanded_cluster(seq_file, annot_file2, sample,mutations_per_expanded_cluster_file,cluster_file,reverse_primer_group)
  #Mutation_selection_per_chain(seq_file, annot_file2, sample, syn_non_syn_mutations,syn_non_syn_mutations_summary,reverse_primer_group)
  #Indel_analysis_per_chain(seq_file, sample, indel_file,annot_file,reverse_primer_group)

  #Get_network_parameters_per_sequence_classification(sample, cluster_file, developmental_classification_file, per_sequence_classification_network_parameters,reverse_primer_group)

  #### Subsampled for clonality
  #Get_network_parameters_per_classificiation_subsample(sample,per_cluster_developmental_classification_file, per_cluster_developmental_classification_network_parameters_subsampled,reverse_primer_group,cluster_file,subsample_depth_file,developmental_classification_file)
  #Isotype_sharing_subsampled(sample,reverse_primer_group, subsample_file,isotype_sharing_SUBSAMPLED,cluster_file,annot_file,isotype_sharing_SUBSAMPLED_summary)

  #Get_cluster_properties_per_isotype(cluster_properties_per_isotype_SUBSAMPLED, sample, per_cluster_developmental_classification_file,subsample_file,cluster_size_file)
  #Proportion_clusters_isotype_overlapping(sample, reverse_primer_group, subsample_file, clusters_isotype_overlapping_SUBSAMPLED)
  #Get_unique_CDR3_regions_per_isotype_group(subsample_file,annot_file, per_cluster_developmental_classification_file, unique_CDR3_regions_per_isotype_group_file,annot_file3,sample,CDR3_length_v_genes_file)
  #Get_public_CDR3_percentages(subsample_file,annot_file, per_cluster_developmental_classification_file, unique_CDR3_regions_per_isotype_group_file,annot_file3,sample,public_CDR3_file)

  #### Unsampled
  #Get_network_parameters_per_classificiation(sample, cluster_file, per_cluster_developmental_classification_file, per_cluster_developmental_classification_network_parameters,reverse_primer_group)
  #Get_class_overlap_subclass(seq_file,class_overlap_subtype_file,sample,developmental_classification_file,cluster_file,reverse_primer_group)
  #Clonality_per_sample(seq_file,sample,clonality_file,cluster_file,reverse_primer_group)
  #Get_proportion_read_mutiple_subtypes(seq_file,per_cluster_developmental_classification_file, shared_isotype_counts,sample,reverse_primer_group)
  #Get_CDR3_lengths_charges_per_classification(sample,seq_file,annot_file, per_cluster_developmental_classification_file,developmental_classification_CDR3_lengths_charges,annot_file_internal,reverse_primer_group)
 # Get_nearest_neighbour_distributions(sample, annot_file, edge_file, nearest_neighbour_file)

########################### DEVELOPEMENT
#Get_probability_connected_nodes(cluster_file, connectivity_probability_file,seq_file,sample)
#Get_class_overlap(seq_file, class_overlap_file,sample)
#Get_mutational_frequencies(seq_file, annot_file_internal,mutational_freq_file,cluster_file,sample,annot_file)
#Get_gene_mutation_frequencies(seq_file, annot_file_internal,mutational_gene_freq_file,cluster_file,sample,annot_file)
#Get_gene_distribution_isotype_specific(seq_file, annot_file_internal,isotype_specific_gene_frequency_file,sample,annot_file)
#IgD_IgM_ratio(annot_file_internal, seq_file,IgD_IgM_ratio_file,cluster_file,sample)
#Clone_size_of_mutated_and_unmutated(annot_file_internal, seq_file,mean_node_size_mutated_unmutated,sample)
#Get_cluster_mutation_sharing_stats(cluster_file, annot_file_internal, seq_file, cluster_mutation_sharing_probability,sample)
#Get_VJ_gene_usage_subtype_specific(annot_file,annot_file_internal, seq_file,subtybe_VJ_gene_usage_file,sample)
#Classify_transitions(seq_file,edge_file, transition_classification_file_directed,cluster_file, children_classification_file_directed,parental_classification_file_directed,annot_file_internal,sample)
#Get_maximum_cluster_information(cluster_file, seq_file, sample, edge_file, annot_file_internal, maximum_cluster_file ,distance_from_central_file,maximum_cluster_file_filtered)
#Get_pre_post_class_switching(cluster_file, sample, annot_file_internal, pre_post_class_switching_file, seq_file)
#Mutation_selection_per_chain(seq_file, annot_file2, sample, syn_non_syn_mutations,syn_non_syn_mutations_summary)
#Get_clonality_above_n_mutations(seq_file, annot_file, clonality_mutation_correlation,cluster_file,sample)
#Get_proportion_non_functional_reads(seq_file, non_functional_rearrangments_count_file, non_ORF_filtered_reads,sample)
#Get_per_cluster_per_isotype_degree(seq_file, edge_file, cluster_file, per_cluster_per_isotype_degree_file,sample,annot_file)
#if(id1.count("unmutated_")==0):
#Get_CDR3_region_composition(sample,seq_file,annot_file, per_cluster_developmental_classification_file,CDR3_region_composition,CDR3_5prime_charge,CDR3_length_file,annot_file4,CDR3_charge_file,CDR23_charge_file,mean_IgD_IgM_ratio,CDR3_5prime_charge_long_short, CDR3_length_vjgenes)
#Clonality_per_sample(seq_file,sample,clonality_file,cluster_file)
#Get_percentage_public_repertoire(sample,seq_unmutated_file,annot_file,public_repertoire)
#Get_secondary_rearrangement_potential(sample, annot_file,annot_file3, secondary_rearrangement_file,seq_file,raw_secondary_rearrangement_file,secondary_rearrangement_cdr3_lengths,mapped_secondary_rearrangements,tmp_file1,v_replacement_transition_counts,secondary_rearrangement_file_SAMPLED,cluster_file,reverse_primer_group)








#Screen_auto_Ab_sequences(annot_file4, sample, seq_file, captured_auto_antibodies_file)
#Get_proportion_read_mutiple_subtypes(seq_file,per_cluster_developmental_classification_file, shared_isotype_counts,sample)
#Isotype_per_cluster_sharing_probabilities(seq_file,per_cluster_developmental_classification_file, per_cluster_inter_isotype_sharing,sample)
#Get_mutations_per_expanded_cluster(seq_file, annot_file2, sample,mutations_per_expanded_cluster_file,cluster_file)
#Number_of_unique_clusters(seq_file,annot_file, unique_cluster_isotypes_file, sample)
#Get_unexpanded_cluster_proportions(seq_file, sample,unexpanded_cluster_file,cluster_file)
#Get_mutations_and_switching_of_expanded_clusters_specific_isotypes(seq_file, expanded_cluster_annotation, cluster_file, annot_file)
#Get_mutational_positions_per_gene(seq_file, developmental_classification_file, mutational_positions_per_gene_file,annot_file7,sample,mutational_positions_per_gene_file_summary)
#Get_example_sequence(seq_file, annot_file)
#Get_mutational_distance_distribution_per_cluster(seq_file, per_cluster_developmental_classification_file, mutational_distance_distribution_per_cluster,tmp_file1,sample)

########### Getting unmutated IGHD/M sequences for some subsets
#Get_unmutated_repertoire(seq_file, annot_file, seq_unmutated_file,clusters_unmutated_file,cluster_file)
#Get_mutational_frequencies(seq_unmutated_file, annot_file_internal,mutational_freq_file,cluster_file,sample,annot_file)
#Classify_sequences_into_developmental_stages_per_cluster(sample, annot_file, clusters_unmutated_file,per_cluster_developmental_classification_file)
#Get_CDR3_lengths_charges_per_classification(sample,seq_unmutated_file,annot_file, per_cluster_developmental_classification_file,developmental_classification_CDR3_lengths_charges,annot_file_internal)

#if(id1.count("unmutated_")!=0):
#  Get_secondary_rearrangement_potential(sample, annot_file,annot_file3, secondary_rearrangement_file,seq_unmutated_file,raw_secondary_rearrangement_file,secondary_rearrangement_cdr3_lengths,mapped_secondary_rearrangements,tmp_file1,v_replacement_transition_counts)
#  Clonality_per_sample(seq_unmutated_file,sample,clonality_file,cluster_file)
#  Get_VJ_gene_usage_subtype_specific(annot_file_internal, seq_unmutated_file,subtybe_VJ_gene_usage_file,sample)
#  Get_percentage_public_repertoire(sample,seq_unmutated_file,annot_file,public_repertoire) 
#  Get_gene_distribution_isotype_specific(seq_unmutated_file, annot_file_internal,isotype_specific_gene_frequency_file,sample,annot_file)




################# Verification tests for annotations
#Compare_mutations_IMGT_internal(annot_file_internal,annot_file,mutation_correlation_file,id,seq_file) 

#########################################
#Identify_example_clusters(sample, annot_file_internal, cluster_file,per_cluster_developmental_classification_file,example_cluster_annot_file,seq_file)


#################Not used yet:
#Get_GL_distance_class_distributions(seq_file, annot_file,GL_distance_class_distribution_file)
#Annotating_clusters_for_analysis(merged_cluster_file, cluster_annotation_file, annot_file,seq_file)
#Make_trees_from_clusters(cluster_annotation_file,merged_cluster_file,seq_file,overall_cluster_large_file,id, output_dir)
#Allele_specific_plots (annot_file, seq_file, cluster_file, allele_specific_file1, allele_specific_file2)
#Cluster_size_connection_sharing_probability(annot_file, seq_file, merged_cluster_nearest_neighbours, cluster_size_connection_sharing_probability_file)
#Get_all_chain_sequences_separately(seq_file, separated_sequence_files)
#Mutational_region_isotypes(cluster_file, annot_file2, seq_file, Mutational_region_isotype_file)
#Get_autoreactive_V4_34_motif_frequency(annot_file4, V4_34_quantification_file, sample,seq_file,reverse_primer_group,annot_file)

## TOP
