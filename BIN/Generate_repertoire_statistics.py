#!/usr/bin/env python
import math
import sys
from collections import defaultdict
sys.path.append('BIN/')
import Align
from Align import *
import Functions
from Functions import *
import os
import re
import time
import commands
import numpy as np
from numpy import outer
from operator import itemgetter, attrgetter, add
import copy
import random
import networkx as nx

def fasta_iterator(fh):
  while True:
    line = fh.readline()
    if line.startswith('>'): break  
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

def Reduce(seq):
  s = seq[0]
  for i in range(1,len(seq)):
    if(seq[i]!=seq[i-1]):
      s = s+seq[i]
  return(s)

def Get_library(ref,end,max_len):
  fh = open (ref, "r")
  seqs=[]
  ids=[]
  all={}
  ident = 5
  length={}
  for header,sequence in fasta_iterator(fh):
    all[header]=1
    length[header]=sequence
    if(end=="end"):
      sequence = sequence[len(sequence)-max_len-ident:len(sequence)-ident]
    else:
      if(end=="start"):
        sequence = sequence[ident:max_len+ident]
      else:
        sequence = sequence[ident:len(sequence)-ident]
    se = Reduce(sequence)
    ids.append(header)
    seqs.append(se)
  fh.close()
  (same,removed)=Get_equivalent(seqs,ids)
  del seqs
  del ids
  del removed
  return(same,length)

def Get_sequence_length_distributions(seq_file,seq_length_ditribution):
  fh=open(seq_file,"r")
  lengths = {}
  for header,sequence in fasta_iterator(fh):
    l = len(sequence)
    f = sum(map(int, header.split("|")[0].split("__")[1].split("_")))
    if(l in lengths):lengths[l] = lengths[l]+f
    else:lengths[l]=f
  fh.close()
  out= "#length\tfrequency\n"
  mean,total = 0,0
  for l in lengths:
    out=out+str(l)+"\t"+str(lengths[l])+"\n"
    mean = mean+(l*lengths[l])
    total = total+lengths[l]
  fh=open(seq_length_ditribution,"w")
  fh.write(out)
  fh.close()
  print seq_length_ditribution,"\t",mean*1.0/total
  return()

def Get_equivalent(seqs,ids):
  same={}
  removed={}
  for i in range(0,len(seqs)):
    if(ids[i] not in same):
      same[ids[i].split("*")[0]]=ids[i].split("*")[0]
      l1 = len(seqs[i])
      for j in range(i,len(seqs)):
        if(i<j):
          l2 =len(seqs[j])
          (s1,s2,p)=Trim_library(seqs[i],seqs[j],l1,l2)
          if(p==1):
            if(s1==s2):
              same[ids[j].split("*")[0]]=ids[i].split("*")[0]
              removed[ids[j]]=1
  return(same,removed)

def Trim_library(s1,s2,l1,l2):
  p=0
  ind = 4
  sample=s1[ind:l1-ind]
  index=s2.find(sample)
  if (index!=-1):
    if(index>ind):
      s2=s2[index-ind:l2]
    else:
      s1=s1[ind-index:l1]
    min_len=min([len(s1),len(s2)])
    p=1
  return(s1,s2,p)

def Blast_match_D(out, refd, tmp_file, seqs,passed, file_D_genes,batch_number):
  fh=open(tmp_file, "w")
  fh.write(out)
  fh.close() 
  command1 = "/software/pubseq/bin/ncbi_blast/blastall -p blastn -a 1 -d "+refd+" -e 2 -i "+tmp_file+" -o "+tmp_file+"_dblast -b 1 -m 8 -W 4"
  os.system(command1)
  fh = open(tmp_file+"_dblast","r")
  out=''
  if(batch_number==0):
    out="#ID\tD_gene\tD_location\tD_sequence\tD_orientation\n"
  genes = {}
  for l in fh:
    l=l.strip().split()
    if(l[0] in genes):
      score_orig = genes[l[0]][1]
      if(score_orig<float(l[11])):
        genes[l[0]]=[l[1], float(l[11]), int(l[6]), int(l[7])]
    else:
      #            D gene, score, start, end
      genes[l[0]]=[l[1], float(l[11]), int(l[6]), int(l[7])]
  fh.close()
  for id in genes:
    a = genes[id]
    start, end = a[2]+seqs[id][1], a[3]+seqs[id][1]
    out=out+id+"\t"+a[0]+"\t"+str(start)+"\t"+str(end)+"\t"+seqs[id][0][a[2]:a[3]]+"\n"
    passed = passed+1
  Write_out(out, file_D_genes)
  out=''
  return(passed)

def Locate_CDR3_start_site_IGH(CDR3_prot,seq, shift,codon,cdr3_ends,cdr3_ends_near):
  score,passed,nn_start = [],0,0
  CDR3_prot_trim,nn ='-','-'
  for b in ["R","V","S","L","K","P"]: 
    if(CDR3_prot.count("CA"+b)==1 and passed ==0):
      ind = CDR3_prot.index("CA"+b)
      score.append([CDR3_prot[ind:len(CDR3_prot)], ind+shift])
      passed = 1
      break
  if(passed==0):
    for b in ["R","S"]:
      if(re.search(r"C."+b, CDR3_prot) and passed ==0):
        ind = CDR3_prot.index(re.search(r"C."+b, CDR3_prot).group(0))
        score.append([CDR3_prot[ind:len(CDR3_prot)], ind+shift])
        passed = 1
        break
  if(passed==0):
    if(CDR3_prot.count("CA")==1 and passed ==0):
      ind = CDR3_prot.index("CA")
      score.append([CDR3_prot[ind:len(CDR3_prot)], ind+shift])
      passed = 1
  if(passed==0):
    for b in ["V","L","K"]:
      if(re.search(r"C."+b, CDR3_prot) and passed ==0):
        ind = CDR3_prot.index(re.search(r"C."+b,CDR3_prot).group(0))
        score.append([CDR3_prot[ind:len(CDR3_prot)], ind+shift])
        passed = 1
        break
  if(passed==0):
    if(CDR3_prot.count("CV")==1 and passed ==0):
      ind = CDR3_prot.index("CV")
      score.append([CDR3_prot[ind:len(CDR3_prot)], ind+shift])
      passed = 1
  if(passed==0):
    for b in ["R","V","S","L","K","T","H","D","M","N"]:
      if(re.search(r".A"+b, CDR3_prot) and passed ==0):
        ind = CDR3_prot.index(re.search(r".A"+b,CDR3_prot).group(0))
        score.append([CDR3_prot[ind:len(CDR3_prot)], ind+shift])
        passed = 1
        break
  if(passed==0):
    for b in ["C"]:
      if(CDR3_prot.count(b)!=0 and passed ==0):
        ind = CDR3_prot.index(b)
        score.append([CDR3_prot[ind:len(CDR3_prot)], ind+shift])
        passed = 1
        break
  if(passed==0):
    for b1 in ["R","L"]:
      for b in ["R","V","S","L","K","T","H","D","M","N"]:
        if(re.search(r"."+b1+b, CDR3_prot) and passed ==0):
          ind = CDR3_prot.index(re.search(r"."+b1+b,CDR3_prot).group(0))
          score.append([CDR3_prot[ind:len(CDR3_prot)], ind+shift])
          passed = 1
          break
      if(passed==1):break
  if(passed==1):
    # get end of CDR3 region
    CDR3_prot = score[0][0]
    for trial in cdr3_ends:
      if(CDR3_prot.count(trial)!=0):
        ind = CDR3_prot.index(trial)+cdr3_ends[trial]
        CDR3_prot_trim = CDR3_prot[0:ind]
        passed=2
        break
      if(passed==2):break
    if(passed!=2):
      for trial in cdr3_ends_near:
        options = ''
        if(re.search(r""+trial, CDR3_prot)):
          ind = CDR3_prot.index(re.search(r""+trial, CDR3_prot).group(0))-cdr3_ends_near[trial]
          CDR3_prot_trim = CDR3_prot[0:ind]
          if(len(CDR3_prot_trim)>0):
            if(len(options)<CDR3_prot_trim):options = CDR3_prot_trim
        if(len(options)>0):
          CDR3_prot_trim = options
          passed=2
          break
    if(passed!=2):
      if(CDR3_prot.count("W")!=0):
        ind = CDR3_prot.index("W")+1
        CDR3_prot_trim = CDR3_prot[0:ind]
        passed=2
    if(passed!=2):CDR3_prot_trim = CDR3_prot
    for i in [0,1,2]:
      prot = Translate(seq[i:len(seq)],codon)
      if(prot.count(score[0][0])!=0):
        ind = i+prot.index(score[0][0])*3
        end = min([ind+((len(CDR3_prot_trim)+1)*3),len(seq)])
        nn = seq[ind:end]
        nn_start = ind
        break
  else:
    print CDR3_prot,"NOT FOUND" 
  return(CDR3_prot_trim, passed,nn_start,nn)

def Blast_match(out, refv,refj, tmp_file, out1,indent,seqs,lengthv,prots,codon,gene,cdr3_ends,cdr3_ends_near):
  fh=open(tmp_file, "w")
  fh.write(out)
  fh.close()
  blast = "/software/pubseq/bin/ncbi_blast/blastall"
  blast = "/usr/bin/blastall"
  command1 = blast+" -p blastn -a 10 -d "+refv+" -e 1e-5 -i "+tmp_file+" -o "+tmp_file+"_vblast -b 1 -m 8 -W 10"
  command2 = blast+" -p blastn -a 10 -d "+refj+" -e 10 -i "+tmp_file+" -o "+tmp_file+"_jblast -b 1 -m 8 -W 4"
  os.system(command1)
  fh=open(tmp_file, "w")
  fh.write(out1)
  fh.close()
  os.system(command2)
  fh = open(tmp_file+"_jblast","r")
  passesv,passesj = {},{}
  for l in fh:
    l=l.strip().split()
    if(int(l[6])<int(l[7]) and int(l[8])<int(l[9])):
      if(l[0] in seqs):
        addition = len(seqs[l[0]])-80 
        passesj[l[0]]= (l[1].upper(), int(l[6])+addition, int(l[7])+addition, int(l[8]), int(l[9]), int(l[3]), int(l[4]), int(l[5]))
      else:
        print l[0]
  fh.close()
  fh = open(tmp_file+"_vblast","r")
  for l in fh:
    l=l.strip().split()
    if(len(l)>=10):
      if(int(l[6])<int(l[7]) and int(l[8])<int(l[9])):
        v_gene = l[1]
        if(v_gene=="IGHV4-B"):v_gene="IGHV4-b*01"
        mut = int(l[4])
        if(l[0] in passesv):
          additional_mutations,v_end_query = passesv[l[0]][11],passesv[l[0]][13]
          if(passesv[l[0]][13]<int(l[7])):
            v_end_query,v_end_ref = int(l[7]), int(l[7])+len(lengthv[l[1]])-int(l[9])
            nn_start = passesv[l[0]][15]
            CDR3_start = passesv[l[0]][9]
            additional_mutations1 = max([nn_start-v_end_query,0])
            if(additional_mutations1<passesv[l[0]][11]):
              additional_mutations = additional_mutations1
          mut,ins = int(l[4])+passesv[l[0]][6], int(l[5])+passesv[l[0]][7]
          passesv[l[0]] = (passesv[l[0]][0],passesv[l[0]][1],passesv[l[0]][2],passesv[l[0]][3],passesv[l[0]][4],passesv[l[0]][5],mut,ins,passesv[l[0]][8],passesv[l[0]][9],passesv[l[0]][10],additional_mutations,passesv[l[0]][12],v_end_query, passesv[l[0]][14],v_end_query)
        else:
          mut = int(l[4])
          v_end_query,v_end_ref = int(l[7]), int(l[7])+len(lengthv[l[1]])-int(l[9])
          end_diff = v_end_ref-v_end_query
          V_end_potential = v_end_query
          if(l[0].split("__")[0] in prots):
            CDR3_prot_potential = prots[l[0].split("__")[0]]
            shift = (V_end_potential/3)-5
            CDR3_prot_potential = CDR3_prot_potential[shift:len(CDR3_prot_potential)]
            if(gene=="IGH"):
              CDR3_p, passed,nn_start,CDR3_nn = Locate_CDR3_start_site_IGH(CDR3_prot_potential, seqs[l[0]], shift,codon,cdr3_ends,cdr3_ends_near)
            if(gene.count("TCR")==1):
              #CDR3_p, passed,nn_start,CDR3_nn = Locate_CDR3_start_site_TCR(CDR3_prot_potential, seqs[l[0]], shift,codon)################################
              CDR3_start, CDR3, additional_mutations, CDR3_found,CDR3_nn,nn_start = v_end_query,"-",v_end_query,"NO","-",v_end_query
              passed = 0
            if(passed!=0):
              CDR3_start, CDR3, additional_mutations, CDR3_found = nn_start, CDR3_p,max([nn_start-v_end_query,0]),"YES"
              if(additional_mutations>0):
                #print "R:",lengthv[l[1]][int(l[8]):int(l[9])]
                #print "Q:",seqs[l[0]][int(l[6]):int(l[7])]
                ref_gene = lengthv[l[1]][int(l[9]):len(lengthv[l[1]])]
                query_seq = seqs[l[0]][int(l[7]):CDR3_start]
                minl = min([len(ref_gene), len(query_seq)])
                ref_gene, query_seq = ref_gene[0:minl], query_seq[0:minl]
                additional_mutations = len([i for i in range(len(ref_gene)) if ref_gene[i]!=query_seq[i]])
                additional_mutations = min([additional_mutations, 20])
            else:CDR3_start, CDR3, additional_mutations, CDR3_found,CDR3_nn,nn_start = v_end_query,"-",v_end_query,"NO","-",v_end_query
            passesv[l[0]]= (v_gene, int(l[6]), int(l[7]), int(l[8]), int(l[9]), int(l[3]), mut, int(l[5]), end_diff, CDR3_start, CDR3, additional_mutations, CDR3_found,v_end_query,CDR3_nn,nn_start)
          else: print l
    else:
      print l
  fh.close()
  return(passesv,passesj)

def Get_codons (codon_file):
  fh = open(codon_file,"r")
  codon = Tree()
  for l in fh:
    l=l.strip()
    l1 = l.split()
    l2 = list(l1[0])
    codon[l2[0]][l2[1]][l2[2]][l1[1]].value=1
  fh.close()
  return (codon)

def Translate(seq,codon):
  p_seq=""
  for cod in range(0, len(seq)/3-1):
    cod = cod*3
    for c in codon[seq[cod+0]][seq[cod+1]][seq[cod+2]]:
      p_seq=p_seq+str(c)
  return(p_seq)

def Get_annotation(passesv,passesj, seqs, out_file,lengthj,regions,failed,CDR3_end,prots, codon,gene):
  #head = "#ID\tIGHV\tFR1s\tFR1e\tCDR1s\tCDR1e\tFR2s\tFR2e\tCDR2s\tCDR2e\tFR3s\tFR3e\tIGHV_end\tIGHJ\tIGHJ_start\tCDR3s\tCDR3e\tCDR3_prot\tCDR3_nn\n"
  out=''
  for id in seqs:
    #o=id.split("__")[0] ################# remove
    o=id
    if(id in passesv):
      v_found = [passesv[id][0],passesv[id][1], passesv[id][2],passesv[id][3],passesv[id][4],passesv[id][5], passesv[id][6], passesv[id][7]]
      o=o+"\t"+v_found[0]
      mmv,mmj,insv,insj=v_found[6],-1,v_found[7],-1
      if(passesv[id][12]=="YES"):mmv = mmv+passesv[id][11]
      if(v_found[0].split("*")[0].upper() in regions): 
        reg = regions[v_found[0].split("*")[0].upper() ]
        for i in range(0,len(reg)):
          if(v_found[1]+reg[i]-v_found[3]>=0):o=o+"\t"+str(v_found[1]+reg[i]-v_found[3])
          else:
            if(i<len(reg)-1 and v_found[1]+reg[i+1]-v_found[3]>=0):o=o+"\t0"
            else:o=o+"\t-1"
      else: 
        for i in range(0,10):
          o=o+"\t-1" 
      o=o+"\t"+str(v_found[2])
      if(id in passesj):
        j_found = [passesj[id][0],passesj[id][1], passesj[id][2],passesj[id][3],passesj[id][4],passesj[id][5], passesj[id][6], passesj[id][7]]
        mmj,insj = j_found[6],j_found[7]
        o=o+"\t"+j_found[0]+"\t"+str(j_found[1])
      else:
        o=o+"\tUNKNOWN\t-1"
        mmj,insj = 0,0
      CDR3_p,CDR3_nn = passesv[id][10], passesv[id][14]
      o=o+"\t"+CDR3_p+"\t"+CDR3_nn
      o=o+"\t"+str(mmv)+"\t"+str(mmj)+"\t"+str(insv)+"\t"+str(insj)
      out=out+o+"\n"
    else:
      o=o+"\tV_not_found"
      failed = failed+1
  Write_out(out, out_file)
  out=''
  return(failed)

def Get_region_boundaries(species,loc):
  regions={}
  if(species=="HOMO_SAPIENS"):
    fh=open(loc+"LIBRARY/human.ndm.imgt2","r")
    J_ref = loc+"LIBRARY/human.ndm.imgt.CDR3.end"
  elif(species=="MUS_MUSCULUS"):
    fh=open(loc+"LIBRARY/mouse.ndm.imgt","r")
    J_ref = loc+"LIBRARY/mouse.ndm.imgt.CDR3.end"
  elif(species=="LLAMA_GLAMA"):
    fh=open(loc+"LIBRARY/llama.ndm.imgt","r")
    J_ref = loc+"LIBRARY/llama.ndm.imgt.CDR3.end"
  else:print "Reference kabat not found"
  for l in fh:
    l=l.strip().split()
    id = l[0]
    if(id.count("VH")==1):id.replace("VH","IGHV")
    id=id.split("*")[0].upper()
    regions[id]=(int(l[1]),int(l[2]),int(l[3]),int(l[4]),int(l[5]),int(l[6]),int(l[7]),int(l[8]),int(l[9]),int(l[10]))
  fh.close()
  fh=open(J_ref,"r")
  CDR3_end={}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id,CDR3_end_nn, CDR3_end_aa = l[0].split("*")[0], int(l[3]), int(l[2])
      CDR3_end[id]=[CDR3_end_nn, CDR3_end_aa]
  fh.close()
  return(regions, CDR3_end)

def Get_protein_sequences(protein_file):
  seqs={}
  fh=open(protein_file, "r")
  for header,seq in fasta_iterator(fh):
    seqs[header.split("__")[0]]=seq.upper()
  fh.close()
  return(seqs)

def Identify_D_genes(annot_file, file_D_genes, refd,tmp_file,seq_file):
  fh=open(file_D_genes,"w")
  fh.close()
  fh=open(annot_file,"r")
  out, ind,seqs,passed, total,batch,batch_number = '',0,{},0,0,1000,0
  D_approx_region = {}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(len(l)>19):
        if(l[15].count("CDR3_")==0):
          V_end, J_start = int(l[12]),int(l[14])
          CDR3_start = int(l[15])
          D_approx_region[l[0]]=[V_end, J_start,V_end]
  fh.close()
  fh=open(seq_file,"r")
  for header,seq in fasta_iterator(fh):
    if(header in D_approx_region):
      seq_d_poss = seq[D_approx_region[header][0]:D_approx_region[header][1]]
      if(len(seq_d_poss)>5):
        out = out+">"+header+"\n"+seq_d_poss+"\n"
        seqs[header]=[seq_d_poss, D_approx_region[header][2]]
        ind,total = ind+1,total+1
        if(ind>batch):
          passed = Blast_match_D(out, refd, tmp_file,seqs,passed, file_D_genes,batch_number)
          batch_number=batch_number+1
          del seqs, out
          seqs, out,ind={},'',0
          print "\r",batch_number, 
          #if(batch_number>30):break
  fh.close()
  print "\nTotal:",total,"Passed:", passed, "% Passed:", passed*100.0/total
  os.system("rm "+tmp_file+" "+tmp_file+"_dblast")
  return()

def Get_CDR3_ends(CDR3_end_file,gene):
  cdr3_ends,cdr3_ends_near = {},{}
  fh=open(CDR3_end_file,"r")
  for l in fh:
    l=l.strip().split()
    if(l[0].count(gene)!=0):
      start = l[4][0:3]
      ind = l[1].index(l[4])
      end = l[1][ind-3:ind]
      cdr3_ends[end] = 3
      cdr3_ends[start] = 0
      cdr3_ends_near["."+start[1:3]] = 0
      cdr3_ends_near[start[0]+"."+start[2]] = 0
      cdr3_ends_near[start[0:2]+"."] = 0
      cdr3_ends_near["."+end[1:3]] = 3
      cdr3_ends_near[end[0]+"."+end[2]] = 3
      cdr3_ends_near[end[0:2]+"."] = 3
  fh.close()
  return(cdr3_ends,cdr3_ends_near)

def Get_ids_test(file):
  fh=open(file,"r")
  ids_test = {}
  for l in fh:
    l=l.strip().split()
    ids_test[l[0]] = [l[1],l[2]]
  fh.close()
  return(ids_test)

def Assign_sequences(Seq_file,gene,annot_file, tmp_file,species,protein_file,refv, refj,loc,CDR3_end_file,codon_file):
  #ids_test = Get_ids_test("/lustre/scratch108/viruses/rbr1/MISEQ/AUTOIMMUNITY_AAV_SLE_HEALTHY1/ORIENTATED_SEQUENCES/ANNOTATIONS/TMP/Incorrect_MERGE_ISO_50_H1.txt")
  head = "#ID\tIGHV\tFR1s\tFR1e\tCDR1s\tCDR1e\tFR2s\tFR2e\tCDR2s\tCDR2e\tFR3s\tFR3e\tIGHV_end\tIGHJ\tIGHJ_start\tCDR3s\tCDR3e\tCDR3_prot\tCDR3_nn\tV_mm\tJ_mm\tV_ins\tJ_ins\n"
  (prots)=Get_protein_sequences(protein_file)
  cdr3_ends,cdr3_ends_near = Get_CDR3_ends(CDR3_end_file,gene)
  fh=open(annot_file, "w")
  fh.write(head)
  fh.close()
  word = 5
  (samev,lengthv)=Get_library(refv,"end",50)
  (samej,lengthj)=Get_library(refj,"start",30)
  (regions,CDR3_end)=Get_region_boundaries(species,loc)
  codon=Get_codons (codon_file)
  fh = open(Seq_file,"r")
  out,out1,t,ind,batch_number,batch, t1='','',0,0,0, 100, time.time()
  seqs,failed,indent={},0,80
  for header,seq in fasta_iterator(fh):
    #if(header.split("__")[0] in prots):
    if(1==1):
    #if(header.split("__")[0] in ids_test):
      (t,ind) = (t+1,ind+1)
      header = header.split("__")[0]
      out = out+">"+header+"\n"+seq+"\n"
      out1=out1+">"+header+"\n"+seq[len(seq)-indent:len(seq)]+"\n"
      seqs[header] = seq
      if(ind >= batch):
        batch_number = batch_number+1
        (passesv,passesj) = Blast_match(out, refv,refj, tmp_file, out1, indent,seqs,lengthv,prots,codon,gene,cdr3_ends,cdr3_ends_near)
        (failed)=Get_annotation(passesv,passesj, seqs, annot_file,lengthj,regions,failed,CDR3_end,prots,codon,gene)
        del seqs,passesv,passesj
        seqs,passesv,passesj = {},{},{}
        ind , out,out1=0,'',''
        print "\rBatch number:",batch_number, "\tSequences completed:",t,"\t% accepted", (t-failed)*100.0/t,"%",
  if(len(seqs)>0):
    (passesv,passesj) = Blast_match(out, refv,refj, tmp_file, out1,indent,seqs,lengthv,prots,codon,gene,cdr3_ends,cdr3_ends_near)
    (failed)=Get_annotation(passesv,passesj, seqs, annot_file,lengthj,regions,failed, CDR3_end,prots,codon,gene)
  t2= time.time()
  fh.close()
  print "\n",Seq_file+"\t"+gene+"\tV\t"+str(t)+"\t"+str(round((t2-t1),3))+"\t sec\t"+str(round((t2-t1)*100000.0/(60*t),3))+"\tmins/100,000\t"+str(failed)+"\tFAILED"
  os.system("rm "+tmp_file+" "+tmp_file+"_jblast "+tmp_file+"_vblast")
  return()

def Get_clusters(cluster_file):
  fh=open(cluster_file, "r")
  clusters={}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      clusters[l[2].split("__")[0]]=l[1]
  fh.close()
  return(clusters)

def Get_max_cluster(seqs, max_id,clusters):
  max_seq,max_aa = seqs[max_id][1], seqs[max_id][2]
  max_v, max_j = seqs[max_id][3], seqs[max_id][4]
  max_c = clusters[max_id.split("__")[0]]
  l1 = len(max_seq)
  max_cluster={}
  max_cluster[max_id]=1
  mismatch = 6
  for id in seqs:
    if(id!=max_id):
      if(clusters[id.split("__")[0]] ==max_c):max_cluster[id]=max_c
      else:
        max_cluster[id]=clusters[id.split("__")[0]]
        seq2, seqs2aa, v, j = seqs[id][1], seqs[id][2], seqs[id][3], seqs[id][4]
        l2 = len(seq2)
        if(len(max_aa) in [len(seqs2aa)-1, len(seqs2aa), len(seqs2aa)+1] and v==max_v and j==max_j):
          if(max_seq==seq2):
            max_cluster[id]=max_c
          else:
            (s1, s2, p1,trimmed)=Trim_sequences(max_seq,seq2,l1,l2, mismatch)
            if (p1==1 and len(s1)>l1-6):
              if(s1==s2):
                max_cluster[id]=max_c
              else:
                (p,mm)=Get_diff(s1,s2,mismatch)
                if(p==1):
                  if(mm<=mismatch):
                    max_cluster[id]=max_c
  return(max_cluster)

def Get_diff(s1,s2,mismatch):
  p=1
  (mm)=Do_counting(s1, s2,mismatch)
  if (mm > mismatch):
    p=0
  return (p,mm)

def Do_counting(s1, s2,mismatch):
  i1=0
  i2=0
  mm=0
  l=min([len(s1),len(s2)])
  for i in range(0,l-2):
    if (s1[i1]==s2[i2]):
      i1=i1+1
      i2=i2+1
    else:
      mm=mm+1
      i2=i2+1
      i1=i1+1
    if (mm>mismatch+1):
      break
  return (mm)

def Trim_sequences(s1,s2,l1,l2,mismatch):
  p,indent=0,6
  trimmed = [0,0]
  stretch = (l1-(2*indent))/mismatch
  i = indent
  for ind in range(0,10):
    if(i>l1-stretch-1):
      break
    sample1=s1[i:i+stretch]
    index=s2.find(sample1)
    if (index!=-1):
      if(index>i):
        if(index-i <=20):s2,trimmed[1]=s2[index-i:l2],index-i
        else:
          if(i-index <=20):s1,trimmed[0]=s1[i-index:l1],i-index
        min_len=min([len(s1),len(s2)])
        if ((max([len(s1),len(s2)]) - min_len) <15):
          s1=s1[0:min_len]
          s2=s2[0:min_len]
          p=1
          break
    else:
      p=0
      i=i+stretch
  return (s1, s2, p,trimmed)

def Refine_CDR3(annot_file, refined_file,cluster_file):
  clusters=Get_clusters(cluster_file)
  fh=open(refined_file,"w")
  fh.close()
  fh=open(annot_file,"r")
  annot,seqs={},{}
  max_f, max_id = 0,''
  out, ind = '',0
  for l in fh:
    l=l.strip()
    if(l[0]=="#"):
      head = l
    else:
      l=l.split()
      if(len(l)>=19):
        if(l[15].count("CDR3_UNCERTAIN")!=1):
          if(l[15]!="CDR3_UNDETERMINED"):
            CDR3 = l[17]
            if(len(CDR3)<25 or (CDR3[0]=="C" and CDR3[len(CDR3)-1]=="W")):
              if(len(CDR3)>28):
                if(CDR3.count("CA")>1):
                  CDR3 = CDR3[CDR3.index("CA"):len(CDR3)]
              o=l[0]
              for i in range(1,19):
                o=o+"\t"+l[i]
              annot[l[0]]=o
              freq = int(l[0].split("__")[1])
              seqs[l[0]]=[freq, l[18],l[17],l[1][0:5], l[13][0:4]]
              if(max_f<freq):max_f, max_id =freq, l[0]
  fh.close()
  if(len(seqs)>0):
    (max_cluster) = Get_max_cluster(seqs, max_id,clusters)
    out, ind =head+"\tCluster\tCDR3_length\tsequence_frequency\tcluster_type\n",0
    fh=open(refined_file,"w")
    fh.close()
    for id in annot:
      cluster = "NON_MAX_CLUSTER"
      if(id in max_cluster):
        cluster = "MAX_CLUSTER"
      out=out+annot[id]+"\t"+str(max_cluster[id])+"\t"+str(len(seqs[id][2]))+"\t"+str(id.split("__")[1])+"\t"+cluster+"\n"
      ind = ind+1
      if(ind>500):
        Write_out(out,refined_file)
        out, ind = "",0
    Write_out(out,refined_file)
  return()

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

def Get_network_statistics_per_chain(cluster_file, sample, dir,per_chain_repertoire_statistics_file,isotyper_primer_set):
  fh=open(cluster_file,"r")
  cluster,vertices=Tree(), Tree()
  index, totalc, totalv, totalreads, sizesv, c_sizes,vertices_in_max_cluster = 0,0,0,0,[],{},0
  total_v, total_reads = [],[]
  sizesv, c_sizes= {},{}
  chains_short = []
  t1,t2 = 0,0
  n = 0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      id = l[2]
      chains, freq, id_short = id.split("|")[1].split("_"), map(int, id.split("__")[1].split("|")[0].split("_")), id.split("__")[0]
      t1 = t1+sum(freq)
      n = n+1
      if(len(chains_short)==0):
        for i in range(0,len(chains)):
          c = chains[i].split("*")[0]
          if(isotyper_primer_set =="INNER_DD"):
            if(c in ["IGHA1","IGHA2"]):c = "IGHA"
            elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
          if(c not in chains_short):chains_short.append(c)
          #if(c in chains_short):chains_index[chains[i].split("*")[0]] = chains_short.index(c)
          #else:
          #  chains_short.append(c)
          #  chains_index[chains[i].split("*")[0]] = chains_short.index(c)
      non_zero = [i for i in range(len(freq)) if freq[i]!=0]
      if(len(total_v)==0):
        total_v, total_reads = [0]*len(chains_short),[0]*len(chains_short)
        for c in chains_short:
          sizesv[c], c_sizes[c] = [],[]
      for i in non_zero:
        c = chains[i].split("*")[0]
        if(isotyper_primer_set =="INNER_DD"):
          if(c in ["IGHA1","IGHA2"]):c = "IGHA"
          elif(c in ["IGHG1","IGHG2"]):c = "IGHG1/2"
        cluster[c][l[1]][freq[i]][id_short].value=1
        index = chains_short.index(c)
        total_v[index] = total_v[index]+1
        sizesv[c] =sizesv[c]+[freq[i]]
        total_reads[index] =total_reads[index]+freq[i]
  fh.close()
  print total_reads, t1, n
  if(t1 != sum(total_reads)):print "ERROR IN COUNTING!!"
  out="#Id\tIsotype\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\n"
  #out="#Id\tAnalysis\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\t% Vertices in largest cluster\tVertex Renyi\tCluster Renyi\tGene\tSpecies\n"
  for c1 in chains_short: 
    cluster_sizes_sub = []
    for clus in cluster[c1]:
      f = 0
      for f1 in cluster[c1][clus]:
        f = f+(f1*len(cluster[c1][clus][f1]))
      cluster_sizes_sub = cluster_sizes_sub+[f]
    if(len(cluster_sizes_sub)>0):
      (vpoints,vvdf)=VDF(sizesv[c1])
      (cpoints,cvdf)=VDF(cluster_sizes_sub)
      vgini, cgini=Gini_index(cpoints,cvdf, vpoints,vvdf,sum(cluster_sizes_sub),sum(sizesv[c1]),total_v)
      max_pop, max_1_pop = cpoints[len(cpoints)-1]*100.0/sum(sizesv[c1]), cpoints[len(cpoints)-2]*100.0/sum(sizesv[c1])
      out = out+str(sample)+"\t"+c1+"\t"+str(sum(sizesv[c1]))+"\t"+str(len(sizesv[c1]))+"\t"+str(vgini)+"\t"+str(cgini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\n"
  fh=open(per_chain_repertoire_statistics_file, "w")
  fh.write(out)
  fh.close()
  return()


def Get_cluster_vertex_distributions(cluster_file):
  fh=open(cluster_file,"r")
  cluster,vertices=Tree(), Tree()
  index, totalc, totalv, totalreads, sizesv, c_sizes,vertices_in_max_cluster = 0,0,0,0,[],{},0
  for l in fh:
    index=index+1
    if (index>1):
      l=l.strip()
      l=l.split()
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

def Print_distributions(points,cdf, file_out):
  out = "#size\tfrequency of size\n"
  for i in range(0,len(points)):
    out=out+str(points[i])+"\t"+str(cdf[i])+"\n"
  fh=open(file_out,"w")
  fh.write(out)
  fh.close()
  return()

def Get_network_statistics(cluster_file,id,dir,network_statistics,gene,species,cluster_size_distribution,vertex_size_distribution):
  (cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads,c_sizes,vertices_in_max_cluster)=Get_cluster_vertex_distributions(cluster_file)
  Print_distributions( vpoints,vvdf, vertex_size_distribution)
  Print_distributions(cpoints,cvdf,cluster_size_distribution)
  (vrenyi,crenyi)=Renyi_entropy(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads)
  (vgini, cgini)=Gini_index(cpoints,cvdf, vpoints,vvdf,totalc,totalv,totalreads)
  (max_pop, max_1_pop)=Proportional_measures(c_sizes, totalreads)
  out="#Id\tAnalysis\tN reads\tN vertices\tVertex Gini Index\tCluster Gini Index\tLargest Cluster (%)\t2nd Largest Cluster (%)\t% Vertices in largest cluster\tVertex Renyi\tCluster Renyi\tGene\tSpecies\n"
  out=out+str(id)+"\tOVERALL\t"+str(totalreads)+"\t"+str(totalv)+"\t"+str(vgini)+"\t"+str(cgini)+"\t"+str(max_pop)+"\t"+str(max_1_pop)+"\t"+str(vertices_in_max_cluster*100.0/totalv)+"\t"+str(vrenyi)+"\t"+str(crenyi)+"\t"+gene+"\t"+species+"\n"
  fh=open(network_statistics, "a")
  fh.write(out)
  fh.close()
  return()

def Get_large_clusters(cluster_file,threshold,max_number_of_clusters_to_include):
  fh=open(cluster_file,"r")
  total,clusters,clust_size=0,Tree(),{}
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
  sig_clust_info =sorted(sig_clust_info,key=itemgetter(1),reverse=True)
  if(len(sig_clust_info)>max_number_of_clusters_to_include):
    ids_new,sig_clust_info_new = {},[]
    for i in range(0,max_number_of_clusters_to_include):
      c = sig_clust_info[i][0]
      for id in clusters[c]:
        ids_new[id]=c
      sig_clust_info_new.append(sig_clust_info[i])
    return(sig_clust_info_new, ids_new, clusters)
  else:
    return(sig_clust_info,ids, clusters)

def Get_ids(header_file, seq_file_all,id,seq_file):
  fh=open(header_file,"r")
  index_used = -1
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[1]==id):
        index_used = int(l[0])
  fh.close()
  print index_used
  ids= {}
  if(index_used!=-1):
    fh=open(seq_file_all,"r")
    for l in fh:
      if(l[0]==">"):
        id = l.strip().replace(">","")
        f = id.split("__")[1].split("_")[index_used]
        if(f!="0"):
          ids[id.split(":")[0]]=int(f)
          #print id.split(":")[0]
          #if(len(ids)>10000):break ############## remove
    fh.close()
  fh=open(seq_file,"r")
  for header,seq in fasta_iterator(fh):
    ids[header.split(":")[0]]=int(header.split("__")[1])
  fh.close()
  print "IDS read", len(ids)
  return(ids)

def Get_annotation_for_clusters(annot_file, ids):
  ids_find ={}
  for id in ids:
    ids_find[id.split("__")[0]]=id
  fh=open(annot_file,"r")
  cluster_annot,cluster_mutations,CDR3s=Tree(),{},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(l[0].split("__")[0] in ids_find):
        c = ids[ids_find[l[0].split("__")[0]]]
        if(len(l)>=14):
          v,j = l[1],l[13]
          if(j.count("J")==0):print l
          vj= l[1]+"\t"+l[13]
          cluster_annot[c][vj][l[0]].value=1
          if(len(l)>=20 and l[19].count("CDR")==0):
            CDR3s[l[0].split("__")[0]]=l[17]
            if(c in cluster_mutations):cluster_mutations[c] = copy.deepcopy(cluster_mutations[c])+[int(l[19])]
            else:cluster_mutations[c] =[int(l[19])]
  cluster_annotation = {}
  for c in cluster_annot:
    max_vj,max_f_vj = '',0
    for vj in cluster_annot[c]:
      if(len(cluster_annot[c][vj])>max_f_vj):max_vj,max_f_vj = vj, len(cluster_annot[c][vj])
    if(c in cluster_mutations):
      m = cluster_mutations[c]
      mean_muts = sum(m)*1.0/len(m)
    else: mean_muts = "NA"
    cluster_annotation[c] = [max_vj, mean_muts]
  del cluster_annot, cluster_mutations
  return(cluster_annotation,CDR3s)

def Get_primers_from_reference (gene,species,loc):
  #if(species=="HOMO_SAPIENS"):fh=open(loc+"LIBRARY/Primer_reference_"+species+"_"+gene+"V.txt","r") 
  #if(species=="MUS_MUSCULUS"):fh=open(loc+"LIBRARY/MUS_MUSCULUS_primers.txt","r")
  if(species=="HOMO_SAPIENS"):fh=open(loc+"LIBRARY/Primers_HOMO_SAPIENS_IG_RBR_Constant_region_MPLX_table.txt", "r")
  if(species=="MUS_MUSCULUS"):fh=open(loc+"LIBRARY/Primers_MUS_MUSCULUS_IG_RBR_Constant_region_MPLX_table.txt", "r")
  primer_reference={}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      primer_reference[l[0]]=l[1]
  fh.close()
  return(primer_reference)

def Get_max_sequences(seq_file,ids):
  ids_find= {}
  for id in ids:
    ids_find[id.split("__")[0]] = id
  fh=open(seq_file,"r")
  seqs={}
  for header,sequence in fasta_iterator(fh):
    if(header.split("__")[0] in ids_find):
      seqs[header.split("__")[0]]=sequence
  fh.close()
  return(seqs)

def Get_cluster_stats(annot_file, cluster_file, id, dir,network_statistics,gene,species,seq_file,loc):
  if(gene in "IGH" and species=="HOMO_SAPIENS"):primer_reference = Get_primers_from_reference (gene,species,loc)
  else:primer_reference={}
  threshold = 1.0 ### the minimum % of reads in cluster to be reported
  sig_clust_info,ids,clusters =Get_large_clusters(cluster_file,threshold,100000000000000000000000000000)
  del clusters
  cluster_annotation,CDR3s = Get_annotation_for_clusters(annot_file, ids)
  seqs = Get_max_sequences(seq_file,ids)
  out="#Id\tAnalysis\tCluster ID\t% of Repertoire\tN Reads\tN Vertices\tV gene\tJ gene\tMean mutations in V\tMax Sequence Size (reads)\tMax Sequence ID\tGene\tSpecies\tSpecific Primer\tLargest vertex sequence\n"
  if(len(sig_clust_info)>0):
    for i in range(0,len(sig_clust_info)):
      c,perc_network, number_vertices, number_reads, max_seq, reads_in_max_seq = sig_clust_info[i]
      if(max_seq.split("__")[0] in seqs):seq = seqs[max_seq.split("__")[0]]
      else:seq = max_seq
      if(c in cluster_annotation):max_vj, mean_muts = cluster_annotation[c]
      else:max_vj, mean_muts = "NA\tNA","NA"
      if(max_vj.split("\t")[0] in primer_reference):primer_sequence = primer_reference[max_vj.split("\t")[0]]
      else:primer_sequence="NA"
      out=out+id+"\tCLUSTER_ANALYSIS\t"+str(c)+"\t"+str(perc_network)+"\t"+str(number_reads)+"\t"+str(number_vertices)+"\t"+max_vj+"\t"+str(mean_muts)+"\t"+str(reads_in_max_seq)+"\t"+str(max_seq)+"\t"+gene+"\t"+species+"\t"+primer_sequence+"\t"+seq+"\n"
  else:
    out=out+id+"\tCLUSTER_ANALYSIS\tNO CLUSTERS LARGER THAN "+str(threshold)+"% OF REPERTOIRE\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t"+gene+"\t"+species+"\tNA\tNA\n"
  fh=open(network_statistics, "a")
  fh.write(out)
  fh.close()
  return()

def Get_unique_protein_sequence_number(protein_file,id):
  seqs = {}
  fh=open(protein_file,"r")
  for header,sequence in fasta_iterator(fh):
    seqs[sequence ]=1
  fh.close()
  print id+"\t"+str(len(seqs))
  return()

def Get_CDR3_lengths(annot_file, CDR3_length_file,id,dir):
  lengths = Tree()
  tot,found = 0,0
  fh=open(annot_file,"r")
  for l in fh:
    if(l[0]!="#"):
      tot = tot+1
      l=l.strip().split()
      if(len(l)>=18):
        CDR3 = l[17]
        if(len(CDR3)>4 and CDR3.count("-")==0):
          lengths[l[0].split(":")[0]+"\t"+str(len(CDR3))][CDR3+":"+l[0]].value = 1
          found = found+1
  fh.close()
  out="#CDR3 length\tFrequency\n"
  for l in lengths:
    out=out+str(l)+"\t"+str(len(lengths[l]))+"\n"
  fh=open(CDR3_length_file,"w")
  fh.write(out)
  fh.close()
  os.system("cp "+CDR3_length_file+" ~rbr1/OUTPUT/")
  print tot,found, found*100.0/tot,"%"
  return()

def Get_CDR3_lengths_of_interest(annot_file, CDR3_length_file,sample,dir):
  lengths = Tree()
  fh=open(annot_file,"r")
  tot,found = 0,0
  lengths_of_interest = [10,13,15,22]
  CDR3_print,CDR3_interest_IDS,CDRs = Tree(),{},{}
  for l in fh:
    if(l[0]!="#"):
      tot = tot+1
      l=l.strip().split()
      if(len(l)>=18):
        CDR3 = l[17]
        if(len(CDR3)>4 and CDR3.count("-")==0):
          lengths[len(CDR3)][CDR3+":"+l[0]].value = 1
          found = found+1
          if(len(CDR3) in lengths_of_interest):
            CDR3_print[len(CDR3)][CDR3][l[0]].value = 1
            CDR3_interest_IDS[l[0]] = len(CDR3)
            CDRs[CDR3]=1
  fh.close()
  llama = "46"
  if(sample in ["01","03","05","09"]):llama = "37"
  prot_cluster_file = dir.replace("ANNOTATIONS/","TEMPORAL_PROTEIN/")+"Protein_Cluster_IDs_"+llama+".txt"
  seq_proteins_all = dir.replace("ANNOTATIONS/","TEMPORAL_PROTEIN/")+"Protein_Fully_reduced_"+llama+".fasta"
  fh=open(seq_proteins_all,"r")
  aliases = {}
  for header,sequence in fasta_iterator(fh):
    for CDR in CDRs:
      if(sequence.count(CDR)!=0):
        aliases[header.split("_")[0]] = CDR
        print len(aliases)
        break
  fh.close()

def Get_gene_frequencies(annot_file, gene_freq_file,gene,id):
  fh=open(annot_file,"r")
  genes,total,found = {},0,0
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      if(len(l)>=13):
        if(l[0].count("__")!=0):
          f, v, j = sum(map(int,l[0].split("__")[1].split("|")[0].split("_"))), l[1].split("*")[0], l[13].split("*")[0]
        else:f, v, j = 1, l[1].split("*")[0], l[13].split("*")[0]
        total = total+f
        if(v.count("not_found")==0 and j.count("J")!=0):
          v=v+"|"+j # change back to v+j for complete list
          if(v in genes):genes[v]=genes[v]+f
          else:genes[v]=f
          found = found+f
        else:
          print l
      elif(len(l)==9):
        if(l[0].count("__")!=0):
          f, v, j = sum(map(int,l[0].split("__")[1].split("|")[0].split("_"))), l[1].split("*")[0], l[3].split("*")[0]
        #if(v.count(gene)!=0 and j.count(gene)!=0):
        if(v.count("not_found")==0 and j.count("J")!=0):
          v=v+"|"+j # change back to v+j for complete list
          if(v in genes):genes[v]=genes[v]+f
          else:genes[v]=f
          found = found+f
      else:
        print l
  fh.close()
  print "TOTAL READS:",total,"FOUND READS:", found
  out=''
  for g in genes:
    out=out+id+"\t"+g+"\t"+str(genes[g])+"\t"+g[0:5]+"\n"
  fh=open(gene_freq_file, "w")
  fh.write(out)
  fh.close()
  return()

def Initalise_file(file):
  fh=open(file,"w")
  fh.close()
  return()

def Write_out(out, file):
  fh = open (file,"a")
  fh.write(out)
  fh.close()
  return()

def Intialise_files(dir):
  dirs_to_add=[dir, dir+"TMP/"]
  for d in dirs_to_add:
    c=commands.getoutput("ls "+d)
    if(c.count("No such file or directory")==1):commands.getoutput("mkdir "+d)
  return()

###########################
dir = sys.argv[1]             ### output directory
id = sys.argv[2]              ### identifier for output files
seq_file= sys.argv[3]         ### input fasta
protein_file = sys.argv[4]    ### protein fasta (does not matter if read number index is not correct)
gene = sys.argv[5]            ### either IGH, IGL, IGK, TRA, TRB, TRD, TRG
species = sys.argv[6]         ### HOMO_SAPIENS...
cluster_file = sys.argv[7]    ### cluster file (output from Read_processing_and_quality.py program)
command = sys.argv[8]         ### Comma separated list: ANNOTATE = Generate raw V-J annotation files, STATISTICS = Get population statistics
forward_primer_group = sys.argv[9] ### Forward primers used (basically to disinguish between ISO_DD and IsoTyper)

#seq_file = seq_file.replace("NETWORKS/Fully_reduced","Nucleotide_ORF_filtered_all")
########################### Files for annotation generation
edge_file = dir.replace("ANNOTATIONS/","")+"NETWORKS/Edges_"+id+".txt"
annot_file = dir+"TMP/Annotation_"+id+".txt"
tmp_file =dir+"TMP/Tmp_annotation_"+id
refined_file = dir+"Refined_annotation_"+id+".txt"
network_statistics = dir+"Network_statistics_"+id+".txt"
cluster_statistics = dir+"Cluster_statistics_"+id+".txt"
gene_freq_file = dir+"Gene_frequencies_"+id+".txt"
gene_freq_file2 = dir+"Gene_IMGT_frequencies_"+id+".txt"
CDR3_frequency_file = dir+"TMP/CDR3_frequencies_"+id+".txt"
cluster_properties_file = dir+"TMP/Cluster_properties_file_"+id+".txt"
clonal_selection_file = dir+"TMP/Clonal_selection_values_"+id+".txt"
pre_secondary_file =  dir+"TMP/Pre_pre_secondary_rearrangement_"+id+".txt"
mutational_file = dir+"Mutations_VJ_"+id+".txt"
#secondary_rearrangement_file = dir+"TMP/Secondary_rearrangement_"+id+".txt"
#CDR3_length_and_V_region_genomic_position_file = dir+"TMP/CDR3_length_and_V_region_genomic_position_"+id+".txt"
secondary_rearrangement_file_output = dir+"Secondary_refined_rearrangement_"+id+".txt"
CDR3_length_file = dir+"CDR3_lengths_"+id+".txt"
file_D_genes = dir+"TMP/D_gene_annotations_"+id+".txt"
positive_residue_count_file = dir+"Positive_residue_counts_"+id+".txt"
region_mutational_file = dir+"Region_mutational_counts_"+id+".txt"
syn_to_non_syn_file = dir+"Non_syn_to_syn_counts_"+id+".txt"
low_mutational_vj_gene_usage = dir+"Low_mutational_vj_gene_usage_"+id+".txt"
IMGT_tmp_file = dir+"IMGT_tmp_file_"+id+".txt"
constant_region_count_file = dir+"Constant_region_counts_"+id+".txt"
cluster_size_distribution = dir+"Distribution_cluster_sizes_"+id+".txt"
vertex_size_distribution = dir+"Distribution_vertex_sizes_"+id+".txt"
seq_length_ditribution = dir+"Distribution_sequence_lengths_"+id+".txt"
per_chain_repertoire_statistics_file = dir+"IsoTyper_chain_repertoire_statistics_file_"+id+".txt"
germline_equiv_file =  dir+"Germline_equivalent_file_"+id+".fasta"
largest_cluster_file = dir+"LARGEST_CLUSTER_ANALYSIS/Largest_cluster_"+id
annot_file_IMGT = dir+"IMGT_ISOTYPER/IMGT_"+id+"_1_Summary.txt"
annot_file_IMGT3 = dir+"IMGT_ISOTYPER/IMGT_"+id+"_3_Nt-sequences.txt"
annot_file_IMGT5 = dir+"IMGT_ISOTYPER/IMGT_"+id+"_5_AA-sequences.txt"
IMGT_QC_file = dir+"IMGT_ISOTYPER/QC_IMGT_"+id+".txt"
mutational_distribution_file = dir+"Number_of_mutations_"+id+".txt"
CDR3_file = dir+"CDR3_sequences_"+id+".txt"
CDR3_cluster_file = cluster_file.replace("Cluster_identities_","Cluster_CDR3_defined_identities_")
custom_reference_blast_file = dir+"Custom_reference_blast_file_"+id+".txt"
pre_seq_file = dir.replace("ANNOTATIONS/","")+"TMP/Trimmed_orientated_all_"+id+".fasta"
subsample_file = dir+"Subsample_file_"+id+".txt"
overall_clonal_ids_file = dir+"Cluster_information_total_"+id+".txt"
########################## Reference files
loc = "/nfs/users/nfs_k/kt16/BCRSeq/"
refv = loc+"LIBRARY/Reference_nn_"+species+"_"+gene+"V.fasta"
refd = loc+"LIBRARY/Reference_nn_"+species+"_"+gene+"D.fasta"
refj = loc+"LIBRARY/Reference_nn_"+species+"_"+gene+"J.fasta"
refvp= loc+"LIBRARY/Reference_protein_"+species+"_"+gene+"V.fasta"
refjp= loc+"LIBRARY/Reference_protein_"+species+"_"+gene+"J.fasta"
CDR3_end_file = loc+"LIBRARY/CDR3_end_"+species+"_imgt.txt"
cRSS_file = loc+"LIBRARY/cRSS_IgHV_sequences.txt"
cRSS_footprint_file = loc+"LIBRARY/cRSS_IGHV_footprint_motifs.txt"
codon_file = loc+"LIBRARY/Codon_table2.txt"
########################## Get commands
command = command.split(",")

########################################################################remove when done
if(forward_primer_group.upper().count("ISO")!=0):constant_region = "TRUE"
else:constant_region = "FALSE"
#constant_region = "FALSE"
constant_region = "TRUE"


isotyper_primer_set = "INNER"
if(constant_region=="TRUE"):
  if(forward_primer_group=="ISO_DD"):isotyper_primer_set = "INNER_DD"
  else:isotyper_primer_set = "INNER"

print "constant_region:",constant_region,"isotyper_primer_set:", isotyper_primer_set

######################### Commands
Intialise_files(dir)
if("ANNOTATE" in command):
  if(species in ["HOMO_SAPIENS","MACACA_MULATTA","MUS_MUSCULUS","LLAMA_GLAMA"]):
    Assign_sequences(seq_file,gene,annot_file, tmp_file,species, protein_file,refv, refj,loc,CDR3_end_file,codon_file)
  Get_gene_frequencies(annot_file, gene_freq_file,gene,id)
if("ANNOTATE_B" in command):
  Get_sequence_length_distributions(seq_file,seq_length_ditribution)
  Refine_CDR3(annot_file, refined_file, cluster_file)
  if(gene not in ["TCR"]):Get_CDR3_lengths(annot_file, CDR3_length_file,id,dir)
  Get_CDR3_lengths_of_interest(annot_file, CDR3_length_file,id,dir)
  Get_unique_protein_sequence_number(protein_file,id)
  Identify_D_genes(annot_file, file_D_genes, refd,tmp_file,seq_file)
if("STATISTICS" in command):
  Initalise_file(network_statistics)
  Initalise_file(cluster_statistics)
  Get_network_statistics(cluster_file, id, dir,network_statistics,gene,species,cluster_size_distribution,vertex_size_distribution)
  Get_cluster_stats(annot_file, cluster_file, id, dir,cluster_statistics,gene,species,seq_file,loc)
  if(constant_region == "TRUE"):Get_network_statistics_per_chain(cluster_file, id, dir,per_chain_repertoire_statistics_file,isotyper_primer_set)