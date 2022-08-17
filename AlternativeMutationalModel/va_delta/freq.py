from Bio import AlignIO
import pandas as pd
import subprocess
import os
import numpy as np
import glob
import sys
from Bio import SeqIO

#calculates threshold for nucleotide change based on shannon
#column entropy

def column_entropy_thresh(freq_df):
	e_act = 0
	e_max = 0
	for i in range(len(freq_df.index)):
		p_xi = freq_df.iloc[i]
		e_act += p_xi*np.log(p_xi)

		p_xm = 1/len(freq_df.index)
		e_max += p_xm*np.log(p_xm)

	thresh = 1-(e_act/e_max)

	if np.isnan(thresh) == True:
		thresh = 0


	return thresh

#Determine nulceotide change
def determine_change(thresh):
	change = []
	for threshold in thresh:
		val = np.random.randint(0,100)
		if val < threshold:
			change.append('Yes')
		else:
			change.append('No')
	return change

def commit_change(index, change):
	new_seq = []
	for (nucleotide, change_val) in zip(index, change):
		if change_val == 'Yes':
			new_nucleotide = np.random.randint(1,5)
			if new_nucleotide == 1:
				new_nucleotide = 'A'
				new_seq.append(new_nucleotide)
			if new_nucleotide == 2:
				new_nucleotide = 'T'
				new_seq.append(new_nucleotide)
			if new_nucleotide == 3:
				new_nucleotide = 'G'
				new_seq.append(new_nucleotide)
			if new_nucleotide == 4:
				new_nucleotide = 'C'
				new_seq.append(new_nucleotide)
			if new_nucleotide == 5:
				new_nucleotide ='-'
				new_seq.append(new_nucleotide)
		else:
			new_nucleotide = nucleotide
			new_seq.append(new_nucleotide)

	new_seq = ''.join(new_seq)
	return new_seq

def add_to_fasta(new_seq, node):
	ofile = open("my_fasta.txt", "a")
	ofile.write(">" + str(node) + "\n" +new_seq + "\n")
	ofile.close()

def find_seq(node):
	records = list(SeqIO.parse("my_fasta.txt", "fasta"))
	rec = 0
	for i in range(0,len(records)):
		print(records[i].name)
		print(node)
		if records[i].name == str(node):
			rec = records[i]
	return rec

##########################################################
#read in alignment to pandas dataframe
print('reading alignment file into pandas dataframe.....')
align = AlignIO.read('delta.fa', 'fasta')
name = []
description = []
for record in align:
	name.append(record.name)
	description.append(record.description)
align2 = pd.DataFrame(align)

#create threshold list, each column threshold included
print('calculating entropy and setting the threshold...')
thresh = []
df = pd.DataFrame()
for i in range(len(align2.columns)):
	df1 = align2.iloc[:,i].value_counts()
	df1 = df1.divide(len(align2.index))
	thresh.append(column_entropy_thresh(df1)*100)

##########################################################
import networkx as nx
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt

#Read in network data
print('reading in the network data....')
df = pd.read_csv('output_abridges.csv')

#Create connection dataframe (pid, and contact_pid columns)
connections1 = df.iloc[:,1]
connections2 = df.iloc[:,3]
connections = pd.concat([connections1, connections2], axis = 1, ignore_index = False)

#Create id dataframe (with tick and exit_state columns)
id1 = df.iloc[:,0]
id2 = df.iloc[:,2]
ids = pd.concat([id1, id2], axis = 1)


#Assigning ticks to correct nodes for coloring based on time
timelist = dict(zip(connections1, id1))
pos = list(timelist.keys()).index(476724)
items = list(timelist.items())
items.insert(pos, (-1, 0))
timelist = dict(items)
colors_list = list(timelist.values())


#draw directed graph network
print('creating network digraph .......')
G=nx.from_pandas_edgelist(df, 'pid', 'contact_pid', create_using=nx.DiGraph())

##########################################################

#creating .fasta file to append
open('my_fasta.txt', 'w')

#nodes = keys, values = predecessors

#determine all seeds
print('determining all seed nodes......')
seeds = []
for (node, pred) in zip(connections1, connections2):
	if pred == -1:
		seeds.append(node)

#walk along every seed, determine edges via depth first search
i=0
for seed in seeds:
	print('determining edges list.......')
	edges = list(nx.dfs_edges(G, source = seed))
	print(edges)
	if len(edges) > 1:
		for edge in reversed(edges):
			if edge[1] == -1:
				print('Determining if sequence is already in fasta......')
				record_val = find_seq(edge[0])
				if str(record_val) == str(0):
					print('Adding seed seq to .fasta ........')
					index = align2.iloc[i].values.tolist()
					index = ''.join(index)
					i+=1
					add_to_fasta(index, str(edge[0]))
			if edge[1] != -1:
				record_val = find_seq(edge[0])
				if str(record_val) == str(0):
					print('Mutating sequence, adding to fasta.....')
					seq_to_change = find_seq(edge[1])
					change = determine_change(thresh)
					print(edge[1])
					new_seq = commit_change(seq_to_change, change)
					add_to_fasta(new_seq, edge[0])












