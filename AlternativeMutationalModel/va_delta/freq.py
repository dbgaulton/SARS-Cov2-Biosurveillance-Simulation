from locale import currency
from multiprocessing import current_process
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

def poor_mut_model(index, change):
	# Rotate to a static new letter
	new_seq = []
	for (nucleotide, change_val) in zip(index, change):
		if change_val == 'Yes':
			if nucleotide == 'A':
				new_nucleotide = 'T'
				new_seq.append(new_nucleotide)
			if nucleotide == 'T':
				new_nucleotide = 'G'
				new_seq.append(new_nucleotide)
			if nucleotide == 'G':
				new_nucleotide = 'C'
				new_seq.append(new_nucleotide)
			if nucleotide == 'C':
				new_nucleotide == '-'
				new_seq.append(new_nucleotide)
			if nucleotide == '-':
				new_nucleotide == 'A'
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
	records = list(SeqIO.parse("my_fasta.txt", "fasta")) # dgaulton: going to parse this each time, need to modify this for re-infections, could be two entries
	rec = 0
	for i in range(0,len(records)):
		print(records[i].name)
		print(node)
		if records[i].name == str(node):
			rec = records[i]
	return rec

def dfs_edges_with_ticks(G, source=None, depth_limit=None):
	if source is None:
		# edges for all components
		nodes = G
	else:
		# edges for components with source
		nodes = [source]
	visited = set()
	if depth_limit is None:
		depth_limit = len(G)
	
	for start in nodes:
		tick = -1

		if start in visited:
			continue
		visited.add(start + "_" + tick) # modify visitied to include tick - TODO make sure this is treated as a string and not as an int
		stack = [(start, depth_limit, iter(G[start]))]
		while stack:
			parent, depth_now, children = stack[-1]
			try:
				child = next(children) # need to confirm if iter for children in MultiDiGraph will treat the node all at once or if will go edge by edge between the neighbor - might need to flatten the results of iter(G[start]) to treat a new edge with an old neighbor as new - that would occur above at stack = 
				# tick = child.
				if child not in visited:
					yield parent, child
					visited.add(child + "_" + tick)
					if depth_now > 1:
						stack.append((child, depth_now - 1, iter(G[child])))
			except StopIteration:
				stack.pop()
	

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
df = pd.read_csv('contact_network_va_delta.csv')


#Create connection dataframe (pid, and contact_pid columns) TODO introduce time to edge in order to avoid cycles during DFS
connections1 = df.iloc[:,1] # pid
connections2 = df.iloc[:,3] # contact_pid
connections = pd.concat([connections1, connections2], axis = 1, ignore_index = False)

#Create id dataframe (with tick and exit_state columns) TODO tick here includes time?
id1 = df.iloc[:,0] # tick
id2 = df.iloc[:,2] # exit_state
ids = pd.concat([id1, id2], axis = 1) # dgaulton: looks like this is never used


#Assigning ticks to correct nodes for coloring based on time
timelist = dict(zip(connections1, id1))
pos = list(timelist.keys()).index(476724) # dgaulton: where did this number come from? Is this from 20 tick cutoff? - looks like picking out that item and moving it
items = list(timelist.items())
items.insert(pos, (-1, 0))
timelist = dict(items)
colors_list = list(timelist.values())


#draw directed graph network
print('creating network digraph .......')
G=nx.from_pandas_edgelist(df, 'contact_pid', 'pid', create_using=nx.MultiDiGraph(), edge_key='tick') # TODO need to add tick here to use down below - will be edge[2], also maybe need to use MultiDiGraph if repeat infection

##########################################################

#creating .fasta file to append
open('my_fasta.txt', 'w')

#nodes = keys, values = predecessors

#determine all seeds
# print('determining all seed nodes......')
# seeds = []
# for (node, pred) in zip(connections1, connections2):
# 	if pred == -1:
# 		seeds.append(node)

#walk along every seed, determine edges via depth first search
# i=0
# for seed in seeds: # dgaulton: will cause double traversals, could fix by passing in list of seeds to one call to edge_dfs
# 	print('determining edges list.......')
# 	edges = list(nx.edge_dfs(G, source = seed))
# 	print(edges)
# 	if len(edges) > 1:
# 		for edge in edges:
# 			if edge[0] == -1:
# 				print('Determining if sequence is already in fasta......') # dgaulton: might take too long at scale, better to store in hash table in mem. Maybe do this check to avoid cycle
# 				record_val = find_seq(edge[1])
# 				if str(record_val) == str(0):
# 					print('Adding seed seq to .fasta ........')
# 					index = align2.iloc[i].values.tolist()
# 					index = ''.join(index)
# 					i+=1
# 					add_to_fasta(index, str(edge[1]))
# 			if edge[0] != -1:
# 				record_val = find_seq(edge[1])
# 				if str(record_val) == str(0): # dgaulton: if not already in fasta, won't work for second visit for reinfection
# 					print('Mutating sequence, adding to fasta.....')
# 					seq_to_change = find_seq(edge[0])
# 					change = determine_change(thresh)
# 					print(edge[1])
# 					new_seq = commit_change(seq_to_change, change)
# 					add_to_fasta(new_seq, edge[1])


# dgaulton: might need more of a bfs that goes tick by tick to make sure we don't mutate an seq from a prior infection for a later infection
# during bfs, have the ability mutate the most recent seq that didn't occur after current infection

# mutations = {} # map of nodes -> list of a node's mutations tupled with tick at which that mutation occurred

# def get_most_recent_mutation(node, tick, mutations):
# 	if len(mutations[node]) == 0:
# 		return None
	
# 	i = 0
# 	for i in range(len(mutations[node])):
# 		if mutations[node][i][0] > tick:
# 			return mutations[node][i-1]
# 		i += 1

# 	return mutations[node][i]

# i=0
# print('determining edges list.......')
# edges = list(nx.edge_bfs(G, source = seeds))
# print(edges)

# if len(edges) > 1:
# 	for edge in edges:
# 		if edge[0] == -1: # seed case
# 			print('Adding seed seq........')
			
# 			index = align2.iloc[i].values.tolist()
# 			index = ''.join(index)
# 			i+=1
			
# 			if mutations[edge[1]] is None:
# 				mutations[edge[1]] = [(0, index)]
# 			else:
# 				mutations[edge[1]].append((0, index))
						
# 		if edge[0] != -1: # normal case
# 			print('Mutating sequence, adding to fasta.....')
# 			seq_to_change = get_most_recent_mutation(edge[0], edge[2], mutations)

# 			if seq_to_change is None:
# 				print("Error: unexpectedly didn't have mutation for parent")

# 			change = determine_change(thresh)
# 			print(edge[1])
# 			new_seq = commit_change(seq_to_change, change)

# 			if mutations[edge[1]] is None:
# 				mutations[edge[1]] = [(edge[2], new_seq)]
# 			else:
# 				mutations[edge[1]].append((edge[2], new_seq))

# print(mutations)

# ofile = open("my_fasta.txt", "a")
# for node in mutations.keys():
# 	for mutation in mutations[node]:
# 		ofile.write(">" + str(node) + "\n" + mutation + "\n")

# ofile.close()


# def dfs_with_tick(G, sources): # TODO inject tick check
# 	nodes = list(G.nbunch_iter(sources))

# 	visited_edges = set()
# 	visited_nodes = set()
# 	edges = {}

# 	kwds = {
# 		"data": False,
# 		"keys": True
# 	}

# 	# start DFS
# 	for start_node in nodes:

# 		current_tick = 0

# 		index = align2.iloc[i].values.tolist()
# 		index = ''.join(index)
# 		i+=1

# 		if mutations[start_node] is None:
# 				mutations[start_node] = [(current_tick, index)]
# 		else:
# 			mutations[start_node].append((current_tick, index))

# 		stack = [start_node]
		
# 		while stack:
# 			current_node = stack[-1]


			
# 			if current_node not in visited_nodes:
# 				edges[current_node] = iter(G.edges(current_node, **kwds))
# 				visited_nodes.add(current_node)

# 			try:
# 				edge = next(edges[current_node])
# 			except StopIteration:
# 				# No more edges from the current node.
# 				stack.pop()
# 			else:
# 				edgeid = edge
# 				if edgeid not in visited_edges:
# 					visited_edges.add(edgeid)
# 					# Mark the traversed "to" node as to-be-explored.
# 					stack.append(edge[1])
# 					yield edge

# need two data structures
# one list to just append to to generate MSA of all sequences as we go
# another dict that maps node to it's most recent sequence

# This assumes the data read into df is in chronological order (ascending according to tick)
current_sequences = {}
i=0

for pid, contact_pid, tick in zip(connections1, connections2, id1):
	print(tick)
	print(contact_pid)
	print(pid)
	if contact_pid == -1: # seed case
		# grab a new real sequence
		print('Adding seed seq to .fasta ........')
		index = align2.iloc[i].values.tolist()
		index = ''.join(index)
		i += 1
		add_to_fasta(index, pid)

		# update the node's current sequence
		current_sequences[pid] = index

	else:
		# Get the parent's sequence, mutate it, and append result to fasta
		print('Mutating sequence, adding to fasta.....')
		seq_to_change = current_sequences[contact_pid]
		change = determine_change(thresh)
		print(pid)
		new_seq = commit_change(seq_to_change, change)
		add_to_fasta(new_seq, pid)
		
		current_sequences[pid] = new_seq
