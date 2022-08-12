import networkx as nx
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt

#Read in network data
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
G=nx.from_pandas_edgelist(df, 'pid', 'contact_pid', create_using=nx.DiGraph())

#check if the graph is a directed acyclic graph
print(nx.is_directed_acyclic_graph(G))

#returns cycles
print(sorted(nx.simple_cycles(G)))

#draw graph
nx.draw(G, with_labels=False, node_color = colors_list, cmap='turbo', node_size = 100, font_size = 5, width = .5 )
plt.show()

