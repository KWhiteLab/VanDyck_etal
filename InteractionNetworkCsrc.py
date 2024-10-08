### The required libraries and packages ###
import networkx as nx
import requests
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm

doc= '/Users/papakobinavandyck/Desktop/PROPKA/CSRC/src_Net.csv'

df=pd.read_csv(doc)

print(df)

interactions = df[['from', 'to', 'weight']]  

carac=pd.DataFrame({'ID':
['H201','H492','C185','E524','D365','D493','D518','E159','E178','E489','E517','K195','K200','K203','R155','R156','R175'],
	'type':['shifts','shifts','shifts','shifts','interactors','interactors','interactors','interactors','interactors','interactors','interactors','interactors','interactors','interactors','interactors','interactors','interactors']})


G=nx.Graph(name='Residue Interaction Graph')
interactions = np.array(interactions)
for i in range(len(interactions)):
    interaction = interactions[i]
    a = interaction[0] # residue a node
    b = interaction[1] # residue b node
    w = float(interaction[2]) # score as weighted edge where high scores = low weight
    G.add_weighted_edges_from([(a,b,w)]) # add weighted edge to graph

carac = carac.set_index('ID')
carac = carac.reindex(G.nodes())
 
carac['type'] = pd.Categorical(carac['type'])
carac['type'].cat.codes

cmap = matplotlib.colors.ListedColormap(['#A35AA2','#8BC7CB'])

pos = nx.spring_layout(G) # position the nodes using the spring layout
plt.figure(figsize=(15,8))
nx.draw(G, with_labels=True, node_color=carac['type'].cat.codes, cmap=cmap, 
        node_size=5000, font_size=20, font_weight="bold", width=0.75, 
        edgecolors='gray')
plt.axis('off')
plt.savefig('cSrc',dpi=600, transparent=True)
plt.show()


# function to rescale list of values to range [newmin,newmax]
def rescale(l,newmin,newmax):
    arr = list(l)
    return [(x-min(arr))/(max(arr)-min(arr))*(newmax-newmin)+newmin for x in arr]
# use the matplotlib plasma colormap
graph_colormap = cm.get_cmap('plasma', 12)
# node color varies with Degree
c = rescale([G.degree(v) for v in G],0.0,0.9) 
c = [graph_colormap(i) for i in c]
# node size varies with betweeness centrality - map to range [10,100] 
bc = nx.betweenness_centrality(G) # betweeness centrality
s =  rescale([v for v in bc.values()],1500,7000)
# edge width shows 1-weight to convert cost back to strength of interaction 
ew = rescale([float(G[u][v]['weight']) for u,v in G.edges],0.1,10)
# edge color also shows weight
ec = rescale([float(G[u][v]['weight']) for u,v in G.edges],0.1,1)
ec = [graph_colormap(i) for i in ec]


pos = nx.spring_layout(G)
plt.figure(figsize=(15,8),facecolor=[0.7,0.7,0.7,0.4])
nx.draw_networkx(G, pos=pos, with_labels=True, node_color=c, node_size=s,edge_color= ec,width=ew,
                 font_color='grey',font_weight='bold',font_size='20')
plt.axis('off')
plt.savefig('Csrc',dpi=600, transparent=True)
plt.show()

