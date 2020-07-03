#Import libraries
import numpy as np
import pandas as pd
import random
import seaborn
import shutil
from functools import reduce
import matplotlib.pyplot as plt
import networkx as nx
import collections
import itertools
from collections import Counter
from operator import itemgetter
import matplotlib.pyplot as plt
import os
import sys
import copy
import time
import subprocess
import heapq
import _pickle as cPickle
from  subprocess import Popen, PIPE,run


# Loading multiplex network from file (e.g. list of layers) into a dict of dictionary
def load_dataset_into_dict_of_dict(file,N):
    """
    :param file: name of the file  containing the list of layers of the multiplex
    :return: multiplex_structure: dictionary of dictionary containing the layers in the format:
    {layer_ID: {key:node; values: first_neighbours}}
    """
    multiplex_structure = dict()
    #take the path of the file
    path = os.path.split(file)[0]
    #open the file and takae all the names of the layers that are in the file
    multiplex_list= pd.read_csv(file, sep= ' ',names=['layer_name'])
    #Iterating over the layers name
    for row in multiplex_list.itertuples():
        layer_to_analyse=path+'/'+row.layer_name    
        # Loading the layer and save it in the dictionary of networks
        # e.g. {layersID:{key:node; values:[first_neighbours}}
        multiplex_structure=load_layers(multiplex_structure,layer_to_analyse,row.Index,N)
    return(multiplex_structure)


##Loading the singular networks as dictionary => {key:node; values:[first_neighbours]}
##set undir = False - > the network is directed
def load_layers(multiplex_structure,layer,index_layer,N,undir = True):

    #Loading the layer (from the edge list representation) and put it 
    layer_struct=np.loadtxt(layer,dtype=int,ndmin=2)
    multiplex_structure[index_layer]={}
    for edges in layer_struct:
        node_i = int(edges[0])
        node_j = int(edges[1])
        #If the node i exists in the multiplex, add node j as node_i neighbour
        if node_i in multiplex_structure[index_layer]:
            multiplex_structure[index_layer][node_i].append(node_j)
        else:
            multiplex_structure[index_layer][node_i]=[node_j]
        if undir == True:
            if node_j in multiplex_structure[index_layer]:
                multiplex_structure[index_layer][node_j].append(node_i)
            else:
                multiplex_structure[index_layer][node_j]=[node_i]
    #Put empty list for isolated nodes
    for i in range(0,N):
        if i not in multiplex_structure[index_layer]:
            multiplex_structure[index_layer][i]=[]
    return(multiplex_structure)

#Function that creates the networkx structure from the edge list dictionary. It returns the graph G as output
def create_graph_fromedgelist(graph_structure):
    """
    :param file: dictionary of neighbours  {key: node; values: first_neighbours}
    :return: Graph as a networkx formatt 
    """
    G=nx.Graph()
    for key in graph_structure.keys():
        for l in graph_structure[key]:
            G.add_edge(key,l)
    return(G)


def compute_connected_comp(multiplex_structure,N):
    """
    Compute the connected component of the multiplex network using the networkx library
    :param file: dictionary of dictionaries containing the multiplex networks, and the total number of nodes of the network N 
    :return: dictionary of the list of connected components on each layer, format: {key: layerID; values: list conn. components} 
    """
    nodes=[]
    [nodes.append(list(multiplex_structure[i].keys())) for i in multiplex_structure]
    unique_nodes=list(set([item for sublist in nodes for item in sublist]))
    list_conn_components= dict()
    conn_components=dict()
    total_nodes=[i for i in range(0,N)]
    for layer_ID in multiplex_structure:
        list_conn_components[layer_ID] = dict.fromkeys(total_nodes,-100000)
        G = create_graph_fromedgelist(multiplex_structure[layer_ID])
        ##Sort the connected components of the graph  by size, and assign an ID to each node according to the size of cluster
        ##Components that are isolated have ID equal to -100000
        ##The component ID starts from 1 
        conn_components[layer_ID]=sorted(nx.connected_components(G), key = len, reverse=True)
        s= len(conn_components[layer_ID])
        for number_component in range(1,s+1):
            for l in conn_components[layer_ID][number_component-1]:
                list_conn_components[layer_ID][l]=number_component
    return(list_conn_components)


#Compute the connected component of a specific "layer_ID" layer and update the dictionary list_conn_components
def compute_connected_comp_layer(multiplex_structure,N,layer_ID,list_conn_components):
    """
    Compute the connected component of a specific "layer_ID" layer and update the dictionary list_conn_components
    :param file: dictionary of dictionaries containing the multiplex networks, the number of nodes N, layer_ID for computing the conn.
    component,  dictionary of the list of connected components on each layer, format: {key: layerID; values: list conn. components} 
    :return: updated dictionary of the list of connected components on each layer, format: {key: layerID; values: list conn. components} 
    """
    nodes=[]
    M=len(multiplex_structure)
    [nodes.append(list(multiplex_structure[i].keys())) for i in multiplex_structure]
    unique_nodes=list(set([item for sublist in nodes for item in sublist]))
    conn_components=dict()
    total_nodes=[i for i in range(0,N)]
    G = create_graph_fromedgelist(multiplex_structure[layer_ID])
    #print("computing components")
    conn_components[layer_ID]=sorted(nx.connected_components(G), key = len, reverse=True)
    s= len(conn_components[layer_ID])
    for number_component in range(1,s+1):
        for l in conn_components[layer_ID][number_component-1]:
            list_conn_components[layer_ID][l]=number_component
    return(list_conn_components)

#Compute the product of the degrees of each layer and return it as a dictionary nodeID:Product degrees 
def computing_product_degree(multiplex_structure,N):
    """
    Compute the product of the degrees of each layer and return it as a dictionary nodeID:Product degrees 
    :param file: dictionary of dictionaries containing the multiplex networks, the number of nodes N
    :return: dictionary containing the product of the degrees of each layers, format: nodeID:Product degrees 
    """
    total_nodes=[i for i in range(0,N)]
    total_degree = dict.fromkeys(total_nodes, 0)
    for layer_ID in multiplex_structure:
        for nodes in multiplex_structure[layer_ID]:
            if total_degree[nodes] == 0:
                total_degree[nodes] = len(multiplex_structure[layer_ID][nodes])
            else:
                total_degree[nodes]*=len(multiplex_structure[layer_ID][nodes])
    return(total_degree)


#### Slow implementation for computing the LMCGC of a multiplex networks ####
#### For duplex up to ~ 10^4 nodes is still fine                         #### 
#Compute the largest cluster of the Mutually Connected Giant Component
#Multiplex structure: dict of dict
#N: number of nodes
#Block_target: list of nodes to be removed from the multiplex
def compute_LMCGC_block(multiplex_structure,N,block_target):
    """
    Compute the largest cluster of the Mutually Connected Giant Component
    :param file: dictionary of dictionaries containing the multiplex networks, the number of nodes N, the list of nodes to be removed from the multiplex
    :return: the size of the largest cluster of the Mutually Connected Giant Component, the percolated multiplex 
    """
    multiplex_new=cPickle.loads(cPickle.dumps(multiplex_structure, -1))
    #Removing the connections of the target in all the layers:
    #print(block_target)
    for target in block_target:
        for layer_ID in multiplex_structure:
            if target in multiplex_structure[layer_ID]:
                for j in multiplex_structure[layer_ID][target]:
                    if j in multiplex_new[layer_ID][target]:
                        multiplex_new[layer_ID][target].remove(j)
                    if target in multiplex_new[layer_ID][j]:
                        multiplex_new[layer_ID][j].remove(target)
    multiplex_next=cPickle.loads(cPickle.dumps(multiplex_new, -1))

    #Now computing the MCGC in the multiplex network [using the algorithm of V.Buldyrev 2010 - Nature -> Not optimised]
    while(1):
        flag = 0
        #multiplex_new = copy.deepcopy(multiplex_next)
        multiplex_new = cPickle.loads(cPickle.dumps(multiplex_next, -1))
        #Computing the connected components of all the layers in the multiplex
        list_connected_components=compute_connected_comp(multiplex_new,N)                   
        for layer_ID in multiplex_new:
            for i in multiplex_new[layer_ID].keys():
                for j in multiplex_new[layer_ID][i]:
                    if i < j:
                        for layer_check in range(0,len(multiplex_new)):
                            #remove the links that don't belong to the same component in the different layers
                            if list_connected_components[layer_check][i] != list_connected_components[layer_check][j] or list_connected_components[layer_check][i] == -100000 or list_connected_components[layer_check][j] == -100000:
                                    multiplex_next[layer_ID][i].remove(j)
                                    multiplex_next[layer_ID][j].remove(i)
                                    flag = 1
            #Computing the connected components of all the layers in the multiplex
            list_connected_components=compute_connected_comp_layer(multiplex_new,N,layer_ID,list_connected_components)
        #If no changes happened in a cycle, then exit!
        if flag == 0:
            break
    G = create_graph_fromedgelist(multiplex_next[0])
    #Check the size of the connected component of one of the layers, after all the link removals, the giant component represents the MCC
    conn_components=sorted(nx.connected_components(G), key = len, reverse=True)
    if not conn_components:
        return(0,0)
    else:
        return(len(conn_components[0]),multiplex_next)

        

# Function computing the total degree of each node of the multiplex structure, in addition
# it also computes the  participation coefficient [See Nicosia-Battiston-Latora 2014 - Structural measures
# for multiplex networks for all the details]
def computing_total_degree_and_participation(multiplex_structure,N):
    """
    :param multiplex_structure: structure containing the multiplex (dictionary of dictionary) 
        as from the function load_dataset_into_dict_of_dict
            N : number of total nodes in the multiplex
    :return:
            total_degree: dictionary {node: total_degree}
            participation_degree: dictionary {node: participation coefficient (of the degree)}
            dgr_part: dictionary: {node: [total_degree[node],participation_coefficient[node]]}
    """
    nodes=[]
    #multiplex_structure_aux=copy.deepcopy(multiplex_structure)
    multiplex_structure_aux = cPickle.loads(cPickle.dumps(multiplex_structure, -1))
    M=len(multiplex_structure)
    [nodes.append(list(multiplex_structure_aux[i].keys())) for i in multiplex_structure_aux]
    #Extract the nodes that are active in the multiplex network
    unique_nodes=list(set([item for sublist in nodes for item in sublist]))
    total_nodes=[i for i in range(0,N)]
    total_degree = dict.fromkeys(total_nodes, 0)
    degree_layers_multiplex = {}
    for layer_ID in multiplex_structure_aux:
        degree_layers_multiplex[layer_ID]=dict.fromkeys(total_nodes,0)
        for nodes in multiplex_structure_aux[layer_ID]:
            total_degree[nodes]+=len(multiplex_structure_aux[layer_ID][nodes])
    #Total degree is a dictionary containing the sum of the degree of each nodes in all the layers
                
    participation_degree = dict.fromkeys(total_nodes, 0)
    for node in unique_nodes:
        temp_sum=0.0
        for layer_ID in multiplex_structure_aux:
            if node in multiplex_structure_aux[layer_ID] and total_degree[node] != 0:
                degree_node= len(multiplex_structure_aux[layer_ID][node])
                degree_layers_multiplex[layer_ID][node] = degree_node
                temp_sum+=((1.0*degree_node)/(1.0*total_degree[node]))**2
        participation_degree[node]=(M/(M-1.0))*(1-temp_sum)
        if node in multiplex_structure_aux[layer_ID] and total_degree[node] == 0:
            participation_degree[node]=0

    dgr_part=dict({key: [total_degree[key], participation_degree[key]] for key in participation_degree})
    return(total_degree,participation_degree,dgr_part,degree_layers_multiplex)


#Compute the intersection graph of the multiplex structure
def find_intersection_graph(multiplex_structure,N):
    intersection_graph={}
    for i in range(N):
        if i in multiplex_structure[0]:
            all_neighbour=multiplex_structure[0][i]
        else:
            continue
        #find all the neighbours of node i in all the layers
        for k in multiplex_structure:
            all_neighbour=set(all_neighbour).intersection(multiplex_structure[k][i])
        intersection_graph[i]=list(all_neighbour)
    return(intersection_graph)

#Compute the union graph of the multiplex structure
def find_union_graph(multiplex_structure,N):
    union_graph={}
    for i in range(N):
        if i in multiplex_structure[0]:
            if multiplex_structure[0][i]!= []:
                all_neighbour=multiplex_structure[0][i]
            else:
                all_neighbour=[]
        else:
            continue
        #find all the neighbours of node i in all the layers
        for k in multiplex_structure:
            all_neighbour=set(all_neighbour).union(multiplex_structure[k][i])
        union_graph[i]=list(all_neighbour)
    return(union_graph)


#Compute the disjunt union graph of the multiplex structure
def find_symmetric_difference_graph(multiplex_structure,N):
    symmetric_graph={}
    for i in range(N):
        if i in multiplex_structure[0]:
            all_neighbour=multiplex_structure[0][i]
        else:
            continue
        for k in multiplex_structure:
            all_neighbour=set(all_neighbour).symmetric_difference(multiplex_structure[k][i])
        symmetric_graph[i]=list(all_neighbour)
    return(symmetric_graph)


# Compute the EMD 1step ranking as in Eq. 4 of the paper "Targeted damage to interdependent networks" 
# by G. J. Baxter, G. Timar, and J. F. F. Mendes, PRE 98, 032307 (2018) 
def EMD1step_ranking(multiplex_structure, N):
    total_nodes=[i for i in range(0,N)]
    emd_ranking=dict.fromkeys(total_nodes,0)
    total_degree,participation_degree,dgr_part,degree_layers_multiplex = computing_total_degree_and_participation(multiplex_structure,N)
    for layerID in multiplex_structure:
        for nodei in multiplex_structure[layerID]:
            temp = 0
            for neigh_of_i in multiplex_structure[layerID][nodei]:
                if degree_layers_multiplex[layerID][neigh_of_i] !=0:
                    activity = np.heaviside(degree_layers_multiplex[layerID][neigh_of_i],0)+np.heaviside(degree_layers_multiplex[(layerID+1)%2][neigh_of_i],0)
                    temp+=((1./activity)*(total_degree[neigh_of_i]/degree_layers_multiplex[layerID][neigh_of_i]))
            emd_ranking[nodei]+=temp
    return(emd_ranking)




# Compute the EMD ranking as in Eq. 1 of the paper "Targeted damage to interdependent networks" 
# by G. J. Baxter, G. Timar, and J. F. F. Mendes, PRE 98, 032307 (2018) 
# The criterion for the solutions to have fully converged is that the largest
# relative difference in all values is less than 10e-7
def EMDfullsteps_ranking(multiplex_structure, N):
    nodes=[]
    M=len(multiplex_structure)
    [nodes.append(list(multiplex_structure[i].keys())) for i in multiplex_structure]
    unique_nodes=list(set([item for sublist in nodes for item in sublist]))
    total_nodes=[i for i in range(0,N)]
    emd_ranking=dict.fromkeys(total_nodes,0)
    total_degree,participation_degree,dgr_part,degree_layers_multiplex = computing_total_degree_and_participation(multiplex_structure,N)
    for layerID in multiplex_structure:
        for nodei in multiplex_structure[layerID]:
            temp = 0
            for neigh_of_i in multiplex_structure[layerID][nodei]:
                activity = np.heaviside(degree_layers_multiplex[layerID][neigh_of_i],0)+np.heaviside(degree_layers_multiplex[(layerID+1)%2][neigh_of_i],0)
                temp+=((1./activity)*(total_degree[neigh_of_i]/degree_layers_multiplex[layerID][neigh_of_i]))
            emd_ranking[nodei]+=temp
    while(1):
        emd_ranking_new=dict.fromkeys(total_nodes,0)
        for layerID in multiplex_structure:
            for nodei in multiplex_structure[layerID]:
                temp = 0
                for neigh_of_i in multiplex_structure[layerID][nodei]:
                    activity = np.heaviside(degree_layers_multiplex[layerID][neigh_of_i],0)+np.heaviside(degree_layers_multiplex[(layerID+1)%2][neigh_of_i],0)
                    temp+=((1./activity)*(emd_ranking[neigh_of_i]/degree_layers_multiplex[layerID][neigh_of_i]))
                emd_ranking_new[nodei]+=temp
        ranking=np.array(list(emd_ranking.values()))
        ranking_new=np.array(list(emd_ranking_new.values()))
        #The elements are ordered in the array 
        diff_toll=np.abs(ranking_new-ranking)
        max_difference=max([diff_toll[i]/ranking[i] for i in range(len(diff_toll)) if ranking[i] !=0])
        if  max_difference < 1e-7:
            emd_ranking=cPickle.loads(cPickle.dumps(emd_ranking_new, -1))
            #emd_ranking =copy.deepcopy(emd_ranking_new)
            break
        #emd_ranking =copy.deepcopy(emd_ranking_new)
        emd_ranking =cPickle.loads(cPickle.dumps(emd_ranking_new, -1))
    return(emd_ranking)


### Compute the betwenneess centrality of a specific layer "layerID" in the multiplex
def compute_betweenness_centrality(multiplex_structure,layerID,N):
    G = create_graph_fromedgelist(multiplex_structure[layerID])
    list_centrality=nx.betweenness_centrality(G)
    for i in range(0,N):
        if i not in list_centrality:
            list_centrality[i]=0
    return(list_centrality)

### K-core ranking of the nodes in the multiplex 
def kcore_ranking(multiplex_structure,N):
    total_nodes=[i for i in range(N)]
    G={}
    G[0]=create_graph_fromedgelist(multiplex_structure[0])
    G[1]=create_graph_fromedgelist(multiplex_structure[1])
    core_numbers={}
    core_numbers[0]=nx.core_number(G[0])
    core_numbers[1]=nx.core_number(G[1])
    max_cores=[max(core_numbers[0].values()),max(core_numbers[1].values())]
    if max(max_cores)>=1:
        ranking_kcore={}
        for layerID in multiplex_structure:
            ranking_kcore[layerID]={}
            for nodei in total_nodes:
                if nodei in core_numbers[layerID]:
                    #if core_numbers[layerID][nodei] == max_cores[layerID] and core_numbers[layerID][nodei] >=2:
                    ranking_kcore[layerID][nodei]=core_numbers[layerID][nodei]
                #    else:
                #       ranking_kcore[layerID][nodei]=0
                else:
                    ranking_kcore[layerID][nodei]=0
    else:
        ranking_kcore={}
        for layerID in multiplex_structure:
            ranking_kcore[layerID]={}
            for nodei in total_nodes:
                ranking_kcore[layerID][nodei]=0
    return(ranking_kcore)


def calc_Pareto_front_fast(objectives,return_mask=True):
    """
    Compute the pareto front of the input (by default, it will maximise each objective function)
    :param objectives: dictionary containing the M-dimensional cost functions. Format  of the dict {nodeID: [f_1, f_2,...,f_M]}
    :return: dictionary containing the Pareto optimality check for each node, i.e. {nodeID: 1}
    """
    a=sorted(objectives.items(), key=lambda x: x[0],reverse=False)
    total_nodes=[i[0] for i in list(a)]
    fronts_dict = dict.fromkeys(total_nodes, 100000)
    costs= np.array([a[i][1] for i in range(0,len(a))])
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index<len(costs):
        nondominated_point_mask = np.any(costs>costs[next_point_index], axis=1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index])+1
    if return_mask:
        is_efficient_mask = np.zeros(n_points, dtype = bool)
        is_efficient_mask[is_efficient] = True
        for val in enumerate(is_efficient_mask):
            if val[1] == True:
                fronts_dict[total_nodes[val[0]]]=1
        return fronts_dict
    else:
        return is_efficient


## Find the 2 core given a networkx graph G ->  return a dictionary containing {nodeID:degree_2core}
def extract_2coreallfromgraph(G,N):
    total_nodes=[i for i in range(N)]
    two_core=nx.k_core(G,k=2)
    two_core_degree=two_core.degree()
    two_core_sorted=sorted(two_core_degree,key=lambda k: k[1],reverse=True)
    ranking_nodes=dict.fromkeys(total_nodes,0)
    for element in two_core_sorted:
        node_i=element[0]
        degree_i=element[1]
        ranking_nodes[node_i]=degree_i
    return(ranking_nodes) 


### Compute the CoreHD ranking as proposed in the paper by
### Zdeborova L, Zhang P, Zhou HJ. "Fast and simple decycling and dismantling of networks."
### Sci Rep. 2016;6:37954 (2016). doi:10.1038/srep37954
def coreHD_ranking(multiplex_structure,N):
    total_nodes=[i for i in range(N)]
    G={}
    G[0]=create_graph_fromedgelist(multiplex_structure[0])
    G[1]=create_graph_fromedgelist(multiplex_structure[1])
    core_numbers={}
    core_numbers[0]=extract_2coreallfromgraph(G[0],N)
    core_numbers[1]=extract_2coreallfromgraph(G[1],N)
    ranking_coreHD={}
    max_core={}
    max_core[0]=max(core_numbers[0].values())
    max_core[1]=max(core_numbers[1].values())
    for layerID in multiplex_structure:
        if max_core[layerID]!=0:
            ranking_coreHD[layerID]=core_numbers[layerID]
        else:
            # If it enters in this case-> The graph is just a forest, the tree breaking procedure used is the one developed in
            # Braunstein, Alfredo, Luca Dall'Asta, Guilhem Semerjian, and Lenka Zdeborova. "Network Dismantling."
            # Proceedings of the National Academy of Sciences, October 18, 2016, 201605083. doi:10.1073/pnas.1605083113
            ranking_coreHD[layerID]=tree_breaker(G[layerID],N)
            ## Also the betweenness centrality is a good proxy for tree breaking but much slower than the method above
            # ranking_coreHD[layerID]=compute_betweenness_centrality(multiplex_structure,layerID,N)
    return(ranking_coreHD) 


## Compute the Euclidean distance between two points given in input
def euclidean_distance(a,b):
    n = len(a)
    distance = 0
    for i in range(n):
        distance += (a[i]-b[i])**2
    return (np.sqrt(distance))



#### Ranking based on the generalisation of CI in the case of duplex network:
###\sum_\alpha \sum_{j in N_i^\alpha} \frac{  k_j^\alpha*k_j^\beta - k_j^{intersect} }{  k_j^{\alpha} }
def compute_ranking_CI_duplex(multiplex_structure,N):
    total_nodes=[i for i in range(N)]
    fitness=dict.fromkeys(total_nodes,0)
    total_degree,participation_degree,dgr_part,degree_layers_multiplex = computing_total_degree_and_participation(multiplex_structure,N)
    aggregate_network=find_union_graph(multiplex_structure,N)
    intersection_graph=find_intersection_graph(multiplex_structure,N)
    ranking={}
    for layerID in multiplex_structure:
        for nodei in multiplex_structure[layerID]:
            temp=0
            for neigh_of_i in multiplex_structure[layerID][nodei]:
                if degree_layers_multiplex[layerID][neigh_of_i] != 0:
                    numerator1=(degree_layers_multiplex[(layerID+1)%2][neigh_of_i])*(degree_layers_multiplex[layerID][neigh_of_i])
                    numerator2=len(intersection_graph[neigh_of_i])
                    numerator=numerator1-numerator2
                    denominator=degree_layers_multiplex[layerID][neigh_of_i]
                    temp+= (numerator/denominator)
            fitness[nodei]+=temp
    for nodei in fitness:
        if len(aggregate_network[nodei])!=0:
            numerator1=(degree_layers_multiplex[0][nodei])*(degree_layers_multiplex[1][nodei])
            numerator2=len(intersection_graph[nodei])
            numerator=numerator1-numerator2
            denominator=len(aggregate_network[nodei])
            fitness[nodei]=fitness[nodei] * (numerator/denominator)
    return(fitness)

#### Ranking based on the generalisation of CI in the case of duplex network
def compute_ranking_DCI(multiplex_structure,N):
    total_nodes=[i for i in range(N)]
    fitness=dict.fromkeys(total_nodes,0)
    total_degree,participation_degree,dgr_part,degree_layers_multiplex = computing_total_degree_and_participation(multiplex_structure,N)
    aggregate_network=find_union_graph(multiplex_structure,N)
    intersection_graph=find_intersection_graph(multiplex_structure,N)
    ranking={}
    for layerID in multiplex_structure:
        for nodei in multiplex_structure[layerID]:
            temp=0
            for neigh_of_i in multiplex_structure[layerID][nodei]:
                if degree_layers_multiplex[layerID][neigh_of_i] != 0:
                    fitness[nodei] += degree_layers_multiplex[(layerID+1)%2][neigh_of_i] - 1
                    ###fitness[nodei]+=temp
    for nodei in fitness:
        if len(aggregate_network[nodei])!=0:
            numerator1=(degree_layers_multiplex[0][nodei])*(degree_layers_multiplex[1][nodei])
            numerator2=len(intersection_graph[nodei])
            numerator=numerator1-numerator2
            denominator=len(aggregate_network[nodei])
            fitness[nodei]=fitness[nodei] * (numerator/denominator)
    return(fitness)


#### Ranking based on the generalisation of CI in the case of duplex network (with the modification for isolated nodes)
def compute_ranking_DCIz(multiplex_structure,N):
    total_nodes=[i for i in range(N)]
    fitness=dict.fromkeys(total_nodes,0)
    total_degree,participation_degree,dgr_part,degree_layers_multiplex = computing_total_degree_and_participation(multiplex_structure,N)
    aggregate_network=find_union_graph(multiplex_structure,N)
    intersection_graph=find_intersection_graph(multiplex_structure,N)
    ranking={}
    for layerID in multiplex_structure:
        for nodei in multiplex_structure[layerID]:
            temp=0
            for neigh_of_i in multiplex_structure[layerID][nodei]:
                if degree_layers_multiplex[layerID][neigh_of_i] != 0:
                    fitness[nodei] += degree_layers_multiplex[(layerID+1)%2][neigh_of_i] - 1
                    ###fitness[nodei]+=temp
    for nodei in fitness:
        if len(aggregate_network[nodei])!=0:
            numerator1=(degree_layers_multiplex[0][nodei]+1)*(degree_layers_multiplex[1][nodei]+1)
            numerator2=3*len(intersection_graph[nodei])
            numerator=numerator1-numerator2-1
            denominator=len(aggregate_network[nodei])
            fitness[nodei]=fitness[nodei] * (numerator/denominator)
    return(fitness)


#### Slow implementation of the Collective Influence algorithm with l=2 (The original code is much faster using a
#### max-heap data structure)
#### presented in the paper: Morone, F., Makse, H. Influence maximization in complex networks
#### through optimal percolation. Nature 524, 65-68 (2015). https://doi.org/10.1038/nature14604
def compute_CI_2_ranking(multiplex_structure,N):
    total_nodes=[i for i in range(N)]
    fitness={}
    total_degree,participation_degree,dgr_part,degree_layers_multiplex = computing_total_degree_and_participation(multiplex_structure,N)
    for layerID in multiplex_structure:
        fitness[layerID]=dict.fromkeys(total_nodes,0)
        for nodei in multiplex_structure[layerID]:
            temp=0
            for neigh_of_i in multiplex_structure[layerID][nodei]:
                for second_neigh in multiplex_structure[layerID][neigh_of_i]:
                    temp+=(degree_layers_multiplex[layerID][second_neigh]-1)
                       
            fitness[layerID][nodei]=temp*(degree_layers_multiplex[layerID][nodei]-1)
    return(fitness)

#### Compute the ranking of the nodes in the multiplex using 
#### the  Believe Propagation algorithm  as presented in
#### "Identifying optimal targets of network attack by belief propagation", Physical Review E  94, 012305 (2016)
#### External C++ code, gently provided by Hai-Jun Zhou - http://lib.itp.ac.cn/html/zhouhj/codes.html
def BPalgorithm_ranking(multiplex_structure,N,target_nodes,jobID=1,iteration=1):
    total_nodes=[i for i in range(N)]
    G={}
    max_nodeID=[0,0]
    numberedges=[0,0]
    list_nodes=[i for i in range(N) if i not in target_nodes]

    total_nodes_map={list_nodes[i]:i+1 for i in range(len(list_nodes))}
    total_nodes_inverse={i+1:list_nodes[i] for i in range(len(list_nodes))}

    for layerID in multiplex_structure:
        filename_BP='network_BP_{0}_job{1}_iteration{2}.txt_new'.format(layerID,jobID,iteration)
        f = open(filename_BP,'w+')
        for nodei in sorted(list(multiplex_structure[layerID].keys())):
            for nodej in multiplex_structure[layerID][nodei]:
                if nodei < nodej:
                    if nodej+1 > max_nodeID[layerID]:
                        max_nodeID[layerID]=nodej+1
                    numberedges[layerID]+=1
                    f.write('%d %d\n' % (total_nodes_map[nodei],total_nodes_map[nodej]))
        f.close()
    os.system("./BP_algorithm network_BP_0_job{0}_iteration{1}.txt_new {2} {3} {4} {5} > network_BP_0_job{6}_iteration{7}.txt_results_new".format(jobID,iteration,max_nodeID[0],numberedges[0],jobID,iteration,jobID,iteration))
    # time.sleep(5.0)
    os.system("./BP_algorithm network_BP_1_job{0}_iteration{1}.txt_new {2} {3} {4} {5} > network_BP_1_job{6}_iteration{7}.txt_results_new".format(jobID,iteration,max_nodeID[1],numberedges[1],jobID,iteration,jobID,iteration))
    # time.sleep(5.0)
    res=collections.OrderedDict()

    with open('network_BP_0_job{0}_iteration{1}.txt_results_new'.format(jobID,iteration)) as f:
        res[0]={total_nodes_inverse[int(k)]: float(v) for line in f for (k, v) in (line.strip().split(' '),)}
    with open('network_BP_1_job{0}_iteration{1}.txt_results_new'.format(jobID,iteration)) as f:
        res[1]={total_nodes_inverse[int(k)]: float(v) for line in f for (k, v) in (line.strip().split(' '),)}
    BP_ranking={}
    for layerID in multiplex_structure:
        BP_ranking[layerID]={}
        flag =0;
        for i in res[layerID]:
            if res[layerID][i] == 100000 and flag ==0:
                flag =1
                BP_ranking[layerID][i]=1
            if res[layerID][i] != 100000:
                BP_ranking[layerID][i]=res[layerID][i]

    for layerID in multiplex_structure:
        for i in total_nodes:
            if i not in BP_ranking[layerID]:
                BP_ranking[layerID][i]=0

    return(BP_ranking)



######################################################################################################
## Modified version of the tree-breaking implementation by Alfredo Braunstein in the
## paper Braunstein, A., Dall'Asta, L., Semerjian, G., Zdeborova, L., 2016.
## Network dismantling. PNAS 201605083. doi:10.1073/pnas.1605083113, arxiv:1603.08883

class Graph:
    def __init__ (self):
        self.V = []
        self.present = []
        self.M = 0
    def size(self):
        return sum(self.present)

    def add_node(self, i):
        if i >= len(self.V):
            delta = i + 1 - len(self.V)
            self.present += [ 1 ] * delta
            self.V += [ [] for j in range(delta) ]
    def add_edge(self, i, j):
        self.add_node(i)
        self.add_node(j)
        self.V[i] += [j]
        self.V[j] += [i]
        self.M += 1
    def remove_node(self, i):
        self.present[i] = 0;
        self.M -= sum(1 for j in self.V[i] if self.present[j])

def size(G, S, i, j):
    if not G.present[i]:
        return 0
    # if S[i]:
    #     print("# the graph is NOT acyclic")
    #     exit()
    S[i] = 1 + sum(size(G, S, k, i) for k in G.V[i] if k != j and G.present[k])
    return S[i]


def tree_breaker(G_layer,N):
    total_nodes=[i for i in range(N)]
    ranking_nodes=dict.fromkeys(total_nodes,0)
    G = Graph()
    n = 0
    for v in G_layer.edges():
        G.add_edge(int(v[0]),int(v[1]))
    N_nodes = G.size()
    S = [0] * len(G.V)
    H = [(-size(G, S, i, None), i) for i in range(len(G.V)) if G.present[i] and not S[i]]
    Ncc = len(H)
    heapq.heapify(H)
    while len(H):
        s,i = heapq.heappop(H)
        scomp = -s;
        sender = None;
        while True:
            sizes = [ (S[k],k) for k in G.V[i] if k != sender and G.present[k]]
            if len(sizes) == 0:
                break
            M, largest = max(sizes)
            if M <= scomp/2:
                for k in G.V[i]:
                    if S[k] > 1 and G.present[k]:
                        heapq.heappush(H, (-S[k], k))
                G.remove_node(i)
                n+=1
                #print("S", i, n, scomp)
                #Assign the ranking to the nodes. Higher values are better (this is done
                #so that in the  Pareto fronts they are selected first)
                ranking_nodes[i]=N-n
                if scomp <= np.sqrt(N):
                    return(ranking_nodes)
                break;
            S[i] = 1 + sum(S[k] for k in G.V[i] if k != largest and G.present[k])
            sender, i = i, largest
    #Return the ranking of the nodes according to the tree breaking greedy procedure
    #Notice: Higher values are better!!!
    return(ranking_nodes)

#####################################################################################################
def compute_MCCfast_fromexternalsource(multiplex_structure,N,block_target):
    multiplex_new=cPickle.loads(cPickle.dumps(multiplex_structure, -1))
    #Removing the connections of the target in all the layers:
    #print(block_target)
    for target in block_target:
        for layer_ID in multiplex_structure:
            if target in multiplex_structure[layer_ID]:
                for j in multiplex_structure[layer_ID][target]:
                    if j in multiplex_new[layer_ID][target]:
                        multiplex_new[layer_ID][target].remove(j)
                    if target in multiplex_new[layer_ID][j]:
                        multiplex_new[layer_ID][j].remove(target)
    #Convert the duplex in the two edge list
    layers={}
    number_links={}
    for layerID in multiplex_new:
        layers[layerID]=[]
        number_links[layerID]=0
        for nodei in multiplex_new[layerID]:
            for nodej in multiplex_new[layerID][nodei]:
                if nodei < nodej:
                    layers[layerID].append([nodei,nodej])
                    number_links[layerID]+=1
    NodeID_max=N
    K_mul= [number_links[x] for x in number_links]
    string1='\n'.join([str(e1)+' '+ str(e2) for e1,e2 in layers[0]])
    string1=string1+'\n#\n'
    string2='\n'.join([str(e1)+' '+ str(e2) for e1,e2 in layers[1]])
    string2=string2+'\n#\n'
    string=string1+string2
    proc = Popen(["./DecrementalGMCC"]+[str(NodeID_max),str(K_mul[0]),str(K_mul[1]),'0'], stdin=PIPE,stdout=PIPE)
    value = bytes(string, 'UTF-8')  # Needed in Python 3.
    proc.stdin.write(value)
    result =int(proc.communicate()[0])
    #proc.terminate()
    LMCC=result
    return(LMCC,multiplex_new)
