#Import libraries
from percolation_basefunctions import *

def targeted_attack_adaptive(multiplex_structure,N,ranking_targeted):
    ranking_nodes = ranking_targeted(multiplex_structure,N)
    ranking_nodes_rank=sorted(ranking_nodes.items(),key= lambda x: x[1],reverse=True)
    ranking_nodes_rank_np=np.array(ranking_nodes_rank)
    ## This is the max element in the ranking
    max_element=ranking_nodes_rank_np[0][1]
    ## Take all the nodes that have the maximum rank
    top_ranking_nodes=ranking_nodes_rank_np[ranking_nodes_rank_np[:,1] == max_element]
    ## Remove one of these nodes at random
    target_ID=random.choice(top_ranking_nodes)[0]
    multiplex_prox = cPickle.loads(cPickle.dumps(multiplex_structure, -1))
    LMCGC,_=compute_LMCGC_block(multiplex_prox,N,[])
    print("Targeted strategy - using the function:",ranking_targeted.__name__)
    print("# of nodes removed, size of LMCGC")
    print(0,LMCGC)
    size_nodes=0
    targeted_steps_attack=[]
    total_target_nodes=[]
    while LMCGC > np.sqrt(N):
        #Remove nodes one by one based on the targeted attack
        target_nodes=[target_ID]
        total_target_nodes.append([target_ID])
        LMCGC,multiplex_prox=compute_LMCGC_block(multiplex_prox,N,target_nodes)
        size_nodes += 1
        targeted_steps_attack.append([size_nodes,LMCGC])
        if multiplex_prox==0:
            break
        ranking_nodes = ranking_targeted(multiplex_prox,N)
        ranking_nodes_rank=sorted(ranking_nodes.items(),key= lambda x: x[1],reverse=True)
        ranking_nodes_rank_np=np.array(ranking_nodes_rank)
        ## This is the max element in the ranking
        max_element=ranking_nodes_rank_np[0][1]
        ## Take all the nodes that have the maximum rank
        top_ranking_nodes=ranking_nodes_rank_np[ranking_nodes_rank_np[:,1] == max_element]
        ## Remove one of these nodes at random
        target_ID=random.choice(top_ranking_nodes)[0]
    print(size_nodes,LMCGC)
    return(total_target_nodes,size_nodes,targeted_steps_attack)



## Pareto Strategy based on the removal of the Pareto point which is closest to the Ideal Point  (3 objectives, 2 single layers and 1 multilayer ranking)
def Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,ranking_singlelayer,ranking_multilayer):
    if ranking_singlelayer.__name__=='BPalgorithm_ranking':
        ranks_single = ranking_singlelayer(multiplex_structure, N,[])
    else:
        ranks_single = ranking_singlelayer(multiplex_structure, N)
    ranks_multi  = ranking_multilayer(multiplex_structure, N)
    width_obj1=(max(list(ranks_single[0].values()))-min(list(ranks_single[0].values())))
    width_obj2=(max(list(ranks_single[1].values()))-min(list(ranks_single[1].values())))
    width_obj3=(max(list(ranks_multi.values()))-min(list(ranks_multi.values())))
    if width_obj1 == 0:
        width_obj1=1
    if width_obj2 == 0:
        width_obj2=1
    if width_obj3 == 0:
        width_obj3=1
    more_objectives=dict({key: [
        ranks_single[0][key]/width_obj1,
        ranks_single[1][key]/width_obj2,
        ranks_multi[key]/width_obj3] for key in ranks_multi})
    pareto_shells=calc_Pareto_front_fast(more_objectives)
    df1 = pd.DataFrame(list(pareto_shells.items()), columns=['nodeID', 'Shell'])
    multiplex_prox = cPickle.loads(cPickle.dumps(multiplex_structure, -1))
    LMCGC,_=compute_LMCGC_block(multiplex_prox,N,[])
    print('''Strategy based on the removal of the Pareto point closer
    to the Ideal Point - Pareto with %s (single) and %s (multi) - 3 objs pareto''' % (ranking_singlelayer.__name__,ranking_multilayer.__name__))
    print("# of nodes removed, size of LMCGC")
    print(0,LMCGC)
    size_nodes=0
    total_target_nodes=[]
    Pareto_targeted_steps_attack=[]
    first_front = list(df1.loc[df1['Shell']==1]['nodeID'])
    while LMCGC > np.sqrt(N):
        first_front = list(df1.loc[df1['Shell']==1]['nodeID'])
        maxvalues = [max(ranks_single[0].values())/width_obj1,
            max(ranks_single[1].values())/width_obj2,
            max(ranks_multi.values())/width_obj3]
        distances_dict={i:euclidean_distance(maxvalues,more_objectives[i]) for i in first_front}
        #Find the minimum of the distances
        minimum=min(distances_dict.values())
        #Find all the nodes that have minimal euclidean distance from the ideal point
        top_ranking_nodes=[[i,distances_dict[i]] for i in distances_dict if distances_dict[i]==minimum]
        #Select one of these nodes at random
        target_ID=random.choice(top_ranking_nodes)[0]
        target_nodes = [target_ID]
        Pareto_targeted_steps_attack.append([size_nodes,LMCGC])
        LMCGC,multiplex_prox=compute_LMCGC_block(multiplex_prox,N,target_nodes)
        size_nodes += len(target_nodes)
        total_target_nodes.append(target_nodes)
        if multiplex_prox==0:
            break
        if ranking_singlelayer.__name__=='BPalgorithm_ranking':
            ranks_single=ranking_singlelayer(multiplex_prox,N,total_target_nodes)
        else:
            ranks_single=ranking_singlelayer(multiplex_prox,N)
        ranks_multi = ranking_multilayer(multiplex_prox,N)
        width_obj1=(max(list(ranks_single[0].values()))-min(list(ranks_single[0].values())))
        width_obj2=(max(list(ranks_single[1].values()))-min(list(ranks_single[1].values())))
        width_obj3=(max(list(ranks_multi.values()))-min(list(ranks_multi.values())))
        if width_obj1 == 0:
            width_obj1=1
        if width_obj2 == 0:
            width_obj2=1
        if width_obj3 == 0:
            width_obj3=1
        more_objectives=dict({key: [
            ranks_single[0][key]/width_obj1,
            ranks_single[1][key]/width_obj2,
            ranks_multi[key]/width_obj3] for key in ranks_multi})
        pareto_shells=calc_Pareto_front_fast(more_objectives)
        df1 = pd.DataFrame(list(pareto_shells.items()), columns=['nodeID', 'Shell'])
        # print(size_nodes,LMCGC)
    print(size_nodes,LMCGC)
    return(total_target_nodes,size_nodes,Pareto_targeted_steps_attack)





## Pareto Strategy based on the removal of the Pareto point which is closest to the Ideal Point (2 objectives -> single layer rankings)
def Pareto_ideal_targeted_attack_2objs(multiplex_structure,N,ranking_singlelayer):
    if ranking_singlelayer.__name__=='BPalgorithm_ranking':
        ranks_single = ranking_singlelayer(multiplex_structure, N,[])
    else:
        ranks_single = ranking_singlelayer(multiplex_structure, N)
    width_obj1=(max(list(ranks_single[0].values()))-min(list(ranks_single[0].values())))
    width_obj2=(max(list(ranks_single[1].values()))-min(list(ranks_single[1].values())))
    if width_obj1 == 0:
        width_obj1=1
    if width_obj2 == 0:
        width_obj2=1
    more_objectives=dict({key: [
        ranks_single[0][key]/width_obj1,
        ranks_single[1][key]/width_obj2] for key in ranks_single[0]})
    pareto_shells=calc_Pareto_front_fast(more_objectives)
    df1 = pd.DataFrame(list(pareto_shells.items()), columns=['nodeID', 'Shell'])
    multiplex_prox = cPickle.loads(cPickle.dumps(multiplex_structure, -1))
    LMCGC,_=compute_LMCGC_block(multiplex_prox,N,[])
    print('''Strategy based on the removal of the Pareto point closer
    to the Ideal Point - Pareto with %s on each layer - 2 objs pareto''' % (ranking_singlelayer.__name__))
    print("# of nodes removed, size of LMCGC")
    print(0,LMCGC)
    size_nodes=0
    total_target_nodes=[]
    Pareto_targeted_steps_attack=[]
    first_front = list(df1.loc[df1['Shell']==1]['nodeID'])
    while LMCGC > np.sqrt(N):
        first_front = list(df1.loc[df1['Shell']==1]['nodeID'])
        maxvalues = [max(ranks_single[0].values())/width_obj1,
            max(ranks_single[1].values())/width_obj2]
        distances_dict={i:euclidean_distance(maxvalues,more_objectives[i]) for i in first_front}
        #Find the minimum of the distances
        minimum=min(distances_dict.values())
        #Find all the nodes that have minimal euclidean distance from the ideal point
        top_ranking_nodes=[[i,distances_dict[i]] for i in distances_dict if distances_dict[i]==minimum]
        #Select one of these nodes at random
        target_ID=random.choice(top_ranking_nodes)[0]
        target_nodes = [target_ID]
        Pareto_targeted_steps_attack.append([size_nodes,LMCGC])
        LMCGC,multiplex_prox=compute_LMCGC_block(multiplex_prox,N,target_nodes)
        size_nodes += len(target_nodes)
        total_target_nodes.append(target_nodes)
        if multiplex_prox==0:
            break
        if ranking_singlelayer.__name__=='BPalgorithm_ranking':
            ranks_single=ranking_singlelayer(multiplex_prox,N,total_target_nodes)
        else:
            ranks_single = ranking_singlelayer(multiplex_prox,N)
        width_obj1=(max(list(ranks_single[0].values()))-min(list(ranks_single[0].values())))
        width_obj2=(max(list(ranks_single[1].values()))-min(list(ranks_single[1].values())))
        if width_obj1 == 0:
            width_obj1=1
        if width_obj2 == 0:
            width_obj2=1
        more_objectives=dict({key: [
            ranks_single[0][key]/width_obj1,
            ranks_single[1][key]/width_obj2] for key in ranks_single[0]})
        pareto_shells=calc_Pareto_front_fast(more_objectives)
        df1 = pd.DataFrame(list(pareto_shells.items()), columns=['nodeID', 'Shell'])
        # print(size_nodes,LMCGC)
        # print(target_nodes,size_nodes,LMCGC)
    print(size_nodes,LMCGC)
    return(total_target_nodes,size_nodes,Pareto_targeted_steps_attack)
