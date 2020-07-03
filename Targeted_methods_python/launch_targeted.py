import sys
from percolation_basefunctions import *
from optimal_percolation_multiplex import *

if len(sys.argv) < 3:
    print("Usage: %s <Number of nodes> <layer_list>" % sys.argv[0])
    sys.exit(1)

N=int(sys.argv[1])
multiplex_filename=sys.argv[2]
multiplex_structure=load_dataset_into_dict_of_dict(multiplex_filename,N)

#################    Multiplex targeted strategies (without Pareto)   #################

HDA_total_target_nodes,HDA_size_nodes,HDA_steps_attack=targeted_attack_adaptive(multiplex_structure,N,computing_product_degree)
print()
EMD1step_total_target_nodes,EMD1step_size_nodes,EMD1step_steps_attack=targeted_attack_adaptive(multiplex_structure,N,EMD1step_ranking)
print()
EMD_total_target_nodes,EMD_size_nodes,EMD_steps_attack=targeted_attack_adaptive(multiplex_structure,N,EMDfullsteps_ranking)
print()
DCI_total_target_nodes,DCI_size_nodes, DCI_steps_attack=targeted_attack_adaptive(multiplex_structure,N,compute_ranking_DCI)
print()
DCIz_total_target_nodes,DCIz_size_nodes, DCIz_steps_attack=targeted_attack_adaptive(multiplex_structure,N,compute_ranking_DCIz)
print()

#################    Pareto Methods  with 3 objs     #################

#k-core with the others

kcore_EMD1step_ideal_total_target_nodes,kcore_EMD1step_ideal_size_nodes,kcore_EMD1step_ideal_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,kcore_ranking,EMD1step_ranking)
print()
kcore_EMD_ideal_total_target_nodes,kcore_EMD_ideal_size_nodes,kcore_EMD_ideal_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,kcore_ranking,EMDfullsteps_ranking)
print()
kcore_DCI_total_target_nodes,kcore_DCI_size_nodes, kcore_DCI_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,kcore_ranking,compute_ranking_DCI)
print()
kcore_DCIz_total_target_nodes,kcore_DCIz_size_nodes, kcore_DCIz_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,kcore_ranking,compute_ranking_DCIz)
print()

#CoreHD with the others

CoreHD_EMD1step_total_target_nodes,CoreHD_EMD1step_size_nodes,CoreHD_EMD1step_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,coreHD_ranking,EMD1step_ranking)
print()
CoreHD_EMD_total_target_nodes,CoreHD_EMD_size_nodes,CoreHD_EMD_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,coreHD_ranking,EMDfullsteps_ranking)
print()
CoreHD_DCI_total_target_nodes,CoreHD_DCI_size_nodes, CoreHD_DCI_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,coreHD_ranking,compute_ranking_DCI)
print()
CoreHD_DCIz_total_target_nodes,CoreHD_DCIz_size_nodes, CoreHD_DCIz_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,coreHD_ranking,compute_ranking_DCIz)
print()

#CI l=2 with the others

CI_2_EMD1step_total_target_nodes,CI_2_EMD1step_size_nodes,CI_2_EMD1step_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,compute_CI_2_ranking,EMD1step_ranking)
print()
CI_2_EMD_total_target_nodes,CI_2_EMD_size_nodes,CI_2_EMD_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,compute_CI_2_ranking,EMDfullsteps_ranking)
print()
CI_2_DCI_total_target_nodes,CI_2_DCI_size_nodes, CI_2_DCI_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,compute_CI_2_ranking,compute_ranking_DCI)
print()
CI_2_DCIz_total_target_nodes,CI_2_DCIz_size_nodes, CI_2_DCIz_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,compute_CI_2_ranking,compute_ranking_DCIz)
print()

#BPD with the others

BPD_EMD1step_total_target_nodes,BPD_EMD1step_size_nodes,BPD_EMD1step_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,BPalgorithm_ranking,EMD1step_ranking)
print()
BPD_EMD_total_target_nodes,BPD_EMD_size_nodes,BPD_EMD_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,BPalgorithm_ranking,EMDfullsteps_ranking)
print()
BPD_DCI_total_target_nodes,BPD_DCI_size_nodes, BPD_DCI_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,BPalgorithm_ranking,compute_ranking_DCI)
print()
BPD_DCIz_total_target_nodes,BPD_DCIz_size_nodes, BPD_DCIz_steps_attack=Pareto_ideal_targeted_attack_3objs(multiplex_structure,N,BPalgorithm_ranking,compute_ranking_DCIz)
print()



################    Pareto Methods  with 2 objs     #################

Pure_CoreHD_total_target_nodes,Pure_CoreHD_size_nodes,Pure_CoreHD_steps_attack=Pareto_ideal_targeted_attack_2objs(multiplex_structure,N,coreHD_ranking)
print()
Pure_CI_2_total_target_nodes,Pure_CI_2_size_nodes,Pure_CI_2_steps_attack=Pareto_ideal_targeted_attack_2objs(multiplex_structure,N,compute_CI_2_ranking)
print()
PureBPD_total_target_nodes,PureBPD_size_nodes,PureBPD_steps_attack=Pareto_ideal_targeted_attack_2objs(multiplex_structure,N,BPalgorithm_ranking)
print()

###### Results
results_dict={
"HDA":HDA_size_nodes,
"EMD1step":EMD1step_size_nodes,
"EMD":EMD_size_nodes,
"DCI":DCI_size_nodes,
"DCIz":DCIz_size_nodes,
"kcore-EMD1step":kcore_EMD1step_ideal_size_nodes,
"kcore-EMD":kcore_EMD_ideal_size_nodes,
"kcore-DCI":kcore_DCI_size_nodes,
"kcore-DCIz":kcore_DCIz_size_nodes,
"CoreHD-EMD1step":CoreHD_EMD1step_size_nodes,
"CoreHD-EMD":CoreHD_EMD_size_nodes,
"CoreHD-DCI":CoreHD_DCI_size_nodes,
"CoreHD-DCIz":CoreHD_DCIz_size_nodes,
"CI2-EMD1step":CI_2_EMD1step_size_nodes,
"CI2-EMD":CI_2_EMD_size_nodes,
"CI2-DCI":CI_2_DCI_size_nodes,
"CI2-DCIz":CI_2_DCIz_size_nodes,
"BPD-EMD1step":BPD_EMD1step_size_nodes,
"BPD-EMD":BPD_EMD_size_nodes,
"BPD-DCI":BPD_DCI_size_nodes,
"BPD-DCIz":BPD_DCIz_size_nodes,
"CoreHD - CoreHD":Pure_CoreHD_size_nodes,
"CI2 - CI2":Pure_CI_2_size_nodes,
"BPD - BPD":PureBPD_size_nodes}

print(sorted(results_dict.items(),key= lambda x: x[1]))
