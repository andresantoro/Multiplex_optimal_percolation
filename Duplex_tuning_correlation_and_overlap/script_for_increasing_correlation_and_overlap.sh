#!/bin/bash
LANG=en_US # This is required for having sequence with dot instead of comma

if [ $# -le 2 ]; then
    echo "Insert <layer1> <layer2> <iteration_number> as input! Exit now"
    exit 1
fi

L1="$1"
L2="$2"
iteration_number=$3

#Tuning the rho correlation from -1 to 1
sh tune_rho_adiab.sh ${L1} ${L2} 0 0.1 1 0.005 0.0000001 $iteration_number


# For each of the duplex with a given correlation value, we increase the overlap to 0.2 and 0.4
counter=1

for i in `seq -f "%.1f" -1.0 0.1 1.0`
do
	if [ $counter -eq 11 ]
      then
        i=0.0
      fi
	ls -v ${L1}_rho_${i} ${L2}_rho_${i}  > list_multiplex.txt
	cat ${L1}_rho_${i}  ${L2}_rho_${i} | awk '{ print $1 "\n" $2}' | sort -n | uniq > nodes_unique.txt
	# sleep 1
	#Rewiring the duplex using biased edge rewiring as explained in the paper 
	#"Optimal percolation in correlated multilayer networks with overlap"
	./Rewiring_overlap/tune_overlap list_multiplex.txt nodes_unique.txt I -m 0.4
	###RENAMING THE FILE WITH OVERLAP AND CORRELATION
	for overlap in `seq -f "%.2f" 0.0 0.2 0.40`
	do
		mv ${L1}_${overlap}* ${L1}_overlap_${overlap}_rho_${i}
		mv ${L2}_${overlap}* ${L2}_overlap_${overlap}_rho_${i}
	done
	counter=$((counter +1))
done
counter=-1

#ONCE ALL THE MULTIPLEX WITH GIVEN OVERLAP AND TUNABLE CORRELATION ARE CREATED, WE MOVE ALL THE FILES IN THE PROPER DIRECTORIES

for overlap in `seq -f "%.2f" 0.00 0.2 0.40`
do
	s_val=`echo "$overlap == 0.0" |bc -l`
	### Case 1: overlap is zero
	if [ $s_val -eq 1 ]
	then
		#Create the folder containing all the multiplex with overlap o=0.0 and variable inter-layer degree correlation
		mkdir multiplex_rho_${iteration_number}
		for  rho in `seq -f "%.1f" -1.0 0.1 1.0`
		do
      		rho=`printf "%.1f" $rho`
      		s_val=`echo "$rho == -0.0" | bc -l`
      		if [ $s_val -eq 1 ]
      		then
        		rho=0.0
      		fi
      		#Move all the file inside the newly created folder
			mv ${L1}_rho_${rho} multiplex_rho_${iteration_number}
			mv ${L2}_rho_${rho} multiplex_rho_${iteration_number}
		done
		
		cd multiplex_rho_${iteration_number}
		counter=0
		for  rho in `seq -f "%.1f" -1.0 0.1 1.0`
		do
      		rho=`printf "%.1f" $rho`
      		s_val=`echo "$rho == -0.0" | bc -l`
      		if [ $s_val -eq 1 ]
      		then
        		rho=0.0
      		fi
			ls -v ${L1}_rho_${rho} ${L2}_rho_${rho} > list_multiplex_${counter}.txt
			counter=$((counter + 1))
		done
		cd ..
	else
		mkdir multiplex_overlap_${overlap}_${iteration_number}
		for  rho in `seq -f "%.1f" -1.0 0.1 1.0`
		do
    		rho=`printf "%.1f" $rho`
      		s_val=`echo "$rho == -0.0" | bc -l`
      		if [ $s_val -eq 1 ]
      		then
        		rho=0.0
      		fi
			mv  ${L1}_overlap_${overlap}_rho_${rho} multiplex_overlap_${overlap}_${iteration_number}
			mv  ${L2}_overlap_${overlap}_rho_${rho} multiplex_overlap_${overlap}_${iteration_number}
		done
		cd multiplex_overlap_${overlap}_${iteration_number}
		counter=0
		for  rho in `seq -f "%.1f" -1.0 0.1 1.0`
		do
      		rho=`printf "%.1f" $rho`
      		s_val=`echo "$rho == -0.0" | bc -l`
      		if [ $s_val -eq 1 ]
      		then
        		rho=0.0
      		fi
			ls -v ${L1}_overlap_${overlap}_rho_${rho} ${L2}_overlap_${overlap}_rho_${rho} > list_multiplex_${counter}.txt
			counter=$((counter + 1))
		done
		cd ..
	fi
done
rm *_bottom_to_top