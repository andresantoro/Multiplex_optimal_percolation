#!/bin/sh

#Launch the code starting from two different layers and gradually increase the overlap
#It will print on file different duplexes according to the value of overlap with steps of 0.1 
#In this case it is not possible to reach perfect overlap (i.e. o=1), hence it will
#stop after 5*10^7 iterations to avoid an infinite loop
./../tune_overlap list_multiplex_different.txt nodes_unique.txt I

