
#!/bin/sh

#Launch the code starting from two different layers and gradually increase the overlap up to 0.45
#The code will stop when the criteria |o- 0.45| < TOLL=0.0005 is true.
#It will print on file different duplexes according to the value of overlap with steps of 0.1 
./../tune_overlap list_multiplex_different.txt nodes_unique.txt I -m 0.45 -t 0.0005

