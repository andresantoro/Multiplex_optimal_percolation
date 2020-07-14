####
##
## Get a netfile and a list containing the mapping
##
##  new_id old_id
##
## and output a relabelled netfile

import sys

if len(sys.argv) < 3:
    print "Usage: %s <netfile> <relabel_map>" % sys.argv[0]
    sys.exit(1)

#### The relabel map
lines = open(sys.argv[2], "r").readlines()

map = {}

for l in lines:
    if l[0] == "#":
        continue
    elems = l.strip(" \n").split(" ")
    new_l, old_l = int(elems[0]), int(elems[1])
    map[old_l] = new_l
### The original graph
lines = open(sys.argv[1], "r").readlines()

for l in lines:
    elems = l.strip(" \n").split(" ")
    s, d = map[int(elems[0])], map[int(elems[1])]
    print s, d, 
    for i in elems[2:]:
        print i,
    print
    

    
