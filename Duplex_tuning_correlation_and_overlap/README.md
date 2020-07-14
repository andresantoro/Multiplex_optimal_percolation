This is a collection of algorithms that allows to construct duplex networks having a tunable amount of interlayer degree correlations (Spearman's rank correlation)
and tunable structural overlap, as explained in the paper: "*Optimal percolation in correlated multilayer networks with overlap*" ArXiv preprint: https://arxiv.org/abs/1910.04783

We provide a simple example with a duplex network consisting of two independent realisations of the Erdos-Renyi random graph models with N=1000 nodes and K=3000 links.

Before running the code:
```./script_for_increasing_correlation_and_overlap.sh layer1.txt layer2.txt 1```,
make sure that you have compiled the C codes in both the two folders using ```make```

Starting from this duplex system with N=1000 nodes (with approximatively interlayer degree correlations equal to zero), we tune the interlayer degree correlations from -1 to 1 with steps of 0.1. For each of these new duplex configurations, we increase the structural overlap *o* from approximatively zero to a value of 0.4 with steps of 0.2. and save all the duplexes created in three distinct folders. 
