# Codes for "Optimal percolation in correlated multilayer networks with overlap"
The Python implementation of the targeted attack strategies presented in the paper can be found in the folder: ```Targeted_methods_python```.

The collection of algorithms for tuning the interlayer degree correlation and structural edge overlap of a duplex network can be found in the folder:```Duplex_tuning_correlation_and_overlap```. Part of this code is based on the libraries developed in the following resources. For more information please visit:

http://www.complex-networks.net/
https://github.com/KatolaZ/mammult






Abstract
----------
Multilayer networks have been found to be prone to abrupt cascading failures under random and targeted attacks, but most of the targeting algorithms proposed so far have been mainly tested on uncorrelated systems. Here we show that the size of the critical percolation set of a multilayer network is substantially affected by the presence of interlayer degree correlations and edge overlap. We provide extensive numerical evidence which confirms that the state-of-the-art optimal percolation strategies consistently fail to identify minimal percolation sets in synthetic and real-world correlated multilayer networks, thus overestimating their robustness. We propose two targeting algorithms, based on the local estimation of path disruptions away from a given node, and a family of Pareto-efficient strategies that take into account both intralayer and interlayer heuristics and can be easily extended to multiplex networks with an arbitrary number of layers. We show that these strategies consistently outperform existing attacking algorithms, on both synthetic and real-world multiplex networks, and provide some interesting insights into the interplay of correlations and overlap in determining the hyperfragility of real-world multilayer networks. Overall, the results presented in the paper suggest that we are still far from having fully identified the salient ingredients determining the robustness of multiplex networks to targeted attacks.


Citation
----------
A. Santoro & V. Nicosia (2019) "*Optimal percolation in correlated multilayer networks with overlap*" ArXiv preprint: https://arxiv.org/abs/1910.04783



## License

This project is licensed under the GNU GPLv3 license - see the [LICENSE](LICENSE) file for details
