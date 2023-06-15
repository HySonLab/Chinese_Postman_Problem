# Implementation in C++ for the load-dependent Chinese postman problem
## Author: Dr. Truong Son Hy
## Copyright 2023

### Data generation
* Data generation is based on paper "The Chinese Postman Problem with Load-Dependent Costs" https://pubsonline.informs.org/doi/abs/10.1287/trsc.2017.0774

* Implementation to generate Eulerian cases: ```gen_eulerian.cpp```.

* Implementation to generate Rural Postmam Problem (RPP) cases: ```gen_rural.cpp```.

* In the ```data/``` direction:
- Prefix `E` denotes the Eulerian cases
- Prefix `C` denotes the cases based on Christofides et al.
- Prefix `H` denotes the cases based on Hertz et al.

### Structure
* Graph object's implementation: ```Graph.h```.

* Meta-heuristics' implementation: ```meta_heuristics.h```.

* Evolutionary Algorithm's implementation: ```evolutionary_algorithm.h```.

* Brute-Force / Back-Tracking's implementation: ```brute_force.h```.


### Running
* Example of how to load / save graph from / to file, use the object and run the Floyd's algorithm: ```test_graph.cpp```.

* Test the dynamic programming for finding the optimal directions of edges: ```test_dp.cpp```.

* Test the Greedy Constructive Heuristics: ```test_greedy.cpp```.

Usage: ```g++ test_greedy.cpp -o test_greedy```, then ```./test_greedy [file name]```.

* Test the Iterative Local Search: ```test_ils.cpp```.

Usage: ```g++ test_ils.cpp -o test_ils```, then ```./test_ils [file name]```.

* Test the Variable Neighborhood Search: ```test_vns.cpp```.

Usage: ```g++ test_vns.cpp -o test_vns```, then ```./test_vns [file name]```.

* Test the Evolutionary Algorithm: ```test_ea.cpp```.

Usage: ```g++ test_ea.cpp -o test_ea```, then ```./test_ea [file name]```.

* Test the Brute-Force / Back-Tracking: ```test_brute_force.cpp```.

Usage: ```g++ test_brute_force.cpp -o test_brute_force```, then ```./test_brute_force [file name]```.

* Example (script) of running all the baselines including (slow) brute-force: ```test_all_baselines.sh```.

Usage: ```sh test_all_baselines.sh```.

* Example (script) of running only meta-heuristics and evolutionary algorithm: ```test_heuristics.sh```.

Usage: ```sh test_heuristics.sh```.

