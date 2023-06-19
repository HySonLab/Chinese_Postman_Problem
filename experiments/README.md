# Implementation in C++ for the load-dependent Chinese postman problem
## Author: Dr. Truong Son Hy
## Copyright 2023


### Experiments (```experiments/```)

**Tests for the basic graph algorithms and data structures:**
* Example of how to load / save graph from / to file, use the object and run the Floyd's algorithm: ```test_graph.cpp```.
* Test the dynamic programming for finding the optimal directions of edges: ```test_dp.cpp```.

**Meta-heuristics:**
* Test the **Greedy Constructive Heuristics**: ```test_greedy.cpp```. Usage: ```g++ test_greedy.cpp -o test_greedy```, then ```./test_greedy [file name]```.
* Test the **Iterative Local Search (ILS)**: ```test_ils.cpp```. Usage: ```g++ test_ils.cpp -o test_ils```, then ```./test_ils [file name]```.
* Test the **Iterative Local Search (ILS)** with multi-threading (much faster): ```test_ils_multithreads.cpp```. Usage: ```g++ test_ils_multithreads.cpp -lpthread```, then ```./a.out [file name]```.
* Test the **Variable Neighborhood Search (VNS)**: ```test_vns.cpp```. Usage: ```g++ test_vns.cpp -o test_vns```, then ```./test_vns [file name]```.

**Evolutionary Algorithm (EA):**
* Test the Evolutionary Algorithm: ```test_ea.cpp```. Usage: ```g++ test_ea.cpp -o test_ea```, then ```./test_ea [file name]```.
* Test the Evolutionary Algorithm with multi-threading (much faster): ```test_ea_multithreads.cpp```. Usage: ```g++ test_ea_multithreads.cpp -lpthread```, then ```./a.out [file name]```.

**Ant Colony Optimization (ACO):**
* Test the Ant Colony Optimization: ```test_aco.cpp```. Usage: ```g++ test_aco.cpp -o test_aco```, then ```./test_aco [file name]```.
* Test the Ant Colony Optimization with multi-threading (much faster): ```test_aco_multithreads.cpp```. Usage: ```g++ test_aco_multithreads.cpp -lpthread```, then ```./a.out [file name]```.

**Brute-Force / Back-Tracking:**
* Test the Brute-Force / Back-Tracking: ```test_brute_force.cpp```. Usage: ```g++ test_brute_force.cpp -o test_brute_force```, then ```./test_brute_force [file name]```.
* Test the Brute-Force / Back-Tracking with multi-threading (much faster): ```test_brute_force_multithreads.cpp```. Usage: ```g++ test_brute_force_multithreads.cpp -lpthread```, then ```./a.out [file name]```.

**Scripts:**
* Example (script) of running all the baselines including brute-force: ```test_all_baselines.sh```. Usage: ```sh test_all_baselines.sh```.
* Example (script) of running only meta-heuristics, evolutionary algorithm and ant colony optimization: ```test_heuristics.sh```. Usage: ```sh test_heuristics.sh```.

