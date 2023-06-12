// Testing dynamic programming for the load-dependent Chinese postman problem
// Author: Dr. Truong Son Hy
// Copyright 2023

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <thread>
#include <assert.h>

#include "Graph.h"
#include "meta_heuristics.cpp"

using namespace std;

// Main program
int main(int argc, char **argv) {
    // Load the input graph
    Graph *graph = new Graph("data/sample_input_1.txt");

    // Run the Floyd's algorithm
    graph -> Floyd_algorithm();

    cout << "All-pair shortest paths:\n" << to_str(graph -> shortest_path) << endl;

    // Test dynamic programming
    vector<Edge> sequence;
    sequence.clear();
    sequence.push_back(graph -> get_edge(0, 1));
    sequence.push_back(graph -> get_edge(2, 0));
    sequence.push_back(graph -> get_edge(1, 3));
    sequence.push_back(graph -> get_edge(3, 2));

    pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sequence);

    cout << "f:\n" << to_str(dp.first) << endl;
    cout << "direction:\n" << to_str(dp.second) << endl;
    cout << "cost: " << compute_cost(graph, sequence, dp.second) << endl;

    return 0;
}

