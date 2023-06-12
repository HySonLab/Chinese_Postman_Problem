// Testing the Iterative Local Search (ILS) for the load-dependent Chinese postman problem
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
    Graph *graph = new Graph("data/sample_input_2.txt");

    // Run the Floyd's algorithm
    graph -> Floyd_algorithm();

    cout << "All-pair shortest paths:\n" << to_str(graph -> shortest_path) << endl;

	// Iterated Local Search
	pair< vector<Edge>, double> result = Iterated_Local_Search(graph);
	const vector<Edge> sigma = result.first;
	const double cost = result.second;

	// Dynamic programming
	pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sigma);

	cout << "Cost (Iterated Local Search): " << dp.first[0][0] << endl;
	assert(abs(cost - dp.first[0][0]) < 1e-6);

	return 0;
}

