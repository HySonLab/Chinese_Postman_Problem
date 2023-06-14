// Testing the Greedy Constructive Heuristics for the load-dependent Chinese postman problem
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
#include "meta_heuristics.h"

using namespace std;

// Main program
int main(int argc, char **argv) {
	cout << "File name:" << argv[1] << endl;

    // Load the input graph
    Graph *graph = new Graph(argv[1]);

    // Run the Floyd's algorithm
    graph -> Floyd_algorithm();

	// cout << "All-pair shortest paths:\n" << to_str(graph -> shortest_path) << endl;

	// Greedy constructive heuristics
	pair< vector<Edge>, double> result = Greedy_Constructive_Heuristic(graph);
	const vector<Edge> sigma = result.first;
	const double cost = result.second;

	// Dynamic programming
	pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sigma);

	cout << "Cost (greedy constructive heuristics): " << dp.first[0][0] << endl;
	assert(abs(cost - dp.first[0][0]) < 1e-6);

	return 0;
}

