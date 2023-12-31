// Testing Brute-Force / Back-Tracking with multi-threading for the load-dependent Chinese postman problem
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

#include "../graph_library/brute_force_multithreads.h"

using namespace std;
using namespace std::chrono;

// Main program
int main(int argc, char **argv) {
	// Starting timepoint
    auto start = high_resolution_clock::now();

    cout << "File name: " << argv[1] << endl;

	// Load the input graph
    Graph *graph = new Graph(argv[1]);

	cout << "Number of nodes: " << graph -> num_nodes << endl;
    cout << "Number of edges: " << graph -> num_edges << endl;
	cout << "Number of deliver-edges (q > 0): " << graph -> num_deliver_edges << endl;

    // Run the Floyd's algorithm
    graph -> Floyd_algorithm();

    // cout << "All-pair shortest paths:\n" << to_str(graph -> shortest_path) << endl;

	// Brute-Force / Back-Tracking
	pair< vector<Edge>, double> result = Brute_Force_MultiThreads(graph);
	const vector<Edge> sigma = result.first;
	const double cost = result.second;

	// Dynamic programming
	pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sigma);

	cout << "Cost (Brute-Force with multi-threading): " << dp.first[0][0] << endl;
	assert(abs(cost - dp.first[0][0]) < 1e-6);

	// Ending timepoint
    auto stop = high_resolution_clock::now();

    // Duration
	auto ms_duration = duration_cast<milliseconds>(stop - start);
    cout << "Running time (milliseconds): " << ms_duration.count() << endl;
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Running time (seconds): " << duration.count() << endl << endl;

	for (int i = 0; i < sigma.size(); ++i) {
		cout << sigma[i].first << " " << sigma[i].second << " " << dp.second[i] << endl;
	}

	return 0;
}

