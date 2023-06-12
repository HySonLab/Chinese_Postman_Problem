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

int main(int argc, char **argv) {
    // Load the input graph
    Graph *graph = new Graph("data/sample_input_1.txt");

    // Run the Floyd's algorithm
    graph -> Floyd_algorithm();

    cout << "All-pair shortest paths:\n" << to_str(graph -> shortest_path) << endl;

	// Greedy constructive heuristics
	vector<Edge> sigma = greedy_constructive_heuristic(graph);

	// Dynamic programming
	pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sigma);

	cout << "Cost (greedy constructive heuristics): " << dp.first[0][0] << endl;

	return 0;
}

