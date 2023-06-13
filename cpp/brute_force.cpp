// Brute-Force / Back-Tracking for the load-dependent Chinese postman problem
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
#include <algorithm>
#include <assert.h>

#include "Graph.h"

using namespace std;


// +-----------------------------+
// | Brute-Force / Back-Tracking |
// +-----------------------------+
void Back_Tracking(Graph *graph, int i, bool *mask, int *perm, int *best_perm, double &best_cost) {
	const int num_edges = graph -> num_edges;
	for (int j = 0; j < num_edges; ++j) {
		if (mask[j] == false) {
			perm[i] = j;
			mask[j] = true;

			if (i + 1 < num_edges) {
				Back_Tracking(graph, i + 1, mask, perm, best_perm, best_cost);
			} else {
				vector<Edge> sequence;
				sequence.clear();
				for (int k = 0; k < num_edges; ++k) {
					sequence.push_back(Edge(graph -> edges[perm[k]]));
				}

				// Dynamic Programming
				pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sequence);

				const double cost = dp.first[0][0];
				if (cost < best_cost) {
					best_cost = cost;
					for (int k = 0; k < num_edges; ++k) {
						best_perm[k] = perm[k];
					}
				}
			}

			mask[j] = false;
		}
	}
}

pair< vector<Edge>, double> Brute_Force(Graph *graph) {
	// Information
	assert(graph -> shortest_path.size() == graph -> num_nodes);

	const int num_nodes = graph -> num_nodes;
	const int num_edges = graph -> num_edges;

	// Mask
	bool *mask = new bool [num_edges];
	for (int i = 0; i < num_edges; ++i) {
		mask[i] = false;
	}

	// Permutation
	int *perm = new int [num_edges];
	int *best_perm = new int [num_edges];
	
	// Best cost
	double best_cost = INF;

	// Back-Tracking
	Back_Tracking(graph, 0, mask, perm, best_perm, best_cost);

	// Solution
	vector<Edge> sequence;
	sequence.clear();
	for (int k = 0; k < num_edges; ++k) {
		sequence.push_back(Edge(graph -> edges[best_perm[k]]));
    }

	// Release memory
	delete[] mask;
	delete[] perm;
	delete[] best_perm;

	return make_pair(sequence, best_cost);
}

