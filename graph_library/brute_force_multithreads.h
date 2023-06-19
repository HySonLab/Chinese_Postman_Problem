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
#include <thread>

#include "Graph.h"

using namespace std;


// +--------------------------------------------------+
// | Brute-Force / Back-Tracking with Multi-Threading |
// +--------------------------------------------------+

static void Back_Tracking(Graph *graph, int i, bool *mask, int *perm, int *best_perm, double &best_cost) {
	const int m = graph -> num_deliver_edges;
	
	for (int j = 0; j < m; ++j) {
		if (mask[j] == false) {
			perm[i] = j;
			mask[j] = true;

			if (i + 1 < m) {
				Back_Tracking(graph, i + 1, mask, perm, best_perm, best_cost);
			} else {
				vector<Edge> sequence;
				sequence.clear();
				for (int k = 0; k < m; ++k) {
					sequence.push_back(Edge(graph -> edges[perm[k]]));
				}

				// Dynamic Programming
				pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sequence);

				const double cost = dp.first[0][0];
				if (cost < best_cost) {
					best_cost = cost;
					for (int k = 0; k < m; ++k) {
						best_perm[k] = perm[k];
					}
				}
			}

			mask[j] = false;
		}
	}
}

// Multi-Threading / Parallel
pair< vector<Edge>, double> Brute_Force_MultiThreads(Graph *graph, const int num_threads = 14) {
	// Information
	assert(graph -> shortest_path.size() == graph -> num_nodes);

	const int num_nodes = graph -> num_nodes;
	const int m = graph -> num_deliver_edges;
	assert(m > 1);

	// Result
	pair< vector<Edge>, double> result;
	result.second = INF;
	
	// Multi-Threading
	std::thread threads[num_threads];

	// Memory for all threads
	vector<bool*> all_mask;
	vector<int*> all_perm;
	vector<int*> all_best_perm;
	vector<double> all_best_cost;

	all_mask.clear();
	all_perm.clear();
	all_best_perm.clear();
	all_best_cost.clear();
	
	for (int i = 0; i < num_threads; ++i) {
		// Mask
		bool *mask = new bool [m];
		all_mask.push_back(mask);

		// Permutation
		int *perm = new int [m];
		int *best_perm = new int [m];

		all_perm.push_back(perm);
		all_best_perm.push_back(best_perm);

		// Best cost
		all_best_cost.push_back(0.0);
	}

	// Process
	int start = 0;
	while (start < m) {
		int finish;
		if (start + num_threads - 1 < m) {
			finish = start + num_threads - 1;
		} else {
			finish = m - 1;
		}

		// Launch threads
		int count = 0;
		for (int j = start; j <= finish; ++j) {
			// Mask
    		bool *mask = all_mask[count];

    		for (int i = 0; i < m; ++i) {
        		mask[i] = false;
    		}

    		// Permutation
    		int *perm = all_perm[count];
    		int *best_perm = all_best_perm[count];

			// Init
			mask[j] = true;
			perm[0] = j;
			all_best_cost[count] = INF;

    		// Back-Tracking
			threads[count] = std::thread(Back_Tracking, graph, 1, mask, perm, best_perm, std::ref(all_best_cost[count]));
			count++;
		}

		// Thread Synchronization
		for (int i = 0; i < count; ++i) {
			threads[i].join();
		}
		
		// Update
		for (int i = 0; i < count; ++i) {
			if (all_best_cost[i] < result.second) {
				result.second = all_best_cost[i];
				vector<Edge> sequence;
    			sequence.clear();
    			for (int k = 0; k < m; ++k) {
        			sequence.push_back(Edge(graph -> edges[all_best_perm[i][k]]));
    			}
				result.first = sequence;
			}
		}

		start = finish + 1;
	}

	// Release memory
	for (int i = 0; i < num_threads; ++i) {
		delete[] all_mask[i];
		delete[] all_perm[i];
		delete[] all_best_perm[i];
	}

	return result;
}

