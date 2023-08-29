// Directed Evolution (DE) for the load-dependent Chinese postman problem
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

#include "../graph_library/Graph.h"
#include "../graph_library/meta_heuristics.h"

using namespace std;


// +-----------------------+
// | k-MOVE (new proposal) |
// +-----------------------+

// For single-thread
pair< vector<Edge>, double> Method_k_MOVE(Graph *graph, const vector<Edge> sigma, const vector<int> k_indices) {
    double best = INF;
    vector<Edge> result;

	// Initialization
	vector<bool> mask;
	mask.clear();
	for (int i = 0; i < sigma.size(); ++i) {
		mask.push_back(true);
	}
	for (int i = 0; i < k_indices.size(); ++i) {
		const int index = k_indices[i];
		assert(index >= 0);
		assert(index < sigma.size());
		mask[index] = false;
	}

	vector<Edge> A;
	A.clear();
	for (int i = 0; i < sigma.size(); ++i) {
		if (mask[i]) {
			A.push_back(Edge(sigma[i]));
		}
	}
	assert(A.size() + k_indices.size() == sigma.size());

	// Put them back
	for (int i = 0; i < k_indices.size(); ++i) {
		const int index = k_indices[i];
		const Edge e = sigma[index];

		// Search for the best place to put the i-th in
		vector<Edge> B;
		B.clear();
		double B_value = INF;

		for (int k = 0; k < A.size(); ++k) {
			vector<Edge> candidate;
			candidate.clear();
			for (int t = 0; t < k; ++t) {
				candidate.push_back(Edge(A[t]));
			}
			candidate.push_back(Edge(e));
			for (int t = k; t < A.size(); ++t) {
				candidate.push_back(Edge(A[t]));
			}

			// Update
			pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, candidate);
			const double cost = dp.first[0][0];
			if (cost < B_value) {
				B_value = cost;
				B = candidate;
			}
		}

		A = B;
		best = B_value;
	}

	result = A;
	assert(result.size() == sigma.size());

	return make_pair(result, best);
}


// +-------------------------+
// | Directed Evolution (DE) |
// +-------------------------+

// Single-thread implementation
pair< vector<Edge>, double> Directed_Evolution(Graph *graph, const int k, const int num_variants = 10, const int early_stop = 5, const int num_iterations = 75, const bool verbose = true) {
    // New Greedy Constructive Heuristic
    pair< vector<Edge>, double > greedy = Greedy_Constructive_Heuristic(graph);
    vector<Edge> sigma_star = greedy.first;
    double best = greedy.second;

	// Number of edges
	const int m = sigma_star.size();

	// Mask
	vector<bool> mask;
	for (int i = 0; i < m; ++i) {
		mask.push_back(false);
	}

    // Iterative
	int count = 0;
    for (int iter = 1; iter <= num_iterations; ++iter) {
        // Check if improved
		bool improve = false;

		// Random indices
		vector< vector<int> > all_indices;
		all_indices.clear();

		for (int variant = 0; variant < num_variants; ++variant) {
			for (int i = 0; i < m; ++i) {
				mask[i] = false;
			}

			vector<int> k_indices;
			k_indices.clear();
			for (int i = 0; i < k; ++i) {
				while (true) {
					const int index = rand() % m;
					if (mask[index] == false) {
						mask[index] = true;
						k_indices.push_back(index);
						break;
					}
				}
			}
			assert(k_indices.size() == k);
			all_indices.push_back(k_indices);
		}
		assert(all_indices.size() == num_variants);

        // k-MOVE
		vector<Edge> sigma;
		sigma.clear();
		for (int i = 0; i < m; ++i) {
			sigma.push_back(Edge(sigma_star[i]));
		}

		for (int variant = 0; variant < num_variants; ++variant) {
       		pair< vector<Edge>, double> Result_k_MOVE = Method_k_MOVE(graph, sigma, all_indices[variant]);

			// Update
			if (Result_k_MOVE.second < best) {
				best = Result_k_MOVE.second;
				sigma_star = Result_k_MOVE.first;
				improve = true;
			}
		}

		/*
		// Random exchange
        sigma = random_exchange(sigma_star);

        // 1-OPT
        pair< vector<Edge>, double> Result_1_OPT = Method_1_OPT(graph, sigma);

        // 2-OPT
        pair< vector<Edge>, double> Result_2_OPT = Method_2_OPT(graph, sigma);

        // 2-EXCHANGE
        pair< vector<Edge>, double> Result_2_EXCHANGE = Method_2_EXCHANGE(graph, sigma);

        // Update
        if (Result_1_OPT.second < best) {
            best = Result_1_OPT.second;
            sigma_star = Result_1_OPT.first;
			improve = true;
        }

        if (Result_2_OPT.second < best) {
            best = Result_2_OPT.second;
            sigma_star = Result_2_OPT.first;
			improve = true;
        }

        if (Result_2_EXCHANGE.second < best) {
            best = Result_2_EXCHANGE.second;
            sigma_star = Result_2_EXCHANGE.first;
			improve = true;
        }
		*/
		
		// Check early stopping
		if (!improve) {
			count += 1;
		} else {
			count = 0;
		}

		if (count == early_stop) {
			cout << "Early stop!" << endl;
			break;
		}

        if (verbose) {
            cout << "Done " << iter << " iterations." << endl;
        }
    }

    return make_pair(sigma_star, best);
}

