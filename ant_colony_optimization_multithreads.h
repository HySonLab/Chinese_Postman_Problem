// Ant colony optimization (ACO) with multi-threading for the load-dependent Chinese postman problem
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
#include "meta_heuristics.h"

using namespace std;

// +----------------------------------------------------+
// | Ant Colony Optimization (ACO) with multi-threading |
// +----------------------------------------------------+

// Search 
static int state_index(const int first, const int second, const vector<Edge> &edges) {
	for (int i = 0; i < edges.size(); ++i) {
		if ((first == edges[i].first) && (second == edges[i].second)) {
			return 2 * i;
		}
		if ((first == edges[i].second) && (second == edges[i].first)) {
            return 2 * i + 1;
        }
	}
	return -1;
}

// Random number from 0 to 1
static double rand_double(const int MAX_RANGE = 1000) {
	double r = rand() % (MAX_RANGE + 1);
	return r / (double)(MAX_RANGE);
}

// Ant's random walk
static void ant_walk(
	Graph *graph, 
	double **tau, 
	double **eta, 
	const int origin_state, 
	int *state_path, 
	Edge *sequence, 
	int *direction, 
	double &cost
) {
	// Information
	const int m = graph -> num_deliver_edges;
	const int start_node = graph -> start_node;

	// Mask if an edge is delivered
	vector<bool> mask;
	mask.clear();
	for (int k = 0; k < m; ++k) {
		mask.push_back(false);
	}

	// Random walk
	int state = origin_state;
	for (int k = 0; k < m; ++k) {
		// List of possible next states to visit
		vector< pair<double, int> > candidate;
		candidate.clear();

		for (int i = 0; i < m; ++i) {
			if (!mask[i]) {
				candidate.push_back(make_pair(- tau[state][2 * i] * eta[state][2 * i], 2 * i));
				candidate.push_back(make_pair(- tau[state][2 * i + 1] * eta[state][2 * i + 1], 2 * i + 1));
			}
		}

		// Sort based on the potential
		sort(candidate.begin(), candidate.end());

		// Reverse the sign
		for (int i = 0; i < candidate.size(); ++i) {
			candidate[i].first = abs(candidate[i].first);
		}

		// Compute the probabilities
		double sum = 0;
		for (int i = 0; i < candidate.size(); ++i) {
			sum += candidate[i].first;
		}

		for (int i = 0; i < candidate.size(); ++i) {
			candidate[i].first /= sum;
		}

		// Random a number from 0 to 1
		const double r = rand_double();

		// Sampling
		int next_state = candidate[0].second;
		double s = 0;
		for (int i = 0; i < candidate.size(); ++i) {
			s += candidate[i].first;
			if (r <= s) {
				next_state = candidate[i].second;
			} else {
				break;
			}
		}

		// Add the corresponding edge to the sequence
		const int edge_index = next_state / 2;
		assert(edge_index >= 0);
		assert(edge_index < m);
		assert(mask[edge_index] == false);

		sequence[k] = Edge(graph -> deliver_edges[edge_index]);

		// Mask it
		mask[edge_index] = true;

		// Update the state
		state = next_state;
		state_path[k] = state;
	}

	// Dynamic programming
	vector<Edge> seq;
	seq.clear();
	for (int i = 0; i < m; ++i) {
		seq.push_back(Edge(sequence[i]));
	}

	pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, seq);

	// Return
	cost = dp.first[0][0];

	for (int i = 0; i < m; ++i) {
		direction[i] = dp.second[i];
	}

	for (int i = 0; i < m; ++i) {
		if (direction[i] == 1) {
			if (state_path[i] % 2 == 0) {
				state_path[i] += 1;
			} else {
				state_path[i] -= 1;
			}
		}
	}
}

// Ant Colony Optimization (ACO) with multi-threading
pair< vector<Edge>, double> Ant_Colony_Optimization_MultiThreads(
	Graph *graph, 
	const int k_max = 1000, 
	const int num_ants = 10, 
	const bool verbose = false,
	const double rho = 0.8, // For updating pheromone
	const int num_threads = 14
) {
	// Epsilon constant
	const double epsilon = 1e-3;

	// Result
	pair< vector<Edge>, double> result;

	// Information
	const double W = graph -> W;
	const int m = graph -> num_deliver_edges;
	const int start_node = graph -> start_node;

	// Sum of loads
	double sum_Q = 0.0;
	for (int i = 0; i < m; ++i) {
		sum_Q += graph -> deliver_edges[i].q;
	}

	// Number of states
	const int num_states = 2 * m + 1;

	// Original state
	const int origin_state = 2 * m;

	// Pheromone
	double **tau = new double* [num_states];
	for (int i = 0; i < num_states; ++i) {
		tau[i] = new double [num_states];
		for (int j = 0; j < num_states; ++j) {
			tau[i][j] = 0.0;
		}
	}

	// Deposit some pheromone
	for (int i = 0; i < m; ++i) {
		// From the original state
		tau[origin_state][2 * i] = epsilon;
		tau[origin_state][2 * i + 1] = epsilon;

		// Between other states
		for (int j = i + 1; j < m; ++j) {
			tau[2 * i][2 * j] = epsilon;
			tau[2 * i + 1][2 * j] = epsilon;
			tau[2 * i][2 * j + 1] = epsilon;
            tau[2 * i + 1][2 * j + 1] = epsilon;
		}
	}

	// Greedy constructive heuristic
    pair< vector<Edge>, double > greedy = Greedy_Constructive_Heuristic(graph);
	assert(greedy.first.size() == m);
	result = greedy;

	// Dynamic programming to determine the direction
	pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, greedy.first);
	assert(dp.second.size() == m);
	assert(abs(dp.first[0][0] - greedy.second) < 1e-6);

	// Initialize the pheromone
	int prev_state = origin_state;
	for (int i = 0; i < m; ++i) {
		int first, second;
		if (dp.second[i] == 0) {
			first = greedy.first[i].first;
			second = greedy.first[i].second;
		} else {
			first = greedy.first[i].second;
            second = greedy.first[i].first;
		}
		
		const int state = state_index(first, second, graph -> deliver_edges);
		assert(state >= 0);
		assert(state < 2 * m);

		// Increase pheromone on this path
		tau[prev_state][state] += epsilon;

		prev_state = state;
	}

	// Priori knowledge
	double **eta = new double* [num_states];
	for (int i = 0; i < num_states; ++i) {
		eta[i] = new double [num_states];
		for (int j = 0; j < num_states; ++j) {
			eta[i][j] = 0.0;
		}
	}

	// Initialize priori knowledge from the original state to other states
	for (int k = 0; k < m; ++k) {
		const int i = graph -> deliver_edges[k].first;
		const int j = graph -> deliver_edges[k].second;
		const double d = graph -> deliver_edges[k].d;
		const double q = graph -> deliver_edges[k].q;

		double c;

		// State 2 * k
		c = (W + sum_Q) * graph -> shortest_path[start_node][i] + (W + sum_Q - q / 2) * d;
		eta[origin_state][2 * k] = 1.0 / sqrt(c);

		// State 2 * k + 1
		c = (W + sum_Q) * graph -> shortest_path[start_node][j] + (W + sum_Q - q / 2) * d;
		eta[origin_state][2 * k + 1] = 1.0 / sqrt(c);
	}

	// Initialize priori knowledge between other states
	for (int k1 = 0; k1 < m; ++k1) {
		const int i1 = graph -> deliver_edges[k1].first;
		const int j1 = graph -> deliver_edges[k1].second;
        
		for (int k2 = 0; k2 < m; ++k2) {
			const int i2 = graph -> deliver_edges[k2].first;
			const int j2 = graph -> deliver_edges[k2].second;
			const double d = graph -> deliver_edges[k2].d;
			const double q = graph -> deliver_edges[k2].q;

			double c;

			// State 2 * k1 to 2 * k2
			// i1 -> j1 -> i2 -> j2
			c = (W + q) * graph -> shortest_path[j1][i2] + (W + q / 2) * d;
			eta[2 * k1][2 * k2] = 1.0 / sqrt(c);

			// State 2 * k1 + 1 to 2 * k2
			// j1 -> i1 -> i2 -> j2
			c = (W + q) * graph -> shortest_path[i1][i2] + (W + q / 2) * d;
			eta[2 * k1 + 1][2 * k2] = 1.0 / sqrt(c);

			// State 2 * k1 to 2 * k2 + 1
			// i1 -> j1 -> j2 -> i2
			c = (W + q) * graph -> shortest_path[j1][j2] + (W + q / 2) * d;
			eta[2 * k1][2 * k2 + 1] = 1.0 / sqrt(c);

			// State 2 * k1 + 1 to 2 * k2 + 1
			// j1 -> i1 -> j2 -> i2
			c = (W + q) * graph -> shortest_path[i1][j2] + (W + q / 2) * d;
			eta[2 * k1 + 1][2 * k2 + 1] = 1.0 / sqrt(c);
		}
	}

	// Memory for multi-threading
	std::thread threads[num_threads];	

	int **state_path = new int* [num_threads];
	Edge **sequence = new Edge* [num_threads];
	int **direction = new int* [num_threads];
	double *cost = new double [num_threads];

	for (int t = 0; t < num_threads; ++t) {
		state_path[t] = new int [m];
		sequence[t] = new Edge [m];
		direction[t] = new int [m];
	}

	// The amount of pheromone to increase
	double **delta_tau = new double* [num_states];
	for (int i = 0; i < num_states; ++i) {
		delta_tau[i] = new double [num_states];
	}

	// Process
	for (int k = 1; k <= k_max; ++k) {
		// The amount of pheromone to increase in this iteration
		for (int i = 0; i < num_states; ++i) {
			for (int j = 0; j < num_states; ++j) {
				delta_tau[i][j] = 0.0;
			}
		}

		// For each ant
		int start = 0;
		while (start < num_ants) {
			// Batch
			int finish;
			if (start + num_threads - 1 < num_ants) {
				finish = start + num_threads - 1;
			} else {
				finish = num_ants - 1;
			}

			// Multi-Threads
			int count = 0;
			for (int ant = start; ant <= finish; ++ant) {
				// Random walk
            	threads[count] = std::thread(ant_walk, graph, tau, eta, origin_state, state_path[count], sequence[count], direction[count], std::ref(cost[count]));
				++count;
			}

			// Threads Synchronization
			for (int t = 0; t < count; ++t) {
				threads[t].join();
			}

			// Update the trace of pheromone
			for (int t = 0; t < count; ++t) {
				int state = origin_state;
				for (int i = 0; i < m; ++i) {
					assert(state_path[t][i] >= 0);
					assert(state_path[t][i] < 2 * m);

					// Increase the pheromone on this path
					delta_tau[state][state_path[t][i]] += 1.0 / sqrt(cost[t]);

					state = state_path[t][i];
				}

				// Update the best solution
				if (cost[t] < result.second) {
					result.first.clear();
					for (int i = 0; i < m; ++i) {
						result.first.push_back(sequence[t][i]);
					}
					result.second = cost[t];
				}

            }

			start = finish + 1;
		}

		// Increase the pheromone
		for (int i = 0; i < num_states; ++i) {
			for (int j = 0; j < num_states; ++j) {
				tau[i][j] = rho * tau[i][j] + delta_tau[i][j];
			}
		}

		if (verbose) {
			if ((k + 1) % 10 == 0) {
				cout << "Completed " << k << " iterations." << endl;
			}
		}
	}

	// Release memory
	for (int i = 0; i < num_states; ++i) {
		delete[] tau[i];
		delete[] eta[i];
		delete[] delta_tau[i];
	}
	delete[] tau;
	delete[] eta;
	delete[] delta_tau;
	delete[] state_path;
	delete[] sequence;
	delete[] direction;

	return result;
}

