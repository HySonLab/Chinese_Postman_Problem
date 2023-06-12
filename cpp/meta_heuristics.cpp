// Meta-heuristics for the load-dependent Chinese postman problem
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


// +---------------------+
// | Dynamic Programming |
// +---------------------+
pair< vector< vector<double> >, vector<int> > dynamic_programming(Graph *graph, const vector<Edge> sequence) {
	// Start node
	const int start_node = graph -> start_node;

	// Check if the graph has computed Floyd's algorithm before
    assert(graph -> shortest_path.size() == graph -> num_nodes);
	const double W = graph -> W;

	// Number of edges in this list
	const int m = sequence.size();

	// Initialize
    vector< vector<double> > f;
    f.clear();
    for (int k = 0; k < m; ++k) {
        vector<double> vect;
        vect.clear();
        vect.push_back(INF);
        vect.push_back(INF);
        f.push_back(vect);
    }

    vector<int> direction;
    direction.clear();
    for (int k = 0; k < m; ++k) {
        direction.push_back(-1);
    }

	// Special case if there is only one edge in the list
	if (m == 1) {
		// Current edge's information
        const int i = sequence[0].first;
        const int j = sequence[0].second;
        const double q = sequence[0].q;
        const double d = sequence[0].d;

		const double option_0 = (W + q / 2) * d + (W + q) * graph -> shortest_path[start_node][i] + W * graph -> shortest_path[j][start_node];
		const double option_1 = (W + q / 2) * d + (W + q) * graph -> shortest_path[start_node][j] + W * graph -> shortest_path[i][start_node];

		if (option_0 <= option_1) {
			f[0][0] = option_0;
			direction[0] = 0;
		} else {
			f[0][0] = option_1;
			direction[0] = 1;
		}

		return make_pair(f, direction);
	}

	// For tracing
	vector< vector<int> > choice;
    choice.clear();
    for (int k = 0; k < m; ++k) {
        vector<int> vect;
        vect.clear();
        vect.push_back(-1);
        vect.push_back(-1);
        choice.push_back(vect);
    }

	// Compute the Q
	vector<double> Q;
	for (int k = 0; k < m; ++k) {
		Q.push_back(0);
	}
	Q[m - 1] = sequence[m - 1].q;
	for (int k = m - 2; k >= 0; --k) {
		Q[k] = Q[k + 1] + sequence[k].q;
	}

	// Dynamic Programming
	for (int k = m - 1; k >= 0; --k) {
		for (int c = 0; c <= 1; ++c) {
			// Skip this case
			if ((k == 0) && (c == 1)) {
				continue;
			}

			// Current edge's information
			const int i = sequence[k].first;
			const int j = sequence[k].second;
			const double q = sequence[k].q;
			const double d = sequence[k].d;

			// Previous node
			int prev = -1;
			if (k == 0) {
				prev = start_node;
			} else {
            	if (c == 0) {
            		// If the previous edge is not flipped
            		prev = sequence[k - 1].second;
				} else {
                	// If the previous edge is flipped
                	prev = sequence[k - 1].first;
            	}
			}

			double option_0;
			double option_1;

			if (k == m - 1) {
				// For the last edge
				option_0 = (W + q / 2) * d + (W + q) * graph -> shortest_path[prev][i] + W * graph -> shortest_path[j][start_node];
				option_1 = (W + q / 2) * d + (W + q) * graph -> shortest_path[prev][j] + W * graph -> shortest_path[i][start_node];
			} else {
				// For all other edges
				option_0 = (W + Q[k] - q / 2) * d + (W + Q[k]) * graph -> shortest_path[prev][i] + f[k + 1][0];
               	option_1 = (W + Q[k] - q / 2) * d + (W + Q[k]) * graph -> shortest_path[prev][j] + f[k + 1][1];
			}

			if (option_0 <= option_1) {
            	f[k][c] = option_0;
                choice[k][c] = 0;
            } else {
                f[k][c] = option_1;
                choice[k][c] = 1;
            }
		}
	}

	// Tracing
	int k = 0;
	int c = 0;
	while (k < m) {
		direction[k] = choice[k][c];
		c = choice[k][c];
		++k;
	}

	return make_pair(f, direction);
}


// +--------------+
// | Compute cost |
// +--------------+
double compute_cost(Graph *graph, const vector<Edge> sequence, const vector<int> direction) {
	// Check if the graph has computed Floyd's algorithm before
    assert(graph -> shortest_path.size() == graph -> num_nodes);
    const double W = graph -> W;

    // Number of edges in this list
    const int m = sequence.size();

	// Compute Q
	double Q = 0;
	for (int k = 0; k < sequence.size(); ++k) {
		Q += sequence[k].q;
	}

	// Compute cost
	double cost = 0;
	
	// Previous node
	int prev = 0;

	for (int k = 0; k < sequence.size(); ++k) {
		// Current edge's information
		const int i = sequence[k].first;
        const int j = sequence[k].second;
        const double q = sequence[k].q;
        const double d = sequence[k].d;

		// First edge
		if (direction[k] == 0) {
			cost += (W + Q) * graph -> shortest_path[prev][i];
		} else {
			cost += (W + Q) * graph -> shortest_path[prev][j];
		}

		// Last edge
		if (k == m - 1) {
			if (direction[k] == 0) {
                cost += W * graph -> shortest_path[j][0];
            } else {
                cost += W * graph -> shortest_path[i][0];
            }
		}

		// Deliver cost
		cost += (W + Q - q / 2) * d;

		// Reduce load
		Q -= q;

		// Update the previous node
		if (direction[k] == 0) {
			prev = j;
		} else {
			prev = i;
		}
	}

	return cost;
}


// +--------------------------------+
// | Greedy Constructive Heuristics |
// +--------------------------------+
pair< vector<Edge>, double> Greedy_Constructive_Heuristic(Graph *graph) {
	// Information
	const int num_nodes = graph -> num_nodes;
	const int num_edges = graph -> num_edges;
		
	vector<Edge> edges;
	edges.clear();
	for (int k = 0; k < num_edges; ++k) {
		edges.push_back(Edge(graph -> edges[k]));
	}

	// Sort the list of edges
	sort(edges.begin(), edges.end());

	// Result
	vector<Edge> sigma_star;
	sigma_star.clear();
	double best;

	// Algorithm
	for (int i = 0; i < num_edges; ++i) {
		double z_min = INF;
		vector<Edge> sigma_prime;

		// The i-th edge
		Edge e = Edge(edges[i]);

		for (int j = 0; j <= i; ++j) {
			// Create another list of edges by adding the i-th edge into the j-th position of sigma_star
			vector<Edge> sigma;
			sigma.clear();
			for (int k = 0; k < j; ++k) {
				sigma.push_back(Edge(sigma_star[k]));
			}
			sigma.push_back(e);
			for (int k = j + 1; k <= i; ++k) {
				sigma.push_back(Edge(sigma_star[k - 1]));
			}

			// Dynamic programming
			pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sigma);

			// Update
			const double z = dp.first[0][0];
			if (z < z_min) {
				z_min = z;
				sigma_prime = sigma;
			}
		}

		// Update sigma_star
		sigma_star = sigma_prime;
		best = z_min;
	}

	return make_pair(sigma_star, best);
}


// +-------------------------------------------+
// | Random exchange with probability 1/factor |
// +-------------------------------------------+
vector<Edge> random_exchange(const vector<Edge> edges, const int factor = 5) {
	// Copy
	vector<Edge> result;
	result.clear();
	for (int i = 0; i < edges.size(); ++i) {
		result.push_back(Edge(edges[i]));
	}

	// Random exchange
	for (int i = 0; i < edges.size(); ++i) {
		if (rand() % factor == 0) {
			const int j = rand() % edges.size();
			swap(result[i], result[j]);
		}
	}

	return result;
}


// +-------+
// | 1-OPT |
// +-------+
pair< vector<Edge>, double> Method_1_OPT(Graph *graph, const vector<Edge> sigma) {
	double best = INF;
	vector<Edge> result;
	
	// Search
	for (int i = 0; i < sigma.size(); ++i) {
		Edge e = sigma[i];

		// Move the i-th edge to the j-th position
		for (int j = 0; j < sigma.size(); ++j) {
			vector<Edge> sequence;
			sequence.clear();
			for (int k = 0; k < j; ++k) {
				if (k != i) {
					sequence.push_back(Edge(sigma[k]));
				}
			}
			sequence.push_back(Edge(e));
			for (int k = j; k < sigma.size(); ++k) {
				if (k != i) {
					sequence.push_back(Edge(sigma[k]));
				}
			}
			assert(sequence.size() == sigma.size());

			// Dynamic programming
			pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sequence);

			// Update
			const double cost = dp.first[0][0];
			if (cost < best) {
				best = cost;
				result = sequence;
			}
		}
	}

	return make_pair(result, best);
}


// +-------+
// | 2-OPT |
// +-------+
pair< vector<Edge>, double> Method_2_OPT(Graph *graph, const vector<Edge> sigma) {
    double best = INF;
    vector<Edge> result;

	// Copy
	vector<Edge> sequence;
	sequence.clear();
	for (int k = 0; k < sigma.size(); ++k) {
		sequence.push_back(Edge(sigma[k]));
	}

	// Search
    for (int i = 0; i < sigma.size(); ++i) {
        Edge e = sigma[i];

        // Swap the i-th and j-th edges
        for (int j = i + 1; j < sigma.size(); ++j) {
			// Swap
			swap(sequence[i], sequence[j]);

            // Dynamic programming
            pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sequence);

            const double cost = dp.first[0][0];
            if (cost < best) {
                best = cost;
                result = sequence;
            }

			// Swap back
			swap(sequence[j], sequence[i]);
        }
    }

    return make_pair(result, best);
}


// +------------+
// | 2-EXCHANGE |
// +------------+
pair< vector<Edge>, double> Method_2_EXCHANGE(Graph *graph, const vector<Edge> sigma) {
    double best = INF;
    vector<Edge> result;

    // Search
    for (int i = 0; i < sigma.size(); ++i) {
        for (int j = i + 1; j < sigma.size(); ++j) {
            // Reverse the order of edges from i-th to j-th
           	vector<Edge> sequence;
			sequence.clear();
			for (int k = 0; k < i; ++k) {
				sequence.push_back(Edge(sigma[k]));
			}
			for (int k = j; k >= i; --k) {
				sequence.push_back(Edge(sigma[k]));
			}
			for (int k = j + 1; k < sigma.size(); ++k) {
				sequence.push_back(Edge(sigma[k]));
			}
			assert(sequence.size() == sigma.size());

            // Dynamic programming
            pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sequence);

            const double cost = dp.first[0][0];
            if (cost < best) {
                best = cost;
                result = sequence;
            }
        }
    }

    return make_pair(result, best);
}


// +--------------------------------------------+
// | Iterated Local Search (ILS) Meta-heuristic |
// +--------------------------------------------+
pair< vector<Edge>, double> Iterated_Local_Search(Graph *graph, const int k_max = 100) {
	// Greedy constructive heuristic
	pair< vector<Edge>, double > greedy  = Greedy_Constructive_Heuristic(graph);
	vector<Edge> sigma_star = greedy.first;
	double best = greedy.second;

	// Iterative
	for (int k = 1; k <= k_max; ++k) {
		// Random exchange
		vector<Edge> sigma = random_exchange(sigma_star);
		
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
		}

		if (Result_2_OPT.second < best) {
            best = Result_2_OPT.second;
            sigma_star = Result_2_OPT.first;
        }

		if (Result_2_EXCHANGE.second < best) {
            best = Result_2_EXCHANGE.second;
            sigma_star = Result_2_EXCHANGE.first;
        }
	}

	return make_pair(sigma_star, best);
}


// +---------------------------------------------------+
// | Variable Neighborhood Search (VNS) Meta-heuristic |
// +---------------------------------------------------+
pair< vector<Edge>, double> Variable_Neighborhood_Search(Graph *graph, const int k_max = 100) {
    // Greedy constructive heuristic
    pair< vector<Edge>, double > greedy  = Greedy_Constructive_Heuristic(graph);
    vector<Edge> sigma_star = greedy.first;
    double best = greedy.second;

    // Iterative
    for (int k = 1; k <= k_max; ++k) {
        // Random exchange
        vector<Edge> sigma = random_exchange(sigma_star);
		
		// Compute the cost of sigma
		pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sigma);
		double cost = dp.first[0][0];

		// 2-EXCHANGE
        pair< vector<Edge>, double> Result_2_EXCHANGE = Method_2_EXCHANGE(graph, sigma);
		
		if (Result_2_EXCHANGE.second < cost) {
			cost = Result_2_EXCHANGE.second;
			sigma = Result_2_EXCHANGE.first;
		} else {
        	// 1-OPT
        	pair< vector<Edge>, double> Result_1_OPT = Method_1_OPT(graph, sigma);

			if (Result_1_OPT.second < cost) {
				cost = Result_1_OPT.second;
				sigma = Result_1_OPT.first;		
			} else {
				// 2-OPT
        		pair< vector<Edge>, double> Result_2_OPT = Method_2_OPT(graph, sigma);

				if (Result_2_OPT.second < cost) {
					cost = Result_2_OPT.second;
					sigma = Result_2_OPT.first;
				}
			}
		}

    	// Update
		if (cost < best) {
			best = cost;
			sigma_star = sigma;
		}
	}

    return make_pair(sigma_star, best);
}

