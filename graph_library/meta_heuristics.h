// Meta-heuristics (e.g., greedy, ILS, and VNS) for the load-dependent Chinese postman problem
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

using namespace std;


// +--------------------------------+
// | Greedy Constructive Heuristics |
// +--------------------------------+
pair< vector<Edge>, double> Greedy_Constructive_Heuristic(Graph *graph) {
	// Information
	const int num_nodes = graph -> num_nodes;
	const int m = graph -> num_deliver_edges;
		
	vector<Edge> edges;
	edges.clear();
	for (int k = 0; k < m; ++k) {
		edges.push_back(Edge(graph -> deliver_edges[k]));
	}

	// Sort the list of edges
	sort(edges.begin(), edges.end());

	// Result
	vector<Edge> sigma_star;
	sigma_star.clear();
	double best;

	// Algorithm
	for (int i = 0; i < m; ++i) {
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


// +-----------------------------------------------+
// | Greedy Constructive Heuristics (new proposal) |
// +-----------------------------------------------+
pair< vector<Edge>, double> Greedy_Constructive_Heuristic_2(Graph *graph) {
    // Information
    const int num_nodes = graph -> num_nodes;
    const int m = graph -> num_deliver_edges;

    vector<Edge> edges;
    edges.clear();
    for (int k = 0; k < m; ++k) {
        edges.push_back(Edge(graph -> deliver_edges[k]));
    }

    // Sort the list of edges
    sort(edges.begin(), edges.end());

	// Mask
	vector<bool> mask;
	mask.clear();
	for (int k = 0; k < m; ++k) {
		mask.push_back(false);
	}

    // Result
    vector<Edge> sigma_star;
    sigma_star.clear();
    double best;

	// Algorithm
    for (int i = 0; i < m; ++i) {
		// If this edge is used already, move on
		if (mask[i] == true) {
			continue;
		}
		mask[i] = true;

		// We need to find where to put it in
        double z_min = INF;
        vector<Edge> sigma_prime;

        // The i-th edge
        Edge e = Edge(edges[i]);

		int position = -1;
		int start = -1;
		int finish = -1;

        for (int j = 0; j <= sigma_star.size(); ++j) {
            // Create another list of edges by adding the i-th edge into the j-th position of sigma_star
            vector<Edge> sigma;
            sigma.clear();
            for (int k = 0; k < j; ++k) {
                sigma.push_back(Edge(sigma_star[k]));
            }

			// Add the edge
            sigma.push_back(e);
            			
			for (int k = j; k < sigma_star.size(); ++k) {
                sigma.push_back(Edge(sigma_star[k]));
            }

            // Dynamic programming
            pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sigma);

            // Update
            const double z = dp.first[0][0];
            if (z < z_min) {
                z_min = z;
                sigma_prime = sigma;
				
				position = j;
				if (position == 0) {
					start = graph -> start_node;
				} else {
					if (dp.second[j - 1] == 0) {
						start = sigma_prime[j - 1].second;
					} else {
						start = sigma_prime[j - 1].first;
					}
				}

				if (dp.second[j] == 0) {
					finish = sigma_prime[j].first;
				} else {
					finish = sigma_prime[j].second;
				}
            }
        }

		assert(position != -1);
	
		// Update sigma_star with the shortest path trick
		sigma_star.clear();
		for (int k = 0; k < position; ++k) {
			sigma_star.push_back(Edge(sigma_prime[k]));
		}

		// Add the shortest path
		assert(position != -1);
		assert(start != -1);
		assert(finish != -1);

		const vector<int> path = graph -> dijkstra_path[start][finish];

		for (int k = 1; k < path.size(); ++k) {
			const int u = path[k - 1];
			const int v = path[k];
			for (int t = i + 1; t < m; ++t) {
				if (!mask[t]) {
					if ((u == edges[t].first) && (v == edges[t].second)) {
						sigma_star.push_back(Edge(edges[t]));
						mask[t] = true;
						break;
					}
					if ((u == edges[t].second) && (v == edges[t].first)) {
                        sigma_star.push_back(Edge(edges[t]));
						mask[t] = true;
                        break;
                    }
				}
			}
		}

		// Add the rest
		for (int k = position; k < sigma_prime.size(); ++k) {
			sigma_star.push_back(Edge(sigma_prime[k]));
		}

		// Dynamic programming
		pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sigma_star);
		best = dp.first[0][0];
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


// +---------------------------+
// | 2-MOVE-OPT (new proposal) |
// +---------------------------+

// For single-thread
pair< vector<Edge>, double> Method_2_MOVE_OPT(Graph *graph, const vector<Edge> sigma) {
    double best = INF;
    vector<Edge> result;

	// Search
	for (int i = 0; i < sigma.size(); ++i) {
		for (int j = 0; j < sigma.size(); ++j) {
			if (i == j) {
				continue;
			}

			// The rest elements except i-th and j-th
			vector<Edge> A;
			A.clear();
			for (int k = 0; k < sigma.size(); ++k) {
				if ((k != i) && (k != j)) {
					A.push_back(Edge(sigma[k]));
				}
			}
			assert(A.size() == sigma.size() - 2);

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
				candidate.push_back(Edge(sigma[i]));
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

			assert(B.size() == sigma.size() - 1);

			// Search for the best place to put the j-th in
			vector<Edge> C;
			C.clear();
			double C_value = INF;

			for (int k = 0; k < B.size(); ++k) {
				vector<Edge> candidate;
				candidate.clear();
				for (int t = 0; t < k; ++t) {
					candidate.push_back(Edge(B[t]));
				}
				candidate.push_back(Edge(sigma[j]));
				for (int t = k; t < B.size(); ++t) {
					candidate.push_back(Edge(B[t]));
				}

				// Update
				pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, candidate);
				const double cost = dp.first[0][0];
				if (cost < C_value) {
					C_value = cost;
					C = candidate;
				}
			}

			assert(C.size() == sigma.size());

			// Update
			if (C_value < best) {
				best = C_value;
				result = C;
			}
		}
	}

    return make_pair(result, best);
}


// +-----------------------+
// | k-MOVE (new proposal) |
// +-----------------------+

// For single-thread
pair< vector<Edge>, double> Method_k_MOVE_OPT(Graph *graph, const vector<Edge> sigma, const vector<int> k_indices) {
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


// +-------+
// | 1-OPT |
// +-------+

// For single-thread
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

// For single-thread
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

// For single-thread
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

// Single-thread implementation
pair< vector<Edge>, double> Iterated_Local_Search(Graph *graph, const int k_max = 75, const bool verbose = false) {
	// Greedy constructive heuristic
	pair< vector<Edge>, double > greedy = Greedy_Constructive_Heuristic(graph);
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

		if (verbose) {
			cout << "Done " << k << " iterations." << endl;
		}
	}

	return make_pair(sigma_star, best);
}


// +-----------------------------------------------------------+
// | New proposal - Iterated Local Search (ILS) Meta-heuristic |
// +-----------------------------------------------------------+

// Single-thread implementation
pair< vector<Edge>, double> Iterated_Local_Search_2(Graph *graph, const int k_max = 75, const bool verbose = false) {
    // New Greedy Constructive Heuristic
    pair< vector<Edge>, double > greedy = Greedy_Constructive_Heuristic_2(graph);
    vector<Edge> sigma_star = greedy.first;
    double best = greedy.second;

    // Iterative
    for (int k = 1; k <= k_max; ++k) {
        // Random exchange
        vector<Edge> sigma = random_exchange(sigma_star);

        // 2-MOVE
        pair< vector<Edge>, double> Result_2_MOVE_OPT = Method_2_MOVE_OPT(graph, sigma);

        // Update
        if (Result_2_MOVE_OPT.second < best) {
            best = Result_2_MOVE_OPT.second;
            sigma_star = Result_2_MOVE_OPT.first;
        }

        if (verbose) {
            cout << "Done " << k << " iterations." << endl;
        }
    }

    return make_pair(sigma_star, best);
}


// +---------------------------------------------------+
// | Variable Neighborhood Search (VNS) Meta-heuristic |
// +---------------------------------------------------+
pair< vector<Edge>, double> Variable_Neighborhood_Search(Graph *graph, const int k_max = 75, const bool verbose = false) {
    // Greedy constructive heuristic
    pair< vector<Edge>, double > greedy = Greedy_Constructive_Heuristic(graph);
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

		if (verbose) {
            cout << "Done " << k << " iterations." << endl;
        }
	}

    return make_pair(sigma_star, best);
}

