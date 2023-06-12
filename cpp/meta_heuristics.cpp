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

using namespace std;

// Dynamic Programming
pair< vector< vector<double> >, vector<int> > dynamic_programming(Graph *graph, const vector<Edge> sequence) {
	// Start node
	const int start_node = graph -> start_node;

	// Check if the graph has computed Floyd's algorithm before
    assert(graph -> shortest_path.size() == graph -> num_nodes);
	const double W = graph -> W;

	// Number of edges in this list
	const int m = sequence.size();

	// Constant
    const double INF = 1e9;

	// Initialize
    vector< vector<double> > f;
    f.clear();
    for (int k = 0; k < graph -> num_nodes; ++k) {
        vector<double> vect;
        vect.clear();
        vect.push_back(INF);
        vect.push_back(INF);
        f.push_back(vect);
    }

    vector<int> direction;
    direction.clear();
    for (int k = 0; k < graph -> num_nodes; ++k) {
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
    for (int k = 0; k < graph -> num_nodes; ++k) {
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

// Compute cost
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

