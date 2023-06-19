// Generate Rural Postman Problem (RPP) cases for the load-dependent Chinese postman problem
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

using namespace std;

// Hyper-parameters
const string directory = "../data/";
const int MAX_X = 10;
const int MAX_Y = 10;
const int MAX_Q = MAX_X * MAX_Y;

// Computing the Euclidean distance
double Euclidean_distance(const pair<double, double> p1, const pair<double, double> p2) {
	const double dx = p1.first - p2.first;
	const double dy = p1.second - p2.second;
	return sqrt(dx * dx + dy * dy);
}

// Generate data
void generate(
	const string file_name, 
	const int num_nodes, 
	const int num_edges, 
	const double W_const, 
	const bool random_q,
	const int factor = 2 // Factor to determine if we need to deliver over an edge or not
){
	vector< pair<double, double> > position;
	position.clear();

	for (int i = 0; i < num_nodes; ++i) {
		while (true) {
			const int x = rand() % (MAX_X + 1);
			const int y = rand() % (MAX_Y + 1);
			bool found = false;
			for (int j = 0; j < i; ++j) {
				if ((position[j].first == x) && (position[j].second == y)) {
					found = true;
					break;
				}
			}
			if (!found) {
				position.push_back(make_pair(x, y));
				break;
			}
		}
	}

	vector< vector<double> > distance;
	distance.clear();
	for (int i = 0; i < num_nodes; ++i) {
		vector<double> vect;
		vect.clear();
		for (int j = 0; j < num_nodes; ++j) {
			vect.push_back(Euclidean_distance(position[i], position[j]));
		}
		distance.push_back(vect);
	}

	vector<int> path;
	path.clear();
	for (int i = 0; i < num_nodes; ++i) {
		path.push_back(i);
	}

	for (int i = 0; i < num_nodes; ++i) {
		const int j = rand() % num_nodes;
		swap(path[i], path[j]);
	}
	
	vector< pair<int, int> > edges;
	edges.clear();
	
	for (int i = 0; i < num_nodes - 1; ++i) {
		edges.push_back(make_pair(path[i], path[i + 1]));
	}

	for (int i = 0; i < num_edges - num_nodes + 1; ++i) {
		while (true) {
			const int first = rand() % num_nodes;
			const int second = rand() % num_nodes;
			if (first == second) {
				continue;
			}
			bool found = false;
			for (int k = 0; k < edges.size(); ++k) {
				if ((first == edges[k].first) && (second == edges[k].second)) {
					found = true;
					break;
				}
				if ((second == edges[k].first) && (first == edges[k].second)) {
                    found = true;
                    break;
                }
			}
			if (!found) {
				edges.push_back(make_pair(first, second));
				break;
			}
		}
	}

	assert(edges.size() == num_edges);

	vector<double> D;
	vector<double> Q;
	D.clear();
	Q.clear();

	for (int i = 0; i < num_edges; ++i) {
		const int u = edges[i].first;
		const int v = edges[i].second;
		const double d = distance[u][v];
		D.push_back(d);
		if (random_q) {
			Q.push_back(rand() % MAX_Q + 1);
		} else {
			Q.push_back(d);
		}
	}

	for (int i = 0; i < num_edges; ++i) {
		// We don't have to deliver over this edge
		if (rand() % factor == 0) {
			Q[i] = 0;
		}
	}

	double sum_Q = 0;
	for (int i = 0; i < num_edges; ++i) {
		sum_Q += Q[i];
	}
	const double W = W_const * sum_Q;

	ofstream file;
	file.open(file_name);
	file << "Number of nodes:" << endl;
	file << num_nodes << endl;
	file << "Number of edges:" << endl;
	file << num_edges << endl;
	file << "W:" << endl;
	file << W << endl;
	file << "Edges (node i, node j, d_ij, q_ij):" << endl;
	for (int i = 0; i < num_edges; ++i) {
		file << edges[i].first << " " << edges[i].second << " " << D[i] << " " << Q[i] << endl;
	}
	file << "Coordinates:" << endl;
	for (int i = 0; i < num_nodes; ++i) {
		file << position[i].first << " " << position[i].second << endl;
	}
	file.close();
}

// Generation
void generate_cases(const vector<int> V, const vector<int> E, const vector<double> W_const, const string prefix) {
	assert(V.size() == E.size());
    int count = 0;

    for (int i = 0; i < V.size(); ++i) {
        const int num_nodes = V[i];
        const int num_edges = E[i];
        for (int j = 0; j < W_const.size(); ++j) {
            count++;
            generate(directory + "/" + prefix + "_" + to_string(count) + ".txt", num_nodes, num_edges, W_const[j], false);

            count++;
            generate(directory + "/" + prefix + "_" + to_string(count) + ".txt", num_nodes, num_edges, W_const[j], true);
        }
    }
}

// Christofides et al.
void Christofides_et_al() {
    const vector<int> V = {11, 14, 17};
    const vector<int> E = {13, 32, 35};
    const vector<double> W_const = {0.0, 0.5, 5.0};
    const string prefix = "C";
	
	generate_cases(V, E, W_const, prefix);
}

// Hertz et al.
void Hertz_et_al() {
	const vector<int> V = {6, 7, 8, 9, 10, 11, 12, 14, 18, 21, 24, 27};
    const vector<int> E = {11, 15, 20, 25, 30, 35, 40, 48, 22, 26, 30, 33};
    const vector<double> W_const = {0.5};
    const string prefix = "H";
	
	generate_cases(V, E, W_const, prefix);
}

// Main program
int main(int argc, char **argv) {
	// Fix random seed
	srand(0);

	// Generate cases based on Christofides et al.
	Christofides_et_al();

	// Generate cases based on Hertz et al.
	Hertz_et_al();

	return 0;
}
