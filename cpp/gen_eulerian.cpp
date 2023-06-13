// Generate Eulerian cases for the load-dependent Chinese postman problem
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
const string directory = "data/";
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
	const bool random_q
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

	const int start_node = path[0];
	
	vector< pair<int, int> > edges;
	edges.clear();
	
	for (int i = 0; i < num_nodes - 1; ++i) {
		edges.push_back(make_pair(path[i], path[i + 1]));
	}

	for (int i = 0; i < num_edges - num_nodes; ++i) {
		const int last_node = path[path.size() - 1];
		while (true) {
			const int next_node = rand() % num_nodes;
			bool found = false;
			if (next_node == last_node) {
				continue;
			}
			if (i == num_edges - num_nodes - 1) {
				if (next_node == start_node) {
					continue;
				}
			}
			for (int j = 0; j < edges.size(); ++j) {
				if ((last_node == edges[j].first) && (next_node == edges[j].second)) {
					found = true;
					break;
				}
				if ((last_node == edges[j].second) && (next_node == edges[j].first)) {
                    found = true;
                    break;
                }
			}
			if (!found) {
				path.push_back(next_node);
				edges.push_back(make_pair(last_node, next_node));
				break;
			}
		}
	}

	const int last_node = path[path.size() - 1];
	assert(last_node != start_node);

	path.push_back(start_node);
	edges.push_back(make_pair(last_node, start_node));

	assert(path.size() == num_edges + 1);
	assert(edges.size() == num_edges);

	vector<double> D;
	vector<double> Q;
	D.clear();
	Q.clear();
	double sum_Q = 0;

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

// Main program
int main(int argc, char **argv) {
	// Fix random seed
	srand(0);

	const vector<int> V = {20};
	const vector<int> E = {75};
	const vector<double> W_const = {0.0, 0.5, 5.0};

	assert(V.size() == E.size());
	int count = 0;

	for (int i = 0; i < V.size(); ++i) {
		const int num_nodes = V[i];
		const int num_edges = E[i];
		for (int j = 0; j < W_const.size(); ++j) {
			count++;
			generate(directory + "/middle_" + to_string(count) + ".txt", num_nodes, num_edges, W_const[j], false);

			count++;
			generate(directory + "/middle_" + to_string(count) + ".txt", num_nodes, num_edges, W_const[j], true);
		}
	}
	return 0;
}
