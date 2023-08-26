// Testing the graph structure for the load-dependent Chinese postman problem
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

#include "../graph_library/Graph.h"

using namespace std;

// Main program
int main(int argc, char **argv) {
	// Load graph from file
	Graph *graph = new Graph("../data/sample_input_1.txt");
	
	// View the graph's content
	cout << graph -> get_content() << endl;
	
	// Save the graph object to file
	graph -> save_to_file("graph.txt");

	// Run the Floyd's algorithm
	graph -> Floyd_algorithm();

	cout << "Floyd's algorithm:" << endl;
	for (int i = 0; i < graph -> num_nodes; ++i) {
		for (int j = 0; j < graph -> num_nodes; ++j) {
			cout << graph -> shortest_path[i][j] << " ";
		}
		cout << endl;
	}

	// Run the Dijkstra's algorihtm
	graph -> Dijkstra_algorithm();

	cout << "\nDijkstra's algorithm:" << endl;
    for (int i = 0; i < graph -> num_nodes; ++i) {
        for (int j = 0; j < graph -> num_nodes; ++j) {
            cout << graph -> dijkstra_d[i][j] << " ";
        }
        cout << endl;
    }

	// Check quality of Dijkstra vs Floyd
	for (int i = 0; i < graph -> num_nodes; ++i) {
        for (int j = 0; j < graph -> num_nodes; ++j) {
			const double diff = abs(graph -> shortest_path[i][j] - graph -> dijkstra_d[i][j]);
			assert(diff < 1e-6);
		}
	}

	// Check paths
	for (int i = 0; i < graph -> num_nodes; ++i) {
        for (int j = 0; j < graph -> num_nodes; ++j) {
			const vector<int> path = graph -> dijkstra_path[i][j];
			
			if (i == j) {
                assert(path.size() == 0);
                continue;
            }
			
			assert(path[0] == i);
			assert(path[path.size() - 1] == j);

			cout << "Shortest path from " << i << " to " << j << ": ";
			for (int k = 0; k < path.size(); ++k) {
				cout << path[k] << " ";
			}
			cout << endl;

			double sum = 0;
			for (int k = 1; k < path.size(); ++k) {
				const int u = path[k - 1];
				const int v = path[k];

				bool found = false;
				for (int e = 0; e < graph -> edges.size(); ++e) {
					Edge edge = graph -> edges[e];
					if ((u == edge.first) && (v == edge.second)) {
						found = true;
						sum += edge.d;
						break;
					}
					if ((u == edge.second) && (v == edge.first)) {
                        found = true;
                        sum += edge.d;
                        break;
                    }
				}
				assert(found == true);
			}
			assert(abs(sum - graph -> dijkstra_d[i][j]) < 1e-6);
		}
	}

	cout << "\nDone" << endl;
	return 0;
}
