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

int main(int argc, char **argv) {
	// Load graph from file
	Graph *graph = new Graph("data/E_10.txt");
	
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
	return 0;
}
