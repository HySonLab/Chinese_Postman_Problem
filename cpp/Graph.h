#ifndef __GRAPH_H_INCLUDED__
#define __GRAPH_H_INCLUDED__

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

// Edge structure
struct Edge {
	Edge(const int first_, const int second_, const double d_, const double q_) {
		first = first_;
		second = second_;
		d = d_;
		q = q_;
	}

	int first;
	int second;
	double d;
	double q;
};

// Coordinate structure
struct Coordinate {
	Coordinate(const double x_, const double y_) {
		x = x_;
		y = y_;
	}

	double x;
	double y;
};

// Convert a string into a vector of strings
vector<string> parse_line(const string str, const char separator = ' ') {
	vector<string> result;
	result.clear();
	int i = 0;
	while (i < str.length()) {
		if (str[i] == separator) {
			++i;
			continue;
		}
		string word = "";
		for (int j = i; j < str.length(); ++j) {
			if (str[j] == separator) {
				break;
			} else {
				word += str[j];
			}
		}
		result.push_back(word);
		i += word.length();
	}
	return result;
}

// Convert a matrix to a string
string to_str(const vector< vector<double> > matrix) {
    string result = "";
    for (int i = 0; i < matrix.size(); ++i) {
        for (int j = 0; j < matrix[i].size(); ++j) {
            result += to_string(matrix[i][j]) + " ";
        }
        if (i + 1 < matrix.size()) {
            result += "\n";
        }
    }
    return result;
}

// Convert a vector to a string
string to_str(const vector<int> vect) {
    string result = "";
    for (int i = 0; i < vect.size(); ++i) {
        result += to_string(vect[i]) + " ";
    }
    return result;
}

// Graph structure
class Graph {
public:
	// Constructor
	Graph(){
	}

	// Constructor given another object
	Graph(Graph *another) {
		this -> num_nodes = another -> num_nodes;
		this -> num_edges = another -> num_edges;
		this -> W = another -> W;

		this -> edges.clear();
		for (int i = 0; i < another -> edges.size(); ++i) {
			this -> edges.push_back(another -> edges[i]);
		}

		this -> position.clear();
		for (int i = 0; i < another -> position.size(); ++i) {
			this -> position.push_back(another -> position[i]);
		}
	}

	// Constructor given the content
	Graph(const int num_nodes, const int num_edges, const double W, const vector<Edge> edges, const vector<Coordinate> position) {
        this -> num_nodes = num_nodes;
        this -> num_edges = num_edges;
		this -> W = W;

        this -> edges.clear();
        for (int i = 0; i < edges.size(); ++i) {
            this -> edges.push_back(edges[i]);
        }

        this -> position.clear();
        for (int i = 0; i < position.size(); ++i) {
            this -> position.push_back(position[i]);
        }
    }

	// Load the content from file
	void load_from_file(const string file_name) {
		ifstream file;
        file.open(file_name);

        string str;
        getline(file, str);
        getline(file, str);
        this -> num_nodes = stoi(str);

        getline(file, str);
        getline(file, str);
        this -> num_edges = stoi(str);

        getline(file, str);
        getline(file, str);
        this -> W = stof(str);

        getline(file, str);
        this -> edges.clear();
        for (int i = 0; i < this -> num_edges; ++i) {
            getline(file, str);
            vector<string> words = parse_line(str);
            assert(words.size() == 4);
            const int first = stoi(words[0]);
            const int second = stoi(words[1]);
            const double d = stof(words[2]);
            const double q = stof(words[3]);
            this -> edges.push_back(Edge(first, second, d, q));
        }

        getline(file, str);
        this -> position.clear();
        for (int i = 0; i < this -> num_nodes; ++i) {
            getline(file, str);
            vector<string> words = parse_line(str);
            assert(words.size() == 2);
            const double x = stof(words[0]);
            const double y = stof(words[1]);
            this -> position.push_back(Coordinate(x, y));
        }

        file.close();
	}

	// Constructor given a file
	Graph(const string file_name) {
		load_from_file(file_name);
	}

	// Convert the content to a string
	string get_content() {
		string result = "";
		result += "Number of nodes:\n";
		result += to_string(this -> num_nodes) + "\n";
		result += "Number of edges:\n";
		result += to_string(this -> num_edges) + "\n";
		result += "W:\n";
		result += to_string(this -> W) + "\n";
		result += "Edges (node i, node j, d_ij, q_ij):\n";
		for (int i = 0; i < this -> edges.size(); ++i) {
			result += to_string(this -> edges[i].first) + " " + to_string(this -> edges[i].second) + " " + to_string(this -> edges[i].d) + " " + to_string(this -> edges[i].q) + "\n";
		}
		result += "Coordinates:\n";
		for (int i = 0; i < this -> position.size(); ++i) {
			result += to_string(this -> position[i].x) + " " + to_string(this -> position[i].y) + "\n";
		}
		return result;
	}

	// Save the object to file
	void save_to_file(const string file_name) {
		ofstream file;
		file.open(file_name);
		file << get_content();
		file.close();
	}

	// Search edge
	Edge get_edge(const int i, const int j) {
		for (int k = 0; k < this -> edges.size(); ++k) {
			const int first = this -> edges[k].first;
			const int second = this -> edges[k].second;
			const double d = this -> edges[k].d;
			const double q = this -> edges[k].q;
			if ((first == i) && (second == j)) {
				return Edge(i, j, d, q);
			}
			if ((first == j) && (second == i)) {
				return Edge(i, j, d, q);
			}
		}
		return Edge(-1, -1, -1, -1);
	}

	// Floyd's algorithm to find all-pair shortest paths
	void Floyd_algorithm() {
		const double INF = 1e9;
		
		this -> shortest_path.clear();
		for (int i = 0; i < this -> num_nodes; ++i) {
			vector<double> vect;
			vect.clear();
			for (int j = 0; j < this -> num_nodes; ++j) {
				vect.push_back(INF);
			}
			this -> shortest_path.push_back(vect);
		}

		for (int i = 0; i < this -> num_nodes; ++i) {
			this -> shortest_path[i][i] = 0;
		}

		for (int i = 0; i < this -> num_edges; ++i) {
			const int first = this -> edges[i].first;
			const int second = this -> edges[i].second;
			const double d = this -> edges[i].d;

			if (d < this -> shortest_path[first][second]) {
				this -> shortest_path[first][second] = d;
			}
			
			if (d < this -> shortest_path[second][first]) {
                this -> shortest_path[second][first] = d;
            }
		}

		for (int k = 0; k < this -> num_nodes; ++k) {
			for (int i = 0; i < this -> num_nodes; ++i) {
				for (int j = 0; j < this -> num_nodes; ++j) {
					if (this -> shortest_path[i][j] > this -> shortest_path[i][k] + this -> shortest_path[k][j]) {
						this -> shortest_path[i][j] = this -> shortest_path[i][k] + this -> shortest_path[k][j];
					}
				}
			}
		}
	}

	// Staring node is always indexed as 0
	const int start_node = 0;

	// Number of nodes
	int num_nodes;

	// Number of edges
	int num_edges;

	// W
	double W;

	// Edges
	vector<Edge> edges;

	// Position of each node
	vector<Coordinate> position;

	// All-pair shortest paths (from Floyd's algorithm)
	vector< vector<double> > shortest_path;
};

#endif
