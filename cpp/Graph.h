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

struct Coordinate {
	Coordinate(const double x_, const double y_) {
		x = x_;
		y = y_;
	}

	double x;
	double y;
};

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

class Graph {
public:
	Graph(){
	}

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

	Graph(const string file_name) {
		load_from_file(file_name);
	}

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

	void save_to_file(const string file_name) {
		ofstream file;
		file.open(file_name);
		file << get_content();
		file.close();
	}

	int num_nodes;
	int num_edges;
	double W;
	vector<Edge> edges;
	vector<Coordinate> position;
};

#endif
