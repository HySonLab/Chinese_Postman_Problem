// Evolutionary Algorithm (EA) for the load-dependent Chinese postman problem
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

#include "../graph_library/Graph.h"
#include "../graph_library/meta_heuristics.h"

using namespace std;

// +------------------------+
// | Evolutionary Algorithm |
// +------------------------+

// Mixing two parents into a child
pair< vector<Edge>, double > mix_parents(Graph *graph, const vector<Edge> father, const vector<Edge> mother) {
	// Information
	const int num_nodes = graph -> num_nodes;
	const int m = graph -> num_deliver_edges;

	assert(m == father.size());
	assert(m == mother.size());

	// Boolean mask
	vector< vector<bool> > mask;
	mask.clear();
	for (int i = 0; i < num_nodes; ++i) {
		vector<bool> vect;
		vect.clear();
		for (int j = 0; j < num_nodes; ++j) {
			vect.push_back(false);
		}
		mask.push_back(vect);
	}
		
	// The child
	vector<Edge> child;
	child.clear();

	// Mixing
	int i = 0;
	int j = 0;
	while ((child.size() < m) && (i < m) && (j < m)) {
		// Random selection
		int choice = rand() % 2;

		// If the father has no more gene, select from the mother's gene
		if (i == m) {
			choice = 1;
		}

		// If the mother has no more gene, select from the father's gene
		if (j == m) {
			choice = 0;
		}

		// Selection
		if (choice == 1) {
			assert(j < m);
			const int x = mother[j].first;
			const int y = mother[j].second;

			// If the edge was not visited then we add to the child's gene
			if ((!mask[x][y]) && (!mask[y][x])) {
				mask[x][y] = true;
				mask[y][x] = true;
				child.push_back(Edge(mother[j]));
			}
			++j;
		} else {
			assert(i < m);
			const int x = father[i].first;
            const int y = father[i].second;

			// If the edge was not visited then we add to the child's gene
            if ((!mask[x][y]) && (!mask[y][x])) {
                mask[x][y] = true;
                mask[y][x] = true;
                child.push_back(Edge(father[i]));
            }
            ++i;
		}
	}

	// If the child doesn't have enough gene, get the rest from the father
	i = 0;
	while ((child.size() < m) && (i < m)) {
		const int x = father[i].first;
        const int y = father[i].second;

        // If the edge was not visited then we add to the child's gene
        if ((!mask[x][y]) && (!mask[y][x])) {
        	mask[x][y] = true;
        	mask[y][x] = true;
        	child.push_back(Edge(father[i]));
		}
		++i;
	}

	// If the child doesn't have enough gene, get the rest from the mother
	j = 0;
	while ((child.size() < m) && (j < m)) {
		const int x = mother[j].first;
        const int y = mother[j].second;

        // If the edge was not visited then we add to the child's gene
        if ((!mask[x][y]) && (!mask[y][x])) {
        	mask[x][y] = true;
        	mask[y][x] = true;
        	child.push_back(Edge(mother[j]));
       	}
       	++j;
	}

	/*
	if (child.size() < m) {
		cout << child.size() << " " << m << endl;
		cout << "Father:" << endl;
		for (int k = 0; k < father.size(); ++k) {
			cout << father[k].first << " " << father[k].second << endl;
		}
		cout << "Mother:" << endl;
        for (int k = 0; k < mother.size(); ++k) {
            cout << mother[k].first << " " << mother[k].second << endl;
        }   
		cout << "Child:" << endl;
        for (int k = 0; k < child.size(); ++k) {
            cout << child[k].first << " " << child[k].second << endl;
        }   
	}
	*/

	assert(child.size() == m);

	// Compute cost
	pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, child);
	const double cost = dp.first[0][0];

	return make_pair(child, cost);
}

// Evolutionary Algorithm
pair< vector<Edge>, double> Evolutionary_Algorithm(
	Graph *graph, 
	const int k_max = 75, 
	const int max_population = 10, 
	const bool verbose = true,
	const bool use_1_OPT = true,
	const bool use_2_OPT = true,
	const bool use_2_EXCHANGE = true
) {
    // Greedy constructive heuristic
    pair< vector<Edge>, double > greedy = Greedy_Constructive_Heuristic(graph);

	// Number of deliver edges
	const int m = greedy.first.size();
	assert(m == graph -> num_deliver_edges);

	// Initialize the population
	vector< pair< vector<Edge>, double > > population;
	population.clear();

	population.push_back(greedy);

	while (population.size() < max_population) {
		// Random exchange
		vector<Edge> sequence = random_exchange(greedy.first);
		pair< vector< vector<double> >, vector<int> > dp = dynamic_programming(graph, sequence);
		population.push_back(make_pair(sequence, dp.first[0][0]));
	}

	if (verbose) {
		cout << "Initial population: " << population.size() << endl;
	}

	// Evolution
	for (int k = 1; k <= k_max; ++k) {
		vector< pair< vector<Edge>, double > > children;
    	children.clear();

		// Create new children via mixing
		for (int i = 0; i < population.size(); ++i) {
			const int j = rand() % population.size();
			pair< vector<Edge>, double > child = mix_parents(graph, population[i].first, population[j].first);
			children.push_back(child);
		}

		// Create new children via mutation
		for (int i = 0; i < population.size(); ++i) {
			if (use_1_OPT) {
				children.push_back(Method_1_OPT(graph, population[i].first));
			}
			if (use_2_OPT) {
				children.push_back(Method_2_OPT(graph, population[i].first));
			}
			if (use_2_EXCHANGE) {
				children.push_back(Method_2_EXCHANGE(graph, population[i].first));
			}
		}

		// Construct the ranking via fitness
		vector< pair<double, int> > rank;
		rank.clear();
		for (int i = 0; i < population.size(); ++i) {
			rank.push_back(make_pair(population[i].second, i));
		}
		for (int i = 0; i < children.size(); ++i) {
			rank.push_back(make_pair(children[i].second, i + population.size()));
		}
		sort(rank.begin(), rank.end());

		// Select the best for the next generation
		// assert(rank.size() >= max_population);
		
		vector< pair< vector<Edge>, double > > next_generation;
		next_generation.clear();

		for (int i = 0; i < min(max_population, int(rank.size())); ++i) {
			const int index = rank[i].second;
			if (index < population.size()) {
				next_generation.push_back(population[index]);
			} else {
				next_generation.push_back(children[index - population.size()]);
			}
		}

		population = next_generation;

		// Prune the population, only keep the unique
		vector<bool> mask;
		mask.clear();
		for (int i = 0; i < population.size(); ++i) {
			mask.push_back(false);
		}

		for (int i = 0; i < population.size(); ++i) {
			if (mask[i]) {
				continue;
			}
			for (int j = i + 1; j < population.size(); ++j) {
				if (mask[j]) {
					continue;
				}

				// Check if these twos are the same
				if (abs(population[i].second - population[j].second) > 1e-3) {
					continue;
				}
				bool same = true;
				for (int k = 0; k < m; ++k) {
					const Edge e1 = population[i].first[k];
					const Edge e2 = population[j].first[k];
					if ((e1.first == e2.first) && (e1.second == e2.second)) {
						continue;
					}
					if ((e1.first == e2.second) && (e1.second == e2.first)) {
						continue;
					}
					same = false;
					break;
				}
				if (same) {
					mask[j] = true;
				}
			}
		}

		vector< pair< vector<Edge>, double > > unique;
		unique.clear();
		for (int i = 0; i < population.size(); ++i) {
			if (!mask[i]) {
				unique.push_back(population[i]);
			}
		}
		population = unique;

		if (verbose) {
            cout << "After " << k << " iterations: Population = " << population.size() << endl;
        }
	}

	// Select the best gene
	vector<Edge> sigma_star;
	double best = INF;

	for (int i = 0; i < population.size(); ++i) {
		if (population[i].second < best) {
			best = population[i].second;
			sigma_star = population[i].first;
		}
	}

	return make_pair(sigma_star, best);
}

