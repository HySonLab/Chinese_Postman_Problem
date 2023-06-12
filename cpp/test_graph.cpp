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
	Graph *graph = new Graph("data/E_10.txt");
	cout << graph -> get_content() << endl;
	graph -> save_to_file("test.txt");
	return 0;
}
