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

string to_str(const vector<int> vect) {
    string result = "";
    for (int i = 0; i < vect.size(); ++i) {
        result += to_string(vect[i]) + " ";
    }
    return result;
}

