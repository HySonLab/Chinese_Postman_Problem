#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include <dirent.h>
#include <sys/types.h>

using namespace std;

vector<string> get_data_names(const char *path) {
	vector<string> result;
	result.clear();
	
	DIR *dr;
	struct dirent *en;
	dr = opendir(path);
	if (dr) {
		while ((en = readdir(dr)) != NULL) {
			const string str = en -> d_name;
			if (str.length() > 4) { 
				const int n = str.length();
				if ((str[n - 1] == 't') && (str[n - 2] == 'x') && (str[n - 3] == 't') && (str[n - 4] == '.')) {
					result.push_back(str);
				}
			}
		}
		closedir(dr);
	}

	sort(result.begin(), result.end());
	return result;
}

int main(void) {
	// Input & output filenames
	const string input_dir = "../data/";
	const string output_dir = "../data/results/";
	vector<string> file_names = get_data_names(input_dir.c_str());
	cout << "Number of examples: " << file_names.size() << endl;

	// Compilation
	system("g++ test_greedy.cpp -o greedy");
	system("g++ test_greedy_2.cpp -o greedy_2");
	system("g++ test_ils_multithreads.cpp -o ils -lpthread");
	system("g++ test_ils_2.cpp -o ils_2");
	system("g++ test_vns.cpp -o vns");
	system("g++ test_de_multithreads.cpp -o de -lpthread");
	system("g++ test_ea_multithreads.cpp -o ea -lpthread");
	system("g++ test_aco_multithreads.cpp -o aco -lpthread");

	// Create new directories
	system("mkdir ../data/results/");
	system("mkdir ../data/results/greedy");
	system("mkdir ../data/results/greedy_2");
	system("mkdir ../data/results/ils");
	system("mkdir ../data/results/ils_2");
	system("mkdir ../data/results/de_k3");
	system("mkdir ../data/results/de_k4");
	system("mkdir ../data/results/de_k5");
	system("mkdir ../data/results/vns");
	system("mkdir ../data/results/ea");
	system("mkdir ../data/results/aco");

	// Run all
    for (int i = 0; i < file_names.size(); ++i) {
        const string input_fn = input_dir + file_names[i];
        
		cout << "------------------------------------------------" << endl;
		cout << "Example " << (i + 1) << " / " << file_names.size() << ":" << endl;
		cout << "Input: " << input_fn << endl;

		string command;
		
		/*
		// Greedy
		command = "./greedy " + input_fn + " > ../data/results/greedy/" + file_names[i];
		system(command.c_str());
		cout << "Done Greedy algorithm" << endl;

		// ILS 
		command = "./ils " + input_fn + " > ../data/results/ils/" + file_names[i];
        system(command.c_str());
        cout << "Done ILS algorithm" << endl;

		// VNS
        command = "./vns " + input_fn + " > ../data/results/vns/" + file_names[i];
        system(command.c_str());
        cout << "Done VNS algorithm" << endl;

		// EA
        command = "./ea " + input_fn + " > ../data/results/ea/" + file_names[i];
        system(command.c_str());
        cout << "Done EA algorithm" << endl;

    	// ACO
        command = "./aco " + input_fn + " > ../data/results/aco/" + file_names[i];
        system(command.c_str());
        cout << "Done ACO algorithm" << endl;

		// New proposal for Greedy
        command = "./greedy_2 " + input_fn + " > ../data/results/greedy_2/" + file_names[i];
        system(command.c_str());
        cout << "Done Greedy_2 algorithm" << endl;

		// New proposal for ILS
        command = "./ils_2 " + input_fn + " > ../data/results/ils_2/" + file_names[i];
        system(command.c_str());
        cout << "Done ILS_2 algorithm" << endl;
		*/

		// DE (k = 3)
		command = "./de " + input_fn + " 3 > ../data/results/de_k3/" + file_names[i];
        system(command.c_str());
        cout << "Done DE algorithm (k = 3)" << endl;

		// DE (k = 4)
        command = "./de " + input_fn + " 4 > ../data/results/de_k4/" + file_names[i];
        system(command.c_str());
        cout << "Done DE algorithm (k = 4)" << endl;

		// DE (k = 5)
        command = "./de " + input_fn + " 5 > ../data/results/de_k5/" + file_names[i];
        system(command.c_str());
        cout << "Done DE algorithm (k = 5)" << endl;
	}

	return(0);
}
