# Import OS module
import os

path = "de_k3"
output_fn = path + ".csv"

output = open(output_fn, 'w')
output.write("Name,Number of nodes,Number of edges,Number of deliver edges,Cost,Time (ms),Time (s)\n")

dir_list = os.listdir(path)
dir_list.sort()

dir_list.remove("aaai_example_1.txt")
dir_list.remove("aaai_example_2.txt")
dir_list.remove("sample_input_1.txt")
dir_list.remove("sample_input_2.txt")

def overlap(s1, s2):
	if len(s1) <= len(s2):
		return False
	return s1[:len(s2)] == s2

for element in dir_list:
	input_fn = path + "/" + element
	with open(input_fn, 'r') as f:
		text = f.read()
	text = text.split('\n')
	
	# Extract information
	name = element
	num_nodes = -1
	num_edges = -1
	num_deliver_edges = -1
	cost = -1
	time_ms = -1
	time_s = -1

	for s in text:
		if overlap(s, "Number of nodes:"):
			assert num_nodes == -1
			s_ = s.split()
			num_nodes = int(s_[len(s_) - 1])
		
		if overlap(s, "Number of edges:"):
			assert num_edges == -1
			s_ = s.split()
			num_edges = int(s_[len(s_) - 1])

		if overlap(s, "Number of deliver"):
			assert num_deliver_edges == -1
			s_ = s.split()
			num_deliver_edges = int(s_[len(s_) - 1])

		if overlap(s, "Cost"):
			assert cost == -1
			s_ = s.split()
			cost = float(s_[len(s_) - 1])

		if overlap(s, "Running time (milliseconds):"):
			assert time_ms == -1
			s_ = s.split()
			time_ms = int(s_[len(s_) - 1])

		if overlap(s, "Running time (seconds):"):
			assert time_s == -1
			s_ = s.split()
			time_s = int(s_[len(s_) - 1])

	assert num_nodes != -1
	assert num_edges != -1
	assert num_deliver_edges != -1
	assert cost != -1
	assert time_ms != -1
	assert time_s != -1
	print(name, num_nodes, num_edges, num_deliver_edges, cost, time_ms, time_s)
	output.write(name + "," + str(num_nodes) + "," + str(num_edges) + "," + str(num_deliver_edges) + "," + str(cost) + "," + str(time_ms) + "," + str(time_s) + "\n")

output.close()
print("Number of files:", len(dir_list))
print("Done")

