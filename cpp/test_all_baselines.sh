rm -rf test_greedy test_ils test_vns test_brute_force
g++ test_greedy.cpp -o test_greedy
g++ test_ils.cpp -o test_ils
g++ test_vns.cpp -o test_vns
g++ test_brute_force.cpp -o test_brute_force

for data_fn in data/sample_input_1.txt data/sample_input_2.txt data/small_1.txt data/small_2.txt data/small_3.txt data/small_4.txt data/small_5.txt data/small_6.txt
do
./test_greedy $data_fn
./test_ils $data_fn
./test_vns $data_fn
./test_brute_force $data_fn
echo "-------------------------------------------------------"
done
