# Script for evaluation some test cases of load-dependent Chinese postman problem
# Author: Dr. Truong Son Hy
# Copyright 2023

rm -rf test_greedy test_ils test_vns test_ea test_brute_force
g++ test_greedy.cpp -o test_greedy
g++ test_ils.cpp -o test_ils
g++ test_vns.cpp -o test_vns
g++ test_ea.cpp -o test_ea
g++ test_brute_force.cpp -o test_brute_force

for data_fn in data/sample_input_1.txt data/sample_input_2.txt
do
./test_greedy $data_fn
./test_ils $data_fn
./test_vns $data_fn
./test_ea $data_fn
./test_brute_force $data_fn
echo "-------------------------------------------------------"
done

for idx in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
do
./test_greedy data/small_$idx.txt
./test_ils data/small_$idx.txt
./test_vns data/small_$idx.txt
./test_ea data/small_$idx.txt
./test_brute_force data/small_$idx.txt
echo "-------------------------------------------------------"
done

echo "Done all baselines"
rm -rf test_greedy test_ils test_vns test_ea test_brute_force
