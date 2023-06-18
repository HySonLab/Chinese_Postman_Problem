# Script for evaluation some test cases of load-dependent Chinese postman problem
# Author: Dr. Truong Son Hy
# Copyright 2023

rm -rf test_greedy test_ils_multithreads test_vns test_ea_multithreads test_brute_force_multithreads
g++ test_greedy.cpp -o test_greedy
g++ test_ils_multithreads.cpp -o test_ils_multithreads -lpthread
g++ test_vns.cpp -o test_vns
g++ test_ea_multithreads.cpp -o test_ea_multithreads -lpthread
g++ test_brute_force_multithreads.cpp -o test_brute_force_multithreads -lpthread

for data_fn in data/sample_input_1.txt data/sample_input_2.txt
do
./test_greedy $data_fn
./test_ils_multithreads $data_fn
./test_vns $data_fn
./test_ea_multithreads $data_fn
./test_brute_force_multithreads $data_fn
echo "-------------------------------------------------------"
done

for idx in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
do
./test_greedy data/small_$idx.txt
./test_ils_multithreads data/small_$idx.txt
./test_vns data/small_$idx.txt
./test_ea_multithreads data/small_$idx.txt
./test_brute_force_multithreads data/small_$idx.txt
echo "-------------------------------------------------------"
done

echo "Done all baselines"
rm -rf test_greedy test_ils_multithreads test_vns test_ea_multithreads test_brute_force_multithreads
