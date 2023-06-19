# Script for evaluation some test cases of load-dependent Chinese postman problem
# Author: Dr. Truong Son Hy
# Copyright 2023

rm -rf test_greedy test_ils_multithreads test_vns test_ea_multithreads test_aco_multithreads
g++ test_greedy.cpp -o test_greedy
g++ test_ils_multithreads.cpp -o test_ils_multithreads -lpthread
g++ test_vns.cpp -o test_vns
g++ test_ea_multithreads.cpp -o test_ea_multithreads -lpthread
g++ test_aco_multithreads.cpp -o test_aco_multithreads -lpthread

prefix=C
for idx in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
do
./test_greedy data/${prefix}_${idx}.txt
./test_ils_multithreads data/${prefix}_${idx}.txt
./test_vns data/${prefix}_${idx}.txt
./test_ea_multithreads data/${prefix}_${idx}.txt
./test_aco_multithreads data/${prefix}_${idx}.txt
echo "-------------------------------------------------------"
done

prefix=H
for idx in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
do
./test_greedy data/${prefix}_${idx}.txt
./test_ils_multithreads data/${prefix}_${idx}.txt
./test_vns data/${prefix}_${idx}.txt
./test_ea_multithreads data/${prefix}_${idx}.txt
./test_aco_multithreads data/${prefix}_${idx}.txt
echo "-------------------------------------------------------"
done

prefix=E
for idx in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
do
./test_greedy data/${prefix}_${idx}.txt
./test_ils_multithreads data/${prefix}_${idx}.txt
./test_vns data/${prefix}_${idx}.txt
./test_ea_multithreads data/${prefix}_${idx}.txt
./test_aco_multithreads data/${prefix}_${idx}.txt
echo "-------------------------------------------------------"
done

echo "Done all baselines"
rm -rf test_greedy test_ils_multithreads test_vns test_ea_multithreads test_aco_multithreads
