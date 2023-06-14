# Script for evaluation some test cases of load-dependent Chinese postman problem
# Author: Dr. Truong Son Hy
# Copyright 2023

rm -rf test_greedy test_ils test_vns test_ea
g++ test_greedy.cpp -o test_greedy
g++ test_ils.cpp -o test_ils
g++ test_vns.cpp -o test_vns
g++ test_ea.cpp -o test_ea

prefix=E
for idx in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
do
./test_greedy data/${prefix}_${idx}.txt
./test_ils data/${prefix}_${idx}.txt
./test_vns data/${prefix}_${idx}.txt
./test_ea data/${prefix}_${idx}.txt
echo "-------------------------------------------------------"
done

echo "Done all baselines"
rm -rf test_greedy test_ils test_vns test_ea
