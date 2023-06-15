# Script for evaluation some test cases of load-dependent Chinese postman problem
# Author: Dr. Truong Son Hy
# Copyright 2023

rm -rf test_greedy test_ils test_vns test_ea
g++ test_greedy.cpp -o test_greedy
g++ test_ils.cpp -o test_ils
g++ test_vns.cpp -o test_vns
g++ test_ea.cpp -o test_ea

prefix=large
for idx in 48
do
./test_greedy data/${prefix}_${idx}.txt
./test_ils data/${prefix}_${idx}.txt
./test_vns data/${prefix}_${idx}.txt
./test_ea data/${prefix}_${idx}.txt
echo "-------------------------------------------------------"
done

echo "Done"
rm -rf test_greedy test_ils test_vns test_ea
