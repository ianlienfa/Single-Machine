#! /usr/local/bin/bash

input=("5_v1" "10_v1" "15_v1" "20_v1" "25_v1" "30_v1")
# input=("5_v1")
# algorithms=("BFSCBB" "localSearch" "DFS" "BFS")
algorithms=("BFS")

# call different algorithms for different instances
for algo in "${algorithms[@]}"; do
    for i in "${input[@]}"; do
        echo "Running test ${i}.in"    
        timeout --signal='SIGINT' 2m ./permSet_outtree $algo < input/${i}.in > out/${algo}/$i.out     
    done
done
