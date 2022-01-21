#! /usr/local/bin/bash
directories=("BFS" "BFSCBB" "DFS" "localSearch")
# echo $directoies
for dir in "${directories[@]}"; do
    outfiles=(`ls ../out/$dir`)
    for out in "${outfiles[@]}"; do
        touch "$dir/$out"
        cat "../out/$dir/$out" | grep min | awk '{print $3}' >> "$dir/$out"        
    done
done
