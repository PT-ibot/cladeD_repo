#!/bin/bash
mkdir phyparts
n=1
cat possible_scenarios.tre | while read -r a
do
	echo ${a} > tree_mod_$n.tre
	for i in {1..50}
	do
		java -jar /path/to/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d gene_trees_$i.nwk -m tree_mod_$n.tre -o phyparts/scenario_output_$i.$n
	done
	n=$((n+1))
done
