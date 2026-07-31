#!/bin/bash -f

clust=$1;

for ((sim=1; sim<=10; sim++)); #loop through sims in cluster
do
	label=$((clust*10 + sim))
  sed "s/XXXXX/$label/g" ion_atmosphere_cluster${clust}.ptraj > tempion_atmosphere_cluster${clust}_sim${sim}.ptraj;
  sed "s/YYYYY/$sim/g" tempion_atmosphere_cluster${clust}_sim${sim}.ptraj > ion_atmosphere_cluster${clust}_sim${sim}.ptraj;
  cpptraj -i ion_atmosphere_cluster${clust}_sim${sim}.ptraj;
  rm tempion_atmosphere_cluster${clust}_sim${sim}.ptraj;
  rm ion_atmosphere_cluster${clust}_sim${sim}.ptraj;
done
