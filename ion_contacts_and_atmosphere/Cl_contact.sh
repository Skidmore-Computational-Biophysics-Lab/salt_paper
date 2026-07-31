#!/bin/bash -f

for ((i=73; i<=256; i++)); #Cl- ions
do
#  cp contact_distance.ptraj ${i}_distance.ptraj;
  sed "s/XXXXX/$i/g" Cl_contact.ptraj > Cl_${i}_distance.ptraj;
 cpptraj -i Cl_${i}_distance.ptraj;
 rm Cl_${i}_distance.ptraj;
done
