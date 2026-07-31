#!/bin/bash -f

declare -a arr=("60" "61" "63" "69" "74" "76")
for i in "${arr[@]}"
do
  echo "NUMBER ${i}"
  sed "s/XXXXX/$i/g" ../ptraj/charged_residue_contact.ptraj > ../ptraj/charged_residue_${i}_distance.ptraj;
 cpptraj -i ../ptraj/charged_residue_${i}_distance.ptraj;
 rm ../ptraj/charged_residue_${i}_distance.ptraj;
done

