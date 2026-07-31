#!/bin/bash -f

for ((i=257; i<=449; i++)); #Na+ ions
do
  sed "s/XXXXX/$i/g" Na_contact.ptraj > Na_${i}_distance.ptraj;
 cpptraj -i Na_${i}_distance.ptraj;
 rm Na_${i}_distance.ptraj;
done

