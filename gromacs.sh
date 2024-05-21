#!/bin/bash
pattern="*.gro"
files=( $pattern )
echo "${files[0]}"  # printf is safer!

gmx grompp -f ../zero.mdp -p ../topol.top -c "${files[0]}" -maxwarn 10 -o zero.tpr &> /dev/null
rm mdout.mdp

for i in *gro
do
echo $i
gmx mdrun -s zero.tpr -rerun $i &> /dev/null
gmx energy -f ener.edr -o $i.xvg <<< "5 7 8 6 9 10 11" &>/dev/null
rm md.log ener.edr
done

rm zero.tpr
