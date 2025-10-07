#!/bin/bash
clear
phase=`echo "scale=2;1"| bc`;
fA=`echo "scale=2;0.4" | bc`;

chiAB=`echo "scale=2;0.14" | bc`; 
dN=`echo "scale=2; 200" | bc`;
ItMax=`echo "scale=2; 10000" | bc`;
Ndeg=`echo "scale=2; 100" | bc`;

if [ ! -d "results-${phase}/" ];
then
        mkdir -p results-${phase}/{M,rho,qu,Data,energy,Diff}
fi

        mpirun -n 4 ./ab ${phase} ${fA} ${chiAB} ${dN} ${ItMax} ${Ndeg}>phase-${phase}_f[${fA}_[${chiAB}].txt &
wait

