#!/bin/bash  
#SBATCH --partition=partMath
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH --array=0-0
#SBATCH --exclusive
#SBATCH --output=./job/%j.out
#j1 is the eta,j2 is the fB, j3=beta(v),j4=fR;
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

         a1=`echo "scale=3;${fA}+0.005*${SLURM_ARRAY_TASK_ID}" | bc`;
        mpirun -n 1 ./ab ${phase} ${a1} ${chiAB} ${dN} ${ItMax} ${Ndeg}>phase-${phase}_f[${a1}_[${chiAB}]-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt &
wait                
