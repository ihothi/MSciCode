#!/bin/tcsh -f
#PBS -N RebinRun
#PBS -l nodes=1:ppn=1,mem=40GB
#PBS -q compute
#PBS -o /home/ihothi/
#PBS -j oe
#PBS -S /bin/tcsh

cd $PBS_O_WORKDIR
source /etc/profile.d/modules.csh
module load dev_tools/oct2017/python-Anaconda-3-5.0.0.1

echo Working directory is $PBS_O_WORKDIR
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo This jobs runs on the following processors:
echo `cat $PBS_NODEFILE`

python3.6 TestSplinterFull.py

echo Time is `date`
