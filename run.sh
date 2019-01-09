#!/bin/sh
#This script will compile all the codes in 0.Cheungcode.Hoa, 1.Preparation, 2.Topology
#Usage
echo "Usage:            sh run.sh     pdb                name   chain  model      map"
echo "Example:          sh run.sh     noLoopsSOD1.pdb    SOD     A      GO        BTmap.dat"
#Excute 0.Cheungcode.Hoa/step1.csh and cp files
cd 0.Cheungcode.Hoa/Multi_Steps/1.CACB
./step1.csh  $1  $2
cp $2.num ../2.HB
cp $2.pdb ../2.HB
cp $2.num ../3.SC
cp $2.pdb ../3.SC
cp $2.ab.crd ../4.PRMTOP
cp INDEX.$2.ATOM ../4.PRMTOP
cp $2.pdb ../../../1.Preparation
cp $2.pdb.CBcrd ../../../1.Preparation
cp $2.num ../../../1.Preparation
cp $2.num ../../../2.Topology
cp INDEX.$2.ATOM ../../../2.Topology
#Excute 0.Cheungcode.Hoa/step2.csh and cp files
cd ../2.HB
./step2.csh $2
cp INDEX.$2.HB ../3.SC
#Excute 0.Cheungcode.Hoa/step3.csh and cp files
cd ../3.SC
./step3.csh $2 $3
cp INDEX.$2.SC ../4.PRMTOP
cp INDEX.$2.HB ../4.PRMTOP
cp INDEX.$2.PAIR ../4.PRMTOP
cp $2.pdb ../4.PRMTOP
cp $2.num ../4.PRMTOP
cp IOab.h ../4.PRMTOP
cp IOab.1.h ../4.PRMTOP
#Excute 0.Cheungcode.Hoa/step4.csh and cp files
cd ../4.PRMTOP
./step4.csh $2 $4 $5
cp BondEnergy.inp ../../../2.Topology
cp pseudo.inp ../../../2.Topology
#Excute 1.Preparation
cd ../../../1.Preparation
./step1.csh $2
cp $2.ab.pdb ../2.Topology
#Excute 2.Topology
cd ../2.Topology
#Remove the first line of "pseudo.inp"
sed -i '1d' pseudo.inp
./step2.csh $2
echo "Finish compiling 0.Cheungcode.Hoa, 1.Preparation, 2.Topology"
echo "In 2.Topology, $2.ab.pdb and $2.top will be used for Gromacs simulations"
