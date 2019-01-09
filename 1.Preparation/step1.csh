#! /bin/csh -f

if ($#argv != 1) then
        echo "Usage:   		./step1.csh  pdb  "
        echo "Example: 		./step1.csh  ci2  "
	goto done 
    endif
set pdb = $argv[1]

###################   (0) COUNT NUMBER OF NATIVE PAIRS #####################
echo $pdb > input
g++ countline.cpp -o count
./count < input
 
###################   (1) DEFINITIONS  #####################################
set NUM = $pdb.num
set NUMAT = $pdb.numat

set numCA = `head -1 $NUM `	
set numAB = `sed -n '2p' $NUM `
set natom = `head -1 $NUMAT `
set numCB = `expr $numAB - $numCA`
	
###################   (2) UPDATE HEADER FILE  ##############################
sed "s/define NUMCA NUMCA/define NUMCA $numCA/g;s/define NATOM NATOM/define NATOM $natom/g;s/define NUMAB NUMAB/define NUMAB $numAB/g;s/define NUMCB NUMCB/define NUMCB $numCB/g" common.header > common.h

###################   (3) CREATE PDB FILE OF CA & CB AND NATIVE MAP  ###########################
g++ pdb.cpp -o process
./process < input

################### CLEANING###########################
rm input
###################   EXIT  ###########################
done:
     exit 0
