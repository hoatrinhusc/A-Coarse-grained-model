#! /bin/csh -f

if ($#argv != 1) then
        echo "Usage:   		./step2.csh  pdb  "
        echo "Example: 		./step2.csh  ci2  "
	goto done 
    endif
set pdb = $argv[1]

################################# Count number of native pairs ###########################

echo $pdb > input
g++ countline.cpp -o count
./count < input
 
###################   (1) DEFINITIONS  ###########################
set NUM = $pdb.num

set numCA = `head -1 $NUM `	
set numAB = `sed -n '2p' $NUM `
set numCB = `expr $numAB - $numCA`
set npair = `head -1 $pdb.npair`
set npseudo = `head -1 $pdb.npseudo`
	
###################   (3) COMPILING beta_coor.c  ##########################
#updating the header file common.h
sed "s/define NUMCA NUMCA/define NUMCA $numCA/g;s/define NUMAB NUMAB/define NUMAB $numAB/g;s/define NUMCB NUMCB/define NUMCB $numCB/g;s/define NPAIR NPAIR/define NPAIR $npair/g;s/define NPSEUDO NPSEUDO/define NPSEUDO $npseudo/g" common.header > common.h

g++ main.cpp -o main

###################   (4) CREATING COARSE GRAINED CRD & INDEX files  ###########################

echo $pdb > input.beta
./main < input.beta

rm input.beta
###################   EXIT  ###########################
done:
     exit 0
