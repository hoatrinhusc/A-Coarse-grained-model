#! /bin/csh -f

if ($#argv != 2) then
        echo "Usage:   		./step1.csh  pdb       name"
        echo "Example: 		./step1.csh  1CFD.pdb  cfd"
	goto done
    endif
set pdb = $argv[1] 
set name = $argv[2]
###################   (1) DEFINITIONS  ###########################

set PDB = $name.pdb
set NATIVE = $name.ab.crd
set ATOM = INDEX.$name.ATOM
set NUM = $name.num
###################   (2) CLEANING UP PDB FILE  ###########################

./removeH_OT.pl $pdb > $PDB	#removing H atoms and OT's
./numbers.pl $PDB > $NUM      #finding natom ntype
sed -i "s/HSE/HIS/g" $PDB
set natom = `head -1 $NUM | awk '{print $1}'`	
set ntype = `tail -1 $NUM | awk '{print $1}'`	

###################   (3) COMPILING beta_coor.c  ##########################
#updating the header file IOab.h
sed "s/define MAXLENGTH  NTYPE/define MAXLENGTH  $ntype/g;s/define NATOM  NATOM/define NATOM  $natom/g;s/define NTYPE  NTYPE/define NTYPE  $ntype/g"  IOab.header > IOab.h

cc beta_coor.c -lm -w -o beta

###################   (4) CREATING COARSE GRAINED CRD & INDEX files  ###########################

echo $PDB > input.beta
./beta < input.beta > output.beta

echo $name > HEAD
echo $ntype >> HEAD
awk '/./ {print $1/3.8, $2/3.8, $3/3.8}'  < $PDB.beta > b.crd
cat HEAD b.crd > a.crd

#formating the crd to AMBER crd
cc format.c -lm -w -o format
./format a.crd > $NATIVE

awk '{print $1,$2, $3/3.8}' $PDB.atom > $ATOM

###################   (5) CLEANING  ###########################

rm -f a.crd b.crd HEAD input.beta output.beta $PDB.beta $PDB.atom

###################   EXIT  ###########################
done:
     exit 0
