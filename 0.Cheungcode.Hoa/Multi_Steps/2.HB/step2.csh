#! /bin/csh -f

if ($#argv != 1) then
        echo "Usage:   		./step2.csh  name"
        echo "Example: 		./step2.csh  cfd"
	goto done
    endif

set name  = $argv[1]

###################   (1) DEFINITIONS  ###########################

set PDB = $name.pdb
set HB = INDEX.$name.HB

###################   (2) HB NATIVE CONTACTS  ###########################

#./dsspcmbi $PDB dssptmp				#using DSSP to generate HB list
./dssp $PDB dssptmp
./hbond-dssp.pl dssptmp > tmplist		#extracting HB pairs
./clean_HBlist.pl tmplist > HBondlist		#Cleaning the list with two rules
						# 1) remove duplicated entry. 2) ensure only a maximum of *2* HB pairs per residue. 
						# If there are more than two, delete the pair that has least bonding energy.

awk '{print NR,"h",$1-1,$2-1}' HBondlist > $HB

###################   (3) CLEANING   ###########################

/bin/rm dssptmp tmplist HBondlist

###################   EXIT  ###########################
done:
     exit 0
