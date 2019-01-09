#! /bin/csh -f

if ($#argv != 3) then
        echo "Usage:   		./step4.csh  name   model   map"
        echo "Example: 		./step4.csh  cfd    GO      BTmap.dat"
	goto done
    endif

set name = $argv[1]
set model = $argv[2]
set map = $argv[3]

###################   (1) DEFINITIONS  ###########################

set PDB = $name.pdb
set NATIVE = $name.ab.crd
set TOP = $name.prmtop
set ATOM = INDEX.$name.ATOM
set SC = INDEX.$name.SC
set HB = INDEX.$name.HB
set PAIR = INDEX.$name.PAIR
set NUM = $name.num
set leapscript = leap_script

set natom = `head -1 $NUM | awk '{print $1}'`	
set ntype = `tail -1 $NUM | awk '{print $1}'`	
set pair = `wc -l $PAIR | awk '{print $1}'`	
set hb = `wc -l $HB | awk '{print $1}'`	
###################   (2) FORCE FIELD  ###########################

#updating the header file IO.h

sed "s/define MAXLENGTH NTYPE/define MAXLENGTH $ntype/g;s/define NATOM NATOM/define NATOM $natom/g;s/define NTYPE NTYPE/define NTYPE $ntype/g;s/define MAXPAIR NPAIR/define MAXPAIR $pair/g;s/define MAXHBPAIR HBPAIR/define MAXHBPAIR $hb/g"  IOablinux.header > IOablinux.h
cc frcfield_ab.xy.amber10.c -lm -w -o frc	#compiling frcfield
./frc $NATIVE $ATOM $PAIR		#creating the force field file
sed -i "s/EP/ZP/g" frc.go               #removing EP because it means "Extra Points" in Amber10


###################   (3) PRMTOP  ###########################

setenv AMBERHOME ~/amber10-PZ
./tlp.cacb.pl $PDB > $leapscript 	#creating tleap script
sed -i "s/EP/ZP/g" $leapscript          #removing EP because it means "Extra Points" in Amber10
$AMBERHOME/exe/tleap -f $leapscript


###################   (4) NONBONDED INTERACTION  ###########################

sed "s/ATOM = TYPE/ATOM = $ntype/g;s/RESNUM = ATOM/RESNUM = $natom/g;s/MAXCONTACT = PAIR/MAXCONTACT = $pair/g" calc6-12out.$model.parm7.tmp > calc6-12out.$model.parm7.pl

#Hoa May 22: to make sure calc6-12out.$model.parm7.pl is excutable
chmod +x calc6-12out.$model.parm7.pl
#Hoa end

./calc6-12out.$model.parm7.pl $ATOM $map $PAIR $NATIVE 
sed -i "s/FORMAT/%FORMAT(5E16.8)/g" 6-12out.dat.att
sed -i "s/FORMAT/%FORMAT(5E16.8)/g" 6-12out.dat.rep

sed '/%FLAG LENNARD_JONES_ACOEF/,/%FLAG LENNARD_JONES_BCOEF/c\%FLAG LENNARD_JONES_ACOEF\n%FLAG LENNARD_JONES_BCOEF' cacb.prmtop > 0.top
sed '/%FLAG LENNARD_JONES_BCOEF/,/%FLAG BONDS_INC_HYDROGEN/c\%FLAG LENNARD_JONES_BCOEF\n%FLAG BONDS_INC_HYDROGEN' 0.top > 1.top
sed '/%FLAG LENNARD_JONES_ACOEF/ r  6-12out.dat.rep'  1.top > 0.top
sed '/%FLAG LENNARD_JONES_BCOEF/ r  6-12out.dat.att'  0.top > $TOP

###################   (5) PRODUCING psuedo.inp trip.inp  ###########################
if($model == "GO")  then
   sed "s/NONNATIVE XXX/NONNATIVE NO/g" countHB.tmp > countHB.c
endif
if ($model == "BT") then 
   sed "s/NONNATIVE XXX/NONNATIVE YES/g" countHB.tmp > countHB.c 
endif
cc countHB.c -lm -w -o countHB.x	#compile countHB.c
cc triple.c -lm -w -o triple	#compile triple.c
sed 2d $NATIVE > vmd.crd
echo $NATIVE $ATOM $PAIR 1 vmd.crd > countHB.in
echo $NATIVE $ATOM $PAIR 20 > triple.in
./countHB.x < countHB.in
./triple < triple.in

###################   (6) CLEANING UP  ###########################

rm -f cacb.prmtop cacb.crd 0.top 1.top
rm -f IO.h $leapscript frc frc.go

###################   EXIT  ###########################
done:
     exit 0
