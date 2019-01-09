#! /bin/csh -f

if ($#argv != 2) then
        echo "Usage:   		./step3.csh  name   chain"
        echo "Example: 		./step3.csh  cfd    A"
	goto done
    endif

set name  = $argv[1]
set chain = $argv[2]

###################   (1) DEFINITIONS  ###########################

set PDB = $name.pdb
set SC = INDEX.$name.SC
set HB = INDEX.$name.HB
set PAIR = INDEX.$name.PAIR
set NUM = $name.num
set natom = `head -1 $NUM | awk '{print $1}'`	
set ntype = `tail -1 $NUM | awk '{print $1}'`	

###################   (2) SC NATIVE CONTACTS  ###########################

./makeSUN 	#compile resc.c

set Max = $natom
@ Max = $Max + 1

@ count = 1 

echo END > tmp
cat $PDB tmp > tmpPDB
while ($count < $Max)
./resc tmpPDB $count $chain  >> a
@ count = $count + 1
end
perl -wnle 'BEGIN {$p=0}; $p=1 if /Table II\b/; $p=0 if /Table III/; print $_ if $p' a >> b
perl -wnle 'BEGIN {$p=0}; $p=1 if /Table III\b/; $p=0 if /Table IV/; print $_ if $p' a >> c
perl -wnle 'BEGIN {$p=0}; $p=1 if /Table IV\b/; $p=0 if /Table I/; print $_ if $p' a >> d

./makemap.ab.pl a

awk '{print NR, "b", $3-1,$4-1}'  MapSC > $SC

###################   (3) JOINING ALL CONTACT PAIRS   ###########################


cat $SC $HB | awk '{print NR,$2,$3,$4}' > $PAIR


###################   (4) UPDATING HEADER FILE   ###########################

set pair = `wc -l $PAIR | awk '{print $1}'`	
set hb = `wc -l $HB | awk '{print $1}'`	

sed "s/define MAXLENGTH  NTYPE/define MAXLENGTH  $ntype/g;s/define NATOM  NATOM/define NATOM  $natom/g;s/define NTYPE  NTYPE/define NTYPE  $ntype/g;s/define MAXPAIR 500/define MAXPAIR $pair/g;s/define MAXHBPAIR 200/define MAXHBPAIR $hb/g"  IOab.header > IOab.h
sed "s/define MAXLENGTH NTYPE/define MAXLENGTH $ntype/g;s/define NATOM NATOM/define NATOM $natom/g;s/define NTYPE NTYPE/define NTYPE $ntype/g;s/define MAXPAIR 500/define MAXPAIR $pair/g;s/define MAXHBPAIR 200/define MAXHBPAIR $hb/g"  IOab.1.header > IOab.1.h

###################   (5) CLEANING   ###########################

/bin/rm a b c d MapSC tmp tmpPDB


###################   EXIT  ###########################
done:
     exit 0
