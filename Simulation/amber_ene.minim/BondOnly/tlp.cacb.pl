#!/usr/bin/perl -w
################################################################################
#		A script to build a Calpha-Cbeta model in amber10              #
#		#######   Antonios Samiotakis 4.13.2010   #######              #
################################################################################

if(@ARGV < 1) {
        print "Usage: tlp.pl file.pdb\n";
        exit;
        }

$pdb = shift @ARGV;

open(PDB,"$pdb") or die "failed to open $pdb, ";


@array = ();
@ltrca = ();
@ltrcb = ();

for(my $j=0;$j<500;$j++){
  my $i= int $j/36;
  my $ja = ($j%36);
  my $i1= 65 + $i;
  if ($ja<10){$i2 = 48 + $ja;}
  else {$i2 = 65 +$ja -10;}
	
  my $c= sprintf("%c%c",$i1,$i2);
#  printf("%c%c\n",$i1,$i2);
  if($c eq 'DU'){$c = 'Z9';}
  push (@ltrca, $c);
}

for(my $j=0;$j<500;$j++){
  my $i= int $j/36;
  my $ja = ($j%36);
  my $i1= 78 + $i;
  if ($ja<10){$i2 = 48 + $ja;}
  else {$i2 = 65 +$ja -10;}

  my $c= sprintf("%c%c",$i1,$i2);
#  printf("%c%c\n",$i1,$i2);
  if($c eq 'DU'){$c = 'Z9';}
  push (@ltrcb, $c);
}



printf "cacb_unit = createUnit CACB\n";




my $count = 0;
my $res = 0;
my $atom = 0;

while (<PDB>) {


           $array[$count][0] = substr($_,0,4);
           $array[$count][1] = substr($_,4,7);
           $array[$count][2] = substr($_,11,5);
           $array[$count][3] = substr($_,16,1);
           $array[$count][4] = substr($_,17,3);
           $array[$count][5] = substr($_,20,3);
           $array[$count][6] = substr($_,23,3);
           $array[$count][7] = substr($_,26,12);
           $array[$count][8] = substr($_,38,8);
           $array[$count][9] = substr($_,46,8);
           $array[$count][10] = substr($_,54,6);
           $array[$count][11] = substr($_,60,6);
           $array[$count][12] = substr($_,66,12);
                
#          printf "$array[$count][4]\n";              

        $count++;
        }



#exit;

my $gly_flag = 0;
for ($l = 0; $l < $count; $l++) {
	

	
	   if(($array[$l][2] eq '  CA ') && ($array[$l][4] ne 'GLY')){
			printf "mon$atom = createAtom $ltrca[$res] $ltrca[$res] 0.0\n";
			printf "set mon$atom position { $array[$l][7]\t $array[$l][8]\t$array[$l][9]}\n";
			printf "cacb_res$res = createResidue R$res\n";
			printf "add cacb_res$res mon$atom\n";
			if ($atom >0){
					if ($gly_flag eq 1){
							$k = $atom-1;
							printf "bond mon$atom mon$k\n";
							}
					else {
						$k = $atom - 2;
						printf "bond mon$atom mon$k\n";
			                     }
					} 
			$atom = $atom +1;
			$gly_flag = 0;
				}


	    if(($array[$l][2] eq '  CA ') && ($array[$l][4] eq 'GLY')){
                        printf "mon$atom = createAtom $ltrca[$res] $ltrca[$res] 0.0\n";
                        printf "set mon$atom position { $array[$l][7]\t $array[$l][8]\t$array[$l][9]}\n";
                        printf "cacb_res$res = createResidue R$res\n";
                        printf "add cacb_res$res mon$atom\n";
			printf "add cacb_unit cacb_res$res\n";
			if ($atom >0){
                                        if ($gly_flag eq 1){
                                                        $k = $atom-1;
                                                        printf "bond mon$atom mon$k\n";
                                                        }
                                        else {
                                                $k = $atom - 2;
                                                printf "bond mon$atom mon$k\n";
                                             }
                                        } 
                        $atom = $atom +1;
			$res = $res +1;
			$gly_flag = 1;
                                }


	   if($array[$l][2] eq '  CB '){
                        printf "mon$atom = createAtom $ltrcb[$res] $ltrcb[$res] 0.0\n";
                        printf "set mon$atom position { $array[$l][7]\t $array[$l][8]\t$array[$l][9]}\n";
                        printf "add cacb_res$res mon$atom\n";
                        printf "add cacb_unit cacb_res$res\n";
			$k = $atom - 1;
                        printf "bond mon$atom mon$k\n";
			$res = $res + 1;
			$atom = $atom + 1;
                                }

	   

}

printf "myforce = loadamberparams frc.go\n";
printf "saveAmberParm cacb_unit cacb.prmtop cacb.crd\n";
printf "quit\n";
