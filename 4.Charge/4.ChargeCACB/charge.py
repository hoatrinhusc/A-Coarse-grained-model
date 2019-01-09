# Written by Hoa Trinh on 2/2/2018
# Might not work with python version > 3
# Check the name of atoms in N and C terminal: usually, Oxygen in C termial will be OXT; N terminal will have 3 hydrogen atoms

import math
fin = open('MopacCharge', 'r')
fout = open('CACB.charge', 'w')

BB = ['N', '1H', '2H', '3H', 'H', 'CA', 'HA', '1HA', '2HA', 'C', 'O', 'OXT']

first = fin.readline()
field = first.split()
resnum = field[6]
chargeCA = chargeCB = 0
scale = math.sqrt(38*0.008314*300)

fin.close()
fin = open('MopacCharge', 'r')
 
for line in fin:
	words = line.split()
	if words[6] == resnum and words[3] in BB:
		chargeCA += float(words[7])
		
 	if words[6] == resnum and words[3] not in BB: 
                chargeCB += float(words[7])
		
	if words[6] != resnum:
		a = list(words[6])
		#fout.write('{0} {1} {2}\n'.format(resnum, 'CA', chargeCA))
		fout.write('{0}\n'.format(chargeCA/scale))
		if chargeCB != 0:
			fout.write('{0}\n'.format(chargeCB/scale))
		#fout.write('{0} {1} {2}\n'.format(resnum, 'CB', chargeCB))
		resnum = words[6]
		chargeCA = float(words[7])
		chargeCB = 0
fin.close()

#Write the charge of CA and CB of the last residue
fout.write('{0}\n'.format(chargeCA/scale))
if chargeCB != 0:
	fout.write('{0}\n'.format(chargeCB/scale))

#fout.write('{0} {1} {2}\n'.format(resnum, 'CA', chargeCA))
#fout.write('{0} {1} {2}\n'.format(resnum, 'CB', chargeCB))
fout.close()



  
		
