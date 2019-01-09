#include "common.h"

int main()
{  
  ATOMIC CA,CB;
  int atnum, resnum, reIndex, count, index, i, j, countCB, atom1[NPAIR], atom2[NPAIR];
  string atom, dummy, residue;
  double dih, radian, bondAB, bondAA, rCB, x[NUMAB], y[NUMAB], z[NUMAB], epsilon[NPAIR];
  double cosA, angA, cosB, angB, dist, sigma, xp, yp, sigij;
  vector<double> ai, aj, ak, al, ij, jk, kl, ik, il, n1, n2, m1;
  int xi0, xi1, xi2, xi3;
  char pdb[20], input1[30], input2[30], input3[30], input4[30], outfile[30]; 

  FILE * TopFile;
        scanf ("%s", pdb);
	sprintf (input1, "%s.ab.pdb", pdb);
	sprintf (input2, "BondEnergy.inp");
 	sprintf (input3, "INDEX.%s.ATOM", pdb);
	sprintf (input4, "pseudo.inp");
	ifstream fp1(input1), fp2(input2), fp3(input3), fp4(input4);
        sprintf(outfile,"%s.top",pdb);
        TopFile = fopen(outfile,"wt");

  index = 1;
  count = 0;
  
  if (fp1.fail())
  {
      cout << "Cannot open PDB file" << endl;
      return 1;
  }
  
  if (fp2.fail())
  {
      cout << "Cannot open native map file" << endl;
      return 1;
  }
  
  if (fp3.fail())
  {
      cout << "Cannot open pdb.atom file" << endl;
      return 1;
  }

  if (fp4.fail())
  {
      cout << "Cannot open pseudo.inp file" << endl;
      return 1;
  }


/* Read coordinates from pdb file input    */
/* Convert unit from angstrom to nanometer */
  for(i = 0; i < NUMAB; i++)
  { 
      fp1 >> dummy >> atnum >> atom >> residue >> dummy >> resnum >> x[i] >> y[i] >> z[i];
      reIndex = resnum - 1;

      if(residue.compare("GLY")!=0 && atom.compare("CA")==0)
      {
	  CA.x[reIndex] = x[i];
	  CA.y[reIndex] = y[i];
	  CA.z[reIndex] = z[i];
	  CA.num[reIndex] = atnum;
      }
      else if(residue.compare("GLY")!=0 && atom.compare("CB")==0)
      {
	  CB.x[reIndex] = x[i];
	  CB.y[reIndex] = y[i];
	  CB.z[reIndex] = z[i];
	  CB.num[reIndex] = atnum;
      }
      else
      {
	  CA.x[reIndex] = x[i];
	  CA.y[reIndex] = y[i];
	  CA.z[reIndex] = z[i];
	  CA.num[reIndex] = atnum;
	  CB.x[reIndex] = CB.y[reIndex] = CB.z[reIndex] = CB.num[reIndex] = 0;
      }        
  }

  fprintf (TopFile, "; Topology file for the coarse-grained model described in the following paper: \n");
  fprintf (TopFile, "; Exploring the Interplay between Topology and Secondary Structural Formation in the Protein Folding Problem \n");
  fprintf (TopFile, "; Cheung et al, J. Phys. Chem. B 2003, 107, 11193-11200 \n \n");
  
/* Defaults Directive */
  fprintf (TopFile, "[ defaults ] \n");
  fprintf (TopFile, "; nbfunc comb-rule gen-pairs \n");
  fprintf (TopFile, "  1      2         no \n \n");
  
/* Atomtype Directive */
  fprintf (TopFile, "[ atomtypes ] \n");
  fprintf (TopFile, "; name    mass   charge    ptype    sigma    epsilon \n ");
  fprintf (TopFile, "%3s CA00   50.000  0.000  A   %11.9e   %11.9e \n ", " ",-2*0.9*0.1*0.5, 0.15);
/* f=0.9: scaling factor; 0.1: convert from angstrom to nm; 0.5: radius of Calpha in reduced unit; m_CA = m_CB = 50 */
/* epsilon = 0.6/4 = 0.15 */

  for (i = 0; i < NUMCA; i++)
  {
      fp3 >> resnum >> residue >> rCB;
  
      if (residue.compare("GLY")!=0)
      {
	  fprintf (TopFile, "%3s CB%02d   50.000  0.000  A   %11.9e   %11.9e \n ", " ", count, rCB*(-2*0.9*0.1), 0.15);
	  count++;
      }
  }
  fp3.close();
  fprintf (TopFile, " \n");
   
/* Molecule Type Directive */
  fprintf (TopFile, "[ moleculetype ] \n");
  fprintf (TopFile, "; name            nrexcl \n");
  fprintf (TopFile, "  Macromolecule   3 \n \n");

/* Atoms Directive */
  fprintf (TopFile, "[ atoms ] \n");
  fprintf (TopFile, "; nr    type   resnr    residue    atom    cgnr  charge \n ");
  fp3.open(input3);

  count = 0;
  for (i = 0; i < NUMCA; i++)
  {
      fp3 >> resnum >> residue >> rCB;
      if (residue.compare("GLY")!=0)
      {
          fprintf (TopFile, "%6d %7s CA00 %7d %6s %6s %6d \n ", index," ", resnum+1, residue.c_str(), "CA", index);
          index ++;
          fprintf (TopFile, "%6d %7s CB%02d %7d %6s %6s %6d \n ", index, " ", count, resnum+1, residue.c_str(), "CB", index);
          count ++;
          index ++;
      }
      else
      {
          fprintf (TopFile, "%6d %7s CA00 %7d %6s %6s %6d \n ", index, " ", resnum+1, residue.c_str(), "CA", index);
          index ++;
      }
  }
  fprintf (TopFile, " \n");

  
/* Bond Directive */
  fprintf (TopFile, "[ bonds ] \n");
  fprintf (TopFile, "; ai     aj      func    r0(nm)  Kb \n");

  for(i = 0; i < NUMCA - 1; i++)
  {
      ai = COOR(&CA.x[i], &CA.y[i], &CA.z[i]);
      aj = COOR(&CB.x[i], &CB.y[i], &CB.z[i]);
      ak = COOR(&CA.x[i+1], &CA.y[i+1], &CA.z[i+1]);
      ik = VECTOR(&ai, &ak); 
    
      bondAA = magnitude(&ik)/10; //Convert Angstrom (PDB) to nm (Gro)
      fprintf (TopFile, "%9d %9d %4d  %11.9e  %11.9e \n", CA.num[i], CA.num[i+1], 1, bondAA, 24000.);
/* Kb = 120*2*100 = 24000 */ 
    
      if (CB.x[i] == 0 && CB.y[i] == 0 && CB.z[i] == 0)
      {
	  continue;
      }
      ij = VECTOR(&ai, &aj);
      bondAB = magnitude(&ij)/10; //Convert Angstrom (PDB) to nm (Gro)
      fprintf (TopFile, "%9d %9d %4d  %11.9e  %11.9e \n", CA.num[i], CB.num[i], 1, bondAB, 24000.);
  }
  fprintf (TopFile, " \n");
  
/* Angle Directive */
  fprintf (TopFile, "[ angles ] \n");
  fprintf (TopFile, "; ai  aj   ak  func  th0(deg)   Ka \n");

/* First atom */
  if (CB.x[0] != 0 || CB.y[0] != 0 || CB.z[0] != 0)
  {
      ai = COOR(&CA.x[0], &CA.y[0], &CA.z[0]);
      aj = COOR(&CB.x[0], &CB.y[0], &CB.z[0]);
      ak = COOR(&CA.x[1], &CA.y[1], &CA.z[1]);
      ik = VECTOR(&ai, &ak);
      ij = VECTOR(&ai, &aj);
      angA = angle(&ij, &ik) ;
      fprintf (TopFile, "%9d %9d %9d %7d  %11.9e  %11.9e \n", CB.num[0], CA.num[0], CA.num[1], 1, angA, 48.);
/* Ka = 24*2 = 48 */
  }

  for (i = 1; i < NUMCA - 1; i++)
  {
      ai = COOR(&CA.x[i], &CA.y[i], &CA.z[i]);
      aj = COOR(&CB.x[i], &CB.y[i], &CB.z[i]);
      ak = COOR(&CA.x[i+1], &CA.y[i+1], &CA.z[i+1]);
      al = COOR(&CA.x[i-1], &CA.y[i-1], &CA.z[i-1]);
      ik = VECTOR(&ai, &ak); 
      il = VECTOR(&ai, &al);
      angA = angle(&ik, &il) ;
      fprintf (TopFile, "%9d %9d %9d %7d  %11.9e  %11.9e \n", CA.num[i-1], CA.num[i], CA.num[i+1], 1, angA, 48.);
   
      if (CB.x[i] == 0 && CB.y[i] == 0 && CB.z[i] == 0)
      {
	  continue;
      }
    
     ij = VECTOR(&ai, &aj);
     angB = angle(&ij, &il);
     fprintf (TopFile, "%9d %9d %9d %7d  %11.9e  %11.9e \n", CA.num[i-1], CA.num[i], CB.num[i], 1, angB, 48.);
    
     angB = angle(&ij, &ik);
     fprintf (TopFile, "%9d %9d %9d %7d  %11.9e  %11.9e \n", CB.num[i], CA.num[i], CA.num[i+1], 1, angB, 48.);
  }

/* Last atom */
  if (CB.x[NUMCA-1] != 0 || CB.y[NUMCA-1] != 0 || CB.z[NUMCA-1] != 0)
  {
      ai = COOR(&CA.x[NUMCA-1], &CA.y[NUMCA-1], &CA.z[NUMCA-1]);
      aj = COOR(&CB.x[NUMCA-1], &CB.y[NUMCA-1], &CB.z[NUMCA-1]);
      ak = COOR(&CA.x[NUMCA-2], &CA.y[NUMCA-2], &CA.z[NUMCA-2]);
      ik = VECTOR(&ai, &ak);
      ij = VECTOR(&ai, &aj);
      angA = angle(&ij, &ik) ;
      fprintf (TopFile, "%9d %9d %9d %7d  %11.9e  %11.9e \n", CA.num[NUMCA-2], CA.num[NUMCA-1], CB.num[NUMCA-1], 1, angA, 48.);
  }
  fprintf (TopFile, " \n");

 
/* Dihedral Directive */
  fprintf (TopFile, "[ dihedrals ] \n");
  fprintf (TopFile, "; ai  aj  ak  al  func  phi0(deg) kd mult \n");

  for(i = 0; i < NUMCA-3; i++)
  {
      ai = COOR(&CA.x[i], &CA.y[i], &CA.z[i]);
      aj = COOR(&CA.x[i+1], &CA.y[i+1], &CA.z[i+1]);
      ak = COOR(&CA.x[i+2], &CA.y[i+2], &CA.z[i+2]);
      al = COOR(&CA.x[i+3], &CA.y[i+3], &CA.z[i+3]);
      
      dih = dih_angle(&ai, &aj, &ak, &al);
      
      fprintf (TopFile, "%9d %9d %9d %9d %4d   %11.9e  %11.9e %4d \n",
			    CA.num[i], CA.num[i+1], CA.num[i+2], CA.num[i+3], 1, dih, 0.6, 1);
      fprintf (TopFile, "%9d %9d %9d %9d %4d   %11.9e  %11.9e %4d \n",
			    CA.num[i], CA.num[i+1], CA.num[i+2], CA.num[i+3], 1, 3*dih, 0.3, 3);
 }


/* Hydrogen bond angular dependence 03/29/2018 */
  for (i=0; i<NPSEUDO; i++)
  {
	fp4 >> xi0 >> xi1 >> xi2 >> xi3;
	ai = COOR(&x[xi1-1], &y[xi1-1], &z[xi1-1]);
        aj = COOR(&x[xi2-1], &y[xi2-1], &z[xi2-1]);
	ij = VECTOR(&ai, &aj);
        sigij = magnitude(&ij)/10;
  
        fprintf (TopFile, "%9d %9d %9d %9d %4d   %11.9e %4d \n", xi0, xi1, xi2, xi3, 10, sigij, 1);
  }

  fp4.close();


/* Improper dihedral angles */
  for(i = 1; i < NUMCA-1; i++)
  { 
      if (CB.x[i] != 0 || CB.y[i] != 0 || CB.z[i] != 0)
      {
	  ai = COOR(&CA.x[i-1], &CA.y[i-1], &CA.z[i-1]);
	  aj = COOR(&CA.x[i+1], &CA.y[i+1], &CA.z[i+1]);
	  ak = COOR(&CA.x[i], &CA.y[i], &CA.z[i]);
	  al = COOR(&CB.x[i], &CB.y[i], &CB.z[i]);

	  dih = idih_angle(&ai, &aj, &ak, &al);

	  fprintf (TopFile, "%9d %9d %9d %9d %4d   %11.9e  %11.9e \n",
				CA.num[i-1], CA.num[i+1], CA.num[i], CB.num[i], 2, dih, fc);
      }
 }
 fprintf (TopFile, " \n");


/* Pairs Directive */
  fprintf (TopFile, "[ pairs ] \n");
  fprintf (TopFile, "; ai aj type, A, B \n");
  
  for(i = 0; i < NPAIR; i++)
  {
      fp2 >> atom1[i] >> atom2[i] >> epsilon[i];
  }

  for(i = 0; i < NPAIR; i++)
  {
      dist = sqrt(pow(x[atom1[i]-1]-x[atom2[i]-1],2.0) + pow(y[atom1[i]-1]-y[atom2[i]-1],2.0) + pow(z[atom1[i]-1]-z[atom2[i]-1],2.0));
      sigma = (dist*0.1)/(1.122);
/* The coefficient 1.122 appears because dist is the native distance, i.e, the distance at which the potential is minimum; */
/* The distance used in Gromacs is minimal distance, i.e, the distance at which the potential goes to 0                    */
      fprintf (TopFile, "%4d %4d %4d  %11.9e  %11.9e \n", atom1[i], atom2[i], 1, sigma, epsilon[i]);
  }
  fprintf (TopFile, " \n");


/* Exclusion Directive */
  fprintf (TopFile, "[ exclusions ] \n");
  for(i = 0; i < NPAIR; i++)
  {
      fprintf (TopFile, "%4d %4d \n", atom1[i], atom2[i]);
  }

  /* 4.14.2018: exclusion 2 central atoms of a pseudo dihedral angle */
  fp4.open(input4);
  for (i=0; i<NPSEUDO; i++)
  {
        fp4 >> xi0 >> xi1 >> xi2 >> xi3;
 	fprintf (TopFile, "%4d %4d \n", xi1, xi3);
  }


  fprintf (TopFile, " \n");

/* Final */
  fprintf (TopFile, "[ system ] \n");
  fprintf (TopFile, "; name \n");
  fprintf (TopFile, "  Macromolecule \n");
  fprintf (TopFile, " \n");
  fprintf (TopFile, "[ molecules ] \n");
  fprintf (TopFile, "; name            #molec \n");
  fprintf (TopFile, "  Macromolecule   1 \n ");
      
  fp1.close();
  fp2.close();
  fp3.close();
  fp4.close();
  fclose(TopFile);
 
  return 0;
}
