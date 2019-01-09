/* Pre-process: Produce the pdb file of only CAs and CBs */

#include "common.h"

int main()
{
  double cb_x[NUMCB], cb_y[NUMCB], cb_z[NUMCB]; /*coordinates x, y, z of CB beads*/
  string dummy, type[NATOM], amin[NATOM], chain[NATOM];
  int atomindex[NATOM], num[NATOM], k, i, j, renum;
  double x[NATOM], y[NATOM], z[NATOM], dummy1;
  char pdb[20], input1[30], input2[30], input3[30], outfile1[30];
  k = 0, j = 0;

        scanf ("%s", pdb);
        sprintf (input1, "%s.pdb.CBcrd", pdb);
        sprintf (input2, "%s.pdb", pdb);
        ifstream File1(input1), File2(input2);
        sprintf(outfile1,"%s.ab.pdb",pdb);
	ofstream fp1 (outfile1);

  if (File1.fail()||File2.fail())
  {
    cout << "Cannot open input files" << endl;
    return 1;
  }

/* Read PDB file */   
  for(i = 0; i < NATOM; i++)
  {
     File2 >> dummy >> atomindex[i] >> type[i] >> amin[i] >> chain[i] >> num[i] >> x[i] >> y[i] >> z[i] >> dummy1 >> dummy1 >> dummy;
  }

/* Renum residue index to start from 1 */
  renum = num[0]-1;

/* Read coordinates of COM of side chain */
  for(i = 0; i < NUMCB; i++)
  {
   File1 >> cb_x[i] >> cb_y[i] >> cb_z[i];
  }

/* Replace the coordinates of CBs in the original PDB file with COM of sidechains */
  for(i = 0; i < NATOM; i++)
  {
    if (type[i].compare("CB")==0)
    {
      x[i] = cb_x[k];
      y[i] = cb_y[k];
      z[i] = cb_z[k];
      k++;
    }
  }

/* Write output of CAs and CBs only in pdb format, reduced unit of length */
  for(i = 0; i < NATOM; i++)
  {
    if(type[i].compare("CA")==0 || type[i].compare("CB")==0)
    {
      j ++;
      fp1 << "ATOM" << setw(7) << j << setw(5) << type[i]  << " " << amin[i] << setw(2) << chain[i] << setw(4) << num[i]-renum << setw(12) << fixed << setprecision(3) << x[i]/3.8 << std::setw(8) << fixed << setprecision(3) << y[i]/3.8 << setw(8) << fixed << setprecision(3) << z[i]/3.8 << "\n";
    }
  }

  File1.close();
  File2.close();
  fp1.close();

return 0;
}
 



