#include<iostream>
#include<fstream>
using namespace std;
 
int main()
{
char pdb[20], input1[30], outfile1[30], input2[30], outfile2[30];

        scanf ("%s", pdb);
        sprintf (input1, "BondEnergy.inp");
        sprintf (input2, "pseudo.inp");
        ifstream fp1 (input1), fp2(input2);
        sprintf (outfile1, "%s.npair", pdb);
	sprintf (outfile2, "%s.npseudo", pdb);
        ofstream fp3 (outfile1), fp4(outfile2);
     
    int count1 = 0, count2 = 0;
    string line;
 
    while (getline(fp1, line))
        count1 ++;
  
    while (getline(fp2, line))
	count2 ++;

    fp3 << count1;
    fp4 << count2;
    fp1.close(); 
    fp2.close();
    fp3.close();
    fp4.close();
  return 0;
}
