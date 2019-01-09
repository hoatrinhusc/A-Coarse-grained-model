#include<iostream>
#include<fstream>
using namespace std;
 
int main()
{
char pdb[20], input[30], outfile[30];

        scanf ("%s", pdb);
        sprintf (input, "%s.pdb", pdb);
        ifstream fp1 (input);
        sprintf(outfile,"%s.numat",pdb);
        ofstream fp2 (outfile);
     
    int count = 0;
    string line;
 
    while (getline(fp1, line))
        count++;
  
    fp2 << count; 
  return 0;
}
