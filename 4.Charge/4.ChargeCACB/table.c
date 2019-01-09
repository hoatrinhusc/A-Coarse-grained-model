/* This code create the tabulated table for the potential and the force of exp(-r/length)*rinv.  the q_iq_j/4piepsilon is included in the gmx already. may 10 2018 */
/* http://www.gromacs.org/@api/deki/files/94/=gromacs_nb.pdf */

#include <stdio.h>
#include <math.h>

int main()
{
    FILE    *fout;
    double  r, rinv;
    double length;

    fout = fopen("table_Debye.xvg", "w");
    fprintf(fout, "#\n# Debye Huckel potential\n#\n");
 
    /* kB = 0.083; T = 300; epsilon0 = 80; I = 0.15M */
    length = sqrt(80*0.008314*300/(4*3.14*138.9*2*0.6022*0.15));
    /* Reduced length */
    length = length/38.;


    for (r=0; r<=30; r+=0.002) {
	rinv = 1/r;
        double f = exp(-r/length)*rinv;
	/* Minus derivative of f with respect to r */
        double fprime = f*(rinv+1/length);

        /* print output */
        if (r<0.04) {
            fprintf(fout, "%12.10e   %12.10e %12.10e   %12.10e %12.10e   %12.10e %12.10e\n", r,0.0,0.0,0.0,0.0,0.0,0.0);
        } else {
            fprintf(fout, "%12.10e   %12.10e %12.10e   %12.10e %12.10e   %12.10e %12.10e\n", r,f,fprime,0,0,0,0);
        }
    }

    fclose(fout);
    return(0);
}
