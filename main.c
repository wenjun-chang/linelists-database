#include <stdio.h>
#include "compute_absorp_exocross_algorithm.h"  /* Include the header here, to obtain the function declaration */

int main(void)
{
    double[:] v = v;
    double T = 1000, p = 0.1, Q = 162879.38910000;
    double[:] v_ij_star, a, elower, g_upper, gamma_p_T;
    
    output = compute_exocross(v, T, p, iso_abundance, iso_mass, Q, v_ij_star, a, elower, g_upper, gamma_p_T);
    printf(output);
    
    int N = sizeof(output) / sizeof(double);
    fp = fopen('/home/toma/Desktop/c_exocross_test.txt',"w");
	for (int i = 0; i < 10; ++i)
    	{
    	fprintf(fp,"%d\n", output[i]);
	    }
	    fclose(fp);
    
    return 0;
}