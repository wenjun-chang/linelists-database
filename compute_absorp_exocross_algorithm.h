//header for the program to compute absorption using ExoCross's algorithm

#ifndef COMPUTE_ABSORP_EXOCROSS_ALGORITHM_H_  /* Include guard */
#define COMPUTE_ABSORP_EXOCROSS_ALGORITHM_H_
/*
double compute_one_wavenum_exocross(double wavenumber, double T, double p, double iso_abundance, double iso_mass, \
                        double Q, double[:] v_ij_star, double[:] a, double[:] elower, double[:] g_upper, \
                        double[:] gamma_p_T);
*/
double[:] compute_exocross(double[:] v, double T, double p, double iso_abundance, double iso_mass, \
                        double Q, double[:] v_ij_star, double[:] a, double[:] elower, double[:] g_upper, \
                        double[:] gamma_p_T);
                        
#endif // COMPUTE_ABSORP_EXOCROSS_ALGORITHM_H_ 