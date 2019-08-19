//program to compute absorption using ExoCross's algorithm

#include "compute_absorp_exocross_algorithm.h"
#include <stdio.h>
#include <math.h>
#define M_PI 3.14159265358979323846

/*
np.polynomial.hermite.hermgauss(30)
(array([-6.86334529, -6.13827922, -5.53314715, -4.98891897, -4.48305536,
        -4.0039086 , -3.54444387, -3.09997053, -2.66713212, -2.24339147,
        -1.82674114, -1.4155278 , -1.00833827, -0.60392106, -0.20112858,
         0.20112858,  0.60392106,  1.00833827,  1.4155278 ,  1.82674114,
         2.24339147,  2.66713212,  3.09997053,  3.54444387,  4.0039086 ,
         4.48305536,  4.98891897,  5.53314715,  6.13827922,  6.86334529]),
 array([2.90825470e-21, 2.81033360e-17, 2.87860708e-14, 8.10618630e-12,
        9.17858042e-10, 5.10852245e-08, 1.57909489e-06, 2.93872523e-05,
        3.48310124e-04, 2.73792247e-03, 1.47038297e-02, 5.51441769e-02,
        1.46735848e-01, 2.80130931e-01, 3.86394890e-01, 3.86394890e-01,
        2.80130931e-01, 1.46735848e-01, 5.51441769e-02, 1.47038297e-02,
        2.73792247e-03, 3.48310124e-04, 2.93872523e-05, 1.57909489e-06,
        5.10852245e-08, 9.17858042e-10, 8.10618630e-12, 2.87860708e-14,
        2.81033360e-17, 2.90825470e-21]))
*/

double[:] hermgauss_points = [-6.86334529, -6.13827922, -5.53314715, -4.98891897, -4.48305536,
                                -4.0039086 , -3.54444387, -3.09997053, -2.66713212, -2.24339147,
                                -1.82674114, -1.4155278 , -1.00833827, -0.60392106, -0.20112858,
                                 0.20112858,  0.60392106,  1.00833827,  1.4155278 ,  1.82674114,
                                 2.24339147,  2.66713212,  3.09997053,  3.54444387,  4.0039086,
                                 4.48305536,  4.98891897,  5.53314715,  6.13827922,  6.86334529];
double[:] hermgauss_weights = [2.90825470e-21, 2.81033360e-17, 2.87860708e-14, 8.10618630e-12,
                                9.17858042e-10, 5.10852245e-08, 1.57909489e-06, 2.93872523e-05,
                                3.48310124e-04, 2.73792247e-03, 1.47038297e-02, 5.51441769e-02,
                                1.46735848e-01, 2.80130931e-01, 3.86394890e-01, 3.86394890e-01,
                                2.80130931e-01, 1.46735848e-01, 5.51441769e-02, 1.47038297e-02,
                                2.73792247e-03, 3.48310124e-04, 2.93872523e-05, 1.57909489e-06,
                                5.10852245e-08, 9.17858042e-10, 8.10618630e-12, 2.87860708e-14,
                                2.81033360e-17, 2.90825470e-21];
/*                              
double compute_one_wavenum_exocross(double wavenumber, double T, double p, double iso_abundance, double iso_mass, \
                        double Q, double[:] v_ij_star, double[:] a, double[:] elower, double[:] g_upper, \
                        double[:] gamma_p_T)
                        {
                            double absorption = 0;
                            int N = sizeof(a) / sizeof(double);
                            for (int i; i < N; ++i)
                            {
                                double sigma_thermal, coefficient, v, sum = 0;
                                dropper_broad = sqrt(2 * k_B * T / (iso_mass * G_TO_AMU * c * c)) * v_ij_star[i];
                                coefficient = gamma_p_T[i] * log(2) / (M_Pi * sigma_thermal) ;
                                v = (wavenumber - v_ij_star[i]) * log(2) / sigma_thermal;
                                
                                for (int j; j < 30; ++j)
                                {
                                    sum += hermgauss_weights[j] / (pow(v - hermgauss_points[j], 2) + pow(gamma_p_T, 2));
                                }
                                absorption += coefficient * sum;
                            }
                            return absorption;
                        }
*/
double[:] compute_exocross(double[:] v, double T, double p, double iso_abundance, double iso_mass, double Q, \
double[:] v_ij_star, double[:] a, double[:] elower, double[:] g_upper, double[:] gamma_p_T)
    {
        double absorption_cross_section[N] = {};
        int N = sizeof(v) / sizeof(double);
        for (int i; i < N; ++i) 
        {
            double absorption = 0;
            int N = sizeof(a) / sizeof(double);
            for (int i; i < N; ++i)
            {
                double sigma_thermal, coefficient, v, sum = 0;
                dropper_broad = sqrt(2 * k_B * T / (iso_mass * G_TO_AMU * c * c)) * v_ij_star[i];
                coefficient = gamma_p_T[i] * log(2) / (M_Pi * sigma_thermal) ;
                v = (wavenumber - v_ij_star[i]) * log(2) / sigma_thermal;
                
                for (int j; j < 30; ++j)
                {
                    sum += hermgauss_weights[j] / (pow(v - hermgauss_points[j], 2) + pow(gamma_p_T, 2));
                }
                absorption += coefficient * sum;
            }
            absorption_cross_section[i] = absorption;
        }
        return absorption_cross_section;
    }

//if do it with the entire fiole isolated from python
double* compute_exocross(int size_v, double* v, double T, double p, double iso_abundance, double iso_mass, double Q, \
double* v_ij_star, double* a, double* elower, double* g_upper, double* gamma_p_T, double* absorption_cross_section)
    {        
        fp = fopen('', "r"); // read mode
        
        for (int i; i < N; ++i) 
        {
            double absorption = 0;
            int N = sizeof(a) / sizeof(double);
            fp = fopen('', "r"); // read mode
            if (fp == NULL) 
                {
                exit(EXIT_FAILURE);
                }
            
            int num_lines 
            for (int i; i < N; ++i)
            {
                FILE* fp;
                char* line = NULL;
                size_t len = 0;
                ssize_t read;
                
                while ((read = getline(&line, &len, fp)) != -1) 
                {
                    fscanf()
                    for (int j; j < 7; ++i)
                    {
                        char
                    }
                }

                double sigma_thermal, coefficient, v, sum = 0;
                dropper_broad = sqrt(2 * k_B * T / (iso_mass * G_TO_AMU * c * c)) * v_ij_star[i];
                coefficient = gamma_p_T[i] * log(2) / (M_Pi * sigma_thermal) ;
                v = (wavenumber - v_ij_star[i]) * log(2) / sigma_thermal;
                
                for (int j; j < 30; ++j)
                {
                    sum += hermgauss_weights[j] / (pow(v - hermgauss_points[j], 2) + pow(gamma_p_T, 2));
                }
                absorption += coefficient * sum;
            }
            absorption_cross_section[i] = absorption;
        }
        return absorption_cross_section;
    }






                        