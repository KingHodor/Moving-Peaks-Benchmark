/* movpeaks.h --- 10/99 */

/* movpeaks.h
 * Copyright (C) 1999 Juergen Branke.
 * This is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License.
 *
 * This file defines the parameters and procedures used in movpeaks.c */



/***** FUNCTIONS *****/
/* initialize all variables at the beginning of the program */
void init_peaks();
/* evaluation function */
double eval_movpeaks(double *gen);

double dummy_eval(double *gen);

double dummy_eval_fitness(double *gen);

/* whenever this function is called, the peaks are changed */
void change_peaks();
/* free disc space at end of program */
void free_peaks();
double get_avg_error();
/* returns the average error of all evaluation calls so far */
double get_offline_performance();
/* returns offline performance */
int get_number_of_evals(); /* returns the number of evaluations so far */

/* the following basis functions are provided :*/
double constant_basis_func(double *gen);
double five_peak_basis_func(double *gen);


/* the following peak functions are provided: */
double peak_function1(double *gen, int peak_number);
double peak_function_cone(double *gen, int peak_number);
double peak_function_hilly(double *gen, int peak_number);
double peak_function_twin(double *gen, int peak_number);

/* the following functions to change the stepsize over time are provided:*/
void change_stepsize_random();
void change_stepsize_linear();

/* allows to set the basis function used */
extern double (*basis_function)(double *gen);


/* defines the form of a single peak */
extern double (*peak_function)(double *gen, int peak_number);

