/* Moving Peaks Function --- 10/99 */



/* movpeaks.c
 * Copyright (C) 1999 Juergen Branke.
 * This is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License.
 *
 * This module generates the Moving Peaks Evaluation Function, a
 * dynamic benchmark problem changing over time.
 */

#include <stdlib.h>
#include <math.h>
#include "movpeaks.h"
#include "Algorithm.h"
/***** PARAMETER SETTINGS *****/



int change_frequency = 0;
/* number of evaluations between changes. change_frequency
		       =0 means that function never changes (or only if function change_peaks is called)*/

unsigned long int movrandseed = 1;
/* seed for built-in random number generator */

int geno_size = 5;
/* number of dimensions, or the number of double
                          valued genes */
double vlength = 1.0;
/* distance by which the peaks are moved, severity */

double height_severity = 7.0;
/* severity of height changes, larger numbers
                               mean larger severity */
double width_severity = 1.0; /* severity of width changes, larger numbers
                               mean larger severity */

/* lambda determines whether there is a direction of the movement, or whether
   they are totally random. For lambda = 1.0 each move has the same direction,
   while for lambda = 0.0, each move has a random direction */
double lambda = 0.5;

int number_of_peaks = 10;
/* number of peaks in the landscape */

int use_basis_function = 0;
/* if set to 1, a static landscape (basis_function) is included in the fitness evaluation */

int calculate_average_error = 1;
/* saves computation time if not needed and set to 0 */
int calculate_offline_performance = 1;
/* saves computation time if not needed and set to 0 */
int calculate_right_peak = 1; /* saves computation time if not needed and set to 0 */

/* minimum and maximum coordinate in each dimension */
double mincoordinate = 0.0, maxcoordinate = 100.0;
/* minimum and maximum height of the peaks          */
/* height chosen randomly when standardheight = 0.0 */
double minheight = 30.0, maxheight = 70.0, standardheight = 50.0;
/* width chosen randomly when standardwidth = 0.0 */
double minwidth = 1.0, maxwidth = 12.0, standardwidth = 1.0;



/* Functions */
/* evaluation function */
//double eval_movpeaks(double *gen);

/* the following basis functions are provided :*/
//double constant_basis_func(double *gen);
//double five_peak_basis_func(double *gen);
/* the following peak functions are provided: */
//double peak_function1(double *gen, int peak_number);
//double peak_function_cone(double *gen, int peak_number);
//double peak_function_hilly(double *gen, int peak_number);
//double peak_function_twin(double *gen, int peak_number);


/* allows to set the basis function used */
double (*basis_function)(double *gen) = constant_basis_func;
/* defines the form of a single peak */
double (*peak_function)(double *gen, int peak_number) = peak_function_cone;
//double (*peak_function)(double *gen, int peak_number) = peak_function1;


/****** END OF PARAMETER SECTION *******************/

//void change_peaks();
/* preliminary declaration of function change_peaks()*/
int recent_change = 1;
/* indicates that a change has just ocurred */
int current_peak;
/* peak on which the current best individual is located */
int maximum_peak;
/* number of highest peak */
double current_maximum;
/* fitness value of currently best individual */
double offline_performance = 0.0;
double offline_error = 0.0;
double avg_error = 0;
/* average error so far */
double current_error = 0;
/* error of the currently best individual */
double global_max;
/* absolute maximum in the fitness landscape */
int evals = 0;
/* number of evaluations so far */
double **peak;
/* data structure to store peak data */
double *shift;
double *coordinates;
int *covered_peaks;
/* which peaks are covered by the population ? */

double **prev_movement;
/* to store every peak's previous movement */


double movrand();
double movnrand();



/* initialize all variables at the beginning of the program */
void init_peaks() {
    int i, j;
    double dummy;

    shift = (double *) calloc(geno_size, sizeof(double));
    coordinates = (double *) calloc(geno_size, sizeof(double));
    covered_peaks = (int *) calloc(number_of_peaks, sizeof(int));
    peak = (double **) calloc(number_of_peaks, sizeof(double *));
    prev_movement = (double **) calloc(number_of_peaks, sizeof(double *));
    for (i = 0; i < number_of_peaks; i++) {
        peak[i] = (double *) calloc(geno_size + 2, sizeof(double));
        prev_movement[i] = (double *) calloc(geno_size, sizeof(double));
    }
    for (i = 0; i < number_of_peaks; i++)
        for (j = 0; j < geno_size; j++) {
// old      peak[i][j] = 100.0*movrand();
            peak[i][j] = (maxcoordinate - mincoordinate) * movrand() + mincoordinate;
            prev_movement[i][j] = movrand() - 0.5;
        }
    if (standardheight <= 0.0)
        for (i = 0; i < number_of_peaks; i++)
            peak[i][geno_size + 1] = (maxheight - minheight) * movrand() + minheight;
    else
        for (i = 0; i < number_of_peaks; i++)
            peak[i][geno_size + 1] = standardheight;
    if (standardwidth <= 0.0)
        for (i = 0; i < number_of_peaks; i++)
            peak[i][geno_size] = (maxwidth - minwidth) * movrand() + minwidth;
    else
        for (i = 0; i < number_of_peaks; i++)
            peak[i][geno_size] = standardwidth;
    if (calculate_average_error) {
        global_max = -100000.0;
        for (i = 0; i < number_of_peaks; i++) {
            for (j = 0; j < geno_size; j++)
                coordinates[j] = peak[i][j];
            dummy = dummy_eval(coordinates);
            if (dummy > global_max)
                global_max = dummy;
        }
    }
}


/* free disc space at end of program */
void free_peaks() {
    int i;

    for (i = 0; i < number_of_peaks; i++) {
        free(peak[i]);
        free(prev_movement[i]);
    }
}

/* current_peak_calc determines the peak of the current best individual */
void current_peak_calc(double *gen) {
    int i;
    double maximum = -100000.0, dummy;

    current_peak = 0;
    maximum = peak_function(gen, 0);
    for (i = 1; i < number_of_peaks; i++) {
        dummy = peak_function(gen, i);
        if (dummy > maximum) {
            maximum = dummy;
            current_peak = i;
        }
    }
}


/* evaluation function */
double eval_movpeaks(double *gen) {
    int i;
    double maximum = -100000.0, dummy;

    if ((change_frequency > 0) && (evals % change_frequency == 0))
        change_peaks();

    for (i = 0; i < number_of_peaks; i++) {
        dummy = peak_function(gen, i);
        if (dummy > maximum)
            maximum = dummy;
    }

    if (use_basis_function) {

        dummy = basis_function(gen);
        /* If value of basis function is higher return it */
        if (maximum < dummy)
            maximum = dummy;
    }
    if (calculate_average_error) {
        avg_error += global_max - maximum;
    }
    if (calculate_offline_performance) {
        if (recent_change || (maximum > current_maximum)) {
            current_error = global_max - maximum;
            if (calculate_right_peak)
                current_peak_calc(gen);
            current_maximum = maximum;
            recent_change = 0;
        }
        offline_performance += current_maximum;
        offline_error += current_error;
    }
    evals++;     /* increase the number of evaluations by one */
    return (maximum);
}

/* dummy evaluation function allows to evaluate without being counted */
double dummy_eval(double *gen) {
    int i;
    double maximum = -100000.0, dummy;

    for (i = 0; i < number_of_peaks; i++) {
        dummy = peak_function(gen, i);
        if (dummy > maximum)
            maximum = dummy;
    }

    if (use_basis_function) {

        dummy = basis_function(gen);
        /* If value of basis function is higher return it */
        if (maximum < dummy)
            maximum = dummy;
    }

    return (maximum);
}

/* dummy evaluation function allows to evaluate without being counted */
double dummy_eval_fitness(double *gen) {
    int i;
    double maximum = -100000.0, dummy;

    for (i = 0; i < number_of_peaks; i++) {
        dummy = peak_function(gen, i);
        if (dummy > maximum)
            maximum = dummy;
    }

    if (use_basis_function) {

        dummy = basis_function(gen);
        /* If value of basis function is higher return it */
        if (maximum < dummy)
            maximum = dummy;
    }
    Algorithm::FitnessEvaluationSize++;
    return (maximum);
}
/* simple random number generator solely for the test function */
/* movrand creates random number between 0 and 1 */
/* This RNG is taken from the book by Kernighan/Ritchie, maybe it would */
/* be worth to try a better one. */

double movrand() {
    /*  static unsigned long int next;*/
    movrandseed = movrandseed * 1103515245 + 12345;
    return (double) ((unsigned int) (movrandseed / 65536) % 32768) / 32767;
}

/* this function produces normally distributed random values */
double movnrand() {
    static int backup = 0;
    static double x2;
    double x1, w;

    if (backup) {
        backup = 0;
        return (x2);
    } else {
        do {
            x1 = 2.0 * movrand() - 1.0;
            x2 = 2.0 * movrand() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);
        w = sqrt((-2.0 * log(w)) / w);
        x2 = w * x2;
        backup = 1;
        return (x1 * w);
    }
}


/* whenever this function is called, the peaks are changed */
void change_peaks() {
    int i, j;
    double sum, sum2, offset, dummy;

    for (i = 0; i < number_of_peaks; i++) {
        /* shift peak locations */
        sum = 0.0;
        for (j = 0; j < geno_size; j++) {
            shift[j] = movrand() - 0.5;
            sum += shift[j] * shift[j];
        }
        if (sum > 0.0)
            sum = vlength / sqrt(sum);
        else                           /* only in case of rounding errors */
            sum = 0.0;
        sum2 = 0.0;
        for (j = 0; j < geno_size; j++) {
            shift[j] = sum * (1.0 - lambda) * shift[j] + lambda * prev_movement[i][j];
            sum2 += shift[j] * shift[j];
        }
        if (sum2 > 0.0)
            sum2 = vlength / sqrt(sum2);
        else                           /* only in case of rounding errors */
            sum2 = 0.0;
        for (j = 0; j < geno_size; j++) {
            shift[j] *= sum2;
            prev_movement[i][j] = shift[j];
            if (((peak[i][j] + prev_movement[i][j]) < mincoordinate) ||
                ((peak[i][j] + prev_movement[i][j]) > mincoordinate)) {
                dummy = (peak[i][j] + prev_movement[i][j] - mincoordinate) / (maxcoordinate - mincoordinate);
                if (((int) dummy % 2) == 0) {
                    dummy = fabs(dummy - floor(dummy));
                } else {
                    dummy = 1 - fabs(dummy - floor(dummy));
                }
                peak[i][j] = mincoordinate + (maxcoordinate - mincoordinate) * dummy;
            }
            else
                peak[i][j] += prev_movement[i][j];
        }

        /* change peak width */
        j = geno_size;
        offset = movnrand() * width_severity;
        if ((peak[i][j] + offset) < minwidth)
            peak[i][j] = 2.0 * minwidth - peak[i][j] - offset;
        else if ((peak[i][j] + offset) > maxwidth)
            peak[i][j] = 2.0 * maxwidth - peak[i][j] - offset;
        else
            peak[i][j] += offset;
        /* change peak height */
        j++;
        offset = height_severity * movnrand();
        if ((peak[i][j] + offset) < minheight)
            peak[i][j] = 2.0 * minheight - peak[i][j] - offset;
        else if ((peak[i][j] + offset) > maxheight)
            peak[i][j] = 2.0 * maxheight - peak[i][j] - offset;
        else
            peak[i][j] += offset;
    }
    if (calculate_average_error) {
        global_max = -100000.0;
        for (i = 0; i < number_of_peaks; i++) {
            for (j = 0; j < geno_size; j++)
                coordinates[j] = peak[i][j];
            dummy = dummy_eval(coordinates);
            if (dummy > global_max) {
                global_max = dummy;
                maximum_peak = i;
            }
        }
    }
    recent_change = 1;
}


/* Basis Functions */

/* This gives a constant value back to the eval-function that chooses the max of them */
double constant_basis_func(double *gen) {
    return 0.0;
}

double five_peak_basis_func(double *gen) {
    int i, j;
    double maximum = -100000.0, dummy;
    static double basis_peak[5][7] =
            {
                    8.0, 64.0, 67.0, 55.0, 4.0, 0.1, 50.0,
                    50.0, 13.0, 76.0, 15.0, 7.0, 0.1, 50.0,
                    9.0, 19.0, 27.0, 67.0, 24.0, 0.1, 50.0,
                    66.0, 87.0, 65.0, 19.0, 43.0, 0.1, 50.0,
                    76.0, 32.0, 43.0, 54.0, 65.0, 0.1, 50.0,
            };
    for (i = 0; i < 5; i++) {
        dummy = (gen[0] - basis_peak[i][0]) * (gen[0] - basis_peak[i][0]);
        for (j = 1; j < geno_size; j++)
            dummy += (gen[j] - basis_peak[i][j]) * (gen[j] - basis_peak[i][j]);
        dummy = basis_peak[i][geno_size + 1] - (basis_peak[i][geno_size] * dummy);
        if (dummy > maximum)
            maximum = dummy;
    }
    return maximum;
}



/* Peak Functions */

/* sharp peaks */
double peak_function1(double *gen, int peak_number) {
    int j;
    double dummy;

    dummy = (gen[0] - peak[peak_number][0]) * (gen[0] - peak[peak_number][0]);
    for (j = 1; j < geno_size; j++)
        dummy += (gen[j] - peak[peak_number][j]) * (gen[j] - peak[peak_number][j]);
    return peak[peak_number][geno_size + 1] / (1 + (peak[peak_number][geno_size]) * dummy);
}

double peak_function_cone(double *gen, int peak_number) {
    int j;
    double dummy;

    dummy = (gen[0] - peak[peak_number][0]) * (gen[0] - peak[peak_number][0]);
    for (j = 1; j < geno_size; j++)
        dummy += (gen[j] - peak[peak_number][j]) * (gen[j] - peak[peak_number][j]);
    return peak[peak_number][geno_size + 1] - (peak[peak_number][geno_size] * sqrt(dummy));
}

double peak_function_hilly(double *gen, int peak_number) {
    int j;
    double dummy;

    dummy = (gen[0] - peak[peak_number][0]) * (gen[0] - peak[peak_number][0]);
    for (j = 1; j < geno_size; j++)
        dummy += (gen[j] - peak[peak_number][j]) * (gen[j] - peak[peak_number][j]);
    return peak[peak_number][geno_size + 1] - (peak[peak_number][geno_size] * dummy) - 0.01 * sin(20.0 * dummy);
}

double peak_function_twin(double *gen, int peak_number) /* two twin peaks moving together */
{
    int j;
    double maximum = -100000.0, dummy;
    static double twin_peak[7] = /* difference to first peak */
            {
                    1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0,
            };

    dummy = pow(gen[0] - peak[peak_number][0], 2);
    for (j = 1; j < geno_size; j++)
        dummy += pow(gen[j] - peak[peak_number][j], 2);
    dummy = peak[peak_number][geno_size + 1] - (peak[peak_number][geno_size] * dummy);
    maximum = dummy;
    dummy = pow(gen[j] - (peak[peak_number][0] + twin_peak[0]), 2);
    for (j = 1; j < geno_size; j++)
        dummy += pow(gen[j] - (peak[peak_number][j] + twin_peak[0]), 2);
    dummy = peak[peak_number][geno_size + 1] + twin_peak[geno_size + 1] -
            ((peak[peak_number][geno_size] + twin_peak[geno_size]) * dummy);
    if (dummy > maximum)
        maximum = dummy;

    return maximum;
}

/* The following procedures may be used to change the step size over time */


void change_stepsize_random() /* assigns vlength a value from a normal distribution */
{
    vlength = movnrand();
}

void change_stepsize_linear() /* sinusoidal change of the stepsize, */
{
    static int counter = 1;
    static double frequency = 3.14159 / 20.0;  /* returns to same value after 20 changes */

    vlength = 1 + sin((double) counter * frequency);
    counter++;
}

double get_avg_error() /* returns the average error of all evaluation calls so far */
{
    return (avg_error / (double) evals);
}

double get_current_error() /* returns the error of the best individual evaluated since last change */
/* To use this function, calculate_average_error and calculate_offline_performance must be set */
{
    return current_error;
}

double get_offline_performance() /* returns offline performance */
{
    return (offline_performance / (double) evals);
}

double get_offline_error() /* returns offline error */
{
    return (offline_error / (double) evals);
}


int get_number_of_evals() /* returns the number of evaluations so far */
{
    return evals;
}

int get_right_peak()  /* returns 1 if current best individual is on highest peak, 0 otherwise */
{
    if (current_peak == maximum_peak)
        return 1;
    else
        return 0;
}

