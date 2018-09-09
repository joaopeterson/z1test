#ifndef Z1TEST_INCLUDED
#define Z1TEST_INCLUDED

#include <iostream> // cout
#include <cmath>    // round
#include <vector>

using namespace std;

double z1test(double x[], int size_x);

double mean_abs_diff(double x[], int size_x);

void cumsum(double x[], int size_x, double result[]);

double sum(double a[], int size_a);

double mean(double a[], int size_a);

double corr(double x[], double y[], int size_n);

double median(double x[], int size_x);

void my_sort(double x[], int size_x);

#endif // Z1TEST_INCLUDED