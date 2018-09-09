/* 
Implementation of the 0-1 test for chaos from Gottwald and Melbourne in C++

Details on the method can be seen in:
http://www.maths.usyd.edu.au/u/gottwald/preprints/testforchaos_MPI.pdf

Author:
João Peterson
*/

#include <iostream>
#include <cmath>
#include <vector>
#include "z1test.h"

using namespace std;

double z1test(double x[], int size_x)
{
    int j[size_x];
    for(int i=0; i<size_x; i++){j[i] = i+1;}
    int roundN10 = round(size_x/10);
    double t[roundN10];
    for(int i=0; i<roundN10; i++){t[i] = i+1;}
    double M[roundN10];

    int ns = 300; // number of samples
    double c[ns];
    for(int i=0; i<ns; i++){c[i] = ((double) rand()/(RAND_MAX))*4*atan(1.0);} // random from 0 to pi

    double kcorr[ns];

    double *p, *q;
    p = new double[size_x];
    q = new double[size_x];

    double xcos[size_x], xsin[size_x];

    for(int its=0; its<ns; its++){
        for(int i=0; i<size_x; i++){
            xcos[i] = x[i]*cos(j[i]*c[its]);
            xsin[i] = x[i]*sin(j[i]*c[its]);
        }

        cumsum(xcos, size_x, p);
        cumsum(xsin, size_x, q);

        for(int n=0; n<roundN10; n++){
            double aux[size_x-n+1];
            for(int i=0; i<(size_x-n+1); i++){
                aux[i] = pow((p[i+n+1]-p[i]), 2.0) + pow((q[i+n+1]-q[i]), 2.0);
            }
            M[n] = mean(aux, size_x-n+1) - pow(mean(x, size_x), 2.0)*(1.0-cos(n*c[its]))/(1.0-cos(c[its]));
        }

        kcorr[its] = corr(t, M, roundN10);
    }

    vector<double> kcorr_less;
    vector<double> kcorr_more;

    for(int i=0; i<ns; i++)
    {
        if(c[i] < mean(c, ns))
        {
            kcorr_less.push_back(kcorr[i]);
        }
        else
        {
            kcorr_more.push_back(kcorr[i]);
        }
    }

    double array_k_less[kcorr_less.size()];
    double array_k_more[kcorr_more.size()];

    // checking for oversampling
    if( ( (max(x) - min(x))/(mean_abs_diff(x, size_x)) ) > 10 ||
       median(array_k_less, int(kcorr_less.size())) - median(array_k_more, int(kcorr_more.size())) > 0.5)
    {
        cout << "data is probably oversampled." << endl;
        cout << "Use coarser sampling or reduce the maximum value of c." << endl;
    }


    delete [] p;
    delete [] q;

    return median(kcorr, ns);
}
// ----------------------------------------------------------------------------
double mean_abs_diff(double x[], int size_x)
{
    double y[size_x-1];

    for(int i=0; i<size_x-1; i++)
    {
        y[i] = abs(x[i+1] - x[i]);
    }

    return mean(y, size_x-1);
}
// ----------------------------------------------------------------------------
void cumsum(double x[], int size_x, double result[])
{
    result[0] = x[0];
    for(int i=1; i<size_x; i++){
        result[i] = x[i] + result[i-1];
    }
}
// ----------------------------------------------------------------------------
double sum(double a[], int size_a)
{
    double s = 0;
    for (int i = 0; i < size_a; i++)
    {
        s += a[i];
    }
    return s;
}
// ----------------------------------------------------------------------------
double mean(double a[], int size_a)
{
    return sum(a, size_a) / size_a;
}
// ----------------------------------------------------------------------------
double corr(double x[], double y[], int size_n)
{
    /* pearson coefficient

                    _n
                    \   (x_i - mean(x))(y_i - mean(y))
                    /_1
      r =  ------------------------------------------------------------
            sqrt(sum((x_i - mean(x))^2)) * sqrt(sum((y_i - mean(y))^2))

    */

    double N = 0;
    double D1 = 0;
    double D2 = 0;
    double mean_x = mean(x, size_n);
    double mean_y = mean(y, size_n);

    for(int i=0; i<size_n; i++){
        N += (x[i] - mean_x)*(y[i] - mean_y);
        D1 += pow((x[i] - mean_x), 2.0);
        D2 += pow((y[i] - mean_y), 2.0);
    }
    D1 = sqrt(D1);
    D2 = sqrt(D2);

    return (N/(D1*D2));
}
// ----------------------------------------------------------------------------
double median(double x[], int size_x)
{
    double result;

    my_sort(x, size_x);

    if(size_x % 2 == 0)
        result = (x[(size_x-1)/2] + x[(size_x+1)/2])/2;
    else
        result = x[size_x/2];

    return result;
}
// ----------------------------------------------------------------------------
void my_sort(double x[], int size_x)
{
    double a;
    for(int i=0; i<size_x; i++){
        for(int j=i; j<size_x; j++){
            if(x[j]<x[i]){
                a = x[j];
                x[j] = x[i];
                x[i] = a;
            }
        }
    }
}