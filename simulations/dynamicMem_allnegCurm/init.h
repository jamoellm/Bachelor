// init.h
#pragma once
void init_placement(int oppositeSpinners, double* x_i, double* y_i,
                    int NN, double A1, double A2, double R0);

void init_placement(double* x_i, double* y_i,
                    int NN, double A1, double A2, double R0);

void init_velocity_Angle(double* xt_i, double* yt_i,
                         double* theta_i,
                         int NN, double v0, double pi);

void init_fillZeros(int* clt, int* cc, int* c,
                    double* x_i,  double* y_i,  double* theta_i,
                    double* xt_i, double* yt_i, double* thetat_i,
                    int NN);

void init_memoryHeap(double* heap0, double* sign0, double* om_Old,
                     int NN, int LC, double Om);
