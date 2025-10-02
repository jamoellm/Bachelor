// physics.h
#pragma once
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>

void evolution(double T,double tL, int LC, double s0, int numberUpscaler, double* lambdas, std::string filename);
void spinnerInteraction_v10(double* x,double* y,double* theta,double* xt,double* yt,double* thetat,int NN,double tau0,double R0,double A1,double A2,double m,double I,double ks,double eta,double etaRot,double kax,double kay,std::vector<double> sign0, double* F, double* H, double* R);
void crosscross_v8(double x1,double y1,double theta1,double xt1,double yt1,double thetat1,double x2,double y2,double theta2,double xt2,double yt2,double thetat2, double* H, double* R);
void lineline(double p11x,double p11y,double p12x,double p12y,double p21x,double p21y,double p22x,double p22y, double* R);
double norm(double x, double y);
double cross(double x1, double y1, double x2, double y2);

double weighting(double memoryCell, int position, int personality);
void influenceDecay_v2(double* influences, int NSS, double decayConstant);