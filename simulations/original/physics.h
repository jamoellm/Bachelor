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
using namespace std;

void evolution(double T, int NN, int NSS,double tL, int LC,double s0);
vector<double> spinnerInteraction_v10(vector<double> x,vector<double> y,vector<double> theta,vector<double> xt,vector<double> yt,vector<double> thetat,int NN,double tau0,double R0,double A1,double A2,double m,double I,double ks,double eta,double etaRot,double kax,double kay,vector<double> sign0);
vector<double> crosscross_v8(double x1,double y1,double theta1,double xt1,double yt1,double thetat1,double x2,double y2,double theta2,double xt2,double yt2,double thetat2);
vector<double> lineline(double p11x,double p11y,double p12x,double p12y,double p21x,double p21y,double p22x,double p22y);
double norm(double x, double y);
double cross(double x1, double y1, double x2, double y2);


/*
int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}
 */
