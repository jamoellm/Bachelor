//  main.cpp
//  spinner
//
//  Created by Shengkai Li on 1/29/23.
//
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <string>
#include <list>
#include <algorithm>
#include <chrono>
#include "write.h"
#include "init.h"
#include "physics.h"


int main(int argc, const char * argv[]) {
    //const double pi = 3.14159265358979323846;
    //cout << "Hello, World!\n";
    double T=1800;
    int NN=8;
    //int NB[7]={0,6,18,12,24,30,36};
    double tL=0.010; // loop time in seconds
    int LA_array[5]={5,7,9,15,21};
    //int NB[9]={9,1,3,5,7,11,13,15,17};
    for (int i=0; i<10; i++) {
        for (int j=3;j<5;j++) {
            evolution(T,NN,8,tL,LA_array[j],(double)i);
        }
    }
    return 0;
}

void evolution(double T, int NN, int NSS,double tL, int LC, double s0) {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    auto t1 = high_resolution_clock::now();
    const double pi = 3.14159265358979323846;

    srand(s0);
    string fn = "spinner_NN"+to_string(NN)+"spinner_NSS"+to_string(NSS)+"_T"+to_string((int)T)+"_LC"+to_string(LC)+"_s"+to_string((int)s0)+"_v3_0.dat";
    ofstream myfile(fn);
    //myfile << "TEST2 TEST2 !!!";

    double dt=1e-4;
    int NT=ceil(T/dt); int NL=ceil(tL/dt); int NS=100; int L0=ceil(NT/10);
    int LA=10; // acc detec window
    int LD=10; //detection interval
    int LDM=100; // min collision detection
    double std0_threshold=0.02*9.8; // the unit is g in Arduino, so need to convert here to SI
    double om_threshold=30*pi/180; // the unit is deg/s in Arduino, so need to convert here to SI
    //int NM0=round(100/dt); //round(150/dt)
    double R0=0.038; double Ri=0.030;
    double m=0.0525; double I=0.33*m*R0*R0;
    double r_imu=0.002;
    double gamma0=1.5; double Om=28;
    double tau0=I*gamma0;
    double A1=0.38; double A2=0.28;
    double ks=1e4; double kd=2e5; // there is a repeated definition afterwards in functions
    double eta=(6.5*0.01)*m; double etaRot=I*gamma0/Om;
    double kax=2; double kay=2;
    double v0=0.05;
    double xp,yp,ff,d,tt;
    double *x_i = (double *)malloc(NN * sizeof(double));
    double *y_i = (double *)malloc(NN * sizeof(double));
    double *theta_i = (double *)malloc(NN * sizeof(double));
    double *x_f = (double *)malloc(NN * sizeof(double));
    double *y_f = (double *)malloc(NN * sizeof(double));
    double *theta_f = (double *)malloc(NN * sizeof(double));
    double *xt_i = (double *)malloc(NN * sizeof(double));
    double *yt_i = (double *)malloc(NN * sizeof(double));
    double *thetat_i = (double *)malloc(NN * sizeof(double));
    double *xt_f = (double *)malloc(NN * sizeof(double));
    double *yt_f = (double *)malloc(NN * sizeof(double));
    double *thetat_f = (double *)malloc(NN * sizeof(double));
    double *xt1 = (double *)malloc(NN * sizeof(double));
    double *yt1 = (double *)malloc(NN * sizeof(double));
    double *thetat1 = (double *)malloc(NN * sizeof(double));
    double *heap0 = (double *)malloc(NSS*LC * sizeof(double));
    double *sign0 = (double *)malloc(NN * sizeof(double));
    double *a0 = (double *)malloc(NSS*LA * sizeof(double));
    //double *om0 = (double *)malloc(NSS*LA * sizeof(double));
    int *clt = (int *)malloc(NSS * sizeof(int));
    int *cc = (int *)malloc(NSS * sizeof(int));
    int *c = (int *)malloc(NSS * sizeof(int));
    double *om_Old = (double *)malloc(NN * sizeof(double));
    //double dom,rand0,P0;
    int ccw,cw;
    double avg0,std0,r_om,om_imu;
    int id,idc,s_other;

    //xt1,yt1,thetat1
    //double xt1[NN],yt1[NN],thetat1[NN],xtt1[NN],ytt1[NN],thetatt1[NN],xtt2[NN],ytt2[NN],thetatt2[NN];
    double xtt1,ytt1,thetatt1,xtt2,ytt2,thetatt2;

    init_fillZeros(clt, cc, c,
                   x_i, y_i, theta_i,
                   xt_i, yt_i, thetat_i,
                   NN, NSS);

    // Assign the initial positions randomly
    init_placement(x_i, y_i, NN, A1, A2, R0);
    init_velocity_Angle(xt_i, yt_i, theta_i, NN, v0, pi);

    // Assign the initial heap and torque for smart spinners
    init_memoryHeap(heap0, sign0, om_Old, NSS, LC, Om);


    //*(x_i+0)=0.8*R0;  *(y_i+0)=0;  *(theta_i+0)=0.3;
    //*(xt_i+0)=-0.05;  *(yt_i+0)=0; *(thetat_i+0)=30;
    //*(x_i+1)=-0.8*R0; *(y_i+1)=0;  *(theta_i+1)=0;
    //*(xt_i+1)=0.05;   *(yt_i+1)=0; *(thetat_i+1)=30;

    vector<double> x,y,theta,xt,yt,thetat,ss0;
    // main loop
    vector<double> F;
    for (int i=0; i<NT; i++) {
        // v(t+0.5dt)=v(t)+0.5*a(t)*dt
        x.assign(x_i,x_i+NN);
        y.assign(y_i,y_i+NN);
        theta.assign(theta_i,theta_i+NN);
        xt.assign(xt_i,xt_i+NN);
        yt.assign(yt_i,yt_i+NN);
        thetat.assign(thetat_i,thetat_i+NN);
        ss0.assign(sign0,sign0+NN);
        F=spinnerInteraction_v10(x,y,theta,xt,yt,thetat,NN,tau0,R0,A1,A2,m,I,ks,eta,etaRot,kax,kay,ss0);
        for (int j=0;j<NN;j++) {
            xtt1=F[j];
            ytt1=F[j+NN];
            thetatt1=F[j+2*NN];
            *(xt1+j)=*(xt_i+j) + 0.5*dt*xtt1;
            *(yt1+j)=*(yt_i+j) + 0.5*dt*ytt1;
            *(thetat1+j)=*(thetat_i+j) + 0.5*dt*thetatt1;
        }

        // x(t+dt)=x(t)+v(t+0.5dt)*dt
        for (int j=0;j<NN;j++) {
            *(x_f+j)=*(x_i+j)+dt* *(xt1+j);
            *(y_f+j)=*(y_i+j)+dt* *(yt1+j);
            *(theta_f+j)=*(theta_i+j)+dt* *(thetat1+j);
        }

        // a(t+dt)=a(x(t+dt))
        // v(t+dt)=v(t+0.5dt)+0.5*a(t+dt)*dt
        x.assign(x_f,x_f+NN);
        y.assign(y_f,y_f+NN);
        theta.assign(theta_f,theta_f+NN);
        xt.assign(xt1,xt1+NN);
        yt.assign(yt1,yt1+NN);
        thetat.assign(thetat1,thetat1+NN);

        F=spinnerInteraction_v10(x,y,theta,xt,yt,thetat,NN,tau0,R0,A1,A2,m,I,ks,eta,etaRot,kax,kay,ss0);
        for (int j=0; j<NN; j++) {
            xtt2=F[j];
            ytt2=F[j+NN];
            thetatt2=F[j+2*NN];
            *(xt_f+j)=*(xt1+j)+0.5*dt*xtt2;
            *(yt_f+j)=*(yt1+j)+0.5*dt*ytt2;
            *(thetat_f+j)=*(thetat1+j)+0.5*dt*thetatt2;
        }

        // today becomes yesterday
        for (int j=0;j<NN;j++) {
            *(x_i+j)=*(x_f+j);
            *(y_i+j)=*(y_f+j);
            *(theta_i+j)=*(theta_f+j);
            *(xt_i+j)=*(xt_f+j);
            *(yt_i+j)=*(yt_f+j);
            *(thetat_i+j)=*(thetat_f+j);
        }

        // loop on Arduino if the time interval reaches multiples of 0.01 s (or other number)
        for (int j=0; j<NSS; j++) {
            if (i%NL==0) {
                if (*(clt+j)<1) (*(clt+j))=0;
                else (*(clt+j))-=1;

                id=*(c+j)%LA;
                *(c+j)+=1;

                // save acc into acc heap
                xtt2=F[j];
                ytt2=F[j+NN];
                om_imu=*(thetat_f+j);
                *(a0+j*LA+id)=sqrt(xtt2*xtt2+ytt2*ytt2+om_imu*om_imu*r_imu);

                if ((*(c+j))%LD==0) {
                    if (*(clt+j)==0) {
                        // evaluate acc fluctuation
                        avg0=0;
                        for (int k=0; k<LA; k++) avg0+=*(a0+j*LA+k);
                        avg0=avg0/LA;
                        std0=0;
                        for (int k=0; k<LA; k++) std0+=(*(a0+j*LA+k)-avg0)*(*(a0+j*LA+k)-avg0);
                        std0=sqrt(std0/LA);
                        if (std0>std0_threshold) { // if in consistent state
                            r_om=((*(om_Old+j))-(*(thetat_f+j)))/(*(om_Old+j)); //evaluate spin drop
                            if (abs(*(om_Old+j))>om_threshold) { //if the spin is significant
                                // push heap
                                *(cc+j) +=1;
                                idc=(*(cc+j))%LC;
                                if (*(om_Old+j)>0) {
                                    if (r_om<0.3) *(heap0+j*LC+idc)=-1; // the one I met is opposite
                                    else *(heap0+j*LC+idc)=1; // the one I met is same
                                }
                                if (*(om_Old+j)<0) {
                                    if (r_om<0.3) *(heap0+j*LC+idc)=1; // the one I met is opposite
                                    else *(heap0+j*LC+idc)=-1; // the one I met is same
                                }

                                // evaluate heap
                                ccw=0;
                                cw=0;
                                for (int k=0; k<LC; k++) {
                                    if (*(heap0+j*LC+k)==1) ccw++;
                                    else cw++;
                                }
                                if (ccw>cw) s_other = 1;
                                else s_other = -1;

                                // change state
                                if ((s_other* (*(sign0+j)))<0) *(sign0+j) = - *(sign0+j);
                                *(clt+j) = LDM;
                            }
                        }
                    }
                    *(om_Old+j)=*(thetat_f+j);
                }
            }
        }


        if (i%L0==0) std::cout << ".\n";

        if (i%NS==0) {
            writeArray(theta_f, NN, myfile);
            writeArray(x_f, NN, myfile);
            writeArray(y_f, NN, myfile);

            writeArray(thetat_f, NN, myfile);
            writeArray(xt_f, NN, myfile);
            writeArray(yt_f, NN, myfile);

            write2DArray(a0, NSS, LA, myfile);
            myfile << "\n";
            writeArray(sign0, NSS, myfile);
            writeArray(cc, NSS, myfile);
            write2DArray(heap0, NSS, LC, myfile);
            myfile << "\n" << "\n";
        }
    }

    //double *x0 = (double *)malloc(NN * sizeof(double));
    //for (int j=0; j<NN; j++) {
    //    *(x0 + j)=(double)j+1.5;
    //}
    //vector<double> tt;
    //tt.assign(x0,x0+NN);
    //tt.clear();
    //tt.assign(x0,x0+3);
    //vector<double> F;
    //F=spinnerInteraction_v10(tt,tt,tt,tt,tt,tt,NN,NB);
    //vector<double> R,H;
    //R=lineline(0,0,1,1,1,0,0.5,1);
    //H=crosscross_v8(0,0,0,0.05,0.05,0.1,0.065,0.005,0.1,0,0,-0.2);
    myfile.close();
    auto t2 = high_resolution_clock::now();
    auto ms_int = duration_cast<milliseconds>(t2 - t1);
    std::cout << ms_int.count() << "ms\n";
}

