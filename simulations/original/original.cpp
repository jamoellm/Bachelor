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
void evolution(double T, int NN, int NSS,double tL, int LC,double s0);
vector<double> spinnerInteraction_v10(vector<double> x,vector<double> y,vector<double> theta,vector<double> xt,vector<double> yt,vector<double> thetat,int NN,double tau0,double R0,double A1,double A2,double m,double I,double ks,double eta,double etaRot,double kax,double kay,vector<double> sign0);
vector<double> crosscross_v8(double x1,double y1,double theta1,double xt1,double yt1,double thetat1,double x2,double y2,double theta2,double xt2,double yt2,double thetat2);
vector<double> lineline(double p11x,double p11y,double p12x,double p12y,double p21x,double p21y,double p22x,double p22y);
double norm(double x, double y);
double cross(double x1, double y1, double x2, double y2);

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
    string fn = "spinner_NN"+to_string(NN)+"spinner_NSS"+to_string(NSS)+"_T"+to_string((int)T)+"_LC"+to_string(LC)+"_s"+to_string((int)s0)+".dat";
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

    for (int i=0; i<NSS; i++) {
        *(clt + i)=0; // collision time count down
        *(cc + i)=0; // collision counter
        *(c + i)=0;
    }

    //xt1,yt1,thetat1
    //double xt1[NN],yt1[NN],thetat1[NN],xtt1[NN],ytt1[NN],thetatt1[NN],xtt2[NN],ytt2[NN],thetatt2[NN];
    double xtt1,ytt1,thetatt1,xtt2,ytt2,thetatt2;

    // Zero the initial values
    for (int j=0; j<NN; j++) {
        *(x_i+j)=0;   *(y_i+j)=0;   *(theta_i+j)=0;
        *(xt_i+j)=0;  *(yt_i+j)=0;  *(thetat_i+j)=0;
    }

    // Assign the initial positions randomly
    int j=1;
    while (j<NN) {
        xp=2 * (((double)rand() / RAND_MAX) - 0.5)*(A1-R0);
        yp=2 * (((double)rand() / RAND_MAX) - 0.5)*(A2-R0);
        ff=1;
        for (int k=0; k<j; k++) {
            d=norm(xp- *(x_i+k),yp- *(y_i+k));
            if (d<(2*R0)) {
                ff=0;
                break;
            }
        }
        if (ff==1) {
            *(x_i+j)=xp;
            *(y_i+j)=yp;
            j++;
        }
    }
    for (j=0; j<NN; j++) {
        tt=2*pi*(double)rand() / RAND_MAX;
        *(xt_i+j)=v0*cos(tt);
        *(yt_i+j)=v0*sin(tt);
        *(theta_i+j)=2*pi*(double)rand() / RAND_MAX;
    }

    // Assign the initial heap and torque for smart spinners
    int NSS2=NSS/2;
    for (int j=0; j<NSS2; j++) {
        ccw=0;
        cw=0;
        for (int k=0; k<LC; k++) {
            if ((k%2)==0) *(heap0 + j*LC + k)=-1;
            else *(heap0 + j*LC + k)=1;
        }
        for (int k=0; k<LC; k++) {
            if (*(heap0 + j*LC + k)==1) ccw++;
            else cw++;
        }
        if (ccw>cw) {
            *(sign0+j)=1;
            *(om_Old+j)=Om;
        }
        else {
            *(sign0+j)=-1;
            *(om_Old+j)=-Om;
        }
    }

    for (int j=NSS2; j<NSS; j++) { //int j=NSS2; j<NSS;
        ccw=0;
        cw=0;
        for (int k=0; k<LC; k++) {
            if ((k%2)==0) *(heap0 + j*LC + k)=1;
            else *(heap0 + j*LC + k)=-1;
        }
        for (int k=0; k<LC; k++) {
            if (*(heap0 + j*LC + k)==1) ccw++;
            else cw++;
        }
        if (ccw>cw) {
            *(sign0+j)=1;
            *(om_Old+j)=Om;
        }
        else {
            *(sign0+j)=-1;
            *(om_Old+j)=-Om;
        }
    }


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
            for (int j=0;j<NN;j++) {
                myfile << *(theta_f+j) << " ";
            }
            myfile << "\n";
            for (int j=0;j<NN;j++) {
                myfile << *(x_f+j) << " ";
            }
            myfile << "\n";
            for (int j=0;j<NN;j++) {
                myfile << *(y_f+j) << " ";
            }
            myfile << "\n";
            for (int j=0;j<NN;j++) {
                myfile << *(thetat_f+j) << " ";
            }
            myfile << "\n";
            for (int j=0;j<NN;j++) {
                myfile << *(xt_f+j) << " ";
            }
            myfile << "\n";
            for (int j=0;j<NN;j++) {
                myfile << *(yt_f+j) << " ";
            }
            myfile << "\n";
            for (int j=0;j<NSS;j++) {
                for (int k=0;k<LA;k++) {
                    myfile << *(a0+j*LA+k) << " ";
                }
                myfile << "\n";
            }
            myfile << "\n";
            for (int j=0;j<NSS;j++) {
                myfile << *(sign0+j) << " ";
            }
            myfile << "\n";
            for (int j=0;j<NSS;j++) {
                myfile << *(cc+j) << " ";
            }
            myfile << "\n";
            for (int j=0;j<NSS;j++) {
                for (int k=0;k<LC;k++) {
                    myfile << *(heap0+j*LC+k) << " ";
                }
                myfile << "\n";
            }
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

vector<double> spinnerInteraction_v10(vector<double> x,vector<double> y,vector<double> theta,vector<double> xt,vector<double> yt,vector<double> thetat,int NN,double tau0,double R0,double A1,double A2,double m,double I,double ks,double eta,double etaRot,double kax,double kay,vector<double> sign0){
    vector<double> H;
    vector<double> F;
    F.resize(NN*3, 0);
    double Fx[NN],Fy[NN],tau[NN],xtt[NN],ytt[NN],thetatt[NN];
    double d,d1,d2,d3,d4,vp,N,f;
    double kd=2e5; double mu=0.2;

    for (int j=0; j<NN; j++) {
        Fx[j]=0;  Fy[j]=0;  tau[j]=0;
        xtt[j]=0; ytt[j]=0; thetatt[j]=0;
    }
    // driving torque
    for (int j=0; j<NN; j++) {
        tau[j]=tau0*sign0[j];
    }

    // cross-cross interaction between spinner j and k
    for (int j=0; j<NN; j++) {
        for (int k=0; k<j; k++) {
            d=norm(x[k]-x[j],y[k]-y[j])-2*R0;
            if (d<0) {
                H=crosscross_v8(x[j],y[j],theta[j],xt[j],yt[j],thetat[j],x[k],y[k],theta[k],xt[k],yt[k],thetat[k]);
                Fx[j]+=H[0];
                Fy[j]+=H[1];
                tau[j]+=H[2];
                Fx[k]+=H[3];
                Fy[k]+=H[4];
                tau[k]+=H[5];
            }
        }
    }

    for (int j=0; j<NN; j++) {
        // collision with the wall
        // right wall
        d1=x[j]-(A1-R0);
        if (xt[j]>0) vp=xt[j];
        else vp=0;
        if (d1>0) {
            N=-d1*ks*(1+(kd/ks)*vp);
            f=-mu*N;
            Fx[j]+=N;
            if ((yt[j]+R0*thetat[j])>0) Fy[j]-=f;
            else Fy[j]+=f;
            tau[j]+=R0*Fy[j];
        }

        // left wall
        d2=(-A1+R0)-x[j];
        if (xt[j]<0) vp=0-xt[j];
        else vp=0;
        if (d2>0) {
            N=d2*ks*(1+(kd/ks)*vp);
            f=mu*N;
            Fx[j]+=N;
            if ((yt[j]-R0*thetat[j])>0) Fy[j]-=f;
            else Fy[j]+=f;
            tau[j]-=R0*Fy[j];
        }

        // top wall
        d3=y[j]-(A2-R0);
        if (yt[j]>0) vp=yt[j];
        else vp=0;
        if (d3>0) {
            N=-d3*ks*(1+(kd/ks)*vp);
            f=-mu*N;
            Fy[j]+=N;
            if (xt[j]-R0*thetat[j]>0) Fx[j]-=f;
            else Fx[j]+=f;
            tau[j]-=R0*Fx[j];
        }

        // bottom wall
        d4=(-A2+R0)-y[j];
        if (yt[j]<0) vp=-yt[j];
        else vp=0;
        if (d4>0) {
            N=d4*ks*(1+(kd/ks)*vp);
            f=mu*N;
            Fy[j]+=N;
            if (xt[j]+R0*thetat[j]>0) Fx[j]-=f;
            else Fx[j]+=f;
            tau[j]+=R0*Fx[j];
        }

        // translational drag
        Fx[j]-=eta*xt[j];
        Fy[j]-=eta*yt[j];

        // rotational drag
        tau[j]-=etaRot*thetat[j];

        // inward air current force
        Fx[j]-=m*kax*x[j]*x[j]*x[j];
        Fy[j]-=m*kay*y[j]*y[j]*y[j];

        // from force/torque to acceleration/ang. acc.
        xtt[j]=Fx[j]/m;
        ytt[j]=Fy[j]/m;
        thetatt[j]=tau[j]/I;
    }

    for (int j=0; j<NN; j++) {
        F[j]=xtt[j];
        F[j+NN]=ytt[j];
        F[j+2*NN]=thetatt[j];
    }
    return F;
}

vector<double> crosscross_v8(double x1,double y1,double theta1,double xt1,double yt1,double thetat1,double x2,double y2,double theta2,double xt2,double yt2,double thetat2){
    double Ri=0.030;  double R0=0.038;
    double ks=1e4; double kd=2e5;
    const double pi = 3.14159265358979323846;
    double xs1[97],ys1[97],xs2[97],ys2[97];
    vector<double> R;
    vector<double> H;
    H.resize(6, 0);
    // H0&1: f1, H2: tau1, H3&4: f2. H5:tau2
    for (int j=0; j<6; j++) {
        H[j]=0;
    }

    // Find the positions of the key points for both spinners
    double Dtheta=2*pi/24;
    double t1=0;      double t2=Dtheta/2;
    double t,r,dd,p11x,p11y,p12x,p12y,p21x,p21y,p22x,p22y,pcx,pcy,d1x,d1y,d2x,d2y,co1x,co1y,co2x,co2y,vc1x,vc1y,vc2x,vc2y,vx,vy,vp,F0;
    for (int j=0; j<24; j++) {
        t=t1+Dtheta*j;
        r=Ri;
        xs1[j*2]=x1+r*cos(t+theta1);
        ys1[j*2]=y1+r*sin(t+theta1);
        xs2[j*2]=x2+r*cos(t+theta2);
        ys2[j*2]=y2+r*sin(t+theta2);

        t=t2+Dtheta*j;
        r=R0;
        xs1[1+j*2]=x1+r*cos(t+theta1);
        ys1[1+j*2]=y1+r*sin(t+theta1);
        xs2[1+j*2]=x2+r*cos(t+theta2);
        ys2[1+j*2]=y2+r*sin(t+theta2);

    }
    xs1[48]=xs1[0]; ys1[48]=ys1[0];
    xs2[48]=xs2[0]; ys2[48]=ys2[0];

    for (int j=0; j<48; j++) {
        for (int k=0; k<48; k++) {
            p11x=xs1[j];    p11y=ys1[j];
            p12x=xs1[j+1];  p12y=ys1[j+1];
            p21x=xs2[k];    p21y=ys2[k];
            p22x=xs2[k+1];  p22y=ys2[k+1];
            dd=norm(0.5*(p11x+p12x)-0.5*(p21x+p22x),0.5*(p11y+p12y)-0.5*(p21y+p22y));
            if (dd<0.05) {
                R=lineline(p11x,p11y,p12x,p12y,p21x,p21y,p22x,p22y);
                d1x=R[0]; d1y=R[1];
                d2x=R[2]; d2y=R[3];
                pcx=R[4]; pcy=R[5];

                if (R[6]==1) {
                    co1x=pcx-x1; co1y=pcy-y1;
                    co2x=pcx-x2; co2y=pcy-y2;
                    vc1x=-thetat1*co1y+xt1; vc1y=thetat1*co1x+yt1;
                    vc2x=-thetat2*co2y+xt2; vc2y=thetat2*co2x+yt2;
                    vx=vc1x-vc2x;
                    vy=vc1y-vc2y;
                    vp=-(vx*d1x+vy*d1y)/norm(d1x,d1y);
                    if (vp<0) vp=0;
                    F0=ks*(1+(kd/ks)*vp);
                    H[0]+=F0*d1x; // f1x
                    H[1]+=F0*d1y; // f1y
                    H[3]+=F0*d2x; // f2x
                    H[4]+=F0*d2y; // f2y

                    H[2]+=F0*cross(co1x,co1y,d1x,d1y); // tau1
                    H[5]+=F0*cross(co2x,co2y,d2x,d2y); // tau2
                }
            }
        }
    }
    return H;
}

vector<double> lineline(double p11x,double p11y,double p12x,double p12y,double p21x,double p21y,double p22x,double p22y){
    vector<double> R;
    R.resize(7, 0);
    // R0&1:d1, R2&3:d2, R4&5:pc, R6:collision or not
    double a=p12y; double b=p12x; double c=p11y; double d=p11x;
    double e=p22y; double f=p22x; double g=p21y; double h=p21x;
    double pcx=(d*(-f*g+a*(f-h)+e*h)+b*(f*g-e*h+c*(-f+h)))/(d*(e-g)+b*(-e+g)+(a-c)*(f-h));
    double pcy=(a*(d*e-d*g-e*h+f*g)+b*c*(g-e)+c*e*h-c*f*g)/((a-c)*(f-h)+b*(g-e)+d*(e-g));
    double d1=(pcx-p11x)*(pcx-p12x)+(pcy-p11y)*(pcy-p12y);
    double d2=(pcx-p21x)*(pcx-p22x)+(pcy-p21y)*(pcy-p22y);
    double rx,ry,ex,ey,d0,d1x,d1y,d2x,d2y;
    d1x=0; d1y=0; d2x=0; d2y=0;
    if (!((d1>0)||(d2>0))) {
        R[6]=1;
        double dd1=norm(pcx-p11x,pcy-p11y);
        double dd2=norm(pcx-p12x,pcy-p12y);
        double dd3=norm(pcx-p21x,pcy-p21y);
        double dd4=norm(pcx-p22x,pcy-p22y);
        double M=min({dd1,dd2,dd3,dd4});
        if (dd1==M) {
            rx=p11x-pcx;
            ry=p11y-pcy;
            ex=p22x-p21x;
            ey=p22y-p21y;
            d0=rx*ex+ry*ey;
            d2x=rx-d0*ex/pow(norm(ex,ey),2);
            d2y=ry-d0*ey/pow(norm(ex,ey),2);
            d1x=-d2x;
            d1y=-d2y;
        }
        if (dd2==M) {
            rx=p12x-pcx;
            ry=p12y-pcy;
            ex=p22x-p21x;
            ey=p22y-p21y;
            d0=rx*ex+ry*ey;
            d2x=rx-d0*ex/pow(norm(ex,ey),2);
            d2y=ry-d0*ey/pow(norm(ex,ey),2);
            d1x=-d2x;
            d1y=-d2y;
        }
        if (dd3==M) {
            rx=p21x-pcx;
            ry=p21y-pcy;
            ex=p12x-p11x;
            ey=p12y-p11y;
            d0=rx*ex+ry*ey;
            d1x=rx-d0*ex/pow(norm(ex,ey),2);
            d1y=ry-d0*ey/pow(norm(ex,ey),2);
            d2x=-d1x;
            d2y=-d1y;
        }
        if (dd4==M) {
            rx=p22x-pcx;
            ry=p22y-pcy;
            ex=p12x-p11x;
            ey=p12y-p11y;
            d0=rx*ex+ry*ey;
            d1x=rx-d0*ex/pow(norm(ex,ey),2);
            d1y=ry-d0*ey/pow(norm(ex,ey),2);
            d2x=-d1x;
            d2y=-d1y;
        }
    }
    R[0]=d1x; R[1]=d1y; R[2]=d2x; R[3]=d2y;
    R[4]=pcx; R[5]=pcy;
    return R;
};

double norm(double x, double y)
{
    double n = sqrt(x*x + y*y);
    return n;
}

double cross(double x1, double y1, double x2, double y2)
{
    double c = x1*y2 - x2*y1;
    return c;
}

/*
int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}
 */
