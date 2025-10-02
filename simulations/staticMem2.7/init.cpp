#include "physics.h"

void init_placement(int oppositeSpinners, double* x_i, double* y_i,
    int NN, double A1, double A2, double R0) {
    double xp, yp, ff, d;
    int j0 = oppositeSpinners == 0 ? 1 : 0;
    int j = j0;

    while (j < oppositeSpinners) {
        xp = (0.8 * ((double)rand() / RAND_MAX) + 0.1) * (A1 - R0);
        yp = (0.8 * ((double)rand() / RAND_MAX) + 0.1) * (A2 - R0);
        ff = 1;
        for (int k = 0; k < j; k++) {
            d = norm(xp - x_i[k], yp - y_i[k]);
            if (d < (2 * R0)) {
                ff = 0;
                break;
            }
        }
        if (ff == 1) {
            x_i[j] = xp;
            y_i[j] = yp;
            j++;
        }
    }

    j = j0;
    while (j < oppositeSpinners) {
        xp = -(0.8 * ((double)rand() / RAND_MAX) + 0.1) * (A1 - R0);
        yp = -(0.8 * ((double)rand() / RAND_MAX) + 0.1) * (A2 - R0);
        ff = 1;
        for (int k = oppositeSpinners; k > oppositeSpinners - j; k--) {
            d = norm(xp - x_i[NN - k], yp - y_i[NN - k]);
            if (d < (2 * R0)) {
                ff = 0;
                break;
            }
        }
        if (ff == 1) {
            x_i[NN - oppositeSpinners + j] = xp;
            y_i[NN - oppositeSpinners + j] = yp;
            j++;
        }
    }

    j = oppositeSpinners == 0 ? 1 : oppositeSpinners;
    int upperBound = oppositeSpinners == 0 ? NN : (NN - oppositeSpinners);
    while (j < upperBound) {
        xp = 2 * (((double)rand() / RAND_MAX) - 0.5) * (A1 - R0);
        yp = 2 * (((double)rand() / RAND_MAX) - 0.5) * (A2 - R0);
        ff = 1;
        for (int k = 0; k < j; k++) {
            d = norm(xp - x_i[k], yp - y_i[k]);
            if (d < (2 * R0)) {
                ff = 0;
                break;
            }
        }
        if (ff == 1) {
            x_i[j] = xp;
            y_i[j] = yp;
            j++;
        }
    }
}

void init_placement(double* x_i, double* y_i,
    int NN, double A1, double A2, double R0) {
    double xp, yp, ff, d;

    int j = 1;
    while (j < NN) {
        xp = 2 * (((double)rand() / RAND_MAX) - 0.5) * (A1 - R0);
        yp = 2 * (((double)rand() / RAND_MAX) - 0.5) * (A2 - R0);
        ff = 1;
        for (int k = 0; k < j; k++) {
            d = norm(xp - x_i[k], yp - y_i[k]);
            if (d < (2 * R0)) {
                ff = 0;
                break;
            }
        }
        if (ff == 1) {
            x_i[j] = xp;
            y_i[j] = yp;
            j++;
        }
    }
}

void init_velocity_Angle(double* xt_i, double* yt_i,
    double* theta_i,
    int NN, double v0, double pi) {
    double tt;
    for (int j = 0; j < NN; j++) {
        tt = 2 * pi * (double)rand() / RAND_MAX;
        xt_i[j] = v0 * cos(tt);
        yt_i[j] = v0 * sin(tt);
        theta_i[j] = 2 * pi * (double)rand() / RAND_MAX;
    }
}

void init_fillZeros(int* clt, int* cc, int* c,
    double* x_i, double* y_i, double* theta_i,
    double* xt_i, double* yt_i, double* thetat_i,
    int NN, int NSS) {
    for (int i = 0; i < NSS; i++) {
        clt[i] = 0; // collision time count down
        cc[i] = 0; // collision counter
        c[i] = 0;
    }

    // Zero the initial values
    for (int j = 0; j < NN; j++) {
        x_i[j] = 0;   y_i[j] = 0;   theta_i[j] = 0;
        xt_i[j] = 0;  yt_i[j] = 0;  thetat_i[j] = 0;
    }
}

void init_memoryHeap(double* heap0, double* sign0, double* om_Old,
    int NSS, int LC, double Om) {
    int ccw, cw;
    int NSS2 = NSS / 2;
    for (int j = 0; j < NSS2; j++) {
        ccw = 0;
        cw = 0;
        for (int k = 0; k < LC; k++) {
            if ((k % 2) == 0) heap0[j * LC + k] = -1;
            else heap0[j * LC + k] = 1;
        }
        for (int k = 0; k < LC; k++) {
            if (heap0[j * LC + k] == 1) ccw++;
            else cw++;
        }
        if (ccw > cw) {
            sign0[j] = 1;
            om_Old[j] = Om;
        }
        else {
            sign0[j] = -1;
            om_Old[j] = -Om;
        }
    }

    for (int j = NSS2; j < NSS; j++) { //int j=NSS2; j<NSS;
        ccw = 0;
        cw = 0;
        for (int k = 0; k < LC; k++) {
            if ((k % 2) == 0) heap0[j * LC + k] = 1;
            else heap0[j * LC + k] = -1;
        }
        for (int k = 0; k < LC; k++) {
            if (heap0[j * LC + k] == 1) ccw++;
            else cw++;
        }
        if (ccw > cw) {
            sign0[j] = 1;
            om_Old[j] = Om;
        }
        else {
            sign0[j] = -1;
            om_Old[j] = -Om;
        }
    }
}
