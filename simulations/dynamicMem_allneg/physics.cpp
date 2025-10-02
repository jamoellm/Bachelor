#include "physics.h"

void spinnerInteraction_v10(double* x, double* y, double* theta, double* xt, double* yt, double* thetat, int NN, double tau0, double R0, double A1, double A2, double m, double I, double ks, double eta, double etaRot, double kax, double kay, std::vector<double> sign0, double* F, double* H, double* R) {
    for (int i = 0; i < 3 * NN; i++) F[i] = 0;
    double Fx[NN], Fy[NN], tau[NN], xtt[NN], ytt[NN], thetatt[NN];
    double d, d1, d2, d3, d4, vp, N, f;
    double kd = 2e5; double mu = 0.2;

    for (int j = 0; j < NN; j++) {
        Fx[j] = 0;  Fy[j] = 0;  tau[j] = 0;
        xtt[j] = 0; ytt[j] = 0; thetatt[j] = 0;
    }
    // driving torque
    for (int j = 0; j < NN; j++) {
        tau[j] = tau0 * sign0[j];
    }

    // cross-cross interaction between spinner j and k
    for (int j = 0; j < NN; j++) {
        for (int k = 0; k < j; k++) {
            d = norm(x[k] - x[j], y[k] - y[j]) - 2 * R0;
            if (d < 0) {
                crosscross_v8(x[j], y[j], theta[j], xt[j], yt[j], thetat[j], x[k], y[k], theta[k], xt[k], yt[k], thetat[k], H, R);
                Fx[j] += H[0];
                Fy[j] += H[1];
                tau[j] += H[2];
                Fx[k] += H[3];
                Fy[k] += H[4];
                tau[k] += H[5];
            }
        }
    }

    for (int j = 0; j < NN; j++) {
        // collision with the wall
        // right wall
        d1 = x[j] - (A1 - R0);
        if (xt[j] > 0) vp = xt[j];
        else vp = 0;
        if (d1 > 0) {
            N = -d1 * ks * (1 + (kd / ks) * vp);
            f = -mu * N;
            Fx[j] += N;
            if ((yt[j] + R0 * thetat[j]) > 0) Fy[j] -= f;
            else Fy[j] += f;
            tau[j] += R0 * Fy[j];
        }

        // left wall
        d2 = (-A1 + R0) - x[j];
        if (xt[j] < 0) vp = 0 - xt[j];
        else vp = 0;
        if (d2 > 0) {
            N = d2 * ks * (1 + (kd / ks) * vp);
            f = mu * N;
            Fx[j] += N;
            if ((yt[j] - R0 * thetat[j]) > 0) Fy[j] -= f;
            else Fy[j] += f;
            tau[j] -= R0 * Fy[j];
        }

        // top wall
        d3 = y[j] - (A2 - R0);
        if (yt[j] > 0) vp = yt[j];
        else vp = 0;
        if (d3 > 0) {
            N = -d3 * ks * (1 + (kd / ks) * vp);
            f = -mu * N;
            Fy[j] += N;
            if (xt[j] - R0 * thetat[j] > 0) Fx[j] -= f;
            else Fx[j] += f;
            tau[j] -= R0 * Fx[j];
        }

        // bottom wall
        d4 = (-A2 + R0) - y[j];
        if (yt[j] < 0) vp = -yt[j];
        else vp = 0;
        if (d4 > 0) {
            N = d4 * ks * (1 + (kd / ks) * vp);
            f = mu * N;
            Fy[j] += N;
            if (xt[j] + R0 * thetat[j] > 0) Fx[j] -= f;
            else Fx[j] += f;
            tau[j] += R0 * Fx[j];
        }

        // translational drag
        Fx[j] -= eta * xt[j];
        Fy[j] -= eta * yt[j];

        // rotational drag
        tau[j] -= etaRot * thetat[j];

        // inward air current force
        // Fx[j] -= m * kax * x[j] * x[j] * x[j];
        // Fy[j] -= m * kay * y[j] * y[j] * y[j];

        // from force/torque to acceleration/ang. acc.
        xtt[j] = Fx[j] / m;
        ytt[j] = Fy[j] / m;
        thetatt[j] = tau[j] / I;
    }

    for (int j = 0; j < NN; j++) {
        F[j] = xtt[j];
        F[j + NN] = ytt[j];
        F[j + 2 * NN] = thetatt[j];
    }
}

void crosscross_v8(double x1, double y1, double theta1, double xt1, double yt1, double thetat1, double x2, double y2, double theta2, double xt2, double yt2, double thetat2, double* H, double* R) {
    double Ri = 0.030;  double R0 = 0.038;
    double ks = 1e4; double kd = 2e5;
    const double pi = 3.14159265358979323846;
    double xs1[97], ys1[97], xs2[97], ys2[97];
    // H0&1: f1, H2: tau1, H3&4: f2. H5:tau2
    for (int j = 0; j < 6; j++) H[j] = 0;
    for (int j = 0; j < 7; j++) R[j] = 0;


    // Find the positions of the key points for both spinners
    double Dtheta = 2 * pi / 24;
    double t1 = 0;      double t2 = Dtheta / 2;
    double t, r, dd, p11x, p11y, p12x, p12y, p21x, p21y, p22x, p22y, pcx, pcy, d1x, d1y, d2x, d2y, co1x, co1y, co2x, co2y, vc1x, vc1y, vc2x, vc2y, vx, vy, vp, F0;
    for (int j = 0; j < 24; j++) {
        t = t1 + Dtheta * j;
        r = Ri;
        xs1[j * 2] = x1 + r * cos(t + theta1);
        ys1[j * 2] = y1 + r * sin(t + theta1);
        xs2[j * 2] = x2 + r * cos(t + theta2);
        ys2[j * 2] = y2 + r * sin(t + theta2);

        t = t2 + Dtheta * j;
        r = R0;
        xs1[1 + j * 2] = x1 + r * cos(t + theta1);
        ys1[1 + j * 2] = y1 + r * sin(t + theta1);
        xs2[1 + j * 2] = x2 + r * cos(t + theta2);
        ys2[1 + j * 2] = y2 + r * sin(t + theta2);

    }
    xs1[48] = xs1[0]; ys1[48] = ys1[0];
    xs2[48] = xs2[0]; ys2[48] = ys2[0];

    for (int j = 0; j < 48; j++) {
        for (int k = 0; k < 48; k++) {
            p11x = xs1[j];    p11y = ys1[j];
            p12x = xs1[j + 1];  p12y = ys1[j + 1];
            p21x = xs2[k];    p21y = ys2[k];
            p22x = xs2[k + 1];  p22y = ys2[k + 1];
            dd = norm(0.5 * (p11x + p12x) - 0.5 * (p21x + p22x), 0.5 * (p11y + p12y) - 0.5 * (p21y + p22y));
            if (dd < 0.05) {
                lineline(p11x, p11y, p12x, p12y, p21x, p21y, p22x, p22y, R);
                d1x = R[0]; d1y = R[1];
                d2x = R[2]; d2y = R[3];
                pcx = R[4]; pcy = R[5];

                if (R[6] == 1) {
                    co1x = pcx - x1; co1y = pcy - y1;
                    co2x = pcx - x2; co2y = pcy - y2;
                    vc1x = -thetat1 * co1y + xt1; vc1y = thetat1 * co1x + yt1;
                    vc2x = -thetat2 * co2y + xt2; vc2y = thetat2 * co2x + yt2;
                    vx = vc1x - vc2x;
                    vy = vc1y - vc2y;
                    vp = -(vx * d1x + vy * d1y) / norm(d1x, d1y);
                    if (vp < 0) vp = 0;
                    F0 = ks * (1 + (kd / ks) * vp);
                    H[0] += F0 * d1x; // f1x
                    H[1] += F0 * d1y; // f1y
                    H[3] += F0 * d2x; // f2x
                    H[4] += F0 * d2y; // f2y

                    H[2] += F0 * cross(co1x, co1y, d1x, d1y); // tau1
                    H[5] += F0 * cross(co2x, co2y, d2x, d2y); // tau2
                }
            }
        }
    }
    // return H;
}

void lineline(double p11x, double p11y, double p12x, double p12y, double p21x, double p21y, double p22x, double p22y, double* R) {
    for (int j = 0; j < 7; j++) R[j] = 0;

    // R0&1:d1, R2&3:d2, R4&5:pc, R6:collision or not
    double a = p12y; double b = p12x; double c = p11y; double d = p11x;
    double e = p22y; double f = p22x; double g = p21y; double h = p21x;
    double pcx = (d * (-f * g + a * (f - h) + e * h) + b * (f * g - e * h + c * (-f + h))) / (d * (e - g) + b * (-e + g) + (a - c) * (f - h));
    double pcy = (a * (d * e - d * g - e * h + f * g) + b * c * (g - e) + c * e * h - c * f * g) / ((a - c) * (f - h) + b * (g - e) + d * (e - g));
    double d1 = (pcx - p11x) * (pcx - p12x) + (pcy - p11y) * (pcy - p12y);
    double d2 = (pcx - p21x) * (pcx - p22x) + (pcy - p21y) * (pcy - p22y);
    double rx, ry, ex, ey, d0, d1x, d1y, d2x, d2y;
    d1x = 0; d1y = 0; d2x = 0; d2y = 0;
    if (!((d1 > 0) || (d2 > 0))) {
        R[6] = 1;
        double dd1 = norm(pcx - p11x, pcy - p11y);
        double dd2 = norm(pcx - p12x, pcy - p12y);
        double dd3 = norm(pcx - p21x, pcy - p21y);
        double dd4 = norm(pcx - p22x, pcy - p22y);
        double M = std::min({ dd1,dd2,dd3,dd4 });
        if (dd1 == M) {
            rx = p11x - pcx;
            ry = p11y - pcy;
            ex = p22x - p21x;
            ey = p22y - p21y;
            d0 = rx * ex + ry * ey;
            d2x = rx - d0 * ex / pow(norm(ex, ey), 2);
            d2y = ry - d0 * ey / pow(norm(ex, ey), 2);
            d1x = -d2x;
            d1y = -d2y;
        }
        if (dd2 == M) {
            rx = p12x - pcx;
            ry = p12y - pcy;
            ex = p22x - p21x;
            ey = p22y - p21y;
            d0 = rx * ex + ry * ey;
            d2x = rx - d0 * ex / pow(norm(ex, ey), 2);
            d2y = ry - d0 * ey / pow(norm(ex, ey), 2);
            d1x = -d2x;
            d1y = -d2y;
        }
        if (dd3 == M) {
            rx = p21x - pcx;
            ry = p21y - pcy;
            ex = p12x - p11x;
            ey = p12y - p11y;
            d0 = rx * ex + ry * ey;
            d1x = rx - d0 * ex / pow(norm(ex, ey), 2);
            d1y = ry - d0 * ey / pow(norm(ex, ey), 2);
            d2x = -d1x;
            d2y = -d1y;
        }
        if (dd4 == M) {
            rx = p22x - pcx;
            ry = p22y - pcy;
            ex = p12x - p11x;
            ey = p12y - p11y;
            d0 = rx * ex + ry * ey;
            d1x = rx - d0 * ex / pow(norm(ex, ey), 2);
            d1y = ry - d0 * ey / pow(norm(ex, ey), 2);
            d2x = -d1x;
            d2y = -d1y;
        }
    }
    R[0] = d1x; R[1] = d1y; R[2] = d2x; R[3] = d2y;
    R[4] = pcx; R[5] = pcy;
};

double norm(double x, double y)
{
    return sqrt(x * x + y * y);
}

double cross(double x1, double y1, double x2, double y2)
{
    return x1 * y2 - x2 * y1;
}

double weighting(double memoryCell, int position, int personality) {
    if (personality == 0) { // pushover
        return memoryCell;
    }
    else if (personality == 1) { // opportunist
        return memoryCell * pow(2, position);
    }
    else if (personality == 2) { // traditionalist
        return memoryCell * pow(2, -position);
    }
    else if (personality == 3) { // contrarian
        return -memoryCell;
    }
    else if (personality == 4) { // curmudgeon
        return 1;
    }
    else if (personality == -4) { // curmudgeon
        return -1;
    }
    return 0.0;
};
