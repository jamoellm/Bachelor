//  main.cpp
//  spinner
//
//  Created by Shengkai Li on 1/29/23.
//
#include <filesystem>
#include <string>
#include "omp.h"

#include "write2buffer.h"
#include "init.h"
#include "physics.h"

const double T = 1800;
const int numberUpscaler = 2;
constexpr int NN = 8 * numberUpscaler * numberUpscaler;

const int TOTAL_SIMSETS = 1;
int threadUnsaveCounter = 0;

std::string folder;

int main(int argc, const char* argv[]) {
    double tL = 0.010; // loop time in seconds
    folder = "./outData/staticMem32/";
    if (!std::filesystem::exists(folder)) {
        std::filesystem::create_directories(folder);
        std::cout << "Ordner erstellt: " + folder + "\n";
    }

    int LA_array[5] = { 5,7,9,15,21 };
    double density_factors[7] = { 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0 }; // {0.5, 1.5, 2.5};

    omp_set_num_threads(1);
#pragma omp parallel for collapse(3)
    for (int i = 0; i < 1; i++) {// run different seeds...
        for (int j = 3; j < 4; j++) {// for different M sizes...
            for (int k = 2; k < 3; k++) {// for different lambdas...
                int tid = omp_get_thread_num();
                std::string s1 = "Thread " + std::to_string(tid) + " edits (" + std::to_string(j) + ", " + std::to_string(i) + ") - ";
                std::string s2 = "M=" + std::to_string(LA_array[j]) + " @seed:" + std::to_string(i) + " - density:" + std::to_string(density_factors[k]) + "\n";
                std::cout << s1 + s2;

                std::string filename = "spinner_T" + std::to_string((int)T)
                    + "_LC" + std::to_string(LA_array[j])
                    + "_density" + std::to_string(density_factors[k])
                    + "_s" + std::to_string((int)i);

                evolution(T, tL, LA_array[j], i, numberUpscaler, density_factors[k], filename);
            }
        }
    }
    return 0;
}

/**
* @param T:  Simulation Time in Seconds
* @param tL: loop time in seconds
* @param LC: Size of Memory -> 'M' in paper
* @param numberUpscaler scales the number of spinners while keeping the density constant
+ @param lambda free parameter
* @param s0: random seed
*/
void evolution(double T, double tL, int LC, double s0, int numberUpscaler, double lambda, std::string filename) {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::microseconds;

    auto t1 = high_resolution_clock::now();
    const double pi = 3.14159265358979323846;

    uint64_t fileIdx = 0;
    int timePerFile = 600; // seconds of data per file
    int numberOfFiles = ceil(T / timePerFile);
    srand(s0);

    std::ofstream myfile(folder + filename + "_" + std::to_string(fileIdx) + ".dat");

    double dt = 1e-4;
    int stepsPerFile = ceil(timePerFile / dt);
    int NT = ceil(T / dt); // Number of Time Steps
    int NL = ceil(tL / dt); // Number of Loops before ardiuino does stuff
    int NS = 100; // Number of Steps between Logs, dt_log = dt*NS
    int NLogs = 10; // how many Updates are printed into Console
    int L0 = ceil(NT / NLogs); // Steps to take before a Log into Console
    int percentageTracker = 0;
    int LA = 10; // acc detec window
    int LD = 10; //detection interval
    int LDM = 100; // min collision detection
    double std0_threshold = 0.02 * 9.8; // the unit is g in Arduino, so need to convert here to SI
    double om_threshold = 30 * pi / 180; // the unit is deg/s in Arduino, so need to convert here to SI
    //int NM0=round(100/dt); //round(150/dt)
    double R0 = 0.038; double Ri = 0.030;
    double m = 0.0525; double I = 0.33 * m * R0 * R0;
    double r_imu = 0.002;
    double gamma0 = 1.5; double Om = 28;
    double tau0 = I * gamma0;
    double A1 = 0.38 * numberUpscaler / sqrt(lambda); double A2 = 0.28 * numberUpscaler / sqrt(lambda); // width and height of arena. kw: walls
    double ks = 1e4; double kd = 2e5; // there is a repeated definition afterwards in functions
    double eta = (6.5 * 0.01) * m; double etaRot = I * gamma0 / Om;
    double kax = 2; double kay = 2;
    double v0 = 0.05;
    double xp, yp, ff, d, tt;
    double x_i[NN], y_i[NN], theta_i[NN];
    double x_f[NN], y_f[NN], theta_f[NN];
    double xt_i[NN], yt_i[NN], thetat_i[NN];
    double xt_f[NN], yt_f[NN], thetat_f[NN];
    double xt1[NN], yt1[NN], thetat1[NN];
    double sign0[NN]; // list of phenotypical signs
    int clt[NN]; // collision time count down
    int cc[NN]; // collision counter
    int c[NN];
    double om_Old[NN];

    // personality types encoded as integers
    // 0: pushover
    // 1: opportunist
    // 2: traditionalist
    // 3: contrarian
    // +4: curmudgeon +
    // -4: curmudgeon -
    int personalities[NN];

    double* heap0 = (double*)malloc(NN * LC * sizeof(double)); // list of all remembered signs
    double* a0 = (double*)malloc(NN * LA * sizeof(double));

    double xtt1, ytt1, thetatt1, xtt2, ytt2, thetatt2;
    double avg0, std0, r_om, om_imu;
    int id, idc, s_other;

    init_fillZeros(clt, cc, c,
        x_i, y_i, theta_i,
        xt_i, yt_i, thetat_i,
        NN, NN);

    // Assign the initial positions randomly
    int oppositeSpinners = 0;
    init_placement(oppositeSpinners, x_i, y_i, NN, A1, A2, R0);
    init_velocity_Angle(xt_i, yt_i, theta_i, NN, v0, pi);

    // Assign the initial heap and torque for smart spinners
    init_memoryHeap(heap0, sign0, om_Old, NN, LC, Om);

    // write personalities of the smart spinners
    for (int j = 0; j < NN; j++) {
        personalities[j] = 0; // only pushovers
        // dynamicInfluenceValues[j] = 1; // + delta_dynamicInfluenceValues[j];
    }

    // write meta data / base parameters into first line of file
    // int NN, int NN, double tL, int LC, double s0
    // ':' to seperate 'meta', ',' to seperate quantities, ' ' to seperate keys from values
    std::string meta = "meta:T " + std::to_string(T)
        + ",dt_log " + std::to_string(NS * dt) + ",n_files " + std::to_string(numberOfFiles)
        + ",r_out " + std::to_string(R0) + ",r_inner " + std::to_string(Ri)
        + ",A1 " + std::to_string(A1) + ",A2 " + std::to_string(A2)
        + ",NN " + std::to_string(NN)
        + ",tL " + std::to_string(tL) + ",LC " + std::to_string(LC);
        // + ",curmudgeons " + std::to_string(curmudgeons);
    myfile << meta << "\n";

    std::vector<double> ss0(NN);
    double F[3 * NN], H[6], R[7];
    // main loop
    for (int i = 0; i < NT; i++) {
        ss0.assign(sign0, sign0 + NN);
        // v(t+0.5dt)=v(t)+0.5*a(t)*dt
        spinnerInteraction_v10(x_i, y_i, theta_i, xt_i, yt_i, thetat_i, NN, tau0, R0, A1, A2, m, I, ks, eta, etaRot, kax, kay, ss0, F, H, R);
        for (int j = 0;j < NN;j++) {
            xtt1 = F[j];
            ytt1 = F[j + NN];
            thetatt1 = F[j + 2 * NN];
            xt1[j] = xt_i[j] + 0.5 * dt * xtt1;
            yt1[j] = yt_i[j] + 0.5 * dt * ytt1;
            thetat1[j] = thetat_i[j] + 0.5 * dt * thetatt1;
        }

        // x(t+dt)=x(t)+v(t+0.5dt)*dt
        for (int j = 0;j < NN;j++) {
            x_f[j] = x_i[j] + dt * xt1[j];
            y_f[j] = y_i[j] + dt * yt1[j];
            theta_f[j] = theta_i[j] + dt * thetat1[j];
        }

        // a(t+dt)=a(x(t+dt))
        // v(t+dt)=v(t+0.5dt)+0.5*a(t+dt)*dt
        spinnerInteraction_v10(x_f, y_f, theta_f, xt1, yt1, thetat1, NN, tau0, R0, A1, A2, m, I, ks, eta, etaRot, kax, kay, ss0, F, H, R);
        for (int j = 0; j < NN; j++) {
            xtt2 = F[j];
            ytt2 = F[j + NN];
            thetatt2 = F[j + 2 * NN];
            xt_f[j] = xt1[j] + 0.5 * dt * xtt2;
            yt_f[j] = yt1[j] + 0.5 * dt * ytt2;
            thetat_f[j] = thetat1[j] + 0.5 * dt * thetatt2;
        }

        // today becomes yesterday
        for (int j = 0;j < NN;j++) {
            x_i[j] = x_f[j];
            y_i[j] = y_f[j];
            theta_i[j] = theta_f[j];
            xt_i[j] = xt_f[j];
            yt_i[j] = yt_f[j];
            thetat_i[j] = thetat_f[j];
        }

        // loop on Arduino if the time interval reaches multiples of 0.01 s (or other number)
        for (int j = 0; j < NN; j++) {
            if (i % NL == 0) {
                if (clt[j] < 1) (clt[j]) = 0;
                else (clt[j]) -= 1;

                id = c[j] % LA;
                c[j] += 1;

                // save acc into acc heap
                xtt2 = F[j];
                ytt2 = F[j + NN];
                om_imu = thetat_f[j];
                a0[j * LA + id] = sqrt(xtt2 * xtt2 + ytt2 * ytt2 + om_imu * om_imu * r_imu);

                if ((c[j]) % LD == 0) {
                    if (clt[j] == 0) {
                        // evaluate acc fluctuation
                        avg0 = 0;
                        for (int k = 0; k < LA; k++) avg0 += a0[j * LA + k];
                        avg0 = avg0 / LA;
                        std0 = 0;
                        for (int k = 0; k < LA; k++) std0 += (a0[j * LA + k] - avg0) * (a0[j * LA + k] - avg0);
                        std0 = sqrt(std0 / LA);
                        if (std0 > std0_threshold) { // if in consistent state
                            r_om = ((om_Old[j]) - (thetat_f[j])) / (om_Old[j]); //evaluate spin drop
                            if (abs(om_Old[j]) > om_threshold) { //if the spin is significant
                                // push heap
                                cc[j] += 1;
                                idc = (cc[j]) % LC;

                                double w = 1;
                                if (om_Old[j] > 0) {
                                    if (r_om < 0.3) heap0[j * LC + idc] = -1 * w; // the one I met is opposite
                                    else heap0[j * LC + idc] = 1; // the one I met is same
                                }
                                if (om_Old[j] < 0) {
                                    if (r_om < 0.3) heap0[j * LC + idc] = 1 * w; // the one I met is opposite
                                    else heap0[j * LC + idc] = -1; // the one I met is same
                                }

                                // evaluate heap
                                double res = 0;
                                for (int k = 0; k < LC; k++) { // iterates over the memory heap,
                                    int K = (idc - k + LC) % LC; // from most to least recently changed.
                                    res += weighting(heap0[j * LC + K], k, personalities[j]);
                                }
                                if (res > 0) s_other = 1;
                                else s_other = -1;

                                // change state
                                if ((s_other * (sign0[j])) < 0) sign0[j] = -sign0[j];
                                clt[j] = LDM;
                            }
                        }
                    }
                    om_Old[j] = thetat_f[j];
                }
            }
        }

        if (i % stepsPerFile == 0 && i != 0 && i != NT) {
            std::cout << fileIdx << "closed\n";
            myfile.close();
            fileIdx += 1;
            myfile.open(folder + filename + "_" + std::to_string(fileIdx) + ".dat");
        }

        if (i % NS == 0) {
            std::string buffer;
            buffer.reserve(6 * NN + NN * (3 + LC));

            writeArray(theta_f, NN, buffer, " ");
            writeArray(x_f, NN, buffer, " ");
            writeArray(y_f, NN, buffer, " ");

            // writeArray(thetat_f, NN, buffer, " ");
            // writeArray(xt_f, NN, buffer, " ");
            // writeArray(yt_f, NN, buffer, " ");

            // write2DArray(a0, NN, LA, buffer, " ");
            buffer += "\n";
            writeArray(sign0, NN, buffer, " ");
            writeArray(cc, NN, buffer, " ");
            // write2DArray(heap0, NN, LC, buffer, " ");
            buffer += "\n\n";

            myfile << buffer;
        }

        if (i % L0 == 0) {
            std::string s = "";
            s += "[";
            for (int j = 0; j < 10; j++) {
                if (j < percentageTracker) s += "#";
                else if (j == percentageTracker)s += ">";
                else s += " ";
            };
            s += "]";
            s += " - (" + std::to_string(omp_get_thread_num()) + ") M=" + std::to_string(LC) + " @seed:" + std::to_string((int)s0) + "\n";
            if (percentageTracker == 10) {
                s += "Done (" + std::to_string(omp_get_thread_num()) + ")\n";
            }
            percentageTracker += 1;
            std::cout << s;
            std::cout.flush();
        }
    }

    int threadSaveCounter = ++threadUnsaveCounter;
    myfile.close();
    auto t2 = high_resolution_clock::now();
    auto ms_int = duration_cast<std::chrono::seconds>(t2 - t1);
    std::string finalOut = std::to_string(ms_int.count()) + " s" + " - " + std::to_string(threadSaveCounter) + "/" + std::to_string(TOTAL_SIMSETS) + " files";
    std::cout << finalOut + "\n";
}

