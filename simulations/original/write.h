// write.h
#pragma once
#include <vector>
#include <array>
#include <fstream>
#include <string>

inline void writeArray(const double* array, int size, std::ofstream& myfile) {
    for (int j=0;j<size;j++) {
        myfile << *(array+j) << " ";
    }
    myfile << "\n";
}

inline void writeArray(const int* array, int size, std::ofstream& myfile) {
    for (int j=0;j<size;j++) {
        myfile << *(array+j) << " ";
    }
    myfile << "\n";
}

inline void write2DArray(const double* array, int size1, int size2, std::ofstream& myfile) {
    for (int j=0;j<size1;j++) {
        for (int k=0;k<size2;k++) {
            myfile << *(array + j * size2 + k) << " ";
        }
        myfile << "\n";
    }
}
