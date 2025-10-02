// write.h
#pragma once
#include <array>
#include <fstream>
#include <string>
#include <vector>

inline void writeArray(const double* array, int size, std::string& buffer, const std::string& delimiter) {
    for (int i = 0; i < size; ++i) {
        buffer += std::to_string(array[i]) + delimiter;
    }
    buffer += "\n";
}

inline void writeArray(const int* array, int size, std::string& buffer, const std::string& delimiter) {
    for (int i = 0; i < size; ++i) {
        buffer += std::to_string(array[i]) + delimiter;
    }
    buffer += "\n";
}

inline void write2DArray(const double* array, int rows, int cols, std::string& buffer, const std::string& delimiter) {
    for (int j = 0; j < rows; ++j) {
        for (int k = 0; k < cols; ++k) {
            buffer += std::to_string(array[j * cols + k]) + delimiter;
        }
        buffer += "\n";
    }
}