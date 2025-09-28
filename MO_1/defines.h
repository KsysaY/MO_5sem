#pragma once
#ifndef __DEFINES_H__
#define __DEFINES_H__
#include <functional>
#include <iostream>
#include <cstdint>

constexpr double PSI = 0.61803398874989484820;
constexpr double PHI = 1.61803398874989484820;
constexpr size_t MAX_ITERATIONS = 1000;
constexpr double ACCURACY_LOW    = 1e-3;
constexpr double ACCURACY_MIDDLE = 1e-6;
constexpr double ACCURACY_HEIGHT = 1e-9;

using function_1d = std::function<double(double)>;

template<typename... params_t>
const char* fstring(const char* format, const params_t&... params) {
    static constexpr size_t buffer_size = 1024;
    static char buffer[buffer_size];
    sprintf_s(buffer, format, params...);
    return buffer;
}
#endif //  __DEFINES_H__
