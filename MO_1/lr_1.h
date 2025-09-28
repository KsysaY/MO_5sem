#pragma once
#include "defines.h"
// using function_id = std::function<double(double)>;

// template<typename... params_t>
// const char* fstring(const char* format, const params_t&... params) {
//     static constexpr size_t buffer_size = 1024;
//     static char buffer[buffer_size];
//     sprintf_s(buffer, format, params...);
//     return buffer;
// }

enum class methods_types : uint8_t {
    bisect,
    golden_ratio,
    fibonacci
};

struct search_result {
    friend std::ostream& operator<<(std::ostream& stream, const search_result& res);

    const methods_types method;
    const size_t iterations;
    const size_t function_calls;
    const double accuracy;
    const double result;

    search_result(methods_types m, size_t iters, size_t calls, double acc, double res)
        :method        (    m)
        ,iterations    (iters)
        ,function_calls(calls)
        ,accuracy      (  acc)
        ,result        (  res)
    {}
};

// Прототипы функций
search_result bisect      (function_1d func, double lhs, double rhs, double eps = ACCURACY_MIDDLE, size_t max_iterations = MAX_ITERATIONS);
search_result golden_ratio(function_1d func, double lhs, double rhs, double eps = ACCURACY_MIDDLE, size_t max_iterations = MAX_ITERATIONS);
search_result fibonacci   (function_1d func, double lhs, double rhs, double eps = ACCURACY_MIDDLE);

/// std::ostream& operator<<(std::ostream& stream, const search_result& result);
