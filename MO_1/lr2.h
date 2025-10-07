#pragma once
#include "defines.h"
#include "numerics/numerics.h"

enum class methods_types : uint8_t {
    bisect,
    golden_ratio,
    fibonacci
};

struct search_result_n {
    friend std::ostream& operator<<(std::ostream& stream, const search_result_n& res);

    const methods_types method;
    const size_t iterations;
    const size_t function_calls;
    const double accuracy;
    const double result;

    search_result_n(methods_types m, size_t iters, size_t calls, double acc, double res)
        :method(m)
        , iterations(iters)
        , function_calls(calls)
        , accuracy(acc)
        , result(res)
    {
    }
};