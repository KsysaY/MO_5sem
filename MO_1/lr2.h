#pragma once
#include "defines.h"
#include "numerics/numerics.h"

enum class methods_types_n : uint8_t {
    bisect,
    golden_ratio,
    fibonacci
};

struct search_result_n {
    friend std::ostream& operator<<(std::ostream& stream, const search_result_n& res);

    const methods_types_n method;
    const UI64 iterations;
    const UI64 function_calls;
    const F64 accuracy;
    const numerics::vector_f64 result;

    search_result_n(methods_types_n m, UI64 iters, UI64 calls, F64 acc, numerics::vector_f64& res)
        :method(m)
        , iterations(iters)
        , function_calls(calls)
        , accuracy(acc)
        , result(res)
    {}
};

search_result_n bisect(function_nd func, const numerics::vector_f64& lhs, const numerics::vector_f64& rhs, const F64 eps = ACCURACY_MIDDLE, const I32 max_iterations = MAX_ITERATIONS);