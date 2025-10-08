#include "lr2.h"

std::ostream& operator<<(std::ostream& stream, const search_result_n& result) {
    switch (result.method) {
    case methods_types_n::bisect:
        stream << fstring("method: bisect\n");
        break;
    case methods_types_n::golden_ratio:
        stream << fstring("method: golden_ratio\n");
        break;
    case methods_types_n::fibonacci:
        stream << fstring("method: fibonacci\n");
        break;
    default:
        stream << fstring("method: Unknown\n");
    }
    stream << fstring("iterations:     %ld\n", result.iterations);
    stream << fstring("function_calls: %ld\n", result.function_calls);
    stream << fstring("accuracy:       %.10f\n", result.accuracy);

    stream << "result: ";
    stream << "[";
    for (size_t i = 0; i < result.result.size(); ++i) {
        stream << std::setprecision(10) << result.result[i];
        if (i < result.result.size() - 1) {
            stream << ", ";
        }
    }
    stream << "]\n";

    return stream;
}

search_result_n bisect(function_nd function, const numerics::vector_f64& left, const numerics::vector_f64& right, const F64 eps, const I32 max_iterations)
{
    numerics::vector_f64 dir, lhs(left), rhs(right);
    dir = numerics::vector_f64::direction(lhs, rhs) * eps;
    UI64 iterations = 0;
    UI64 function_iter = 0;
    numerics::vector_f64 result;

    for (; iterations != max_iterations && (numerics::vector_f64::distance(lhs, rhs)) > 2 * eps; iterations++)
    {
        result = (lhs + rhs) * 0.5;
        function_iter += 2;
        if (function(result - dir) > function(result + dir))
            lhs = result;
        else
            rhs = result;
    }

    numerics::vector_f64 final_result = (lhs + rhs) * 0.5;
    F64 final_accuracy = numerics::vector_f64::distance(lhs, rhs);
    return search_result_n(methods_types_n::bisect, iterations, function_iter, final_accuracy, final_result);
}