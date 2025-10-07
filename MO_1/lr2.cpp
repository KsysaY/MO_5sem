#include "lr2.h"

std::ostream& operator<<(std::ostream& stream, const search_result_n& result) {
    switch (result.method) {
    case methods_types::bisect:
        stream << fstring("method: bisect\n");
        break;
    case methods_types::golden_ratio:
        stream << fstring("method: golden_ratio\n");
        break;
    case methods_types::fibonacci:
        stream << fstring("method: fibonacci\n");
        break;
    default:
        stream << fstring("method: Unknown\n");
    }
    stream << fstring("iterations:     %ld\n", result.iterations);
    stream << fstring("function_calls: %ld\n", result.function_calls);
    stream << fstring("accuracy:       %.10f\n", result.accuracy);
    stream << fstring("result:         %.10f\n", result.result);
    return stream;
}

//static void bisect(function_nd function, const numerics::vector_f64& left, const numerics::vector_f64& right, const F64 eps, const I32 max_iterations)
//{
//    result.clear();
//    result.type = search_method_type_nd::bisection;
//    numerics::vector_f64 dir, lhs(left), rhs(right);
//    dir = numerics::vector_f64::direction(lhs, rhs) * eps;
//    for (; result.iterations != max_iterations && (result.accuracy = numerics::vector_f64::distance(lhs, rhs)) > 2 * eps; result.iterations++, result.function_probes += 2)
//    {
//        result.result = (lhs + rhs) * 0.5;
//        if (function(result.result - dir) > function(result.result + dir))
//            lhs = result.result;
//        else
//            rhs = result.result;
//    }
//}