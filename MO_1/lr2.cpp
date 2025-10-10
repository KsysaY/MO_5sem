#include "lr2.h"
#include "lr_1.h"

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
    case methods_types_n::per_coordinate_descend:
        stream << fstring("method: per_coordinate_descend\n");
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

search_result_n golden_ratio(function_nd function, const numerics::vector_f64& left, const numerics::vector_f64& right, const F64 eps, const I32 max_iterations) 
{
    numerics::vector_f64 lhs(left), rhs(right);
    numerics::vector_f64 x_l = rhs - (rhs - lhs) * PSI;
    numerics::vector_f64 x_r = lhs + (rhs - lhs) * PSI;
    F64 f_l = function(x_l);
    F64 f_r = function(x_r);
    UI64 iterations = 0;
    UI64 function_iter = 2;
    for (; iterations != max_iterations && (numerics::vector_f64::distance(lhs, rhs)) > 2 * eps; iterations++, function_iter++)
    {
        if (f_l > f_r)
        {
            lhs = x_l;
            x_l = x_r;
            f_l = f_r;
            x_r = lhs + (rhs - lhs) * PSI;
            f_r = function(x_r);
        }
        else
        {
            rhs = x_r;
            x_r = x_l;
            f_r = f_l;
            x_l = rhs - (rhs - lhs) * PSI;
            f_l = function(x_l);
        }
    }

    numerics::vector_f64 final_result = (lhs + rhs) * 0.5;
    F64 final_accuracy = numerics::vector_f64::distance(lhs, rhs);
    return search_result_n(methods_types_n::golden_ratio, iterations, function_iter, final_accuracy, final_result);
}

search_result_n fibonacci(function_nd function, const numerics::vector_f64& left, const numerics::vector_f64& right, const F64 eps)
{
    numerics::vector_f64 lhs(left), rhs(right);
    F64 condition = numerics::vector_f64::distance(lhs, rhs) / eps;
    F64 fib_t{ 0.0 }, fib_1{ 1.0 }, fib_2{ 1.0 };
    UI64 iterations = 0;
    UI64 function_iter = 2;
    while (fib_2 < condition)
    {
        fib_t = fib_1;
        fib_1 = fib_2;
        fib_2 += fib_t;
        iterations++;
        function_iter++;
    }
    numerics::vector_f64 x_l = lhs + (rhs - lhs) * ((fib_2 - fib_1) / fib_2);
    numerics::vector_f64 x_r = lhs + (rhs - lhs) * (fib_1 / fib_2);
    F64 f_l = function(x_l);
    F64 f_r = function(x_r);
    for (I32 index = iterations; index > 0; index--)
    {
        fib_t = fib_2 - fib_1;
        fib_2 = fib_1;
        fib_1 = fib_t;
        if (f_l > f_r)
        {
            lhs = x_l;
            f_l = f_r;
            x_l = x_r;
            x_r = lhs + (rhs - lhs) * (fib_1 / fib_2);
            f_r = function(x_r);
        }
        else
        {
            rhs = x_r;
            x_r = x_l;
            f_r = f_l;
            x_l = lhs + (rhs - lhs) * ((fib_2 - fib_1) / fib_2);
            f_l = function(x_l);
        }
    }

    numerics::vector_f64 final_result = (lhs + rhs) * 0.5;
    F64 final_accuracy = numerics::vector_f64::distance(lhs, rhs);
    return search_result_n(methods_types_n::fibonacci, iterations, function_iter, final_accuracy, final_result);
}

search_result_n per_coord_descend(function_nd function, const numerics::vector_f64& x_start, const F64 eps, const I32 max_iters)
{
    UI64 total_probes = 0;
    numerics::vector_f64 x_0(x_start);
    numerics::vector_f64 x_1(x_start);
    numerics::vector_f64 result;
    F64 step = 1.0;
    F64 x_i, y_1, y_0;
    I32 opt_coord_n = 0, coord_id;
    UI64 iterations;
    for (iterations = 0; iterations < max_iters; ++iterations)
    {
        coord_id = iterations % x_0.size();
        x_1[coord_id] -= eps;
        y_0 = function(x_1);
        x_1[coord_id] += 2.0 * eps;
        y_1 = function(x_1);
        x_1[coord_id] -= eps;
        x_1[coord_id] = y_0 > y_1 ? x_1[coord_id] += step : x_1[coord_id] -= step;
        x_i = x_0[coord_id];
        search_result_n result = fibonacci(function, x_0, x_1, eps);
        x_0 = result.result;
        total_probes += result.function_calls + 2;
        if (std::abs(x_1[coord_id] - x_i) < 2 * eps)
        {
            opt_coord_n++;
            if (opt_coord_n == x_1.size())
            {
                break;
            }
            continue;
        }
        opt_coord_n = 0;
    }

    return search_result_n(methods_types_n::per_coordinate_descend, iterations, total_probes, , x_0);
}