#include "lr2.h"
#include "lr_1.h"
#include "numerics/linalg/numeric_vector.h"

search_result_n::search_result_n(methods_types_n m, UI64 iters, UI64 calls, F64 acc, numerics::vector_f64&& res)
    :method(m)
    , iterations(iters)
    , function_calls(calls)
    , accuracy(acc)
    , result(std::move(res))
{
}


std::ostream& operator<<(std::ostream& stream, const search_result_n& result) {
    switch (result.method) {
    case methods_types_n::bisect:
        stream << fstring(                "method: bisect\n");
        break;
    case methods_types_n::golden_ratio:
        stream << fstring(          "method: golden_ratio\n");
        break;
    case methods_types_n::fibonacci:
        stream << fstring(             "method: fibonacci\n");
        break;
    case methods_types_n::per_coordinate_descend:
        stream << fstring("method: per_coordinate_descend\n");
        break;
    case methods_types_n::gradient_descend:
        stream << fstring(      "method: gradient_descend\n");
        break;
    case methods_types_n::conj_gradient_descend:
        stream << fstring( "method: conj_gradient_descend\n");
        break;
    case methods_types_n::newton_raphson:
        stream << fstring(        "method: newton_raphson\n");
        break;
    case methods_types_n::internal_penalty:
        stream << fstring(      "method: internal_penalty\n");
        break;
    case methods_types_n::external_penalty:
        stream << fstring(      "method: external_penalty\n");
        break;
    default:
        stream << fstring("method: Unknown\n");
    }
    stream << fstring("iterations:     %ld\n", result.iterations);
    stream << fstring("function_calls: %ld\n", result.function_calls);
    stream << fstring("accuracy:       %10.4e\n", result.accuracy);

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
    numerics::vector_f64 dir(numerics::vector_f64::direction(left, right) *= eps), lhs(left), rhs(right);
    
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

    // numerics::vector_f64 final_result = (lhs + rhs) * 0.5;
    return search_result_n(
        methods_types_n::bisect,                 //
        iterations,                              // 
        function_iter,                           //
        numerics::vector_f64::distance(lhs, rhs),//
        std::move(result));                      //
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

    return search_result_n(
        methods_types_n::golden_ratio,           //
        iterations,                              //
        function_iter,                           //
        numerics::vector_f64::distance(lhs, rhs),//
        (lhs + rhs) * 0.5);                      //

}

search_result_n fibonacci(function_nd function, const numerics::vector_f64& left, const numerics::vector_f64& right, const F64 eps)
{
    numerics::vector_f64 lhs(left), rhs(right);
    F64 condition = numerics::vector_f64::distance(lhs, rhs) / eps;
    F64 fib_t{ 0.0 }, fib_1{ 1.0 }, fib_2{ 1.0 };
    UI64 iterations = 0;
    while (fib_2 < condition)
    {
        fib_t = fib_1;
        fib_1 = fib_2;
        fib_2 += fib_t;
        iterations++;
    }
    const UI64 function_iter = iterations+2;
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

    return search_result_n(
        methods_types_n::fibonacci,              // Вызываемый метод
        iterations,                              // Количество итераций
        function_iter,                           // Количество вызовов функции
        numerics::vector_f64::distance(lhs, rhs),// Точность
        (lhs + rhs) * 0.5);                      // Экстремум
}

search_result_n per_coord_descend(function_nd function, const numerics::vector_f64& x_start, const F64 eps, const I32 max_iters)
{
    UI64 total_probes = 0;
    F64 step = 1.0;
    I32 opt_coord_n = 0;
    UI64 iterations;
    numerics::vector_f64 x_curr(x_start), x_next(x_start);
    F64 accuracy = std::numeric_limits<double>::infinity();

    for (iterations = 0; iterations < max_iters; ++iterations)
    {
        I32 coord_id = iterations % x_curr.size();
        F64 original_value = x_curr[coord_id];

        x_curr[coord_id] -= eps;
        const F64 f_left = function(x_curr);
        x_curr[coord_id] += 2 * eps;
        const F64 f_right = function(x_curr);
        x_curr[coord_id] = original_value;
        total_probes += 2;

        const F64 search_direction(f_left > f_right ? step : -step);
        x_next[coord_id] = x_curr[coord_id] + search_direction;

        search_result_n line_search_result = fibonacci(function, x_curr, x_next, eps);
        accuracy = std::min(accuracy, line_search_result.accuracy);
        total_probes += line_search_result.function_calls;

        numerics::vector_f64 x_prev = x_curr;
        x_curr = line_search_result.result;
        x_next = x_curr;

        if (numerics::vector_f64::distance(x_curr, x_prev) < 2 * eps)
        {
            ++opt_coord_n;
            if (opt_coord_n == x_curr.size())
                break;
        }
        else
            opt_coord_n = 0;
    }

    return search_result_n(
        methods_types_n::per_coordinate_descend,// Вызываемый метод 
        iterations,                             // Количество итераций
        total_probes,                           // Количество вызовов функции
        accuracy,                               // Точность
        std::move(x_curr));                     // Экстремум
}

search_result_n gradient_descend(function_nd function, const numerics::vector_f64& x_start, const F64 eps, const I32 max_iters) {
    numerics::vector_f64 x_i(x_start), x_i_1, gradient;
    F64 step = 1.0;
    UI64 iterations = 0;
    UI64 total_probes = 0;
    F64 accuracy = std::numeric_limits<double>::infinity();

    for (; iterations <= max_iters; iterations++) {
        gradient = numerics::vector_f64::gradient(function, x_i, eps);
        x_i_1 = x_i - step * gradient;
        total_probes += 2 * x_i.size();

        search_result_n line_search_result = fibonacci(function, x_i, x_i_1, eps);
        accuracy = std::min(accuracy, line_search_result.accuracy);
        total_probes += line_search_result.function_calls;

        x_i_1 = line_search_result.result;

        if (numerics::vector_f64::distance(x_i_1, x_i) < 2 * eps)
            break;

        x_i = x_i_1;
    }
    return search_result_n(
        methods_types_n::gradient_descend,      // Вызываемый метод 
        iterations,                             // Количество итераций
        total_probes,                           // Количество вызовов функции
        accuracy,                               // Точность
        std::move(x_i_1));                      // Экстремум
}

search_result_n conj_gradient_descend(function_nd function, const numerics::vector_f64& x_start, const F64 eps, const I32 max_iters) {
    numerics::vector_f64 x_i(x_start), x_i_1;
    numerics::vector_f64 gradient = numerics::vector_f64::gradient(function, x_i, eps);
    numerics::vector_f64 gradient_old = (- 1.0) * gradient, gradient_new;
    F64 w;
    F64 step = 1.0;
    UI64 iterations = 0;
    UI64 total_probes = 2 * x_i.size();
    F64 accuracy = std::numeric_limits<double>::infinity();

    for (; iterations <= max_iters; iterations++) {
        x_i_1 = x_i + step * gradient_old;

        search_result_n line_search_result = fibonacci(function, x_i, x_i_1, eps);
        accuracy = std::min(accuracy, line_search_result.accuracy);
        total_probes += line_search_result.function_calls;

        x_i_1 = line_search_result.result;
        gradient_new = numerics::vector_f64::gradient(function, x_i_1, eps);
        total_probes += 2 * x_i.size();

        w = std::pow(gradient_new.magnitude(), 2) / std::pow(gradient_old.magnitude(), 2);
        gradient_old = gradient_old * w - gradient_new;

        if (numerics::vector_f64::distance(x_i_1, x_i) < 2 * eps) 
            break;
   
        x_i = x_i_1;
        gradient_old = gradient_new;
    }
    return search_result_n(
        methods_types_n::conj_gradient_descend, // Вызываемый метод 
        iterations,                             // Количество итераций
        total_probes,                           // Количество вызовов функции
        accuracy,                               // Точность
        std::move(x_i_1));                      // Экстремум
}

search_result_n newton_raphson(function_nd function, const numerics::vector_f64& x_start, const F64 eps, const I32 max_iters) {
    numerics::vector_f64 x_i(x_start), x_i_1, gradient;
    numerics::matrix_f64 hessian, hessian_1;
    F64 step = 1.0;
    UI64 iterations = 0;
    UI64 total_probes = 0;

    for (; iterations <= max_iters; iterations++) {
        gradient = numerics::vector_f64::gradient(function, x_i, eps);
        hessian = numerics::matrix_f64::hessian(function, x_i, eps);
        hessian_1 = numerics::matrix_f64::invert(hessian);
        x_i_1 = x_i - hessian_1 * gradient;
        total_probes += 2 * x_i.size() + x_i.size() * x_i.size();

        if (numerics::vector_f64::distance(x_i_1, x_i) < 2 * eps)
            break;

        x_i = x_i_1;
    }
    return search_result_n(                         
        methods_types_n::newton_raphson,            // Вызываемый метод 
        iterations,                                 // Количество итераций
        total_probes,                               // Количество вызовов функции
        numerics::vector_f64::distance(x_i, x_i_1), // Точность
        std::move(x_i_1));                          // Экстремум
}

search_result_n internal_penalty(function_nd function, const std::vector<function_nd>& constraints, const numerics::vector_f64& x_start, const F64 eps, const I32 max_iters) {
    numerics::vector_f64 x_i(x_start), x_i_1, gradient, search_point;
    F64 step = 1.0; 
    UI64 iterations = 0;
    UI64 total_probes = 0;

    for (; iterations < max_iters; iterations++) {
        auto penalty_function = [&](const numerics::vector_f64& x) -> double {
            double penalty = 0.0;

            for (const auto& constraint : constraints) {
                double fi = constraint(x);
                penalty += step / (-fi);
            }

            return function(x) + penalty;
            };

        numerics::vector_f64 gradient = numerics::vector_f64::gradient(penalty_function, x_i, eps);
        total_probes += 2 * x_i.size();

        numerics::vector_f64 search_point = x_i - gradient * 1.0;

        search_result_n result = fibonacci(penalty_function, x_i, search_point, eps);
        x_i_1 = result.result;
        total_probes += result.function_calls;

        step *= 0.5;

        F64 current_accuracy = numerics::vector_f64::distance(x_i_1, x_i);

        if (numerics::vector_f64::distance(x_i_1, x_i) < 2 * eps)
            break;

        x_i = x_i_1;
    }
    return search_result_n(
        methods_types_n::internal_penalty,          // Вызываемый метод 
        iterations,                                 // Количество итераций
        total_probes,                               // Количество вызовов функции
        numerics::vector_f64::distance(x_i_1, x_i), // Точность
        std::move(x_i_1));                          // Экстремум
}