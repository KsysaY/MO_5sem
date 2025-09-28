#include "lr_1.h"


std::ostream& operator<<(std::ostream& stream, const search_result& result) {
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

search_result bisect(function_1d func, double lhs, double rhs, double eps, size_t max_iterations) {
    size_t iter = 0;
    double center;
    size_t func_iter = 0;

    if (lhs > rhs)
        std::swap(lhs, rhs);

    for (; iter < max_iterations && std::abs(rhs - lhs) > 2 * eps; ++iter) {
        center = (lhs + rhs) * 0.5;
        func_iter+=2;
        if (func(center + eps) > func(center - eps)) {
            rhs = center;
        }
        else {
            lhs = center;
        }
    }

    return search_result(methods_types::bisect, iter, func_iter, std::abs(lhs - rhs), (lhs + rhs) * 0.5);
}

search_result golden_ratio(function_1d func, double lhs, double rhs, const double eps, const size_t max_iterations) {
    if (lhs > rhs)
        std::swap(lhs, rhs);
    int iteration = 0;
    double x_l = rhs - (rhs - lhs) * PSI;
    double x_r = lhs + (rhs - lhs) * PSI;
    double f_l = func(x_l);
    double f_r = func(x_r);
    size_t func_calls = 2;
    for (; iteration != max_iterations && std::abs(rhs - lhs) > 2 * eps; iteration++) {
        if (f_l > f_r) {
            lhs = x_l;
            x_l = x_r;
            f_l = f_r;
            x_r = lhs + (rhs - lhs) * PSI;
            f_r = func(x_r);
            func_calls++;
        }
        else {
            rhs = x_r;
            x_r = x_l;
            f_r = f_l;
            x_l = rhs - (rhs - lhs) * PSI;
            f_l = func(x_l);
            func_calls++;
        }
    }
    return search_result(methods_types::golden_ratio, iteration, func_calls, std::abs(rhs - lhs), (rhs + lhs) * 0.5);
}

search_result fibonacci(function_1d func, double lhs, double rhs, const double eps) {
    if (lhs > rhs)
        std::swap(lhs, rhs);

    double acc = (rhs - lhs) / eps;
    double fib, fib_n_1 = 1, fib_n = 1;
    int iter = 0;

    while (fib_n < acc) {
        fib = fib_n_1;
        fib_n_1 = fib_n;
        fib_n += fib;
        iter++;
    }

    double x_r = lhs + (fib_n_1 / fib_n) * (rhs - lhs);
    double x_l = lhs + ((fib_n - fib_n_1) / fib_n) * (rhs - lhs);

    double f_l = func(x_l);
    double f_r = func(x_r);
    int func_calls = 2;

    double delta = (rhs - lhs) / 100.0;
    for (int i = iter; i > 0; i--) {
        fib = fib_n - fib_n_1;
        fib_n = fib_n_1;
        fib_n_1 = fib;

        if (f_l > f_r) {
            lhs = x_l;
            f_l = f_r;
            x_l = x_r;
            x_r = lhs + (fib_n_1 / fib_n) * (rhs - lhs);

            /* if (std::abs(x_r - x_l) < delta) {
                 x_r = x_r + delta;
             }*/
            f_r = func(x_r);
            func_calls++;
        }
        else {
            rhs = x_r;
            x_r = x_l;
            f_r = f_l;
            x_l = lhs + ((fib_n - fib_n_1) / fib_n) * (rhs - lhs);

            /*if (std::abs(x_r - x_l) < delta) {
                x_l = x_l - delta;
            }*/
            f_l = func(x_l);
            func_calls++;
        }
    }

    return search_result(methods_types::fibonacci, iter, func_calls, std::abs(rhs - lhs) * 0.5, (lhs + rhs) * 0.5);
}