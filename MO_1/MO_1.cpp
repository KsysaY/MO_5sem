//#include <iostream>
//#include <cmath>
//#include <functional>
//#include <iomanip>
//
//using function_id = std::function<double(double)>;
//constexpr double PSI = 0.61803398874989484820;
//
//template<typename... params_t>
//const char* fstring(const char* format, const params_t&... params) {
//    static constexpr size_t buffer_size = 1024;
//    static char buffer[buffer_size];
//    sprintf_s(buffer, format, params...);
//    return buffer;
//}
//
//enum class methods_types : uint8_t {
//    bisect,
//    golden_ratio,
//    fibonacci
//};
//
//struct search_result {
//    friend std::ostream& operator<<(std::ostream& stream, const search_result& res);
//
//    const methods_types method;
//    const size_t iterations;
//    const size_t function_calls;
//    const double accuracy;
//    const double result;
//
//    search_result(methods_types m, size_t iters, size_t calls, double acc, double res)
//        : method(m),
//        iterations(iters),
//        function_calls(calls),
//        accuracy(acc),
//        result(res) {
//    }
//};
//
//std::ostream& operator<<(std::ostream& stream, const search_result& result) {
//    switch (result.method) {
//    case methods_types::bisect:
//        stream << fstring("method: bisect\n");
//        break;
//    case methods_types::golden_ratio:
//        stream << fstring("method: golden_ratio\n");
//        break;
//    case methods_types::fibonacci:
//        stream << fstring("method: fibonacci\n");
//        break;
//    default:
//        stream << fstring("method: Unknown\n");
//    }
//    stream << fstring("iterations: %ld\n", result.iterations);
//    stream << fstring("function_calls: %ld\n", result.function_calls);
//    stream << fstring("accuracy: %.10f\n", result.accuracy);
//    stream << fstring("result: %.10f\n", result.result);
//    return stream;
//}
//
//search_result bisect(function_id func, double lhs, double rhs, double eps = 1e-9, size_t max_iters = 1000);
//
//search_result bisect(function_id func, double lhs, double rhs, double eps, size_t max_iters) {
//    size_t iter = 0;
//    double center;
//    size_t func_iter = 0;
//
//    if (lhs > rhs) 
//        std::swap(lhs, rhs);
//
//    for (; iter < max_iters && std::abs(rhs - lhs) > 2 * eps; ++iter) {
//        center = (lhs + rhs) * 0.5;
//        func_iter++;
//        if (func(center + eps) > func(center - eps)) {
//            rhs = center;
//        }
//        else {
//            lhs = center;
//        }
//    }
//
//    return search_result(methods_types::bisect, iter, func_iter, std::abs(lhs - rhs), (lhs + rhs) * 0.5);
//}
//
//search_result golden_ratio(function_id func, double left, double right, const double eps, const int max_iterations) {
//    if (left > right) 
//        std::swap(left, right);
//    int iteration = 0;
//    double x_l = right - (right - left) * PSI;
//    double x_r = left + (right - left) * PSI;
//    double f_l = func(x_l);
//    double f_r = func(x_r);
//    size_t func_calls = 2;
//    for (; iteration != max_iterations && std::abs(right - left) > 2 * eps; iteration++) {
//        if (f_l > f_r) {
//            left = x_l;
//            x_l = x_r;
//            f_l = f_r;
//            x_r = left + (right - left) * PSI;
//            f_r = func(x_r);
//            func_calls++; 
//        }
//        else {
//            right = x_r;
//            x_r = x_l;
//            f_r = f_l;
//            x_l = right - (right - left) * PSI;
//            f_l = func(x_l);
//            func_calls++;
//        }
//    }
//    return search_result(methods_types::golden_ratio, iteration, func_calls, std::abs(right - left), (right + left) * 0.5);
//}
//
//search_result fibonacci(function_id func, double left, double right, const double eps) {
//    if (left > right)
//        std::swap(left, right);
//
//    double acc = (right - left) / eps;
//    double fib, fib_n_1 = 1, fib_n = 1;
//    int iter = 0;
//
//    while (fib_n < acc) {
//        fib = fib_n_1;
//        fib_n_1 = fib_n;
//        fib_n += fib;
//        iter++;
//    }
//
//    double x_r = left + (fib_n_1 / fib_n) * (right - left);
//    double x_l = left + ((fib_n - fib_n_1) / fib_n) * (right - left);
//
//    double f_l = func(x_l);
//    double f_r = func(x_r);
//    int func_calls = 2;
//
//    double delta = (right - left) / 100.0;
//    for (int i = iter; i > 0; i--) {
//        fib = fib_n - fib_n_1;
//        fib_n = fib_n_1;
//        fib_n_1 = fib;
//
//        if (f_l > f_r) {
//            left = x_l;
//            f_l = f_r;
//            x_l = x_r;
//            x_r = left + (fib_n_1 / fib_n) * (right - left);
//
//           /* if (std::abs(x_r - x_l) < delta) {
//                x_r = x_r + delta;
//            }*/
//            f_r = func(x_r);
//            func_calls++;
//        }
//        else {
//            right = x_r;
//            x_r = x_l;
//            f_r = f_l;
//            x_l = left + ((fib_n - fib_n_1) / fib_n) * (right - left);
//
//            /*if (std::abs(x_r - x_l) < delta) {
//                x_l = x_l - delta;
//            }*/
//            f_l = func(x_l);
//            func_calls++;
//        }
//    }
//
//    double x_min = (left + right) * 0.5;
//    double accuracy = std::abs(right - left) * 0.5;
//    return search_result(methods_types::fibonacci, iter, func_calls, accuracy, x_min);
//}
//
//int main() {
//    function_id my_function = [](double x) { return (x - 1) * (x - 3); };
//    double left = 0.0;
//    double right = 3.0;
//    double eps = 1e-6;
//    int max_iters = 100;
//
//    search_result bisect_result = bisect(my_function, left, right, eps, max_iters);
//    std::cout << "Bisect Result:\n" << bisect_result << std::endl;
//
//    search_result golden_ratio_result = golden_ratio(my_function, left, right, eps, max_iters);
//    std::cout << "Golden Ratio Result:\n" << golden_ratio_result << std::endl;
//
//    search_result fibonacci_result = fibonacci(my_function, left, right, eps);
//    std::cout << "Fibonacci Result:\n" << fibonacci_result << std::endl;
//
//    return 0;
//}