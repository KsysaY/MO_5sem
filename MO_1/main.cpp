#include "lr_1.h"
#include "lr2.h"
#include <iostream>

int lab1() {
    function_1d my_function = [](double x) { return (x - 1) * (x - 2); };
    double left = 0.0;
    double right = 3.0;
    double eps = 1e-6;
    int max_iters = 100;


    search_result bisect_result = bisect(my_function, left, right, eps, max_iters);
    std::cout << "Bisect Result:\n" << bisect_result << std::endl;

    search_result golden_ratio_result = golden_ratio(my_function, left, right, eps, max_iters);
    std::cout << "Golden Ratio Result:\n" << golden_ratio_result << std::endl;

    search_result fibonacci_result = fibonacci(my_function, left, right, eps );//офигеть, чтобы было одинаковое количество функций: eps*2.0
    std::cout << "Fibonacci Result:\n" << fibonacci_result << std::endl;

    return 0;
}

int lab2() {
    function_nd my_function = [](const numerics::vector_f64& args) -> double {
        return (args[0] - 5) * args[0] + (args[1] - 3) * args[1];
    };
    numerics::vector_f64 left = { 0,0 };
    numerics::vector_f64 right = { 5,3 };
    double eps = 1e-6;
    int max_iters = 100;


    search_result_n bisect_result = bisect(my_function, left, right, eps, max_iters);
    std::cout << "Bisect Result:\n" << bisect_result << std::endl;

    search_result_n golden_ratio_result = golden_ratio(my_function, left, right, eps, max_iters);
    std::cout << "Golden Ratio Result:\n" << golden_ratio_result << std::endl;

    search_result_n fibonacci_result = fibonacci(my_function, left, right, eps);
    std::cout << "Fibonacci Result:\n" << fibonacci_result << std::endl;

    search_result_n per_coord_descend_result = per_coord_descend(my_function, left, eps, max_iters);
    std::cout << "Per Coord Descend Result:\n" << per_coord_descend_result << std::endl;

    return 0;
}

int lab3() {
    function_nd my_function = [](const numerics::vector_f64& args) -> double {
        return (args[0] - 5) * args[0] + (args[1] - 3) * args[1];
        };
    numerics::vector_f64 left = { 0,0 };
    numerics::vector_f64 right = { 5,3 };
    double eps = 1e-6;
    int max_iters = 100;

    search_result_n gradient_result = gradient_descend(my_function, left, eps, max_iters);
    std::cout << "Gradient Descend Result:\n" << gradient_result << std::endl;

    search_result_n conj_gradient_result = conj_gradient_descend(my_function, left, eps, max_iters);
    std::cout << "Conjugate Gradient Result:\n" << conj_gradient_result << std::endl;

    return 0;
}

int main()
{
    //lab1();
    lab2();
    lab3();
    return 0;
}