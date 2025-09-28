#include "lr_1.h"
#include <iostream>

int main() {
    function_1d my_function = [](double x) { return (x - 1) * (x - 3); };
    function_1d my_function1 = [](double x) { return (x - 2)*(x-2)+1; };
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