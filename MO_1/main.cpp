#include "lr_1.h"
#include "lr2.h"
#include "simplex.h"
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

    search_result_n newton_raphson_result = newton_raphson(my_function, left, eps, max_iters);
    std::cout << "Newton Raphson Result:\n" << newton_raphson_result << std::endl;

    std::vector<function_nd> constraints;
    std::vector<function_nd> equality_constraints;

    // x ≥ 2
    constraints.push_back([](const numerics::vector_f64& x) -> double {
        return 2 - x[0];
        });
    // y ≥ 1  
    constraints.push_back([](const numerics::vector_f64& x) -> double {
        return 1 - x[1];
        });
    // x + y ≤ 7
    constraints.push_back([](const numerics::vector_f64& x) -> double {
        return x[0] + x[1] - 7;
        });
    // x + y = 6
    equality_constraints.push_back([](const numerics::vector_f64& x) -> double {
        return x[0] + x[1] - 6; 
        });

    numerics::vector_f64 start = { 3.0, 2.0 };

    search_result_n internal_result = internal_penalty(my_function, constraints, start, eps, max_iters);
    std::cout << "Internal Penalty Result:\n" << internal_result << std::endl;

    search_result_n external_result = external_penalty(my_function, constraints, equality_constraints, start, eps, max_iters);
    std::cout << "External Penalty Result:\n" << external_result << std::endl;
    
    return 0;
}

int lab4() {
    //numerics::vector_f64 c = { 2.0, -1.0 };  //2x - y
    //bool seek = true;
    //numerics::vector_f64 b = { 12.0, 8.0, 5.0, 1.0 };
    //int max_iters = 100;

    //numerics::matrix_f64 A({ -1, 3, 1, 1, 1, 0, 1, -2}, 4, 2);

    //numerics::vector_f64 c = { 1.0, -1.0 };  //x - y
    //bool seek = false;
    //numerics::vector_f64 b = { 2.0, 9.0, 27.0, 3.0 };
    //int max_iters = 100;

    //numerics::matrix_f64 A({ -3, 1, -1, 2, 4, 1, 1, -1 }, 4, 2);

    numerics::vector_f64 c = { 1.0, 0.0, 0.0 };  //x
    bool seek = true;
    numerics::vector_f64 b = { -2.0, 20.0 };
    int max_iters = 100;

    numerics::matrix_f64 A({ -1, -1, -1, 4, 5, 5 }, 2, 3);

    simplex_result sx_result = minimize(seek, c, A, b, max_iters);
    std::cout << "Simplex Result:\n" << sx_result << std::endl;

    return 0;
}

int main()
{
    //lab1();
    //lab2();
    //lab3();
    lab4();
    return 0;
}