#pragma once
#include "numerics/numerics.h"
#include "defines.h"

struct simplex_result {
    friend std::ostream& operator<<(std::ostream& stream, const simplex_result& res);

    F64 value;                       //Значение целевой функции
    numerics::vector_f64 solution;   //Решение
    UI64 iterations;                 //Количество итераций

    simplex_result(F64 val, const numerics::vector_f64& sol, UI64 iters)
        : value(val), solution(sol), iterations(iters) { }
};

simplex_result minimize(bool seek, const numerics::vector_f64& c, const numerics::matrix_f64& A, const numerics::vector_f64& b, const I32 max_iterations);