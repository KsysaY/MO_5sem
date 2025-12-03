#include "lr2.h"
#include "lr_1.h"
#include "simplex.h"
#include "numerics/linalg/numeric_vector.h"

std::ostream& operator<<(std::ostream& stream, const simplex_result& result) {

    stream << fstring("value:          %f\n", result.value);
    stream << fstring("iterations:     %ld\n", result.iterations);

    stream << "result: ";
    stream << "[";
    for (size_t i = 0; i < result.solution.size(); ++i) {
        stream << std::setprecision(10) << result.solution[i];
        if (i < result.solution.size() - 1) {
            stream << ", ";
        }
    }
    stream << "]\n";

    return stream;
}

simplex_result minimize(bool seek, const numerics::vector_f64& c, const numerics::matrix_f64& A, const numerics::vector_f64& b, const I32 max_iterations) 
{
    UI64 m = A.rows_count(); //строка
    UI64 n = A.cols_count(); //столбец
    UI64 iterations = 0;
    numerics::vector_f64 res(c.size(), 0.0);
    numerics::vector_f64 basis(c.size());
    F64 optimal_value = 0.0;
    numerics::vector_f64 c_copy(c);
    UI64 ind_i = 0;
    UI64 ind_j = 0;

    if (!seek) {
        for (UI64 i = 0; i < c.size(); ++i) {
            c_copy[i] = -c[i];
        }
    }

    //Симплекс-таблица
    numerics::matrix_f64 table = numerics::matrix_f64::zeros(m + 1, n + 1);

    table(0, 0) = 0.0;

    for (UI64 j = 0; j < n; ++j) {
        table(0, j + 1) = c_copy[j];
    }

    for (UI64 i = 0; i < m; ++i) { 
        table(i + 1, 0) = b[i];

        for (UI64 j = 0; j < n; ++j) {
            table(i + 1, j + 1) = A(i, j);
        }
    }

    for (iterations = 0; iterations < max_iterations; ++iterations) {

        F64 min_coeff = 0.0;
        bool optimal = true;

        for (UI64 j = 1; j <= n ; ++j) {
            if (table(0, j) > 1e-10) {
                optimal = false;
                if (table(0, j) > min_coeff) {
                    min_coeff = table(0, j);
                    ind_j = j;
                }
            }
        }

        if (optimal) {
            break;
        }

        if (ind_j == 0) {
            break;
        }

        F64 min_ratio = std::numeric_limits<F64>::max();

        for (UI64 i = 1; i <= m; ++i) {
            if (table(i, ind_j) > 1e-10) {
                F64 ratio = table(i, 0) / table(i, ind_j);
                if (ratio >= 0 && ratio < min_ratio) {
                    min_ratio = ratio;
                    ind_i = i;
                }
            }
        }

        basis[ind_j - 1] = ind_i;

        for (UI64 i = 0; i <= m; ++i) {
            if (i != ind_i) {
                for (UI64 j = 0; j <= n; ++j) {
                    if (j != ind_j) {
                        table(i, j) = table(i, j) - table(i, ind_j) * table(ind_i, j) / table(ind_i, ind_j);
                    }
                }
            }
        }

        for (UI64 i = 0; i <= m; ++i) {
            if (i != ind_i) {
                table(i, ind_j) = - table(i, ind_j) / table(ind_i, ind_j);
            }
        }

        for (UI64 j = 0; j < n; ++j) {
            if (j != ind_j) {
                table(ind_i, j) = table(ind_i, j) / table(ind_i, ind_j);
            }
        }

        table(ind_i, ind_j) = 1 / table(ind_i, ind_j);

        //for (UI64 i = 0; i <= m; ++i) {
        //    // Метка строки
        //    if (i == 0) {
        //        std::cout << std::setw(10) << "Z";
        //    }
        //    else {
        //        std::cout << std::setw(10) << "x" << (n + i);
        //    }

        //    // Значения в строке
        //    for (UI64 j = 0; j <= n; ++j) {
        //        std::cout << std::setw(12) << std::setprecision(4) << std::fixed << table(i, j);
        //    }
        //    std::cout << "\n";
        //}      
    }

    if (seek)
        optimal_value = -table(0, 0);

    for (UI64 j = 0; j < c.size(); ++j) {
        for (UI64 i = 1; i <= m; ++i) {
            if (basis[j] == i) {
                res[j] = table(i, 0);
                std::cout << res[j];
            }
        }
    }

    return simplex_result(optimal_value, res, iterations);
}
