#include "lr2.h"
#include "lr_1.h"
#include "simplex.h"
#include "numerics/linalg/numeric_vector.h"

std::ostream& operator<<(std::ostream& stream, const simplex_result& result) {

    stream << fstring("value:          %ld\n", result.value);
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
    numerics::vector_f64 res;
    F64 optimal_value = 0.0;
    numerics::vector_f64 b_copy(b);

    if (!seek) {
        for (UI64 i = 0; i < b.size(); ++i) {
            b_copy[i] = -b[i];
        }
    }

    //Симплекс-таблица
    numerics::matrix_f64 table = numerics::matrix_f64::zeros(m + 1, n + m + 1);

    table(0, 0) = 0.0;

    for (UI64 j = 0; j < n; ++j) {
        table(0, j + 1) = -c[j];
    }

    for (UI64 i = 0; i < m; ++i) {
        for (UI64 j = 1; j < n; ++j) {
            table(i + 1, j) = A(i, j);
        }
        table(i + 1, 0) = b[i];
    }

    for (iterations = 0; iterations < max_iterations; ++iterations) {

        UI64 ind_j = 0;
        F64 min_coeff = 0.0;

        for (UI64 j = 0; j < n ; ++j) {
            if (table(0, j + 1) > min_coeff) {
                min_coeff = table(0, j + 1);
                ind_j = j + 1;
            }
        }

        if (ind_j == 0) {
            break;
        }

        UI64 ind_i = m + 1;
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

        res[ind_j] = table(0, ind_j);

        table(ind_i, ind_j) = 1 / table(ind_i, ind_j);

        for (UI64 i = 0; i <= m; ++i) {
            if (i != ind_i) {
                table(i, ind_j) = - table(i, ind_j) / table(ind_i, ind_j);
            }
        }

        for (UI64 j = 0; j < n; ++j) {
            if (j != ind_j) {
                table(ind_i, j) = -table(ind_i, j) / table(ind_i, ind_j);
            }
        }

        for (UI64 i = 0; i < m; ++i) {
            for (UI64 j = 1; j < n; ++j) {
                if (i != ind_i or j != ind_j) {
                    table(i, j) = table(i, j) - (table(i, ind_j) * table(ind_i, j) / table(ind_i, ind_j));
                }
            }
        }
        if (!seek)
            optimal_value = table(0, 0);
        else
            optimal_value = -table(0, 0);
    }

    return simplex_result(optimal_value, res, iterations);
}
