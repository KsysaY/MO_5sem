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

numerics::matrix_f64 reference_plan(const numerics::vector_f64& c, const numerics::matrix_f64& A, const numerics::vector_f64& b, const I32 max_iterations, UI64& iterations) {
    UI64 m = A.rows_count(); //строка
    UI64 n = A.cols_count(); //столбец
    numerics::vector_f64 ksi(b.size());
    UI64 ind_i = 0;
    UI64 ind_j = 0;

    numerics::matrix_f64 table = numerics::matrix_f64::zeros(m + 2, n + 1 + b.size());
    numerics::matrix_f64 new_table = numerics::matrix_f64::zeros(m + 1, n + 1);

    table(1, 0) = 0.0;

    for (UI64 i = 0; i < m; ++i) {
        if (b[i] < 0) {
            table(i + 2, 0) = -b[i];

            for (UI64 j = 0; j < n; ++j) {
                table(i + 2, j + 1) = -A(i, j);
            }
        }
        else {
            table(i + 2, 0) = b[i];

            for (UI64 j = 0; j < n; ++j) {
                table(i + 2, j + 1) = A(i, j);
            }
        }
    }

    for (UI64 i = 0; i < m; ++i) {
        if (b[i] < 0)
            table(i + 2, i + n + 1) = -1;
        else
            table(i + 2, i + n + 1) = 1;
    }

    for (UI64 j = 0; j < c.size(); ++j) {
        table(1, j + 1) = c[j];
    }

    for (UI64 j = 0; j < n + 1 + b.size(); ++j) {
        for (UI64 i = 0; i < m; ++i) {
            table(0, j) += table(i + 2, j);
        }
    }

    for (iterations = 0; iterations < max_iterations; ++iterations) {

        F64 min_coeff = 0.0;
        bool optimal = true;

        for (UI64 j = 1; j < n + 1 + b.size(); ++j) {
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

        ksi[iterations] = ind_j;

        F64 min_ratio = std::numeric_limits<F64>::max();

        for (UI64 i = 2; i < m + 2; ++i) {
            if (table(i, ind_j) > 1e-10) {
                F64 ratio = table(i, 0) / table(i, ind_j);
                if (ratio >= 0 && ratio < min_ratio) {
                    min_ratio = ratio;
                    ind_i = i;
                }
            }
        }

        for (UI64 i = 0; i < m + 2 ; ++i) {
            if (i != ind_i) {
                for (UI64 j = 0; j < n + 1 + b.size(); ++j) {
                    if (j != ind_j) {
                        table(i, j) = table(i, j) - table(i, ind_j) * table(ind_i, j) / table(ind_i, ind_j);
                    }
                }
            }
        }

        for (UI64 i = 0; i < m + 2; ++i) {
            if (i != ind_i) {
                table(i, ind_j) = -table(i, ind_j) / table(ind_i, ind_j);
            }
        }

        for (UI64 j = 0; j < n + 1 + b.size(); ++j) {
            if (j != ind_j) {
                table(ind_i, j) = table(ind_i, j) / table(ind_i, ind_j);
            }
        }

        table(ind_i, ind_j) = 1 / table(ind_i, ind_j);
    }
    // удалить лишние строки + добавить базис + проверить интерации


    return new_table;
}

simplex_result minimize(bool seek, const numerics::vector_f64& c, const numerics::matrix_f64& A, const numerics::vector_f64& b, const I32 max_iterations) 
{
    UI64 m = A.rows_count(); //строка
    UI64 n = A.cols_count(); //столбец
    UI64 iterations = 0;
    numerics::vector_f64 res(c.size(), 0.0);
    numerics::vector_f64 basis(c.size());
    F64 value = 0.0;
    numerics::vector_f64 c_copy(c);
    UI64 ind_i = 0;
    UI64 ind_j = 0;
    bool neg = true;

    if (!seek) {
        for (UI64 i = 0; i < c.size(); ++i) {
            c_copy[i] = -c[i];
        }
    }

    //Симплекс-таблица
    numerics::matrix_f64 table = numerics::matrix_f64::zeros(m + 1, n + 1);

    //Опорный план
    for (UI64 i = 0; i < b.size(); ++i) {
        if (b[i] < 0) {
            neg = false;
        }
    }

    if (!neg) {
        table = reference_plan(c_copy, A, b, max_iterations, iterations);
    }
    else {
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
    }

    if (seek)
        value = -table(0, 0);
    else
        value = table(0, 0);

    for (UI64 j = 0; j < c.size(); ++j) {
        for (UI64 i = 1; i <= m; ++i) {
            if (basis[j] == i) {
                res[j] = table(i, 0);
            }
        }
    }

    return simplex_result(value, res, iterations);
}
