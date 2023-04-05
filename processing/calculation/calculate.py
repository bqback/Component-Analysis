import scipy as sp
import numpy as np
import classes
import plotting
from .ica import run_ica
from typing import List, Dict


def value_at_point(gasses, coefficients, data_point):
    value = 0
    for i in range(len(gasses)):
        print(f'\t{gasses[i].name}: {coefficients[i] * gasses[i].KNa[data_point]}')
        value += coefficients[i] * gasses[i].KNa[data_point]
    return value


def normalize(coefficients, nh3_pos):
    norm_coefficients = [coeff / (coefficients[nh3_pos]*1e7) for coeff in coefficients]
    return norm_coefficients


def get_a(spectra):
    array = []
    for spectrum in spectra:
        array.append(spectrum.KNa)

    a_horizontal = np.array(array)
    a = np.transpose(a_horizontal)

    return a


def y_vals(values, coefficients):
    return values.dot(coefficients)


def calculate(data: Dict[str, classes.Spectrum | List[classes.Spectrum]]):
    a = get_a(data["gasses"])

    run_ica(a)

    b = np.array(data["sample"].KNa)

    res = sp.optimize.lsq_linear(a, b, bounds=(0, np.inf))
    solution_nnls, residuals_nnls = sp.optimize.nnls(a, b)

    nh3_idx = next((index for (index, d) in enumerate(data["gasses"]) if "NH3" in d.name), None)

    lsql = normalize(res.x, nh3_idx)
    lsql_ub = normalize(res.unbounded_sol[0], nh3_idx)
    nnls = normalize(solution_nnls, nh3_idx)

    print('Solution for coefficients')
    for i in range(len(data["gasses"])):
        print(f'\n{data["gasses"][i].name}:'
              f'\n\tlsq_linear [0, inf): {res.x[i]}, normalized around NH3 {lsql[i]}'
              f'\n\tlsq_linear unbounded: {res.unbounded_sol[0][i]}, normalized around NH3 {lsql_ub[i]}'
              f'\n\tnon-negative least squares: {solution_nnls[i]}, normalized around NH3 {nnls[i]}')

    print(f'\n\nSolving for 255th data point (highest peak)')

    print(f'\nlsq_linear solution')
    result_lsql = value_at_point(data["gasses"], res.x, 254)

    print(f'\nlsq_linear solution, normalized')
    result_lsql_norm = value_at_point(data["gasses"], lsql, 254)

    print(f'\nlsq_linear unbounded solution')
    result_lsql_ub = value_at_point(data["gasses"], res.unbounded_sol[0], 254)

    print(f'\nlsq_linear unbounded solution, normalized')
    result_lsql_ub_norm = value_at_point(data["gasses"], lsql_ub, 254)

    print(f'\nNNLS solution')
    result_nnls = value_at_point(data["gasses"], solution_nnls, 254)

    print(f'\nNNLS solution, normalized')
    result_nnls_norm = value_at_point(data["gasses"], nnls, 254)

    print(f'\nSample:\t\t\t\t{data["sample"].KNa[254]}')
    print(f'LSQ_L result:\t\t{result_lsql}\t\t normalized\t\t{result_lsql_norm}')
    print(f'LSQ_L_UB result:\t{result_lsql_ub}\t\t normalized\t\t{result_lsql_ub_norm}')
    print(f'NNLS result:\t\t{result_nnls}\t\t normalized\t\t{result_nnls_norm}')

    plotting.plot_data(data["sample"],
                       data["gasses"],
                       [res.x, res.unbounded_sol[0], solution_nnls],
                       ["least_squares linear", "least_squares unbounded", "non-negative least squares"])

    # test_sum = 0
    # for i in range(len(test)):
    #     test_sum += test[i]*solution[i]
    # print(test_sum)
    # print(test_solution)

