import numpy as np
from classes import Spectrum
from typing import List


def calculate(data: List[Spectrum]):
    array = []
    test = []
    test_solution = data[0].KNa[0]
    for gas in data[1:]:
        array.append(gas.KNa)
        test.append(gas.KNa[0])
    a_horizontal = np.array(array)
    a = np.transpose(a_horizontal)
    b = data[0].KNa
    print(test)
    print(test_solution)
    solution, sol_error, _rank, _singvalues = np.linalg.lstsq(a, b)
    for i in range(len(data[1:])):
        print(f'{data[i+1].name}: {solution[i]}')
    print(solution)
    print(sol_error)

