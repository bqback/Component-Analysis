from typing import List, Tuple
import numpy as np


def convert(data_lines: List[str]) -> np.ndarray:
    floats = list(list(map(float, line.split())) for line in data_lines)
    floats_transposed = list(zip(*floats))
    return np.array(floats_transposed)
