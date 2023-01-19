from typing import List, Tuple


def convert(data_lines: List[str]) -> List[Tuple[float]]:
    floats = list(list(map(float, line.split())) for line in data_lines)
    floats_transposed = list(zip(*floats))
    return floats_transposed
