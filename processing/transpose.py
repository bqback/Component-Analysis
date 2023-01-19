from typing import List


def transpose(lines: List) -> List[List]:
    floats = list(map(float, (line.split() for line in lines)))
    return floats
