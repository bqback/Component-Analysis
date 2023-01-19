from typing import List
from pathlib import Path


def parse(path: str, name: str) -> List:
    dir_path = Path(path)
    file_path = dir_path / name
    with file_path.open('r') as file:
        lines = list(filter(None, (line.lstrip().rstrip() for line in file)))
        trimmed_lines = list(filter(lambda s: s[0] != "#", lines))
    return trimmed_lines
