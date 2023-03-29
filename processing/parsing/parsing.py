from typing import List, Tuple
from pathlib import Path


def parse(filename: Path) -> Tuple[str, List]:
    with filename.open('r') as file:
        lines = list(filter(None, (line.lstrip().rstrip() for line in file)))
        name = lines[1].lstrip("# ")
        data = list(filter(lambda s: s[0] != "#", lines))
    return name, data

