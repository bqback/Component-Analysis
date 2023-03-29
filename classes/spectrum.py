from dataclasses import dataclass, field, InitVar


@dataclass
class Spectrum:
    name: str
    data: InitVar[list[tuple[float]]]
    WN: list[float] = field(init=False)
    KNa: list[float] = field(init=False)
    Ka: list[float] = field(init=False)
    WN_min: float = field(init=False)
    WN_max: float = field(init=False)
    WN_step: float = field(init=False)
    num_points: float = field(init=False)

    def __post_init__(self, data):
        self.WN, self.KNa, self.Ka = data
        self.WN_min = min(self.WN)
        self.WN_max = max(self.WN)
        self.WN_step = self.WN[1] - self.WN[0]
        self.num_points = len(self.WN)

    def __str__(self):
        header = f"{self.name}\nWN\t\tKNa\t\tKa\n"
        # data = "\n".join(
        #     (f"{self.WN[idx]}\t{self.KNa[idx]}\t\t{self.Ka[idx]}" for idx, _ in enumerate(self.WN))
        # )
        return header
