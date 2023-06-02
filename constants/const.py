DATA_PATH = "data/0-4000_1"
GAS_LIST = ["H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2", "SO2", "NH3"]
MODEL = "model"
REQUIRED_SAMPLES = 15
SNR = 25
TRUE_SOLUTION = {
    "[Non-airy] IAO model, high latitude, winter, H=0": {
        "Pure H2O": 0.000654836697,
        "Pure CO2": 0.000334962550,
        "Pure O3":  0.000000020118,
        "Pure N2O": 0.000000281650,
        "Pure CO":  0.000000472770,
        "Pure CH4": 0.000001488722,
        "Pure O2":  0.208219963045,
        "Pure SO2": 0.000000080471,
        "Pure NH3": 0.000000010059,
        "Pure N2":  0.790787883918
    },
    "[Non-airy] IAO model, high latitude, summer, H=0": {
        "Pure H2O": 0.006510000210,
        "Pure CO2": 0.000333000004,
        "Pure O3":  0.000000020000,
        "Pure N2O": 0.000000280000,
        "Pure CO":  0.000000470000,
        "Pure CH4": 0.000001480000,
        "Pure O2":  0.207000002265,
        "Pure SO2": 0.000000080000,
        "Pure NH3": 0.000000010000,
        "Pure N2":  0.786154657521
    },
    "[Non-airy] IAO model, mean latitude, winter, H=0": {
        "Pure H2O": 0.005880000070,
        "Pure CO2": 0.000333000004,
        "Pure O3":  0.000000020000,
        "Pure N2O": 0.000000280000,
        "Pure CO":  0.000000470000,
        "Pure CH4": 0.000001480000,
        "Pure O2":  0.207000002265,
        "Pure SO2": 0.000000080000,
        "Pure NH3": 0.000000010000,
        "Pure N2":  0.786784657661
    },
    "[Non-airy] IAO model, mean latitude, summer, H=0": {
        "Pure H2O": 0.015599999577,
        "Pure CO2": 0.000333000004,
        "Pure O3":  0.000000020000,
        "Pure N2O": 0.000000280000,
        "Pure CO":  0.000000470000,
        "Pure CH4": 0.000001480000,
        "Pure O2":  0.207000002264,
        "Pure SO2": 0.000000080000,
        "Pure NH3": 0.000000010000,
        "Pure N2":  0.777064658155
    },
    "[Non-airy] IAO model, tropics, H=0": {
        "Pure H2O": 0.023800000000,
        "Pure CO2": 0.000333000000,
        "Pure O3":  0.000000020000,
        "Pure N2O": 0.000000280000,
        "Pure CO":  0.000000470000,
        "Pure CH4": 0.000001480000,
        "Pure O2":  0.207000000000,
        "Pure SO2": 0.000000080000,
        "Pure NH3": 0.000000010000,
        "Pure N2":  0.768864660000
    }
}
CONSTRAINT_TYPES = {
    0: "unconstrained",
    1: "constrained from above",
    -1: "constrained from below"
}
