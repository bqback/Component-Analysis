import processing
import constants
import classes
from typing import List


def gather_data(name: str) -> classes.Spectrum:
    raw_data = processing.parse(constants.DATA_PATH, f'{name}.txt')
    converted_data = processing.convert(raw_data)
    return classes.Spectrum(name, converted_data)


def check_data(data: List[classes.Spectrum]) -> bool:
    valid = True
    sample = data[0]
    for gas in data[1:]:
        if sample.WN_min != gas.WN_min:
            print(f'Minimal wave number for spectra {gas.name} does not match \\'
                  f'the minimal wave number for the air sample \\'
                  f'(expected {sample.WN_min}, got {gas.WN_min}')
            valid = False
            # raise ValueError(f'Minimal wave number for spectra {gas.name} does not match \\'
            #                  f'the minimal wave number for the air sample \\'
            #                  f'(expected {sample.WN_min}, got {gas.WN_min}')
        if sample.WN_max != gas.WN_max:
            print(f'Maximum wave number for spectra {gas.name} does not match \\'
                  f'the maximum wave number for the air sample \\'
                  f'(expected {sample.WN_max}, got {gas.WN_max}')
            valid = False
            # raise ValueError(f'Maximum wave number for spectra {gas.name} does not match \\'
            #                  f'the maximum wave number for the air sample \\'
            #                  f'(expected {sample.WN_max}, got {gas.WN_max}')
        if sample.num_points != gas.num_points:
            print(f'Number of data points for spectra {gas.name} does not match \\'
                  f'the number of data points for the air sample \\'
                  f'(expected {sample.num_points}, got {gas.num_points}')
            valid = False
            # raise ValueError(f'Number of data points for spectra {gas.name} does not match \\'
            #                  f'the number of data points for the air sample \\'
            #                  f'(expected {sample.num_points}, got {gas.num_points}')
    return valid


def main() -> None:
    spectra = []
    air_model = gather_data(constants.MODEL)
    spectra.append(air_model)
    for gas in constants.GAS_LIST:
        spectrum = gather_data(gas)
        spectra.append(spectrum)
    if not check_data(spectra):
        raise ValueError('Unable to proceed due to data mismatch')
    else:
        processing.calculate(spectra)


if __name__ == '__main__':
    main()
