import processing
import constants
import classes
from typing import List, Dict
import os
import pathlib


def gather_data() -> Dict[str, classes.Spectrum | List[classes.Spectrum]]:
    spectra = []
    print(constants.DATA_PATH)
    for child in pathlib.Path(constants.DATA_PATH).iterdir():
        if child.suffix == ".txt":
            sample_name, raw_data = processing.parse(child)
            converted_data = processing.convert(raw_data)
            spectra.append(classes.Spectrum(sample_name, converted_data))
    sample = list(filter(lambda s: "model" in s.name, spectra))
    if not sample:
        raise Exception(f'Data directory /{constants.DATA_PATH} did not contain a valid model spectrum file')
    if len(sample) > 1:
        raise Exception(f'Multiple model spectra files found in data directory /{constants.DATA_PATH}')
    sample = sample[0]
    spectra.remove(sample)
    return {"sample": sample, "gasses": spectra}


def check_data(data: Dict[str, classes.Spectrum | List[classes.Spectrum]]) -> bool:
    valid = True
    for gas in data["gasses"]:
        if data["sample"].WN_min != gas.WN_min:
            print(f'Minimal wave number for spectra {gas.name} does not match \\'
                  f'the minimal wave number for the air sample \\'
                  f'(expected {data["sample"].WN_min}, got {gas.WN_min}')
            valid = False
            # raise ValueError(f'Minimal wave number for spectra {gas.name} does not match \\'
            #                  f'the minimal wave number for the air sample \\'
            #                  f'(expected {sample.WN_min}, got {gas.WN_min}')
        if data["sample"].WN_max != gas.WN_max:
            print(f'Maximum wave number for spectra {gas.name} does not match \\'
                  f'the maximum wave number for the air sample \\'
                  f'(expected {data["sample"].WN_max}, got {gas.WN_max}')
            valid = False
            # raise ValueError(f'Maximum wave number for spectra {gas.name} does not match \\'
            #                  f'the maximum wave number for the air sample \\'
            #                  f'(expected {sample.WN_max}, got {gas.WN_max}')
        if data["sample"].num_points != gas.num_points:
            print(f'Number of data points for spectra {gas.name} does not match \\'
                  f'the number of data points for the air sample \\'
                  f'(expected {data["sample"].num_points}, got {gas.num_points}')
            valid = False
            # raise ValueError(f'Number of data points for spectra {gas.name} does not match \\'
            #                  f'the number of data points for the air sample \\'
            #                  f'(expected {sample.num_points}, got {gas.num_points}')
    return valid


def main() -> None:
    spectra = gather_data()
    if not check_data(spectra):
        raise ValueError('Unable to proceed due to data mismatch')
    else:
        processing.calculate(spectra)


if __name__ == '__main__':
    main()
