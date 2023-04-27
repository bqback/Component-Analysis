import processing
import constants
import classes
from typing import List, Dict
import os
import pathlib


def gather_data() -> Dict[str, classes.Spectrum | List[classes.Spectrum]]:
    sources = []
    samples = []
    print(f"Data sourced from {constants.DATA_PATH}")
    for child in pathlib.Path(constants.DATA_PATH).iterdir():
        if child.suffix == ".txt":
            sample_name, raw_data = processing.parse(child)
            converted_data = processing.convert(raw_data)
            if "model" in sample_name:
                samples.append(classes.Spectrum(sample_name, converted_data))
            else:
                sources.append(classes.Spectrum(sample_name, converted_data))
            if child.name != sample_name + ".txt":
                child.rename(child.with_name(sample_name + ".txt"))
    if len(samples) == 0:
        raise Exception(f'Data directory /{constants.DATA_PATH} did not contain a valid model spectrum file')
    return {"samples": samples, "sources": sources}


def check_data(data: Dict[str, classes.Spectrum | List[classes.Spectrum]]) -> bool:
    valid = True
    for gas in data["sources"]:
        if data["samples"][0].WN_min != gas.WN_min:
            print(f'Minimal wave number for spectra {gas.name} does not match \\'
                  f'the minimal wave number for the air sample \\'
                  f'(expected {data["sample"].WN_min}, got {gas.WN_min}')
            valid = False
            # raise ValueError(f'Minimal wave number for spectra {gas.name} does not match \\'
            #                  f'the minimal wave number for the air sample \\'
            #                  f'(expected {sample.WN_min}, got {gas.WN_min}')
        if data["samples"][0].WN_max != gas.WN_max:
            print(f'Maximum wave number for spectra {gas.name} does not match \\'
                  f'the maximum wave number for the air sample \\'
                  f'(expected {data["sample"].WN_max}, got {gas.WN_max}')
            valid = False
            # raise ValueError(f'Maximum wave number for spectra {gas.name} does not match \\'
            #                  f'the maximum wave number for the air sample \\'
            #                  f'(expected {sample.WN_max}, got {gas.WN_max}')
        if data["samples"][0].num_points != gas.num_points:
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
    # fetch_info()
    # H2O: 0.588000
    # CO2: 0.033300
    # O3: 0.000002
    # N2O: 0.000028
    # CO: 0.000047
    # CH4: 0.000148   O2: 20.700000   SO2: 0.000008   NH3: 0.000001   N2: 78.678466
    # H2O_conc = 588e-3
    # CO2_conc = 333e-4
    # O3_conc = 2e-6
    # N2O_conc = 28e-6
    # CO_conc = 47e-6
    # CH4_conc = 148e-6
    # O2_conc = 20.7
    # SO2_conc = 8e-6
    # NH3_conc = 1e-6
    # N2_conc = 100 - sum(
    #     [H2O_conc, CO2_conc, O3_conc, N2O_conc, CO_conc, CH4_conc, O2_conc, SO2_conc, NH3_conc]
    # )
    #
    # db_begin('data')
    # nu_H2O, coef_H2O = absorptionCoefficient_Voigt(
    #     SourceTables='H2O', Diluent={'self': 1.0}, File='absCoefH2O.txt'
    # )
    # nu_CO2, coef_CO2 = absorptionCoefficient_Voigt(
    #     SourceTables='CO2', Diluent={'self': 1.0}, File='absCoefCO2.txt'
    # )
    # nu_O3, coef_O3 = absorptionCoefficient_Voigt(
    #     SourceTables='O3', Diluent={'self': 1.0}, File='absCoefO3.txt'
    # )
    # nu_N2O, coef_N2O = absorptionCoefficient_Voigt(
    #     SourceTables='N2O', Diluent={'self': 1.0}, File='absCoefN2O.txt'
    # )
    # nu_CO, coef_CO = absorptionCoefficient_Voigt(
    #     SourceTables='CO', Diluent={'self': 1.0}, File='absCoefCO.txt'
    # )
    # nu_CH4, coef_CH4 = absorptionCoefficient_Voigt(
    #     SourceTables='CH4', Diluent={'self': 1.0}, File='absCoefCH4.txt'
    # )
    # nu_O2, coef_O2 = absorptionCoefficient_Voigt(
    #     SourceTables='O2', Diluent={'self': 1.0}, File='absCoefO2.txt'
    # )
    # nu_SO2, coef_SO2 = absorptionCoefficient_Voigt(
    #     SourceTables='SO2', Diluent={'self': 1.0}, File='absCoefSO2.txt'
    # )
    # nu_NH3, coef_NH3 = absorptionCoefficient_Voigt(
    #     SourceTables='NH3', Diluent={'self': 1.0}, File='absCoefNH3.txt'
    # )
    # nu_N2, coef_N2 = absorptionCoefficient_Voigt(
    #     SourceTables='N2', Diluent={'self': 1.0}, File='absCoefN2.txt'
    # )
    # nu_air, coef_air = absorptionCoefficient_Voigt(
    #     SourceTables=['H2O', 'CO2', 'O3', 'N2O', 'CO', 'CH4', 'O2', 'SO2', 'NH3', 'N2'],
    #     WavenumberRange=[0, 7500],
    #     WavenumberWingHW=50,
    #     Diluent={'H2O': H2O_conc*1e-2, 'CO2': CO2_conc*1e-2, 'O3': O3_conc*1e-2,
    #              'N2O': N2O_conc*1e-2, 'CO': CO_conc*1e-2, 'CH4': CH4_conc*1e-2,
    #              'O2': O2_conc*1e-2, 'SO2': SO2_conc*1e-2, 'NH3': NH3_conc*1e-2,
    #              'N2': N2_conc*1e-2},
    #     File='absCoefMixture.txt'
    # )
    # plot(nu_air, coef_air)
    # show()


if __name__ == '__main__':
    main()
