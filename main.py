import processing
import constants
import classes


def gather_data():
    raw_data = processing.parse(constants.DATA_PATH, 'model.txt')
    converted_data = processing.convert(raw_data)
    spectrum = classes.Spectrum("Model", converted_data)


if __name__ == '__main__':
    gather_data()
