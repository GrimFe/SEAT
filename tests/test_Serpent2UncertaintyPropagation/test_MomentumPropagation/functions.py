from uncertainties import ufloat

def check_ufloat_equality(a: ufloat, b: ufloat) -> bool:
    test1 = a.nominal_value == b.nominal_value
    test2 = a.std_dev == b.std_dev
    return test1 and test2