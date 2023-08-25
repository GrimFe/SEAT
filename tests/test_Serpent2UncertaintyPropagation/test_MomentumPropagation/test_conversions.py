import uncertainties
from uncertainties import ufloat

from SEAT.Serpent2UncertaintyPropagation.MomentumPropagation.conversions import *

from functions import check_ufloat_equality

class Test_Functions:
    def test_AffineScalarFun_2_Variable(self):
        x = ufloat(1, 0.001)
        y = x * 2
        assert isinstance(x, uncertainties.core.Variable)
        assert isinstance(y, uncertainties.core.AffineScalarFunc)
        assert isinstance(AffineScalarFun_2_Variable(y),
                          uncertainties.core.Variable)
        assert check_ufloat_equality(AffineScalarFun_2_Variable(y),
                                     ufloat(2, 0.002))