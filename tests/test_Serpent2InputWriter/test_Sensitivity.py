import SEAT.Serpent2InputWriter.sensitivity as sensitivity

from SEAT.Serpent2InputWriter.detectors import Detector
from SEAT.Serpent2InputWriter.composition import Material


D1 = Detector("D1")
D2 = Detector("D2")

class Test_Detratio:
    def test_str(self):
        assert sensitivity.Detratio('dr', D1, D2).__str__() == "dr D1 D2"
        assert sensitivity.Detratio('dr', D1, D2, wd=2).__str__() == "dr D1 D2 1 2"


class Test_Response:
    BASE_RESPONSE = {'beff': True, 'leff': True}
    DETRATIOS = [sensitivity.Detratio('dr', D1, D2)]
    VOIDS = [Material("TEST_MATERIAL")]

    def test_base(self):
        target = {'keff': False,
                  'beff': True,
                  'leff': True,
                  'lambda': False}
        assert sensitivity.Response(base=self.BASE_RESPONSE)._base == target

    def test_str(self):
        target = 'sens resp beff\n'+\
                'sens resp leff\n'+\
                'sens resp detratiodr D1 D2\n'+\
                'sens resp void TEST_MATERIAL\n'
        assert sensitivity.Response(base=self.BASE_RESPONSE,
                                    detratios=self.DETRATIOS,
                                    voids=self.VOIDS).__str__() == target

    def test_lambda_total(self):
        assert sensitivity.Response(base={'lambda': True},
                                    lambda_total=True).__str__() == 'sens resp lambda 1 0\n'

    def test_lambda_groupwise(self):
        assert sensitivity.Response(base={'lambda': True},
                                    lambda_groupwise=True).__str__() == 'sens resp lambda 0 1\n'


class Test_BapasePerturbation:
    def test_xs(self):
        assert sensitivity._BasePerturbation(realist=['allmt'])._xs == ' allmt\n'
        assert sensitivity._BasePerturbation(realist=['all'])._xs == ' all\n'
        assert sensitivity._BasePerturbation(
            realist=['capt', 'fiss'])._xs == ' realist capt fiss\n'
        assert sensitivity._BasePerturbation(
            mtlist=[102, 18])._xs == ' mtlist 102 18\n'

    def test_zailist(self):
        target = [922350, 922380]
        assert sensitivity._BasePerturbation(
            zailist=[922350, 922380])._zailist == target
        assert sensitivity._BasePerturbation(
            zailist=['U235', 'U238'])._zailist == target
        assert sensitivity._BasePerturbation(
            zailist=[922350, 'U238'])._zailist == target


class Test_Perturbation:
    BASE_PERT = {'chi': True, 'elamu': True}

    def test_base(self):
        target = {'xs': False,
                  'chi': True,
                  'nubar': False,
                  'elamu': True,
                  'inlmu': False,
                  'eleg': False,
                  'temperature': False}
        assert sensitivity.Perturbation(base=self.BASE_PERT)._base == target

    def test_str(self):
        base_p = {'xs': True,
                  'chi': True,
                  'nubar': True,
                  'elamu': True,
                  'inlmu': True,
                  'eleg': True,
                  'temperature': True}
        nleg = 2
        p = sensitivity.Perturbation(base=base_p, nleg=nleg, mtlist=[102, 18],
                                     zailist=[922380, 922350],
                                     matlist=[Material("TEST_MATERIAL")])
        target = 'sens pert xs  mtlist 102 18\n'+\
            'sens pert chi\n'+\
            'sens pert nubar\n'+\
            'sens pert elamu\n'+\
            'sens pert inlmu\n'+\
            'sens pert eleg 2\n'+\
            'sens pert temperature\n'+\
            'sens pert zailist 922380 922350\n'+\
            'sens pert matlist TEST_MATERIAL\n'
        assert p.__str__() == target


class Test_CustomPerturbation:
    def test_str(self):
        cp = sensitivity.CustomPerturbation(name='TEST_NAME',
                                            efunc='TEST_FILE',
                                            realist=['fiss', 'capt'],
                                            zailist=[922350, 922380],
                                            matlist=[Material("TEST_MATERIAL")])
        target = 'sens pert custom TEST_NAME TEST_FILE zailist 922350 922380 '+\
            'matlist TEST_MATERIAL realist fiss capt\n'
        assert cp.__str__() == target


class Test_Sensitivity:
    def test_str(self):
        r = sensitivity.Response(base={'keff': True})
        p = sensitivity.Perturbation(realist=['capt', 'fiss'])
        target = '--- Responses ---\n'+\
            'sens resp keff\n\n'+\
            '--- Perturbations ---\n'+\
            'sens pert xs  realist capt fiss\n'+\
            'sens pert zailist all\n'+\
            'sens pert matlist all\n'
        assert sensitivity.Sensitivity(r, p).__str__() == target