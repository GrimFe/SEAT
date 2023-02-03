import pkg_resources

import SEAT._names
import SEAT.natural
import SEAT.composites
from SEAT.Serpent2InputWriter import composition

TEST_NAME = "Test"
COMPONENTS = {'He3': 2.46744e-10, 'He4': 0.000123372}
COMPOSITION_STR = "20030 2.46744e-10\n20040 0.000123372"

TEMPERATURE = 1100
MATERIAL_REPRESENTATION = composition.MaterialRepresentation(COMPONENTS, atomic=True, tmp=TEMPERATURE)


class Test_MaterialComposition:
    ## still not working
    # REPRESENTATION_FILE = pkg_resources.resource_filename(SEAT._names.PAKAGE_NAME,
    #                                                             'tests/test_Serpent2InputWriter/composition_file.txt')
    atomic = True
    enriched_uranium = composition.MaterialComposition({'U238': 0.9699411482790922, 'U235': 0.03, 'U234': 4.885172090774484e-05}, atomic=True)
    polluted_uranium = composition.MaterialComposition({'Pu239': 0.01, 'U238': 0.9828125999999999, 'U235': 0.007128, 'U234': 4.9500000000000004e-05}, atomic=True)
    UO2_enriched_str = '922380 0.30331472335539233\n922350 0.03\n922340 1.5276644607620947e-05\n80160 0.6650466666666666\n80170 0.00025333333333333333\n80180 0.0013666666666666666'
    H2O_polluted_str = '942390 0.01\n10010 0.6599241\n10020 7.59e-05\n80160 0.32919809999999994\n80170 0.00012539999999999999\n80180 0.0006765'

    def test_str(self) -> None:
        comp0 = composition.MaterialComposition(COMPONENTS, self.atomic)
        comp1 = composition.MaterialComposition({'HE3': 2.46744e-10, 'HE4': 0.000123372}, self.atomic)
        comp2 = composition.MaterialComposition({'he3': 2.46744e-10, 'he4': 0.000123372}, self.atomic)
        assert comp0.__str__() == comp1.__str__() == comp2.__str__() == COMPOSITION_STR

    def test_atomic(self) -> None:
        comp = composition.MaterialComposition(COMPONENTS, atomic=False)
        assert comp.__str__() == "20030 -2.46744e-10\n20040 -0.000123372"
        cmop1 = composition.MaterialComposition.from_zam({20030: 2.46744e-10, 20040: 0.000123372}, atomic=False)
        assert cmop1.__str__() == "20030 -2.46744e-10\n20040 -0.000123372"

    # def test_from_file(self) -> None:
    #     comp = composition.MaterialComposition.from_file(self.REPRESENTATION_FILE)
    #     assert comp.__str__() == COMPOSITION_STR

    def test_parse(self):
        assert composition.MaterialComposition.parse(COMPOSITION_STR, atomic=self.atomic).__str__() == COMPOSITION_STR

    def test_from_list(self):
        lst = ["20030 2.46744e-10", "20040 0.000123372"]
        assert composition.MaterialComposition._from_list(lst, atomic=self.atomic).__str__() == COMPOSITION_STR
        lst = ["20030***2.46744e-10", "20040***0.000123372"]
        assert composition.MaterialComposition._from_list(lst, separator='***', atomic=self.atomic).__str__() == COMPOSITION_STR

    def test_enriched_custom(self) -> None:
        # test enriched_custom(nuclides with nuclides)
        c = composition.MaterialComposition.enriched_custom(SEAT.natural.U, {'U235': 0.03}, self.atomic)
        assert c.__str__() == self.enriched_uranium.__str__()
        # test enriched_custom(elements with nuclides)
        e = composition.MaterialComposition.enriched_custom(SEAT.composites.UO2, {'U235': 0.03}, True)
        assert e.__str__() == self.UO2_enriched_str

    def test_enriched(self) -> None:
        # test enriched from natural abundances
        e = composition.MaterialComposition.enriched('U', {'U235': 0.03})
        assert e.__str__() == self.enriched_uranium.__str__()
        # test enriched from composites
        e = composition.MaterialComposition.enriched('UO2', {'U235': 0.03})
        assert e.__str__() == self.UO2_enriched_str

    def test_polluted_custom(self) -> None:
        # test polluted_custom(nuclides with nuclides)
        p = composition.MaterialComposition.polluted_custom(SEAT.natural.U, {'Pu239': 0.01}, self.atomic)
        assert p.__str__() == self.polluted_uranium.__str__()
        # test polluted_custom(elements with nuclides)
        p = composition.MaterialComposition.polluted_custom(SEAT.composites.H2O, {'Pu239': 0.01}, self.atomic)
        assert p.__str__() == self.H2O_polluted_str
        # test polluted_custom(nuclides with elements)
        p = composition.MaterialComposition.polluted_custom(SEAT.natural.Pb, {'Sb': 0.5}, self.atomic)
        assert p.__str__() == '511210 0.28605\n511230 0.21395\n822040 0.00689\n822060 0.11978\n822070 0.11037\n822080 0.26296'
        # test polluted_custom(elements with elements)
        p = composition.MaterialComposition.polluted_custom(SEAT.composites.H2O, {'Sb': 0.5}, self.atomic)
        assert p.__str__() == '511210 0.28605\n511230 0.21395\n10010 0.333295\n10020 3.8333333333333334e-05\n80160 0.16626166666666664\n80170 6.333333333333333e-05\n80180 0.00034166666666666666'

    def test_polluted(self) -> None:
        # test polluted from natural abundances
        p = composition.MaterialComposition.polluted('U', {'Pu239': 0.01})
        assert p.__str__() == self.polluted_uranium.__str__()
        # test polluted from composites
        p = composition.MaterialComposition.polluted('H2O', {'Pu239': 0.01})
        assert p.__str__() == self.H2O_polluted_str

    def test_natural(self) -> None:
        assert composition.MaterialComposition.natural('Be').__str__() == '40090 1.0'

    def test_predefined(self) -> None:
        assert composition.MaterialComposition.predefined('Bi').__str__() == '832090 1.0'

    def test_composite(self) -> None:
        c = composition.MaterialComposition.composite({'H': 0.5, 'B': 0.5}, self.atomic)
        assert c.__str__() == '10010 0.4999425\n10020 5.75e-05\n50100 0.0995\n50110 0.4005'

    def test_from_zam(self) -> None:
        comp = composition.MaterialComposition.from_zam({20030: 2.46744e-10, 20040: 0.000123372}, self.atomic)
        assert comp.__str__() == '20030 2.46744e-10\n20040 0.000123372'

    def test_from_za(self) -> None:
        comp = composition.MaterialComposition.from_za({2003: 2.46744e-10, 2004: 0.000123372}, self.atomic)
        assert comp.__str__() == '2003 2.46744e-10\n2004 0.000123372'

    def test_from_sep_nuclides(self) -> None:
        comp = composition.MaterialComposition.from_sep_nuclides({'He-3': 2.46744e-10, 'He-4': 1.23372e-4}, sep='-', atomic=self.atomic)
        assert comp.__str__() == COMPOSITION_STR

    def test_from_output_file(self) -> None:
        pass

    def test_adjust(self) -> None:
        c = composition.MaterialComposition.composite(SEAT.composites.H2O, atomic=self.atomic)
        assert c.components == {'H1': 0.66659, 'H2': 7.666666666666667e-05, 'O16': 0.3325233333333333, 'O17': 0.00012666666666666666, 'O18': 0.0006833333333333333}
        assert sum(c.components.values()) == 1
        assert (c.components['H1'] + c.components['H2']) / sum(c.components.values()) == 0.6666666666666666
        assert c.adjust('jendl40u', verbose=False).components == {'H1': 0.66659, 'H2': 7.666666666666667e-05, 'O16': 0.33333333333333326}
        assert sum(c.components.values()) == 0.9999999999999999
        assert (c.components['H1'] + c.components['H2']) / sum(c.components.values()) == 0.6666666666666667
        # test adjust when no adjustment is needed (e.g. adjusted.adjust))
        assert c.adjust('jendl40u', verbose=False) == c
        assert c.adjust('jendl40u', verbose=False) is c


class Test_MaterialRepresentation:
    REPRESENTATION_STRING = '2003.11c 2.46744e-10\n2004.11c 0.000123372\n'

    def test_str(self) -> None:
        representation = MATERIAL_REPRESENTATION
        assert representation.__str__() == self.REPRESENTATION_STRING
        representation1 = composition.MaterialRepresentation(COMPONENTS, atomic=True, tmp=TEMPERATURE, data_type='b')
        assert representation1.__str__() == self.REPRESENTATION_STRING.replace('c', 'b')
        # test zam builder
        r1 = composition.MaterialRepresentation.from_zam({20030: 2.46744e-10, 20040: 0.000123372}, atomic=True, tmp=TEMPERATURE)
        assert r1.__str__() == self.REPRESENTATION_STRING
        # test za builder
        r2 = composition.MaterialRepresentation.from_za({2003: 2.46744e-10, 2004: 0.000123372}, atomic=True, tmp=TEMPERATURE)
        assert r2.__str__() == self.REPRESENTATION_STRING

    def test_get_temperature(self) -> None:
        assert MATERIAL_REPRESENTATION.get_temperature() == '{:.0f}'.format(TEMPERATURE / 100)


class Test_Material:
    TEMPERATURE = 1200
    DENSITY = -10.4
    REPRESENTATION = composition.MaterialRepresentation(COMPONENTS, atomic=True, tmp=1200)
    BURN = 1
    RGB = '255 255 255'
    VOL = 1
    MASS = 1
    FIX = '09c 900'
    MODER = 'lwe70 2004'
    TARGET_STRING_ALL = f"mat {TEST_NAME} {DENSITY} tmp {TEMPERATURE} rgb {RGB} vol {VOL} mass {MASS} burn {BURN} fix {FIX} moder {MODER}\n" + \
                        f"2003.{TEMPERATURE / 100:.0f}c {list(COMPONENTS.items())[0][1]}\n" + \
                        f"2004.{TEMPERATURE / 100:.0f}c {list(COMPONENTS.items())[1][1]}\n\n"

    def _target_string(self, temperature_type: str, burn: int) -> str:
        if temperature_type != 'tft':
            string = f"mat {TEST_NAME} {self.DENSITY} {temperature_type} {self.TEMPERATURE} burn {burn}\n" + \
                     f"2003.{self.TEMPERATURE / 100:.0f}c {list(COMPONENTS.items())[0][1]}\n" + \
                     f"2004.{self.TEMPERATURE / 100:.0f}c {list(COMPONENTS.items())[1][1]}\n\n"
        else:
            string = f"mat {TEST_NAME} {self.DENSITY} {temperature_type} {self.TEMPERATURE} {self.TEMPERATURE} burn {burn}\n" + \
                     f"2003.{self.TEMPERATURE / 100:.0f}c {list(COMPONENTS.items())[0][1]}\n" + \
                     f"2004.{self.TEMPERATURE / 100:.0f}c {list(COMPONENTS.items())[1][1]}\n\n"
        return string

    def test_str_tmp(self) -> None:
        mat = composition.Material(name=TEST_NAME, dens=self.DENSITY, representation=self.REPRESENTATION, burn=True,
                                   tmp=self.TEMPERATURE)
        assert mat.__str__() == self._target_string('tmp', 1)

    def test_str_tms(self) -> None:
        mat = composition.Material(name=TEST_NAME, dens=self.DENSITY, representation=self.REPRESENTATION, burn=True,
                                   tms=self.TEMPERATURE)
        assert mat.__str__() == self._target_string('tms', 1)

    def test_str_tft(self) -> None:
        mat = composition.Material(name=TEST_NAME, dens=self.DENSITY, representation=self.REPRESENTATION, burn=True,
                                   tft=(self.TEMPERATURE, self.TEMPERATURE))
        assert mat.__str__() == self._target_string('tft', 1)

    def test_other_parameters(self) -> None:
        mat = composition.Material(name=TEST_NAME, dens=self.DENSITY,
                                   representation=self.REPRESENTATION,
                                   burn=True, tmp=self.TEMPERATURE,
                                   rgb=(255, 255, 255), vol=1, mass=1,
                                   fix={"09c": 900}, moder={'lwe70': 2004})
        assert mat.__str__() == self.TARGET_STRING_ALL

    def test_mix(self) -> None:
        mat1 = composition.Material(name='dummy1')
        c1 = 0.8
        mat2 = composition.Material(name='dummy2')
        c2 = 0.2
        mat = composition.Material.mix(name=TEST_NAME, materials=[(mat1, c1), (mat2, c2)])
        assert mat.__str__() == f"""mix {TEST_NAME}\n{mat1.name} {c1}\n{mat2.name} {c2}\n\n\n"""

    def test_get_temperature(self) -> None:
        mat = composition.Material(name=TEST_NAME, tms=self.TEMPERATURE)
        assert mat.get_temperature() == str(self.TEMPERATURE)
        mat = composition.Material(name=TEST_NAME, tmp=self.TEMPERATURE)
        assert mat.get_temperature() == str(self.TEMPERATURE)
        mat = composition.Material(name=TEST_NAME, tft=(self.TEMPERATURE, self.TEMPERATURE))
        assert mat.get_temperature() == f'{self.TEMPERATURE} {self.TEMPERATURE}'

    def test_get_temperature_kind(self) -> None:
        mat = composition.Material(name=TEST_NAME, dens=self.DENSITY, representation=self.REPRESENTATION, burn=True,
                                   tft=(self.TEMPERATURE, self.TEMPERATURE))
        assert mat.get_temperature_kind() == 'tft'

    def test_divide(self) -> None:
        pass


class Test_Division:
    VOLUME = 1

    def test_str(self) -> None:
        material = composition.Material(name=TEST_NAME)
        division = composition.Division(material, self.VOLUME)
        assert division.__str__() == f"{TEST_NAME} {self.VOLUME}"


class Test_Composition:
    TMP = 1200
    rep = MATERIAL_REPRESENTATION
    mat = composition.Material(name='TestMaterial4Composition', dens=1, representation=rep, burn=False, tmp=TMP)

    lib = "lib1"
    path = "/path/to/lib/1"

    comp = composition.Composition(
        materials=[mat],
        libraries={lib: path},
        scattering_name='lwe70',
        scattering_type='',
        scattering_tmp=-1,
        scattering_libs=['lwtr.16t'])

    STR = f"{lib} '{path}'\n\ntherm lwe70  lwtr.16t\n\n{mat}\n"

    def test_str(self) -> None:
        assert self.comp.__str__() == self.STR

    def test_get_subdivisions(self) -> None:
        pass
