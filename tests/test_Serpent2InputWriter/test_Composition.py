from SEAT.Serpent2InputWriter import composition

TEST_NAME = "Test"
COMPOSITION_DCT = {'He-3': 2.46744e-10, 'He-4': 1.23372e-4}
NUCLIDE_COMPOSITION = composition.MaterialComposition.from_nuclides(COMPOSITION_DCT)
COMPOSITION_STR = "He-3 2.46744e-10\nHe-4 0.000123372"


class Test_MaterialComposition:
    # REPRESENTATION_FILE_PATH = r'tests\test_Serpent2InputWriter\composition_file.txt'

    def test_str(self) -> None:
        comp = composition.MaterialComposition(components=COMPOSITION_STR)
        assert comp.__str__() == COMPOSITION_STR

    # def test_from_file(self) -> None:
    #    comp = composition.MaterialComposition.from_file(self.REPRESENTATION_FILE_PATH)
    #    assert comp.__str__() == COMPOSITION_STR

    def test_from_za(self) -> None:
        comp = composition.MaterialComposition.from_za({2003: 2.46744e-10, 2004: 0.000123372})
        assert comp.__str__() == '2003 2.46744e-10\n2004 0.000123372'

    def test_from_za_mass(self) -> None:
        comp = composition.MaterialComposition.from_za({2003: 2.46744e-10, 2004: 0.000123372}, atomic=False)
        assert comp.__str__() == '2003 -2.46744e-10\n2004 -0.000123372'

    def test_from_zam(self) -> None:
        comp = composition.MaterialComposition.from_za({20030: 2.46744e-10, 20040: 0.000123372})
        assert comp.__str__() == '20030 2.46744e-10\n20040 0.000123372'

    def test_from_nuclides(self) -> None:
        comp1 = composition.MaterialComposition.from_nuclides({'HE3': 2.46744e-10, 'HE4': 0.000123372})
        comp2 = composition.MaterialComposition.from_nuclides({'He3': 2.46744e-10, 'He4': 0.000123372})
        comp3 = composition.MaterialComposition.from_nuclides({'he3': 2.46744e-10, 'he4': 0.000123372})
        assert comp1.__str__() == comp2.__str__() == comp3.__str__() == '2003 2.46744e-10\n2004 0.000123372'
        assert NUCLIDE_COMPOSITION.__str__() == COMPOSITION_STR

    def test_from_output_file(self) -> None:
        pass


class Test_MaterialRepresentation:
    TEMPERATURE = 1100
    REPRESENTATION_STRING = 'He-3.11c 2.46744e-10\nHe-4.11c 0.000123372\n'

    def test_from_material_composition(self) -> None:
        representation = composition.MaterialRepresentation.from_material_composition(NUCLIDE_COMPOSITION,
                                                                                      tmp=self.TEMPERATURE)
        assert representation.__str__() == self.REPRESENTATION_STRING

    def test_from_string(self) -> None:
        representation = composition.MaterialRepresentation.from_string(COMPOSITION_STR, tmp=self.TEMPERATURE)
        assert representation.__str__() == self.REPRESENTATION_STRING

    def test_get_temperature(self) -> None:
        representation = composition.MaterialRepresentation.from_string(COMPOSITION_STR, tmp=self.TEMPERATURE)
        assert representation.get_temperature() == '{:.0f}'.format(self.TEMPERATURE / 100)


class Test_Material:
    TEMPERATURE = 1200
    DENSITY = -10.4
    REPRESENTATION = composition.MaterialRepresentation.from_material_composition(NUCLIDE_COMPOSITION, tmp=TEMPERATURE)
    BURN = 1
    RGB = '255 255 255'
    VOL = 1
    MASS = 1
    FIX = '1 1'
    MODER = 'lwe70 2004'
    TARGET_STRING_ALL = f"mat {TEST_NAME} {DENSITY} tmp {TEMPERATURE} rgb {RGB} vol {VOL} mass {MASS} burn {BURN} fix {FIX} moder {MODER}\n" + \
                        f"{list(COMPOSITION_DCT.items())[0][0]}.{TEMPERATURE / 100:.0f}c {list(COMPOSITION_DCT.items())[0][1]}\n" + \
                        f"{list(COMPOSITION_DCT.items())[1][0]}.{TEMPERATURE / 100:.0f}c {list(COMPOSITION_DCT.items())[1][1]}\n\n"

    def _target_string(self, temperature_type: str, burn: int) -> str:
        if temperature_type != 'tft':
            string = f"mat {TEST_NAME} {self.DENSITY} {temperature_type} {self.TEMPERATURE} burn {burn}\n" + \
                     f"{list(COMPOSITION_DCT.items())[0][0]}.{self.TEMPERATURE / 100:.0f}c {list(COMPOSITION_DCT.items())[0][1]}\n" + \
                     f"{list(COMPOSITION_DCT.items())[1][0]}.{self.TEMPERATURE / 100:.0f}c {list(COMPOSITION_DCT.items())[1][1]}\n\n"
        else:
            string = f"mat {TEST_NAME} {self.DENSITY} {temperature_type} {self.TEMPERATURE} {self.TEMPERATURE} burn {burn}\n" + \
                     f"{list(COMPOSITION_DCT.items())[0][0]}.{self.TEMPERATURE / 100:.0f}c {list(COMPOSITION_DCT.items())[0][1]}\n" + \
                     f"{list(COMPOSITION_DCT.items())[1][0]}.{self.TEMPERATURE / 100:.0f}c {list(COMPOSITION_DCT.items())[1][1]}\n\n"
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
        mat = composition.Material(name=TEST_NAME, dens=self.DENSITY, representation=self.REPRESENTATION, burn=True,
                                   tmp=self.TEMPERATURE,
                                   rgb=(255, 255, 255), vol=1, mass=1, fix=(1, 1), moder=[('lwe70', 2004)])
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
    mat_comp = composition.MaterialComposition.from_nuclides(COMPOSITION_DCT)
    rep = composition.MaterialRepresentation.from_material_composition(mat_comp, TMP)
    mat = composition.Material(name='TestMaterial4Composition', dens=1, representation=rep, burn=False, tmp=TMP)

    lib = "lib1"
    path = "/path/to/lib/1"

    comp = composition.Composition(
        materials=[mat],
        libraries={lib: path},
        scattering={('', 'lwe70'): (-1, ['lwtr.16t'])})

    STR = f"{lib} '{path}'\n\ntherm lwe70  lwtr.16t\n\n{mat}\n"

    def test_str(self) -> None:
        assert self.comp.__str__() == self.STR

    def test_get_subdivisions(self) -> None:
        pass
