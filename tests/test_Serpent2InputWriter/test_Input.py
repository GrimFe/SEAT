from SEAT.Serpent2InputWriter import Simulation
from SEAT.Serpent2InputWriter import Comment, StandaloneComment
from SEAT.Serpent2InputWriter.base import reformat
from SEAT.Serpent2InputWriter.input import check_allowed_key


class test_Functions:
    def test_check_allowed_key(self) -> None:
        assert check_allowed_key({1: 'a'}, {1: 'a', 2: 'b'})
        assert not check_allowed_key({3: 'c'}, {1: 'a', 2: 'b'})

class test_DivisionWrapper:
    def test_format_sort(self) -> None:
        pass


class Test_DivisionWriter:
    def test_to_string(self) -> None:
        pass


class Test_Simulation:
    #                       Importing the simulation test elements
    from test_Composition import Test_Composition
    from test_Geometry import Test_Geometry
    from test_Depletion import Test_Depletion

    composition = Test_Composition.comp
    geometry = Test_Geometry.geom
    depletion = Test_Depletion.dep

    COMMENTS = {'Intro': StandaloneComment('Introductiory comment'),
                'Geometry': Comment('Geometry comment'),
                'Steps': Comment('Depletion comment'),
                'Materials': Comment('Composition comment'),
                'Others': Comment('Extra comment'),
                'Concluding': Comment('Concluding comment')}
    simulation = Simulation(geometry, depletion, composition, comments=COMMENTS)

    A, B = 1, 2

    def test_seed(self):
        assert self.simulation._seed == ''
        self.simulation.seed = 1
        assert self.simulation._seed == 'set seed 1\n'

    def test_dix(self):
        assert self.simulation._dix == 'set dix 0\n'
        self.simulation.dix = True
        assert self.simulation._dix == 'set dix 1\n'

    def test_fpcut(self):
        assert self.simulation._fpcut == 'set fpcut 0\n\n'
        self.simulation.fpcut = 1
        assert self.simulation._fpcut == 'set fpcut 1\n\n'

    def test_xscalc(self):
        assert self.simulation._xscalc == ''
        self.simulation.xscalc = 1
        assert self.simulation._xscalc == 'set xscalc 1'

    def test_printm_fraction(self):
        assert self.simulation._printm_fraction == 'set printm 0 1\n'
        self.simulation.printm_fraction = (True, 0)
        assert self.simulation._printm_fraction == 'set printm 1 0\n'

    def test_pcc(self) -> None:
        self.simulation.pcc = 'CE'
        self.simulation.pcc_param = (self.A, self.B)
        assert self.simulation._pcc == f'set pcc {0} {self.A} {self.B}\n'

    def test_bumode(self) -> None:
        self.simulation.bumode = 'CRAM'
        self.simulation.cram_param = (self.A, self.B)
        assert self.simulation._bumode == f'set bumode 2 {self.A} {self.B}\n'

    def test_gcu(self) -> None:
        from SEAT.Serpent2InputWriter import Universe
        universe = Universe('gcuTestUniverse')
        assert self.simulation._gcu == 'set gcu -1\n\n'
        self.simulation.gcu = [universe]
        assert self.simulation._gcu == f'set gcu {universe.name}\n\n'

    def test_population(self) -> None:
        assert self.simulation._population == f'set pop {100000} {250} {30} {1} {1} {1}\n'
        self.simulation.population['particles'] = 200000
        assert self.simulation._population == f'set pop {200000} {250} {30} {1} {1} {1}\n'

    def test_ures(self) -> None:
        assert self.simulation._ures == 'set ures 0\n\n'
        rp = self.composition.materials[0].representation
        self.simulation.ures = {"active": True, "representation": rp, "dilution_cut": 1e-9}
        assert self.simulation._ures == 'set ures 1 ' +\
                                        reformat(str(rp.get_nuclide_representation()), "[],'") +\
                                        ' ' + str(1e-9) + '\n\n'

    def test_bc(self) -> None:
        assert self.simulation._bc == 'set bc 1\n'
        self.simulation.bc = ('vacuum', 2, 'P')
        self.simulation.albedo = 0.5
        assert self.simulation._bc == 'set bc 1 2 3 0.5\n'

    def test_inventory(self) -> None:
        assert self.simulation._inventory == 'set inventory all\n'
        self.simulation.inventory = (3, 'ng')
        assert self.simulation._inventory == 'set inventory top 3 ng\n'
        self.simulation.inventory = ['U235', 'U238']
        assert self.simulation._inventory == 'set inventory 922350\n922380\n'

    def test_divisions(self) -> None:
        pass

    def test_str(self) -> None:
        pass
