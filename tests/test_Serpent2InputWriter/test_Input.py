import SEAT.Serpent2InputWriter.input as inpt

from SEAT.Serpent2InputWriter import Comment, StandaloneComment
from SEAT.Serpent2InputWriter.base import reformat
from SEAT.Serpent2InputWriter.input import check_allowed_key

from test_Composition import Test_Composition
from test_Geometry import Test_Geometry
from test_Depletion import Test_Depletion



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


class Test_GeometryPlot:
    xpix = 10
    ypix = 100
    kind = 'yz'
    pos = 5
    min1 = -2
    min2 = -3
    max1 = 2
    max2 = 3
    plt = inpt.GeometryPlot(xpix, ypix, kind, pos, min1, min2, max1, max2)

    def test_importance(self):
        assert self.plt._importance is None

    def test_str(self):
        test = f'1 {self.xpix} {self.ypix} {self.pos} {self.min1} {self.max1}' +\
            f' {self.min2} {self.max2}\n'
        assert self.plt.__str__() == test


class Test_ImportancePlot:
    xpix = 10
    ypix = 100
    kind = 'yz'
    fmin = 20
    fmax = 30
    e = 300
    plt = inpt.ImportancePlot(xpix, ypix, kind, fmin=fmin, fmax=fmax, energy=e)

    def test_importance(self):
        assert self.plt._importance == f'{self.fmin} {self.fmax} {self.e} '


class Test_MeshPlot:
    pass


class Test_Simulation:
    composition = Test_Composition.comp
    geometry = Test_Geometry.geom

    COMMENTS = {'Intro': StandaloneComment('Introductiory comment'),
                'Geometry': Comment('Geometry comment'),
                'Depletion': Comment('Depletion comment'),
                'Materials': Comment('Composition comment'),
                'Others': Comment('Extra comment'),
                'Concluding': Comment('Concluding comment')}
    simulation = inpt.Simulation(geometry, composition, comments=COMMENTS)

    def test_seed(self):
        assert self.simulation._seed == ''
        self.simulation.seed = 1
        assert self.simulation._seed == 'set seed 1\n'

    def test_dix(self):
        assert self.simulation._dix == 'set dix 0\n'
        self.simulation.dix = True
        assert self.simulation._dix == 'set dix 1\n'

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

    def test_depletion(self) -> None:
        assert self.simulation._depletion is None

    def test_str(self) -> None:
        pass


class Test_DepletionSimulation:
    composition = Test_Composition.comp
    geometry = Test_Geometry.geom
    depletion = Test_Depletion.dep

    COMMENTS = {'Intro': StandaloneComment('Introductiory comment'),
                'Geometry': Comment('Geometry comment'),
                'Depletion': Comment('Depletion comment'),
                'Materials': Comment('Composition comment'),
                'Others': Comment('Extra comment'),
                'Concluding': Comment('Concluding comment')}
    simulation = inpt.DepletionSimulation(geometry, composition, comments=COMMENTS,
                                          depletion=depletion)

    A, B = 1, 2

    fpcut = 1
    xscalc = 1
    printm_fraction = (True, 0)
    pcc = 'CE'
    pcc_param = (A, B)
    bumode = 'CRAM'
    cram_param = (A, B)
    inventory = (3, 'ng')

    def test_fpcut(self):
        assert self.simulation._fpcut == 'set fpcut 0\n'
        self.simulation.fpcut = 1
        assert self.simulation._fpcut == 'set fpcut 1\n'

    def test_xscalc(self):
        assert self.simulation._xscalc == ''
        self.simulation.xscalc = 1
        assert self.simulation._xscalc == 'set xscalc 1\n'

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

    def test_inventory(self) -> None:
        assert self.simulation._inventory == 'set inventory all\n'
        self.simulation.inventory = (3, 'ng')
        assert self.simulation._inventory == 'set inventory top 3 ng\n'
        self.simulation.inventory = ['U235', 'U238']
        assert self.simulation._inventory == 'set inventory 922350\n922380\n'

    def test_divisions(self) -> None:
        pass

    def test_depletion(self):
        self.simulation.fpcut = self.fpcut
        self.simulation.xscalc = self.xscalc
        self.simulation.printm_fraction = self.printm_fraction
        self.simulation.pcc = self.pcc
        self.simulation.pcc_param = self.pcc_param
        self.simulation.bumode = self.bumode
        self.simulation.cram_param = self.cram_param
        self.simulation.inventory = self.inventory

        test = self.simulation._inventory + self.simulation._fpcut +\
                self.simulation._bumode + self.simulation._pcc +\
                self.simulation._xscalc + self.simulation._printm_fraction +\
                "/* Depletion comment */\n" + self.depletion.__str__() +\
                self.simulation.divisions().to_string() + '\n'

        assert self.simulation._depletion == test
