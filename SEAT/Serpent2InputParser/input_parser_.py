from collections.abc import Iterator
from dataclasses import dataclass, field

import SEAT.Serpent2InputWriter as siw


@dataclass(slots=True, frozen=False)
class Creator:
    """
    Prototype for the creators: creates the name
    """
    name: str = field(init=False)
    _done: bool = False
    _next_keyword: str = field(default_factory=str)

    def create(self) -> str:
        """
        Iteratively fills the fields Serpent2 object identities.
        If the end of the file is met, then the '__END__OF__FILE__' dummy keyword is used to finish the iteration in
        a clean way.

        :return: string - The next keyword. This is ment to be used in the match-case.
        """
        self.name = next(word_iterator)
        while not self._done:
            try:
                self.fill(next(word_iterator))
            except StopIteration:
                self.fill('__END__OF__FILE__')
        return self._next_keyword

    def fill(self, word: str):
        pass


@dataclass(slots=True, frozen=False)
class MaterialCreator(Creator):
    """
    Creator of Pin object instances.

    :param dens: string or None - the material density. Default is 'None'
    :param representation: string or None - the material representation. Default is 'None'
    :param tmp: string or None - the material temperature for doppler preprocessing. Default is 'None'
    :param tms: string or None - the material temperature for on-the-fly broadening. Default is 'None'
    :param tft: tuple or None - the material temperature limits for coupled multiphysics calculations. Default is 'None'
    :param rgb: tuple or None - the material color. Default is 'None'
    :param vol: string or None - the material volume. Default is 'None'
    :param mass: string or None - the material mass. Default is 'None'
    :param burn: string or None - the material subregions to burn separately. Default is 'None'
    :param fix: string or None - the library information for decay nuclides. Default is 'None'
    :param moder: string or None - the use of thermal scattering data for a nuclide. Default is 'None'
    """
    dens: str | None = None
    representation: str | None = None
    tmp: str | None = None
    tms: str | None = None
    tft: tuple[str, str] | None = None
    rgb: tuple[str, str, str] | None = None
    vol: str | None = None
    mass: str | None = None
    burn: str | None = None
    fix: tuple | None = None
    moder: list | None = None
    _moderated: bool = False

    def fill(self, word: str):
        """
        Fills the fields of the Material.

        :param word: string - the word to fill the object fields with.
        """
        global material_creators
        self.dens = next(word_iterator)
        # parse material options
        for w in word_iterator:
            if '.' in w:
                break
            else:
                match w:
                    case 'tmp':
                        self.tmp = next(word_iterator)
                    case 'tms':
                        self.tms = next(word_iterator)
                    case 'tft':
                        self.tft = (next(word_iterator), next(word_iterator))
                    case 'rgb':
                        self.rgb = (next(word_iterator), next(word_iterator), next(word_iterator))
                    case 'vol':
                        self.vol = next(word_iterator)
                    case 'mass':
                        self.mass = next(word_iterator)
                    case 'burn':
                        self.burn = next(word_iterator)
                    case 'fix':
                        self.fix = (next(word_iterator), next(word_iterator))
                    case 'moder':
                        self.moder = [] if self.moder is None else self.moder
                        self._moderated = True
                        done = False
                        while not done:
                            self.moder.append((next(word_iterator), next(word_iterator)))
                            done = True  # include here the possibility to have more moderators.
                            # One should stop iterating when a '.' is found in the string.
                    case _:
                        representation_elements = []
                        for ww in word_iterator:
                            if ww in KEYS:
                                self.representation = ''.join(representation_elements)
                                material_creators[self.name] = self
                                break
                            else:
                                representation_elements.append(f"{next(word_iterator)} {next(word_iterator)}\n")

    def initialize(self) -> siw.Material:
        tft = (float(self.tft[0]), float(self.tft[1])) if self.tft is not None else self.tft
        rgb = (int(self.rgb[0]), int(self.rgb[1]), int(self.rgb[2])) if self.rgb is not None else self.rgb
        fix = (self.fix[0], self.fix[1]) if self.fix is not None else self.fix
        moder = [(library, int(za)) for library, za in self.moder]
        return siw.Material(name=self.name, dens=float(self.dens),
                            representation=siw.MaterialRepresentation.from_string(self.representation),
                            tmp=float(self.tmp), tms=float(self.tms), tft=tft, rgb=rgb, vol=float(self.vol),
                            mass=float(self.mass), burn=int(self.burn), fix=fix, moder=moder)


@dataclass(slots=True, frozen=False)
class SurfaceCreator(Creator):
    """
    Creator of Surface object instances.

    """
    pass


@dataclass(slots=True, frozen=False)
class PinCreator(Creator):
    """
    Creator of Pin object instances.
    Separately creates materials and radi.

    :param materials: list - list of materials to use in the Pin initialisation. Default is '[]'
    :param radi: list - list of radi to use in the Pin initialisation. Default is '[]'
    """
    materials: list[str] = field(default_factory=list)
    radi: list[float] = field(default_factory=list)

    def fill(self, word: str):
        """
        Fills the fields of the Pin.
        Checks the passed string and appends it either to the pins or to the materials.
        If a new keyword is met, the iteration is stopped updating the self._done and the self._next_keyword parameters.

        :param word: string - the word to fill the object fields with.
        """
        global pin_creators, universe_creators
        if word in KEYS:  # When a keyword is met, the loop is stopped
            pin_creators[self.name] = self
            universe_creators[self.name] = self
            self._next_keyword = word
            self._done = True
        elif word.replace('.', '').isnumeric():  # radi will be numeric
            self.radi.append(float(word))
            self._done = False
        else:
            self.materials.append(word)  # else it is going to be a material
            self._done = False

    def initialize(self) -> siw.Pin:
        return siw.Pin(name=self.name, radi=[(materials[m], r) for m, r in zip(self.materials, self.radi)])


@dataclass(slots=True, frozen=False)
class CellCreator(Creator):
    """
    Creator of Cell object instances.
    Creates everything in one go in the fill method.

    :param father: string - the identity of the cell father universe. Default is "''"
    :param kind: string - the type of cell. It can either be 'fill', 'material' or 'outside'. Default is ''
    :param filler: string or None - the identity of the universe to fill the cell with. Default is 'None'
    :param material: string or None - the identity of the material to fill the cell with. Default is 'None'
    :param materials: list - list of identities of the surfaces delimiting the cell. Default is '[]'
    """
    father: str = field(default_factory=str)
    kind: str = field(default_factory=str)
    filler: str | None = None
    material: str | None = None
    delimiters: list[str] = field(default_factory=list)

    def fill(self, word: str):
        """
        Fills the fields of the Cell.
        Checks the passed string and initialises the proper parameters.
        Parses the 'kind' parameter from the string.
        When a keyword is met, the iteration is stopped.

        :param word: string - the word to fill the object fields with.
        """
        global cell_creators, universe_creators
        self.father = next(word_iterator)
        option = next(word_iterator)
        match option:
            case 'fill':  # the cell is filled with a universe
                self.kind = option
                self.filler = next(word_iterator)
            case 'outside':  # the cell is external
                self.kind = option
            case _:  # the cell is filled with a material
                self.kind = 'material'
                self.material = next(word_iterator)
        for w in word_iterator:
            if w in KEYS:  # if a keyword is found: conclude than break
                cell_creators[self.name] = self
                universe_creators[self.name] = self
                self._next_keyword = w
                self._done = True
                break
            self.delimiters.append(w)  # append delimiting sourfaces identities

    def initialize(self) -> siw.Cell:
        filler = universes[self.filler] if self.filler is not None else self.filler
        material = materials[self.material] if self.material is not None else self.material
        return siw.Cell(name=self.name, delimiters=[surfaces[s] for s in self.delimiters],
                        father=universes[self.father], kind=self.kind, filler=filler, material=material)


@dataclass(slots=True, frozen=False)
class LatticeCreator:
    pass


@dataclass(slots=True, frozen=False)
class DivisionCreator:
    pass


@dataclass(slots=True, frozen=False)
class NormalisationCreator:
    pass


@dataclass(slots=True, frozen=False)
class IntervalCreator:
    pass


def create_creator(word: str):
    option = None
    match word:
        case '/*':
            option = None
        case '%':
            option = None
        case 'mat':
            option = None
        case 'pin':
            option = PinCreator().create()
        case 'cell':
            option = None
        case 'lat':
            option = None
        case 'surf':
            option = None
        case 'dep':
            option = None
    if option is None:
        raise Exception(f"Keyword '{word}' does not meet any case.")
    return option


def initialize() -> siw.Simulation:
    pass


def parse(file: str) -> siw.Simulation:
    global word_iterator
    first: bool = True
    option: str = field(default_factory=str)
    with open(file, 'r') as f:
        words = f.read().split()
    word_iterator = iter(words)
    for word in word_iterator:
        option = create_creator(word) if first else create_creator(option)
        first = False
    return initialize()


def key_check(word: str):
    return word in KEYS


KEYS = ['mat', 'pin', 'cell', '__END__OF__FILE__', '%', '/*']

material_creators: dict[str, MaterialCreator] = {}
materials: dict[str, siw.Material] = {}

surface_creators: dict[str, SurfaceCreator] = {}
surfaces: dict[str, siw.Surface] = {}

cell_creators: dict[str, CellCreator] = {}
cells: dict[str, siw.Cell] = {}

pin_creators: dict[str, PinCreator] = {}
pins: dict[str, siw.Pin] = {}

universe_creators: dict[str, Creator] = {}
universes: dict[str, siw.NestedUniverse] = {}
# word_iterator initialised here and modified in the following to be the iterator over the file
word_iterator: Iterator[str] = field(default_factory=Iterator)

if __name__ == '__main__':
    filename = 'parser_test.txt'
    parse(filename)
