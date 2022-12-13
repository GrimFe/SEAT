from collections.abc import Iterator
from dataclasses import dataclass

import SEAT.Serpent2InputWriter as siw

KEYS = ['mat', 'pin', 'cell', '__END_OF_FILE__']

materials = {}
surfaces = {}
cells = {}
pins = {}

word_iterator: Iterator['str'] = None  # initialised here and modified in the following to be the iterator over the file


@dataclass(slots=True, frozen=False)
class Creator:
    """
    Prototype for the creators: creates the name
    """
    word: str
    name: str

    def get_name(self):
        self.name = self.word


@dataclass(slots=True, frozen=False)
class MaterialCreator:
    data: str

    def from_line(self) -> siw.Material:
        data_ = self.data.split('\n')[0]
        comp = '\n'.join(self.data.split('\n')[1:])
        name = data_[1]
        dens = float(data_[2])
        # ... how to continue?

        return siw.Material(name=name, dens=dens, representation=siw.MaterialRepresentation(comp))


@dataclass(slots=True, frozen=False)
class SurfaceCreator:
    data: str

    def from_line(self) -> siw.Surface:
        return siw.Surface()


@dataclass(slots=True, frozen=False)
class PinCreator:
    data: str

    def from_line(self) -> siw.Pin:
        data_ = self.data.split()
        name = data_[1]
        radi_builder = zip(data_[2::2], data_[3::2])
        radi = [(materials[m].copy(), float(r)) for m, r in radi_builder]
        if radi[-1][1] != data_[-1]:  # i.e. if external material (is this optional?)
            radi.append((materials[data_[-1]], 0))
        return siw.Pin(name=name, radi=radi)

    def from_word(self):
        global word_iterator
        for word in word_iterator:
            self.update(word)

    def update(self, word):
        if key_check(word):
            self.radi_builder = zip(self.materials, self.radi)
            create_creator_from_word(word)
        elif word.replace('.', '').isnumeric():
            self.radi.append(float(word))
        else:
            self.materials.append(word)


@dataclass(slots=True, frozen=False)
class CellCreator:
    data: str

    def from_line(self) -> siw.Cell:
        data_ = self.data.split()
        name = data_[1]
        father = siw.NestedUniverse(name=data_[2])
        kind = data_[3]
        filler, material = None, None
        if kind == 'fill':
            filler = siw.Universe(name=data_[4])
        elif kind == 'material':
            meterial = materials[data_[4]].copy()
        elif kind == 'outside':
            pass
        else:
            raise Exception(f"Allowed cell types are: 'filler', 'material', 'outside': {data_[4]} was passed.")
        delimiters = []
        operator = ''
        for string in data_[5:]:
            if string.isalnum():
                surf_name = string
                s = surfaces[surf_name].copy()
                s._operator = operator + ' '  # the extra space might be not needed
                delimiters.append(s)
        else:
            operator += operator
        return siw.Cell(name=name, father=father, kind=kind, filler=filler, material=material, delimiters=delimiters)


@dataclass(slots=True, frozen=False)
class LatticeCreator:
    data: str

    def from_line(self) -> siw.Lattice:
        return siw.Lattice()


@dataclass(slots=True, frozen=False)
class DivisionCreator:
    pass


@dataclass(slots=True, frozen=False)
class NormalisationCreator:
    data: str

    def from_line(self) -> siw.Normalization:
        return siw.Normalization()


@dataclass(slots=True, frozen=False)
class IntervalCreator:
    data: str

    def from_line(self) -> siw.Interval:
        return siw.Interval()


def create_creator_from_word(word: str):
    print('------------------')
    match word:
        case '/*':
            pass
        case '%':
            pass
        case 'mat':
            pass
        case 'pin':
            PinCreator(word).from_word()
        case 'cell':
            pass
        case 'lat':
            pass
        case 'surf':
            pass
        case 'dep':
            pass


def create_creator_from_line(builder: str):
    print('------------------')
    if is_comment(builder):
        print('Skipping comment')
    elif builder.lstrip().startswith('mat'):
        # print(f"Creating material creator for {builder.lstrip().split(NEWLINE)[0].split(' ')[0:2]}")
        pass
    elif builder.lstrip().startswith('pin'):
        # print(f"Creating pin creator for {builder.lstrip().split(NEWLINE)[0].split(' ')[0:2]}")
        PinCreator(builder).from_line()
    elif builder.lstrip().startswith('cell'):
        # print(f"Creating cell creator for {builder.lstrip().split(NEWLINE)[0].split(' ')[0:2]}")
        pass
    elif builder.lstrip().startswith('lat'):
        # print(f"Creating lattice creator for {builder.lstrip().split(NEWLINE)[0].split(' ')[0:2]}")
        pass
    elif builder.lstrip().startswith('surf'):
        # print(f"Creating surface creator for {builder.lstrip().split(NEWLINE)[0].split(' ')[0:2]}")
        pass
    elif builder.lstrip().startswith('dep'):
        # print(f"Creating interval creator for {builder.lstrip().split(NEWLINE)[0].split(' ')[0:2]}")
        pass


def is_comment(st: str):
    return st.lstrip().startswith('/*') or st.strip().startswith('%')


def key_check(ln: str):
    begins_with_key = ln.lstrip().split(' ')[0] in KEYS
    return begins_with_key or is_comment(ln)


if __name__ == '__main__':
    string = """
    pin 1   pin1
            pin1
    
    cell 1  cell1
    
    mat 1   mat1
            mat1
            mat1
            mat1
    /*one
    thing
    */
    pin 2   pin2
            pin2
    
    cell 2  cell2
    pin 3   pin3
            pin3
    mat 2   mat2
            mat2
            mat2
            mat2
    cell 3  cell3
    % one other thing
    mat 3   mat3
            mat3
            mat3
            mat3
    __END_OF_FILE__
    """
    concat = ''
    for line in string.split('\n'):
        if line != '':
            if key_check(line) and concat != '':
                create_creator_from_line(concat)
                concat = line
            else:
                concat = concat + '\n' + line
