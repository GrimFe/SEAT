from SEAT.Serpent2InputWriter import Surface
from SEAT.Serpent2InputWriter import base

TEST_NAME = "Test"

parameters = [0, 0, 1]
s1 = Surface(name=TEST_NAME + '1', parameters=parameters)
s2 = Surface(name=TEST_NAME + '2', parameters=parameters)
s3 = Surface(name=TEST_NAME + '3', parameters=parameters)


class Test_Comment:
    comment = base.Comment(txt=TEST_NAME)

    def test_str(self) -> None:
        assert self.comment.__str__() == f'/* {TEST_NAME} */\n'

    def test_sum(self) -> None:
        comment2 = base.Comment(txt=TEST_NAME)
        assert (self.comment + comment2).__str__() == f'/* {TEST_NAME}\n   {TEST_NAME} */\n'


class Test_StandaloneComment:
    def test_str(self) -> None:
        comment = base.StandaloneComment(txt=TEST_NAME)
        assert comment.__str__() == f'/* {TEST_NAME} */\n\n'


class Test_InlineComment:
    def test_str(self) -> None:
        comment = base.InlineComment(txt=TEST_NAME)
        assert comment.__str__() == f'  % {TEST_NAME}\n'


class Test_Entity:
    def test_str(self) -> None:
        entity = base.Entity(name=TEST_NAME)
        assert entity.__str__() == TEST_NAME

    def test_numerical_name(self):
        entity = base.Entity(name=1)
        assert entity.__str__() == '1'


class Test_Functions:
    def test_reformat(self) -> None:
        assert base.reformat('a,b;c?d', ',;?') == 'abcd'

    def test_get_size(self) -> None:
        assert base.get_size([1, 2, [3, 4, 5]]) == 5

    def test_flatten(self) -> None:
        assert base.get_size([1, 2, 3, 4, 5]) == 5

    def test_intersect(self) -> None:
        intersect = base.intersect([s1, s2, s3])
        assert [item._operator for item in intersect] == ['', '', '']

    def test_unite(self) -> None:
        intersect = base.unite([s1, s2, s3])
        assert [item._operator for item in intersect] == ['', ' : ', ' : ']

    def test_surface_complement(self) -> None:
        intersect = base.surface_complement([s1, s2, s3])
        assert [item._operator for item in intersect] == ['-', '-', '-']


class Test_Other:
    def test_str(self) -> None:
        other = base.Other(TEST_NAME + ' string')
        assert other.__str__() == TEST_NAME + ' string'
