import SEAT.nuclides

class test_Functions:
    def test_expand_za(self) -> None:
        assert SEAT.nuclides.expand_za(92235) == (92, 235, 0)
        assert SEAT.nuclides.expand_za(95641) == (95, 241, 1)
        assert SEAT.nuclides.expand_za(95241, method=None, meta=1) == (95, 241, 1)

    def test_expand_zam(self) -> None:
        assert SEAT.nuclides.expand_zam(922350) == (92, 235, 0)
        assert SEAT.nuclides.expand_zam(952411) == (95, 241, 1)

    def test_get_zam(self) -> None:
        assert SEAT.nuclides.get_zam(92, 235, 0) == 922350
        assert SEAT.nuclides.get_zam(95, 241, 1) == 952411 

    def test_get_za(self) -> None:
        assert SEAT.nuclides.get_za(92, 235, 0) == (92235, 0)
        assert SEAT.nuclides.get_za(95, 241, 1) == (95641, 1)
        assert SEAT.nuclides.get_za(95, 241, 0, method=None) == (95241, 0)

    def test_za2zam(self) -> None:
        assert SEAT.nuclides.za2zam(92235) == 922350
        assert SEAT.nuclides.za2zam(95641) == 952411

    def test_zam2za(self) -> None:
        assert SEAT.nuclides.zam2za(922350) == (92235, 0)
        assert SEAT.nuclides.zam2za(952411) == (95641, 1)

    def test_to_latex(self) -> None:
        assert SEAT.nuclides.to_latex(92, 235, 0) == '$^{235}$U'
        assert SEAT.nuclides.to_latex(95, 241, 1) == '$^{241\text{m}}$Am'

    def test_za2latex(self) -> None:
        assert SEAT.nuclides.za2latex(92235) == '$^{235}$U'
        assert SEAT.nuclides.za2latex(95641) == '$^{241\text{m}}$Am'

    def test_zam2latex(self) -> None:
        SEAT.nuclides.zam2latex(922350) == '$^{235}$U'
        assert SEAT.nuclides.za2latex(952411) == '$^{241\text{m}}$Am'

    def test_zam2nuclide(self) -> None:
        assert SEAT.nuclides.zam2nuclide(952411) == 'Am241m'
        assert SEAT.nuclides.zam2nuclide(952411, atomic_number=True) == '95Am241m'
        assert SEAT.nuclides.zam2nuclide(952411, sep='-') == 'Am-241m'

    def test_nuclide2zam(self) -> None:
        assert SEAT.nuclides.nuclide2zam('U235') == 922350
        assert SEAT.nuclides.nuclide2zam('Am241m') == 952411

    def test_nuclide2za(self) -> None:
        assert SEAT.nuclides.nuclide2za('U235') == (92235, 0)
        assert SEAT.nuclides.nuclide2za('Am241m') == (95641, 1)

    def test_get_meta_letter(self) -> None:
        assert SEAT.nuclides.get_meta_letter(0) == 'g'
        assert SEAT.nuclides.get_meta_letter(0, skip_ground=True) == ''

    def test_nuclide2element(self) -> None:
        assert SEAT.nuclides.nuclide2element('Am241m') == 'Am'