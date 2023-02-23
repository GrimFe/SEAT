from SEAT.Serpent2InputWriter import depletion
from SEAT.Serpent2InputWriter import composition

TEST_MATERIAL = composition.Material(name='TestMaterial')


class Test_Interval:
    interval = depletion.Interval(1)

    def test_step(self) -> None:
        interval = depletion.Interval(1, total=True)
        assert interval.step == 'tot'
        assert self.interval.step == 'step'

    def test_str(self) -> None:
        assert self.interval.__str__() == f'dep {self.interval.kind}{self.interval.step} {self.interval.span}\n'


class Test_Normalisation:
    POWER = 1

    def test_str(self) -> None:
        normalization = depletion.Normalization(value=self.POWER, material=TEST_MATERIAL, kind='power')
        assert normalization.__str__() == f'set {normalization.kind} {normalization.value} {normalization.material.name}\n'


class Test_Depletion:
    POWERS = [1, 2]
    TIMES = [10, 20]

    normalisations = [depletion.Normalization(power, TEST_MATERIAL) for power in POWERS]
    intervals = [depletion.Interval(t) for t in TIMES]

    dep = depletion.Depletion(list(zip(normalisations, intervals)))

    def test_str(self) -> None:
        assert self.dep.__str__() == f"{self.normalisations[0]}{self.intervals[0]}{self.normalisations[1]}{self.intervals[1]}\n"
