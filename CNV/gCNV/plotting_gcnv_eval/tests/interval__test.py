##
# Testing interval module
#
from interval import Interval


def test_interval_class():
    chrom = "1"
    start = 1000
    end = 2000
    test_interval = Interval(chrom=chrom, start=start, end=end)

    test_interval_correct_string = "1\t1000\t2000"

    assert (test_interval.__str__() == test_interval_correct_string)


def test_interval_lt():
    lesserIntervalChrom = Interval(chrom="1", start=500, end=600)
    greaterIntervalChrom = Interval(chrom="2", start=1000, end=2000)
    assert (greaterIntervalChrom > lesserIntervalChrom)

    lesserIntervalPos = Interval(chrom="1", start=1000, end=2000)
    greaterIntervalPos = Interval(chrom="1", start = 3000, end=4000)
    assert (greaterIntervalPos > lesserIntervalPos)

    overlapping_interval1 = Interval(chrom="1", start=1000, end=2000)
    overlapping_interval2 = Interval(chrom="1", start=1500, end=2500)
    assert not (overlapping_interval2 < overlapping_interval1)


def main():
    test_interval_class()
    test_interval_lt()


if __name__ == '__main__':
    main()