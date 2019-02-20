##
# Testing interval module
#
from interval import Interval


def test_interval_class():
    chrom = "1"
    start = 1000
    end = 2000
    test_interval = Interval(chrom=chrom, start=start, end=end)

    test_interval_correct_string = "Interval(1:1000-2000)"

    assert (test_interval.__str__() == test_interval_correct_string)


def test_interval_lt():
    lesser_interval_chrom = Interval(chrom="1", start=500, end=600)
    greater_interval_chrom = Interval(chrom="2", start=1000, end=2000)
    assert (greater_interval_chrom > lesser_interval_chrom)

    lesser_interval_pos = Interval(chrom="1", start=1000, end=2000)
    greater_interval_pos = Interval(chrom="1", start = 3000, end=4000)
    assert (greater_interval_pos > lesser_interval_pos)

    overlapping_interval1 = Interval(chrom="1", start=1000, end=2000)
    overlapping_interval2 = Interval(chrom="1", start=1500, end=2500)
    assert not (overlapping_interval2 < overlapping_interval1)


def main():
    test_interval_class()
    test_interval_lt()


if __name__ == '__main__':
    main()
