class Interval:
    """Stores a Genomic interval"""

    def __init__(self, chrom: str, start: int, end: int):
        self.chrom = chrom
        assert (start <= end), "Cannot create an interval where start %d is greater than end %d (on chromosome %s)"\
                               % (start, end, chrom)
        self.start = start
        self.end = end

    def length(self):
        return self.start - self.end + 1

    def intersects_with(self, interval):
        """
        Check whether two intervals overlap by at least one base

        Args:
            interval: interval to check for intersection with

        Returns: boolean

        """
        if self.chrom != interval.chrom:
            return False
        return ((self.start > interval.start) and (self.start < interval.end)) or \
               ((interval.start > self.start) and (interval.start < self.end))

    def __str__(self):
        return "Interval(" + self.chrom + ":" + str(self.start) + "-" + str(self.end) + ")"

    def __hash__(self):
        return (((hash(self.chrom) * 31) + hash(self.start)) * 31) + hash(self.end)

    def __eq__(self, other):
        return isinstance(other, Interval) and (self.chrom == other.chrom) \
               and (self.start == other.start) and (self.end == other.end)

    # This apparently is no longer necessary in Python 3
    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if self.chrom == other.chrom:
            return self.end <= other.start
        elif str.isdigit(self.chrom) and str.isdigit(other.chrom):
            return int(self.chrom) < int(other.chrom)
        elif str.isdigit(self.chrom) and not str.isdigit(other.chrom):
            return True
        elif not str.isdigit(self.chrom) and str.isdigit(other.chrom):
            return False
        else:
            assert (self.chrom == 'X' or self.chrom == 'Y') and (other.chrom == 'X' or other.chrom == 'Y')
            return self.chrom == 'X'
