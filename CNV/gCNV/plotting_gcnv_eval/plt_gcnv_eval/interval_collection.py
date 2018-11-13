from interval import Interval
import io_plt
import pandas as pd
import bisect

class IntervalCollection:
    """Collection of genomic intervals"""

    def __init__(self, interval_list: list, header: str=None):
        self.assert_interval_list_sorted(interval_list)
        self.header = header
        self.interval_list = interval_list
        self.last_searched_index = None

    @classmethod
    def read_interval_list(cls, interval_list_file):
        header = io_plt.read_comments_lines(interval_list_file)
        intervals_df = pd.read_csv(open(interval_list_file, 'r'), comment="@", delimiter="\t", usecols=[0,1,2], header=None, names=["CONTIG","START","END"])
        intervals_series = intervals_df.apply(lambda x: Interval(str(x.CONTIG), int(x.START), int(x.END)), axis=1)
        interval_list = intervals_series.tolist()
        return cls(interval_list=interval_list, header=header)

    @staticmethod
    def assert_interval_list_sorted(interval_list):
        for index in range(len(interval_list) - 1):
            assert (interval_list[index] < interval_list[index+1]), \
                     ("Interval list is not sorted for intervals Interval(%s) and Interval(%s)" % (str(interval_list[index]), str(interval_list[index+1])))

    def write_interval_list(self, output_file_path: str, output_file_name):
        with open(output_file_path + output_file_name, 'w') as output:
            if (self.header != None):
                output.write(self.header)
                output.write("\n")
            for interval in self.interval_list:
                output.write(interval.to_interval_file_string() + "\n")

    def find_intersecting_interval_indices(self, interval: Interval):
        if (self.last_searched_index == None):
            return self.__perform_binary_search(interval)
        else:
            query_result = self.__extend_left_right(interval, self.last_searched_index)
            if (not query_result):
                return self.__perform_binary_search(interval)
            else:
                return query_result

    def __perform_binary_search(self, interval: Interval):
        index = bisect.bisect_left(self.interval_list, interval)
        if (index < len(self.interval_list) and self.interval_list[index].intersects_with(interval)):
            last_searched_index = index
            return self.__extend_left_right(interval, index)
        else:
            return []

    def __extend_left_right(self, interval: Interval, intersecting_interval_index: int):
        indices = []
        current_index = intersecting_interval_index
        while (current_index >= 0 and self.interval_list[current_index].intersects_with(interval)):
            indices.insert(0, current_index)
            current_index -= 1
        current_index = intersecting_interval_index + 1
        while (current_index < len(self.interval_list) and self.interval_list[current_index].intersects_with(interval)):
            indices.append(current_index)
            current_index += 1
        return indices

    def __eq__(self, other):
        if (self.header != other.header):
            return False

        if (len(self.interval_list) != len(other.interval_list)):
            return False

        for index, interval in enumerate(self.interval_list):
            if (interval != other.interval_list[index]):
                return False
        return True

    # This apparently is no longer necessary in Python 3
    def __ne__(self, other):
        return not self.__eq__(other)