##
# Testing interval collection module
#
from interval_collection import IntervalCollection
from interval import Interval

import argparse

def test_interval_collection_class(temp_dir: str):
    test_interval_list = [Interval("1", 1001, 2000), Interval("1", 2001, 3000),
                          Interval("1", 3001, 4000), Interval("1", 10000, 20000),
                          Interval("2", 1001, 2000), Interval("2", 2001, 3000)]
    test_interval_collection_header = "@TESTHEADER"
    test_interval_collection = IntervalCollection(
        interval_list=test_interval_list, header=test_interval_collection_header)
    test_interval_collection.write_interval_list(temp_dir, "test_intervals.tsv")
    parsed_interval_collection = IntervalCollection.read_interval_list(temp_dir + "test_intervals.tsv")
    assert parsed_interval_collection == test_interval_collection

    #Test search
    assert(test_interval_collection.find_intersecting_interval_indices(Interval("1", 500, 1500)) == [0])

    assert(test_interval_collection.find_intersecting_interval_indices(Interval("1", 2500, 3500)) == [1, 2])

    assert(test_interval_collection.find_intersecting_interval_indices(Interval("1", 100000, 100100)) == [])
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('temp_dir', metavar='TemporaryDirectory', type=str,
                    help='temporary directory for tests')
    args = parser.parse_args()
    temp_dir = args.temp_dir
    test_interval_collection_class(temp_dir)


if __name__ == '__main__':
    main()