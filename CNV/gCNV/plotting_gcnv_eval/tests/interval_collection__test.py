##
# Testing interval collection module
#
from interval_collection import IntervalCollection
from interval import Interval

import argparse


def test_interval_collection_class(temp_dir: str):
    interval1_1001_2000 = Interval("1", 1001, 2000)
    interval1_2001_3000 = Interval("1", 2001, 3000)
    interval1_3001_4000 = Interval("1", 3001, 4000)
    interval1_10000_20000 = Interval("1", 10000, 20000)
    interval2_1001_2000 = Interval("2", 1001, 2000)
    interval2_2001_3000 = Interval("2", 2001, 3000)
    test_interval_list = [interval1_1001_2000, interval1_2001_3000, interval1_3001_4000, interval1_10000_20000,
                          interval2_1001_2000, interval2_2001_3000]
    test_interval_collection_header = "@TESTHEADER"
    test_interval_collection = IntervalCollection(
        interval_list=test_interval_list, header=test_interval_collection_header)

    # Test search method
    assert test_interval_collection.find_intersection(Interval("1", 500, 1500)) == [interval1_1001_2000]
    assert test_interval_collection.find_intersection(Interval("1", 2500, 3500)) == [interval1_2001_3000,
                                                                                     interval1_3001_4000]
    assert test_interval_collection.find_intersection(Interval("1", 100000, 100100)) == []

    # Test find_intersection_with_interval_and_truncate
    assert test_interval_collection.find_intersection_with_interval_and_truncate(Interval("1", 2500, 3500)) == [
        Interval("1", 2500, 3000), Interval("1", 3001, 3500)]
    # Test subtraction operator
    interval_collection_to_subtract = IntervalCollection([Interval("1", 1001, 2000), Interval("1", 2500, 2501), Interval("3", 100, 200)])
    assert test_interval_collection - interval_collection_to_subtract == IntervalCollection(
        interval_list=[interval1_3001_4000, interval1_10000_20000, interval2_1001_2000, interval2_2001_3000],
        header=test_interval_collection_header)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('temp_dir', metavar='TemporaryDirectory', type=str,
                        help='temporary directory for tests')
    args = parser.parse_args()
    temp_dir = args.temp_dir
    test_interval_collection_class(temp_dir)


if __name__ == '__main__':
    main()
