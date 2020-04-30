from collections import OrderedDict

from interval import Interval


class ReferenceDictionary:

    SQ_FIELD_START = "@SQ"
    SEQUENCE_NAME_TAG = "SN"
    SEQUENCE_LENGTH = "LN"
    FIELD_SEPARATOR = ":"

    def __init__(self, contigs: list, contig_to_interval_map: dict):
        self.contigs = contigs
        self.contig_to_interval_map = contig_to_interval_map

    @classmethod
    def read_in(cls, dict_file_path: str):
        contigs = [str(i + 1) for i in range(21)]
        contigs.extend(['X', 'Y'])
        contig_to_interval_map = OrderedDict()
        with open(dict_file_path, 'r') as dict_file:
            for line in dict_file:
                if line.startswith(ReferenceDictionary.SQ_FIELD_START):
                    fields = line.split(sep='\t')[1:]
                    attr_to_value_map = {}
                    for field in fields:
                        split_field = field.split(sep=ReferenceDictionary.FIELD_SEPARATOR)
                        attr = split_field[0]
                        value = split_field[1]
                        attr_to_value_map[attr] = value
                    contig = attr_to_value_map[ReferenceDictionary.SEQUENCE_NAME_TAG]
                    length = attr_to_value_map[ReferenceDictionary.SEQUENCE_LENGTH]
                    if contig in contigs:
                        contig_to_interval_map[contig] = Interval(contig, 1, int(length))
        return cls(contigs, contig_to_interval_map)

    def get_contig_interval_for_chrom_name(self, contig_name: str)->Interval:
        return self.contig_to_interval_map.get(contig_name)
