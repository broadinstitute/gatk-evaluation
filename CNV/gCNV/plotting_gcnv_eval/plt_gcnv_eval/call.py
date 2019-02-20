import copy
from enum import Enum

from interval import Interval


class EventType(Enum):
    """Enumeration of possible alleles"""

    REF = 0
    DEL = 1
    DUP = 2
    NO_CALL = 3

    @classmethod
    def gcnv_genotype_to_event_type(cls, gcnv_call: int):
        return cls(gcnv_call)


class Call:
    """Stores an event type and call qualities for a single interval and single sample"""

    def __init__(self, interval: Interval, sample: str, event_type: EventType, call_attributes: dict):
        self.interval = interval
        self.sample = sample
        self.event_type = event_type
        self.call_attributes = call_attributes

    @classmethod
    def deep_copy(cls, call):
        interval = Interval(call.interval.chrom, call.interval.start, call.interval.end)
        return cls(interval=interval, sample=call.sample,
                   event_type=call.event_type, call_attributes=copy.deepcopy(call.call_attributes))

    def __eq__(self, other):
        return self.interval == other.interval and self.sample == other.sample \
               and self.event_type == other.event_type and self.call_attributes == other.call_attributes

    def __hash__(self) -> int:
        return super().__hash__()

    def __str__(self) -> str:
        return "Call(" + str(self.interval) + ", call_type: " + str(self.event_type) + ")"


