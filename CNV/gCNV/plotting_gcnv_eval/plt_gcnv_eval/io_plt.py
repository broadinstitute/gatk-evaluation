# NEVER name your file io.py in the root module directory!!!
import vcf
import time
import datetime


def read_comments_lines(file_path):
    comment_prefixes = ["@", "#"]
    comment_string = ""
    with open(file_path, 'r') as file:
        for line in file.readlines():
            if not line[0] in comment_prefixes:
                # assume comment lines are contiguous
                break
            comment_string += line
    return comment_string.strip("\n")


def read_sample_names_from_vcfs(gcnv_vcfs: list):
    samples_names = []
    for vcf_file in gcnv_vcfs:
        vcf_reader = vcf.Reader(open(vcf_file, 'r'))
        assert len(vcf_reader.samples) == 1
        sample_name = vcf_reader.samples[0]
        samples_names.append(sample_name)
    return samples_names


def log(message: str):
    ts = time.time()
    st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
    print(st + " " + message)
