from distutils.core import setup
import sys
import re

assert sys.version_info >= (3, 4), "germline_cnv_evaluation requires Python 3.4.x or later"

def get_version_string():
    version_file = "germline_cnv_evaluation/_version.py"
    version_str_line = open(version_file, "rt").read()
    version_regexp = r"^__version__ = ['\"]([^'\"]*)['\"]"
    re_out = re.search(version_regexp, version_str_line, re.M)
    if re_out is not None:
        return re_out.group(1)
    else:
        raise RuntimeError("Unable to find version string in %s." % (version_file,))

setup(
    name='germline_cnv_evaluation',
    author='Mehrtash Babadi',
    version=get_version_string(),
    author_email='mehrtash@broadinstitute.org',
    packages=['germline_cnv_evaluation'],
    description='Germline CNV Evaluation Suite',
    install_requires=[
        "gcnvkernel >= 0.5.1",
        "intervaltree_bio",
        "PyVCF",
        "pandas",
        "numpy"
    ])
