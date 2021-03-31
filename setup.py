import setuptools
import codecs
import os.path

# Suggestion #1 from
# https://packaging.python.org/guides/single-sourcing-package-version/
# for handling versioning
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

DISTNAME = 'skipguide'
VERSION = get_version(DISTNAME + '/__init__.py')
DESCRIPTION = 'Prediction of CRISPR-Cas9 mediated Exon Skipping'

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name=DISTNAME,
    version=VERSION,
    author="Wilson Louie",
    author_email="wilsonlouie1@gmail.com",
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gifford-lab/skipguide",
    packages=setuptools.find_packages(),
    classifiers=[
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
    ],
    license='MIT',
    python_requires=">=3.5",
    include_package_data=True,
    install_requires=[
        'cyvcf2',
        'pyranges<=0.0.67',
        'mmsplice==1.0.3',
        'spliceai'
    ]
)