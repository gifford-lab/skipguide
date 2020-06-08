import setuptools
# import skipguide


with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="skipguide",
    version="0.0.1",
    author="Wilson Louie",
    author_email="wilsonlouie1@gmail.com",
    description="Exon Skipping Prediction from CRISPR gRNA",
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
    # license='MIT',
    python_requires=">=3.5,<=3.7"
)