#! /usr/bin/env python3
import os
from setuptools import setup, find_packages


def get_description():
    README = os.path.join(os.path.abspath(os.path.dirname(__file__)), "README.md")
    with open(README, "r") as f:
        return f.read()


setup(
    name="altanalyze3",
    description="AltAnalyze3",
    long_description=get_description(),
    long_description_content_type="text/markdown",
    version="0.0.1",
    url="https://github.com/SalomonisLab/altanalyze3.git",
    download_url="https://github.com/SalomonisLab/altanalyze3.git",
    author="",
    author_email="",
    license="Apache-2.0",
    include_package_data=True,
    packages=find_packages(
        exclude=["docs", "tests", "cwls"]
    ),
    install_requires=[
        "pysam"
    ],
    zip_safe=False,
    scripts=["altanalyze3/bin/altanalyze3"],
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Healthcare Industry",
        "License :: OSI Approved :: Apache Software License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Medical Science Apps."
    ]
)
