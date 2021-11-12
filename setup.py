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
        exclude=["doc", "test"]
    ),
    install_requires=[],
    zip_safe=False,
    classifiers=[]
)
