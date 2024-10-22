from setuptools import setup, find_packages
import os

def parse_requirements():
    with open('requirements.txt') as f:
        return f.read().splitlines()

setup(
    name="altanalyze3",
    version="0.0.1",
    description="AltAnalyze3",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="AltAnalyze Development Team",
    author_email="altanalyze@gmail.com",
    url="https://github.com/SalomonisLab/altanalyze3.git",
    license="Apache-2.0",
    include_package_data=True,
    packages=find_packages(exclude=["docs", "tests", "cwls", "tmp"]),
    install_requires=parse_requirements(),
    scripts=["altanalyze3/bin/altanalyze3"],
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    zip_safe=False,
)
