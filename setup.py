from setuptools import setup, find_packages
import os

def parse_requirements():
    with open('requirements.txt') as f:
        return f.read().splitlines()

setup(
    name="altanalyze3",
    version="0.1.3",
    description="AltAnalyze3",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="AltAnalyze Development Team",
    author_email="altanalyze@gmail.com",
    url="https://github.com/SalomonisLab/altanalyze3.git",
    license="Apache-2.0",
    include_package_data=True,
    package_data={
        "altanalyze3.components.visualization.isv_web": ["templates/*.html", "static/*.js", "static/*.css", "*.md"],
        "altanalyze3.components.cellHarmony.webapp": ["templates/*.html", "static/*.js", "static/*.css", "static/*.gif", "static/*.png", "*.md"],
        "altanalyze3.components.cellHarmony": ["resources/biology/*.json", "resources/biology/*.md", "resources/biology/*.py", "discover_resources/source_manifest.json"],
        "altanalyze3.components.model_registry": ["catalog.json"],
        "altanalyze3.components.rna2metabolite": ["resources/*.json", "missingness_release.json", "artifacts/rna2metabolite_aml_bundle.pkl.gz", "artifacts/rna2metabolite_aml_NA30_20261006_bundle.pkl.gz"],
        "altanalyze3.components.snaf": ["docs/*.md", "surface/data/*.json"],
        "altanalyze3.components.fastCNV": ["resources/*.tsv"],
        "altanalyze3.components.fastComm": ["resources/*.tsv", "resources/*/*.tsv", "resources/*/*.json", "configs/*.json"],
        "altanalyze3.components.long_read": ["resources/*.gz"],
    },
    packages=find_packages(exclude=["docs", "tests", "cwls", "tmp"]),
    install_requires=parse_requirements(),
    extras_require={"iso2function": ["igraph>=0.11,<1.1", "adjustText>=0.8,<2"]},
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
