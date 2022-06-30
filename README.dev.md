# Installation

1. Install `virtualenv`
   ```
   pip3 install virtualenv
   ```
2. Create and activate virtual environment
   ```
   mkdir venv
   cd venv
   virtualenv .
   source ./bin/activate
   ```
3. Clone current repository and install AltAnalyze3 from there
   ```
   git clone https://github.com/SalomonisLab/altanalyze3.git
   cd altanalyze3
   pip3 install -e .
   ```
4. Test installation
   ```
   altanalyze3 --help
   ```

Note, for proper AltAnalyze functionality, you need to have [samtools](https://formulae.brew.sh/formula/samtools) intalled.
