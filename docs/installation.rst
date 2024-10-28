Installation Guide
==================

Follow the steps below to install AltAnalyze3 on your system.

**From PyPi**

.. code-block:: python

   pip install altanalyze3

**From GitHub**

Clone the GitHub repository:

.. code-block:: python

   git clone https://github.com/SalomonisLab/altanalyze3.git

Navigate to the project directory and install:

.. code-block:: python

   cd altanalyze3
   pip install -e .

Install dependencies:

.. code-block:: python

   pip install -r docs/requirements.txt

Verify the installation:

.. code-block:: python

   import altanalyze3
   print("AltAnalyze3 installed successfully!")
   
**On a Compute-Cluster**

.. code-block:: python

   proxy_on
   python -m venv /path_to_target/altanalyze3
   source /path_to_target/altanalyze3/bin/activate
   pip install setuptools
   cd /path_to_target/altanalyze3
   git clone https://github.com/SalomonisLab/altanalyze3.git
   module load gcc/12.2.0
   pip install -e .
   pip list
