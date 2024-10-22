Installation Guide
==================

Follow the steps below to install AltAnalyze3 on your system.

1. Clone the GitHub repository:
   .. code-block:: bash

      git clone https://github.com/SalomonisLab/altanalyze3.git

2. Navigate to the project directory and install:
   .. code-block:: bash

      cd altanalyze3
      pip install -e .

3. Install dependencies:
   .. code-block:: bash

      pip install -r docs/requirements.txt

4. Verify the installation:
   .. code-block:: python

      import altanalyze3
      print("AltAnalyze3 installed successfully!")
