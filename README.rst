.. image:: https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg
  :width: 200
  :alt: Contributor Covenant
  :target: https://www.contributor-covenant.org/version/2/0/code_of_conduct/


pwem
====

**pwem** is a Python module of Scipion framework for image processing in Electron Microscopy

The entire collection is licensed under the terms of the GNU Public License,
version 3 (GPLv3).

Development
-----------

To install **pwem** for development purposes, one can do:

.. code-block:: bash

    # Create a clean virtual environment
    conda create -n scipion python=3.8
    conda activate
    git clone https://github.com/scipion-em/scipion-em.git
    cd scipion-em
    pip install -e .

Running tests
-------------

.. code-block:: bash

    conda activate scipion
    cd scipion-em
    export SCIPION_DOMAIN="pwem"
    python -m unittest discover

API documentation
-----------------

https://scipion-em.github.io/docs/release-3.0.0/api/pwem/pwem.html
