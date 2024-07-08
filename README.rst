.. image:: https://img.shields.io/pypi/v/scipion-em.svg
        :target: https://pypi.python.org/pypi/scipion-em
        :alt: PyPI release

.. image:: https://sonarcloud.io/api/project_badges/measure?project=scipion-em_scipion-em&metric=alert_status
        :alt: Quality Gate Status
        :target: https://sonarcloud.io/summary/new_code?id=scipion-em_scipion-em

.. image:: https://img.shields.io/sonar/quality_gate/scipion-em_scipion-em?server=https%3A%2F%2Fsonarcloud.io
        :target: https://sonarcloud.io/dashboard?id=scipion-em_scipion-em
        :alt: SonarCloud quality gate

.. image:: https://img.shields.io/pypi/l/scipion-em.svg
        :target: https://pypi.python.org/pypi/scipion-em
        :alt: License

.. image:: https://img.shields.io/pypi/pyversions/scipion-em.svg
        :target: https://pypi.python.org/pypi/scipion-em
        :alt: Supported Python versions

.. image:: https://img.shields.io/pypi/dm/scipion-em
        :target: https://pypi.python.org/pypi/scipion-em
        :alt: Downloads

.. image:: https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg
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
