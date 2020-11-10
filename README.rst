====
pwem
====

**pwem** is a Python module of Scipion framework for image processing in Electron Microscopy


The entire collection is licensed under the terms of the GNU Public License,
version 3 (GPLv3).

-------------
Development
-------------

To install **pwem** for development purposes, one can do:

::

    # Create a clean virtual environment
    python -m venv ~/myenv
    source ~/myenv/bin/activate
    git clone git@github.com:scipion-em/scipion-em.git
    cd scipion-em
    python -m pip install -e .  # Install in the environment as development

-------------
Running tests
-------------

First make sure that **pwem** is available as a Python module in your
current Python environment. During development, I tend to set the PYTHONPATH:

::

    cd scipion-em
    # Either you have installed as mentioned above, or modify the PYTHONPATH
    export PYTHONPATH=$PYTHONPATH:$PWD
    # After pyworkflow is accessible as a module, then:
    cd pwem/tests

    python -m unittest discover

