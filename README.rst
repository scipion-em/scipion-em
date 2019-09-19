pwem
===========

Development
-------------

We are now fully going toward Python 3!

To install **pwem** for development purposes, one can do:

.. code-block:: bash

    # Create a clean virtual environment
    python -m venv ~/myenv
    source ~/myenv/bin/activate
    git clone git@github.com:scipion-em/scipion-em.git
    cd scipion-em
    python -m pip install -e .  # Install in the environment as development

Running tests
.............
First make sure that **pwem** is available as a Python module in your
current Python environment. During development, I tend to set the PYTHONPATH:

.. code-block:: bash

    cd scipion-em
    # Either you have installed as mentioned above, or modify the PYTHONPATH
    export PYTHONPATH=$PYTHONPATH:$PWD
    # After pyworkflow is accesible as a module, then:
    cd pwem/tests

    python -m unittest discover

