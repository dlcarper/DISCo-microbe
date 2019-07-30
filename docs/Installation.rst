==========================
DISCo-microbe Installation
==========================

DISCo-microbe is a python based package that can easily be installed with ``pip`` or directly from our github_ page.

pip
===

If you have MacOSX or you can use:

.. code-block:: bash

    >>> pip install disco-microbe


git
===

The newest development version of disco-microbe may be installed from our github_ page. You can use git directly, or download a zipfile.
    .. _github: https://github.com/dlcarper/Disco-microbe

.. code-block:: bash

    >>> git clone https://github.com/dlcarper/Disco-microbe.git
    >>> cd Disco-microbe
    >>> python setup.py install


anaconda
========
Anaconda_ provides the ``conda`` environment and is the recommended way to install DISCO-microbe if you are operating a Windows machine without a native ``python``. Once you have ``anaconda`` installed, create a ``conda`` environment and then use ``source`` and ``activate`` to launch a ``python`` environment.
    .. _Anaconda: https://www.anaconda.com

.. code-block:: bash

    >>> conda create -n disco-microbe-env python=3.7
    >>> source activate disco-microbe-env

You can now proceed with a pip install:

.. code-block:: bash

    >>> pip install disco-microbe

Test your installation
======================

You can test your installation by running the DISCo-microbe help command:

.. code-block:: bash

    >>> disco -h
