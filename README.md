# DISCO-microbe

Design of an Identifiable Synthetic Community of Microbes (DISCo-microbe) is an easy-to-use command-line program, for creation of diverse communities of organisms that can be distinguished through next-generation sequencing technology. DISCo-microbe consists of two modules, create and subsample. The create module constructs a highly diverse community at a specified sequence difference from an input of aligned DNA/RNA sequences, e.g., 16S sequence. The module can either design a de novo community or design a community that includes targeted organisms. create solves problem (1) by easily generating a diverse community of members through an easily documentable method, ensuring reproducibility. The subsample module provides options for dividing the community into subsets, according to either the number of members or the proportions of a grouping variable, both of which can be specified by the user. subsample module solves problem (2) by allowing the user to subsample an already distinguishable community of members based on attributes of interest. Although this software was designed for construction of microbial communities, any DNA/RNA alignment can be used as input; consequently, users are not restricted to any particular organismal group or marker gene. This program is implemented in Python and is available through GitHub and PYPI. 

## Installation


DISCo-microbe is a python based package that can easily be installed with ``pip`` or directly from our github_ page.

## pip

If you have a native installation of ``python`` you can use:

    >>> pip install disco-microbe


## git

The newest development version of disco-microbe may be installed from our github_ page. You can use git directly, or download a zipfile.

 https://github.com/dlcarper/Disco-microbe

    >>> git clone https://github.com/dlcarper/Disco-microbe.git
    >>> cd Disco-microbe
    >>> python setup.py install


## anaconda

Anaconda_ provides the ``conda`` environment and is the recommended way to install DISCO-microbe if you are operating a Windows machine without a native ``python``. Once you have ``anaconda`` installed, create a ``conda`` environment and then use ``source`` and ``activate`` to launch a ``python`` environment.

 https://www.anaconda.com


    >>> conda create -n disco-microbe-env python=3.7
    >>> source activate disco-microbe-env

You can now proceed with a pip install:


    >>> pip install disco-microbe

## Test your installation

You can test your installation by running the DISCo-microbe help command:


    >>> disco -h

## For full tutorial go to https://disco-design-of-an-identifiable-synthetic-community.readthedocs.io/en/latest/Tutorial.html
