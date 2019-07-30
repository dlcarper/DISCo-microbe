DISCo-microbe Tutorial
======================

This tutorial covers a DISCo workflow to generate a list of organisms that are identifiable at a user-specified edit distance (or minimum number of nucleotide differences) for a targeted region of DNA. The goal is to become familiarized with the general workflow, data files, and parameter settings in DISCo. We will use a referenced-based alignment of 16S rRNA gene sequences trimmed to the V4 region as an example, but the workflow applies to any other marker gene or organisms (e.g. fungal or other organism). Follow along by copy/pasting the
code-blocks into a command line terminal.

.. note::

    If you haven't installed DISCo yet, go here first:
    Installation_

    .. _Installation: https://github.com/dlcarper/DISCo-microbe/blob/master/docs/Installation.rst


Getting the data
----------------

Use the commands below to download and extract the data in a command-line interface. This will create a new directory called ``TUTORIAL_FILES/`` located in your current directory.

.. code-block:: bash

    >>> wget -L https://github.com/dlcarper/DISCo-microbe/raw/master/TUTORIAL_FILES.tar.gz
    >>> tar -xvzf TUTORIAL_FILES.tar.gz

Use the command ``cd`` to navigate and ``ls`` to look inside the ``TUTORIAL_FILES/`` directory. You will see that it contains different files to input in each module as well as examples of outputs from the program.

.. code-block:: bash

    >>> ls TUTORIAL_FILES/

.. parsed-literal::
    Tutorial_Metdata_file.txt
    Tutorial_alignment.fasta
    Tutorial_proportions_file.txt
    Tutorial_starter_community_file.txt

Create Module
-------------

The create module has two required arguments, an alignment of DNA or RNA sequences in FASTA format (``--i-alignment``) and a user-specified minimum sequence distance between community members (``--p-editdistance``). For this tutorial we will use an alignment of 16S rRNA genes that were subsampled from RDP, aligned with SINA, and columns with only gaps were removed. We use an edit distance of 3 so that our final community list will contain community members that have a minimum of 3 nucleotide differences. We will also include a seed for reproducibility purposes (``--p-seed``). To have taxonomic information in the final community list, we also input a tab-delimited metadata file (``--i-metadata``) that contains the sequence identifier followed by the taxonomic identification for each sequence.

.. code-block:: bash

    >>> disco create --i-alignment Tutorial_alignment.fasta --p-editdistance 3 --p-seed 10 --i-metadata Tutorial_Metdata_file.txt --o-community-list community_ED3_with_taxonomy.txt

This command should generate two files: a text file containing a list of members that differ by at least 3 nucleotides and a distance dictionary. The distance dictionary is a database of the pairwise sequence similarities for the provided community alignment. This distance database can be used as a starting point for the create module if you were to add community members to your existing alignment and only needed to calculate the pairwise distances for the new members.

Option to include specific strains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If you wish to generate a community that includes specific strains, you may add an additional argument (``--p-include-strains``) of a list of strains that must be included in the final community.

.. note::

    The list of strains must have an edit distance greater than or equal to the edit distance you provide for the final community. If your input list of strains have less nucleotide differences than the edit distance you specify, you will receive a conflict warning and must address the conflicts before proceeding.

.. code-block:: bash

        >>> disco create --i-alignment Tutorial_alignment.fasta --p-editdistance 3 --p-seed 10 --i-metadata Tutorial_Metdata_file.txt --o-community-list community_ED3_with_taxonomy_specific.txt  --p-include-strains Tutorial_starter_community_file.txt


Subsample Module
----------------

The subsample module is designed to take the final output community from the create module and provide a subsample of the community. The subsample module requires the input of the community generated from the create module. From here, the community can be subsampled to either include a specific number of strains (``--p-num-taxa``) or to represent specific proportions of a grouping variable (``--p-proportion``).

Subsample by total number of members to include
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To subsample by number of members to include, we need to provide the output community from the create module (``--i-input``) and the number of strains to include in the final community (``--p-num-taxa``). We will limit our community to 100 members and also include a seed number for reproducibility.

.. code-block:: bash

    >>> disco subsample --i-input-community community_ED3_with_taxonomy.txt --p-num-taxa 100 --p-seed 10

The above command should generate a tab delimited file that contains a list with only 100 community members that have a minimum of 3 nucleotide differences.

Subsample by proportions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To subsample by proportions of a grouping variable, we need to provide the output community from the create module (``--i-input``) and a file containing proportions of each group you wish to include (``--p-proportion``). We will subsample our community to reflect taxonomic proportions at the class level, of a natural microbiome and also include a seed number for reproducibility. We also need to indicate the column of the input community that we want to group by (here we use class).

.. code-block:: bash

    >>> disco subsample --i-input-community community_ED3_with_taxonomy.txt --p-proportion Tutorial_proportions_file.txt --p-seed 10 --p-group-by "Class"


Tutorial Completed
------------------
Congratulations! You have created a list of microbial strains that differ by at least 3 nucleotides. You then subsampled that list to either contain a specified number of strains or to reflect a specified proportion of groups. Please use the help option to view all options for the create and subsample modules.

.. code-block:: bash

    >>> disco create -h
    >>> disco subsample -h
