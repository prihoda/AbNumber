.. AbNumber documentation master file, created by
   sphinx-quickstart on Thu Oct 29 22:07:08 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AbNumber's documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. automodule:: abnumber
   :members:

Antibody chain representation, alignment and numbering using ANARCI

Getting Started
===============

Install AbNumber using Bioconda:

.. code-block:: bash

   conda install -c bioconda abnumber

Credits
-------

This tool is based on `ANARCI
<https://github.com/oxpig/ANARCI>`_, please cite the ANARCI paper:
`ANARCI: antigen receptor numbering and receptor classification
<https://doi.org/10.1093/bioinformatics/btv552>`_


Examples
========

See the `Example Jupyter Notebook
<https://github.com/prihoda/AbNumber/tree/master/examples/AbNumber_getting_started.ipynb>`_ for usage examples.

Chain
=====

.. autoclass:: Chain
   :members:

Alignment
=========

See the `Example Jupyter Notebook
<https://github.com/prihoda/AbNumber/tree/master/examples/AbNumber_getting_started.ipynb>`_ for usage examples.

.. autoclass:: Alignment
   :members:

Position
========

See the `Example Jupyter Notebook
<https://github.com/prihoda/AbNumber/tree/master/examples/AbNumber_getting_started.ipynb>`_ for usage examples.

.. autoclass:: Position   
   :members: