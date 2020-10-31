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

Getting Started
===============

1) Set up Bioconda

Follow [Bioconda Getting Started](https://bioconda.github.io/user/install.html) to enable Bioconda

2) Install AbNumber using Bioconda:

.. code-block:: bash

   conda install abnumber


Examples
========

>>> from abnumber import Chain
>>>
>>> seq = 'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSSAKTTAPSVYPLA'
>>> chain = Chain(seq, scheme='imgt')
>>>
>>> print(chain.format())
QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
                       ^^^^^^^^                 ^^^^^^^^                                      ^^^^^^^^^^^^
>>> chain.chain_type
'H'
>>> chain.seq
'QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS'
>>> chain.tail
'AKTTAPSVYPLA'
>>> chain.cdr3_seq
'ARYYDDHYCLDY'
>>> print(chain.tall_format())
fw1 H1    Q
fw1 H2    V
fw1 H3    Q
fw1 H4    L
fw1 H5    Q
fw1 H6    Q
fw1 H7    S
...
>>> seq2 = 'QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYDDYLDRWGQGTTLTVSSAKTTAP'
>>> chain2 = Chain(seq2, scheme='imgt')
>>>
>>> alignment = chain.align(chain2)
>>> print(alignment.format())
QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPS-RGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSS
||||.||||||.||||+|||||||||||.||||||||||||||||+||||||||.|.||||||||||||||||||||||||||.+|||||||||||||||||....||.|||||||||||
QVQLVQSGAELDRPGATVKMSCKASGYTTTRYTMHWVKQRPGQGLDWIGYINPSDRSYTNYNQKFKDKATLTTDKSSSTAYMQKTSLTSEDSAVYYCARYYD--DYLDRWGQGTTLTVSS
                        ^^^^^^^^                 ^^^^^^^^^                                      ^^^^^^^^^^^^


Chain
=====

.. autoclass:: Chain
   :members:

Alignment
=========

.. autoclass:: Alignment
   :members:

Position
========

.. autoclass:: Position   
   :members:
