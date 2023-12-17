
The Geno class
==============

The main object of the package is the :class:`~genal.Geno` class that contains the SNP-level data and manipulates it through its methods.

.. autoclass:: genal.Geno
   :members:
   
Clumping function
=================

Clumping is performed with the :func:`~genal.clump.clump_data` function:

.. autofunction:: genal.clump.clump_data

Extract and PRS functions
=========================

The SNP extraction from genetic files is done with :func:`~genal.extract_prs.extract_snps_func`:

.. autofunction:: genal.extract_prs.extract_snps_func

...