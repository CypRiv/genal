==============
The Geno class
==============

The main object of the package is the :class:`~genal.Geno` class that contains the SNP-level data and manipulates it through its methods.

.. autoclass:: genal.Geno

==============
Main functions
==============

Preprocessing
=============

The preprocessing of the SNP-level data is performed with the :func:`~genal.Geno.preprocess_data` method:

.. automethod:: genal.Geno.preprocess_data


Clumping
========

Clumping is performed with the :func:`~genal.Geno.clump` method:

.. automethod:: genal.Geno.clump

Polygenic Risk Scoring
======================

The computation of a polygenic risk score in a target population is performed with the :func:`~genal.Geno.prs` method:

.. automethod:: genal.Geno.prs

Querying outcome data
=====================

Before running Mendelian Randomization, the extraction of the genetic instruments from the :class:`~genal.Geno` object containing the SNP-outcome association data is done with :func:`~genal.Geno.query_outcome` method:

.. automethod:: genal.Geno.query_outcome

Mendelian Randomization
=======================

Various Mendelian Randomization methods are computed with the :func:`~genal.Geno.MR` method:

.. automethod:: genal.Geno.MR

MR-PRESSO
=========

The MR-PRESSO algorithm to detect and correct horizontal pleiotropy is executed with :func:`~genal.Geno.MRpresso` method:

.. automethod:: genal.Geno.MRpresso

Phenotype assignment
====================

Before running SNP-association tests, assigning a dataframe with phenotypic data to the :class:`~genal.Geno` object is done with :func:`~genal.Geno.set_phenotype` method:

.. automethod:: genal.Geno.set_phenotype

SNP-association tests
=====================

SNP-association testing is conducted with :func:`~genal.Geno.association_test` method:

.. automethod:: genal.Geno.association_test

Genetic lifting
===============

Lifting the SNP data to another genetic build is done with :func:`~genal.Geno.lift` method:

.. automethod:: genal.Geno.lift

GWAS Catalog
============

Querying the GWAS Catalog to extract traits associated with the SNPs is done with :func:`~genal.Geno.query_gwas_catalog` method:

.. automethod:: genal.Geno.query_gwas_catalog