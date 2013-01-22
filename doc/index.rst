Scrimer
=======

Scrimer is a pipeline for designing PCR and genotyping primers from 454 transcriptome data.

Pipeline workflow
-----------------

.. toctree::
   :maxdepth: 2
   
   prepare-reference
   remove-adaptors
   reference-assembly
   map-contigs
   map-reads
   choose-variants
   design-primers

Component documentation
-----------------------

.. toctree::
   :maxdepth: 3

   components

Dataflow
--------

Dataflow diagram of the pipeline. Inputs are in green, processing steps in yellow and results in red.

.. figure:: images/dataflow.png
   
Software used
-------------

.. store all software references here, and cite them throughout the documents

.. [BioPython] http://biopython.org/
.. [lastz] http://www.bx.psu.edu/~rsharris/lastz/
