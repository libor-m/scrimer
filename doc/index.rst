Wellcome to Scrimer
===================

Scrimer is a GNU/Linux pipeline for designing PCR and genotyping primers from 454 transcriptomic data. 

Installation
------------
Scrimer itself is a convenience wrapper that uses several other tools to accomplish its goals. 
You need a default installation of [Python]_ with [pysam]_, [pybedtools]_ and [PyVCF]_. We recommend
to use [virtualenv]_ to use and manage different versions of Python packages in your Python installation. 
We provide the ? `frozen env file` to ease your installation. It's basically used as:
  
    some sample
    code lines that install virtualenv frozen env

Refer to `virtualenv refrence manual`_ for further information. 

``5prime_stats.py`` that is not needed in basic use of the pipeline requires [BioPython]_ and [numpy]_.)

Apart from the python modules, the scrimer pipeline relies on some other tools that should be installed 
in your PATH. Default installations - following the instructions supplied with package are sufficient.

- [bedtools]_ is a dependency of pybedtools, for manipulating with gff and bed files
- [samtools]_ is a dependency of pysam, for manipulating with short read alignments, calling variants
- [tabix]_ creates compressed and indexed verisions of annotation files
- [gmap]_ produces a spliced mapping of your contigs to the reference genome
- [smalt]_ maps short reads to consensus contigs to discover variants
- [GNU parallel]_ is used throughout the pipeline to speed up some lenghty calculations

Additional tools can be installed to provide some more options.

- [FastQC]_ can be used to check the tag cleaning process
- [agrep] and [tre-agrep] can be used to check the tag cleaning process
- [IGV] is great at visualizing the data when checking the results
- [newbler] is the best option for assembling 454 mRNA data
- [sim4db]_ can be used as alternative spliced mapper
- [Pipe view]_ can be used to display progress of longer operations

Pipeline workflow
-----------------
The workflow is divided into several steps. Intermediate results of these steps can be 
inspected and checked for validity. It is usually not meaningful to continue with the process 
if something went wrong, so the whole pipeline is a series of pre-made shell commands, 
that are supposed to be executed one by one by pasting them into the console.

Following pages describe these steps in detail, provide some explanations for the choices
we made and suggest ways how to check the validity of the intermediate results.

.. toctree::
   :maxdepth: 1
   
   project-sh
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
   :maxdepth: 2

   components

Dataflow
--------

Dataflow diagram of the pipeline.

.. figure:: images/dataflow.png
   
   Inputs are in green, processing steps in yellow and results in red.

License
-------
Scrimer is licensed under `GNU Affero General Public License <http://www.gnu.org/licenses/agpl.html>`_. 
Contact the author if you're interested in other licensing terms.

Author
------
| Libor Morkovsky
| Department of Zoology
| Charles University in Prague
| Czech Republic
| morkovsk[at]natur.cuni.cz

Software used
-------------

.. store all software references here, and cite them throughout the documents

.. [BioPython] http://biopython.org/
.. [Python] http://www.python.org/
.. [lastz] http://www.bx.psu.edu/~rsharris/lastz/
.. [bedtools] http://code.google.com/p/bedtools/
.. [tabix] http://samtools.sourceforge.net/tabix.shtml
.. [sort-alt] https://github.com/lh3/foreign/tree/master/sort
.. [gmap] http://research-pub.gene.com/gmap/
.. [samtools] http://samtools.sourceforge.net/
.. [smalt] http://www.sanger.ac.uk/resources/software/smalt/
.. [fastqc] http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. [pyvcf] https://github.com/jamescasbon/PyVCF
