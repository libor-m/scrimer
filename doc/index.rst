Wellcome to Scrimer
===================

Scrimer is a GNU/Linux pipeline for designing PCR and genotyping primers from 454 transcriptomic data. 

Installation
------------
Scrimer itself is a convenience wrapper that uses several other tools to accomplish its goals. 
You need a default installation of Python 2.7 [#Python]_ with pysam [#pysam]_, pybedtools [#pybedtools]_ and PyVCF [#PyVCF]_. We recommend
to use virtualenv [#virtualenv]_ to use and manage different versions of Python packages in your Python installation. 
We provide the ``requirements.txt`` file to ease your installation. It's basically used as::
  
    # create and activate new python virtual environment for scrimer
    virtualenv ~/scrimer-env
    . ~/scrimer-env/activate
    
    # install required python modules
    pip install -r requirements.txt

Apart from the python modules, the scrimer pipeline relies on some other tools that should be installed 
in your PATH. Default installations - following the instructions supplied with package are sufficient.

- bedtools [#bedtools]_ is a dependency of pybedtools, used for manipulating with gff and bed files
- samtools [#samtools]_ is a dependency of pysam, used for manipulating with short read alignments, ad for calling variants
- LASTZ [#lastz]_ is used to find the longest isotigs
- tabix [#tabix]_ creates compressed and indexed verisions of annotation files
- GMAP [#gmap]_ produces a spliced mapping of your contigs to the reference genome
- smalt [#smalt]_ maps short reads to consensus contigs to discover variants
- GNU parallel [#parallel]_ is used throughout the pipeline to speed up some lenghty calculations

Additional tools can be installed to provide some more options.

- FastQC [#FastQC]_ can be used to check the tag cleaning process
- agrep [#agrep]_ and tre-agrep [#tre-agrep]_ can be used to check the tag cleaning process
- sort-alt [#sort-alt] provides alphanumeric sorting of chromosome names
- IGV [#IGV]_ is great at visualizing the data when checking the results
- newbler [#newbler]_ is the best option for assembling 454 mRNA data [#mundry]_
- MIRA [#mira]_ does well on 454 transcriptome assembly as well [#mundry]_
- sim4db [#sim4db]_ can be used as alternative spliced mapper
- Pipe Viewer [#pv]_ can be used to display progress of longer operations
- BioPython [#BioPython]_ and NumPy [#numpy]_ are required for running ``5prime_stats.py``

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

Python packages
***************

.. [#Python] Python http://www.python.org/
.. [#virtualenv] virtualenv http://www.virtualenv.org/en/latest/
.. [#pysam] pysam http://code.google.com/p/pysam/
.. [#pybedtools] pybedtools http://pythonhosted.org/pybedtools/
.. [#PyVCF] PyVCF https://github.com/jamescasbon/PyVCF

Other software
**************

.. [#lastz] lastz http://www.bx.psu.edu/~rsharris/lastz/
.. [#bedtools] bedtools http://code.google.com/p/bedtools/
.. [#tabix] tabix http://samtools.sourceforge.net/tabix.shtml
.. [#sort-alt] sort-alt https://github.com/lh3/foreign/tree/master/sort
.. [#gmap] gmap http://research-pub.gene.com/gmap/
.. [#samtools] samtools http://samtools.sourceforge.net/
.. [#smalt] smalt http://www.sanger.ac.uk/resources/software/smalt/
.. [#parallel] GNU parallel http://www.gnu.org/software/parallel/

Optional software
*****************

.. [#BioPython] BioPython http://biopython.org/
.. [#numpy] numpy http://www.numpy.org/

.. [#FastQC] FastQC http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. [#agrep] agrep https://github.com/Wikinaut/agrep
.. [#tre-agrep] tre-agrep http://laurikari.net/tre/
.. [#IGV] IGV http://www.broadinstitute.org/igv/
.. [#newbler] newbler http://454.com/products/analysis-software/index.asp
.. [#mira] MIRA http://www.chevreux.org/projects_mira.html
.. [#sim4db] sim4db http://sourceforge.net/apps/mediawiki/kmer/index.php?title=Main_Page
.. [#pv] Pipe Viewer http://www.ivarch.com/programs/pv.shtml

Papers
******
.. [#mundry] Mundry,M. et al. (2012) Evaluating Characteristics of De Novo Assembly Software on 454 Transcriptome Data: A Simulation Approach. PLoS ONE, 7, e31410.
