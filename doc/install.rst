Scrimer installation
====================

You need a default installation of **Python 2.7** [#Python]_ with **virtualenv** [#virtualenv]_.
If you're in an environment where you're not able to install virtualenv, we recommend to use 
the technique described at http://eli.thegreenplace.net/2013/04/20/bootstrapping-virtualenv/.

.. code-block:: bash
  
    # create and activate new python virtual environment for scrimer
    # in home directory of current user
    virtualenv ~/scrimer-env
    . ~/scrimer-env/bin/activate
    
    # install cython in advance because of pybedtools
    # and distribute because of pyvcf
    pip install cython distribute
    
    # now install scrimer from pypi
    # with it's additional dependencies (pyvcf, pysam, pybedtools)
    pip install scrimer
	
	# if you're in a grid environment, this can help
	virtualnev --relocatable ~/scrimer-env

Scrimer depends on several python modules, that should be installed automatically using the above procedue.

- **pysam** [#pysam]_ is used to manipulate with indexed fasta and bam files
- **pybedtools** [#pybedtools]_ is used to read and write the annotations
- **PyVCF** [#PyVCF]_ is used to access variants data

Apart from the python modules, the scrimer pipeline relies on some other tools that should be installed 
in your PATH. Default installations - following the included instructions - are sufficient.

- **bedtools** [#bedtools]_ is a dependency of pybedtools, used for manipulating with gff and bed files
- **samtools** [#samtools]_ is used for manipulating short read alignments, and for calling variants
- **LASTZ** [#lastz]_ is used to find the longest isotigs
- **tabix** [#tabix]_ creates compressed and indexed verisions of annotation files
- **GMAP** [#gmap]_ produces a spliced mapping of your contigs to the reference genome
- **smalt** [#smalt]_ maps short reads to consensus contigs to discover variants
- **GNU parallel** [#parallel]_ is used throughout the pipeline to speed up some lenghty calculations [#tange]_
- **blat** and **isPcr** [#blat]_ are used to check the designed primers

Additional tools can be installed to provide some more options.

- **FastQC** [#FastQC]_ can be used to check the tag cleaning process
- **agrep** [#agrep]_ and **tre-agrep** [#tre-agrep]_ can be used to check the tag cleaning process, 
  compile by ``make`` command
- **sort-alt** [#sortalt]_ provides alphanumeric sorting of chromosome names, 
  compile by ``make`` command
- **IGV** [#IGV]_ is great at visualizing the data when checking the results
- **newbler** [#newbler]_ is the best option for assembling 454 mRNA data [#mundry]_ [#kumar]_
- **MIRA** [#mira]_ does well on 454 transcriptome assembly as well [#mundry]_ [#kumar]_
- **sim4db** [#sim4db]_ can be used as alternative spliced mapper, 
  part of kmer suite, apply our patch [#sim4db-patch]_ to get standard conformant output
- **Pipe Viewer** [#pv]_ can be used to display progress of longer operations
- **BioPython** [#BioPython]_ and **NumPy** [#numpy]_ are required for running ``5prime_stats.py``
- **mawk** [#mawk]_, awk is often used in the pipeline, and mawk is usually an order of magnitude faster

.. _path:

Add installed tools to your PATH
--------------------------------
To easily manage locations of tools that you're using with Scrimer pipeline, create a text file
containing paths to directories, where binaries of your tools are located.
The format is one path per line, for example::

  /opt/bedtools/bin
  /opt/samtools-0.1.18
  /opt/lastz/bin
  /opt/tabix
  /opt/gmap/bin
  /data/samba/liborm/sw_testbed/smalt-0.7.4
  /data/samba/liborm/sw_testbed/FastQC
  /data/samba/liborm/sw_testbed/kmer/sim4db

Put this file to your virtual environment directory, e.g. ``~/scrimer-env/paths``.
You can run the following snippet when starting your work session:

.. code-block:: bash

  export PATH=$( cat ~/scrimer-env/paths | tr "\n" ":" ):$PATH

.. _software:

References
----------

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
.. [#tabix] tabix http://samtools.sourceforge.net/tabix.shtml, http://sourceforge.net/projects/samtools/files/
.. [#sortalt] sort-alt https://github.com/lh3/foreign/tree/master/sort
.. [#gmap] gmap http://research-pub.gene.com/gmap/
.. [#samtools] samtools http://samtools.sourceforge.net/, http://sourceforge.net/projects/samtools/files/
.. [#smalt] smalt http://www.sanger.ac.uk/resources/software/smalt/, 
   we used 0.7.0.1, because the latest version (0.7.3) crashes
.. [#parallel] GNU parallel http://www.gnu.org/software/parallel/
.. [#blat] http://users.soe.ucsc.edu/~kent/src/, get ``blatSrc35.zip`` and  ``isPcr33.zip``, 
   before ``make`` do ``export MACHTYPE`` and ``export BINDIR=<dir>``

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
.. [#sim4db-patch] patch for sim4db gff output http://sourceforge.net/p/kmer/patches/2/
.. [#pv] Pipe Viewer http://www.ivarch.com/programs/pv.shtml
.. [#mawk] mawk http://invisible-island.net/mawk/
.. [#yed] yEd http://www.yworks.com/en/products_yed_about.html

Papers
******
.. [#mundry] Mundry,M. et al. (2012) Evaluating Characteristics of De Novo Assembly Software on 454 Transcriptome Data: A Simulation Approach. PLoS ONE, 7, e31410. DOI: http://dx.doi.org/10.1371/journal.pone.0031410
.. [#kumar] Kumar,S. and Blaxter,M.L. (2010) Comparing de novo assemblers for 454 transcriptome data. BMC Genomics, 11, 571. DOI: http://dx.doi.org/10.1186/1471-2164-11-571
.. [#tange] Tange,O. (2011) GNU Parallel - The Command-Line Power Tool. ;login: The USENIX Magazine, 36, 42-47.
