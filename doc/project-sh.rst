Set up project dependent settings
=================================

All commands in scrimer scripts and manual suppose that you will set some environment 
variables that define your project and that you add the required tools into your ``PATH``. 

Directory layout
----------------
.. note::

    The way of organizing your data presented here is just our suggestion. Python scripts 
    doing most of the work are not dependent on any directory structure.

Genomes directory
^^^^^^^^^^^^^^^^^
We suppose that genome data can be shared among different projects and different people
on the same machine. Thus we place it in location different from project specific data.

Project directory
^^^^^^^^^^^^^^^^^
A directory containing files specific for one input dataset. Various setps can be run with
various settings in the same project directory. We organize our files in a *waterfall*
structure of directories, where each directory name is prefixed with two digit number.
The directory name is some short meaningful description of the step, first digit in the 
prefix corresponds to part of the process (read mapping, variant calling etc.), and the 
second digit distinguishes substeps or runs with different settings.

To start a new project, create a new direcotry. Scrimer is working with ``.fastq`` 
data. Put your ``.fastq`` data in a subdirecotry called ``00-raw``. 

project.sh
----------
Create a file called ``project.sh`` in your project directory. It will consist of ``KEY=VALUE``
lines that will define your project specific settings, and each time you want to use scrimer
you'll start by:

.. code-block:: bash

    cd my/project/directory
    . project.sh

Example ``project.sh`` :

.. code-block:: bash

  # number of cores you want to use for parallel calculations
  CPUS=8

  # location of genome data in your system
  # you need write access to add a new reference genome to that location
  GENOMES=/data/genomes
  
  # reference genome used
  GENOME=taeGut1
  GENOMEDIR=$GENOMES/$GENOME
  GENOMEFA=$GENOMEDIR/$GENOME.fa
  
  # genome in blat format
  GENOME2BIT=$GENOMEDIR/$GENOME.2bit

  # gmap index location
  GMAP_IDX_DIR=$GENOMEDIR
  GMAP_IDX=gmap_${GENOME}
  
  # smalt index
  SMALT_IDX=$GENOMEDIR/smalt/${GENOME}k13s4

  # primers used to synthetize cDNA
  # (sequences were found in .pdf report from the company that did the normalization)
  PRIMER1=AAGCAGTGGTATCAACGCAGAGTTTTTGTTTTTTTCTTTTTTTTTTNN  
  PRIMER2=AAGCAGTGGTATCAACGCAGAGTACGCGGG
  PRIMER3=AAGCAGTGGTATCAACGCAGAGT
  
Adding the tools to your PATH
-----------------------------
For each scrimer session, all the executables that are used have to be in one of 
the directories mentioned in ``PATH``.
You can set up your ``PATH`` easily by using the file you created during :ref:`installation <path>`:

.. code-block:: bash

    export PATH=$( cat ~/scrimer-env/paths | tr "\n" ":" ):$PATH
    
Such line can be at the end of your ``project.sh`` file, so everythig is set up at once.

Alternatively you can copy all the tool executables into your virtual environment.
