Set up project dependent settings
=================================

All commands in scrimer scripts and manual suppose that you will set some environment 
variables that define your project and that you put the required tools into your path. 

.. note::

    The way of organizing your data presented here is just our suggestion. Python scripts 
    doing most of the work are not dependent on any directory structure.

Project directory
-----------------
To start a new project, create a new direcotry. Scrimer is working with ``.fastq`` 
data. Put your ``.fastq`` data in a subdirecotry called ``00-raw``. 

project.sh
----------
Create a file called ``project.sh`` in your project directory. It will consist of ``KEY=VALUE``
lines that will define your project specific settings, and each time you want to use scrimer
you'll start by::

    cd my/project/directory
    . project.sh

Example ``project.sh`` ::

  # location of genome data in your system
  # you need write access to add a new reference genome to that location
  GENOMES=/data/genomes
  
  # reference genome used
  GENOME=taeGut1
  GENOMEDIR=$GENOMES/$GENOME
  GENOMEFA=$GENOMEDIR/$GENOME.fa

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
For each scrimer session, you need to have all used executables in ``PATH``.
Using the file with paths according to your system you created during installation,
you can set up your ``PATH`` easily by::

    export PATH=$( cat ~/scrimer-env/paths | tr "\n" ":" ):$PATH
    
Such line can be at the end of your ``project.sh`` file, so everythig is set up at once.
