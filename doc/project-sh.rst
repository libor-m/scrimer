Set up project dependent settings
=================================

All commands in scrimer scripts and manual suppose that you will set some environment 
variables that define your project and that you put the required tools into your path. 

NOTE: The way of organizing your data presented here is just our suggestion.

Project directory
-----------------
To start a new project, create a new direcotry. For now, scrimer is working with ``.fastq`` 
data. Put your ``.fastq`` data in a subdirecotry called ``00-raw``. 

Project.sh
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

  # primers used to synthetize cDNA
  # (sequences were found in .pdf report from the company that did the normalization)
  PRIMER1=AAGCAGTGGTATCAACGCAGAGTTTTTGTTTTTTTCTTTTTTTTTTNN  
  PRIMER2=AAGCAGTGGTATCAACGCAGAGTACGCGGG
  PRIMER3=AAGCAGTGGTATCAACGCAGAGT
  
Adding the tools to your PATH
-----------------------------
To add a tool to your ``PATH``, you need to know in which directory the executable resides.
Say ``tabix`` is located in ``/opt/tabix/bin``. Then you add ``tabix`` to your ``PATH``
by::

    export PATH=/opt/tabix/bin:$PATH
    
Such line, containing paths to all of the tools used by scrimer in your system can be at 
the end of your ``project.sh`` file, so everythig is set up at once when you source the file.
