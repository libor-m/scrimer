Wellcome to Scrimer
===================

Scrimer is a GNU/Linux pipeline for designing PCR and genotyping primers from 454 transcriptomic data. 

.. note::

    We present scrimer as an *end-to-end solution*, from raw reads to usable primers. This is intended to help
    novice users, so they can see the whole picture. But it has an important downside as well. The whole pipeline 
    is a series of steps that are considered to be complex on their own (contig assembly, read mapping, variant calling etc.).
    
    Our approach was to choose the most common solution for given step, use some reasonable default settings, and give the user 
    an option to choose another program - we're using standart formats for input and output. The fine tuning of each step is up
    to the users.

Installation
------------
Scrimer is a set of Python scripts that serve as a glue for several external programs. 
Bash commands used to run the external programs can be found in this documentation.

A detailed description how to install Scrimer can be found here:

.. toctree::
    :maxdepth: 1
    
    install

Pipeline workflow
-----------------
The workflow is divided into several steps. Intermediate results of these steps can be 
visually inspected and checked for validity. It is usually not meaningful to continue with the process 
if something went wrong. Whole pipeline is thus a series of pre-made shell commands, 
that are supposed to be executed one by one by pasting them into the console.

Following pages describe these steps in detail, provide some explanations for the choices
we made and suggest ways how to check validity of the intermediate results.

It takes about a day to push your data through the pipeline, if everything goes well.

.. toctree::
   :maxdepth: 1
   
   project-sh
   0-prepare-reference
   1-remove-adaptors
   2-reference-assembly
   3-map-contigs
   4-map-reads
   5-choose-variants
   6-design-primers
   igv

Scrimer dataflow
----------------

Dataflow diagram of the pipeline. Inputs are in green, processing steps in yellow and results in red. 
Edges are labelled with data format. Image was created with :ref:`yEd <yed>`.

.. figure:: images/dataflow.png
   

Scrimer components
------------------

.. toctree::
   :maxdepth: 2

   components

Source code and reporting bugs
------------------------------
Source code and a bugtracker can be found at https://github.com/libor-m/scrimer.

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
