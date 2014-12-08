Components of the Scrimer pipeline
##################################

Components bound to the pipeline
================================

Transfer annotation from the reference genome to the contigs
------------------------------------------------------------
File: ``scripts/liftover.py``

.. automodule:: liftover

Create a contig scaffold according to contig hits in the reference genome
-------------------------------------------------------------------------
File: ``scripts/scaffold.py``

.. automodule:: scaffold

Design primers using sequences, annotations and selected variants
-----------------------------------------------------------------
File: ``scripts/design_primers.py``

.. automodule:: design_primers

Export primers from gff3 to the FASTA, tabular or isPcr format
--------------------------------------------------------------
File: ``scripts/extract_primers.py``

.. automodule:: extract_primers


Tools for FASTA/FASTQ manipulation
==================================

Given pairs of matching sequences, create clusters and find the longest representative
--------------------------------------------------------------------------------------
File: ``scripts/find_redundant_sequences.py``

.. automodule:: find_redundant_sequences

Filter sequences in a FASTQ file based on their position in the file
--------------------------------------------------------------------
File: ``scripts/fastq_kill_lines.py``

.. automodule:: fastq_kill_lines

Filter sequences in FASTQ, FASTA based on their identifier
----------------------------------------------------------
File: ``scripts/seq_filter_by_id.py``

Taken from :ref:`BioPython <BioPython>`. FASTA and FASTQ readers are pasted in the file,
so the program is standalone.

.. automodule:: seq_filter_by_id

Break sequences in FASTA file into fragments
--------------------------------------------
File: ``scripts/fasta_fragments.py``

Taken from :ref:`lastz <lastz>`. Break sequences in a fasta file into fragments, so 
some kind of short read aligner can be used for further processing.

.. automodule:: fasta_fragments

General purpose tools
=====================

Primer3 wrapper
---------------
File: ``scrimer/primer3_connector.py``

.. automodule:: primer3_connector
   :members:

Convert CIGAR matches to sim4db 'script'
----------------------------------------
File: ``scripts/cigar_to_sim4db_scr.py``

.. automodule:: cigar_to_sim4db_scr

Count different bases in 5' end of reads in FASTQ
-------------------------------------------------
File: ``scripts/5prime_stats.py``

.. automodule:: 5prime_stats

Variant filters for PyVCF vcf_filter.py
---------------------------------------
File: ``scrimer/pyvcf_filters.py``

.. automodule:: pyvcf_filters
   :members:
   :undoc-members:
