Component documentation
#######################

Components bound to the pipeline
================================

Transfer annotation from reference genome to the mapped reads
-------------------------------------------------------------
File: ``liftover.py``

.. automodule:: liftover

Create a contig scaffold with help of spliced mapping
-----------------------------------------------------
File: ``scaffold.py``

.. automodule:: scaffold

Design primers using sequences, annotations and selected variants
-----------------------------------------------------------------
File: ``design_primers.py``

.. automodule:: design_primers

Export primers from gff3 to FASTA or isPcr format
-------------------------------------------------
File: ``extract_primers.py``

.. automodule:: extract_primers


Tools for FASTA/FASTQ manipulation
==================================

Given pairs of matching sequences, find the longest one
-------------------------------------------------------
File: ``find_redundant_sequences.py``

.. automodule:: find_redundant_sequences

Filter sequences in FASTQ file based on their order in the file
---------------------------------------------------------------
File: ``fastq_kill_lines.py``

.. automodule:: fastq_kill_lines

Filter sequences in FASTQ, FASTA based on their order in the file
-----------------------------------------------------------------
File: ``seq_filter_by_id.py``

Taken from [BioPython]_. FASTA and FASTQ readers are pasted in the file.

.. automodule:: seq_filter_by_id

Break sequences in FASTA file into fragments
--------------------------------------------
File: ``fastq_kill_lines.py``

Taken from [lastz]_.

.. automodule:: fasta_fragments

Other tools
===========

Primer3 wrapper
---------------
File: ``primer3_connector.py``

.. automodule:: primer3_connector
   :members:

Convert CIGAR matches to sim4db 'script'
----------------------------------------
File: ``cigar_to_sim4db_scr.py``

.. automodule:: cigar_to_sim4db_scr

Count different bases in 5' end of reads in FASTQ
-------------------------------------------------
File: ``5prime_stats.py``

.. automodule:: 5prime_stats

Variant filters for PyVCF vcf_filter.py
---------------------------------------
File: ``pyvcf_filters.py``

.. automodule:: pyvcf_filters
   :members:
   :undoc-members:
