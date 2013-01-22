#!/usr/bin/env python
"""Filter a FASTA, FASTQ or SSF file with IDs from a tabular file.

Takes six command line options, tabular filename, ID column numbers (comma
separated list using one based counting), input filename, input type (e.g.
FASTA or SFF) and two output filenames (for records with and without the
given IDs, same format as input sequence file).

If either output filename is just a minus sign, that file is not created.
This is intended to allow output for just the matched (or just the non-matched)
records.

When filtering an SFF file, any Roche XML manifest in the input file is
preserved in both output files.

Note in the default NCBI BLAST+ tabular output, the query sequence ID is
in column one, and the ID of the match from the database is in column two.
Here sensible values for the column numbers would therefore be "1" or "2".

This tool is a short Python script which requires Biopython 1.54 or later
for SFF file support. If you use this tool in scientific work leading to a
publication, please cite the Biopython application note:

Cock et al 2009. Biopython: freely available Python tools for computational
molecular biology and bioinformatics. Bioinformatics 25(11) 1422-3.
http://dx.doi.org/10.1093/bioinformatics/btp163 pmid:19304878.

This script is copyright 2010 by Peter Cock, SCRI, UK. All rights reserved.
See accompanying text file for licence details (MIT/BSD style).

This is version 0.0.1 of the script.
"""
import sys

#Dan Blankenberg

class fastaSequence( object ):
    def __init__( self ):
        self.identifier = None
        self.sequence = '' #holds raw sequence string: no whitespace
    def __len__( self ):
        return len( self.sequence )
    def __str__( self ):
        return "%s\n%s\n" % ( self.identifier, self.sequence )

class fastaReader( object ):
    def __init__( self, fh ):
        self.file = fh
    def close( self ):
        return self.file.close()
    def next( self ):
        line = self.file.readline()
        #remove header comment lines
        while line and line.startswith( '#' ):
            line = self.file.readline()
        if not line:
            raise StopIteration 
        assert line.startswith( '>' ), "FASTA headers must start with >"
        rval = fastaSequence()
        rval.identifier = line.strip()
        offset = self.file.tell()
        while True:
            line = self.file.readline()
            if not line or line.startswith( '>' ):
                if line:
                    self.file.seek( offset ) #this causes sequence id lines to be read twice, once to determine previous sequence end and again when getting actual sequence; can we cache this to prevent it from being re-read?
                return rval
            #454 qual test data that was used has decimal scores that don't have trailing spaces 
            #so we'll need to parse and build these sequences not based upon de facto standards
            #i.e. in a less than ideal fashion
            line = line.rstrip()
            if ' ' in rval.sequence or ' ' in line:
                rval.sequence = "%s%s " % ( rval.sequence, line )
            else:
                rval.sequence += line
            offset = self.file.tell()
    def __iter__( self ):
        while True:
            yield self.next()

class fastaWriter( object ):
    def __init__( self, fh ):
        self.file = fh
    def write( self, fastq_read ):
        #this will include color space adapter base if applicable
        self.file.write( ">%s\n%s\n" % ( fastq_read.identifier[1:], fastq_read.sequence ) )
    def close( self ):
        return self.file.close() 

def stop_err(msg, err=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(err)

def main():
    #Parse Command Line
    try:
        tabular_file, cols_arg, in_file, seq_format, out_positive_file, out_negative_file = sys.argv[1:]
    except ValueError:
        stop_err("Expected six arguments, got %i:\n%s" % (len(sys.argv)-1, " ".join(sys.argv)))
    try:
        columns = [int(arg)-1 for arg in cols_arg.split(",")]
    except ValueError:
        stop_err("Expected list of columns (comma separated integers), got %s" % cols_arg)

    if out_positive_file == "-" and out_negative_file == "-":
        stop_err("Neither output file requested")


    #Read tabular file and record all specified identifiers
    ids = set()
    handle = open(tabular_file, "rU")
    if len(columns)>1:
        #General case of many columns
        for line in handle:
            if line.startswith("#"):
                #Ignore comments
                continue
            parts = line.rstrip("\n").split("\t")
            for col in columns:
                ids.add(parts[col])
        print "Using %i IDs from %i columns of tabular file" % (len(ids), len(columns))
    else:
        #Single column, special case speed up
        col = columns[0]
        for line in handle:
            if not line.startswith("#"):
                ids.add(line.rstrip("\n").split("\t")[col])
        print "Using %i IDs from tabular file" % (len(ids))
    handle.close()


    if seq_format.lower()=="sff":
        #Now write filtered SFF file based on IDs from BLAST file
        try:
            from Bio.SeqIO.SffIO import SffIterator, SffWriter
        except ImportError:
            stop_err("Requires Biopython 1.54 or later")

        try:
            from Bio.SeqIO.SffIO import ReadRocheXmlManifest
        except ImportError:
            #Prior to Biopython 1.56 this was a private function
            from Bio.SeqIO.SffIO import _sff_read_roche_index_xml as ReadRocheXmlManifest
        in_handle = open(in_file, "rb") #must be binary mode!
        try:
            manifest = ReadRocheXmlManifest(in_handle)
        except ValueError:
            manifest = None
        #This makes two passes though the SFF file with isn't so efficient,
        #but this makes the code simple.
        if out_positive_file != "-":
            out_handle = open(out_positive_file, "wb")
            writer = SffWriter(out_handle, xml=manifest)
            in_handle.seek(0) #start again after getting manifest
            pos_count = writer.write_file(rec for rec in SffIterator(in_handle) if rec.id in ids)
            out_handle.close()
        if out_negative_file != "-":
            out_handle = open(out_negative_file, "wb")
            writer = SffWriter(out_handle, xml=manifest)
            in_handle.seek(0) #start again
            neg_count = writer.write_file(rec for rec in SffIterator(in_handle) if rec.id not in ids)
            out_handle.close()
        #And we're done
        in_handle.close()
        #At the time of writing, Galaxy doesn't show SFF file read counts,
        #so it is useful to put them in stdout and thus shown in job info.
        if out_positive_file != "-" and out_negative_file != "-":
            print "%i with and %i without specified IDs" % (pos_count, neg_count)
        elif out_positive_file != "-":
            print "%i with specified IDs" % pos_count
        elif out_negative_file != "-":
            print "%i without specified IDs" % neg_count
    elif seq_format.lower()=="fasta":
        #Write filtered FASTA file based on IDs from tabular file
        reader = fastaReader(open(in_file, "rU"))
        if out_positive_file != "-" and out_negative_file != "-":
            print "Generating two FASTA files"
            positive_writer = fastaWriter(open(out_positive_file, "w"))
            negative_writer = fastaWriter(open(out_negative_file, "w"))
            for record in reader:
                #The [1:] is because the fastaReader leaves the > on the identifer.
                if record.identifier and record.identifier.split()[0][1:] in ids:
                    positive_writer.write(record)
                else:
                    negative_writer.write(record)
            positive_writer.close()
            negative_writer.close()
        elif out_positive_file != "-":
            print "Generating matching FASTA file"
            positive_writer = fastaWriter(open(out_positive_file, "w"))
            for record in reader:
                #The [1:] is because the fastaReader leaves the > on the identifer.
                if record.identifier and record.identifier.split()[0][1:] in ids:
                    positive_writer.write(record)
            positive_writer.close()
        elif out_negative_file != "-":
            print "Generating non-matching FASTA file"
            negative_writer = fastaWriter(open(out_negative_file, "w"))
            for record in reader:
                #The [1:] is because the fastaReader leaves the > on the identifer.
                if not record.identifier or record.identifier.split()[0][1:] not in ids:
                    negative_writer.write(record)
            negative_writer.close()
    elif seq_format.lower().startswith("fastq"):
        #Write filtered FASTQ file based on IDs from tabular file
        from galaxy_utils.sequence.fastq import fastqReader, fastqWriter
        reader = fastqReader(open(in_file, "rU"))
        if out_positive_file != "-" and out_negative_file != "-":
            print "Generating two FASTQ files"
            positive_writer = fastqWriter(open(out_positive_file, "w"))
            negative_writer = fastqWriter(open(out_negative_file, "w"))
            for record in reader:
                #The [1:] is because the fastaReader leaves the @ on the identifer.
                if record.identifier and record.identifier.split()[0][1:] in ids:
                    positive_writer.write(record)
                else:
                    negative_writer.write(record)
            positive_writer.close()
            negative_writer.close()
        elif out_positive_file != "-":
            print "Generating matching FASTQ file"
            positive_writer = fastqWriter(open(out_positive_file, "w"))
            for record in reader:
                #The [1:] is because the fastaReader leaves the @ on the identifer.
                if record.identifier and record.identifier.split()[0][1:] in ids:
                    positive_writer.write(record)
            positive_writer.close()
        elif out_negative_file != "-":
            print "Generating non-matching FASTQ file"
            negative_writer = fastqWriter(open(out_negative_file, "w"))
            for record in reader:
                #The [1:] is because the fastaReader leaves the @ on the identifer.
                if not record.identifier or record.identifier.split()[0][1:] not in ids:
                    negative_writer.write(record)
            negative_writer.close()
    else:
        stop_err("Unsupported file type %r" % seq_format)

if __name__ == "__main__": main()
