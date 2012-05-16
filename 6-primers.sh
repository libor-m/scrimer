TOOLS=tabix:bedtools:primer3

GFF=60-gff-primers/lx3-primers.gff3

# for all selected variants design pcr and genotyping primers
# takes about a minute for 1000 selected variants, 5 MB gzipped vcf, 26 MB uncompressed genome, 5 MB gzipped gff
./design_primers.py 33-virtual-genome/lx3.fasta 33-virtual-genome/lx3.sorted.gff3.gz 51-variants-parallel/lx3-variants-selected.vcf.gz > $GFF

# it's quite vital to sort and index the file to use it in IGV
sortBed -i $GFF | bgzip > ${GFF%.*}.sorted.gff3.gz
tabix -f -p gff ${GFF%.*}.sorted.gff3.gz

# first, check primers with sole blat
# the next step is to set up a gfServer and use gfPcr