TOOLS=tabix:bedtools:primer3

GFF=60-gff-primers/lx3-primers.gff3

# for all selected variants design pcr and genotyping primers
# takes about a minute for 1000 selected variants, 5 MB gzipped vcf, 26 MB uncompressed genome, 5 MB gzipped gff
./design_primers.py 33-virtual-genome/lx3.fasta 33-virtual-genome/lx3.sorted.gff3.gz 51-variants-parallel/lx3-variants-selected.vcf.gz > $GFF

# it's quite vital to sort and index the file to use it in IGV
sortBed -i $GFF | bgzip > ${GFF%.*}.sorted.gff3.gz
tabix -f -p gff ${GFF%.*}.sorted.gff3.gz

###
# use blat and isPcr to map previously manually designed primers
###

# convert fa to 2bit
faToTwoBit 33-virtual-genome/lx3.fasta 61-check-primers/lx3.2bit

# use isPcr to check the products
isPcr 61-check-primers/lx3.2bit 61-check-primers/manual_pcr.ispcr 61-check-primers/manual_pcr.ispcr.fa

# isPcr with psl output to be easily loaded into IGV
isPcr -out=psl 61-check-primers/lx3.2bit 61-check-primers/manual_pcr.ispcr 61-check-primers/manual_pcr.ispcr.psl

# blat the manually designed primers againts the virtual genome
# primer sequences are short, so the blat default parameters
# have to be changed a bit
blat -minScore=15 -tileSize=6 -maxIntron=0 61-check-primers/lx3.2bit 61-check-primers/manual_pcr.fa 61-check-primers/manual_pcr.psl

# agrep is quite enough for simple checks on assemblies of this size
SEQ=GCACATTTCATGGTCTCCAA
agrep $SEQ 0a-jp-newbler-contigs/lu??_contigs.fasta|grep $SEQ

