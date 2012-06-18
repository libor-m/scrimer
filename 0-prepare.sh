# merge the same individuals from different runs to a single file
# should be done for sam files using the appropriate @RG headers
# http://www.broadinstitute.org/gsa/wiki/index.php/Frequently_Asked_Questions#My_BAM_file_doesn.27t_have_read_group_and_sample_information.__Do_I_really_need_it.3F
# http://seqanswers.com/forums/showthread.php?t=4180
# http://www.broadinstitute.org/gsa/wiki/index.php/Input_files_for_the_GATK#Sequencing_Reads
# .. explained: 
# READ GROUP is unique id for instrument/library/run/sample combination
# the metadata are given in the header, the assignment to a readgroup is given for each alignment
# resolution 
# >> map the MID separated samples one by original file
# then samtools merge -rh (into one file per species?) with @RG:
# ID:original filename
# PL:LS454
# LB:can be misused to distinguish LM and LL? (library == SM in our data) -- used by samtools rmdup command, has to be set to the correct value
# SM:luXX

# download and prepare chicken genome
wget -m -np ftp://hgdownload.cse.ucsc.edu/goldenPath/galGal3/bigZips/
find hgdownload.cse.ucsc.edu -type f -exec mv {} . \;
rm -r hgdownload.cse.ucsc.edu
# check data integrity
md5sum -c md5sum.txt
cat *.md5|md5sum -c



faToTwoBit taeGut1.fa twoBit/taeGut1.2bit