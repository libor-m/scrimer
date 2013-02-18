# Scrimer: 454 transcriptome data to SNaPShot primers
This project is a pipeline of bash commands driving custom Python scripts and 'industry standard' 
bioinformatics applications. After finding variable sites in the data, 
scrimer uses a (distant) reference genome to predict exons in the transcripts.
Predicted exons are used to find sequences suitable as PCR products that do not cross 
exon intron boundaries and contain some variable site. PCR and genotyping primers for 
SNaPshot screening of such site are created.

# Documentation
http://scrimer.rtfd.org

Documentation in sphinx format is in the ``doc`` directory.

# Status
Publication ready. Working on documentation. Publication is in review.


# License
Scrimer is licensed under GNU Affero General Public License.

If you use any part of the pipeline, please cite the appropriate publication (in prep now).

# Contact
    Libor Morkovsky
    Department of Zoology, Charles University in Prague
    morkovsk@natur.cuni.cz
