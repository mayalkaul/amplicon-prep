# amplicon-prep
Merges, trims, and creates sorted bam files to view paired-end amplicon sequencing data

#### Requires
-conda
-pandaseq (merges paired end reads)
-cutadapt (trims linked primers)
-bwa (maps reads onto their target amplicon)
-samtools (creates the bam files that can be opened by a viewer -- Tablet and IGV, among many other, are free and decent)

#### Data format
user must create a design.csv file which tells nextflow where to find your files. As an example:

<img width="1278" alt="image" src="https://user-images.githubusercontent.com/20071084/213905715-40ff1092-3726-4e81-92c7-0076cb897e2e.png">

In this case, I have a folder with Genewiz sequencing (30-802983634/00_fastq/), which has compressed fastq files for the forward and reverse sequencing results for three samples (in this case, amplicon sequencing for CRISPR-Cas9 efficiency using three different gRNAs). "Fwd" and "Rev" take the paths to those reads. 

In that folder, I have made a fasta file with my target sequence.  In this case, they all have the same target. "Target" takes the path to that file. 

"Up" and "Down" are the sequences that bracket the amplicon of interest upstream and downstream, respectively.  This will generally be your priming sequences (forward primer, and the reverse complement of your reverse primer), but I'd recommend lengthening those sequences, so cutadapt can throw out any sequences that are the result of off-target binding.

"Name" is user defined. 

#### Running from terminal and expected output

after activating cutadapt environment (conda activate cutadaptenv), run nextflow run <location1>/grna_screen.nf -profile conda --directory <location2>, where the <location>s are the path to this script, and the directory containing your design.csv, respectively.  

I'd recommend making a folder for your design.csv, since this code will create an "out" folder which contains fastq files of intermediate steps (merged sequences, trimmed sequences, mapped sequences, bam files, and sorted bam files)

sorted bam files will be openable in various sequencing viewers. 
