# amplicon-prep
Merges, trims, and creates sorted bam files to view paired-end amplicon sequencing data

#### Requirements (and command line commands to get them)
-conda

-sequence viewer: Tablet and IGV are free

all of the following can be installed by command line via "conda install -c bioconda <package name>"

-nextflow (builds the pipeline)

-pandaseq (merges paired end reads)

-cutadapt (trims linked primers)

-bwa (maps reads onto their target amplicon)

-samtools (creates the bam files that can be opened by a viewer)

#### Data format
user must create a design.csv file which tells nextflow where to find your files. As an example:

<img width="1278" alt="image" src="https://user-images.githubusercontent.com/20071084/213905715-40ff1092-3726-4e81-92c7-0076cb897e2e.png">

In this case, I have a folder with Genewiz sequencing (30-802983634/00_fastq/), which has compressed fastq files for the forward and reverse sequencing results for three samples (in this case, amplicon sequencing for CRISPR-Cas9 efficiency using three different gRNAs). "Fwd" and "Rev" take the paths to those reads. 

In that folder, I have made a fasta file with my target sequence.  In this case, they all have the same target. "Target" takes the path to that file. 

"Up" and "Down" are the sequences that bracket the amplicon of interest upstream and downstream, respectively.  This will generally be your priming sequences (forward primer, and the reverse complement of your reverse primer), but I'd recommend lengthening those sequences, so cutadapt can throw out any sequences that are the result of off-target binding.

"Name" is user defined. 

#### Running from terminal and expected output

after activating cutadapt environment (conda activate cutadaptenv), run nextflow run location1/grna_screen.nf -profile conda --directory location2, where the locations are the path to this script, and the directory containing your design.csv, respectively.  

I'd recommend making a folder for your design.csv, since this code will create an "out" folder which contains fastq files of intermediate steps (merged sequences, trimmed sequences, mapped sequences, bam files, and sorted bam files)

sorted bam files will be openable in various sequencing viewers. 

in getDoc there will be a spreadsheet that gives information on a base-by-base basis for a maximum depth of 8000 (to keep runtime reasonable: can be changed.  Check out samtools mpileup command for instructions)
  <img width="852" alt="image" src="https://user-images.githubusercontent.com/20071084/216056443-30364a72-9ae7-4b72-acf8-9a75ec7d571e.png">

In the above example, my sample is named "rosaO1O2". In the first position, my reference has an "A", and there are 19 reads over that area.  6 of those match "C", with no deletions or insertions.

Caveats
  -counts are off by 1-2. (the base-row will always be 1) 
  -currently does not account for inverted segments ("," in mpileup column). 

