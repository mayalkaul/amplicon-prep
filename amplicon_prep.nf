
//////////////////////////////////////////////////////////
/*						Parameters						*/
//////////////////////////////////////////////////////////
// params.directory='/Users/maya/Desktop/30-802983634/00_fastq/' // <-- alternative to passing directory in from command line: replace with your location

//////////////////////////////////////////////////////////
/*						Design file						*/
//////////////////////////////////////////////////////////
params.design=params.directory+'design.csv'
Channel
	.fromPath(params.design)
	.splitCsv(header:true,sep:',')
	.map{ row->[row.Name,
				row.Target,
				file(row.Fwd),
				file(row.Rev),
				row.Up,
				row.Down]}
	.set{ch_nabbed}

//////////////////////////////////////////////////////////
/*						Merge 							*/
//////////////////////////////////////////////////////////
process merge {
	publishDir params.directory+'out/merge'
	input:
	tuple val(name),path(target),path(reads1),path(reads2),val(up),val(down) from ch_nabbed
	output: 
	tuple val(name),path(target),path("${name}_merged.fasta"),val(up),val(down) into ch_merged
	script:
	"""
	pandaseq -f "${reads1}" -r "${reads2}" >> "${name}_merged.fasta"
	"""
}

//////////////////////////////////////////////////////////
/*						Trim/filter						*/
//////////////////////////////////////////////////////////
process trim{
	publishDir params.directory+'out/trim'
	input:
	tuple val(name),path(target),path(reads),val(up),val(down) from ch_merged
	output:
	tuple val(name),path(target),path("${name}_trim.fasta") into ch_trimmed
	script: 
	// because using -a, it's important that the upstream and downstream sequences are anchored
	"""
	cutadapt -a $up...$down -o ${name}_trim.fasta $reads --discard-untrimmed
	"""
}

//////////////////////////////////////////////////////////
/*						aligning						*/
//////////////////////////////////////////////////////////

process map{
	publishDir params.directory+'out/map'
	input:
	tuple val(name),path(target), path(reads) from ch_trimmed
	output:
	tuple val(name),path("${name}_aln_pe.sam") into ch_sam
	script:
	"""
    bwa index '${target}'
    bwa mem '${target}' $reads > ${name}_aln_pe.sam
	"""
}

process makeBam {
    publishDir params.directory+'out/makeBam'
    input:
    tuple val(name), path(aln_pe) from ch_sam
    output:
    tuple val(name), path("${name}_aln_pe.bam")into ch_bam
    script:
    """
    samtools view -h -b -S $aln_pe > ${name}_aln_pe.bam
    """

}
process sortBam {
    publishDir params.directory+'out/sortBam'
    input:
    tuple val(name), path(aln_pe_bam) from ch_bam
    output:
    tuple val(name), path("${name}_sorted.bam"), path("${name}_sorted.bam.bai") into ch_sorted_bam
    script:
    """
    samtools sort $aln_pe_bam -o ${name}_sorted.bam
    samtools index ${name}_sorted.bam ${name}_sorted.bam.bai
    """
}
