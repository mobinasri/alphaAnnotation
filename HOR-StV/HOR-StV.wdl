version 1.0

workflow hor_stv_workflow {

	call hor_stv

	output {
		File hor_stv_bed     = hor_stv.hor_stv_bed
		File hor_stv_row_bed = hor_stv.hor_stv_row_bed
		File hor_stv_stats   = hor_stv.hor_stv_stats
	}

	meta {
		author: "Julian Lucas"
		email: "juklucas@ucsc.edu"
		description: "Calls [StV](https://github.com/fedorrik/stv) for Alpha Sat. HOR Structrual Variant (StV) prediction using HOR-monomer annotation"
	}
}

task hor_stv {
	input {
		File asm_monomer_bed
		String sample_name
		String asm_tag

		Int memSizeGB   = 4
		Int threadCount = 1
		Int diskSizeGB  = 32
	}

	parameter_meta {
		asm_monomer_bed: "BED file with HOR monomers. Can be generated by [HumAS-HMMER](https://github.com/kmiga/alphaAnnotation/blob/main/alphaSat-HMMER/alphaSat-HMMER.wdl) for example"
		sample_name: "sample name to attach to output files"
		asm_tag: "tag to attach to output files"
	}

	command <<<
		set -eux -o pipefail

		## Sort just in case we didn't pass a sorted bed file.
		bedtools sort -i ~{asm_monomer_bed} > sorted.bed


		## Find HOR Structural Variants
		stv.sh sorted.bed

		## Rename files so they are easier to keep track of
		mv stv.bed ~{sample_name}_~{asm_tag}_HOR_StV.bed
		mv stv_row.bed ~{sample_name}_~{asm_tag}_HOR_StV_row.bed
		mv stv_stats.tsv ~{sample_name}_~{asm_tag}_HOR_StV_stats.txt
	>>>

	output {
		File hor_stv_bed     = "~{sample_name}_~{asm_tag}_HOR_StV.bed"
		File hor_stv_row_bed = "~{sample_name}_~{asm_tag}_HOR_StV_row.bed"
		File hor_stv_stats   = "~{sample_name}_~{asm_tag}_HOR_StV_stats.txt"
	}

	runtime {
		memory: memSizeGB + " GB"
		cpu: threadCount
		disks: "local-disk " + diskSizeGB + " SSD"
		docker: "juklucas/hor_stv@sha256:839cc33b2ae8f19a67f50ffbe07882f3937f6056c4ffca51eefdfbbc1f4580fe"
		preemptible: 1
	}
}