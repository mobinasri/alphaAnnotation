version 1.0


import "./tasks/rDNA_annotation.wdl" as rDNA_annotation
import "./tasks/CenSatAnnotation.wdl" as finalizeCenSat
import "./tasks/gap_Annotation.wdl" as gapWorkflow
import "../alphaSat-HMMER/alphaSat-HMMER.wdl" as alphaSat

workflow centromereAnnotation {
    input {
        File fasta 
        File rDNAhmm_profile="../utilities/rDNA1.0.hmm"
        File hor_and_sf_bed
        File RMOut
        File HSatBed
        String fName=sub(basename(fasta), "\.(fa|fasta)(\.gz)?$", "")
        Boolean fix_sequence_ids = false
    }

    call formatAssembly {
        input:
            fasta=fasta,
            fName=fName,
            fix_sequence_ids=fix_sequence_ids
    }

    call rDNA_annotation.annotateRDNA as annotateRDNA {
        input:
            fasta=formatAssembly.formattedFasta,
            hmm_profile=rDNAhmm_profile

    }

    call gapWorkflow.annotateGaps as annotateGaps {
        input:
            fasta=formatAssembly.formattedFasta
    }

    call separate_hor_and_sf_bed {
        input:
            hor_and_sf_bed = hor_and_sf_bed
    }

    call createFile as emptyBed{
        input:
            content = "",
            filename = "empty.bed"
    }
    call alphaSat.summarize_alpha {
        input:
            as_hor_bed   = separate_hor_and_sf_bed.as_hor_bed,
            as_sf_bed    = separate_hor_and_sf_bed.as_sf_bed,
            assembly_id  = fName
    }

    call finalizeCenSat.cenSatAnnotation as cenSatAnnotation {
        input:
            RMOut=RMOut,
            aSatBed=summarize_alpha.as_summary_bed,
            aSatStrand=emptyBed.outFile,
            HSatBed=HSatBed,
            rDNABed=annotateRDNA.rDNAbed,
            gapBed=annotateGaps.gapBed
    }

    
    output {
        File cenSatAnnotations = cenSatAnnotation.cenSatAnnotations
    }

}


task createFile{
     input {
        String content = ""
        String filename = "mock.txt"
        # runtime configurations
        Int memSize=2
        Int threadCount=2
        Int diskSize=2
        String dockerImage="ubuntu:18.04"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mkdir output
        echo ~{content} > output/~{filename}
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
    }
    output {
        File outFile = glob("output/*")[0] 
    }
} 


task separate_hor_and_sf_bed {
     input {
        File hor_and_sf_bed
        String fName
        # runtime configurations
        Int memSize=2
        Int threadCount=2
        Int diskSize=2
        String dockerImage="ubuntu:18.04"
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mkdir output
        cat ~{hor_and_sf_bed} | grep -E 'S[0-9]+.*H[0-9]+' > output/~{fName}.as_hor.bed
        cat ~{hor_and_sf_bed} | grep -v -E 'S[0-9]+.*H[0-9]+' > output/~{fName}.as_sf.bed
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
    }
    output {
        File as_hor_bed = glob("output/*.as_hor.bed")[0]
        File as_sf_bed = glob("output/*.as_sf.bed")[0]
    }
}


task formatAssembly {
    input{
        File fasta
        String fName
        Boolean fix_sequence_ids
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # check that file is not empty and is in the correct format 
        if [ ! -s ~{fasta} ]; then echo "Fasta file is empty" ; exit 2 ; fi

        # unzip the fasta if it is zipped 
        if [[ ~{fasta} =~ \.gz$ ]] ; then
            gunzip -fc ~{fasta} > ~{fName}.fa 
        else 
            cat ~{fasta} > ~{fName}.fa 
        fi
        
        # check that file is nucleotide sequences and not proteins 
        
        #make sure there are no duplicate sequence_ids
        awk '/^>/ && seen[$1]++ {print "Error: Duplicate sequence IDs found:", $1; exit 1}' ~{fName}.fa

        ## Fix sequence IDs (if requested):
        ## Replace all instances of :,*,;, and - to underscores in sequence names 
        ## and record the correspondence to the unmodified names
        if [ ~{fix_sequence_ids} == true ]; then
            
            tee >( awk '{if (substr($1,1,1) == ">") { orig=substr($1,2); gsub(/[:*;-]/, "_", $1); print orig, substr($1,2) } }'  > ~{fName}.sequence_id_key.txt ) < ~{fName}.fa | awk '{ if (substr($1,1,1) == ">") { gsub(/[:*;-]/, "_", $1) } print $1 }'  > ~{fName}.formatted.fa

            #make sure there are no sequence name conflicts after renaming
            #this is so that when we will be renaming back, there will only be one correct option
            awk 'NR==FNR{map[$1]=$2; next} $2 in map && $1 != map[$2]{print "Error: "$2" found on a different row than "$1; exit 1}' ~{fName}.sequence_id_key.txt ~{fName}.sequence_id_key.txt
        else
            ## no renaming of sequence IDs, create empty file of original/new correspondance
            touch ~{fName}.sequence_id_key.txt

            ## rename file for output
            cp ~{fName}.fa ~{fName}.formatted.fa
        fi

    >>>

    output {
        File formattedFasta="~{fName}.formatted.fa"
        File sequence_id_key="~{fName}.sequence_id_key.txt"
    }

    runtime {
        docker: "ubuntu:18.04"
    }

}

task renameFinalOutputs {
    input{
        String fName
        File sequence_id_key
        File RMBed
        File RMOut
        File RMrmskBed
        File RMrmskAlignBed
        File RMMaskedFasta
        File as_hor_sf_bed
        File as_strand_bed
        File as_hor_bed
        File as_sf_bed
        File cenSatAnnotations
        File cenSatStrand
        File centromeres
        Boolean fix_sequence_ids
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        if [ ~{fix_sequence_ids} == true ]; then
            ## aggretate/summarized censat annotations
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{cenSatAnnotations} > ~{fName}.cenSat.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{cenSatStrand} > ~{fName}.SatelliteStrand.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{centromeres} > ~{fName}.active.centromeres.bed

            ## RepeatMasker: Not all RM outputs can be renamed in this way -- for example per-chrom files in the tar.gz and bigbed outputs
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{RMBed} > ~{fName}.RepeatMasker.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $5; next} {if ($5 in map) $5 = map[$1]} 1' ~{sequence_id_key} ~{RMOut} > ~{fName}.RepeatMasker.out
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{RMrmskBed} > ~{fName}.RepeatMasker.rmsk.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{RMrmskAlignBed} > ~{fName}.RepeatMasker.rmskAlign.bed
            zcat ~{RMMaskedFasta} \
                | awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} /^>/ {split($0,a,">"); if (a[2] in map) print ">"map[a[2]]; else print $0; next} {print}' ~{sequence_id_key} - \
                | gzip > ~{fName}.RepeatMasker.masked.fasta.gz
    
            ## ASat
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{as_hor_sf_bed} > ~{fName}.as_hor_sf.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{as_strand_bed} > ~{fName}.as_strand.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{as_hor_bed} > ~{fName}.as_hor.bed
            awk 'BEGIN {OFS="\t"} NR==FNR {map[$2] = $1; next} {if ($1 in map) $1 = map[$1]} 1' ~{sequence_id_key} ~{as_sf_bed} > ~{fName}.as_sf.bed            
        else 
            ## aggretate/summarized censat annotations
            cp ~{cenSatAnnotations} ~{fName}.cenSat.bed
            cp  ~{cenSatStrand} ~{fName}.SatelliteStrand.bed
            cp  ~{centromeres} ~{fName}.active.centromeres.bed

            ## RepeatMasker
            cp  ~{RMBed} ~{fName}.RepeatMasker.bed
            cp  ~{RMOut} ~{fName}.RepeatMasker.out
            cp  ~{RMrmskBed} ~{fName}.RepeatMasker.rmsk.bed
            cp  ~{RMrmskAlignBed} ~{fName}.RepeatMasker.rmskAlign.bed
            cp  ~{RMMaskedFasta} ~{fName}.RepeatMasker.masked.fasta.gz

            ## ASat
            cp  ~{as_hor_sf_bed} ~{fName}.as_hor_sf.bed
            cp  ~{as_strand_bed} ~{fName}.as_strand.bed
            cp  ~{as_hor_bed} ~{fName}.as_hor.bed
            cp  ~{as_sf_bed} ~{fName}.as_sf.bed            
        fi
    >>>

    output {
        File final_cenSatAnnotations="~{fName}.cenSat.bed"
        File final_cenSatStrand="~{fName}.SatelliteStrand.bed"
        File final_centromeres="~{fName}.active.centromeres.bed"

        File final_repeatMaskerBed="~{fName}.RepeatMasker.bed"
        File final_repeatMaskerOut="~{fName}.RepeatMasker.out"
        File final_rmskBed="~{fName}.RepeatMasker.rmsk.bed"
        File final_rmskAlignBed="~{fName}.RepeatMasker.rmskAlign.bed"
        File final_rmMaskedFasta="~{fName}.RepeatMasker.masked.fasta.gz"

        File final_as_hor_sf_bed="~{fName}.as_hor_sf.bed"
        File final_as_strand_bed="~{fName}.as_strand.bed"
        File final_as_hor_bed="~{fName}.as_hor.bed"
        File final_as_sf_bed="~{fName}.as_sf.bed"
    }

    runtime {
        docker: "ubuntu:18.04"
    }

}
