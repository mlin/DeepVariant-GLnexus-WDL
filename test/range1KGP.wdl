# Slice a genomic range from the (deep resequencing of the) 1000 Genomes Project cohort CRAM files,
# generate per-sample gVCF using DeepVariant, then merge them to a pVCF using GLnexus.

version 1.0

import "../wdl/DeepVariant_GLnexus.wdl" as dvglx

workflow range1KGP {
    input {
        # ERR ID and sample, for constructing CRAM file URI
        # example: "ERR3240114/HG00096"
        # table from which these can be derived: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index
        # e.g. grep -v ^\# 1000G_2504_high_coverage.sequence.index | cut -f3,10 --output-delimiter / | awk '{printf("\"%s\",\n",$0)}'
        Array[String]+ ERR_slash_sample

        # GRCh38 reference
        File ref_fasta
        File ref_fasta_idx

        # desired genomic range e.g. "chr17" or "chr12:111760000-111820000"
        String range
        String range_name = range

        # used in output filenames
        String cohort_name = "1KGP"
    }

    call samtools_slice_1000G_bam as slicer {
        input:
            ERR_slash_sample = ERR_slash_sample,
            range = range,
            range_name = range_name,
            ref_fasta = ref_fasta,
            ref_fasta_idx = ref_fasta_idx
    }

    call dvglx.DeepVariant_GLnexus {
        input:
            bam = slicer.bams,
            range = range,
            ref_fasta = ref_fasta,
            ref_fasta_idx = ref_fasta_idx,
            output_name = "~{cohort_name}.~{range_name}"
    }

    output {
        Array[File]+ bams = slicer.bams
        Array[File]+ gvcfs = DeepVariant_GLnexus.gvcf_gz
        File pvcf_gz = DeepVariant_GLnexus.pvcf_gz
    }
}

task samtools_slice_1000G_bam {
    # Slice genomic range BAMs out of the CRAMs in the 1000genomes public S3 bucket
    #
    # We parallelize this internally to amortize the CRAM refget operation

    input {
        Array[String]+ ERR_slash_sample
        String range
        String range_name

        File ref_fasta
        File ref_fasta_idx

        String timeout = "1h"
        Int cpu = 32
    }

    # oversubscribe CPUs to anticipate networking bottlenecks
    Int P = cpu*2

    command <<<
        set -euxo pipefail
        apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y samtools parallel

        cat << "EOF" > do1
        sample=$(basename "$1")
        rm -f "$sample.final.cram.crai"
        outfn="$sample.~{range_name}.bam"
        timeout "~{timeout}" samtools view -b -T "~{ref_fasta}" \
            -o "$outfn" \
            "https://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/data/$1.final.cram" \
            "~{range}"
        echo "$outfn"
        EOF
        cat "~{write_lines(ERR_slash_sample)}" | \
            parallel -t -k -P ~{P} --retries 3 --halt 2 \
            bash -eux -o pipefail do1 {} > slice_filenames
    >>>

    output {
        Array[File]+ bams = read_lines("slice_filenames")
    }

    runtime {
        docker: "ubuntu:19.10"
        cpu: cpu
    }
}
