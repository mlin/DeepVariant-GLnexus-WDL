# Use DeepVariant to generate VCF & gVCF for one sample
# ref: https://github.com/google/deepvariant/blob/r0.9/docs/deepvariant-gvcf-support.md
#
#              +----------------------------------------------------------------------------+
#              |                                                                            |
#              |  DeepVariant.wdl                                                           |
#              |                                                                            |
#              |  +-----------------+    +-----------------+    +------------------------+  |
# sample.bam   |  |                 |    |                 |    |                        |  |
#  genome.fa ----->  make_examples  |---->  call_variants  |---->  postprocess_variants  |-----> gVCF
#      range   |  |                 |    |                 |    |                        |  |
#              |  +-----------------+    +--------^--------+    +------------------------+  |
#              |                                  |                                         |
#              |                                  |                                         |
#              +----------------------------------|-----------------------------------------+
#                                                 |
#                                        DeepVariant Model
version 1.0

workflow DeepVariant {
    input {
        # reference genome & index (if available; generated otherwise)
        File ref_fasta
        File? ref_fasta_idx

        # Genomic range(s) to call. Provide at most one of range or ranges_bed.
        # If neither is provided, calls the whole reference genome.
        String? range       # e.g. chr12:111760000-111820000
        File? ranges_bed

        # DeepVariant model type -- wgs, wes, or pacbio
        String model_type = "wgs"

        # Read alignments - bam & bai (bai generated if omitted)
        # The output vcf/gvcf filename is derived from the bam's.
        File bam
        File? bai

        # gVCF advanced setting
        Int? gvcf_gq_binsize
    }

    if (!defined(ref_fasta_idx)) {
        # index ref_fasta if needed
        call samtools_faidx {
            input:
                fasta = ref_fasta
        }
    }
    File ref_fasta_idx2 = select_first([ref_fasta_idx, samtools_faidx.fai])

    if (!defined(bai)) {
        call samtools_index {
            input:
                bam = bam
        }
    }
    File bai2 = select_first([bai, samtools_index.bai])

    call make_examples {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_idx = ref_fasta_idx2,
            range = range,
            ranges_bed = ranges_bed,
            bam = bam,
            bai = bai2,
            gvcf_gq_binsize = gvcf_gq_binsize
    }

    call call_variants {
        input:
            examples = make_examples.examples,
            model_type = model_type
    }

    call postprocess_variants {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_idx = ref_fasta_idx2,
            call_variants_output = call_variants.call_variants_output,
            gvcf_tfrecords = make_examples.gvcf_tfrecords
    }

    call bgzip as bgzip_vcf {
        input:
            file = postprocess_variants.vcf
    }

    call bgzip as bgzip_gvcf {
        input:
            file = postprocess_variants.gvcf
    }

    output {
        File vcf_gz = bgzip_vcf.file_gz
        File gvcf_gz = bgzip_gvcf.file_gz
    }
}

task samtools_faidx {
    input {
        File fasta
    }
    command <<<
        set -euxo pipefail
        samtools faidx "~{fasta}"
    >>>
    output {
        File fai = "~{fasta}.fai"
    }
    runtime {
        docker: "biocontainers/samtools:v1.9-4-deb_cv1"
    }
}

task samtools_index {
    input {
        File bam
    }
    command <<<
        set -euxo pipefail
        samtools index -@ 4 "~{bam}"
    >>>
    output {
        File bai = "~{bam}.bai"
    }
    runtime {
        docker: "biocontainers/samtools:v1.9-4-deb_cv1"
        cpu: 4
    }
}

# DeepVariant make_examples
task make_examples {
    input {
        # reference genome
        File ref_fasta
        File ref_fasta_idx

        # bam & bai
        File bam
        File bai

        # Genomic range(s) to run on. Provide at most one of range or ranges_bed.
        # If neither is provided, calls the whole reference genome.
        String? range       # e.g. chr12:111766922-111817529
        File? ranges_bed

        # advanced gVCF setting
        Int? gvcf_gq_binsize
    }

    command <<<
        range_arg=""
        if [ -n "~{range}" ]; then
            range_arg="--regions ~{range}"
        fi
        if [ -n "~{ranges_bed}" ]; then
            range_arg="--regions ~{ranges_bed}"
        fi
        binsize_arg=""
        if [ -n "~{gvcf_gq_binsize}" ]; then
            binsize_arg="--gvcf_gq-binsize ~{gvcf_gq_binsize}"
        fi

        set -euxo pipefail
        export SHELL=/bin/bash

        output_name=$(basename "~{bam}" .bam)
        mkdir examples/ gvcf/ logs/
        output_fn="examples/$output_name.tfrecord@$(nproc).gz"
        gvcf_fn="gvcf/$output_name.gvcf.tfrecord@$(nproc).gz"

        seq 0 $(( `nproc` - 1 )) | NO_GCE_CHECK=True parallel --halt 2 -t \
            "/opt/deepvariant/bin/make_examples --mode calling --ref '~{ref_fasta}' --reads '~{bam}' --examples '$output_fn' --gvcf '$gvcf_fn' --task {} $range_arg $binsize_arg 2>&1" > /dev/null
    >>>

    runtime {
        docker: "gcr.io/deepvariant-docker/deepvariant:0.9.0"
        cpu: 64
    }

    output {
        Array[File]+ examples = glob("examples/*.gz")
        Array[File]+ gvcf_tfrecords = glob("gvcf/*")
    }
}

# DeepVariant call_variants (CPU)
task call_variants {
    input {
        Array[File]+ examples
        String model_type
    }

    command <<<
        set -euxo pipefail

        n_examples=~{length(examples)}
        examples_dir=$(dirname "~{examples[0]}")
        output_name=$(basename "~{examples[0]}" | grep -oP '.+(?=\.tfrecord-\d+-of-\d+.gz)')

        mkdir output
        NO_GCE_CHECK=True /opt/deepvariant/bin/call_variants \
            --outfile "output/$output_name.call_variants.tfrecord.gz" \
            --examples "$examples_dir/$output_name.tfrecord@$n_examples.gz" \
            --checkpoint "/opt/models/~{model_type}/model.ckpt"
    >>>

    runtime {
        docker: "gcr.io/deepvariant-docker/deepvariant:0.9.0"
        cpu: 64
    }

    output {
        File call_variants_output = glob("output/*.gz")[0]
    }
}

# DeepVariant postprocess_variants
task postprocess_variants {
    input {
        File ref_fasta
        File ref_fasta_idx

        Array[File]+ gvcf_tfrecords
        
        File call_variants_output
    }

    command <<<
        set -euxo pipefail

        mkdir gvcf output
        n_gvcf_tfrecords="~{length(gvcf_tfrecords)}"
        gvcf_tfrecords_dir=$(dirname "~{gvcf_tfrecords[0]}")
        output_name=$(basename "~{call_variants_output}" .call_variants.tfrecord.gz)
        NO_GCE_CHECK=True /opt/deepvariant/bin/postprocess_variants \
            --ref "~{ref_fasta}" --infile "~{call_variants_output}" \
            --nonvariant_site_tfrecord_path "$gvcf_tfrecords_dir/$output_name.gvcf.tfrecord@$n_gvcf_tfrecords.gz" \
            --outfile "output/$output_name.vcf" \
            --gvcf_outfile "output/$output_name.gvcf"
    >>>

    runtime {
        docker: "gcr.io/deepvariant-docker/deepvariant:0.9.0"
        cpu: 2
    }

    output {
        File vcf = glob("output/*.vcf")[0]
        File gvcf = glob("output/*.gvcf")[0]
    }
}

task bgzip {
    input {
        File file
    }
    String filename = basename(file)
    command <<<
        bgzip -@ 4 -c "~{file}" > "~{filename}.gz"
    >>>
    output {
        File file_gz = "~{filename}.gz"
    }
    runtime {
        docker: "biocontainers/tabix:v1.9-11-deb_cv1"
        cpu: 4
    }
}
