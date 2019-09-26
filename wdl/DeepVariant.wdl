# Use DeepVariant to generate VCF & gVCF for one sample
# ref: https://github.com/google/deepvariant/blob/r0.8/docs/deepvariant-gvcf-support.md
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
        # reference genome & index (if available, generated otherwise)
        File ref_fasta
        File? ref_fasta_idx

        # Genomic range(s) to call. Provide at most one of range or ranges_bed.
        # If neither is provided, calls the whole reference genome.
        String? range       # e.g. chr12:111766922-111817529
        File? ranges_bed

        # Read alignments - bam & bai (bai auto-generated if omitted)
        # The output vcf/gvcf filename is derived from the bam's.
        File bam
        File? bai
    }

    if (!defined(ref_fasta_idx)) {
        # index ref_fasta if needed
        call samtools_faidx {
            input:
                fasta = ref_fasta
        }
    }
    File ref_fasta_idx2 = select_first([ref_fasta_idx, samtools_faidx.fai])

    call make_examples {
        input:
            ref_fasta = ref_fasta,
            ref_fasta_idx = ref_fasta_idx2,
            range = range,
            ranges_bed = ranges_bed,
            bam = bam,
            bai = bai
    }

    call call_variants {
        input:
            examples = make_examples.examples
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
        curl -LSs --retry -O http://mirrors.kernel.org/ubuntu/pool/universe/s/samtools/samtools_1.9-4_amd64.deb
        dpkg -i samtools*.deb
        samtools faidx "~{fasta}"
    >>>

    output {
        File fai = "~{fasta}.fai"
    }

    runtime {
        docker: "ubuntu:disco"
    }
}

# DeepVariant make_examples
task make_examples {
    input {
        # reference genome
        File ref_fasta
        File ref_fasta_idx

        # bam & bai (bai auto-generated if omitted)
        File bam
        File? bai

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

        if [ ! -f "~{ref_fasta}.fai" ]; then
            cp "~{ref_fasta_idx}" "~{ref_fasta}.fai"
        fi

        if [ -z "~{bai}" ]; then
            # samtools deb from xenial base image of deepvariant-docker
            curl -LSs --retry -O http://mirrors.kernel.org/ubuntu/pool/universe/s/samtools/samtools_0.1.19-1ubuntu1_amd64.deb > samtools.deb
            dpkg -i samtools.deb
            samtools index "~{bam}"
        fi

        output_name=$(basename "~{bam}" .bam)
        mkdir examples/ gvcf/ logs/
        output_fn="examples/$output_name.tfrecord@$(nproc).gz"
        gvcf_fn="gvcf/$output_name.gvcf.tfrecord@$(nproc).gz"

        seq 0 $(( `nproc` - 1 )) | NO_GCE_CHECK=True parallel --halt 2 -t --results logs/ \
            "/opt/deepvariant/bin/make_examples --mode calling --ref ref.fasta --reads '~{bam}' --examples '$output_fn' --gvcf '$gvcf_fn' --task {} $range_arg $binsize_arg 2>&1" > /dev/null
    >>>

    runtime {
        docker: "gcr.io/deepvariant-docker/deepvariant:0.8.0"
        cpu: "64"
    }

    output {
        Array[File]+ examples = glob("examples/*")
        Array[File]+ gvcf_tfrecords = glob("gvcf/*")
        Array[File]+ logs = glob("logs/*")
    }
}

# DeepVariant call_variants (CPU)
task call_variants {
    input {
        Array[File]+ examples
        String model_type = "wgs" # "wgs", "wes", or "pacbio"
    }

    command <<<
        set -euxo pipefail

        n_examples=~{length(examples)}
        examples_dir=$(dirname "~{examples[0]}")
        output_name=$(basename "~{examples[0]}" .tfrecord) # FIXME

        NO_GCE_CHECK=True /opt/deepvariant/bin/call_variants \
            --outfile "output/$output_name.call_variants.tfrecord.gz" \
            --examples "$examples_dir/$output_name.tfrecord@$n_examples.gz" \
            --checkpoint "/opt/models/~{model_type}/model.ckpt"
    >>>

    runtime {
        docker: "gcr.io/deepvariant-docker/deepvariant:0.8.0"
        cpu: "64"
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

        if [ ! -f "~{ref_fasta}.fai" ]; then
            cp "~{ref_fasta_idx}" "~{ref_fasta}.fai"
        fi

        mkdir gvcf output
        n_gvcf_tfrecords="~{length(gvcf_tfrecords)}"
        output_name=$(basename "~{call_variants_output}" .call_variants.tfrecord.gz)
        NO_GCE_CHECK=True /opt/deepvariant/bin/postprocess_variants \
            --ref "~{ref_fasta}" --infile "~{call_variants_output}" \
            --nonvariant_site_tfrecord_path "gvcf/$output_name.gvcf.tfrecord@$n_gvcf_tfrecords.gz" \ FIXME
            --outfile "output/$output_name.vcf" \
            --gvcf_outfile "output/$output_name.gvcf"
    >>>

    runtime {
        docker: "gcr.io/deepvariant-docker/deepvariant:0.8.0"
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

    command <<<
        set -euxo pipefail
        curl -LSs --retry -O http://mirrors.kernel.org/ubuntu/pool/universe/h/htslib/tabix_1.9-10_amd64.deb
        dpkg -i tabix*.deb
        bgzip -@ 4 -c "~{file}" > "~{file}.gz"
    >>>

    output {
        File file_gz = "~{file}.gz"
    }

    runtime {
        docker: "ubuntu:disco"
        cpu: 4
    }
}
