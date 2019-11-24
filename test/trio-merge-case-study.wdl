# Inspired by https://github.com/google/deepvariant/blob/r0.9/docs/trio-merge-case-study.md

version 1.0

import "../wdl/DeepVariant_GLnexus.wdl" as dvglx

workflow trio_merge_case_study {
    input {
        Array[File]+ bams
        Array[File]+ bais

        File ref_fasta_gz

        File targets_bed
    }

    call gunzip as ref {
        input:
            gz = ref_fasta_gz
    }

    call dvglx.DeepVariant_GLnexus {
        input:
            bam = bams,
            bai = bais,
            ranges_bed = targets_bed,
            ref_fasta = ref.uncompressed,
            model_type="WES",
            glnexus_config="DeepVariantWES",
            output_name="AshkenaziTrio"
    }

    output {
        Array[File] gvcf_gz = DeepVariant_GLnexus.gvcf_gz
        File pvcf_gz = DeepVariant_GLnexus.pvcf_gz
    }
}

task gunzip {
    input {
        File gz
    }
    command {
        mkdir out/
        gunzip -c "~{gz}" > "out/$(basename "~{gz}" .gz)"
    }
    output {
        File uncompressed = glob("out/*")[0]
    }
}
