# Slice a genomic range from the DeepVariant gVCF files for the 1000 Genomes Project, then use
# GLnexus to merge them into a pVCF.

version 1.0

import "../wdl/GLnexus.wdl" as lib

workflow GLnexus_range1KGP {
    input {
        Array[String]+ samples

        # desired genomic range e.g. "chr17" or "chr12:111760000-111820000"
        String range
        String range_name = range

        # used in output filenames
        String cohort_name = "1KGP"
    }

    scatter (sample in samples) {
        call tabix_slice_1KGP_gVCF as slicer {
            input:
                sample = sample,
                range = range,
                range_name = range_name
        }
    }

    call lib.GLnexus {
        input:
            gvcf = slicer.gvcf_gz,
            range = range,
            output_name = "~{cohort_name}.~{range_name}"
    }
}

task tabix_slice_1KGP_gVCF {
    # Slice genomic range out of the DeepVariant gVCFs stored under gs://brain-genomics-public/

    input {
        String sample
        String range
        String range_name

        String timeout = "1h"
    }

    String outfile = "~{sample}.~{range_name}.gvcf.gz"

    command <<<
        set -euxo pipefail
        timeout "~{timeout}" tabix -h \
            "gs://brain-genomics-public/research/cohort/1KGP/dv_vcf/v1/~{sample}.dv0.8.0.g.vcf.gz" "~{range}" \
            | bgzip -c > "~{outfile}"
    >>>

    output {
        File gvcf_gz = outfile
    }

    runtime {
        docker: "biocontainers/tabix:v1.9-11-deb_cv1"
        cpu: 1
        maxRetries: 2
    }
}
