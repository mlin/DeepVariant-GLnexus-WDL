version 1.0

task GLnexus {
    input {
        Array[File]+ gvcf
        String output_name

        String config = "DeepVariantWGS"
        File? config_yml

        String? range
        File? ranges_bed

        Boolean squeeze = false

        Int cpu = 16
        Int memoryGB = cpu*4
        Int diskGB = 3*floor(size(gvcf, "GB"))+1
    }

    command <<<
        set -euxo pipefail

        bed_arg=""
        if [ -n "~{ranges_bed}" ]; then
            bed_arg="--bed ~{ranges_bed}"
        elif [ -n "~{range}" ]; then
            echo "~{range}" | tr :- '\t' | tr -d , > range.bed
            if [ -z "$(cut -sf 2,3 range.bed)" ]; then
                echo -e "~{range}\t0\t999999999" > range.bed
            fi
            bed_arg="--bed range.bed"
        fi

        squeeze_cmd="cat"
        squeeze_arg=""
        if [ "~{squeeze}" == "true" ]; then
            squeeze_cmd="spvcf squeeze"
            squeeze_arg="--squeeze"
        fi

        glnexus_cli \
            --config "~{if defined(config_yml) then config_yml else config}" \
            --list $bed_arg $squeeze_arg "~{write_lines(gvcf)}" \
            | bcftools view - | $squeeze_cmd | bgzip -@ 4 -c > "~{output_name}.vcf.gz"
    >>>

    runtime {
        docker: "quay.io/mlin/glnexus:v1.2.2-7-g843c67e"
        cpu: cpu
        memory: "~{memoryGB}G"
        disks: "local-disk ~{diskGB} HDD"
    }

    output {
        File pvcf_gz = "${output_name}.vcf.gz"
    }
}
