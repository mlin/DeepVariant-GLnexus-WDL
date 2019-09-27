version 1.0

workflow fetch {
    # one 1000G sample per population (lexicographically first sample id in each)
    # table with sample and ERR IDs: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index
    Array[String] samples = [
        "ERR3240114/HG00096",
        "ERR3240160/HG00171",
        "ERR3241665/HG00403",
        "ERR3241804/HG00551",
        "ERR3242128/HG00759",
        "ERR3241828/HG01112",
        "ERR3241916/HG01500",
        "ERR3241984/HG01565",
        "ERR3243045/HG01583",
        "ERR3242062/HG01595",
        "ERR3242292/HG01879",
        "ERR3242357/HG02461",
        "ERR3242482/HG02922",
        "ERR3242642/HG03006",
        "ERR3242467/HG03052",
        "ERR3242669/HG03642",
        "ERR3242698/HG03713",
        "ERR3239458/NA06984",
        "ERR3239335/NA18486",
        "ERR3239646/NA18525",
        "ERR3239557/NA18939",
        "ERR3239683/NA19017",
        "ERR3239893/NA19625",
        "ERR3239894/NA19648",
        "ERR3239785/NA20502",
        "ERR3239997/NA20845"
    ]

    call samtools_slice_1000G_bam {
        input:
            ERR_slash_sample = samples,
            # ["ERR3239334/NA12878"],
            region_name = "ALDH2"
    }

    call curl as curl_ref_fasta {
        input:
            url = "https://s3.amazonaws.com/1000genomes/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    }

    output {
        File ref_fasta = curl_ref_fasta.file
        Array[File] bams = samtools_slice_1000G_bam.bams
    }
}

task curl {
    input {
        String url
    }

    command <<<
        apt-get update && apt-get install -y curl
        mkdir file
        cd file
        curl --retry 3 -Lss -O "~{url}"
    >>>

    output {
        File file = glob("file/*")[0]
    }

    runtime {
        docker: "ubuntu:disco"
    }
}

task samtools_slice_1000G_bam {
    input {
        Array[String] ERR_slash_sample
        String region = "chr12:111760000-111820000"
        String? region_name
    }

    command <<<
        set -euxo pipefail
        apt-get update && apt-get install -y samtools parallel

        cat << "EOF" > do1
        sample=$(basename "$1")
        outfn="$sample.~{select_first([region_name, region])}.bam"
        samtools view -@ 4 -b \
            -o "$outfn" \
            "https://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/data/$1.final.cram" \
            ~{region}
        echo "$outfn"
        EOF
        cat "~{write_lines(ERR_slash_sample)}" | \
            parallel -t -k -P 4 --retries 3 --halt 2 \
            bash -eux -o pipefail do1 {} > slice_filenames
    >>>

    output {
        Array[File] bams = read_lines("slice_filenames")
    }

    runtime {
        docker: "ubuntu:disco"
        cpu: 7
    }
}
