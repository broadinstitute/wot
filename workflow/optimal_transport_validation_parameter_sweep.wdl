import "https://api.firecloud.org/ga4gh/v1/tools/regev:optimal_transport/versions/11/plain-WDL/descriptor" as wot


workflow optimal_transport_validation_parameter_sweep {
    File matrix
    File cell_days
    Boolean? transpose = false
    File? cell_growth_rates
    File? gene_filter
    File? cell_filter
    String? cell_day_filter
    String? format = "loom"
    String? out = "wot"
    Int? num_cpu = 64
    String? memory = "57.6GB"
    Int? preemptible = 2
    File? covariate

    Int? ot_validation_disk_space = 150
    Boolean? full_distances = true
    # first column contains parameter name, 2nd column contains comma separated list of values
    File parameter_file
	String? zones = "us-east1-d us-west1-a us-west1-b"

    call create_inputs {
        input:
            parameter_file=parameter_file,
            preemptible=preemptible,
            zones=zones
    }

    scatter (p in create_inputs.parameters) {
        call wot.optimal_transport {
            input:
                matrix=matrix,
                cell_days=cell_days,
                transpose=transpose,
                cell_growth_rates=cell_growth_rates,
                gene_filter=gene_filter,
                cell_filter=cell_filter,
                cell_day_filter=cell_day_filter,
                format = format,
                out = out,
                num_cpu = num_cpu,
                memory = memory,
                preemptible = preemptible,
                covariate=covariate,
                full_distances=full_distances,
                run_validation=true,
                run_ot=false,
                parameters=p,
                zones=zones,
                ot_validation_disk_space=ot_validation_disk_space
        }
    }
 }


task create_inputs {
    # first column contains parameter name, 2nd column contains comma separated list of values
    File parameter_file
    Int preemptible
    String zones

    command {
        set -e

        python <<CODE

        import itertools

        parameter_names = []
        parameter_values = []
        with open('${parameter_file}', 'rt') as f:
            for line in f:
                tokens = line.strip().split('\t')
                parameter_names.append(tokens[0].strip())
                value_tokens = tokens[1].split(',')
                parameter_values.append(value_tokens)

        parameter_values_product = list(itertools.product(*parameter_values))

        for i in range(len(parameter_values_product)):
            with open('parameters-' + str(i) + '.txt', 'w') as w:
                params = parameter_values_product[i]
                for j in range(len(params)):
                    w.write(parameter_names[j])
                    w.write('\t')
                    w.write(str(params[j]))
                    w.write('\n')

        CODE
    }

    output {
        Array[File] parameters = glob('parameters-*.txt')
    }

    runtime {
        docker: "continuumio/miniconda3:latest"
        memory: "1 GB"
        zones: "${zones}"
        bootDiskSizeGb: 10
        disks: "local-disk " + sub((1 + size(parameter_file, 'GB')), "\\..*", "") + " HDD"
        cpu: 1
        preemptible: "${preemptible}"
    }
}
