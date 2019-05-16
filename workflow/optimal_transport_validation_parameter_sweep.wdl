import "https://api.firecloud.org/ga4gh/v1/tools/wot:optimal_transport/versions/1/plain-WDL/descriptor" as wot


workflow optimal_transport_validation_parameter_sweep {
    # first column contains parameter name, 2nd column contains comma separated list of values
    File parameter_file

	File matrix
	File cell_days
	File? cell_growth_rates
	File? parameters
	#	File? config
	Boolean? transpose
	Int? local_pca
	Int? growth_iters
	File? gene_filter
	File? cell_filter
	String? cell_day_filter
	Int? scaling_iter
	Int? inner_iter_max
	Float? epsilon
	Float? lambda1
	Float? lambda2
	Int? max_iter
	Int? batch_size
	Int? tolerance
	Float? epsilon0
	Float? tau
	Int? ncells
	Int? ncounts
	String? solver
	String? cell_days_field
	String? cell_growth_rates_field
	Boolean? verbose
	String? format = "h5ad"
	#   	Boolean? no_overwrite
	String? out = "wot"

    # validation parameters
   	File? day_triplets
   	Int? interp_size
   	File? covariate
   	String? covariate_field
   	Boolean? full_distances

   	Int? num_cpu = 2
   	String? memory = "52GB"
   	Int? preemptible = 2
   	String? zones = "us-east1-d us-west1-a us-west1-b"
   	Int? ot_validation_disk_space = 150


    call create_inputs {
        input:
            parameter_file=parameter_file,
            preemptible=preemptible,
            zones=zones
    }

    scatter (p in create_inputs.parameters) {
        call wot.optimal_transport {
            input:
            	parameters=p,
               	matrix=matrix,
				cell_days=cell_days,
				cell_growth_rates=cell_growth_rates,
				parameters=parameters,
				#	 config
				transpose=transpose,
				local_pca=local_pca,
				growth_iters=growth_iters,
				gene_filter=gene_filter,
				cell_filter=cell_filter,
				cell_day_filter=cell_day_filter,
				scaling_iter=scaling_iter,
				inner_iter_max=inner_iter_max,
				epsilon=epsilon,
				lambda1=lambda1,
				lambda2=lambda2,
				max_iter=max_iter,
				batch_size=batch_size,
				tolerance=tolerance,
				epsilon0=epsilon0,
				tau=tau,
				ncells=ncells,
				ncounts=ncounts,
				solver=solver,
				cell_days_field=cell_days_field,
				cell_growth_rates_field=cell_growth_rates_field,
				verbose=verbose,
				format=format,
				out=out,

				day_triplets=day_triplets,
				interp_size=interp_size,
				covariate=covariate,
				covariate_field=covariate_field,
				full_distances=full_distances,

				num_cpu=num_cpu,
				memory=memory,
				preemptible=preemptible,
				zones=zones,
				disk_space=ot_validation_disk_space


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
