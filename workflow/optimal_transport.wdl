workflow optimal_transport {
	Boolean run_validation = false
    Boolean run_ot = true

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
	Int? ot_disk_space = 150
	String wot_version = "1.0.2"

	if(run_ot) {

		call ot {
			input:
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
				#   	 no_overwrite
				out=out,

				num_cpu=num_cpu,
				memory=memory,
				preemptible=preemptible,
				zones=zones,
				disk_space=ot_disk_space,
				wot_version=wot_version
		}

    }

    if(run_validation) {
		call ot_validation {
			input:
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
				disk_space=ot_validation_disk_space,
				wot_version=wot_version
		}

	}
}

task create_day_pairs {
	File cell_days
	File? cell_filter
	Int preemptible
	String zones

    command {
        set -e

        python <<CODE
        import pandas as pd
        import numpy as np
        dtype = dict()
        cell_filter = '${cell_filter}'
        dtype['day'] = np.float64
        df = pd.read_table('${cell_days}', index_col='id', engine='python', sep=None, dtype=dtype)
        if cell_filter is not '':
			cell_ids = pd.read_table(cell_filter, index_col=0, header=None).index.values
            df = df[df.index.isin(cell_ids)]
        t = sorted(df['day'].unique())
        day_pairs = [(t[i], t[i + 1]) for i in range(len(t) - 1)]
        with open('day_pairs.txt', 'w') as w:
        	for p in day_pairs:
        		w.write(str(p[0]))
        		w.write(',')
        		w.write(str(p[1]))
        		w.write('\n')
        CODE

    }

    output {
    	Array[String] day_pairs = read_lines('day_pairs.txt')
    }

    runtime {
		docker: "continuumio/miniconda3:latest"
       	memory: "1 GB"
       	zones: "${zones}"
		bootDiskSizeGb: 10
		disks: "local-disk " + sub((1 + size(cell_days, 'GB')), "\\..*", "") + " HDD"
		cpu: 1
		preemptible: "${preemptible}"
    }
}

task create_day_triplets {
	File cell_days
	File? cell_filter
	Int preemptible
	String zones

    command {
        set -e

        python <<CODE
        import pandas as pd
        import numpy as np
        dtype = dict()
        cell_filter = '${cell_filter}'
        dtype['day'] = np.float64
        df = pd.read_table('${cell_days}', index_col='id', engine='python', sep=None, dtype=dtype)
        if cell_filter is not '':
			cell_ids = pd.read_table(cell_filter, index_col=0, header=None).index.values
            df = df[df.index.isin(cell_ids)]
        unique_times = sorted(df['day'].unique())

		with open('day_triplets.txt', 'w') as w:
			for i in range(len(unique_times) - 2):
				t0 = unique_times[i]
				t05 = unique_times[i + 1]
				t1 = unique_times[i + 2]
				w.write(str(t0))
				w.write(',')
				w.write(str(t05))
				w.write(',')
				w.write(str(t1))
				w.write('\n')
		CODE
    }

    output {
    	Array[String] triplets = read_lines('day_triplets.txt')
    }

    runtime {
		docker: "continuumio/miniconda3:latest"
		zones: "${zones}"
       	memory: "1 GB"
		bootDiskSizeGb: 10
		disks: "local-disk " + sub((1 + size(cell_days, 'GB')), "\\..*", "") + " HDD"
		cpu: 1
		preemptible: "${preemptible}"
    }
}

task ot {
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
	String? format
	#   	Boolean? no_overwrite
	String? out


	Int num_cpu
	String memory
	Int preemptible
	String zones
	Int disk_space
	String wot_version

    command {
        set -e
        wot optimal_transport \
        --matrix ${matrix} \
        --cell_days ${cell_days} \
		${"--cell_growth_rates " + cell_growth_rates} \
		${"--parameters " + parameters} \
		${true="--transpose" false="" transpose} \
		${"--local_pca " + local_pca} \
		${"--growth_iters " + growth_iters} \
		${"--gene_filter " + gene_filter} \
		${"--cell_filter " + cell_filter} \
		${"--cell_day_filter " + cell_day_filter} \
		${"--scaling_iter " + scaling_iter} \
		${"--inner_iter_max " + inner_iter_max} \
		${"--epsilon " + epsilon} \
		${"--lambda1 " + lambda1} \
		${"--lambda2 " + lambda2} \
		${"--max_iter " + max_iter} \
		${"--batch_size " + batch_size} \
		${"--tolerance " + tolerance} \
		${"--epsilon0 " + epsilon0} \
		${"--tau " + tau} \
		${"--ncells " + ncells} \
		${"--ncounts " + ncounts} \
		${"--cell_days_field " + cell_days_field} \
		${"--cell_growth_rates_field " + cell_growth_rates_field} \
		${true="--verbose" false="" verbose} \
		${"--format " + format} \
		${"--out " + out}
    }

    output {
    	Array[File] transport_maps = glob("*.${format}")
    }

    runtime {
		docker: "regevlab/wot-${wot_version}"
       	memory: "${memory}"
       	zones: "${zones}"
		bootDiskSizeGb: 10
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: "${preemptible}"
    }
}


task ot_validation {
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
	String? format
#   	Boolean? no_overwrite
	String? out

	# validation specific parameters
	File? day_triplets
	Int? interp_size
	File? covariate
	String? covariate_field
	Boolean? full_distances


	Int num_cpu
	String memory
	Int preemptible
	String zones
	Int disk_space
	String wot_version

    command {
        set -e

        wot optimal_transport_validation \
		--matrix ${matrix} \
		--cell_days ${cell_days} \
		${"--cell_growth_rates " + cell_growth_rates} \
		${"--parameters " + parameters} \
		${true="--transpose" false="" transpose} \
		${"--local_pca " + local_pca} \
		${"--growth_iters " + growth_iters} \
		${"--gene_filter " + gene_filter} \
		${"--cell_filter " + cell_filter} \
		${"--cell_day_filter " + cell_day_filter} \
		${"--scaling_iter " + scaling_iter} \
		${"--inner_iter_max " + inner_iter_max} \
		${"--epsilon " + epsilon} \
		${"--lambda1 " + lambda1} \
		${"--lambda2 " + lambda2} \
		${"--max_iter " + max_iter} \
		${"--batch_size " + batch_size} \
		${"--tolerance " + tolerance} \
		${"--epsilon0 " + epsilon0} \
		${"--tau " + tau} \
		${"--ncells " + ncells} \
		${"--ncounts " + ncounts} \
		${"--cell_days_field " + cell_days_field} \
		${"--cell_growth_rates_field " + cell_growth_rates_field} \
		${true="--verbose" false="" verbose} \
		${"--format " + format} \
		${"--out " + out}

		${"--day_triplets " + day_triplets} \
		${"--interp_size " + interp_size} \
		${"--covariate " + covariate} \
		${"--covariate_field " + covariate_field} \
		${true="--full_distances" false="" full_distances}
    }

    output {
    	File summary =  "${out}_validation_summary.txt"
    	Array[File] cv_outputs =  glob("${out}_cv_validation_summary*")
    	Array[File] full_outputs =  glob("${out}_full_validation_summary*")
    }

    runtime {
		docker: "regevlab/wot-${wot_version}"
		zones: "${zones}"
       	memory: "${memory}"
		bootDiskSizeGb: 10
		disks: "local-disk ${disk_space} HDD"
		cpu: "${num_cpu}"
		preemptible: "${preemptible}"
    }
}





