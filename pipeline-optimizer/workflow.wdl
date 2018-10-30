workflow ExampleMetricsWorkflow {
    Float variable
    Float constant

    call ExampleMetricsTask {
        input:
            variable = variable,
	    constant = constant
    }

    output {
        Float objective_value = ExampleMetricsTask.objective_value
    }
}

task ExampleMetricsTask {
    Float variable
    Float constant

    command <<<
        set -e
        python2 -c "print ${variable} * ${variable} - 2 * ${variable} + ${constant}"
    >>>

    runtime {
        docker: "gatk-evaluation/pipeline-optimizer:test"
        memory: "2000 MB"
        disks: "local-disk 100 HDD"
        preemptible: 3
        cpu: 1
    }

    output {
        Float objective_value = read_float(stdout())
    }
}

