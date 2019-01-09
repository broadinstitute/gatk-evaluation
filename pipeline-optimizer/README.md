This directory contains a simple example of a pipeline that can be optimized using Advisor and Cromwell.  The ingredients are:

1) A WDL workflow (`workflow.wdl`) that returns the objective function as one of the outputs (`Float objective_value`).  The example here (`ExampleMetricsWorkflow`) calculates a simple quadratic function `x^2 - 2x + c` of the variable `x` and the constant `c`.

2) A JSON template (`template.json`) that specifies the parameters of the WDL that will be held constant over the optimization.  In this case, we fix `c = ExampleMetricsWorkflow.constant = 1.0`.

3) A JSON (`scan.json`) specifying the parameters over which to optimize.  In this case, we scan over `x = ExampleMetricsWorkflow.variable` in the range `[-10, 10]`.  We also specify other parameters for the optimization, such as whether to minimize or maximize the objective function, the maximum number of trials for the run, and the number of random trials with which to initialize our optimization algorithm.

4) A pipeline-optimizer Docker (`us.gcr.io/broad-dsde-methods/gatk-evaluation:pipeline-optimizer`) with the Advisor and Cromwell APIs installed.  Running the command
    
    `docker run --network host -it us.gcr.io/broad-dsde-methods/gatk-evaluation:pipeline-optimizer`

    locally will spin up an Advisor server at `http://127.0.0.1:8001` to which you can submit studies and trials.

5) A python script (`launch_study.py`) for generating Advisor trials, submitting jobs to a Cromwell server, and waiting for returned objective functions.  See comments in the script for details about behavior upon encountering Cromwell exceptions or failed runs.

6) A bash wrapper script (`launch_study.sh`) that starts the optimization run.  An example invocation might be:
    
   `./launch_study.sh -d us.gcr.io/broad-dsde-methods/gatk-evaluation:pipeline-optimizer -a http://127.0.0.1:8001 -c http://cromwell-v34.dsde-methods.broadinstitute.org -n ExampleStudy -l BayesianOptimization -w workflow.wdl -t template.json -s scan.json`
   
   This launches a study named `ExampleStudy` that optimizes the cost function using the `BayesianOptimization` algorithm.  The `launch_study.py` script should be relatively robust to connection interruptions and workflow failures, but note that the Advisor framework also allows for studies to be stopped and restarted.

Users should feel free to modify these basic ingredients as appropriate.  
