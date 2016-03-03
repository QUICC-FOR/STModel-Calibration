### Some notes on how the sampler works

The first phase is adaptation. This usually progresses pretty quickly (less than 10,000 iterations), but it's not necessarily a huge problem if it takes longer. FYI, on my laptop running on 2 cores, it took 18 hours and 3500 iterations to complete adaptation in a test run of the intercept only model (so no climate signal -- much faster than the full model). Running on 20+ cores on the cluster should be at least 10x as fast, but sometimes the sampler gets very slow for reasons I don't fully understand. All of this is logged, so you can inspect the log files to see how things are going (you can see from the screenshot `logging_example.png` what this looks like on my machine). Adaptation will complete when the acceptance rate for all parameters is between around 0.15 and 0.35.

The second phase is sampling. This will continue until the number of iterations you asked for is complete or you hit the wall time (no worries about this, you can resume the sampler). This will log every 500 iterations (the first step is shown in the screenshot).

### Starting the jobs
You will be running 8 total jobs; 4 chains for the full model and 4 chains for the intercept only model. You should take care that each job runs on a different node, and it's best to wait a few seconds before starting each one (to avoid RNG problems due to clock synchronization between the nodes). To get it going, use the following steps:

1. Unzip the mcmc_4state directory. This should be **outside** the git repo for the model code itself. The instructions below are going to assume you have everything in your home directory, so the model code will be in `~/STModel-MCMC` and the run files will be in `~/mcmc_4state`.
2. Copy the binary from the source code directory to the run directory, renaming the binary in the process so that we have a good chain of reproducability (we know what version of the sampler we used):

        cp ~/STModel-MCMC/bin/stm4_mcmc ~/mcmc_4state/bin/stm4_mcmc1.5
3. Use your job sumission system to launch the samplers. See the one-liners in `start_mcmc` as examples. **Each script in `start_mcmc` should be a different job on a different node.** I'll explain a bit what everything means so you can change them if needed. The first one looks like this:
 
         ~/mcmc_4state/bin/stm4_mcmc1.5 -d -p dat/mcmc_inits1.txt -t dat/mcmc_calib.txt -o res/ch1 -i 200000 -c 20 -l 5 -v 2 2>res/ch1/log.txt &

There is just one option you might/will have to change, depending on the configuration of the cluster.

* `-c 20` is the number of cores to use. Each job has to be on its own node, so this can be at maximum the number of cores in one node.


Most of the options you don't need to change:

* `-d` tells the sampler to also record DIC. Leave this alone.
* `-p dat/mcmc_inits1.txt` gives the location of the inits file. This is chain 1, so we use mcmc_inits1.txt. There are 4 of these for each of the 4 chains in the full model, plus one `mcmc_inits_int.txt` for the intercept-only model (all 4 chains use the same inits for that model).
* `-t dat/mcmc_calib.txt` is the calibration data formatted for the sampler. Its the same for all 8 runs.
* `-o res/ch1` is the output directory. `ch1` is chain 1, `ch1_int` is chain 1 for the intercept only model, etc. It's already set up correctly, but just remember, **these have to be different for each chain, or files will be overwritten.**
* `-i 200000` the number of iterations to run. Leave it like this for now; if we need more we will do it by resuming the samplers.
* `-l 5` The number of years between transitions. This is fixed by the data; leave it alone.
* `-v 2` The number specifies how much logging to do. No real reason to change this unless you want to do some troubleshooting. Anything higher than 2 produces a lot of text, so be ready.
* `2>res/ch1/log.txt &` Directs the logging into a file and runs the program in the background.

### Resuming jobs
You have resume scripts in `resume_mcmc`. These look a lot like the start scripts, but with slightly different options. However, you will have to do some customization to get the right number of iterations. When you are ready to resume a job (for example, the job in `ch1`--just change the name to do the others), do the following:

1. check the length of the posterior to see how many iterations completed. The easiest way is to use `wc`:

        wc -l ~/mcmc_4state/res/ch1/posterior.csv       
This number will be one larger than the number of iterations (for example, you might get a value of `10001`, meaning 10000 iterations have been done).
2. Use this value to change the `-i 200000` flag in the resume script for that chain. For for our example, with chain 1, we want 200.000 iterations and have completed 10.000, so open up `mcmc_resume/resume1.sh` and change `-i 200000` to `-i 190000`.
3. All the other options can stay the same (you will see that they are much simpler for resuming -- most of the information is saved in `resumeData.txt`. All you have to do is launch the resume script from your job scheduler and it will automatically start adding new rows to the `posterior.csv` file.

You can resume the scripts as many times as you need. Just make sure the old job is finished before starting a new one with the same chain!