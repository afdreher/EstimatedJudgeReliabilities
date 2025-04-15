# Estimated Judge Reliabilities


This GitHub repository is the artifact for our KDD 2024 paper.

> Estimated Judge Reliabilities for Weighted Bradley-Terry-Luce Are Not Reliable.  Andrew F. Dreher, Etienne Vouga, Donald S. Fussell.
> In KDD '24: The 30th ACM SIGKDD Conference on Knowledge Discovery and Data Mining Proceedings. https://doi.org/10.1145/3637528.3671907

## Overview

The included files contain the code and data necessary to run the examples in Creating Dictators (Section 4) and to run and generate the data in Simulation (Section 6).  

One can download the Fonts Attribute dataset from [O'Donovan et al.](https://www.dgp.toronto.edu/~donovan/font/index.html) to reproduce some of the Section 5 and Table 4 results.


## Installation

The examples in this repository were developed using the following software versions on OS X version 11.7.8 (20G1351) for x86. Simulations were performed using both Mac OS and [Ubuntu 20.04.6 LTS](https://releases.ubuntu.com/focal/)  (GNU/Linux 5.15.0-76-generic x86_64).

### Software

The following software versions were used for the examples.

[Mathematica](https://www.wolfram.com/mathematica/): 12.2.0.0

[Julia](https://julialang.org): 1.7.2

[Python](https://www.python.org) 3.11.4
  

#### Julia packages

| Hash        | Name               | Version     |
| ----------- | ------------------ | ----------- |
  [c7e460c6]  | ArgParse           | 1.1.4       |
  [69666777]  | Arrow              | 2.4.1       |
  [6e4b80f9]  | BenchmarkTools     | 1.3.2       |
  [336ed68f]  | CSV                | 0.8.5       |
  [a93c6f00]  | DataFrames         | 1.4.4       |
  [1313f7d8]  | DataFramesMeta     | 0.13.0      |
  [31c24e10]  | Distributions      | 0.25.79     |
  [59287772]  | Formatting         | 0.4.2       |
  [60bf3e95]  | GLPK               | 1.1.0       |
  [09f84164]  | HypothesisTests    | 0.10.11     |
  [7073ff75]  | IJulia             | 1.24.0      |
  [6deec6e2]  | IndexedTables      | 1.0.0       |
  [b6b21f68]  | Ipopt              | 1.1.0       |
  [c8e1da08]  | IterTools          | 1.4.0       |
  [682c06a0]  | JSON               | 0.21.3      |
  [4076af6c]  | JuMP               | 1.5.0       |
  [b4fcebef]  | Lasso              | 0.7.0       |
  [b8f27783]  | MathOptInterface   | 1.11.1      |
  [76087f3c]  | NLopt              | 0.6.5       |
  [bac558e1]  | OrderedCollections | 1.4.1       |
  [91a5bcdd]  | Plots              | 1.38.0      |
  [92933f4c]  | ProgressMeter      | 1.7.2       |
  [438e738f]  | PyCall             | 1.94.1      |
  [1a8c2f83]  | Query              | 1.0.0       |
  [276daf66]  | SpecialFunctions   | 2.1.7       |
  [2913bbd2]  | StatsBase          | 0.33.21     |
  [bd369af6]  | Tables             | 1.10.0      |
  [ddb6d928]  | YAML               | 0.4.8       |
  [56ddb016]  | Logging            |             |
  [de0858da]  | Printf             |             |
  [9a3f8284]  | Random             |             |
  [2f01184e]  | SparseArrays       |             |


#### Python packages

| Name               | Version     |
| ------------------ | ----------- |
| arrow              | 1.3.0       |
| pandas             | 2.2.1       |


## Mathematica Examples

In `Supplemental.nb`, we provide implementations of the functions described in the paper using Mathematica, paralleling Sections 3 and 4.  A few results are expanded beyond the paper versions to show in more detail.


## Research Code

The folder `julia/` contains the code, input, and settings files necessary to run the examples using our [Julia](https://julialang.org)-based implementation. We have provided additional details and notes in the README files in many folders.

- `/src`

    This folder contains the research code written in Julia.  Additional packages will need to be installed, as shown above. These have been tested on both Linux ([Ubuntu 20.04.6 LTS](https://releases.ubuntu.com/focal/) (GNU/Linux 5.4.0-153-generic x86_64)) and macOS Big Sur (macOS 11.7.8 (20G1351)).

    Note that the paper uses the term *JUDGES*, but the source files may refer to the decision makers as either *USERS* or *PARTICIPANTS*.  These terms were incorporated early in the development cycle and have not been updated to match the paper terminology. 

    There are a few different types of files:

    1. Method files

       - `crowdbt.jl` 
       
          This is our implementation of Chen et al.'s Crowd-BT algorithm.

        - `odonovan.jl` 
        
          This is our implementation of O'Donovan et al.'s algorithm.  Because O'Donovan et al.'s method is equivalent to Bradley-Terry-Luce when the weights are ignored, we use this same code for BTL.

          There are many different options available.  We used these during early experimentation and have left them in case others wish to run their own experiments with different parameters.

          The flag `--estimate-user-reliability` indicates that the code should use Eq. 5; without it, the code will use Eq. 1.

    2. Evaluation files

        - `nll_of_scale.jl` 
        
          This file contains the code to calculate the negative log-likelihood of a scale with an optional user weight file, which is used to compute Table 6 in section A.2 of the Appendix.

    3. Helper files

        - `data_loading.jl` 
        
          This file contains helper methods to load the data for both `crowdbt.jl` and `odonovan.jl`.

        - `negative_log.jl` 
        
          This file contains helper methods for calculating the sigmoids and negative log values.


- `/input`
    
    These contain the examples in the CSV format used by O'Donovan et al.

    NOTE: To reproduce some of the values from Section 5, one needs to create a new folder, called `fonts`, and place the file  `userChoices.csv` inside.

    The subfolders with input provided are:

    - `positive`

        This folder contains a single example that implements an example from Section 3.

    - `dictators`

        This folder contains the input necessary to reproduce the results in Section 4. The folder `symmetric` contains examples that verify and expand upon Sections 4.1 and 4.2, while `broken` contains the files for Section 4.3.

    - `simulation`

        This folder contains example input extracted from our simulation work (Section 6) and can be used to verify the reproduction of the input files, as described below.


- `/settings`

    This folder contains various experimental settings for running the examples or reproducing the results in the paper.  These files replace the command line options; one can use the command line instead of creating additional files.

- `/results`

    This folder contains the scales used for Section 5 (Fonts) and a limited number of zipped simulation results.  We provide the entire results for $\omega = 0.1$ and 8 judges per pair.  We can make additional zipped archives available if requested. We are not publishing the entire sweep here because the zipped files  use around 1.5 GB, but the different methods exhibited only small differences from the BTL order (Figure 5).

    Results for the simulation are provided for both "same" and "convenience" allocation methods.

- `/transcripts`

    This folder contains the output from IPOPT for most of the worked examples.  These can be used, for example, to compare live runs of the code against previous output. 


## Reproduction

To reproduce the experiments, use:

```
julia [path_to_source] [path_to_input] --experiment [path_to_settings]
```

### Core Instability (Section 3)

For example, to reproduce the positive example in Section 3 using O'Donovan et al.'s method, with a user limit of 10 and no lambda, use:

```
julia src/odonovan.jl input/positive/positive.csv --experiment settings/odonovan_scale_limit_100_and_user_limit_10_unregularized.yaml
```


### Creating Dictators (Sections 4.1 and 4.2)

The input files `input/dictators/symmetric/symmetric_4.csv` and `input/dictators/symmetric/symmetric_5.csv` correspond to the input pattern described in Section 4.  Additional input files with 6, 10, and 25 judges are provided.

The following two inputs demonstrate the change from preferring the dictatorial solution to the equality one in Section 4.2.

```
julia src/odonovan.jl input/dictators/symmetric/symmetric_4.csv --experiment settings/odonovan_scale_limit_100_and_user_limit_10_unregularized.yaml 

julia src/odonovan.jl input/dictators/symmetric/symmetric_5.csv --experiment settings/odonovan_scale_limit_100_and_user_limit_10_unregularized.yaml 
```

Running the files with Crowd-BT shows that the equality solution will not prevail.  The file with 25 judges shows that the gap between the minimum and maximum weights narrows but is not eliminated.

```
julia src/crowdbt.jl input/dictators/symmetric/symmetric_4.csv --experiment settings/crowdbt_scale_limit_100_unregularized.yaml 

julia src/crowdbt.jl input/dictators/symmetric/symmetric_5.csv --experiment settings/crowdbt_scale_limit_100_unregularized.yaml 
```

### Breaking the Symmetry (Section 4.3)

Although Table 1 was computed using Mathematica, the following commands can be used to reproduce the results using the Julia code.

```
julia src/crowdbt.jl input/dictators/broken/broken_5_extra_error_in_u1_to_agree_with_u4.csv --experiment settings/crowdbt_scale_limit_100_unregularized.yaml  

julia src/crowdbt.jl input/dictators/broken/broken_5_extra_error_in_u1_to_agree_with_u4.csv --experiment settings/crowdbt_scale_limit_100_regularized.yaml  

julia src/odonovan.jl input/dictators/broken/broken_5_extra_error_in_u1_to_agree_with_u4.csv --experiment settings/odonovan_scale_limit_100_and_user_limit_10_unregularized.yaml  

julia src/odonovan.jl input/dictators/broken/broken_5_extra_error_in_u1_to_agree_with_u4.csv --experiment settings/odonovan_scale_limit_100_and_user_limit_10_regularized.yaml  
```

The values may not match the Mathematica calculations exactly because the solution is parameter sensitive.  Starting with random values of $\mu_{\text{initial}}$ can yeild better negative log-likelihood scores.  In the early draft of the paper, Table 1 appears slightly different because we did not enforce $\mu_{\text{initial}} = 1$.

The Supplemental.pdf file contains the full table of N x (N-1) solutions.  These files are available as `broken_5_extra_error_in_u[N]_to_agree_with_u[M].csv`


### Fonts (Section 5)

One can recreate Table 4 using the `userChoices.csv` data provided by O'Donovan et al., which is available from their [project page](https://www.dgp.toronto.edu/~donovan/font/index.html).  The necessary file is available inside the Attribute Dataset.


```
julia src/odonovan.jl input/fonts/userChoices.csv --experiment [experiment_settings]
```

Use the following experiment settings:

- btl_scale_limit_100_fonts.yaml
- odonovan_scale_limit_100_and_mu_1_fonts.yaml
- odonovan_scale_limit_100_and_mu_6_146_fonts.yaml
- odonovan_scale_limit_100_and_mu_10_fonts.yaml
- odonovan_scale_limit_100_and_mu_50_fonts.yaml

These will generate files in the current working directory.  One can change the file name in the experiment_settings yaml files or use the corresponding command line flags.

Once the results are computed, use `nll_of_scale.jl` to compute the negative log-likelihoods.

Running the following two commands computes the BTL negative log-likelihood values (the far left column of Table 4)

```
julia src/odonovan.jl input/fonts/userChoices.csv --experiment settings/btl_scale_limit_100_fonts.yaml

julia src/nll_of_scale.jl input/fonts/userChoices.csv btl_results.csv
```

One can use the included Python script to compute the $\tau$ between the font results.  The command below computes $\tau$ between the BTL order and our result using $\mu^+ = 6.146$, which corresponds to $\tau (\clubsuit, \spadesuit)$ in Table 4.

```
python src/compute_tau_for_fonts.py results/fonts/btl/results.csv results/fonts/odonovan/original/mu_pm_6_146/scales.csv
```


### Simulation (Section 6)

Because of the size of the simulation data, we provide a limited number of sample input files. The provided code can generate the remaining files, which can be compared against the included sample to verify the recreation. When the space-saving [Arrow](https://arrow.apache.org) format is used, one should expect the simulation input data to use around 85 GB; the storage requirements will be much larger with CSV files. When saving the transcripts of the IPOPT output, one should expect between 3 and 5 GB for the output.

We performed the simulation using an [HTCondor](https://research.cs.wisc.edu/htcondor/htc.html) cluster with more than 1500 cores, mostly Intel Xenon (e.g., E5530, X5355, X5440).  While the computation time will depend on several factors, one should expect several CPU years to recreate the results shown in Section 6. The memory requirements, likewise, vary.  Crowd-BT with a large number of judges and pairs uses substantial resources; we occasionally had to re-run these inputs isolated on their own machine to avoid running out of system resources.

#### Recreating the input data

Recreating the input data consists of 3 primary steps. An optional fourth step enriches the analysis, but it is not necessary for running the simulation.

1. Create the decision data using `create_same_users.jl`. This generates the input for the "same" allocation.

    For example, the following command will generate files with 800 decisions each for $\omega = 0.1$ using the seed $713$.  It iterates over the number of judges per pair, so this will generate files for [8, 16, 24, 32] judges per pair.  One can add flags (e.g., "--users-per-pair") to restrict the numbers generated.  Additionally, one can leave off all of the flags to generate the full simulation set, although we caution that generating it in stages can be desirable in case unexpected interruptions occur.

    ```
    julia generate_same_users.jl --seed 713 --pair-counts 800 --width 10
    ```

    By default, the Julia script automatically creates a new folder `data` inside the code directory.

    The summary data provided for verification in the `input/simulation` folder contains the complete summary data for the given parameters. Unless one generates all of the seeds and pairs, the files will not match exactly, so quickly running `diff` may initially produce unexpected results.

2. Generate the "convenience" sample from the "same" sample to simulate crowdsourced judges.

      ```
      julia generate_convenience_from_same.jl --seeds 713 --pair-counts 800 --width 10
      ```

      This script creates a new folder containing the "convenience" samples derived from the "same" samples.  
      
      Warning: this script may take some time to process and requires the corresponding "same" sample to be created first. It will crash if even a single file is missing.

      The summaries for the "convenience" sample are not generated by default because the summary should be identical to those for the corresponding "same" sample since the votes do not change. However, one can compute the summary using `create_summaries.jl` to verify the result.

      ```
      julia create_summaries.jl --seeds 713 --pair-counts 800 --width 10 --series convenience
      ```
      

3. Run the simulation using `odonovan.jl` or `crowdbt.jl` with the desired settings.

4. (Optional) Create user summaries to better understand the vote distribution in the convenience sample.

   Generating summaries of the judges is mostly useful for the "convenience" sample since it gives the number of votes each judge makes as well as the disagreement from the s* order.  One can also create user summaries for the series using the "same" assignment method, but it is less interesting since every judge rates every pair.

    ```
    julia compute_user_summaries.jl  --seeds 713 --per-pair 8 --pair-counts 800 --width 10 --series convenience
    ```
  
    Please note: for this script, unlike the others, one needs to specify the per-pair, or it will default to just 8.  It does not iterate over all of the default per-pair options.


## Additional Results

In `Supplemental.pdf`, we show the full table created by adding a different error in every position for every judge (Section 4.3) and additional analysis of the Fonts results (Section 5).  These are provided for transparency and are not necessary to understand the paper's results.


## A remark about discontinuity and LLM alignment:

At KDD and afterward, several perceptive people asked how the paper results could impact the alignment (post-training) task with large language models. The answer is that the unreliability results shown in the paper do not depend on having a single, connected scale. The behavior is essentially unchanged. Thus, any system using weighted BTL as the core operation where the weights are estimated from the population's answers should have similar issues.

As an extreme example, consider the case where every pair embodies a scale with two stimuli. In this case, one stimulus will always be at the lower limit and the other at the positive limit. Which stimulus occupies which endpoint and what the values are for the estimated judge weights follow that of Section 3. The same set-coverage-like behavior dominates.
