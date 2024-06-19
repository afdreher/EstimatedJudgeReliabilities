We have provided input files in the format of O'Donovan et al.'s input that 
both correspond to the examples shown in the Mathematica notebook. Additional
examples are also included to show what happens for input that is larger than
can reasonably be shown in the Mathematica workbook.


===== FORMAT =====

All of the files are CSV files that adhere to the following format used by 
O'Donovan et al.:

  [ATTRIBUTE], [HIT_ID], [USER_ID], [STIMULUS_A], [STIMULUS_B], [DECISION]

  - The [ATTRIBUTE] field will always be "example", since we only use a single
    scale here.
  - The [HIT_ID] corresponds to the pair of stimuli being compared.  
  - The [USER_ID] corresponds to u[N].
  - The stimuli, [STIMULUS_A] and [STIMULUS_B], will be characters 
    corresponding to their position: 'A'=1, 'B'=2, and so forth.
  - The [DECISION] field is either 'more' or 'less' and means that STIMULUS_A 
    is either preferred to STIMULUS_B or vice versa, respectively.  It is 
    equivalent to having all of the decisions be either 'more' or 'less' and 
    swapping the stimulus positions, but we choose not to do so here.


===== FILES ====

There are three main directories of files:

positive/
  - This folder provides an example of case where the weighted methods identify
    the least reliable judge. It is not included in the paper for space, but
    is included here as a demonstration.  There is only the one file, 
    "positive.csv".

dictators/
  - This folder corresponds to the inputs used in Creating Dictators 
    (Section 4)

  symmetric/
    - This folder contains examples corresponding to symmetric judges

    symmetric_[N].csv
      - These files contain N judges with the basic pattern of 1 deviation from
        the expected order each.  The most important files are 4 and 5 since
        these are used to show the boundary between dictatorship and equality
        with O'Donovan et al.'s method.  Inputs 6, 10, 25 show that even for
        increasing numbers of judges, Chen et al.'s method prefers creating
        dictators.

  broken/
    - This folder contains examples of breaking the symmetry with different 
      changes to the preferences.

    broken_5_extra_error_in_u[A]_to_agree_with_u[B].csv
      - These files build on the basic pattern of symmetric_5.csv, but add in an
        extra deviation for user A to match user B's deviation.

  simulation/
    - This folder contains the sample data extracted from the simulation. It
      contains 1 seed (713) and 1 width (omega = 0.10) and for each of the 
      simulated number of judges (8, 16, 24, 32).

      The folder structure is [ALLOCATION]/[DISTANCE]/[JUDGE_COUNT]/[WIDTH]/[SEED]

      This example contains one distance method, even, which is the variety 
      described in the paper. The code references a method "expectation",
      which perturbs the distances, but renormalizes them so that in
      expectation, they are separated by 1 unit.

      The sample method can be "convenience" and "same", which are as described
      in the paper. The same sample means every judge participates in every
      pair. Convenience emulates crowdsourcing.

      Lastly, the files are named as:

        seed_[SEED]_count_[COUNT]_per_pair_[JUDGE_COUNT].arrow.lz4
        
      Here, all the seeds will be 713 and the JUDGE_COUNT will be one of 8, 16,
      24, or 32.  When running the generation scripts (see the src/ folder), 
      the seeds, counts, widths, etc. can all be different.

      The user_summary.csv files contain quick summaries of the judge votes.
      These show the number of times that a judge "disagrees" with the s*
      order and the total number of decisions that judge rendered.


==== FONTS ==== 

To run analysis of Fonts, one should download O'Donovan et al.'s fonts
attribute dataset from https://www.dgp.toronto.edu/~donovan/font/index.html