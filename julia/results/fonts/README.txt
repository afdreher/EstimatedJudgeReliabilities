These folders contain the results from running the weighted methods on
O'Donovan et al.'s Fonts dataset.

Inside each folder are 3 files:

  [judges.csv] - This is the resulting judge weight file (mu_i values for
      O'Donovan; eta_i values for Crowd-BT)

  [scales.csv] - This is the resulting latent scale file (sigma_i values)

  [transcript.txt] - This shows the IPOPT solver transcript


For the primary folders:

  - [crowdbt] contains the results from running src/crowdbt.jl. These results
        are computed by treating each attribute as  an independent scale, which
        is the behavior from Chen et al.'s work.  The parameters used are:

    ===== PARAMETERS FOR CROWD-BT =====

      lambda = 0
      scale = [-100, 100]
      judge starting value = 0.9
      solver tolerance = 1E-6


  - [odonovan] contains the results from src/odonovan.jl with various limits of
        mu_i.  There are two folders:
   
      - [independent] Contains a run of src/odonovan.jl with 
          -6.146 <= mu <= 6.146 where each attribute is as an independent
          scale such that the mu_i values are not shared. This mimics the
          Crowd-BT formulation, but differs from the version in 
          O'Donovan et al.

      - [original] These are the results used to compute the NLL scores used
          in Table 6, shown in section A.2. This version has the shared judge
          weights as described in O'Donovan et al.  

          Use the Julia script src/nll_of_scale.jl to find the values shown in
          the table.  For example, the following command should produce the 
          mu <= |10| results: 

            julia src/nll_of_scale.jl input/fonts/userChoices.csv results/fonts/mu_pm_10/scales.csv --judges results/fonts/mu_pm_10/judges.csv

          There are 4 different limits on mu_i:

            - [mu_pm_1] means -1 <= mu <= 1

            - [mu_pm_6_146] means -6.146 <= mu <= 6.146

            - [mu_pm_10] means -10 <= mu <= 10

            - [mu_pm_50] means -50 <= mu <= 50


        For all, the same parameters are used:

        ===== PARAMETERS FOR O'DONOVAN =====

          lambda = 0
          scale = [-100, 100]
          judge starting value = 1.0
          solver tolerance = 1E-6