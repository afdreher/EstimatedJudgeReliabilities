These YAML files contain the options settings for running the examples provided
in the paper.

1) btl_scale_limit_100_unregularized.yaml
  - This file contains the settings for processing the decisions using the 
    vanilla Bradley-Terry-Luce form where scale values must lie on the interval 
    [-100, 100].

2) btl_scale_limit_100_regularized.yaml
  - This file contains the settings for processing the decisions using the
    Bradley-Terry-Luce form where scale values must lie on the interval 
    [-100, 100], but with a lambda penalty applied.  When lambda = 1, this is
    the Mease formulation.

1) btl_scale_limit_100_fonts.yaml
  - This file contains the settings for processing the decisions using the 
    vanilla Bradley-Terry-Luce form where scale values must lie on the interval 
    [-100, 100]. It has an output file for fonts already provided.

3) crowdbt_scale_limit_100_unregularized.yaml
  - This file contains the settings for processing the decisions using the 
    Chen et al.'s Crowd-BT form where scale values must lie on the interval 
    [-100, 100].

4) crowdbt_scale_limit_100_regularized.yaml
  - This file contains the settings for processing the decisions using the 
    Chen et al.'s Crowd-BT form where scale values must lie on the interval 
    [-100, 100] and where we have added a scaled version of the lambda 
    penalty from Mease. 

5) odonovan_scale_limit_100_and_mu_1_fonts.yaml
  - This file contains the settings for processing the decisions using the 
    O'Donovan et al.'s form where scale values must lie on the interval 
    [-100, 100], and the judge weight must lie on the interval [-1, 1].
    Results are written to odonovan_mu_1_results.csv.

6) odonovan_scale_limit_100_and_user_limit_6_146_fonts.yaml
  - This file contains the settings for processing the decisions using the 
    O'Donovan et al.'s form where scale values must lie on the interval 
    [-100, 100], and the judge weight must lie on the interval [-6.146, 6.146].
    Results are written to odonovan_mu_6_146_results.csv.

7) odonovan_scale_limit_100_and_mu_10_fonts.yaml
  - This file contains the settings for processing the decisions using the 
    O'Donovan et al.'s form where scale values must lie on the interval 
    [-100, 100], and the judge weight must lie on the interval [-10, 10].
    Results are written to odonovan_mu_10_results.csv.

8) odonovan_scale_limit_100_and_mu_50_fonts.yaml
  - This file contains the settings for processing the decisions using the 
    O'Donovan et al.'s form where scale values must lie on the interval 
    [-100, 100], and the judge weight must lie on the interval [-50, 50].
    Results are written to odonovan_mu_50_results.csv.

9) odonovan_scale_limit_100_and_user_limit_10_unregularized.yaml
  - This file contains the settings for processing the decisions using the 
    O'Donovan et al.'s form where scale values must lie on the interval 
    [-100, 100], and the judge weight must lie on the interval [-10, 10].

10) odonovan_scale_limit_500_and_user_limit_50_unregularized.yaml
  - This file contains the settings for processing the decisions using the 
    O'Donovan et al.'s form where scale values must lie on the interval 
    [-500, 500], and the judge weight must lie on the interval [-50, 50].

11) odonovan_scale_limit_100_and_user_limit_10_regularized.yaml
  - This file contains the settings for processing the decisions using the 
    O'Donovan et al.'s form where scale values must lie on the interval 
    [-100, 100], and the judge weight must lie on the interval [-10, 10], but to 
    which we have added a scaled version of the lambda penalty from Mease. 