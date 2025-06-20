# exG

    Code
      make_data(p_vector, designSSexG, n_trials = 10, functions = list(SSD = mySSD_function))
    Output
         trials subjects     S     R  SSD        rt
      1       1        1  left  left  Inf 0.7679116
      3       2        1  left  left  Inf 0.6162440
      5       3        1  left  left  Inf 0.6845308
      7       4        1  left  left  Inf 0.6260193
      9       5        1  left  left  Inf 0.5996385
      11      6        1  left right  Inf 0.6858717
      13      7        1  left right  Inf 0.6631416
      15      8        1  left  left  Inf 0.6100059
      17      9        1  left  left  Inf 0.5582757
      19     10        1  left  left  Inf 0.5187092
      21     11        1 right right  Inf 0.4792897
      23     12        1 right right 0.26 0.4665539
      25     13        1 right right 0.35 0.5336244
      27     14        1 right right  Inf 0.4747087
      29     15        1 right right  Inf 0.6470309
      31     16        1 right right  Inf 0.5723853
      33     17        1 right right  Inf 0.6253725
      35     18        1 right right  Inf 0.5772106
      37     19        1 right right 0.46 0.5168540
      39     20        1 right right  Inf 0.5263251

---

    Code
      init_chains(emc, particles = 10, cores_for_chains = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                          1
      mu        -0.90210035
      mu_lMTRUE  0.25671609
      sigma     -1.56090286
      tau       -1.63296389
      muS       -0.05540344
      sigmaS     1.27858701
      tauS       2.24313992
      gf        -0.42248920
      tf        -0.18398258
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -105.9188
      
      $idx
      [1] 1
      

