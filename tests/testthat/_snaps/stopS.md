# exG

    Code
      make_data(p_vector, designSSexG, n_trials = 10)
    Output
         trials subjects     S     R        rt
      1       1        1  left  left 0.4389288
      3       2        1  left  left 0.4607559
      5       3        1  left  left 1.0076918
      7       4        1  left  <NA>        NA
      9       5        1  left  left 0.4538172
      11      6        1  left  <NA>        NA
      13      7        1  left  left 0.5854384
      15      8        1  left  left 0.5193409
      17      9        1  left  left 0.5449355
      19     10        1  left  left 0.7572228
      21     11        1 right right 0.5042637
      23     12        1 right right 0.5649602
      25     13        1 right  <NA>        NA
      27     14        1 right  <NA>        NA
      29     15        1 right right 0.6147173
      31     16        1 right  left 0.6463509
      33     17        1 right right 0.6529462
      35     18        1 right right 0.6385749
      37     19        1 right  left 0.6280618
      39     20        1 right right 0.6454608

---

    Code
      init_chains(emc, particles = 10, cores_for_chains = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                          1
      mu        -1.91781543
      mu_lMTRUE -2.03490167
      sigma      0.14331586
      tau        0.01536863
      muS        0.42753549
      sigmaS     0.84325328
      tauS      -0.16526771
      gf        -0.06399883
      tf        -0.66880302
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -151.505
      
      $idx
      [1] 1
      

