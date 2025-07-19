# DDMGNG

    Code
      init_chains(emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                        1
      v         0.7732297
      v_Sright  0.6825030
      a        -0.2179361
      Z        -0.6387820
      t0       -1.6570450
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -13.6954
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_vector, dGNG, n_trials = 10)
    Output
         trials subjects     S     R        rt
      1       1        1  left  left        NA
      2       2        1  left  left        NA
      3       3        1  left  left        NA
      4       4        1  left  left        NA
      5       5        1  left  left        NA
      6       6        1  left  left        NA
      7       7        1  left  left        NA
      8       8        1  left  left        NA
      9       9        1  left  left        NA
      10     10        1  left  left        NA
      11     11        1 right right 1.2853401
      12     12        1 right right 0.5337337
      13     13        1 right right 0.5287350
      14     14        1 right right 0.4991414
      15     15        1 right right 0.5803559
      16     16        1 right  left        NA
      17     17        1 right  left        NA
      18     18        1 right  left        NA
      19     19        1 right  left        NA
      20     20        1 right  left        NA

# probit

    Code
      init_chains(emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                           1
      mean_Sleft  -0.5126331
      mean_Sright  0.8401962
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -9.45741
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_vector, dprobit, n_trials = 10)
    Output
         trials subjects     S     R rt
      1       1        1  left  left NA
      3       2        1  left right NA
      5       3        1  left  left NA
      7       4        1  left  left NA
      9       5        1  left  left NA
      11      6        1  left  left NA
      13      7        1  left  left NA
      15      8        1  left right NA
      17      9        1  left  left NA
      19     10        1  left  left NA
      21     11        1 right right NA
      23     12        1 right right NA
      25     13        1 right right NA
      27     14        1 right right NA
      29     15        1 right right NA
      31     16        1 right right NA
      33     17        1 right  left NA
      35     18        1 right right NA
      37     19        1 right right NA
      39     20        1 right right NA

