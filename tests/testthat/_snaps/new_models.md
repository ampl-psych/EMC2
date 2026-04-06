# DDMGNG

    Code
      init_chains(emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                         1
      v        -1.23162378
      v_Sright -0.43601513
      a        -1.23465746
      Z         0.04221836
      t0       -0.79198362
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -149.5918
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_vector, dGNG, n_trials = 10)
    Output
         trials subjects     S     R        rt
      1       1        1  left right 0.5699457
      2       2        1  left right 1.0712502
      3       3        1  left right 0.9057347
      4       4        1  left right 0.8566880
      5       5        1  left right 0.6253342
      6       6        1  left right 1.0906341
      7       7        1  left right 0.5305895
      8       8        1  left right 0.5611776
      9       9        1  left right 0.5721548
      10     10        1  left right 0.7099338
      11     11        1 right right 0.5117904
      12     12        1 right right 0.5066511
      13     13        1 right right 0.7898381
      14     14        1 right  left        NA
      15     15        1 right right 0.9470190
      16     16        1 right  left        NA
      17     17        1 right  left        NA
      18     18        1 right  left        NA
      19     19        1 right right 0.4737746
      20     20        1 right  left        NA

# probit

    Code
      init_chains(emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                           1
      mean_Sleft  -0.2455619
      mean_Sright  0.5282036
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -11.19952
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_vector, dprobit, n_trials = 10)
    Output
         trials subjects     S     R rt
      1       1        1  left  left NA
      3       2        1  left  left NA
      5       3        1  left  left NA
      7       4        1  left right NA
      9       5        1  left  left NA
      11      6        1  left  left NA
      13      7        1  left  left NA
      15      8        1  left right NA
      17      9        1  left  left NA
      19     10        1  left  left NA
      21     11        1 right right NA
      23     12        1 right right NA
      25     13        1 right  left NA
      27     14        1 right right NA
      29     15        1 right right NA
      31     16        1 right right NA
      33     17        1 right right NA
      35     18        1 right  left NA
      37     19        1 right right NA
      39     20        1 right right NA

