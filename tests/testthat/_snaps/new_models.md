# DDMGNG

    Code
      init_chains(emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                        1
      v         0.7031960
      v_Sright  1.7179854
      a        -0.3277478
      Z         0.3962836
      t0       -0.5348966
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -187.1293
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_vector, dGNG, n_trials = 10)
    Output
         trials subjects     S     R        rt
      1       1        1  left right 0.7424707
      2       2        1  left right 0.9788611
      3       3        1  left right 0.7264890
      4       4        1  left right 0.6121781
      5       5        1  left right 0.4952376
      6       6        1  left right 0.6591371
      7       7        1  left right 0.5383684
      8       8        1  left right 0.6319298
      9       9        1  left right 0.4657021
      10     10        1  left right 0.4827050
      11     11        1 right right 0.5888895
      12     12        1 right  left        NA
      13     13        1 right  left        NA
      14     14        1 right  left        NA
      15     15        1 right  left        NA
      16     16        1 right right 0.7068593
      17     17        1 right right 0.5447677
      18     18        1 right  left        NA
      19     19        1 right  left        NA
      20     20        1 right right 0.9042361

