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
      1 -34.14854
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_vector, dGNG, n_trials = 10)
    Output
         subjects     S trials     R        rt
      1         1  left      1  left        NA
      2         1 right      1  left        NA
      3         1  left      2  left        NA
      4         1 right      2 right 0.7486671
      5         1  left      3  left        NA
      6         1 right      3 right 1.0350038
      7         1  left      4  left        NA
      8         1 right      4 right 0.4898115
      9         1  left      5  left        NA
      10        1 right      5 right 0.7117752
      11        1  left      6  left        NA
      12        1 right      6  left        NA
      13        1  left      7  left        NA
      14        1 right      7 right 1.2646380
      15        1  left      8  left        NA
      16        1 right      8 right 0.4848544
      17        1  left      9  left        NA
      18        1 right      9 right 0.4545908
      19        1  left     10  left        NA
      20        1 right     10 right 0.5352639

# probit

    Code
      init_chains(emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                           1
      mean_Sleft  -1.8345288
      mean_Sright  0.5068369
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -9.816787
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_vector, dprobit, n_trials = 10)
    Output
         subjects     S trials     R rt
      1         1  left      1  left NA
      3         1 right      1 right NA
      5         1  left      2 right NA
      7         1 right      2  left NA
      9         1  left      3 right NA
      11        1 right      3  left NA
      13        1  left      4  left NA
      15        1 right      4 right NA
      17        1  left      5 right NA
      19        1 right      5 right NA
      21        1  left      6  left NA
      23        1 right      6  left NA
      25        1  left      7  left NA
      27        1 right      7 right NA
      29        1  left      8  left NA
      31        1 right      8  left NA
      33        1  left      9  left NA
      35        1 right      9 right NA
      37        1  left     10  left NA
      39        1 right     10 right NA

