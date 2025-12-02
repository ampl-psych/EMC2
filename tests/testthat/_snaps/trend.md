# trend

    Code
      init_chains(LNR2cov, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                       1
      m       0.67512673
      m_lMd  -0.93638351
      s       0.36913678
      t0     -1.35711929
      m.w    -0.69011792
      m.d_ei  0.02412719
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -181.0193
      
      $idx
      [1] 1
      

# trend_shared

    Code
      init_chains(LNR2cov_shared, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                          1
      m           0.7365984
      m_lMd       1.0261898
      s           1.5774684
      t0         -1.7425644
      m.w        -0.8788452
      m.d_ei     -1.0144008
      m_lMd.w     0.7455618
      m_lMd.d_pd -0.3036067
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -20.76214
      
      $idx
      [1] 1
      

# premap trend works

    Code
      init_chains(LNR_premap, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                          1
      m          -1.5190914
      m_lMd       1.8129215
      s           0.4789562
      t0          0.6542443
      lMd.d1     -1.1511910
      lMd.d1_lR2  0.4772103
      m.w        -1.3819492
      m.d_ei     -0.1335546
      lMd.d2      1.7386078
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -460.517
      
      $idx
      [1] 1
      

# pretransform trend works

    Code
      init_chains(LNR_pretrans, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                        1
      m       -0.23862315
      m_lMd   -1.60968099
      s       -0.07479868
      t0      -0.62918350
      s.w     -0.46940990
      s.w_lR2  1.63913718
      m.w      0.31843950
      m.q0     2.58531781
      m.alpha -0.27882665
      s.d_ed   0.30316343
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -420.7097
      
      $idx
      [1] 1
      

# posttransform trend works

    Code
      init_chains(LNR_posttrans, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                        1
      m       -1.01882175
      m_lMd    2.29672893
      s        0.17880129
      t0       0.41984658
      s.w     -1.29129696
      s.w_lR2  1.68757362
      m.w      0.04419662
      m.d_pd  -1.65495809
      s.d_pi   0.68423590
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -460.517
      
      $idx
      [1] 1
      

# different trend base functions work

    Code
      init_chains(LNR_bases, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      m      -1.8337446
      m_lMd   1.4686824
      s       0.2724623
      t0     -1.8894146
      m.w    -2.0431978
      m.d_ei  1.1556991
      s.w    -0.7800323
      s.d_ed -0.2839388
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -166.7394
      
      $idx
      [1] 1
      

# polynomial trends work

    Code
      init_chains(LNR_poly, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                     1
      m     -1.1074205
      m_lMd -0.9857907
      s      0.7015643
      t0    -1.5078429
      m.d1   0.5197721
      m.d2   2.5500148
      m.d3  -0.3549074
      s.d1   0.2384021
      s.d2  -0.4849534
      s.d3  -1.7699108
      s.d4  -0.6364685
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -132.6838
      
      $idx
      [1] 1
      

# phase-specific trends work (premap, pretransform, posttransform)

    Code
      init_chains(LNR_phases, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                        1
      m       -0.07196167
      m_lMd    0.06587958
      s       -0.45122368
      t0      -2.20935035
      m.w     -1.29131004
      s.w      0.13365561
      s.d_ed  -0.29009503
      t0.w     0.22045064
      t0.d_pi  1.64918797
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -287.738
      
      $idx
      [1] 1
      

# par_input trend uses t0 as input to m trend

    Code
      init_chains(LNR_par_input, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                     1
      m     -1.6593198
      m_lMd  1.4048088
      s      3.3523195
      t0    -1.7750359
      m.w   -0.3235253
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -62.60216
      
      $idx
      [1] 1
      

# share works premap trend

    Code
      init_chains(LNR_shared_premap, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      m     -0.03723144
      m_lMd  0.88997977
      s     -2.26240896
      t0     0.58246707
      m.d2   0.22518542
      m.d3  -1.05376927
      s.d2  -0.01292946
      s.d3  -2.43473304
      s.d4   0.36702330
      shrd  -1.00831782
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -460.517
      
      $idx
      [1] 1
      

# share works posttransform trend

    Code
      init_chains(LNR_shared_posttransform, particles = 3, cores_per_chain = 1)[[1]]$
        samples
    Output
      $alpha
      , , 1
      
                      1
      m      0.27150933
      m_lMd  0.84470019
      s      0.05222749
      t0    -0.40192418
      m.d2   0.01704472
      m.d3  -1.27332214
      s.d2   0.34399709
      s.d3   0.60230341
      s.d4  -0.53128283
      shrd   0.35946573
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -436.8101
      
      $idx
      [1] 1
      

# trend_conditional

    Code
      attributes(make_data(p_vec, design_cond, n_trials = n_trials,
        conditional_on_data = FALSE, return_trialwise_parameters = TRUE))
    Message
      Imputing trial2 with random values
    Output
      $names
      [1] "trials"   "subjects" "S"        "R"        "trial2"   "rt"      
      
      $class
      [1] "data.frame"
      
      $row.names
       [1]  1  3  5  7  9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39
      
      $p_vector
           m m_lMd         s        t0 m.w m.q0    m.alpha
      1 -0.5     0 -1.203973 -1.609438 0.5    0 -0.8416212
      
      $trialwise_parameters
                     m   s  t0    m_trial2
       [1,] -0.5000000 0.3 0.2  0.00000000
       [2,] -0.5310232 0.3 0.2 -0.06204634
       [3,] -0.5558417 0.3 0.2 -0.11168340
       [4,] -0.7207814 0.3 0.2 -0.44156289
       [5,] -0.8527332 0.3 0.2 -0.70546648
       [6,] -0.7869319 0.3 0.2 -0.57386379
       [7,] -0.7342908 0.3 0.2 -0.46858164
       [8,] -0.4527176 0.3 0.2  0.09456474
       [9,] -0.2274591 0.3 0.2  0.54508183
      [10,] -0.3193801 0.3 0.2  0.36123978
      [11,] -0.3929169 0.3 0.2  0.21416614
      [12,] -0.3679439 0.3 0.2  0.26411213
      [13,] -0.3479655 0.3 0.2  0.30406892
      [14,] -0.4457397 0.3 0.2  0.10852068
      [15,] -0.5239590 0.3 0.2 -0.04791791
      [16,] -0.6879987 0.3 0.2 -0.37599739
      [17,] -0.8192305 0.3 0.2 -0.63846098
      [18,] -0.7815806 0.3 0.2 -0.56316113
      [19,] -0.7514606 0.3 0.2 -0.50292125
      [20,] -0.8051583 0.3 0.2 -0.61031658
      [21,] -0.8481164 0.3 0.2 -0.69623284
      [22,] -0.6298986 0.3 0.2 -0.25979727
      [23,] -0.4553244 0.3 0.2  0.08935118
      [24,] -0.6243666 0.3 0.2 -0.24873323
      [25,] -0.7596004 0.3 0.2 -0.51920076
      [26,] -0.6728986 0.3 0.2 -0.34579726
      [27,] -0.6035372 0.3 0.2 -0.20707446
      [28,] -0.5273251 0.3 0.2 -0.05465028
      [29,] -0.4663555 0.3 0.2  0.06728906
      [30,] -0.3391646 0.3 0.2  0.32167084
      [31,] -0.2374119 0.3 0.2  0.52517625
      [32,] -0.2973069 0.3 0.2  0.40538612
      [33,] -0.3452230 0.3 0.2  0.30955402
      [34,] -0.4089034 0.3 0.2  0.18219320
      [35,] -0.4598477 0.3 0.2  0.08030454
      [36,] -0.6396310 0.3 0.2 -0.27926193
      [37,] -0.7834575 0.3 0.2 -0.56691510
      [38,] -0.7735402 0.3 0.2 -0.54708049
      [39,] -0.7656064 0.3 0.2 -0.53121280
      [40,] -0.7973888 0.3 0.2 -0.59477754
      

# trend_multiple

    Code
      init_chains(LNR_multi, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                         1
      m          0.1595770
      m_lMd     -0.1311673
      s          1.0856025
      t0        -0.8015625
      m.w       -0.1148203
      m.d_ei    -0.1297780
      m.w.delta -1.2439554
      m.q0       0.1056588
      m.alpha   -1.3316796
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -371.9567
      
      $idx
      [1] 1
      

