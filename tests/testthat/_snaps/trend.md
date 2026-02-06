# trend

    Code
      init_chains(LNR2cov, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      m      -2.7444171
      m_lMd   0.5476119
      s       1.9769066
      t0     -2.0011088
      m.w     1.2153540
      m.d_ei -0.2860677
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -48.81096
      
      $idx
      [1] 1
      

# trend_shared

    Code
      init_chains(LNR2cov_shared, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                           1
      m           3.06176217
      m_lMd      -1.29148010
      s           1.71412324
      t0         -1.18214962
      m.w         1.38416964
      m.d_ei     -0.49761490
      m_lMd.w     0.05966031
      m_lMd.d_pd  0.83571136
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -192.4246
      
      $idx
      [1] 1
      

# premap trend works

    Code
      init_chains(LNR_premap, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                            1
      m           0.138841321
      m_lMd       0.597885432
      s           0.795426900
      t0         -1.182857399
      lMd.d1      0.427554507
      lMd.d1_lR2  0.004832919
      m.w        -0.013016616
      m.d_ei      0.996688307
      lMd.d2      0.803537699
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -185.5338
      
      $idx
      [1] 1
      

# pretransform trend works

    Code
      init_chains(LNR_pretrans, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                        1
      m       -0.02764523
      m_lMd    0.56222037
      s       -1.04587870
      t0      -2.57222702
      s.w     -0.29067254
      s.w_lR2  2.13028966
      m.w      0.07270155
      m.q0    -0.51286198
      m.alpha -1.05489862
      s.d_ed  -0.39551142
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -121.3379
      
      $idx
      [1] 1
      

# posttransform trend works

    Code
      init_chains(LNR_posttrans, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                       1
      m       -0.8858629
      m_lMd   -1.1231660
      s        0.7205207
      t0      -1.4218538
      s.w      0.1719652
      s.w_lR2  0.4364415
      m.w     -0.3847409
      m.d_pd  -0.1073068
      s.d_pi  -0.3799550
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -190.2867
      
      $idx
      [1] 1
      

# different trend base functions work

    Code
      init_chains(LNR_bases, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                       1
      m       0.09801042
      m_lMd   1.36758960
      s       1.28261336
      t0     -0.96209016
      m.w     1.32060356
      m.d_ei -0.19936759
      s.w     1.02552088
      s.d_ed -0.71163103
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -365.708
      
      $idx
      [1] 1
      

# polynomial trends work

    Code
      init_chains(LNR_poly, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                     1
      m     -0.2507031
      m_lMd  0.6991970
      s      0.1287192
      t0    -0.9348798
      m.d1   1.9555753
      m.d2  -2.3501194
      m.d3  -1.1840015
      s.d1   0.2815442
      s.d2   0.2861484
      s.d3   0.4585218
      s.d4   0.5578078
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -358.0446
      
      $idx
      [1] 1
      

# phase-specific trends work (premap, pretransform, posttransform)

    Code
      init_chains(LNR_phases, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                       1
      m       -0.3927859
      m_lMd    0.5309645
      s       -0.3249141
      t0      -1.9579093
      m.w      0.8567414
      s.w     -0.6980340
      s.d_ed   1.6384842
      t0.w    -0.2191840
      t0.d_pi -0.8801203
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -289.8456
      
      $idx
      [1] 1
      

# par_input trend uses t0 as input to m trend

    Code
      init_chains(LNR_par_input, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      m      1.26745820
      m_lMd -1.30562857
      s     -0.03611131
      t0    -0.59117615
      m.w   -0.64113399
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -419.1373
      
      $idx
      [1] 1
      

# share works premap trend

    Code
      init_chains(LNR_shared_premap, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                     1
      m     -0.3998853
      m_lMd -0.4771568
      s      0.4612157
      t0    -0.8711165
      m.d2  -0.8065226
      m.d3   0.8246624
      s.d2   0.4571500
      s.d3   1.5296289
      s.d4   1.0298274
      shrd  -0.4164044
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -387.7981
      
      $idx
      [1] 1
      

# share works posttransform trend

    Code
      init_chains(LNR_shared_posttransform, particles = 10, cores_per_chain = 1)[[1]]$
        samples
    Output
      $alpha
      , , 1
      
                      1
      m      1.25661377
      m_lMd -0.77424945
      s      0.11498310
      t0    -1.13396388
      m.d2  -0.30681825
      m.d3   0.01042179
      s.d2   0.70392462
      s.d3   0.20343661
      s.d4   1.20844566
      shrd  -1.32706234
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -336.1057
      
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
       [2,] -0.3292343 0.3 0.2  0.34153135
       [3,] -0.1926218 0.3 0.2  0.61475643
       [4,] -0.2611014 0.3 0.2  0.47779724
       [5,] -0.3158851 0.3 0.2  0.36822988
       [6,] -0.5090088 0.3 0.2 -0.01801753
       [7,] -0.6635077 0.3 0.2 -0.32701545
       [8,] -0.4735738 0.3 0.2  0.05285234
       [9,] -0.3216267 0.3 0.2  0.35674657
      [10,] -0.4052861 0.3 0.2  0.18942781
      [11,] -0.4722136 0.3 0.2  0.05557281
      [12,] -0.3985552 0.3 0.2  0.20288953
      [13,] -0.3396285 0.3 0.2  0.32074291
      [14,] -0.3828055 0.3 0.2  0.23438900
      [15,] -0.4173471 0.3 0.2  0.16530586
      [16,] -0.5662962 0.3 0.2 -0.13259233
      [17,] -0.6854554 0.3 0.2 -0.37091088
      [18,] -0.7516159 0.3 0.2 -0.50323173
      [19,] -0.8045442 0.3 0.2 -0.60908841
      [20,] -0.6854462 0.3 0.2 -0.37089233
      [21,] -0.5901677 0.3 0.2 -0.18033547
      [22,] -0.5586274 0.3 0.2 -0.11725475
      [23,] -0.5333951 0.3 0.2 -0.06679017
      [24,] -0.6457328 0.3 0.2 -0.29146570
      [25,] -0.7356031 0.3 0.2 -0.47120611
      [26,] -0.5838372 0.3 0.2 -0.16767447
      [27,] -0.4624246 0.3 0.2  0.07515085
      [28,] -0.6102719 0.3 0.2 -0.22054383
      [29,] -0.7285498 0.3 0.2 -0.45709958
      [30,] -0.8146464 0.3 0.2 -0.62929277
      [31,] -0.8835237 0.3 0.2 -0.76704731
      [32,] -0.8359223 0.3 0.2 -0.67184470
      [33,] -0.7978413 0.3 0.2 -0.59568260
      [34,] -0.5917613 0.3 0.2 -0.18352259
      [35,] -0.4268973 0.3 0.2  0.14620542
      [36,] -0.5130261 0.3 0.2 -0.02605215
      [37,] -0.5819291 0.3 0.2 -0.16385822
      [38,] -0.7302152 0.3 0.2 -0.46043034
      [39,] -0.8488440 0.3 0.2 -0.69768803
      [40,] -0.8350868 0.3 0.2 -0.67017356
      

# trend_multiple

    Code
      init_chains(LNR_multi, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                         1
      m         -1.9716301
      m_lMd     -0.7343946
      s          0.2149761
      t0        -1.1131262
      m.w       -0.9148022
      m.d_ei     1.2214175
      m.w.delta -0.5846742
      m.q0       0.2086729
      m.alpha   -0.4205229
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -382.9031
      
      $idx
      [1] 1
      

