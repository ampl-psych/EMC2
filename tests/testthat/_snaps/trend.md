# trend

    Code
      init_chains(LNR2cov, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      m      -0.9685927
      m_lMd   0.7061091
      s       1.4890213
      t0     -1.8150926
      m.w     0.3304096
      m.d_ei -1.1421557
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -31.24098
      
      $idx
      [1] 1
      

# trend_shared

    Code
      init_chains(LNR2cov_shared, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                          1
      m          -0.9685927
      m_lMd       0.7061091
      s           1.4890213
      t0         -1.8150926
      m.w         0.3304096
      m.d_ei     -1.1421557
      m_lMd.w     0.1571934
      m_lMd.d_pd -2.0654072
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -31.03245
      
      $idx
      [1] 1
      

# premap trend works

    Code
      init_chains(LNR_premap, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                            1
      m            -0.9685927
      m_lMd         0.7061091
      s             1.4890213
      t0           -1.8150926
      m_lMd.d1      0.3304096
      m_lMd.d1_lR2 -1.1421557
      m.w           0.1571934
      m.d_ei       -2.0654072
      m_lMd.d2     -0.4405469
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -32.26135
      
      $idx
      [1] 1
      

# pretransform trend works

    Code
      init_chains(LNR_pretrans, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                        1
      m       -0.96859273
      m_lMd    0.70610908
      s        1.48902132
      t0      -1.81509255
      s.w      0.33040958
      s.w_lR2 -1.14215571
      m.w      0.15719342
      m.q0    -2.06540724
      m.alpha -0.44054688
      s.d_ed   0.00395328
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -65.4863
      
      $idx
      [1] 1
      

# posttransform trend works

    Code
      init_chains(LNR_posttrans, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                       1
      m       -0.9685927
      m_lMd    0.7061091
      s        1.4890213
      t0      -1.8150926
      s.w      0.3304096
      s.w_lR2 -1.1421557
      m.w      0.1571934
      m.d_pd  -2.0654072
      s.d_pi  -0.4405469
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -30.67787
      
      $idx
      [1] 1
      

# different trend base functions work

    Code
      init_chains(LNR_bases, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      m      -0.9685927
      m_lMd   0.7061091
      s       1.4890213
      t0     -1.8150926
      m.w     0.3304096
      m.d_ei -1.1421557
      s.w     0.1571934
      s.d_ed -2.0654072
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -34.8246
      
      $idx
      [1] 1
      

# polynomial trends work

    Code
      init_chains(LNR_poly, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      m     -0.96859273
      m_lMd  0.70610908
      s      1.48902132
      t0    -1.81509255
      m.d1   0.33040958
      m.d2  -1.14215571
      m.d3   0.15719342
      s.d1  -2.06540724
      s.d2  -0.44054688
      s.d3   0.00395328
      s.d4  -0.53711605
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -167.7708
      
      $idx
      [1] 1
      

# phase-specific trends work (premap, pretransform, posttransform)

    Code
      init_chains(LNR_phases, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                       1
      m       -0.9685927
      m_lMd    0.7061091
      s        1.4890213
      t0      -1.8150926
      m.w      0.3304096
      s.w     -1.1421557
      s.d_ed   0.1571934
      t0.w    -2.0654072
      t0.d_pi -0.4405469
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -280.2161
      
      $idx
      [1] 1
      

# par_input trend uses t0 as input to m trend

    Code
      init_chains(LNR_par_input, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                     1
      m     -0.9685927
      m_lMd  0.7061091
      s      1.4890213
      t0    -1.8150926
      m.w    0.3304096
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -32.93429
      
      $idx
      [1] 1
      

# share works premap trend

    Code
      init_chains(LNR_shared_premap, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      m     -0.96859273
      m_lMd  0.70610908
      s      1.48902132
      t0    -1.81509255
      m.d2   0.33040958
      m.d3  -1.14215571
      s.d2   0.15719342
      s.d3  -2.06540724
      s.d4  -0.44054688
      shrd   0.00395328
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -149.5046
      
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
      m     -0.96859273
      m_lMd  0.70610908
      s      1.48902132
      t0    -1.81509255
      m.d2   0.33040958
      m.d3  -1.14215571
      s.d2   0.15719342
      s.d3  -2.06540724
      s.d4  -0.44054688
      shrd   0.00395328
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -102.7324
      
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
                  m m_lMd   s  t0 m.w m.q0 m.alpha m.trial2.Qvalue subject trial
      1  -0.5000000     0 0.3 0.2 0.5    0     0.2      0.00000000       1     1
      2  -0.4669590     0 0.3 0.2 0.5    0     0.2      0.06608192       1     1
      3  -0.4405263     0 0.3 0.2 0.5    0     0.2      0.11894745       1     2
      4  -0.5666366     0 0.3 0.2 0.5    0     0.2     -0.13327318       1     2
      5  -0.6675248     0 0.3 0.2 0.5    0     0.2     -0.33504969       1     3
      6  -0.6183005     0 0.3 0.2 0.5    0     0.2     -0.23660107       1     3
      7  -0.5789211     0 0.3 0.2 0.5    0     0.2     -0.15784217       1     4
      8  -0.7696776     0 0.3 0.2 0.5    0     0.2     -0.53935518       1     4
      9  -0.9222828     0 0.3 0.2 0.5    0     0.2     -0.84456560       1     5
      10 -0.8818809     0 0.3 0.2 0.5    0     0.2     -0.76376185       1     5
      11 -0.8495594     0 0.3 0.2 0.5    0     0.2     -0.69911886       1     6
      12 -0.7792522     0 0.3 0.2 0.5    0     0.2     -0.55850443       1     6
      13 -0.7230064     0 0.3 0.2 0.5    0     0.2     -0.44601289       1     7
      14 -0.7321168     0 0.3 0.2 0.5    0     0.2     -0.46423352       1     7
      15 -0.7394050     0 0.3 0.2 0.5    0     0.2     -0.47881003       1     8
      16 -0.6927845     0 0.3 0.2 0.5    0     0.2     -0.38556893       1     8
      17 -0.6554880     0 0.3 0.2 0.5    0     0.2     -0.31097605       1     9
      18 -0.6067650     0 0.3 0.2 0.5    0     0.2     -0.21352996       1     9
      19 -0.5677865     0 0.3 0.2 0.5    0     0.2     -0.13557310       1    10
      20 -0.7125973     0 0.3 0.2 0.5    0     0.2     -0.42519453       1    10
      21 -0.8284458     0 0.3 0.2 0.5    0     0.2     -0.65689168       1    11
      22 -0.7159775     0 0.3 0.2 0.5    0     0.2     -0.43195498       1    11
      23 -0.6260028     0 0.3 0.2 0.5    0     0.2     -0.25200563       1    12
      24 -0.4813411     0 0.3 0.2 0.5    0     0.2      0.03731785       1    12
      25 -0.3656117     0 0.3 0.2 0.5    0     0.2      0.26877663       1    13
      26 -0.3150312     0 0.3 0.2 0.5    0     0.2      0.36993769       1    13
      27 -0.2745667     0 0.3 0.2 0.5    0     0.2      0.45086654       1    14
      28 -0.3109429     0 0.3 0.2 0.5    0     0.2      0.37811412       1    14
      29 -0.3400439     0 0.3 0.2 0.5    0     0.2      0.31991219       1    15
      30 -0.3251240     0 0.3 0.2 0.5    0     0.2      0.34975190       1    15
      31 -0.3131882     0 0.3 0.2 0.5    0     0.2      0.37362367       1    16
      32 -0.2285526     0 0.3 0.2 0.5    0     0.2      0.54289473       1    16
      33 -0.1608442     0 0.3 0.2 0.5    0     0.2      0.67831158       1    17
      34 -0.2552787     0 0.3 0.2 0.5    0     0.2      0.48944270       1    17
      35 -0.3308262     0 0.3 0.2 0.5    0     0.2      0.33834759       1    18
      36 -0.3325889     0 0.3 0.2 0.5    0     0.2      0.33482215       1    18
      37 -0.3339991     0 0.3 0.2 0.5    0     0.2      0.33200180       1    19
      38 -0.4419046     0 0.3 0.2 0.5    0     0.2      0.11619070       1    19
      39 -0.5282291     0 0.3 0.2 0.5    0     0.2     -0.05645818       1    20
      40 -0.5443814     0 0.3 0.2 0.5    0     0.2     -0.08876279       1    20
      

# trend_multiple

    Code
      init_chains(LNR_multi, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                         1
      m         -0.9685927
      m_lMd      0.7061091
      s          1.4890213
      t0        -1.8150926
      m.w        0.3304096
      m.d_ei    -1.1421557
      m.w.delta  0.1571934
      m.q0      -2.0654072
      m.alpha   -0.4405469
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -32.51614
      
      $idx
      [1] 1
      

# trend_covmap

    Code
      init_chains(LNR_covmap, particles = 3, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                        1
      m        -0.9685927
      m_lMd     0.7061091
      s         1.4890213
      t0       -1.8150926
      m.w_map1  0.3304096
      m.w_map2 -1.1421557
      m.q0      0.1571934
      m.alpha  -2.0654072
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -15.04154
      
      $idx
      [1] 1
      

