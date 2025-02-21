# trend

    Code
      init_chains(LNR2cov, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      m      -1.1979009
      m_lMd   1.0470309
      s      -0.9388293
      t0     -0.9865591
      m.B0    0.8184259
      m.d_ei -0.7822744
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -365.4276
      
      $idx
      [1] 1
      

# trend_shared

    Code
      init_chains(LNR2cov_shared, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                          1
      m           1.0738577
      m_lMd      -2.2657574
      s           1.2122577
      t0         -1.7407933
      m.B0       -1.8507349
      m.d_ei     -0.6956514
      m_lMd.B0    2.6699944
      m_lMd.d_pd  0.2087916
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -18.85132
      
      $idx
      [1] 1
      

# premap trend works

    Code
      init_chains(LNR_premap, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                          1
      m          -0.6235522
      m_lMd       1.0999322
      s           0.2105472
      t0         -2.1280977
      lMd.d1     -0.4869198
      lMd.d1_lR2 -0.9021004
      m.B0        0.2567161
      m.d_ei     -1.5609029
      lMd.d2     -1.6329639
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -1.691976
      
      $idx
      [1] 1
      

# pretransform trend works

    Code
      init_chains(LNR_pretrans, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                        1
      m         0.3340588
      m_lMd     0.5141329
      s         1.7853760
      t0        0.1481788
      s.B0      0.1730763
      s.B0_lR2 -0.1974613
      m.B0      2.4566301
      m.q0     -0.3629998
      m.alpha   0.2239200
      s.d_ed    0.2086578
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -432.8909
      
      $idx
      [1] 1
      

# posttransform trend works

    Code
      init_chains(LNR_posttrans, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                         1
      m        -0.88157151
      m_lMd     0.02317448
      s        -0.56063160
      t0       -0.33560090
      s.B0     -0.02423430
      s.B0_lR2 -0.03752685
      m.B0     -0.58364760
      m.d_pd    0.65977460
      s.d_pi    0.28870637
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -439.7668
      
      $idx
      [1] 1
      

# different trend base functions work

    Code
      init_chains(LNR_bases, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                        1
      m      -0.931458191
      m_lMd  -0.742610532
      s       0.860133296
      t0     -1.028127504
      m.B0    0.865253657
      m.d_ei -0.002536259
      s.B0    0.931652459
      s.d_ed -1.429944707
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -307.9806
      
      $idx
      [1] 1
      

# polynomial trends work

    Code
      init_chains(LNR_poly, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      m     -0.54913508
      m_lMd  0.03106616
      s      0.19307917
      t0    -1.14773228
      m.d1  -2.02743395
      m.d2   0.79001905
      m.d3   1.01595861
      s.d1  -0.73848808
      s.d2  -1.95925906
      s.d3  -0.03969343
      s.d4   1.45342098
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -208.2988
      
      $idx
      [1] 1
      

# share works premap trend

    Code
      init_chains(LNR_shared_premap, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                     1
      m     -0.3996138
      m_lMd  2.5265480
      s     -0.7102682
      t0    -2.0213681
      shrd  -0.1659889
      m.d2   0.3035319
      m.d3   2.2508028
      s.d2   1.0873114
      s.d3  -0.5137314
      s.d4  -1.0852663
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -5.104573
      
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
      m     -0.7060564
      m_lMd  1.1075357
      s      0.2666760
      t0    -1.0135051
      shrd   0.2408547
      m.d2  -1.0385293
      m.d3  -0.4287139
      s.d2   0.2257823
      s.d3  -0.1467424
      s.d4   0.2546146
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -301.6662
      
      $idx
      [1] 1
      

