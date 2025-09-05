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
      1 -49.95395
      
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
      1 -194.7257
      
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
      1 -186.3582
      
      $idx
      [1] 1
      

# pretransform trend works

    Code
      init_chains(LNR_pretrans, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                         1
      m        -0.02764523
      m_lMd     0.56222037
      s        -1.04587870
      t0       -2.57222702
      s.w      -0.29067254
      s.w_lR2   2.13028966
      m.w       0.07270155
      m.q0     -0.51286198
      m.alpha  -1.05489862
      s.d_ed   -0.39551142
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -417.3681
      
      $idx
      [1] 1
      

# posttransform trend works

    Code
      init_chains(LNR_posttrans, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                        1
      m        -0.8858629
      m_lMd    -1.1231660
      s         0.7205207
      t0       -1.4218538
      s.w       0.1719652
      s.w_lR2   0.4364415
      m.w      -0.3847409
      m.d_pd   -0.1073068
      s.d_pi   -0.3799550
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -193.8656
      
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
      1 -365.9303
      
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
      1 -358.1376
      
      $idx
      [1] 1
      

# share works premap trend

    Code
      init_chains(LNR_shared_premap, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      m     -3.04036083
      m_lMd -0.03862094
      s      1.05595193
      t0    -1.35353617
      shrd  -0.16071726
      m.d2   1.43433607
      m.d3   0.29033594
      s.d2  -0.15329406
      s.d3   0.13855736
      s.d4  -1.01216385
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -167.8266
      
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
      m      0.23662324
      m_lMd  0.10949662
      s      0.18978810
      t0    -2.35492413
      shrd   0.49788000
      m.d2  -0.99634927
      m.d3   0.08045588
      s.d2   0.07067752
      s.d3  -1.02057618
      s.d4   0.52312880
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -36.2619
      
      $idx
      [1] 1
      

