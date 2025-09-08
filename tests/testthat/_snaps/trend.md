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
      m.B0    1.2153540
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
      m.B0        1.38416964
      m.d_ei     -0.49761490
      m_lMd.B0    0.05966031
      m_lMd.d_pd  0.83571136
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -192.7699
      
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
      m.B0       -0.013016616
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
      m        -0.02764523
      m_lMd     0.56222037
      s        -1.04587870
      t0       -2.57222702
      s.B0     -0.29067254
      s.B0_lR2  2.13028966
      m.B0      0.07270155
      m.q0     -0.51286198
      m.alpha  -1.05489862
      s.d_ed   -0.39551142
      
      
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
      m        -0.8858629
      m_lMd    -1.1231660
      s         0.7205207
      t0       -1.4218538
      s.B0      0.1719652
      s.B0_lR2  0.4364415
      m.B0     -0.3847409
      m.d_pd   -0.1073068
      s.d_pi   -0.3799550
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -33.59888
      
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
      m.B0    1.32060356
      m.d_ei -0.19936759
      s.B0    1.02552088
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
      m.B0     0.8567414
      s.B0    -0.6980340
      s.d_ed   1.6384842
      t0.B0   -0.2191840
      t0.d_pi -0.8801203
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -225.4447
      
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
      m.B0  -0.64113399
      
      
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
      shrd  -0.8065226
      m.d2   0.8246624
      m.d3   0.4571500
      s.d2   1.5296289
      s.d3   1.0298274
      s.d4  -0.4164044
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -393.5949
      
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
      shrd  -0.30681825
      m.d2   0.01042179
      m.d3   0.70392462
      s.d2   0.20343661
      s.d3   1.20844566
      s.d4  -1.32706234
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -337.8635
      
      $idx
      [1] 1
      

