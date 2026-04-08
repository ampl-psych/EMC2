# CDM simulates and initializes

    Code
      cdm_data
    Output
        trials subjects          R        rt
      1      1        1 -0.3133472 1.0517062
      2      2        1  0.6478522 0.5303746
      3      3        1 -1.4022045 0.3531808
      4      4        1  0.7872438 0.9028722

---

    Code
      init_chains(cdm_emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                     1
      v      0.6730084
      theta  0.1647118
      a     -0.2286852
      t0    -0.4854718
      s     -0.3349777
      sv     0.3168479
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -47.61282
      
      $idx
      [1] 1
      

# SDM simulates and initializes

    Code
      sdm_data
    Output
        trials subjects         R        rt        R2
      1      1        1 0.8859462 0.6333035 2.6895510
      2      2        1 1.2606809 0.3453745 0.4989501
      3      3        1 1.1413576 0.5874518 2.2878413
      4      4        1 2.1124296 0.4434894 1.9875732

---

    Code
      init_chains(sdm_emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      v       0.2905103
      theta1  0.4004717
      theta2 -1.9158174
      a       0.5236415
      t0      2.6656195
      s       1.0530416
      sv     -1.4946516
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -92.1034
      
      $idx
      [1] 1
      

# HSDM simulates and initializes

    Code
      hsdm_data
    Output
        trials subjects        R        rt        R2         R3
      1      1        1 1.821304 0.3140800 1.7436259 -0.1250342
      2      2        1 1.672926 0.4395042 0.9739164  0.2621289
      3      3        1 2.110135 0.4217394 1.7526926  1.6678593
      4      4        1 1.280640 0.3519383 1.8401769  0.8536089

---

    Code
      init_chains(hsdm_emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      v       1.3253250
      theta1  0.2410887
      theta2 -0.2696486
      theta3 -1.3040563
      a      -1.4577938
      t0     -1.4830809
      s      -0.4383813
      sv     -2.1095788
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -38.49926
      
      $idx
      [1] 1
      

# PSDM simulates and initializes

    Code
      psdm_data
    Output
        trials subjects        R        rt
      1      1        1 1.130801 0.4864837
      2      2        1 1.415597 0.3112879
      3      3        1 2.157332 0.7243085
      4      4        1 1.840693 0.2807131

---

    Code
      init_chains(psdm_emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                       1
      v       1.09900070
      theta1  2.46578576
      a       0.01986150
      t0     -0.99358606
      s      -0.34020721
      sv      0.04050694
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -51.05512
      
      $idx
      [1] 1
      

# PHSDM simulates and initializes

    Code
      phsdm_data
    Output
        trials subjects        R        rt       R2
      1      1        1 2.185560 0.3386887 1.976002
      2      2        1 2.073026 0.4544274 1.273449
      3      3        1 1.329916 0.4246930 2.197339
      4      4        1 1.162111 0.6650170 2.742994

---

    Code
      init_chains(phsdm_emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                      1
      v       0.7677234
      theta1 -0.4616785
      theta2  0.5979320
      a      -0.7161771
      t0     -0.7645011
      s      -0.9509113
      sv     -1.7777507
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -67.40101
      
      $idx
      [1] 1
      

