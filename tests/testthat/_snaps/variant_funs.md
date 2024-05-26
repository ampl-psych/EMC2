# run_factor

    Code
      LNR_factor[[1]]$samples$theta_lambda[, , idx]
    Output
                  [,1]       [,2]
      m      0.8819422 0.00000000
      m_lMd  1.2451270 0.35839878
      s      0.7179031 0.26796354
      t0    -0.8450668 0.09326284

---

    Code
      LNR_factor[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.7661381 -1.0144493
      m_lMd -0.3214950 -0.6941293
      s     -0.8375086 -0.4176126
      t0    -2.2546042 -1.5120190

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
              m     m_lMd         s        t0 
      -1.654374 -1.227006 -1.072185 -1.433432 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                     m      m_lMd          s         t0
      m      0.8224309  1.0981301  0.6331491 -0.7453001
      m_lMd  1.0981301  1.7007127  0.9899184 -1.0187902
      s      0.6331491  0.9899184  0.7882633 -0.5816850
      t0    -0.7453001 -1.0187902 -0.5816850  0.7811450

# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8918293 -1.0112046
      m_lMd -0.3032468 -0.6554517
      s     -0.9048326 -0.4559385
      t0    -2.0056110 -1.5367747

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8838214 -0.6478713 -0.6885587 -1.9623001 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                     m      m_lMd          s        t0
      m     0.01634551 0.00000000 0.00000000 0.0000000
      m_lMd 0.00000000 0.01606372 0.00000000 0.0000000
      s     0.00000000 0.00000000 0.03936784 0.0000000
      t0    0.00000000 0.00000000 0.00000000 0.3949877

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8200062 -0.9867921
      m_lMd -0.3466458 -0.5887587
      s     -0.9256529 -0.5087412
      t0    -2.0937707 -1.5618025

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.3975078 -0.6297665 -0.9090594 -2.2586168 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                   m    m_lMd         s        t0
      m     1.053811 0.000000 0.0000000 0.0000000
      m_lMd 0.000000 2.539353 0.0000000 0.0000000
      s     0.000000 0.000000 0.3434857 0.3333133
      t0    0.000000 0.000000 0.3333133 1.5763158

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9700347 -1.1965591
      m_lMd -0.2607282 -0.4144615
      s     -0.9501442 -0.3004955
      t0    -1.7970237 -1.4804571

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor),
      filter = "preburn", cores_for_props = 1)
    Output
               MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -245   0   91    0  374     0        283  -193   -82 -476
      diag   -484   1 -354    1 -238     1        116  -470  -571 -586
      factor -307   0  118    0  458     0        340  -223  -491 -563

