# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8378073 -0.9872863
      m_lMd -0.2876366 -0.5898393
      s     -1.0249672 -0.6047458
      t0    -2.3783442 -1.7427491

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9743222 -0.6134234 -0.7679532 -2.2539005 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                     m      m_lMd          s       t0
      m     0.02531788 0.00000000 0.00000000 0.000000
      m_lMd 0.00000000 0.09894308 0.00000000 0.000000
      s     0.00000000 0.00000000 0.03923812 0.000000
      t0    0.00000000 0.00000000 0.00000000 0.106451

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9722393 -0.9902785
      m_lMd -0.3535471 -0.6376234
      s     -0.8799747 -0.5240860
      t0    -1.9649581 -1.7118019

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9862506 -0.7445712 -0.1982309 -2.1982936 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                       m      m_lMd          s         t0
      m     0.0002440168 0.00000000  0.0000000  0.0000000
      m_lMd 0.0000000000 0.05788441  0.0000000  0.0000000
      s     0.0000000000 0.00000000  0.2532656 -0.1998045
      t0    0.0000000000 0.00000000 -0.1998045  0.1911231

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.0126139 -0.5952710
      m_lMd -0.6268864 -0.5299123
      s     -0.5883318 -0.9124103
      t0    -1.7063384 -2.8790704

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, blocked = LNR_blocked),
      stage = "preburn", cores_for_props = 1)
    Output
                MD   wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single  -243 0.000  248    0  649     0        401  -153  -555 -555
      diag    -535 0.055 1066    0 1916     0        851   215  -607 -635
      blocked -541 0.945 -144    1   79     1        223  -367  -564 -590

