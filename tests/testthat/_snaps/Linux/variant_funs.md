# run_factor

    Code
      LNR_factor[[1]]$samples$theta_lambda[, , idx]
    Output
                    F1         F2
      m      0.3871466  0.0000000
      m_lMd -0.2639166  0.4511165
      s     -0.1509166 -0.1631589
      t0     0.1987681 -0.5779386

---

    Code
      LNR_factor[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8617966 -0.9699954
      m_lMd -0.3472919 -0.6453849
      s     -0.9796450 -0.6128228
      t0    -2.1063684 -1.5453802

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.4538775 -0.7512564 -0.4593961 -1.0559969 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd           s          t0
      m      0.25332810 -0.1021744 -0.05842685  0.07695239
      m_lMd -0.10217439  0.3060136 -0.03377430 -0.31317584
      s     -0.05842685 -0.0337743  0.14837942  0.06429844
      t0     0.07695239 -0.3131758  0.06429844  0.43902223

# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8905103 -0.9378623
      m_lMd -0.2973458 -0.6303210
      s     -0.9700930 -0.5016096
      t0    -1.9440951 -1.5285791

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9254360 -0.5878139 -0.8392987 -1.8542595 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                    m      m_lMd          s         t0
      m     0.0023799 0.00000000 0.00000000 0.00000000
      m_lMd 0.0000000 0.03864044 0.00000000 0.00000000
      s     0.0000000 0.00000000 0.05257501 0.00000000
      t0    0.0000000 0.00000000 0.00000000 0.01781726

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.7703327 -0.9199820
      m_lMd -0.2865068 -0.6444077
      s     -1.1229272 -0.5667622
      t0    -2.4548518 -1.6286011

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -1.2238423 -0.3295895 -0.2308937 -0.2395329 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                  m     m_lMd         s        t0
      m     1.60083 0.0000000 0.0000000 0.0000000
      m_lMd 0.00000 0.1998595 0.0000000 0.0000000
      s     0.00000 0.0000000 0.1699233 0.5070567
      t0    0.00000 0.0000000 0.5070567 1.8337766

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.0396038 -0.6343280
      m_lMd -0.1583645 -0.3631572
      s     -0.8743594 -0.8204171
      t0    -1.9011606 -2.1694764

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor,
        blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD   wMD  DIC  wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single  -312 0.000  -34 0.000  135 0.000        169  -202  -371 -371
      diag    -496 0.003 -407 0.997 -319 0.996         88  -495  -563 -583
      factor  -508 0.997 -247 0.000 -107 0.000        140  -387  -515 -527
      blocked -468 0.000 -395 0.003 -307 0.004         88  -483  -550 -571

