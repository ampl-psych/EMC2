# run_factor

    Code
      LNR_factor[[1]]$samples$theta_lambda[, , idx]
    Output
                    F1          F2
      m      0.4042440 0.000000000
      m_lMd -0.4318224 0.379480463
      s      0.6085597 0.101328100
      t0     0.1192645 0.006457266

---

    Code
      LNR_factor[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9538521 -0.9353932
      m_lMd -0.3830433 -0.6787215
      s     -0.8030984 -0.4866918
      t0    -1.7952026 -1.6199259

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.5966560 -1.0288748 -0.2052301 -1.6511029 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                      m       m_lMd           s          t0
      m      0.19259767 -0.17456160  0.24600658  0.04821194
      m_lMd -0.17456160  0.35581441 -0.22433767 -0.04905066
      s      0.24600658 -0.22433767  0.45313739  0.07323385
      t0     0.04821194 -0.04905066  0.07323385  0.03879326

# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8746630 -0.8803158
      m_lMd -0.3238564 -0.5286320
      s     -0.9305809 -0.5595418
      t0    -2.0076934 -1.7339672

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8781100 -0.4785104 -0.6912904 -1.9212876 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd          s         t0
      m     0.001665839 0.00000000 0.00000000 0.00000000
      m_lMd 0.000000000 0.01593168 0.00000000 0.00000000
      s     0.000000000 0.00000000 0.03541001 0.00000000
      t0    0.000000000 0.00000000 0.00000000 0.02246014

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9277165 -0.8420786
      m_lMd -0.3919777 -0.5142153
      s     -0.8024551 -0.6875488
      t0    -1.8154290 -1.8660141

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8107844 -0.4890125 -0.5082068 -1.7735883 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd          s         t0
      m     0.003418992 0.00000000 0.00000000 0.00000000
      m_lMd 0.000000000 0.03850907 0.00000000 0.00000000
      s     0.000000000 0.00000000 0.05772351 0.03310210
      t0    0.000000000 0.00000000 0.03310210 0.05980684

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8379934 -0.9074886
      m_lMd -0.4239260 -0.5827497
      s     -0.9029937 -0.6085757
      t0    -2.1009790 -1.7736820

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor,
        blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single  -495   0 -417    1 -360     1         57  -474  -531 -531
      diag    -503   0 -349    0 -230     0        119  -468  -568 -586
      factor  -423   0 2584    0 4080     0       1496  1089  -392 -407
      blocked -520   1 -201    0  -10     0        191  -392  -568 -582

