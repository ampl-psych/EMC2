# run_factor

    Code
      LNR_factor[[1]]$samples$theta_lambda[, , idx]
    Output
      NULL

---

    Code
      LNR_factor[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8295289 -0.9638242
      m_lMd -0.2831784 -0.5197682
      s     -1.0057495 -0.5451441
      t0    -2.1806935 -1.6467459

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -1.0028059 -0.5132872 -0.5078411 -1.8032623 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                        m         m_lMd             s            t0
      m      2.744149e-02  1.769231e-14  1.189351e-14 -1.158945e-14
      m_lMd  1.769231e-14  6.802777e-02  2.881939e-02 -2.808263e-02
      s      1.189351e-14  2.881939e-02  1.180234e-01 -1.927401e-02
      t0    -1.158945e-14 -2.808263e-02 -1.927401e-02  4.215373e-02

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

