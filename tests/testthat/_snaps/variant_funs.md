# run_factor

    Code
      LNR_factor[[1]]$samples$theta_lambda[, , idx]
    Output
                  [,1]      [,2]
      m      0.7908027 0.0000000
      m_lMd  3.0296415 0.6608763
      s     -1.4390380 1.0627351
      t0     0.8435664 1.0414702

---

    Code
      LNR_factor[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8151099 -0.8161671
      m_lMd -0.2378029 -0.5595275
      s     -0.9552971 -0.6504969
      t0    -2.1748102 -1.9464431

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.5169707  0.9304177 -0.6020754 -0.9499170 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                     m     m_lMd          s         t0
      m      0.6592341  2.395849 -1.1379950  0.6670945
      m_lMd  2.3958486  9.673562 -3.6574326  3.2439868
      s     -1.1379950 -3.657433  3.2804565 -0.1071171
      t0     0.6670945  3.243987 -0.1071171  1.8186457

# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9020895 -0.8695611
      m_lMd -0.3212180 -0.5936720
      s     -0.9453282 -0.5776115
      t0    -1.9817381 -1.7298590

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8303137 -0.6372378 -0.6819422 -1.8187986 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                     m     m_lMd         s         t0
      m     0.02535829 0.0000000 0.0000000 0.00000000
      m_lMd 0.00000000 0.1408963 0.0000000 0.00000000
      s     0.00000000 0.0000000 0.4420488 0.00000000
      t0    0.00000000 0.0000000 0.0000000 0.07646437

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.6168358 -0.5524390
      m_lMd -0.2882936 -0.3364735
      s     -1.0836746 -1.0307488
      t0    -3.2216830 -4.4020933

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
                m       m_lMd           s          t0 
      -0.47825755 -0.07236939 -0.91027040 -1.74705647 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                    m    m_lMd          s         t0
      m     0.2139126 0.000000  0.0000000  0.0000000
      m_lMd 0.0000000 0.404989  0.0000000  0.0000000
      s     0.0000000 0.000000  0.2229343 -0.1561151
      t0    0.0000000 0.000000 -0.1561151  3.1891020

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t        bd6t
      m     -0.5144656 -0.97423083
      m_lMd -1.6789264  0.04761233
      s     -0.2610040  0.38310109
      t0    -1.4550129 -1.20733178

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor),
      filter = "preburn", cores_for_props = 1)
    Output
               MD   wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single  512 0.000 2448    0 3568     0       1120  1329   363  209
      diag   -372 0.996 -364    1 -257     1        107  -471    14 -579
      factor -361 0.004 1563    0 2628     0       1065   498  -500 -568

