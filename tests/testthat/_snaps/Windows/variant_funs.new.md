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
      m     -0.8043590 -0.8742376
      m_lMd -0.2906632 -0.5648409
      s     -0.9885715 -0.6207059
      t0    -2.4729997 -1.8756629

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8011021 -0.3580180 -1.7864820 -0.8808069 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                        m         m_lMd             s            t0
      m      2.374857e-02  7.838620e-12 -8.550733e-12 -6.097784e-14
      m_lMd  7.838620e-12  2.897825e-01 -2.718986e-01 -1.939002e-03
      s     -8.550733e-12 -2.718986e-01  4.664389e-01 -2.000179e-01
      t0    -6.097784e-14 -1.939002e-03 -2.000179e-01  3.142439e-01

# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9191356 -1.0814948
      m_lMd -0.3330559 -0.7065455
      s     -0.8900481 -0.3831875
      t0    -2.0810165 -1.5745563

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
                m       m_lMd           s          t0 
      -0.98592856 -0.61764476  0.03152585 -1.89599243 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                   m      m_lMd         s        t0
      m     1.037155 0.00000000 0.0000000 0.0000000
      m_lMd 0.000000 0.01391776 0.0000000 0.0000000
      s     0.000000 0.00000000 0.1824788 0.0000000
      t0    0.000000 0.00000000 0.0000000 0.0814452

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8372551 -0.9715821
      m_lMd -0.2886091 -0.6724534
      s     -1.0111942 -0.4885203
      t0    -2.3428754 -1.7048039

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8936678 -0.4589840 -0.9086171 -1.3074862 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                    m      m_lMd         s        t0
      m     0.0166629 0.00000000 0.0000000 0.0000000
      m_lMd 0.0000000 0.02815875 0.0000000 0.0000000
      s     0.0000000 0.00000000 0.4840285 0.1140295
      t0    0.0000000 0.00000000 0.1140295 0.2475800

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8392820 -1.0584339
      m_lMd -0.2453826 -0.5590586
      s     -0.9814889 -0.3720923
      t0    -2.3615214 -1.6510395

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor,
        blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single  -522   0 -357    1 -221     1        136  -494  -630 -630
      diag    -412   0 -175    0   22     0        197  -372  -539 -568
      factor  -622   1  478    0 1025     0        547   -69  -606 -617
      blocked -523   0  -93    0  183     0        276  -368  -613 -644

