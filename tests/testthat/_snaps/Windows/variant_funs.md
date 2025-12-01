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
      m     -0.9160650 -1.0464409
      m_lMd -0.3688805 -0.6324767
      s     -0.8697713 -0.3858096
      t0    -2.0085773 -1.5940102

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9395881 -0.3717887 -0.4903179 -2.0308361 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                        m         m_lMd             s            t0
      m      4.704915e-02 -1.349783e-11  1.733364e-11 -2.469458e-11
      m_lMd -1.349783e-11  2.201783e-02 -4.920711e-03  7.010347e-03
      s      1.733364e-11 -4.920711e-03  1.167946e-01  3.195008e-02
      t0    -2.469458e-11  7.010347e-03  3.195008e-02  1.423015e-01

# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9470492 -1.0637712
      m_lMd -0.3317728 -0.6817155
      s     -0.8732219 -0.4168069
      t0    -1.9612423 -1.6056028

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9728177 -0.6768029 -0.6949762 -1.6464152 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m     m_lMd        s        t0
      m     0.004275191 0.0000000 0.000000 0.0000000
      m_lMd 0.000000000 0.1365253 0.000000 0.0000000
      s     0.000000000 0.0000000 0.116714 0.0000000
      t0    0.000000000 0.0000000 0.000000 0.1347646

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8228285 -0.9371903
      m_lMd -0.2817573 -0.6832554
      s     -1.0462369 -0.5569486
      t0    -2.4475033 -1.7908100

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -1.1214937 -0.5197986 -0.4495438 -2.0889890 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                     m     m_lMd           s          t0
      m     0.07375605 0.0000000  0.00000000  0.00000000
      m_lMd 0.00000000 0.0893964  0.00000000  0.00000000
      s     0.00000000 0.0000000  0.17229799 -0.07513489
      t0    0.00000000 0.0000000 -0.07513489  0.26119822

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8211001 -0.7845790
      m_lMd -0.3076343 -0.4735092
      s     -1.0085658 -0.7406569
      t0    -2.3460615 -2.3117218

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor,
        blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single  -525   0 -383    0 -262 0.004        121  -505  -626 -626
      diag    -497   0 -283    0 -100 0.000        184  -467  -628 -651
      factor  -627   1 -109    0  124 0.000        234  -343  -570 -577
      blocked -553   0 -403    1 -273 0.996        130  -534  -642 -664

