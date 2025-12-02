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
      m     -0.8479736 -0.8926952
      m_lMd -0.3226738 -0.4977847
      s     -1.0147345 -0.6747540
      t0    -2.3176046 -1.9300261

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9780874 -0.5403146 -0.5657009 -1.9942174 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                        m         m_lMd             s            t0
      m      2.648496e-02  8.175122e-14  7.057745e-14 -3.281909e-14
      m_lMd  8.175122e-14  5.888972e-02  2.957135e-02 -1.375092e-02
      s      7.057745e-14  2.957135e-02  1.135393e-01 -1.253234e-02
      t0    -3.281909e-14 -1.375092e-02 -1.253234e-02  2.721841e-02

# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9079896 -0.9327237
      m_lMd -0.3490792 -0.5356950
      s     -0.9134089 -0.5523713
      t0    -2.0958893 -1.7942577

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8754524 -0.4624233 -0.6844776 -2.0209660 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m       m_lMd          s         t0
      m     0.001787335 0.000000000 0.00000000 0.00000000
      m_lMd 0.000000000 0.009365169 0.00000000 0.00000000
      s     0.000000000 0.000000000 0.03372024 0.00000000
      t0    0.000000000 0.000000000 0.00000000 0.01951607

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                 as1t       bd6t
      m     -0.831102 -0.9959333
      m_lMd -0.301806 -0.6283231
      s     -1.008355 -0.4343118
      t0    -2.391361 -1.6418025

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8974612 -0.6581013 -0.1735832 -1.2915573 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                     m      m_lMd         s        t0
      m     0.01761997 0.00000000 0.0000000 0.0000000
      m_lMd 0.00000000 0.07205873 0.0000000 0.0000000
      s     0.00000000 0.00000000 0.1130024 0.1307153
      t0    0.00000000 0.00000000 0.1307153 0.2134335

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9897073 -0.9585400
      m_lMd -0.3399877 -0.4782692
      s     -0.8737060 -0.5818668
      t0    -1.8503583 -1.7458396

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor,
        blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single  -541   0 -440    1 -358     1         83  -523  -606 -606
      diag    -557   0 -415    0 -294     0        121  -536  -628 -658
      factor  -641   1 3126    0 4995     0       1868  1258  -605 -610
      blocked -562   0 -245    0  -30     0        215  -461  -643 -676

