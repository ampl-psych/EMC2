# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9693333 -0.9640234
      m_lMd -0.4409364 -0.6000380
      s     -0.7190573 -0.5559688
      t0    -1.8426007 -1.7529013

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9377841 -0.6067421 -0.7638091 -1.6031017 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd          s         t0
      m     0.001234336 0.00000000 0.00000000 0.00000000
      m_lMd 0.000000000 0.03327357 0.00000000 0.00000000
      s     0.000000000 0.00000000 0.01251535 0.00000000
      t0    0.000000000 0.00000000 0.00000000 0.04757134

---

    Code
      compare(list(diag = LNR_diag), stage = "preburn", cores_for_props = 1)
    Output
             MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      diag -536   1  277     1 -85    1  201     1        286  -370  -610 -656
    Message
      p_waic produced warnings, use waic_warnings() to inspect

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.7416180 -0.9346588
      m_lMd -0.3533509 -0.7934335
      s     -0.9803318 -0.3123977
      t0    -2.6572651 -1.6084508

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.7849105 -0.3747564 -0.1156958 -1.4576668 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                     m      m_lMd         s        t0
      m     0.02201408 0.00000000 0.0000000 0.0000000
      m_lMd 0.00000000 0.07367456 0.0000000 0.0000000
      s     0.00000000 0.00000000 1.2417562 0.1181575
      t0    0.00000000 0.00000000 0.1181575 0.2302009

---

    Code
      compare(list(blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      blocked -528   1 4625     1 476    1 1040     1        564   -87  -631 -651
    Message
      p_waic produced warnings, use waic_warnings() to inspect

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.7569884 -0.8593122
      m_lMd -0.3453346 -0.6544262
      s     -0.9959459 -0.6699200
      t0    -2.2834271 -1.7788640

---

    Code
      compare(list(single = LNR_single), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -336   1  464     1  79    1  423     1        344  -265  -505 -609
    Message
      p_waic produced warnings, use waic_warnings() to inspect

