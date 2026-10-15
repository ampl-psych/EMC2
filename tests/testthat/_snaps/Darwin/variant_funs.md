# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9743351 -0.9833168
      m_lMd -0.4386883 -0.5367333
      s     -0.7338663 -0.5722065
      t0    -1.8065456 -1.7854034

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -1.0061445 -0.5730070 -0.5592815 -1.7125241 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                     m     m_lMd          s         t0
      m     0.01563927 0.0000000 0.00000000 0.00000000
      m_lMd 0.00000000 0.0125815 0.00000000 0.00000000
      s     0.00000000 0.0000000 0.02627957 0.00000000
      t0    0.00000000 0.0000000 0.00000000 0.02260007

---

    Code
      compare(list(diag = LNR_diag), stage = "preburn", cores_for_props = 1)
    Output
             MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      diag -400   1 5075     1 491    1 1064     1        573   -82  -562 -655
    Message
      p_waic produced warnings, use waic_warnings() to inspect

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8034791 -1.0238329
      m_lMd -0.3277423 -0.4290092
      s     -1.0259435 -0.4402202
      t0    -2.5499149 -1.6199128

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8421076 -0.2267080 -0.7557078 -1.5410847 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                    m      m_lMd          s         t0
      m     0.1579914 0.00000000 0.00000000 0.00000000
      m_lMd 0.0000000 0.01308866 0.00000000 0.00000000
      s     0.0000000 0.00000000 0.05507848 0.02764013
      t0    0.0000000 0.00000000 0.02764013 0.33836474

---

    Code
      compare(list(blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      blocked -368   1 3530     1 722    1 1415     1        693    29  -546 -663
    Message
      p_waic produced warnings, use waic_warnings() to inspect

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.0103237 -1.0523847
      m_lMd -0.3734309 -0.5476712
      s     -0.3097086 -0.6119045
      t0    -1.6697117 -1.7132085

---

    Code
      compare(list(single = LNR_single), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD WAIC wWAIC DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -135   1  882     1 485    1 1007     1        522   -36  -318 -558
    Message
      p_waic produced warnings, use waic_warnings() to inspect

