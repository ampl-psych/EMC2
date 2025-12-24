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

---

    Code
      compare(list(diag = LNR_diag), stage = "preburn", cores_for_props = 1)
    Output
             MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      diag -579   1 1066    1 1916     1        851   215  -607 -635

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9668625 -1.0003825
      m_lMd -0.2819924 -0.5445420
      s     -0.8501189 -0.5963007
      t0    -1.9733013 -1.7466873

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -1.0052941 -1.5321301 -0.9540201 -1.8960213 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                      m     m_lMd          s         t0
      m     0.003426314 0.0000000 0.00000000 0.00000000
      m_lMd 0.000000000 0.7972991 0.00000000 0.00000000
      s     0.000000000 0.0000000 0.02009918 0.02058433
      t0    0.000000000 0.0000000 0.02058433 0.04993745

---

    Code
      compare(list(blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      blocked -527   1 6444    1 9984     1       3540  2904  -611 -636

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -1.0064461 -0.7239915
      m_lMd -0.3903504 -0.8186065
      s     -0.5314815 -0.2610311
      t0    -1.6904537 -1.8390487

---

    Code
      compare(list(single = LNR_single), stage = "preburn", cores_for_props = 1)
    Output
               MD wMD DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -390   1 245    1  599     1        354  -109  -463 -463

