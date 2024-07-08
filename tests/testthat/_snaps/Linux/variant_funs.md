# run_factor

    Code
      LNR_factor[[1]]$samples$theta_lambda[, , idx]
    Output
                      F1         F2
      m      0.361174688  0.0000000
      m_lMd -0.444228520  0.5179415
      s     -0.295217709 -0.2682314
      t0    -0.002400008 -0.7269956

---

    Code
      LNR_factor[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8449149 -0.9391686
      m_lMd -0.3082034 -0.6155118
      s     -0.9379678 -0.4871466
      t0    -1.9409439 -1.5628295

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.5392683 -1.1564255 -0.1857980 -0.7779520 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                       m        m_lMd            s           t0
      m      0.271892650 -0.160444097 -0.106625164 -0.000866822
      m_lMd -0.160444097  0.502904983 -0.007784044 -0.375475004
      s     -0.106625164 -0.007784044  0.271083750  0.195711573
      t0    -0.000866822 -0.375475004  0.195711573  0.598383441

# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9378169 -0.9489708
      m_lMd -0.3076695 -0.6012228
      s     -0.8798202 -0.4407736
      t0    -1.8683764 -1.5734354

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9534561 -0.6384648 -0.7878486 -1.7974908 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd          s        t0
      m     0.000257467 0.00000000 0.00000000 0.0000000
      m_lMd 0.000000000 0.06233054 0.00000000 0.0000000
      s     0.000000000 0.00000000 0.03554651 0.0000000
      t0    0.000000000 0.00000000 0.00000000 0.0132861

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                   as1t        bd6t
      m      -0.6064332  -0.5266625
      m_lMd  -0.2040174  -0.3446261
      s      -1.2660982  -1.1346625
      t0    -10.4778607 -31.0242957

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -1.0458184 -0.1057472 -1.3440932  0.8481241 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                   m   m_lMd           s         t0
      m     6.905368 0.00000  0.00000000   0.000000
      m_lMd 0.000000 1.03098  0.00000000   0.000000
      s     0.000000 0.00000  0.01855142  -2.198102
      t0    0.000000 0.00000 -2.19810177 294.561468

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.6630457 -1.1170701
      m_lMd -0.3022816 -0.3745601
      s     -1.1542119 -0.3384082
      t0    -3.0194702 -1.5177944

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor), stage = "preburn",
      cores_for_props = 1)
    Output
               MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single -312   0  -62    0   97     0        159  -220  -379 -379
      diag   -477   1 -353    1 -234     1        119  -472  -568 -592
      factor  365   0 -112    0  114     0        225  -337  -551 -563

