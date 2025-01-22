# run_factor

    Code
      LNR_factor[[1]]$samples$theta_lambda[, , idx]
    Output
                    F1         F2
      m     0.65189145  0.0000000
      m_lMd 0.43051042  0.2739173
      s     0.38441388 -0.3794117
      t0    0.06447518  0.2381119

---

    Code
      LNR_factor[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9326169 -0.8480190
      m_lMd -0.2972457 -0.6236316
      s     -0.8387106 -0.5695598
      t0    -1.8680927 -1.7047142

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.4512140 -0.2588467 -0.6315056 -1.1288822 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                     m      m_lMd           s          t0
      m     0.44855967 0.28064606  0.25059612  0.04203082
      m_lMd 0.28064606 0.27654533  0.06156676  0.09298022
      s     0.25059612 0.06156676  0.35488555 -0.06555729
      t0    0.04203082 0.09298022 -0.06555729  0.18663637

# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8478104 -0.9921555
      m_lMd -0.3865287 -0.6148483
      s     -0.9142786 -0.4746343
      t0    -2.0645389 -1.5692885

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9378624 -0.3402953 -0.6335480 -1.6310249 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m     m_lMd          s        t0
      m     0.003458138 0.0000000 0.00000000 0.0000000
      m_lMd 0.000000000 0.1977962 0.00000000 0.0000000
      s     0.000000000 0.0000000 0.01151789 0.0000000
      t0    0.000000000 0.0000000 0.00000000 0.2529093

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9260897 -0.8587017
      m_lMd -0.2945346 -0.4345684
      s     -0.9381094 -0.7072942
      t0    -1.9666953 -1.8735804

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.9364989 -0.3636221 -0.5077303 -2.0114817 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd           s          t0
      m     0.001114941 0.00000000  0.00000000  0.00000000
      m_lMd 0.000000000 0.00211203  0.00000000  0.00000000
      s     0.000000000 0.00000000  0.26341869 -0.05385428
      t0    0.000000000 0.00000000 -0.05385428  0.06735875

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8234659 -0.9020162
      m_lMd -0.3214201 -0.6333021
      s     -0.9439186 -0.5232201
      t0    -2.2166120 -1.6597359

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor,
        blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD   wMD  DIC  wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single  -494 0.000 -395 0.000 -324     0         71  -465  -536 -536
      diag    -517 0.947 -450 0.985 -380     1         70  -521  -577 -591
      factor  -512 0.053 -349 0.000 -242     0        108  -457  -550 -565
      blocked -484 0.000 -442 0.015 -360     0         82  -524  -576 -605

