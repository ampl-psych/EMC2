# run_factor

    Code
      LNR_factor[[1]]$samples$theta_lambda[, , idx]
    Output
                    F1          F2
      m      0.4042440 0.000000000
      m_lMd -0.4318224 0.379480463
      s      0.6085597 0.101328100
      t0     0.1192645 0.006457266

---

    Code
      LNR_factor[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.9538521 -0.9353932
      m_lMd -0.3830433 -0.6787215
      s     -0.8030984 -0.4866918
      t0    -1.7952026 -1.6199259

---

    Code
      LNR_factor[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.5966560 -1.0288748 -0.2052301 -1.6511029 

---

    Code
      LNR_factor[[1]]$samples$theta_var[, , idx]
    Output
                      m       m_lMd           s          t0
      m      0.19259767 -0.17456160  0.24600658  0.04821194
      m_lMd -0.17456160  0.35581441 -0.22433767 -0.04905066
      s      0.24600658 -0.22433767  0.45313739  0.07323385
      t0     0.04821194 -0.04905066  0.07323385  0.03879326

# run_diag

    Code
      LNR_diag[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.8697654 -0.8967558
      m_lMd -0.3358620 -0.5448376
      s     -0.8933652 -0.5917299
      t0    -1.9745522 -1.6971386

---

    Code
      LNR_diag[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.8247087 -0.4798699 -0.6709188 -1.9290857 

---

    Code
      LNR_diag[[1]]$samples$theta_var[, , idx]
    Output
                      m      m_lMd          s         t0
      m     0.001956766 0.00000000 0.00000000 0.00000000
      m_lMd 0.000000000 0.01533777 0.00000000 0.00000000
      s     0.000000000 0.00000000 0.03421786 0.00000000
      t0    0.000000000 0.00000000 0.00000000 0.01902554

# run_blocked

    Code
      LNR_blocked[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.7045390 -0.9539209
      m_lMd -0.2739985 -0.6802633
      s     -1.1179201 -0.4835457
      t0    -2.7406658 -1.5598814

---

    Code
      LNR_blocked[[1]]$samples$theta_mu[, idx]
    Output
               m      m_lMd          s         t0 
      -0.7859169 -0.3300567  0.7986264 -0.6736691 

---

    Code
      LNR_blocked[[1]]$samples$theta_var[, , idx]
    Output
                     m     m_lMd         s       t0
      m     0.04326159 0.0000000 0.0000000 0.000000
      m_lMd 0.00000000 0.2575511 0.0000000 0.000000
      s     0.00000000 0.0000000 0.9349592 1.020365
      t0    0.00000000 0.0000000 1.0203646 2.492648

# run_single

    Code
      LNR_single[[1]]$samples$alpha[, , idx]
    Output
                  as1t       bd6t
      m     -0.7546871 -0.9552423
      m_lMd -0.2646319 -0.6581083
      s     -0.8971468 -0.3139218
      t0    -2.2014758 -1.5023459

# run_bridge

    Code
      compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor,
        blocked = LNR_blocked), stage = "preburn", cores_for_props = 1)
    Output
                MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      single  -476   0 -313    1 -191     1        122  -436  -558 -558
      diag    -501   1    9    0  302     0        294  -285  -557 -579
      factor  -422   0 2584    0 4080     0       1496  1089  -392 -407
      blocked -152   0   60    0  343     0        284  -224  -478 -507

