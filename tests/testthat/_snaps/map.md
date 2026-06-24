# mapped_pars

    Code
      mapped_pars(des)
    Output
      $v 
        covariate           E         
       -0.968592726552943  speed     : v - 0.969 * v_covariate
       0.706109077785223   neutral   : v + 0.706 * v_covariate + v_Eneutral + 0.706 * v_covariate:Eneutral
       1.48902132366997    accuracy  : v + 1.49 * v_covariate + v_Eaccuracy + 1.49 * v_covariate:Eaccuracy
       -1.81509255136049   speed     : v - 1.82 * v_covariate
       0.330409581768214   neutral   : v + 0.33 * v_covariate + v_Eneutral + 0.33 * v_covariate:Eneutral
       -1.14215571117096   accuracy  : v - 1.14 * v_covariate + v_Eaccuracy - 1.14 * v_covariate:Eaccuracy
      
      $B 
        E         
       speed     : exp(B)
       neutral   : exp(B + B_Eneutral)
       accuracy  : exp(B + B_Eaccuracy)
      
      $t0 
        S      
       left   : exp(t0)
       right  : exp(t0 + t0_Sright)
      

---

    Code
      mapped_pars(des, p_vector = rnorm(length(sampled_pars(des))))
    Output
                E     S  covariate    lR      v sv     B A    t0     b
      1     speed  left -0.2721696  left  0.029  1 0.721 0 1.013 0.721
      2     speed  left -0.2721696 right  0.029  1 0.721 0 1.013 0.721
      3   neutral  left  1.2311057  left -3.677  1 1.071 0 1.013 1.071
      4   neutral  left  1.2311057 right -3.677  1 1.071 0 1.013 1.071
      5  accuracy  left -1.3606844  left  2.111  1 0.422 0 1.013 0.422
      6  accuracy  left -1.3606844 right  2.111  1 0.422 0 1.013 0.422
      7     speed right -0.3208495  left  0.098  1 0.721 0 0.535 0.721
      8     speed right -0.3208495 right  0.098  1 0.721 0 0.535 0.721
      9   neutral right -1.1235776  left -2.022  1 1.071 0 0.535 1.071
      10  neutral right -1.1235776 right -2.022  1 1.071 0 0.535 1.071
      11 accuracy right  1.0520209  left  2.864  1 0.422 0 0.535 0.422
      12 accuracy right  1.0520209 right  2.864  1 0.422 0 0.535 0.422

---

    Code
      mapped_pars(prior(des, mu_mean = c(v_covariate = 1)))
    Output
                E     S  covariate    lR      v sv B A t0 b
      1     speed  left -0.9865591  left -0.987  1 1 0  1 1
      2     speed  left -0.9865591 right -0.987  1 1 0  1 1
      3   neutral  left  0.8184259  left  0.818  1 1 0  1 1
      4   neutral  left  0.8184259 right  0.818  1 1 0  1 1
      5  accuracy  left -0.7822744  left -0.782  1 1 0  1 1
      6  accuracy  left -0.7822744 right -0.782  1 1 0  1 1
      7     speed right -0.5559056  left -0.556  1 1 0  1 1
      8     speed right -0.5559056 right -0.556  1 1 0  1 1
      9   neutral right -0.8314149  left -0.831  1 1 0  1 1
      10  neutral right -0.8314149 right -0.831  1 1 0  1 1
      11 accuracy right -1.9219137  left -1.922  1 1 0  1 1
      12 accuracy right -1.9219137 right -1.922  1 1 0  1 1

---

    Code
      mapped_pars(samples_LNR)
    Output
                E     S    lR    lM      m     s    t0
      1     speed  left  left  TRUE -1.225 0.587 0.197
      2     speed  left right FALSE -0.711 0.587 0.197
      3   neutral  left  left  TRUE -1.225 0.587 0.197
      4   neutral  left right FALSE -0.711 0.587 0.197
      5  accuracy  left  left  TRUE -1.225 0.587 0.197
      6  accuracy  left right FALSE -0.711 0.587 0.197
      7     speed right  left FALSE -0.711 0.587 0.197
      8     speed right right  TRUE -1.225 0.587 0.197
      9   neutral right  left FALSE -0.711 0.587 0.197
      10  neutral right right  TRUE -1.225 0.587 0.197
      11 accuracy right  left FALSE -0.711 0.587 0.197
      12 accuracy right right  TRUE -1.225 0.587 0.197

---

    Code
      mapped_pars(get_prior(samples_LNR))
    Output
                E     S    lR    lM m s t0
      1     speed  left  left  TRUE 0 1  1
      2     speed  left right FALSE 0 1  1
      3   neutral  left  left  TRUE 0 1  1
      4   neutral  left right FALSE 0 1  1
      5  accuracy  left  left  TRUE 0 1  1
      6  accuracy  left right FALSE 0 1  1
      7     speed right  left FALSE 0 1  1
      8     speed right right  TRUE 0 1  1
      9   neutral right  left FALSE 0 1  1
      10  neutral right right  TRUE 0 1  1
      11 accuracy right  left FALSE 0 1  1
      12 accuracy right right  TRUE 0 1  1

---

    Code
      mapped_pars(get_design(samples_LNR))
    Output
      $m 
        lM     
       TRUE   : m + 0.5 * m_lMd
       FALSE  : m - 0.5 * m_lMd
      

# map

    Code
      credint(samples_LNR, selection = "mu", map = "E")
    Output
      $mu
                     2.5%    50%  97.5%
      m_Eaccuracy  -1.008 -0.964 -0.919
      m_Eneutral   -1.008 -0.964 -0.919
      m_Espeed     -1.008 -0.964 -0.919
      s_Eaccuracy   0.554  0.583  0.619
      s_Eneutral    0.554  0.583  0.619
      s_Espeed      0.554  0.583  0.619
      t0_Eaccuracy  0.177  0.193  0.206
      t0_Eneutral   0.177  0.193  0.206
      t0_Espeed     0.177  0.193  0.206
      

---

    Code
      credint(samples_LNR, selection = "mu", map = list(~ E * S))
    Output
      $mu
                            2.5%    50%  97.5%
      m_(Intercept)       -1.008 -0.964 -0.919
      m_Eneutral          -1.008 -0.964 -0.919
      m_Eaccuracy         -1.008 -0.964 -0.919
      m_Sright            -1.008 -0.964 -0.919
      m_Eneutral:Sright   -1.008 -0.964 -0.919
      m_Eaccuracy:Sright  -1.008 -0.964 -0.919
      s_(Intercept)        0.554  0.583  0.619
      s_Eneutral           0.554  0.583  0.619
      s_Eaccuracy          0.554  0.583  0.619
      s_Sright             0.554  0.583  0.619
      s_Eneutral:Sright    0.554  0.583  0.619
      s_Eaccuracy:Sright   0.554  0.583  0.619
      t0_(Intercept)       0.177  0.193  0.206
      t0_Eneutral          0.177  0.193  0.206
      t0_Eaccuracy         0.177  0.193  0.206
      t0_Sright            0.177  0.193  0.206
      t0_Eneutral:Sright   0.177  0.193  0.206
      t0_Eaccuracy:Sright  0.177  0.193  0.206
      

---

    Code
      credint(samples_LNR, selection = "mu", map = TRUE)
    Output
      $mu
                  2.5%    50%  97.5%
      m_lMFALSE -0.752 -0.700 -0.655
      m_lMTRUE  -1.271 -1.228 -1.167
      s          0.554  0.583  0.619
      t0         0.177  0.193  0.206
      

