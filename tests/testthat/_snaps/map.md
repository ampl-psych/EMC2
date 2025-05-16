# map

    Code
      mapped_pars(des)
    Output
      $v 
        covariate           E         
       -0.513716908811104  speed     : v - 0.514 * v_covariate
       0.601813121611908   neutral   : v + 0.602 * v_covariate + v_Eneutral + 0.602 * v_covariate:Eneutral
       -1.54982186331112   accuracy  : v - 1.55 * v_covariate + v_Eaccuracy - 1.55 * v_covariate:Eaccuracy
       -1.70962282511647   speed     : v - 1.71 * v_covariate
       0.7701488026301     neutral   : v + 0.77 * v_covariate + v_Eneutral + 0.77 * v_covariate:Eneutral
       -0.716873043602934  accuracy  : v - 0.717 * v_covariate + v_Eaccuracy - 0.717 * v_covariate:Eaccuracy
      
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
                E     S   covariate    lR      v sv     B A    t0     b
      1     speed  left -0.29797771  left  1.637  1 3.048 0 0.292 3.048
      2     speed  left -0.29797771 right  1.637  1 3.048 0 0.292 3.048
      3   neutral  left  0.34063680  left  0.805  1 1.793 0 0.292 1.793
      4   neutral  left  0.34063680 right  0.805  1 1.793 0 0.292 1.793
      5  accuracy  left  0.13373072  left -0.213  1 3.051 0 0.292 3.051
      6  accuracy  left  0.13373072 right -0.213  1 3.051 0 0.292 3.051
      7     speed right  0.86266186  left  0.057  1 3.048 0 0.189 3.048
      8     speed right  0.86266186 right  0.057  1 3.048 0 0.189 3.048
      9   neutral right  0.05063779  left  0.895  1 1.793 0 0.189 1.793
      10  neutral right  0.05063779 right  0.895  1 1.793 0 0.189 1.793
      11 accuracy right  1.22458533  left -2.828  1 3.051 0 0.189 3.051
      12 accuracy right  1.22458533 right -2.828  1 3.051 0 0.189 3.051

---

    Code
      mapped_pars(prior(des, mu_mean = c(v_covariate = 1)))
    Output
                E     S   covariate    lR      v sv B A t0 b
      1     speed  left -0.09534976  left -0.095  1 1 0  1 1
      2     speed  left -0.09534976 right -0.095  1 1 0  1 1
      3   neutral  left -0.54328498  left -0.543  1 1 0  1 1
      4   neutral  left -0.54328498 right -0.543  1 1 0  1 1
      5  accuracy  left -2.15286166  left -2.153  1 1 0  1 1
      6  accuracy  left -2.15286166 right -2.153  1 1 0  1 1
      7     speed right -1.12658056  left -1.127  1 1 0  1 1
      8     speed right -1.12658056 right -1.127  1 1 0  1 1
      9   neutral right  0.67546458  left  0.675  1 1 0  1 1
      10  neutral right  0.67546458 right  0.675  1 1 0  1 1
      11 accuracy right -0.08025789  left -0.080  1 1 0  1 1
      12 accuracy right -0.08025789 right -0.080  1 1 0  1 1

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
      

