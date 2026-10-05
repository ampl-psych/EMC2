# MTLNR snapshots: chain init, simulated data, natural-scale mapping

    Code
      init_chains(emc, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                          1
      m          -0.9685927
      m_lMTRUE    0.7061091
      s           1.4890213
      t0         -1.8150926
      rho         0.3304096
      c1         -1.1421557
      c1_lRright  0.1571934
      c2         -2.0654072
      c2_lRright -0.4405469
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -193.0817
      
      $idx
      [1] 1
      

---

    Code
      make_data(mtlnr_p, des, n_trials = 10)
    Output
         trials subjects     S     R        rt RR
      1       1        1  left  left 0.4533743  2
      3       2        1  left  left 0.8856033  2
      5       3        1  left right 0.9557567  2
      7       4        1  left  left 0.3779199  3
      9       5        1  left  left 0.7335829  1
      11      6        1  left  left 0.4334906  3
      13      7        1  left right 0.5207177  2
      15      8        1  left  left 0.3637794  1
      17      9        1  left  left 0.5339992  3
      19     10        1  left  left 0.6339255  1
      21     11        1 right right 0.5485097  2
      23     12        1 right  left 0.9004454  2
      25     13        1 right right 0.7672100  2
      27     14        1 right  left 0.4708543  2
      29     15        1 right  left 1.1818207  1
      31     16        1 right right 1.1531434  2
      33     17        1 right right 0.6551928  3
      35     18        1 right right 0.4301364  3
      37     19        1 right right 0.3731337  3
      39     20        1 right  left 1.9095959  2

---

    Code
      mapped_pars(des, mtlnr_p)
    Output
            S    lR    lM    m   s  t0 rho    c1    c2    d1    d2
      1  left  left  TRUE -1.1 0.8 0.3 0.5 0.400 0.600 0.368 0.670
      2  left right FALSE -0.5 0.8 0.3 0.5 0.442 0.491 0.393 0.643
      3 right  left FALSE -0.5 0.8 0.3 0.5 0.400 0.600 0.368 0.670
      4 right right  TRUE -1.1 0.8 0.3 0.5 0.442 0.491 0.393 0.643

# rating_summary snapshot

    Code
      rating_summary(dat, factors = "S")
    Output
         postn  cell    resp     p       q10       q50       q90
      1      1  left  left.3 0.325 0.3901944 0.5369874 0.7840765
      3      1  left  left.2 0.295 0.4172763 0.6404236 1.2040471
      5      1  left  left.1 0.175 0.4621618 0.6046966 0.9281359
      7      1  left right.1 0.125 0.5199095 0.7233067 0.8573340
      9      1  left right.2 0.050 0.3887281 0.6882900 1.0080987
      11     1  left right.3 0.030 0.3823335 0.6338415 0.6986068
      2      1 right  left.3 0.025 0.5149544 0.5771045 0.6678162
      4      1 right  left.2 0.090 0.4495120 0.5840303 0.8612179
      6      1 right  left.1 0.140 0.5278857 0.7171173 1.2413839
      8      1 right right.1 0.155 0.4863773 0.7785415 1.4929806
      10     1 right right.2 0.285 0.4522445 0.6361729 1.1644833
      12     1 right right.3 0.305 0.3849295 0.4896386 0.9770985

