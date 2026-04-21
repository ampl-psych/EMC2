# ordered_probit

    Code
      init_chains(ord_probit_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                              1
      location_Slow  -0.7422707
      location_Smid   0.7732297
      location_Shigh  0.6825030
      cut_lRlow      -0.2179361
      cut_lRmid      -0.6387820
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -35.53273
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_ord_probit, d_ord_probit, n_trials = 10)
    Output
         trials subjects    S    R
      1       1        1  low  low
      4       2        1  low  low
      7       3        1  low  low
      10      4        1  low  mid
      13      5        1  low  mid
      16      6        1  low  mid
      19      7        1  low  low
      22      8        1  low  low
      25      9        1  low  mid
      28     10        1  low  mid
      31     11        1  mid  mid
      34     12        1  mid  mid
      37     13        1  mid high
      40     14        1  mid  mid
      43     15        1  mid  mid
      46     16        1  mid  mid
      49     17        1  mid  low
      52     18        1  mid high
      55     19        1  mid high
      58     20        1  mid  mid
      61     21        1 high  mid
      64     22        1 high  mid
      67     23        1 high high
      70     24        1 high high
      73     25        1 high high
      76     26        1 high high
      79     27        1 high high
      82     28        1 high high
      85     29        1 high  mid
      88     30        1 high high

# ordered_logit

    Code
      init_chains(ord_logit_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                               1
      location_Sleft  -0.5126331
      location_Sright  0.8401962
      cut              0.6442838
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -13.18434
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_ord_logit, d_ord_logit, n_trials = 10)
    Output
         trials subjects     S     R
      1       1        1  left  left
      3       2        1  left right
      5       3        1  left right
      7       4        1  left  left
      9       5        1  left  left
      11      6        1  left right
      13      7        1  left  left
      15      8        1  left  left
      17      9        1  left  left
      19     10        1  left right
      21     11        1 right right
      23     12        1 right right
      25     13        1 right right
      27     14        1 right right
      29     15        1 right right
      31     16        1 right right
      33     17        1 right  left
      35     18        1 right right
      37     19        1 right right
      39     20        1 right right

# multinomial_logit

    Code
      init_chains(mnl_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                               1
      utility_Sleft   0.02412719
      utility_Sright -2.74441709
      utility_Sup     0.54761193
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -32.95837
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_mnl, d_mnl, n_trials = 10)
    Output
         trials subjects     S     R
      1       1        1  left    up
      4       2        1  left right
      7       3        1  left right
      10      4        1  left  left
      13      5        1  left    up
      16      6        1  left    up
      19      7        1  left right
      22      8        1  left    up
      25      9        1  left    up
      28     10        1  left  left
      31     11        1 right    up
      34     12        1 right  left
      37     13        1 right right
      40     14        1 right  left
      43     15        1 right right
      46     16        1 right  left
      49     17        1 right right
      52     18        1 right right
      55     19        1 right  left
      58     20        1 right right
      61     21        1    up right
      64     22        1    up    up
      67     23        1    up  left
      70     24        1    up right
      73     25        1    up    up
      76     26        1    up    up
      79     27        1    up right
      82     28        1    up    up
      85     29        1    up    up
      88     30        1    up  left

# multinomial_probit

    Code
      init_chains(mnp_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                              1
      utility_Sleft   0.1730763
      utility_Sright -0.1974613
      utility_Sup     2.4566301
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -32.95837
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_mnp, d_mnp, n_trials = 10)
    Output
         trials subjects     S     R
      1       1        1  left  left
      4       2        1  left    up
      7       3        1  left right
      10      4        1  left    up
      13      5        1  left  left
      16      6        1  left  left
      19      7        1  left    up
      22      8        1  left  left
      25      9        1  left  left
      28     10        1  left  left
      31     11        1 right    up
      34     12        1 right  left
      37     13        1 right    up
      40     14        1 right  left
      43     15        1 right    up
      46     16        1 right    up
      49     17        1 right right
      52     18        1 right  left
      55     19        1 right right
      58     20        1 right right
      61     21        1    up  left
      64     22        1    up  left
      67     23        1    up right
      70     24        1    up    up
      73     25        1    up    up
      76     26        1    up  left
      79     27        1    up    up
      82     28        1    up right
      85     29        1    up right
      88     30        1    up  left

