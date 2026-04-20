# ordered_probit

    Code
      init_chains(ord_probit_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                              1
      location_Slow  -0.2660328
      location_Smid   0.3207204
      location_Shigh -0.7470537
      cut_lRlow      -0.2179812
      cut_lRmid      -0.5137169
      
      
      $stage
      [1] "init"
      
      $subj_ll
            [,1]
      1 -36.7689
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_ord_probit, d_ord_probit, n_trials = 10)
    Output
         trials subjects    S    R
      1       1        1  low  low
      4       2        1  low  mid
      7       3        1  low  mid
      10      4        1  low  low
      13      5        1  low  low
      16      6        1  low  low
      19      7        1  low  low
      22      8        1  low  low
      25      9        1  low  low
      28     10        1  low  low
      31     11        1  mid  low
      34     12        1  mid  mid
      37     13        1  mid  mid
      40     14        1  mid  low
      43     15        1  mid  mid
      46     16        1  mid high
      49     17        1  mid  mid
      52     18        1  mid  mid
      55     19        1  mid  mid
      58     20        1  mid high
      61     21        1 high  mid
      64     22        1 high high
      67     23        1 high  mid
      70     24        1 high  mid
      73     25        1 high  mid
      76     26        1 high high
      79     27        1 high  low
      82     28        1 high  low
      85     29        1 high high
      88     30        1 high  mid

# ordered_logit

    Code
      init_chains(ord_logit_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $alpha
      , , 1
      
                               1
      location_Sleft  -1.8150926
      location_Sright  0.3304096
      cut             -1.1421557
      
      
      $stage
      [1] "init"
      
      $subj_ll
             [,1]
      1 -13.42335
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_ord_logit, d_ord_logit, n_trials = 10)
    Output
         trials subjects     S     R
      1       1        1  left  left
      3       2        1  left  left
      5       3        1  left right
      7       4        1  left  left
      9       5        1  left right
      11      6        1  left  left
      13      7        1  left  left
      15      8        1  left  left
      17      9        1  left  left
      19     10        1  left  left
      21     11        1 right  left
      23     12        1 right  left
      25     13        1 right right
      27     14        1 right  left
      29     15        1 right  left
      31     16        1 right right
      33     17        1 right right
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
      utility_Sleft  1.19461175
      utility_Sright 0.77458193
      utility_Sup    0.08710445
      
      
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
      1       1        1  left  left
      4       2        1  left right
      7       3        1  left    up
      10      4        1  left right
      13      5        1  left    up
      16      6        1  left right
      19      7        1  left  left
      22      8        1  left  left
      25      9        1  left right
      28     10        1  left  left
      31     11        1 right  left
      34     12        1 right  left
      37     13        1 right right
      40     14        1 right  left
      43     15        1 right  left
      46     16        1 right right
      49     17        1 right  left
      52     18        1 right    up
      55     19        1 right right
      58     20        1 right  left
      61     21        1    up  left
      64     22        1    up right
      67     23        1    up right
      70     24        1    up    up
      73     25        1    up right
      76     26        1    up right
      79     27        1    up  left
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
      utility_Sleft  1.19461175
      utility_Sright 0.77458193
      utility_Sup    0.08710445
      
      
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
      1       1        1  left    up
      4       2        1  left right
      7       3        1  left  left
      10      4        1  left  left
      13      5        1  left    up
      16      6        1  left  left
      19      7        1  left right
      22      8        1  left  left
      25      9        1  left right
      28     10        1  left right
      31     11        1 right right
      34     12        1 right right
      37     13        1 right  left
      40     14        1 right  left
      43     15        1 right right
      46     16        1 right    up
      49     17        1 right  left
      52     18        1 right    up
      55     19        1 right right
      58     20        1 right  left
      61     21        1    up  left
      64     22        1    up right
      67     23        1    up  left
      70     24        1    up right
      73     25        1    up  left
      76     26        1    up right
      79     27        1    up  left
      82     28        1    up    up
      85     29        1    up    up
      88     30        1    up right

