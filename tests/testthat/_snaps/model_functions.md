# LNR

    Code
      init_chains(LNR_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                        [,1]
      m           -0.9685927
      m_lMd        0.7061091
      m_Eneutral   1.4890213
      m_Eaccuracy -1.8150926
      s            0.3304096
      t0          -1.1421557
      
      $theta_var
      , , 1
      
                            m        m_lMd   m_Eneutral m_Eaccuracy           s
      m            0.15548530 -0.120772081  0.018078294  0.04633176 -0.04484492
      m_lMd       -0.12077208  0.184889995 -0.009471693 -0.04763271  0.02012377
      m_Eneutral   0.01807829 -0.009471693  0.058583553 -0.01291727 -0.02728142
      m_Eaccuracy  0.04633176 -0.047632707 -0.012917266  0.09477040 -0.01673469
      s           -0.04484492  0.020123770 -0.027281415 -0.01673469  0.07945096
      t0           0.03785971 -0.018083101 -0.008237056  0.01937868 -0.01285052
                            t0
      m            0.037859710
      m_lMd       -0.018083101
      m_Eneutral  -0.008237056
      m_Eaccuracy  0.019378684
      s           -0.012850519
      t0           0.076944459
      
      
      $a_half
                       [,1]
      m           0.4056687
      m_lMd       1.9319610
      m_Eneutral  0.2242516
      m_Eaccuracy 0.4969591
      s           0.3334446
      t0          0.2464756
      
      $epsilon
           [,1]
      as1t  1.5
      bd6t  1.5
      
      $origin
           [,1]
      as1t    2
      bd6t    2
      
      $alpha
      , , 1
      
                        as1t       bd6t
      m           -0.6599300 -0.5836701
      m_lMd        0.3772913  0.1802469
      m_Eneutral   1.6147295  1.7008123
      m_Eaccuracy -2.1286744 -1.3967109
      s            0.5289465  0.4253892
      t0          -1.3065202 -1.3823828
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -891.6900
      bd6t -861.8978
      
      $last_theta_var_inv
                [,1]       [,2]       [,3]      [,4]      [,5]      [,6]
      [1,] 17.928942  9.9292394 -2.6234713 -2.086301  5.400555 -5.341685
      [2,]  9.929239 11.9481645  0.3450765  2.123976  2.803122 -2.107418
      [3,] -2.623471  0.3450765 23.5046084  5.314472  8.257526  3.928791
      [4,] -2.086301  2.1239763  5.3144716 14.075925  2.918336 -0.963037
      [5,]  5.400555  2.8031216  8.2575264  2.918336 18.577442  1.253107
      [6,] -5.341685 -2.1074180  3.9287913 -0.963037  1.253107 16.001843
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_LNR, design_LNR, n_trials = 10)
    Output
          subjects        E     S trials     R        rt
      1       as1t    speed  left      1 right 0.6524027
      5       as1t  neutral  left      1 right 0.1966150
      9       as1t accuracy  left      1 right 0.2066419
      13      as1t    speed right      1  left 0.2051926
      17      as1t  neutral right      1  left 0.2110602
      21      as1t accuracy right      1  left 0.2865263
      25      as1t    speed  left      2 right 1.3896176
      29      as1t  neutral  left      2 right 0.2241148
      33      as1t accuracy  left      2  left 1.2610216
      37      as1t    speed right      2  left 0.2250366
      41      as1t  neutral right      2  left 1.9593321
      45      as1t accuracy right      2  left 0.1975113
      49      as1t    speed  left      3 right 0.2191481
      53      as1t  neutral  left      3 right 0.2062084
      57      as1t accuracy  left      3 right 0.3458841
      61      as1t    speed right      3  left 0.1917125
      65      as1t  neutral right      3 right 0.2310378
      69      as1t accuracy right      3 right 0.8798443
      73      as1t    speed  left      4 right 0.8798883
      77      as1t  neutral  left      4 right 0.2743619
      81      as1t accuracy  left      4 right 0.1984744
      85      as1t    speed right      4  left 0.3659923
      89      as1t  neutral right      4  left 0.3305002
      93      as1t accuracy right      4  left 0.3787522
      97      as1t    speed  left      5  left 0.3171107
      101     as1t  neutral  left      5 right 0.2833863
      105     as1t accuracy  left      5 right 0.5533500
      109     as1t    speed right      5  left 0.6168620
      113     as1t  neutral right      5  left 0.2478463
      117     as1t accuracy right      5  left 0.2627446
      121     as1t    speed  left      6 right 0.4980723
      125     as1t  neutral  left      6  left 1.1084587
      129     as1t accuracy  left      6 right 0.2053519
      133     as1t    speed right      6 right 0.3226604
      137     as1t  neutral right      6  left 0.2795828
      141     as1t accuracy right      6  left 0.6148091
      145     as1t    speed  left      7 right 0.3093358
      149     as1t  neutral  left      7  left 0.5855643
      153     as1t accuracy  left      7 right 0.1947236
      157     as1t    speed right      7  left 0.2508862
      161     as1t  neutral right      7 right 0.2745357
      165     as1t accuracy right      7  left 0.2003754
      169     as1t    speed  left      8 right 5.6549638
      173     as1t  neutral  left      8 right 0.3006800
      177     as1t accuracy  left      8 right 0.2202169
      181     as1t    speed right      8  left 0.2055304
      185     as1t  neutral right      8 right 0.2493889
      189     as1t accuracy right      8  left 1.0355901
      193     as1t    speed  left      9 right 0.2702954
      197     as1t  neutral  left      9 right 0.1973118
      201     as1t accuracy  left      9 right 0.5270859
      205     as1t    speed right      9  left 0.2119288
      209     as1t  neutral right      9 right 0.3214212
      213     as1t accuracy right      9 right 0.3095826
      217     as1t    speed  left     10 right 0.2301976
      221     as1t  neutral  left     10 right 0.2422833
      225     as1t accuracy  left     10  left 1.8152301
      229     as1t    speed right     10 right 0.7174467
      233     as1t  neutral right     10  left 0.7369551
      237     as1t accuracy right     10 right 0.2188424
      3       bd6t    speed  left      1 right 0.3424278
      7       bd6t  neutral  left      1  left 2.8141979
      11      bd6t accuracy  left      1  left 0.2396563
      15      bd6t    speed right      1  left 0.6497393
      19      bd6t  neutral right      1  left 0.3641709
      23      bd6t accuracy right      1  left 0.3796578
      27      bd6t    speed  left      2 right 0.4302404
      31      bd6t  neutral  left      2 right 0.2140387
      35      bd6t accuracy  left      2 right 0.2013967
      39      bd6t    speed right      2 right 0.2889413
      43      bd6t  neutral right      2  left 0.4042207
      47      bd6t accuracy right      2  left 0.2136079
      51      bd6t    speed  left      3 right 0.3795840
      55      bd6t  neutral  left      3  left 0.5196810
      59      bd6t accuracy  left      3 right 0.2046349
      63      bd6t    speed right      3  left 0.3174439
      67      bd6t  neutral right      3 right 0.2892084
      71      bd6t accuracy right      3  left 2.6703282
      75      bd6t    speed  left      4 right 0.4029771
      79      bd6t  neutral  left      4 right 0.4149557
      83      bd6t accuracy  left      4  left 0.3110277
      87      bd6t    speed right      4 right 0.4503768
      91      bd6t  neutral right      4  left 0.2765337
      95      bd6t accuracy right      4  left 0.2469175
      99      bd6t    speed  left      5 right 0.4186540
      103     bd6t  neutral  left      5 right 0.2202781
      107     bd6t accuracy  left      5 right 0.2318917
      111     bd6t    speed right      5  left 0.2393652
      115     bd6t  neutral right      5  left 0.2978355
      119     bd6t accuracy right      5  left 0.2598228
      123     bd6t    speed  left      6 right 1.1033846
      127     bd6t  neutral  left      6 right 0.2024674
      131     bd6t accuracy  left      6 right 0.5196946
      135     bd6t    speed right      6 right 0.3401150
      139     bd6t  neutral right      6 right 0.3034629
      143     bd6t accuracy right      6  left 0.3102542
      147     bd6t    speed  left      7 right 0.5393920
      151     bd6t  neutral  left      7  left 0.2367627
      155     bd6t accuracy  left      7  left 0.3067202
      159     bd6t    speed right      7  left 0.3964789
      163     bd6t  neutral right      7  left 0.2423918
      167     bd6t accuracy right      7 right 0.5639506
      171     bd6t    speed  left      8 right 0.3099389
      175     bd6t  neutral  left      8 right 1.2686087
      179     bd6t accuracy  left      8 right 0.8795366
      183     bd6t    speed right      8  left 0.6813552
      187     bd6t  neutral right      8  left 0.2541792
      191     bd6t accuracy right      8  left 0.2873379
      195     bd6t    speed  left      9 right 0.4286412
      199     bd6t  neutral  left      9 right 2.1251756
      203     bd6t accuracy  left      9 right 0.6344836
      207     bd6t    speed right      9 right 0.2514761
      211     bd6t  neutral right      9 right 1.3794866
      215     bd6t accuracy right      9  left 0.2744479
      219     bd6t    speed  left     10 right 0.2113637
      223     bd6t  neutral  left     10 right 0.3336251
      227     bd6t accuracy  left     10  left 0.7346564
      231     bd6t    speed right     10  left 0.3938712
      235     bd6t  neutral right     10  left 0.3524264
      239     bd6t accuracy right     10  left 0.3341685

# LBA

    Code
      init_chains(LBA_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                        [,1]
      v           -0.9685927
      v_lMd        0.7061091
      sv_lMTRUE    1.4890213
      B           -1.8150926
      B_Eneutral   0.3304096
      B_Eaccuracy -1.1421557
      B_lRright    0.1571934
      A           -2.0654072
      t0          -0.4405469
      
      $theta_var
      , , 1
      
                             v         v_lMd    sv_lMTRUE             B   B_Eneutral
      v            0.059329585 -0.0104497238 -0.010774234  0.0021450832 -0.018128601
      v_lMd       -0.010449724  0.0559162133 -0.002872401  0.0014127637  0.004728529
      sv_lMTRUE   -0.010774234 -0.0028724006  0.059773323  0.0087906529  0.022576320
      B            0.002145083  0.0014127637  0.008790653  0.0447136854  0.001458968
      B_Eneutral  -0.018128601  0.0047285288  0.022576320  0.0014589675  0.038522063
      B_Eaccuracy  0.018325623  0.0019110134 -0.037503336 -0.0100477306 -0.023423360
      B_lRright   -0.012850266 -0.0030964417  0.020333128  0.0098693855  0.019288722
      A           -0.009137082 -0.0181065189  0.007641907 -0.0029721272  0.008676374
      t0          -0.019762600  0.0006819411  0.001329633  0.0006671515  0.007687211
                   B_Eaccuracy    B_lRright            A            t0
      v            0.018325623 -0.012850266 -0.009137082 -0.0197625997
      v_lMd        0.001911013 -0.003096442 -0.018106519  0.0006819411
      sv_lMTRUE   -0.037503336  0.020333128  0.007641907  0.0013296327
      B           -0.010047731  0.009869385 -0.002972127  0.0006671515
      B_Eneutral  -0.023423360  0.019288722  0.008676374  0.0076872108
      B_Eaccuracy  0.078051779 -0.028702651 -0.004381484  0.0022878527
      B_lRright   -0.028702651  0.042243108  0.005224348 -0.0090916441
      A           -0.004381484  0.005224348  0.052137019  0.0067800486
      t0           0.002287853 -0.009091644  0.006780049  0.0477904171
      
      
      $a_half
                       [,1]
      v           0.4941231
      v_lMd       2.5931756
      sv_lMTRUE   1.4684537
      B           0.8437531
      B_Eneutral  0.6261293
      B_Eaccuracy 0.8155694
      B_lRright   2.2752573
      A           0.7726374
      t0          0.3149199
      
      $epsilon
           [,1]
      as1t  1.5
      bd6t  1.5
      
      $origin
           [,1]
      as1t    2
      bd6t    2
      
      $alpha
      , , 1
      
                        as1t       bd6t
      v           -0.8659366 -0.8274847
      v_lMd        0.8828966  0.5761732
      sv_lMTRUE    1.4190220  1.7831332
      B           -1.5914170 -1.6746325
      B_Eneutral   0.2736373  0.5774843
      B_Eaccuracy -0.7775480 -1.8565900
      B_lRright    0.1379297  0.5908926
      A           -1.9607101 -2.4159174
      t0          -0.6633743 -0.9088031
      
      
      $stage
      [1] "init"
      
      $subj_ll
                 [,1]
      as1t -12583.929
      bd6t  -6090.795
      
      $last_theta_var_inv
                 [,1]      [,2]       [,3]      [,4]       [,5]      [,6]
       [1,] 26.028515  6.224107 -1.3572503 -3.811330   3.070339 -3.657354
       [2,]  6.224107 22.834612  1.1572385 -1.918507  -6.426290 -1.390907
       [3,] -1.357250  1.157238 27.2973400 -2.546630  -9.560950  9.393468
       [4,] -3.811330 -1.918507 -2.5466296 25.241309   4.568563  1.562002
       [5,]  3.070339 -6.426290 -9.5609496  4.568563  46.584345  3.445213
       [6,] -3.657354 -1.390907  9.3934682  1.562002   3.445213 22.028224
       [7,]  8.055800  5.016142 -2.0467420 -8.084231 -16.170739  7.104309
       [8,]  3.587576  8.795323 -1.3430047  1.224921  -4.769682 -1.611160
       [9,] 11.470411  3.049327 -0.4123182 -4.351590  -8.494190 -1.804338
                    [,7]        [,8]       [,9]
       [1,]   8.05580010  3.58757606 11.4704105
       [2,]   5.01614235  8.79532339  3.0493269
       [3,]  -2.04674198 -1.34300467 -0.4123182
       [4,]  -8.08423140  1.22492096 -4.3515897
       [5,] -16.17073884 -4.76968243 -8.4941901
       [6,]   7.10430863 -1.61116037 -1.8043382
       [7,]  44.63330665 -0.03626134 14.1866918
       [8,]  -0.03626134 23.94575450 -1.1814287
       [9,]  14.18669184 -1.18142869 30.0158974
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_LBA, design_LBA, n_trials = 10)
    Output
          subjects        E     S trials     R        rt
      1       as1t    speed  left      1  left 1.4536265
      5       as1t  neutral  left      1  left 2.1604799
      9       as1t accuracy  left      1  left 1.4754803
      13      as1t    speed right      1 right 0.9325650
      17      as1t  neutral right      1 right 1.6152186
      21      as1t accuracy right      1  left 1.6191709
      25      as1t    speed  left      2  left 0.9661781
      29      as1t  neutral  left      2  left 2.2083721
      33      as1t accuracy  left      2  left 1.5309142
      37      as1t    speed right      2 right 1.3213516
      41      as1t  neutral right      2 right 1.1851660
      45      as1t accuracy right      2 right 1.3201569
      49      as1t    speed  left      3  left 1.2544948
      53      as1t  neutral  left      3  left 2.0300890
      57      as1t accuracy  left      3  left 1.9658718
      61      as1t    speed right      3 right 2.1462549
      65      as1t  neutral right      3 right 2.0698275
      69      as1t accuracy right      3  left 1.4650098
      73      as1t    speed  left      4 right 1.0085966
      77      as1t  neutral  left      4  left 2.0493509
      81      as1t accuracy  left      4  left 2.0201509
      85      as1t    speed right      4 right 1.7347319
      89      as1t  neutral right      4 right 1.6125577
      93      as1t accuracy right      4 right 1.2566897
      97      as1t    speed  left      5 right 0.8382184
      101     as1t  neutral  left      5  left 2.5712633
      105     as1t accuracy  left      5  left 2.0562983
      109     as1t    speed right      5 right 1.0963985
      113     as1t  neutral right      5  left 1.5559065
      117     as1t accuracy right      5 right 1.1315209
      121     as1t    speed  left      6  left 1.2454877
      125     as1t  neutral  left      6  left 1.9603786
      129     as1t accuracy  left      6  left 1.2872130
      133     as1t    speed right      6 right 0.8971890
      137     as1t  neutral right      6 right 1.9179483
      141     as1t accuracy right      6 right 1.3965247
      145     as1t    speed  left      7 right 1.7369445
      149     as1t  neutral  left      7  left 1.9432939
      153     as1t accuracy  left      7  left 1.3417835
      157     as1t    speed right      7 right 1.0805940
      161     as1t  neutral right      7 right 2.0428886
      165     as1t accuracy right      7 right 1.1145099
      169     as1t    speed  left      8  left 1.7776611
      173     as1t  neutral  left      8 right 1.3409702
      177     as1t accuracy  left      8  left 1.7187141
      181     as1t    speed right      8 right 1.5532235
      185     as1t  neutral right      8 right 1.9012293
      189     as1t accuracy right      8 right 1.7998848
      193     as1t    speed  left      9  left 1.1197374
      197     as1t  neutral  left      9 right 1.3913653
      201     as1t accuracy  left      9 right 1.7521075
      205     as1t    speed right      9 right 0.9537609
      209     as1t  neutral right      9  left 0.9867002
      213     as1t accuracy right      9 right 2.2070464
      217     as1t    speed  left     10  left 1.2442672
      221     as1t  neutral  left     10  left 1.7141112
      225     as1t accuracy  left     10  left 1.3119795
      229     as1t    speed right     10 right 1.1031988
      233     as1t  neutral right     10 right 1.3273941
      237     as1t accuracy right     10 right 1.4690808
      3       bd6t    speed  left      1  left 1.1333569
      7       bd6t  neutral  left      1 right 1.7604053
      11      bd6t accuracy  left      1  left 1.4731090
      15      bd6t    speed right      1  left 0.8214118
      19      bd6t  neutral right      1 right 1.7271047
      23      bd6t accuracy right      1  left 1.1430470
      27      bd6t    speed  left      2  left 1.1887117
      31      bd6t  neutral  left      2 right 1.6976806
      35      bd6t accuracy  left      2  left 1.7163586
      39      bd6t    speed right      2  left 0.9247969
      43      bd6t  neutral right      2 right 1.7739821
      47      bd6t accuracy right      2 right 1.3156174
      51      bd6t    speed  left      3  left 0.9284531
      55      bd6t  neutral  left      3  left 1.6511994
      59      bd6t accuracy  left      3  left 1.2253081
      63      bd6t    speed right      3 right 1.2457512
      67      bd6t  neutral right      3 right 2.0139053
      71      bd6t accuracy right      3 right 1.8417664
      75      bd6t    speed  left      4  left 1.0972716
      79      bd6t  neutral  left      4  left 1.7693049
      83      bd6t accuracy  left      4  left 1.4958240
      87      bd6t    speed right      4 right 1.2045694
      91      bd6t  neutral right      4  left 1.5566072
      95      bd6t accuracy right      4 right 1.4363973
      99      bd6t    speed  left      5 right 1.1963245
      103     bd6t  neutral  left      5  left 1.4490702
      107     bd6t accuracy  left      5  left 1.3768588
      111     bd6t    speed right      5 right 1.6477356
      115     bd6t  neutral right      5  left 1.2180042
      119     bd6t accuracy right      5  left 2.1264373
      123     bd6t    speed  left      6  left 1.0088861
      127     bd6t  neutral  left      6 right 2.2457813
      131     bd6t accuracy  left      6 right 1.3415882
      135     bd6t    speed right      6 right 1.1266162
      139     bd6t  neutral right      6 right 2.2688797
      143     bd6t accuracy right      6 right 1.7039162
      147     bd6t    speed  left      7 right 1.1217635
      151     bd6t  neutral  left      7 right 2.8032292
      155     bd6t accuracy  left      7  left 1.8483821
      159     bd6t    speed right      7 right 1.4166983
      163     bd6t  neutral right      7  left 1.9382056
      167     bd6t accuracy right      7 right 1.3040573
      171     bd6t    speed  left      8 right 1.2803645
      175     bd6t  neutral  left      8  left 2.1822864
      179     bd6t accuracy  left      8  left 1.4169982
      183     bd6t    speed right      8  left 0.9983440
      187     bd6t  neutral right      8  left 1.9582264
      191     bd6t accuracy right      8 right 1.5126236
      195     bd6t    speed  left      9 right 1.3098783
      199     bd6t  neutral  left      9  left 1.4615988
      203     bd6t accuracy  left      9  left 2.2226335
      207     bd6t    speed right      9 right 1.1978558
      211     bd6t  neutral right      9 right 1.6111760
      215     bd6t accuracy right      9 right 1.6594657
      219     bd6t    speed  left     10 right 0.8109034
      223     bd6t  neutral  left     10  left 1.3506600
      227     bd6t accuracy  left     10  left 1.5620442
      231     bd6t    speed right     10  left 0.9063811
      235     bd6t  neutral right     10 right 1.3352358
      239     bd6t accuracy right     10 right 1.3213549

# RDM

    Code
      init_chains(RDM_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                        [,1]
      v           -0.9685927
      v_lMd        0.7061091
      s_lMTRUE     1.4890213
      B           -1.8150926
      B_Eneutral   0.3304096
      B_Eaccuracy -1.1421557
      B_lRright    0.1571934
      A           -2.0654072
      t0          -0.4405469
      
      $theta_var
      , , 1
      
                             v         v_lMd     s_lMTRUE             B   B_Eneutral
      v            0.059329585 -0.0104497238 -0.010774234  0.0021450832 -0.018128601
      v_lMd       -0.010449724  0.0559162133 -0.002872401  0.0014127637  0.004728529
      s_lMTRUE    -0.010774234 -0.0028724006  0.059773323  0.0087906529  0.022576320
      B            0.002145083  0.0014127637  0.008790653  0.0447136854  0.001458968
      B_Eneutral  -0.018128601  0.0047285288  0.022576320  0.0014589675  0.038522063
      B_Eaccuracy  0.018325623  0.0019110134 -0.037503336 -0.0100477306 -0.023423360
      B_lRright   -0.012850266 -0.0030964417  0.020333128  0.0098693855  0.019288722
      A           -0.009137082 -0.0181065189  0.007641907 -0.0029721272  0.008676374
      t0          -0.019762600  0.0006819411  0.001329633  0.0006671515  0.007687211
                   B_Eaccuracy    B_lRright            A            t0
      v            0.018325623 -0.012850266 -0.009137082 -0.0197625997
      v_lMd        0.001911013 -0.003096442 -0.018106519  0.0006819411
      s_lMTRUE    -0.037503336  0.020333128  0.007641907  0.0013296327
      B           -0.010047731  0.009869385 -0.002972127  0.0006671515
      B_Eneutral  -0.023423360  0.019288722  0.008676374  0.0076872108
      B_Eaccuracy  0.078051779 -0.028702651 -0.004381484  0.0022878527
      B_lRright   -0.028702651  0.042243108  0.005224348 -0.0090916441
      A           -0.004381484  0.005224348  0.052137019  0.0067800486
      t0           0.002287853 -0.009091644  0.006780049  0.0477904171
      
      
      $a_half
                       [,1]
      v           0.4941231
      v_lMd       2.5931756
      s_lMTRUE    1.4684537
      B           0.8437531
      B_Eneutral  0.6261293
      B_Eaccuracy 0.8155694
      B_lRright   2.2752573
      A           0.7726374
      t0          0.3149199
      
      $epsilon
           [,1]
      as1t  1.5
      bd6t  1.5
      
      $origin
           [,1]
      as1t    2
      bd6t    2
      
      $alpha
      , , 1
      
                        as1t       bd6t
      v           -0.8659366 -0.8274847
      v_lMd        0.8828966  0.5761732
      s_lMTRUE     1.4190220  1.7831332
      B           -1.5914170 -1.6746325
      B_Eneutral   0.2736373  0.5774843
      B_Eaccuracy -0.7775480 -1.8565900
      B_lRright    0.1379297  0.5908926
      A           -1.9607101 -2.4159174
      t0          -0.6633743 -0.9088031
      
      
      $stage
      [1] "init"
      
      $subj_ll
                 [,1]
      as1t -11878.967
      bd6t  -6641.485
      
      $last_theta_var_inv
                 [,1]      [,2]       [,3]      [,4]       [,5]      [,6]
       [1,] 26.028515  6.224107 -1.3572503 -3.811330   3.070339 -3.657354
       [2,]  6.224107 22.834612  1.1572385 -1.918507  -6.426290 -1.390907
       [3,] -1.357250  1.157238 27.2973400 -2.546630  -9.560950  9.393468
       [4,] -3.811330 -1.918507 -2.5466296 25.241309   4.568563  1.562002
       [5,]  3.070339 -6.426290 -9.5609496  4.568563  46.584345  3.445213
       [6,] -3.657354 -1.390907  9.3934682  1.562002   3.445213 22.028224
       [7,]  8.055800  5.016142 -2.0467420 -8.084231 -16.170739  7.104309
       [8,]  3.587576  8.795323 -1.3430047  1.224921  -4.769682 -1.611160
       [9,] 11.470411  3.049327 -0.4123182 -4.351590  -8.494190 -1.804338
                    [,7]        [,8]       [,9]
       [1,]   8.05580010  3.58757606 11.4704105
       [2,]   5.01614235  8.79532339  3.0493269
       [3,]  -2.04674198 -1.34300467 -0.4123182
       [4,]  -8.08423140  1.22492096 -4.3515897
       [5,] -16.17073884 -4.76968243 -8.4941901
       [6,]   7.10430863 -1.61116037 -1.8043382
       [7,]  44.63330665 -0.03626134 14.1866918
       [8,]  -0.03626134 23.94575450 -1.1814287
       [9,]  14.18669184 -1.18142869 30.0158974
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_RDM, design_RDM, n_trials = 10)
    Output
          subjects        E     S trials     R        rt
      1       as1t    speed  left      1  left 0.9182619
      5       as1t  neutral  left      1  left 1.3568147
      9       as1t accuracy  left      1  left 1.2535025
      13      as1t    speed right      1 right 1.3595822
      17      as1t  neutral right      1 right 1.6907934
      21      as1t accuracy right      1 right 1.5532023
      25      as1t    speed  left      2  left 0.8444087
      29      as1t  neutral  left      2  left 1.3393144
      33      as1t accuracy  left      2  left 1.3358297
      37      as1t    speed right      2 right 0.9862341
      41      as1t  neutral right      2 right 0.9799580
      45      as1t accuracy right      2 right 0.9630661
      49      as1t    speed  left      3  left 1.0245223
      53      as1t  neutral  left      3  left 1.4129907
      57      as1t accuracy  left      3  left 1.1492037
      61      as1t    speed right      3 right 0.7804297
      65      as1t  neutral right      3 right 1.4033613
      69      as1t accuracy right      3 right 0.9903596
      73      as1t    speed  left      4  left 0.9493701
      77      as1t  neutral  left      4  left 1.4143090
      81      as1t accuracy  left      4  left 1.1163854
      85      as1t    speed right      4 right 0.8354452
      89      as1t  neutral right      4 right 1.2111328
      93      as1t accuracy right      4 right 0.9976115
      97      as1t    speed  left      5 right 0.7434855
      101     as1t  neutral  left      5  left 1.2962211
      105     as1t accuracy  left      5  left 1.0891922
      109     as1t    speed right      5 right 0.7598303
      113     as1t  neutral right      5 right 1.3872515
      117     as1t accuracy right      5  left 1.6969225
      121     as1t    speed  left      6  left 1.0320553
      125     as1t  neutral  left      6  left 1.4913216
      129     as1t accuracy  left      6  left 0.9299625
      133     as1t    speed right      6 right 0.8011501
      137     as1t  neutral right      6 right 1.4199055
      141     as1t accuracy right      6 right 1.4221061
      145     as1t    speed  left      7  left 0.8308089
      149     as1t  neutral  left      7  left 1.4731482
      153     as1t accuracy  left      7  left 1.4730117
      157     as1t    speed right      7 right 0.8832382
      161     as1t  neutral right      7 right 1.5159172
      165     as1t accuracy right      7 right 1.7815679
      169     as1t    speed  left      8  left 1.2151523
      173     as1t  neutral  left      8 right 1.2979840
      177     as1t accuracy  left      8  left 1.2123689
      181     as1t    speed right      8 right 1.0418076
      185     as1t  neutral right      8 right 1.5071605
      189     as1t accuracy right      8 right 1.2041380
      193     as1t    speed  left      9  left 0.9093217
      197     as1t  neutral  left      9  left 1.5675691
      201     as1t accuracy  left      9  left 1.3694578
      205     as1t    speed right      9 right 0.8695812
      209     as1t  neutral right      9  left 0.9882402
      213     as1t accuracy right      9  left 1.0016315
      217     as1t    speed  left     10  left 1.0235876
      221     as1t  neutral  left     10  left 1.4560284
      225     as1t accuracy  left     10  left 1.1606575
      229     as1t    speed right     10 right 1.1462631
      233     as1t  neutral right     10 right 2.0040101
      237     as1t accuracy right     10 right 1.3591102
      3       bd6t    speed  left      1  left 1.0111386
      7       bd6t  neutral  left      1  left 1.1874563
      11      bd6t accuracy  left      1  left 1.3713312
      15      bd6t    speed right      1  left 0.7966629
      19      bd6t  neutral right      1 right 1.3892778
      23      bd6t accuracy right      1 right 1.5511251
      27      bd6t    speed  left      2 right 1.0028353
      31      bd6t  neutral  left      2  left 1.5231490
      35      bd6t accuracy  left      2  left 1.1938073
      39      bd6t    speed right      2 right 1.0097182
      43      bd6t  neutral right      2 right 1.2939758
      47      bd6t accuracy right      2 right 1.5001656
      51      bd6t    speed  left      3  left 1.4068085
      55      bd6t  neutral  left      3  left 1.2890259
      59      bd6t accuracy  left      3  left 1.6137723
      63      bd6t    speed right      3 right 0.9729191
      67      bd6t  neutral right      3 right 1.4360370
      71      bd6t accuracy right      3 right 1.1661584
      75      bd6t    speed  left      4  left 0.7799093
      79      bd6t  neutral  left      4  left 1.5564834
      83      bd6t accuracy  left      4  left 1.2279961
      87      bd6t    speed right      4 right 1.0165810
      91      bd6t  neutral right      4  left 1.6377343
      95      bd6t accuracy right      4 right 1.3911158
      99      bd6t    speed  left      5  left 1.1696183
      103     bd6t  neutral  left      5  left 1.2463580
      107     bd6t accuracy  left      5  left 1.0579920
      111     bd6t    speed right      5 right 1.1977988
      115     bd6t  neutral right      5 right 1.5030490
      119     bd6t accuracy right      5 right 1.6128328
      123     bd6t    speed  left      6  left 1.2759126
      127     bd6t  neutral  left      6  left 1.7988705
      131     bd6t accuracy  left      6  left 1.4098915
      135     bd6t    speed right      6  left 0.9339767
      139     bd6t  neutral right      6 right 1.6111036
      143     bd6t accuracy right      6 right 1.2261954
      147     bd6t    speed  left      7 right 0.8399248
      151     bd6t  neutral  left      7  left 1.7582247
      155     bd6t accuracy  left      7  left 1.1626884
      159     bd6t    speed right      7 right 0.9489773
      163     bd6t  neutral right      7 right 1.2652466
      167     bd6t accuracy right      7 right 1.5137636
      171     bd6t    speed  left      8  left 0.9131276
      175     bd6t  neutral  left      8  left 1.3736444
      179     bd6t accuracy  left      8  left 1.0281685
      183     bd6t    speed right      8 right 1.0191257
      187     bd6t  neutral right      8 right 1.6076562
      191     bd6t accuracy right      8 right 1.2395016
      195     bd6t    speed  left      9  left 0.9003089
      199     bd6t  neutral  left      9  left 1.2843223
      203     bd6t accuracy  left      9  left 1.5685607
      207     bd6t    speed right      9 right 1.0582006
      211     bd6t  neutral right      9 right 1.3990724
      215     bd6t accuracy right      9 right 1.2194096
      219     bd6t    speed  left     10  left 1.3899057
      223     bd6t  neutral  left     10  left 1.2065876
      227     bd6t accuracy  left     10  left 1.2870212
      231     bd6t    speed right     10 right 0.8595019
      235     bd6t  neutral right     10 right 1.0803544
      239     bd6t accuracy right     10  left 1.4520010

# DDM

    Code
      init_chains(DDM_s, particles = 10, cores_per_chain = 1)[[1]]$samples
    Output
      $theta_mu
                        [,1]
      v_Sleft     -0.9685927
      v_Sright     0.7061091
      a            1.4890213
      a_Eneutral  -1.8150926
      a_Eaccuracy  0.3304096
      t0          -1.1421557
      Z            0.1571934
      sv          -2.0654072
      SZ          -0.4405469
      
      $theta_var
      , , 1
      
                       v_Sleft      v_Sright            a    a_Eneutral  a_Eaccuracy
      v_Sleft      0.059329585 -0.0104497238 -0.010774234  0.0021450832 -0.018128601
      v_Sright    -0.010449724  0.0559162133 -0.002872401  0.0014127637  0.004728529
      a           -0.010774234 -0.0028724006  0.059773323  0.0087906529  0.022576320
      a_Eneutral   0.002145083  0.0014127637  0.008790653  0.0447136854  0.001458968
      a_Eaccuracy -0.018128601  0.0047285288  0.022576320  0.0014589675  0.038522063
      t0           0.018325623  0.0019110134 -0.037503336 -0.0100477306 -0.023423360
      Z           -0.012850266 -0.0030964417  0.020333128  0.0098693855  0.019288722
      sv          -0.009137082 -0.0181065189  0.007641907 -0.0029721272  0.008676374
      SZ          -0.019762600  0.0006819411  0.001329633  0.0006671515  0.007687211
                            t0            Z           sv            SZ
      v_Sleft      0.018325623 -0.012850266 -0.009137082 -0.0197625997
      v_Sright     0.001911013 -0.003096442 -0.018106519  0.0006819411
      a           -0.037503336  0.020333128  0.007641907  0.0013296327
      a_Eneutral  -0.010047731  0.009869385 -0.002972127  0.0006671515
      a_Eaccuracy -0.023423360  0.019288722  0.008676374  0.0076872108
      t0           0.078051779 -0.028702651 -0.004381484  0.0022878527
      Z           -0.028702651  0.042243108  0.005224348 -0.0090916441
      sv          -0.004381484  0.005224348  0.052137019  0.0067800486
      SZ           0.002287853 -0.009091644  0.006780049  0.0477904171
      
      
      $a_half
                       [,1]
      v_Sleft     0.4941231
      v_Sright    2.5931756
      a           1.4684537
      a_Eneutral  0.8437531
      a_Eaccuracy 0.6261293
      t0          0.8155694
      Z           2.2752573
      sv          0.7726374
      SZ          0.3149199
      
      $epsilon
           [,1]
      as1t  1.5
      bd6t  1.5
      
      $origin
           [,1]
      as1t    2
      bd6t    2
      
      $alpha
      , , 1
      
                          as1t        bd6t
      v_Sleft     -0.795745564 -0.61298014
      v_Sright     0.831296796  0.53982132
      a            1.424407847  1.19546709
      a_Eneutral  -1.922987294 -1.62502932
      a_Eaccuracy  0.049542401  0.03054263
      t0          -1.396100362 -1.38269638
      Z           -0.005054203 -0.21436007
      sv          -2.197044667 -2.08693786
      SZ          -0.436668005 -0.46579217
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -3255.245
      bd6t -2695.119
      
      $last_theta_var_inv
                 [,1]      [,2]       [,3]      [,4]       [,5]      [,6]
       [1,] 26.028515  6.224107 -1.3572503 -3.811330   3.070339 -3.657354
       [2,]  6.224107 22.834612  1.1572385 -1.918507  -6.426290 -1.390907
       [3,] -1.357250  1.157238 27.2973400 -2.546630  -9.560950  9.393468
       [4,] -3.811330 -1.918507 -2.5466296 25.241309   4.568563  1.562002
       [5,]  3.070339 -6.426290 -9.5609496  4.568563  46.584345  3.445213
       [6,] -3.657354 -1.390907  9.3934682  1.562002   3.445213 22.028224
       [7,]  8.055800  5.016142 -2.0467420 -8.084231 -16.170739  7.104309
       [8,]  3.587576  8.795323 -1.3430047  1.224921  -4.769682 -1.611160
       [9,] 11.470411  3.049327 -0.4123182 -4.351590  -8.494190 -1.804338
                    [,7]        [,8]       [,9]
       [1,]   8.05580010  3.58757606 11.4704105
       [2,]   5.01614235  8.79532339  3.0493269
       [3,]  -2.04674198 -1.34300467 -0.4123182
       [4,]  -8.08423140  1.22492096 -4.3515897
       [5,] -16.17073884 -4.76968243 -8.4941901
       [6,]   7.10430863 -1.61116037 -1.8043382
       [7,]  44.63330665 -0.03626134 14.1866918
       [8,]  -0.03626134 23.94575450 -1.1814287
       [9,]  14.18669184 -1.18142869 30.0158974
      
      $idx
      [1] 1
      

---

    Code
      make_data(p_DDM, design_DDM, n_trials = 10)
    Output
          subjects        E     S trials     R        rt
      1       as1t    speed  left      1  left 0.6704792
      3       as1t  neutral  left      1  left 0.9907978
      5       as1t accuracy  left      1 right 2.2457532
      7       as1t    speed right      1 right 0.4853054
      9       as1t  neutral right      1 right 0.3740953
      11      as1t accuracy right      1 right 0.8545565
      13      as1t    speed  left      2 right 0.3483702
      15      as1t  neutral  left      2  left 0.7892018
      17      as1t accuracy  left      2  left 0.5317843
      19      as1t    speed right      2 right 0.3927556
      21      as1t  neutral right      2 right 0.7891189
      23      as1t accuracy right      2 right 1.1108798
      25      as1t    speed  left      3 right 0.8734347
      27      as1t  neutral  left      3  left 0.3471626
      29      as1t accuracy  left      3  left 0.5490293
      31      as1t    speed right      3 right 0.5949942
      33      as1t  neutral right      3 right 2.1445549
      35      as1t accuracy right      3 right 0.4910811
      37      as1t    speed  left      4  left 0.4401786
      39      as1t  neutral  left      4  left 1.3214967
      41      as1t accuracy  left      4  left 0.3699592
      43      as1t    speed right      4 right 0.3747106
      45      as1t  neutral right      4 right 0.8956814
      47      as1t accuracy right      4  left 1.3078173
      49      as1t    speed  left      5 right 0.3206229
      51      as1t  neutral  left      5 right 0.5461436
      53      as1t accuracy  left      5  left 0.5585949
      55      as1t    speed right      5 right 0.3047198
      57      as1t  neutral right      5 right 0.3249469
      59      as1t accuracy right      5 right 0.4450077
      61      as1t    speed  left      6  left 0.3549161
      63      as1t  neutral  left      6  left 0.8488631
      65      as1t accuracy  left      6  left 0.5425459
      67      as1t    speed right      6 right 0.3786981
      69      as1t  neutral right      6  left 0.5512905
      71      as1t accuracy right      6 right 0.5610937
      73      as1t    speed  left      7  left 0.6292160
      75      as1t  neutral  left      7  left 0.7881852
      77      as1t accuracy  left      7  left 0.6770581
      79      as1t    speed right      7  left 0.4285827
      81      as1t  neutral right      7  left 1.0912666
      83      as1t accuracy right      7  left 2.2791974
      85      as1t    speed  left      8  left 0.4986977
      87      as1t  neutral  left      8  left 0.8194816
      89      as1t accuracy  left      8  left 0.3988712
      91      as1t    speed right      8  left 0.3570129
      93      as1t  neutral right      8 right 0.7220455
      95      as1t accuracy right      8  left 4.9173973
      97      as1t    speed  left      9  left 0.4635973
      99      as1t  neutral  left      9 right 2.3802724
      101     as1t accuracy  left      9 right 0.6767897
      103     as1t    speed right      9  left 0.6815723
      105     as1t  neutral right      9 right 0.2817931
      107     as1t accuracy right      9 right 0.4405847
      109     as1t    speed  left     10  left 0.2586958
      111     as1t  neutral  left     10 right 0.9079815
      113     as1t accuracy  left     10  left 0.9043604
      115     as1t    speed right     10 right 0.3517752
      117     as1t  neutral right     10  left 1.6029021
      119     as1t accuracy right     10  left 0.4164281
      2       bd6t    speed  left      1 right 0.3558950
      4       bd6t  neutral  left      1  left 0.4383693
      6       bd6t accuracy  left      1  left 2.6004300
      8       bd6t    speed right      1  left 0.9677727
      10      bd6t  neutral right      1 right 0.3963958
      12      bd6t accuracy right      1  left 0.4102032
      14      bd6t    speed  left      2  left 0.3596023
      16      bd6t  neutral  left      2  left 0.4515037
      18      bd6t accuracy  left      2  left 0.3322909
      20      bd6t    speed right      2  left 1.1910805
      22      bd6t  neutral right      2  left 0.9744928
      24      bd6t accuracy right      2 right 0.7199995
      26      bd6t    speed  left      3  left 0.7249741
      28      bd6t  neutral  left      3  left 0.3614126
      30      bd6t accuracy  left      3 right 0.4578476
      32      bd6t    speed right      3 right 0.6345737
      34      bd6t  neutral right      3  left 0.5119013
      36      bd6t accuracy right      3 right 2.4335010
      38      bd6t    speed  left      4  left 0.4028299
      40      bd6t  neutral  left      4  left 0.3371061
      42      bd6t accuracy  left      4  left 0.9431179
      44      bd6t    speed right      4  left 0.5372962
      46      bd6t  neutral right      4 right 0.5863872
      48      bd6t accuracy right      4 right 0.5787578
      50      bd6t    speed  left      5 right 1.6795988
      52      bd6t  neutral  left      5  left 0.4688466
      54      bd6t accuracy  left      5  left 0.4813896
      56      bd6t    speed right      5 right 0.9888955
      58      bd6t  neutral right      5 right 0.5648297
      60      bd6t accuracy right      5 right 0.3892264
      62      bd6t    speed  left      6  left 1.0679385
      64      bd6t  neutral  left      6 right 0.6659524
      66      bd6t accuracy  left      6  left 0.3454326
      68      bd6t    speed right      6 right 0.3612517
      70      bd6t  neutral right      6 right 1.0799850
      72      bd6t accuracy right      6 right 1.6119733
      74      bd6t    speed  left      7  left 0.3320072
      76      bd6t  neutral  left      7 right 1.0442965
      78      bd6t accuracy  left      7  left 0.3013069
      80      bd6t    speed right      7  left 1.1229204
      82      bd6t  neutral right      7 right 7.5471825
      84      bd6t accuracy right      7  left 0.7837445
      86      bd6t    speed  left      8  left 1.4620092
      88      bd6t  neutral  left      8  left 0.5405186
      90      bd6t accuracy  left      8  left 2.8014449
      92      bd6t    speed right      8 right 0.3434099
      94      bd6t  neutral right      8 right 0.3940605
      96      bd6t accuracy right      8 right 0.5348753
      98      bd6t    speed  left      9  left 0.3384052
      100     bd6t  neutral  left      9 right 0.6003083
      102     bd6t accuracy  left      9  left 1.0441463
      104     bd6t    speed right      9 right 0.3687537
      106     bd6t  neutral right      9 right 0.8069842
      108     bd6t accuracy right      9 right 0.6772552
      110     bd6t    speed  left     10  left 0.6862207
      112     bd6t  neutral  left     10  left 0.5237003
      114     bd6t accuracy  left     10  left 1.1545852
      116     bd6t    speed right     10  left 0.8844238
      118     bd6t  neutral right     10  left 1.4948986
      120     bd6t accuracy right     10  left 0.7836393

