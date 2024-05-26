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
      as1t -1012.1654
      bd6t  -828.2889
      
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
      make_data(p_LNR, design_LNR, trials = 10)
    Output
          subjects        E     S trials     R        rt
      1       as1t    speed  left      1 right 0.6621568
      5       as1t  neutral  left      1 right 0.2063691
      9       as1t accuracy  left      1 right 0.2163960
      13      as1t    speed right      1  left 0.2149467
      17      as1t  neutral right      1  left 0.2208144
      21      as1t accuracy right      1  left 0.2962805
      25      as1t    speed  left      2 right 1.3993717
      29      as1t  neutral  left      2 right 0.2338689
      33      as1t accuracy  left      2  left 1.2707757
      37      as1t    speed right      2  left 0.2347907
      41      as1t  neutral right      2  left 1.9690862
      45      as1t accuracy right      2  left 0.2072654
      49      as1t    speed  left      3 right 0.2289023
      53      as1t  neutral  left      3 right 0.2159625
      57      as1t accuracy  left      3 right 0.3556382
      61      as1t    speed right      3  left 0.2014666
      65      as1t  neutral right      3 right 0.2407919
      69      as1t accuracy right      3 right 0.8895984
      73      as1t    speed  left      4 right 0.8896424
      77      as1t  neutral  left      4 right 0.2841160
      81      as1t accuracy  left      4 right 0.2082286
      85      as1t    speed right      4  left 0.3757464
      89      as1t  neutral right      4  left 0.3402543
      93      as1t accuracy right      4  left 0.3885063
      97      as1t    speed  left      5  left 0.3268648
      101     as1t  neutral  left      5 right 0.2931404
      105     as1t accuracy  left      5 right 0.5631042
      109     as1t    speed right      5  left 0.6266161
      113     as1t  neutral right      5  left 0.2576004
      117     as1t accuracy right      5  left 0.2724988
      121     as1t    speed  left      6 right 0.5078265
      125     as1t  neutral  left      6  left 1.1182128
      129     as1t accuracy  left      6 right 0.2151060
      133     as1t    speed right      6 right 0.3324145
      137     as1t  neutral right      6  left 0.2893369
      141     as1t accuracy right      6  left 0.6245632
      145     as1t    speed  left      7 right 0.3190899
      149     as1t  neutral  left      7  left 0.5953184
      153     as1t accuracy  left      7 right 0.2044777
      157     as1t    speed right      7  left 0.2606403
      161     as1t  neutral right      7 right 0.2842899
      165     as1t accuracy right      7  left 0.2101295
      169     as1t    speed  left      8 right 5.6647179
      173     as1t  neutral  left      8 right 0.3104341
      177     as1t accuracy  left      8 right 0.2299710
      181     as1t    speed right      8  left 0.2152846
      185     as1t  neutral right      8 right 0.2591430
      189     as1t accuracy right      8  left 1.0453442
      193     as1t    speed  left      9 right 0.2800495
      197     as1t  neutral  left      9 right 0.2070659
      201     as1t accuracy  left      9 right 0.5368400
      205     as1t    speed right      9  left 0.2216829
      209     as1t  neutral right      9 right 0.3311753
      213     as1t accuracy right      9 right 0.3193368
      217     as1t    speed  left     10 right 0.2399517
      221     as1t  neutral  left     10 right 0.2520374
      225     as1t accuracy  left     10  left 1.8249842
      229     as1t    speed right     10 right 0.7272009
      233     as1t  neutral right     10  left 0.7467092
      237     as1t accuracy right     10 right 0.2285965
      3       bd6t    speed  left      1 right 0.3521819
      7       bd6t  neutral  left      1  left 2.8239520
      11      bd6t accuracy  left      1  left 0.2494104
      15      bd6t    speed right      1  left 0.6594934
      19      bd6t  neutral right      1  left 0.3739250
      23      bd6t accuracy right      1  left 0.3894119
      27      bd6t    speed  left      2 right 0.4399945
      31      bd6t  neutral  left      2 right 0.2237928
      35      bd6t accuracy  left      2 right 0.2111508
      39      bd6t    speed right      2 right 0.2986954
      43      bd6t  neutral right      2  left 0.4139748
      47      bd6t accuracy right      2  left 0.2233620
      51      bd6t    speed  left      3 right 0.3893381
      55      bd6t  neutral  left      3  left 0.5294351
      59      bd6t accuracy  left      3 right 0.2143891
      63      bd6t    speed right      3  left 0.3271981
      67      bd6t  neutral right      3 right 0.2989625
      71      bd6t accuracy right      3  left 2.6800823
      75      bd6t    speed  left      4 right 0.4127312
      79      bd6t  neutral  left      4 right 0.4247098
      83      bd6t accuracy  left      4  left 0.3207819
      87      bd6t    speed right      4 right 0.4601309
      91      bd6t  neutral right      4  left 0.2862879
      95      bd6t accuracy right      4  left 0.2566716
      99      bd6t    speed  left      5 right 0.4284081
      103     bd6t  neutral  left      5 right 0.2300322
      107     bd6t accuracy  left      5 right 0.2416458
      111     bd6t    speed right      5  left 0.2491193
      115     bd6t  neutral right      5  left 0.3075896
      119     bd6t accuracy right      5  left 0.2695769
      123     bd6t    speed  left      6 right 1.1131387
      127     bd6t  neutral  left      6 right 0.2122215
      131     bd6t accuracy  left      6 right 0.5294487
      135     bd6t    speed right      6 right 0.3498691
      139     bd6t  neutral right      6 right 0.3132170
      143     bd6t accuracy right      6  left 0.3200083
      147     bd6t    speed  left      7 right 0.5491461
      151     bd6t  neutral  left      7  left 0.2465168
      155     bd6t accuracy  left      7  left 0.3164743
      159     bd6t    speed right      7  left 0.4062331
      163     bd6t  neutral right      7  left 0.2521459
      167     bd6t accuracy right      7 right 0.5737047
      171     bd6t    speed  left      8 right 0.3196930
      175     bd6t  neutral  left      8 right 1.2783628
      179     bd6t accuracy  left      8 right 0.8892908
      183     bd6t    speed right      8  left 0.6911094
      187     bd6t  neutral right      8  left 0.2639333
      191     bd6t accuracy right      8  left 0.2970921
      195     bd6t    speed  left      9 right 0.4383953
      199     bd6t  neutral  left      9 right 2.1349298
      203     bd6t accuracy  left      9 right 0.6442377
      207     bd6t    speed right      9 right 0.2612303
      211     bd6t  neutral right      9 right 1.3892407
      215     bd6t accuracy right      9  left 0.2842021
      219     bd6t    speed  left     10 right 0.2211178
      223     bd6t  neutral  left     10 right 0.3433792
      227     bd6t accuracy  left     10  left 0.7444105
      231     bd6t    speed right     10  left 0.4036253
      235     bd6t  neutral right     10  left 0.3621805
      239     bd6t accuracy right     10  left 0.3439226

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
      as1t -13997.387
      bd6t  -7641.792
      
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
      make_data(p_LBA, design_LBA, trials = 10)
    Output
          subjects        E     S trials     R        rt
      1       as1t    speed  left      1 right  3.030501
      5       as1t  neutral  left      1  left  3.539124
      9       as1t accuracy  left      1 right  2.018035
      13      as1t    speed right      1 right  4.039522
      17      as1t  neutral right      1 right  3.096771
      21      as1t accuracy right      1 right  2.831589
      25      as1t    speed  left      2  left  1.582925
      29      as1t  neutral  left      2  left  2.274229
      33      as1t accuracy  left      2 right  1.846171
      37      as1t    speed right      2 right  1.889838
      41      as1t  neutral right      2 right  2.484168
      45      as1t accuracy right      2  left  3.046477
      49      as1t    speed  left      3  left  1.486228
      53      as1t  neutral  left      3  left  3.778649
      57      as1t accuracy  left      3 right  1.576586
      61      as1t    speed right      3 right  2.388316
      65      as1t  neutral right      3 right  3.756077
      69      as1t accuracy right      3  left  6.599845
      73      as1t    speed  left      4  left  1.745405
      77      as1t  neutral  left      4  left  3.839705
      81      as1t accuracy  left      4  left  3.845541
      85      as1t    speed right      4 right  1.060495
      89      as1t  neutral right      4 right  3.446446
      93      as1t accuracy right      4 right  1.908258
      97      as1t    speed  left      5  left  2.236360
      101     as1t  neutral  left      5  left  2.467514
      105     as1t accuracy  left      5 right  2.376394
      109     as1t    speed right      5 right  1.815742
      113     as1t  neutral right      5  left  6.137105
      117     as1t accuracy right      5 right  3.050245
      121     as1t    speed  left      6 right  4.654374
      125     as1t  neutral  left      6 right  3.737266
      129     as1t accuracy  left      6  left  3.261177
      133     as1t    speed right      6 right  2.419888
      137     as1t  neutral right      6 right  2.328470
      141     as1t accuracy right      6  left  1.912428
      145     as1t    speed  left      7 right  1.773037
      149     as1t  neutral  left      7  left  2.060795
      153     as1t accuracy  left      7  left  3.284590
      157     as1t    speed right      7 right  2.278090
      161     as1t  neutral right      7  left  3.814300
      165     as1t accuracy right      7 right  3.006519
      169     as1t    speed  left      8  left  2.730204
      173     as1t  neutral  left      8  left  9.251958
      177     as1t accuracy  left      8 right  5.571195
      181     as1t    speed right      8 right  1.697408
      185     as1t  neutral right      8  left  2.403360
      189     as1t accuracy right      8 right  1.772008
      193     as1t    speed  left      9  left  1.371370
      197     as1t  neutral  left      9  left  1.877307
      201     as1t accuracy  left      9 right  4.181847
      205     as1t    speed right      9 right  5.671736
      209     as1t  neutral right      9 right  6.015465
      213     as1t accuracy right      9  left  1.970632
      217     as1t    speed  left     10  left  1.920388
      221     as1t  neutral  left     10  left  2.253761
      225     as1t accuracy  left     10  left  1.770226
      229     as1t    speed right     10 right  1.801675
      233     as1t  neutral right     10 right  2.370456
      237     as1t accuracy right     10 right  1.745306
      3       bd6t    speed  left      1  left  2.212102
      7       bd6t  neutral  left      1  left  3.743711
      11      bd6t accuracy  left      1  left  3.324826
      15      bd6t    speed right      1  left  2.016301
      19      bd6t  neutral right      1 right  2.106579
      23      bd6t accuracy right      1  left  2.139603
      27      bd6t    speed  left      2  left  3.334463
      31      bd6t  neutral  left      2  left  2.029409
      35      bd6t accuracy  left      2 right  3.957159
      39      bd6t    speed right      2 right  2.117888
      43      bd6t  neutral right      2 right  3.021339
      47      bd6t accuracy right      2 right  2.316572
      51      bd6t    speed  left      3 right  3.056940
      55      bd6t  neutral  left      3  left  3.128322
      59      bd6t accuracy  left      3  left  3.181291
      63      bd6t    speed right      3  left  1.283469
      67      bd6t  neutral right      3 right  2.567409
      71      bd6t accuracy right      3 right  2.278658
      75      bd6t    speed  left      4  left  2.210870
      79      bd6t  neutral  left      4  left  2.744873
      83      bd6t accuracy  left      4  left  2.774626
      87      bd6t    speed right      4  left  1.480785
      91      bd6t  neutral right      4  left  7.277249
      95      bd6t accuracy right      4 right  1.548620
      99      bd6t    speed  left      5  left  5.784232
      103     bd6t  neutral  left      5 right  3.327742
      107     bd6t accuracy  left      5 right  3.104018
      111     bd6t    speed right      5  left  1.203665
      115     bd6t  neutral right      5  left  2.090526
      119     bd6t accuracy right      5  left  4.917617
      123     bd6t    speed  left      6  left  2.501372
      127     bd6t  neutral  left      6  left  5.301814
      131     bd6t accuracy  left      6  left  3.460860
      135     bd6t    speed right      6 right  2.318174
      139     bd6t  neutral right      6 right  3.773084
      143     bd6t accuracy right      6  left  3.505149
      147     bd6t    speed  left      7  left  1.799115
      151     bd6t  neutral  left      7  left  3.401848
      155     bd6t accuracy  left      7  left  2.444221
      159     bd6t    speed right      7 right  1.289936
      163     bd6t  neutral right      7  left  5.230063
      167     bd6t accuracy right      7 right  2.277541
      171     bd6t    speed  left      8  left  1.522399
      175     bd6t  neutral  left      8 right  1.849652
      179     bd6t accuracy  left      8  left  1.379855
      183     bd6t    speed right      8 right  1.790735
      187     bd6t  neutral right      8 right  3.215419
      191     bd6t accuracy right      8  left  2.678309
      195     bd6t    speed  left      9 right  1.449396
      199     bd6t  neutral  left      9  left  2.880126
      203     bd6t accuracy  left      9  left  2.541788
      207     bd6t    speed right      9  left  4.824451
      211     bd6t  neutral right      9 right  2.519785
      215     bd6t accuracy right      9 right  5.377238
      219     bd6t    speed  left     10  left  2.250957
      223     bd6t  neutral  left     10  left  2.622829
      227     bd6t accuracy  left     10 right  2.837247
      231     bd6t    speed right     10 right  1.353242
      235     bd6t  neutral right     10 right 10.014760
      239     bd6t accuracy right     10 right  2.474592

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
      as1t -13015.14
      bd6t  -8124.47
      
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
      make_data(p_RDM, design_RDM, trials = 10)
    Output
          subjects        E     S trials     R        rt
      1       as1t    speed  left      1  left 0.9280160
      5       as1t  neutral  left      1  left 1.3665688
      9       as1t accuracy  left      1  left 1.2632566
      13      as1t    speed right      1 right 1.3693363
      17      as1t  neutral right      1 right 1.7005475
      21      as1t accuracy right      1 right 1.5629564
      25      as1t    speed  left      2  left 0.8541628
      29      as1t  neutral  left      2  left 1.3490685
      33      as1t accuracy  left      2  left 1.3455839
      37      as1t    speed right      2 right 0.9959882
      41      as1t  neutral right      2 right 0.9897121
      45      as1t accuracy right      2 right 0.9728202
      49      as1t    speed  left      3  left 1.0342764
      53      as1t  neutral  left      3  left 1.4227448
      57      as1t accuracy  left      3  left 1.1589579
      61      as1t    speed right      3 right 0.7901838
      65      as1t  neutral right      3 right 1.4131155
      69      as1t accuracy right      3 right 1.0001137
      73      as1t    speed  left      4  left 0.9591242
      77      as1t  neutral  left      4  left 1.4240631
      81      as1t accuracy  left      4  left 1.1261396
      85      as1t    speed right      4 right 0.8451994
      89      as1t  neutral right      4 right 1.2208869
      93      as1t accuracy right      4 right 1.0073656
      97      as1t    speed  left      5 right 0.7532396
      101     as1t  neutral  left      5  left 1.3059752
      105     as1t accuracy  left      5  left 1.0989463
      109     as1t    speed right      5 right 0.7695844
      113     as1t  neutral right      5 right 1.3970057
      117     as1t accuracy right      5  left 1.7066766
      121     as1t    speed  left      6  left 1.0418094
      125     as1t  neutral  left      6  left 1.5010757
      129     as1t accuracy  left      6  left 0.9397166
      133     as1t    speed right      6 right 0.8109042
      137     as1t  neutral right      6 right 1.4296597
      141     as1t accuracy right      6 right 1.4318602
      145     as1t    speed  left      7  left 0.8405631
      149     as1t  neutral  left      7  left 1.4829023
      153     as1t accuracy  left      7  left 1.4827658
      157     as1t    speed right      7 right 0.8929923
      161     as1t  neutral right      7 right 1.5256713
      165     as1t accuracy right      7 right 1.7913220
      169     as1t    speed  left      8  left 1.2249064
      173     as1t  neutral  left      8 right 1.3077382
      177     as1t accuracy  left      8  left 1.2221230
      181     as1t    speed right      8 right 1.0515617
      185     as1t  neutral right      8 right 1.5169147
      189     as1t accuracy right      8 right 1.2138921
      193     as1t    speed  left      9  left 0.9190758
      197     as1t  neutral  left      9  left 1.5773232
      201     as1t accuracy  left      9  left 1.3792119
      205     as1t    speed right      9 right 0.8793353
      209     as1t  neutral right      9  left 0.9979943
      213     as1t accuracy right      9  left 1.0113856
      217     as1t    speed  left     10  left 1.0333417
      221     as1t  neutral  left     10  left 1.4657825
      225     as1t accuracy  left     10  left 1.1704116
      229     as1t    speed right     10 right 1.1560172
      233     as1t  neutral right     10 right 2.0137642
      237     as1t accuracy right     10 right 1.3688643
      3       bd6t    speed  left      1  left 1.0208927
      7       bd6t  neutral  left      1  left 1.1972104
      11      bd6t accuracy  left      1  left 1.3810853
      15      bd6t    speed right      1  left 0.8064170
      19      bd6t  neutral right      1 right 1.3990319
      23      bd6t accuracy right      1 right 1.5608792
      27      bd6t    speed  left      2 right 1.0125895
      31      bd6t  neutral  left      2  left 1.5329031
      35      bd6t accuracy  left      2  left 1.2035614
      39      bd6t    speed right      2 right 1.0194723
      43      bd6t  neutral right      2 right 1.3037299
      47      bd6t accuracy right      2 right 1.5099197
      51      bd6t    speed  left      3  left 1.4165626
      55      bd6t  neutral  left      3  left 1.2987800
      59      bd6t accuracy  left      3  left 1.6235264
      63      bd6t    speed right      3 right 0.9826732
      67      bd6t  neutral right      3 right 1.4457911
      71      bd6t accuracy right      3 right 1.1759125
      75      bd6t    speed  left      4  left 0.7896634
      79      bd6t  neutral  left      4  left 1.5662376
      83      bd6t accuracy  left      4  left 1.2377503
      87      bd6t    speed right      4 right 1.0263352
      91      bd6t  neutral right      4  left 1.6474885
      95      bd6t accuracy right      4 right 1.4008699
      99      bd6t    speed  left      5  left 1.1793724
      103     bd6t  neutral  left      5  left 1.2561121
      107     bd6t accuracy  left      5  left 1.0677462
      111     bd6t    speed right      5 right 1.2075529
      115     bd6t  neutral right      5 right 1.5128032
      119     bd6t accuracy right      5 right 1.6225869
      123     bd6t    speed  left      6  left 1.2856667
      127     bd6t  neutral  left      6  left 1.8086246
      131     bd6t accuracy  left      6  left 1.4196456
      135     bd6t    speed right      6  left 0.9437308
      139     bd6t  neutral right      6 right 1.6208577
      143     bd6t accuracy right      6 right 1.2359495
      147     bd6t    speed  left      7 right 0.8496789
      151     bd6t  neutral  left      7  left 1.7679788
      155     bd6t accuracy  left      7  left 1.1724425
      159     bd6t    speed right      7 right 0.9587314
      163     bd6t  neutral right      7 right 1.2750007
      167     bd6t accuracy right      7 right 1.5235177
      171     bd6t    speed  left      8  left 0.9228817
      175     bd6t  neutral  left      8  left 1.3833986
      179     bd6t accuracy  left      8  left 1.0379226
      183     bd6t    speed right      8 right 1.0288798
      187     bd6t  neutral right      8 right 1.6174103
      191     bd6t accuracy right      8 right 1.2492557
      195     bd6t    speed  left      9  left 0.9100630
      199     bd6t  neutral  left      9  left 1.2940764
      203     bd6t accuracy  left      9  left 1.5783148
      207     bd6t    speed right      9 right 1.0679547
      211     bd6t  neutral right      9 right 1.4088265
      215     bd6t accuracy right      9 right 1.2291637
      219     bd6t    speed  left     10  left 1.3996598
      223     bd6t  neutral  left     10  left 1.2163417
      227     bd6t accuracy  left     10  left 1.2967753
      231     bd6t    speed right     10 right 0.8692560
      235     bd6t  neutral right     10 right 1.0901085
      239     bd6t accuracy right     10  left 1.4617551

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
      as1t -3463.140
      bd6t -2944.067
      
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
      make_data(p_DDM, design_DDM, trials = 10)
    Output
          subjects        E     S trials     R        rt
      1       as1t    speed  left      1  left 0.4051153
      3       as1t  neutral  left      1  left 1.1555614
      5       as1t accuracy  left      1  left 0.5004580
      7       as1t    speed right      1 right 0.5486047
      9       as1t  neutral right      1 right 0.8940213
      11      as1t accuracy right      1 right 0.4702255
      13      as1t    speed  left      2  left 0.3723834
      15      as1t  neutral  left      2  left 1.8486886
      17      as1t accuracy  left      2  left 1.5727215
      19      as1t    speed right      2 right 0.4971404
      21      as1t  neutral right      2 right 0.6204149
      23      as1t accuracy right      2 right 0.5450486
      25      as1t    speed  left      3  left 1.4835021
      27      as1t  neutral  left      3  left 0.4173369
      29      as1t accuracy  left      3  left 2.4369634
      31      as1t    speed right      3 right 1.6178208
      33      as1t  neutral right      3 right 0.4569681
      35      as1t accuracy right      3  left 0.4380294
      37      as1t    speed  left      4 right 0.5229175
      39      as1t  neutral  left      4  left 1.3761103
      41      as1t accuracy  left      4  left 1.8648743
      43      as1t    speed right      4 right 0.5330341
      45      as1t  neutral right      4  left 0.8107023
      47      as1t accuracy right      4  left 0.5772930
      49      as1t    speed  left      5 right 0.5588077
      51      as1t  neutral  left      5 right 0.8228730
      53      as1t accuracy  left      5  left 1.4508758
      55      as1t    speed right      5  left 1.4554709
      57      as1t  neutral right      5 right 0.4077888
      59      as1t accuracy right      5  left 3.5878112
      61      as1t    speed  left      6 right 0.3267569
      63      as1t  neutral  left      6  left 1.4557686
      65      as1t accuracy  left      6 right 0.5694368
      67      as1t    speed right      6 right 2.3809979
      69      as1t  neutral right      6  left 1.5258845
      71      as1t accuracy right      6  left 0.7676759
      73      as1t    speed  left      7 right 0.5964250
      75      as1t  neutral  left      7  left 0.3998947
      77      as1t accuracy  left      7  left 0.8572288
      79      as1t    speed right      7 right 0.3848368
      81      as1t  neutral right      7 right 8.9306756
      83      as1t accuracy right      7 right 0.7719064
      85      as1t    speed  left      8 right 0.5408207
      87      as1t  neutral  left      8  left 0.8330954
      89      as1t accuracy  left      8  left 1.1042810
      91      as1t    speed right      8 right 0.9368667
      93      as1t  neutral right      8  left 1.2216048
      95      as1t accuracy right      8 right 1.2020308
      97      as1t    speed  left      9  left 0.4259109
      99      as1t  neutral  left      9  left 0.5140515
      101     as1t accuracy  left      9  left 0.6659963
      103     as1t    speed right      9 right 0.6366333
      105     as1t  neutral right      9 right 0.5687022
      107     as1t accuracy right      9 right 0.8979096
      109     as1t    speed  left     10 right 1.0555149
      111     as1t  neutral  left     10  left 0.7691387
      113     as1t accuracy  left     10  left 1.7141380
      115     as1t    speed right     10 right 0.3060367
      117     as1t  neutral right     10  left 0.3654133
      119     as1t accuracy right     10 right 1.5808686
      2       bd6t    speed  left      1  left 0.4824442
      4       bd6t  neutral  left      1  left 0.4312960
      6       bd6t accuracy  left      1 right 0.6836209
      8       bd6t    speed right      1 right 0.4036104
      10      bd6t  neutral right      1  left 0.8242776
      12      bd6t accuracy right      1 right 0.3735829
      14      bd6t    speed  left      2  left 1.4424930
      16      bd6t  neutral  left      2 right 1.5389097
      18      bd6t accuracy  left      2  left 0.4368082
      20      bd6t    speed right      2 right 0.4137409
      22      bd6t  neutral right      2 right 0.6206644
      24      bd6t accuracy right      2 right 1.9228466
      26      bd6t    speed  left      3 right 0.4377277
      28      bd6t  neutral  left      3 right 0.9936258
      30      bd6t accuracy  left      3  left 1.2167281
      32      bd6t    speed right      3 right 0.4239670
      34      bd6t  neutral right      3 right 1.2828606
      36      bd6t accuracy right      3 right 0.6546611
      38      bd6t    speed  left      4  left 0.7046626
      40      bd6t  neutral  left      4  left 0.6382084
      42      bd6t accuracy  left      4  left 0.5630979
      44      bd6t    speed right      4 right 0.3930529
      46      bd6t  neutral right      4 right 0.5637251
      48      bd6t accuracy right      4  left 0.7199701
      50      bd6t    speed  left      5  left 1.5805567
      52      bd6t  neutral  left      5 right 0.4444290
      54      bd6t accuracy  left      5  left 0.3881707
      56      bd6t    speed right      5 right 0.3867762
      58      bd6t  neutral right      5  left 0.3660658
      60      bd6t accuracy right      5  left 1.8875954
      62      bd6t    speed  left      6  left 0.7828390
      64      bd6t  neutral  left      6 right 0.5159542
      66      bd6t accuracy  left      6  left 0.3364400
      68      bd6t    speed right      6 right 0.4342419
      70      bd6t  neutral right      6 right 1.3316921
      72      bd6t accuracy right      6 right 0.4684726
      74      bd6t    speed  left      7 right 2.7144179
      76      bd6t  neutral  left      7  left 0.6602649
      78      bd6t accuracy  left      7  left 1.5521877
      80      bd6t    speed right      7 right 0.4362133
      82      bd6t  neutral right      7 right 1.5224105
      84      bd6t accuracy right      7 right 0.7384710
      86      bd6t    speed  left      8 right 0.3355405
      88      bd6t  neutral  left      8  left 3.3384981
      90      bd6t accuracy  left      8  left 1.9218797
      92      bd6t    speed right      8 right 1.4853188
      94      bd6t  neutral right      8  left 1.2472914
      96      bd6t accuracy right      8 right 0.5519437
      98      bd6t    speed  left      9  left 0.3093947
      100     bd6t  neutral  left      9 right 0.4081733
      102     bd6t accuracy  left      9 right 0.8035991
      104     bd6t    speed right      9 right 0.5259909
      106     bd6t  neutral right      9 right 0.4850557
      108     bd6t accuracy right      9 right 0.7247657
      110     bd6t    speed  left     10  left 0.4781965
      112     bd6t  neutral  left     10  left 0.4753573
      114     bd6t accuracy  left     10 right 0.8034570
      116     bd6t    speed right     10 right 0.4628794
      118     bd6t  neutral right     10 right 1.0207296
      120     bd6t accuracy right     10 right 0.3664488

