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
      make_data(p_LNR, design_LNR, n_trials = 10)
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
      make_data(p_LBA, design_LBA, n_trials = 10)
    Output
          subjects        E     S trials     R        rt
      1       as1t    speed  left      1  left 1.4633807
      5       as1t  neutral  left      1  left 2.1702340
      9       as1t accuracy  left      1  left 1.4852344
      13      as1t    speed right      1 right 0.9423191
      17      as1t  neutral right      1 right 1.6249727
      21      as1t accuracy right      1  left 1.6289250
      25      as1t    speed  left      2  left 0.9759322
      29      as1t  neutral  left      2  left 2.2181263
      33      as1t accuracy  left      2  left 1.5406683
      37      as1t    speed right      2 right 1.3311058
      41      as1t  neutral right      2 right 1.1949201
      45      as1t accuracy right      2 right 1.3299111
      49      as1t    speed  left      3  left 1.2642489
      53      as1t  neutral  left      3  left 2.0398431
      57      as1t accuracy  left      3  left 1.9756260
      61      as1t    speed right      3 right 2.1560090
      65      as1t  neutral right      3 right 2.0795816
      69      as1t accuracy right      3  left 1.4747640
      73      as1t    speed  left      4 right 1.0183507
      77      as1t  neutral  left      4  left 2.0591050
      81      as1t accuracy  left      4  left 2.0299050
      85      as1t    speed right      4 right 1.7444860
      89      as1t  neutral right      4 right 1.6223118
      93      as1t accuracy right      4 right 1.2664438
      97      as1t    speed  left      5 right 0.8479725
      101     as1t  neutral  left      5  left 2.5810174
      105     as1t accuracy  left      5  left 2.0660525
      109     as1t    speed right      5 right 1.1061526
      113     as1t  neutral right      5  left 1.5656606
      117     as1t accuracy right      5 right 1.1412750
      121     as1t    speed  left      6  left 1.2552418
      125     as1t  neutral  left      6  left 1.9701327
      129     as1t accuracy  left      6  left 1.2969671
      133     as1t    speed right      6 right 0.9069431
      137     as1t  neutral right      6 right 1.9277024
      141     as1t accuracy right      6 right 1.4062788
      145     as1t    speed  left      7 right 1.7466986
      149     as1t  neutral  left      7  left 1.9530480
      153     as1t accuracy  left      7  left 1.3515376
      157     as1t    speed right      7 right 1.0903481
      161     as1t  neutral right      7 right 2.0526427
      165     as1t accuracy right      7 right 1.1242640
      169     as1t    speed  left      8  left 1.7874152
      173     as1t  neutral  left      8 right 1.3507243
      177     as1t accuracy  left      8  left 1.7284682
      181     as1t    speed right      8 right 1.5629776
      185     as1t  neutral right      8 right 1.9109834
      189     as1t accuracy right      8 right 1.8096389
      193     as1t    speed  left      9  left 1.1294915
      197     as1t  neutral  left      9 right 1.4011195
      201     as1t accuracy  left      9 right 1.7618616
      205     as1t    speed right      9 right 0.9635151
      209     as1t  neutral right      9  left 0.9964543
      213     as1t accuracy right      9 right 2.2168005
      217     as1t    speed  left     10  left 1.2540213
      221     as1t  neutral  left     10  left 1.7238653
      225     as1t accuracy  left     10  left 1.3217336
      229     as1t    speed right     10 right 1.1129529
      233     as1t  neutral right     10 right 1.3371482
      237     as1t accuracy right     10 right 1.4788349
      3       bd6t    speed  left      1  left 1.1431110
      7       bd6t  neutral  left      1 right 1.7701594
      11      bd6t accuracy  left      1  left 1.4828631
      15      bd6t    speed right      1  left 0.8311659
      19      bd6t  neutral right      1 right 1.7368588
      23      bd6t accuracy right      1  left 1.1528011
      27      bd6t    speed  left      2  left 1.1984658
      31      bd6t  neutral  left      2 right 1.7074347
      35      bd6t accuracy  left      2  left 1.7261127
      39      bd6t    speed right      2  left 0.9345510
      43      bd6t  neutral right      2 right 1.7837362
      47      bd6t accuracy right      2 right 1.3253715
      51      bd6t    speed  left      3  left 0.9382072
      55      bd6t  neutral  left      3  left 1.6609535
      59      bd6t accuracy  left      3  left 1.2350623
      63      bd6t    speed right      3 right 1.2555053
      67      bd6t  neutral right      3 right 2.0236594
      71      bd6t accuracy right      3 right 1.8515205
      75      bd6t    speed  left      4  left 1.1070257
      79      bd6t  neutral  left      4  left 1.7790590
      83      bd6t accuracy  left      4  left 1.5055781
      87      bd6t    speed right      4 right 1.2143236
      91      bd6t  neutral right      4  left 1.5663614
      95      bd6t accuracy right      4 right 1.4461514
      99      bd6t    speed  left      5 right 1.2060786
      103     bd6t  neutral  left      5  left 1.4588243
      107     bd6t accuracy  left      5  left 1.3866129
      111     bd6t    speed right      5 right 1.6574897
      115     bd6t  neutral right      5  left 1.2277583
      119     bd6t accuracy right      5  left 2.1361914
      123     bd6t    speed  left      6  left 1.0186402
      127     bd6t  neutral  left      6 right 2.2555354
      131     bd6t accuracy  left      6 right 1.3513423
      135     bd6t    speed right      6 right 1.1363703
      139     bd6t  neutral right      6 right 2.2786338
      143     bd6t accuracy right      6 right 1.7136703
      147     bd6t    speed  left      7 right 1.1315176
      151     bd6t  neutral  left      7 right 2.8129833
      155     bd6t accuracy  left      7  left 1.8581362
      159     bd6t    speed right      7 right 1.4264524
      163     bd6t  neutral right      7  left 1.9479597
      167     bd6t accuracy right      7 right 1.3138114
      171     bd6t    speed  left      8 right 1.2901186
      175     bd6t  neutral  left      8  left 2.1920405
      179     bd6t accuracy  left      8  left 1.4267523
      183     bd6t    speed right      8  left 1.0080981
      187     bd6t  neutral right      8  left 1.9679805
      191     bd6t accuracy right      8 right 1.5223778
      195     bd6t    speed  left      9 right 1.3196324
      199     bd6t  neutral  left      9  left 1.4713529
      203     bd6t accuracy  left      9  left 2.2323876
      207     bd6t    speed right      9 right 1.2076100
      211     bd6t  neutral right      9 right 1.6209301
      215     bd6t accuracy right      9 right 1.6692198
      219     bd6t    speed  left     10 right 0.8206576
      223     bd6t  neutral  left     10  left 1.3604141
      227     bd6t accuracy  left     10  left 1.5717983
      231     bd6t    speed right     10  left 0.9161352
      235     bd6t  neutral right     10 right 1.3449899
      239     bd6t accuracy right     10 right 1.3311090

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
      make_data(p_RDM, design_RDM, n_trials = 10)
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
      as1t -3444.491
      bd6t -2932.240
      
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
      1       as1t    speed  left      1  left 0.6802333
      3       as1t  neutral  left      1  left 1.0005519
      5       as1t accuracy  left      1 right 2.2555073
      7       as1t    speed right      1 right 0.4950596
      9       as1t  neutral right      1 right 0.3838494
      11      as1t accuracy right      1 right 0.8643106
      13      as1t    speed  left      2 right 0.3581243
      15      as1t  neutral  left      2  left 0.7989559
      17      as1t accuracy  left      2  left 0.5415385
      19      as1t    speed right      2 right 0.4025097
      21      as1t  neutral right      2 right 0.7988730
      23      as1t accuracy right      2 right 1.1206340
      25      as1t    speed  left      3 right 0.8831889
      27      as1t  neutral  left      3  left 0.3569167
      29      as1t accuracy  left      3  left 0.5587834
      31      as1t    speed right      3 right 0.6047484
      33      as1t  neutral right      3 right 2.1543090
      35      as1t accuracy right      3 right 0.5008352
      37      as1t    speed  left      4  left 0.4499327
      39      as1t  neutral  left      4  left 1.3312508
      41      as1t accuracy  left      4  left 0.3797134
      43      as1t    speed right      4 right 0.3844647
      45      as1t  neutral right      4 right 0.9054355
      47      as1t accuracy right      4  left 1.3175715
      49      as1t    speed  left      5 right 0.3303770
      51      as1t  neutral  left      5 right 0.5558977
      53      as1t accuracy  left      5  left 0.5683490
      55      as1t    speed right      5 right 0.3144739
      57      as1t  neutral right      5 right 0.3347010
      59      as1t accuracy right      5 right 0.4547618
      61      as1t    speed  left      6  left 0.3646702
      63      as1t  neutral  left      6  left 0.8586173
      65      as1t accuracy  left      6  left 0.5523000
      67      as1t    speed right      6 right 0.3884522
      69      as1t  neutral right      6  left 0.5610446
      71      as1t accuracy right      6 right 0.5708478
      73      as1t    speed  left      7  left 0.6389701
      75      as1t  neutral  left      7  left 0.7979393
      77      as1t accuracy  left      7  left 0.6868122
      79      as1t    speed right      7  left 0.4383368
      81      as1t  neutral right      7  left 1.1010207
      83      as1t accuracy right      7  left 2.2889515
      85      as1t    speed  left      8  left 0.5084518
      87      as1t  neutral  left      8  left 0.8292358
      89      as1t accuracy  left      8  left 0.4086253
      91      as1t    speed right      8  left 0.3667670
      93      as1t  neutral right      8 right 0.7317996
      95      as1t accuracy right      8  left 4.9271514
      97      as1t    speed  left      9  left 0.4733514
      99      as1t  neutral  left      9 right 2.3900265
      101     as1t accuracy  left      9 right 0.6865439
      103     as1t    speed right      9  left 0.6913264
      105     as1t  neutral right      9 right 0.2915473
      107     as1t accuracy right      9 right 0.4503388
      109     as1t    speed  left     10  left 0.2684499
      111     as1t  neutral  left     10 right 0.9177357
      113     as1t accuracy  left     10  left 0.9141145
      115     as1t    speed right     10 right 0.3615293
      117     as1t  neutral right     10  left 1.6126562
      119     as1t accuracy right     10  left 0.4261822
      2       bd6t    speed  left      1 right 0.3656491
      4       bd6t  neutral  left      1  left 0.4481234
      6       bd6t accuracy  left      1  left 2.6101841
      8       bd6t    speed right      1  left 0.9775268
      10      bd6t  neutral right      1 right 0.4061499
      12      bd6t accuracy right      1  left 0.4199573
      14      bd6t    speed  left      2  left 0.3693564
      16      bd6t  neutral  left      2  left 0.4612578
      18      bd6t accuracy  left      2  left 0.3420450
      20      bd6t    speed right      2  left 1.2008347
      22      bd6t  neutral right      2  left 0.9842469
      24      bd6t accuracy right      2 right 0.7297536
      26      bd6t    speed  left      3  left 0.7347282
      28      bd6t  neutral  left      3  left 0.3711668
      30      bd6t accuracy  left      3 right 0.4676017
      32      bd6t    speed right      3 right 0.6443278
      34      bd6t  neutral right      3  left 0.5216554
      36      bd6t accuracy right      3 right 2.4432551
      38      bd6t    speed  left      4  left 0.4125840
      40      bd6t  neutral  left      4  left 0.3468602
      42      bd6t accuracy  left      4  left 0.9528720
      44      bd6t    speed right      4  left 0.5470503
      46      bd6t  neutral right      4 right 0.5961414
      48      bd6t accuracy right      4 right 0.5885120
      50      bd6t    speed  left      5 right 1.6893529
      52      bd6t  neutral  left      5  left 0.4786007
      54      bd6t accuracy  left      5  left 0.4911437
      56      bd6t    speed right      5 right 0.9986496
      58      bd6t  neutral right      5 right 0.5745838
      60      bd6t accuracy right      5 right 0.3989805
      62      bd6t    speed  left      6  left 1.0776926
      64      bd6t  neutral  left      6 right 0.6757066
      66      bd6t accuracy  left      6  left 0.3551867
      68      bd6t    speed right      6 right 0.3710058
      70      bd6t  neutral right      6 right 1.0897391
      72      bd6t accuracy right      6 right 1.6217274
      74      bd6t    speed  left      7  left 0.3417613
      76      bd6t  neutral  left      7 right 1.0540506
      78      bd6t accuracy  left      7  left 0.3110610
      80      bd6t    speed right      7  left 1.1326745
      82      bd6t  neutral right      7 right 7.5569366
      84      bd6t accuracy right      7  left 0.7934986
      86      bd6t    speed  left      8  left 1.4717633
      88      bd6t  neutral  left      8  left 0.5502727
      90      bd6t accuracy  left      8  left 2.8111990
      92      bd6t    speed right      8 right 0.3531640
      94      bd6t  neutral right      8 right 0.4038146
      96      bd6t accuracy right      8 right 0.5446294
      98      bd6t    speed  left      9  left 0.3481593
      100     bd6t  neutral  left      9 right 0.6100624
      102     bd6t accuracy  left      9  left 1.0539004
      104     bd6t    speed right      9 right 0.3785078
      106     bd6t  neutral right      9 right 0.8167383
      108     bd6t accuracy right      9 right 0.6870093
      110     bd6t    speed  left     10  left 0.6959749
      112     bd6t  neutral  left     10  left 0.5334544
      114     bd6t accuracy  left     10  left 1.1643393
      116     bd6t    speed right     10  left 0.8941779
      118     bd6t  neutral right     10  left 1.5046528
      120     bd6t accuracy right     10  left 0.7933934

