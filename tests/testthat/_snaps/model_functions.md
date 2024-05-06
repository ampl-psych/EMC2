# LNR

    Code
      init_chains(LNR_s, particles = 10, cores_per_chain = 4)[[1]]$samples
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
      
                        as1t         bd6t
      m           -1.7437564 -0.451325041
      m_lMd        0.8894434 -0.004573125
      m_Eneutral   1.3453151  1.767460316
      m_Eaccuracy -1.6999044 -2.108929606
      s            0.8516853  0.442019143
      t0          -1.5816652 -1.285432743
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -1325.892
      bd6t -1094.095
      
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
      1       as1t    speed  left      1 right 0.2755576
      5       as1t  neutral  left      1 right 0.2025982
      9       as1t accuracy  left      1 right 0.6768155
      13      as1t    speed right      1 right 0.7825357
      17      as1t  neutral right      1  left 0.4874694
      21      as1t accuracy right      1  left 0.3528860
      25      as1t    speed  left      2 right 0.3509746
      29      as1t  neutral  left      2 right 0.2619079
      33      as1t accuracy  left      2  left 0.2993013
      37      as1t    speed right      2 right 0.5555909
      41      as1t  neutral right      2 right 0.3509435
      45      as1t accuracy right      2 right 0.5820264
      49      as1t    speed  left      3  left 0.8384373
      53      as1t  neutral  left      3 right 0.2881257
      57      as1t accuracy  left      3  left 0.3733459
      61      as1t    speed right      3  left 0.2245143
      65      as1t  neutral right      3  left 0.2623760
      69      as1t accuracy right      3  left 0.2196318
      73      as1t    speed  left      4  left 0.3306013
      77      as1t  neutral  left      4 right 0.2787778
      81      as1t accuracy  left      4  left 0.5250431
      85      as1t    speed right      4  left 0.6169692
      89      as1t  neutral right      4  left 0.2521739
      93      as1t accuracy right      4  left 0.2112857
      97      as1t    speed  left      5 right 0.4373107
      101     as1t  neutral  left      5 right 0.8202014
      105     as1t accuracy  left      5 right 1.3263698
      109     as1t    speed right      5  left 0.3392336
      113     as1t  neutral right      5  left 0.3516762
      117     as1t accuracy right      5  left 0.2223888
      121     as1t    speed  left      6 right 0.2460490
      125     as1t  neutral  left      6 right 0.5752245
      129     as1t accuracy  left      6 right 0.2814323
      133     as1t    speed right      6 right 0.8500997
      137     as1t  neutral right      6  left 0.3181316
      141     as1t accuracy right      6 right 0.3544758
      145     as1t    speed  left      7 right 0.2570443
      149     as1t  neutral  left      7  left 2.1088579
      153     as1t accuracy  left      7 right 0.6006545
      157     as1t    speed right      7 right 0.4709772
      161     as1t  neutral right      7 right 0.3533543
      165     as1t accuracy right      7 right 0.4350710
      169     as1t    speed  left      8  left 0.3966046
      173     as1t  neutral  left      8 right 0.2598129
      177     as1t accuracy  left      8 right 0.2062906
      181     as1t    speed right      8  left 0.2118133
      185     as1t  neutral right      8  left 0.3757369
      189     as1t accuracy right      8  left 1.7571616
      193     as1t    speed  left      9 right 0.2211219
      197     as1t  neutral  left      9 right 0.2430565
      201     as1t accuracy  left      9 right 0.2394999
      205     as1t    speed right      9  left 0.2794567
      209     as1t  neutral right      9  left 0.3372966
      213     as1t accuracy right      9  left 0.2787819
      217     as1t    speed  left     10 right 0.2552589
      221     as1t  neutral  left     10 right 0.2233440
      225     as1t accuracy  left     10 right 0.3310308
      229     as1t    speed right     10  left 0.7254711
      233     as1t  neutral right     10  left 0.5264426
      237     as1t accuracy right     10  left 0.2192021
      3       bd6t    speed  left      1 right 0.2790124
      7       bd6t  neutral  left      1 right 1.0933827
      11      bd6t accuracy  left      1  left 0.2470619
      15      bd6t    speed right      1  left 0.2148357
      19      bd6t  neutral right      1 right 0.3537985
      23      bd6t accuracy right      1  left 0.7661086
      27      bd6t    speed  left      2 right 0.4673085
      31      bd6t  neutral  left      2 right 0.2393454
      35      bd6t accuracy  left      2 right 0.3449137
      39      bd6t    speed right      2 right 0.4570706
      43      bd6t  neutral right      2  left 0.4675106
      47      bd6t accuracy right      2 right 0.2536847
      51      bd6t    speed  left      3 right 0.6119290
      55      bd6t  neutral  left      3  left 0.4360208
      59      bd6t accuracy  left      3 right 0.2479397
      63      bd6t    speed right      3 right 0.2108371
      67      bd6t  neutral right      3  left 0.5689329
      71      bd6t accuracy right      3  left 0.2055206
      75      bd6t    speed  left      4  left 0.8239745
      79      bd6t  neutral  left      4 right 0.5916102
      83      bd6t accuracy  left      4 right 0.3083552
      87      bd6t    speed right      4  left 0.2074455
      91      bd6t  neutral right      4  left 1.0291990
      95      bd6t accuracy right      4  left 0.2387766
      99      bd6t    speed  left      5 right 0.2171693
      103     bd6t  neutral  left      5 right 0.2477152
      107     bd6t accuracy  left      5 right 0.2589259
      111     bd6t    speed right      5  left 0.5285502
      115     bd6t  neutral right      5  left 0.2650557
      119     bd6t accuracy right      5  left 0.2311833
      123     bd6t    speed  left      6 right 0.3348222
      127     bd6t  neutral  left      6 right 0.3242666
      131     bd6t accuracy  left      6 right 2.9628795
      135     bd6t    speed right      6  left 0.2136685
      139     bd6t  neutral right      6  left 0.2198051
      143     bd6t accuracy right      6  left 0.2291927
      147     bd6t    speed  left      7 right 0.2735157
      151     bd6t  neutral  left      7 right 0.2535005
      155     bd6t accuracy  left      7 right 0.4116405
      159     bd6t    speed right      7  left 0.2040516
      163     bd6t  neutral right      7  left 0.2337994
      167     bd6t accuracy right      7  left 0.2114074
      171     bd6t    speed  left      8  left 1.1127029
      175     bd6t  neutral  left      8 right 0.2745290
      179     bd6t accuracy  left      8 right 0.6393959
      183     bd6t    speed right      8  left 0.2271189
      187     bd6t  neutral right      8  left 0.2758488
      191     bd6t accuracy right      8  left 0.2080041
      195     bd6t    speed  left      9 right 0.3147874
      199     bd6t  neutral  left      9 right 0.8385318
      203     bd6t accuracy  left      9  left 0.2522106
      207     bd6t    speed right      9 right 0.3602160
      211     bd6t  neutral right      9  left 0.2108524
      215     bd6t accuracy right      9  left 0.2177526
      219     bd6t    speed  left     10  left 0.3079804
      223     bd6t  neutral  left     10 right 0.3360226
      227     bd6t accuracy  left     10 right 0.4199177
      231     bd6t    speed right     10  left 2.0354719
      235     bd6t  neutral right     10  left 0.2788527
      239     bd6t accuracy right     10 right 0.7382237

# LBA

    Code
      init_chains(LBA_s, particles = 10, cores_per_chain = 4)[[1]]$samples
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
      v           -1.2748892 -0.8494871
      v_lMd        0.5956117  0.4824282
      sv_lMTRUE    1.3952474  1.2472335
      B           -1.9260735 -2.0320307
      B_Eneutral   0.2458292  0.3543701
      B_Eaccuracy -1.3690513 -0.8626310
      B_lRright    0.3166578  0.3197739
      A           -2.1157142 -1.9422547
      t0          -0.7838674 -0.8355143
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -9920.128
      bd6t -9421.163
      
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
      1       as1t    speed  left      1  left  1.990572
      5       as1t  neutral  left      1  left  1.958785
      9       as1t accuracy  left      1  left  2.253859
      13      as1t    speed right      1 right  1.630186
      17      as1t  neutral right      1 right  4.074374
      21      as1t accuracy right      1 right  1.788465
      25      as1t    speed  left      2  left  1.749243
      29      as1t  neutral  left      2  left  2.954065
      33      as1t accuracy  left      2 right  1.976452
      37      as1t    speed right      2 right  2.565688
      41      as1t  neutral right      2  left  5.816528
      45      as1t accuracy right      2  left  1.429333
      49      as1t    speed  left      3  left  1.710729
      53      as1t  neutral  left      3 right  3.450728
      57      as1t accuracy  left      3 right  2.736339
      61      as1t    speed right      3 right  2.432148
      65      as1t  neutral right      3  left  3.352541
      69      as1t accuracy right      3 right  7.496598
      73      as1t    speed  left      4  left  1.942640
      77      as1t  neutral  left      4  left  5.364524
      81      as1t accuracy  left      4 right  1.824391
      85      as1t    speed right      4  left  1.909663
      89      as1t  neutral right      4 right  2.890292
      93      as1t accuracy right      4 right  2.406236
      97      as1t    speed  left      5 right  1.454888
      101     as1t  neutral  left      5 right  2.176252
      105     as1t accuracy  left      5 right  2.045797
      109     as1t    speed right      5  left  1.588279
      113     as1t  neutral right      5  left  3.011128
      117     as1t accuracy right      5 right  3.226346
      121     as1t    speed  left      6  left  1.951459
      125     as1t  neutral  left      6  left  3.816063
      129     as1t accuracy  left      6  left  1.971351
      133     as1t    speed right      6 right  1.443986
      137     as1t  neutral right      6  left  6.356136
      141     as1t accuracy right      6 right  1.794648
      145     as1t    speed  left      7  left  5.945754
      149     as1t  neutral  left      7 right  3.256723
      153     as1t accuracy  left      7 right  1.485119
      157     as1t    speed right      7 right  1.755943
      161     as1t  neutral right      7 right  3.638573
      165     as1t accuracy right      7 right  3.966576
      169     as1t    speed  left      8  left  2.434721
      173     as1t  neutral  left      8 right  4.435659
      177     as1t accuracy  left      8  left  2.058331
      181     as1t    speed right      8 right  2.244097
      185     as1t  neutral right      8 right  2.294874
      189     as1t accuracy right      8 right  4.038249
      193     as1t    speed  left      9 right  2.195047
      197     as1t  neutral  left      9  left  3.435898
      201     as1t accuracy  left      9 right  4.878960
      205     as1t    speed right      9  left  1.858387
      209     as1t  neutral right      9  left  2.860723
      213     as1t accuracy right      9 right  2.045361
      217     as1t    speed  left     10  left  1.290320
      221     as1t  neutral  left     10 right  2.845883
      225     as1t accuracy  left     10  left  3.524384
      229     as1t    speed right     10 right  2.003214
      233     as1t  neutral right     10  left  4.397111
      237     as1t accuracy right     10 right  8.332840
      3       bd6t    speed  left      1  left  2.950885
      7       bd6t  neutral  left      1  left  3.513701
      11      bd6t accuracy  left      1 right  3.504860
      15      bd6t    speed right      1 right  1.703662
      19      bd6t  neutral right      1 right  2.338615
      23      bd6t accuracy right      1 right  3.286494
      27      bd6t    speed  left      2  left  2.675333
      31      bd6t  neutral  left      2 right  2.236844
      35      bd6t accuracy  left      2 right  2.618327
      39      bd6t    speed right      2 right  1.115592
      43      bd6t  neutral right      2 right  2.914574
      47      bd6t accuracy right      2 right  7.050540
      51      bd6t    speed  left      3 right  5.093368
      55      bd6t  neutral  left      3 right  2.491515
      59      bd6t accuracy  left      3  left  8.569944
      63      bd6t    speed right      3 right  1.531315
      67      bd6t  neutral right      3  left  3.829156
      71      bd6t accuracy right      3 right  2.737317
      75      bd6t    speed  left      4 right  3.316738
      79      bd6t  neutral  left      4  left  3.953523
      83      bd6t accuracy  left      4  left  4.318171
      87      bd6t    speed right      4 right  2.015606
      91      bd6t  neutral right      4 right  2.201556
      95      bd6t accuracy right      4 right  1.980273
      99      bd6t    speed  left      5  left  2.839655
      103     bd6t  neutral  left      5  left  2.629991
      107     bd6t accuracy  left      5  left  2.346083
      111     bd6t    speed right      5 right  2.981165
      115     bd6t  neutral right      5 right 12.884714
      119     bd6t accuracy right      5 right  2.252274
      123     bd6t    speed  left      6  left  2.219416
      127     bd6t  neutral  left      6  left  2.856504
      131     bd6t accuracy  left      6  left  1.792302
      135     bd6t    speed right      6  left  2.644760
      139     bd6t  neutral right      6 right  3.735933
      143     bd6t accuracy right      6 right  2.274908
      147     bd6t    speed  left      7  left  2.218290
      151     bd6t  neutral  left      7  left  3.691373
      155     bd6t accuracy  left      7  left  3.100901
      159     bd6t    speed right      7 right  2.076683
      163     bd6t  neutral right      7 right  2.945258
      167     bd6t accuracy right      7  left  1.726307
      171     bd6t    speed  left      8  left  1.935618
      175     bd6t  neutral  left      8  left  4.079779
      179     bd6t accuracy  left      8 right  2.731046
      183     bd6t    speed right      8 right  2.807654
      187     bd6t  neutral right      8  left  4.231869
      191     bd6t accuracy right      8 right  1.175273
      195     bd6t    speed  left      9 right  2.538113
      199     bd6t  neutral  left      9  left  3.198599
      203     bd6t accuracy  left      9 right  3.837472
      207     bd6t    speed right      9 right  5.222315
      211     bd6t  neutral right      9 right  3.163063
      215     bd6t accuracy right      9  left  1.241953
      219     bd6t    speed  left     10  left  2.314831
      223     bd6t  neutral  left     10  left  5.634757
      227     bd6t accuracy  left     10 right  1.402478
      231     bd6t    speed right     10 right  1.439969
      235     bd6t  neutral right     10  left  4.437530
      239     bd6t accuracy right     10 right  3.100985

# RDM

    Code
      init_chains(RDM_s, particles = 10, cores_per_chain = 4)[[1]]$samples
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
      v           -1.2748892 -0.7532267
      v_lMd        0.5956117  0.2242688
      s_lMTRUE     1.3952474  1.3894611
      B           -1.9260735 -1.4897176
      B_Eneutral   0.2458292  0.3957119
      B_Eaccuracy -1.3690513 -1.3320533
      B_lRright    0.3166578  0.3242525
      A           -2.1157142 -2.1016525
      t0          -0.7838674 -0.8287774
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -8613.200
      bd6t -9004.968
      
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
      1       as1t    speed  left      1 right 0.8222103
      5       as1t  neutral  left      1  left 1.1314587
      9       as1t accuracy  left      1  left 1.2634804
      13      as1t    speed right      1 right 0.8840027
      17      as1t  neutral right      1 right 1.5106676
      21      as1t accuracy right      1 right 1.2981546
      25      as1t    speed  left      2  left 1.1538982
      29      as1t  neutral  left      2  left 1.3059854
      33      as1t accuracy  left      2  left 1.2299928
      37      as1t    speed right      2 right 1.1582296
      41      as1t  neutral right      2 right 1.6234553
      45      as1t accuracy right      2 right 1.0841324
      49      as1t    speed  left      3  left 1.0603885
      53      as1t  neutral  left      3  left 1.2417056
      57      as1t accuracy  left      3  left 1.2220239
      61      as1t    speed right      3 right 0.9223162
      65      as1t  neutral right      3 right 1.2235672
      69      as1t accuracy right      3 right 1.0746118
      73      as1t    speed  left      4  left 1.2480860
      77      as1t  neutral  left      4  left 1.7592539
      81      as1t accuracy  left      4  left 1.2760465
      85      as1t    speed right      4 right 1.1822597
      89      as1t  neutral right      4 right 1.5625712
      93      as1t accuracy right      4  left 1.1736538
      97      as1t    speed  left      5  left 1.1735565
      101     as1t  neutral  left      5  left 1.2356905
      105     as1t accuracy  left      5  left 1.2013194
      109     as1t    speed right      5  left 0.7937991
      113     as1t  neutral right      5 right 1.6722165
      117     as1t accuracy right      5 right 1.0033618
      121     as1t    speed  left      6  left 0.9482384
      125     as1t  neutral  left      6  left 1.1923728
      129     as1t accuracy  left      6  left 1.0289639
      133     as1t    speed right      6 right 0.8279774
      137     as1t  neutral right      6 right 1.0370566
      141     as1t accuracy right      6 right 1.5487737
      145     as1t    speed  left      7  left 0.7783782
      149     as1t  neutral  left      7  left 1.7738582
      153     as1t accuracy  left      7  left 1.6109707
      157     as1t    speed right      7 right 1.0553684
      161     as1t  neutral right      7 right 1.4806025
      165     as1t accuracy right      7 right 0.9980931
      169     as1t    speed  left      8  left 1.0960504
      173     as1t  neutral  left      8  left 1.7496989
      177     as1t accuracy  left      8  left 1.7875851
      181     as1t    speed right      8 right 1.0779126
      185     as1t  neutral right      8  left 1.5847181
      189     as1t accuracy right      8 right 1.1334881
      193     as1t    speed  left      9  left 1.0620281
      197     as1t  neutral  left      9  left 1.4466611
      201     as1t accuracy  left      9  left 1.4115570
      205     as1t    speed right      9 right 0.9704942
      209     as1t  neutral right      9 right 1.3560140
      213     as1t accuracy right      9 right 0.7370308
      217     as1t    speed  left     10  left 0.7722566
      221     as1t  neutral  left     10  left 1.3131958
      225     as1t accuracy  left     10  left 1.2898566
      229     as1t    speed right     10 right 1.0640279
      233     as1t  neutral right     10 right 1.8281453
      237     as1t accuracy right     10 right 1.3131207
      3       bd6t    speed  left      1 right 0.9581313
      7       bd6t  neutral  left      1  left 1.9106312
      11      bd6t accuracy  left      1  left 1.1192002
      15      bd6t    speed right      1 right 0.8238558
      19      bd6t  neutral right      1 right 1.2953471
      23      bd6t accuracy right      1 right 1.5123866
      27      bd6t    speed  left      2  left 0.9251321
      31      bd6t  neutral  left      2  left 1.5232809
      35      bd6t accuracy  left      2  left 1.0176998
      39      bd6t    speed right      2 right 0.8758530
      43      bd6t  neutral right      2 right 1.6483671
      47      bd6t accuracy right      2 right 1.3332667
      51      bd6t    speed  left      3  left 0.8251797
      55      bd6t  neutral  left      3  left 1.5676923
      59      bd6t accuracy  left      3  left 1.1502573
      63      bd6t    speed right      3  left 1.1687351
      67      bd6t  neutral right      3 right 1.1901727
      71      bd6t accuracy right      3 right 1.1179527
      75      bd6t    speed  left      4  left 0.7771911
      79      bd6t  neutral  left      4  left 1.1804826
      83      bd6t accuracy  left      4  left 1.0148929
      87      bd6t    speed right      4 right 1.4711346
      91      bd6t  neutral right      4 right 1.2032808
      95      bd6t accuracy right      4 right 1.2201967
      99      bd6t    speed  left      5  left 1.0523840
      103     bd6t  neutral  left      5  left 1.7768784
      107     bd6t accuracy  left      5  left 1.3378905
      111     bd6t    speed right      5 right 1.1390029
      115     bd6t  neutral right      5 right 1.6991595
      119     bd6t accuracy right      5 right 1.4529666
      123     bd6t    speed  left      6  left 0.9987922
      127     bd6t  neutral  left      6  left 1.4966811
      131     bd6t accuracy  left      6  left 1.3825747
      135     bd6t    speed right      6 right 1.1917950
      139     bd6t  neutral right      6  left 1.6583115
      143     bd6t accuracy right      6 right 1.5113941
      147     bd6t    speed  left      7  left 0.8290193
      151     bd6t  neutral  left      7  left 1.4392097
      155     bd6t accuracy  left      7  left 1.2893487
      159     bd6t    speed right      7 right 1.2399665
      163     bd6t  neutral right      7 right 1.4434496
      167     bd6t accuracy right      7  left 0.8043393
      171     bd6t    speed  left      8  left 0.8558665
      175     bd6t  neutral  left      8  left 1.6613360
      179     bd6t accuracy  left      8  left 1.4116084
      183     bd6t    speed right      8 right 0.8360994
      187     bd6t  neutral right      8 right 1.4510405
      191     bd6t accuracy right      8 right 1.0945947
      195     bd6t    speed  left      9  left 1.2783586
      199     bd6t  neutral  left      9  left 1.2208552
      203     bd6t accuracy  left      9  left 1.3174780
      207     bd6t    speed right      9 right 1.0527764
      211     bd6t  neutral right      9 right 1.3167410
      215     bd6t accuracy right      9 right 1.4548632
      219     bd6t    speed  left     10 right 0.8499064
      223     bd6t  neutral  left     10 right 2.0337011
      227     bd6t accuracy  left     10  left 0.9986664
      231     bd6t    speed right     10 right 1.1677391
      235     bd6t  neutral right     10 right 1.1760958
      239     bd6t accuracy right     10 right 1.5289592

# DDM

    Code
      init_chains(DDM_s, particles = 10, cores_per_chain = 4)[[1]]$samples
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
      
                         as1t       bd6t
      v_Sleft     -0.81572253 -1.0867281
      v_Sright     1.03870104  0.7321494
      a            1.22854033  1.3829047
      a_Eneutral  -2.14289844 -1.4698452
      a_Eaccuracy -0.09083814  0.3671716
      t0          -1.13678510 -1.6263327
      Z            0.10742647  0.2940408
      sv          -2.61001680 -2.2368029
      SZ          -0.56886020 -0.2263479
      
      
      $stage
      [1] "init"
      
      $subj_ll
                [,1]
      as1t -4673.532
      bd6t -3862.727
      
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
      1       as1t    speed  left      1  left 0.5507627
      3       as1t  neutral  left      1  left 0.2980996
      5       as1t accuracy  left      1  left 0.4347898
      7       as1t    speed right      1 right 0.4203045
      9       as1t  neutral right      1 right 0.6851127
      11      as1t accuracy right      1 right 0.8796972
      13      as1t    speed  left      2  left 0.4295901
      15      as1t  neutral  left      2  left 0.5689743
      17      as1t accuracy  left      2  left 1.5632208
      19      as1t    speed right      2 right 0.8167386
      21      as1t  neutral right      2 right 0.3568233
      23      as1t accuracy right      2  left 1.2735248
      25      as1t    speed  left      3  left 0.4033563
      27      as1t  neutral  left      3  left 0.7681822
      29      as1t accuracy  left      3  left 0.4684874
      31      as1t    speed right      3 right 0.5608058
      33      as1t  neutral right      3 right 0.6136671
      35      as1t accuracy right      3 right 0.5653802
      37      as1t    speed  left      4  left 1.6801981
      39      as1t  neutral  left      4  left 1.9731549
      41      as1t accuracy  left      4  left 0.7910344
      43      as1t    speed right      4 right 1.0703770
      45      as1t  neutral right      4 right 2.3420273
      47      as1t accuracy right      4 right 3.4174948
      49      as1t    speed  left      5  left 0.9345521
      51      as1t  neutral  left      5  left 1.8056773
      53      as1t accuracy  left      5  left 0.6651166
      55      as1t    speed right      5  left 0.5955325
      57      as1t  neutral right      5 right 0.5349894
      59      as1t accuracy right      5 right 0.4044412
      61      as1t    speed  left      6 right 1.6731333
      63      as1t  neutral  left      6  left 0.7524843
      65      as1t accuracy  left      6  left 2.5218907
      67      as1t    speed right      6 right 0.6827266
      69      as1t  neutral right      6 right 0.7440494
      71      as1t accuracy right      6 right 0.3940918
      73      as1t    speed  left      7  left 0.5183934
      75      as1t  neutral  left      7  left 0.5999088
      77      as1t accuracy  left      7 right 1.6495720
      79      as1t    speed right      7 right 0.3325124
      81      as1t  neutral right      7 right 0.7845693
      83      as1t accuracy right      7 right 0.2984396
      85      as1t    speed  left      8  left 0.6047512
      87      as1t  neutral  left      8  left 0.4041152
      89      as1t accuracy  left      8 right 1.3347864
      91      as1t    speed right      8  left 0.5827979
      93      as1t  neutral right      8  left 0.4515240
      95      as1t accuracy right      8 right 0.4886136
      97      as1t    speed  left      9  left 0.9369213
      99      as1t  neutral  left      9  left 0.4578845
      101     as1t accuracy  left      9 right 0.9257603
      103     as1t    speed right      9 right 0.7714174
      105     as1t  neutral right      9 right 1.1145605
      107     as1t accuracy right      9 right 0.6190086
      109     as1t    speed  left     10 right 1.0599100
      111     as1t  neutral  left     10 right 0.5306243
      113     as1t accuracy  left     10  left 0.8009746
      115     as1t    speed right     10 right 0.3642373
      117     as1t  neutral right     10  left 0.9159588
      119     as1t accuracy right     10 right 0.5240794
      2       bd6t    speed  left      1  left 0.3946918
      4       bd6t  neutral  left      1  left 0.4812603
      6       bd6t accuracy  left      1  left 0.6539427
      8       bd6t    speed right      1 right 0.3521524
      10      bd6t  neutral right      1 right 1.0890507
      12      bd6t accuracy right      1  left 0.4679715
      14      bd6t    speed  left      2 right 0.8749487
      16      bd6t  neutral  left      2 right 0.9518047
      18      bd6t accuracy  left      2 right 1.0540175
      20      bd6t    speed right      2 right 0.6929149
      22      bd6t  neutral right      2  left 0.6114332
      24      bd6t accuracy right      2 right 0.4172701
      26      bd6t    speed  left      3  left 0.6604484
      28      bd6t  neutral  left      3  left 0.4910742
      30      bd6t accuracy  left      3  left 0.3494988
      32      bd6t    speed right      3 right 0.3744440
      34      bd6t  neutral right      3  left 0.6628965
      36      bd6t accuracy right      3  left 1.6557658
      38      bd6t    speed  left      4  left 0.4031482
      40      bd6t  neutral  left      4 right 1.4050725
      42      bd6t accuracy  left      4  left 0.6846832
      44      bd6t    speed right      4 right 1.0583867
      46      bd6t  neutral right      4 right 0.8893607
      48      bd6t accuracy right      4  left 0.5603036
      50      bd6t    speed  left      5  left 0.5034308
      52      bd6t  neutral  left      5  left 0.4936842
      54      bd6t accuracy  left      5  left 0.3963543
      56      bd6t    speed right      5 right 0.3655050
      58      bd6t  neutral right      5  left 0.4395762
      60      bd6t accuracy right      5 right 0.6701831
      62      bd6t    speed  left      6  left 1.7902786
      64      bd6t  neutral  left      6 right 0.4278195
      66      bd6t accuracy  left      6 right 1.0987299
      68      bd6t    speed right      6 right 0.4911527
      70      bd6t  neutral right      6 right 0.3477567
      72      bd6t accuracy right      6 right 0.6382508
      74      bd6t    speed  left      7  left 0.6370293
      76      bd6t  neutral  left      7  left 1.2014549
      78      bd6t accuracy  left      7  left 0.5033713
      80      bd6t    speed right      7 right 0.7577291
      82      bd6t  neutral right      7  left 0.4967881
      84      bd6t accuracy right      7  left 0.7301008
      86      bd6t    speed  left      8  left 2.7821642
      88      bd6t  neutral  left      8 right 0.5305342
      90      bd6t accuracy  left      8  left 0.6244036
      92      bd6t    speed right      8 right 1.6683624
      94      bd6t  neutral right      8  left 2.3045041
      96      bd6t accuracy right      8  left 2.7339417
      98      bd6t    speed  left      9 right 0.3595699
      100     bd6t  neutral  left      9  left 0.3446541
      102     bd6t accuracy  left      9 right 0.3740896
      104     bd6t    speed right      9  left 0.6720155
      106     bd6t  neutral right      9  left 1.3486839
      108     bd6t accuracy right      9 right 0.5775676
      110     bd6t    speed  left     10  left 0.3225284
      112     bd6t  neutral  left     10  left 0.7773451
      114     bd6t accuracy  left     10  left 0.9295459
      116     bd6t    speed right     10 right 0.8894185
      118     bd6t  neutral right     10 right 0.9539637
      120     bd6t accuracy right     10 right 1.1460566

