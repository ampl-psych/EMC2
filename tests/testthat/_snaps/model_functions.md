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
          trials subjects        E     S     R        rt
      1        1     as1t    speed  left right 0.6621568
      3        2     as1t    speed  left right 0.2057630
      5        3     as1t    speed  left right 0.2148357
      7        4     as1t    speed  left  left 0.3104419
      9        5     as1t    speed  left  left 0.3391626
      11       6     as1t    speed  left right 0.3029197
      13       7     as1t    speed  left right 1.3993717
      15       8     as1t    speed  left right 0.2306459
      17       9     as1t    speed  left  left 1.1688779
      19      10     as1t    speed  left right 0.4420536
      41      11     as1t  neutral  left right 0.6309831
      43      12     as1t  neutral  left  left 0.2536847
      45      13     as1t  neutral  left right 0.2319419
      47      14     as1t  neutral  left right 0.2159625
      49      15     as1t  neutral  left right 0.3556382
      51      16     as1t  neutral  left  left 0.2119769
      53      17     as1t  neutral  left right 0.2055206
      55      18     as1t  neutral  left right 0.2933270
      57      19     as1t  neutral  left right 0.9621727
      59      20     as1t  neutral  left right 0.2841160
      81      21     as1t accuracy  left right 0.2082286
      83      22     as1t accuracy  left right 1.0291990
      85      23     as1t accuracy  left right 0.2387766
      87      24     as1t accuracy  left right 0.3454905
      89      25     as1t accuracy  left  left 0.3402073
      91      26     as1t accuracy  left right 0.2931404
      93      27     as1t accuracy  left right 0.5631042
      95      28     as1t accuracy  left right 0.2650557
      97      29     as1t accuracy  left right 0.2311833
      99      30     as1t accuracy  left right 0.3249946
      121     31     as1t    speed right  left 0.3348222
      123     32     as1t    speed right  left 0.3124411
      125     33     as1t    speed right right 0.3009970
      127     34     as1t    speed right right 0.3324145
      129     35     as1t    speed right  left 0.2808354
      131     36     as1t    speed right  left 0.5841607
      133     37     as1t    speed right  left 0.2735157
      135     38     as1t    speed right  left 0.2484093
      137     39     as1t    speed right right 0.2299376
      139     40     as1t    speed right  left 0.2606403
      161     41     as1t  neutral right right 0.2842899
      163     42     as1t  neutral right  left 0.2101295
      165     43     as1t  neutral right  left 1.4312545
      167     44     as1t  neutral right  left 0.2745290
      169     45     as1t  neutral right right 0.4214576
      171     46     as1t  neutral right  left 0.2168920
      173     47     as1t  neutral right right 0.2591430
      175     48     as1t  neutral right  left 1.0453442
      177     49     as1t  neutral right  left 0.3268598
      179     50     as1t  neutral right right 0.2522106
      201     51     as1t accuracy right right 2.6889300
      203     52     as1t accuracy right  left 0.2239633
      205     53     as1t accuracy right right 0.3311753
      207     54     as1t accuracy right right 0.3193368
      209     55     as1t accuracy right  left 0.4696293
      211     56     as1t accuracy right  left 0.3360226
      213     57     as1t accuracy right  left 0.4199177
      215     58     as1t accuracy right right 0.7826471
      217     59     as1t accuracy right  left 0.7467092
      219     60     as1t accuracy right right 0.2285965
      21       1     bd6t    speed  left right 0.3521819
      23       2     bd6t    speed  left  left 2.5742499
      25       3     bd6t    speed  left  left 0.2447084
      27       4     bd6t    speed  left right 0.3353470
      29       5     bd6t    speed  left right 0.3184436
      31       6     bd6t    speed  left right 0.4884841
      33       7     bd6t    speed  left right 0.4399945
      35       8     bd6t    speed  left right 0.2215287
      37       9     bd6t    speed  left right 0.2100897
      39      10     bd6t    speed  left right 0.2133570
      61      11     bd6t  neutral  left right 0.2446073
      63      12     bd6t  neutral  left  left 0.3726230
      65      13     bd6t  neutral  left right 0.4092510
      67      14     bd6t  neutral  left  left 0.5294351
      69      15     bd6t  neutral  left right 0.2143891
      71      16     bd6t  neutral  left right 0.9475786
      73      17     bd6t  neutral  left right 0.2133931
      75      18     bd6t  neutral  left right 1.2623103
      77      19     bd6t  neutral  left right 0.4351044
      79      20     bd6t  neutral  left right 0.4247098
      101     21     bd6t accuracy  left  left 0.3207819
      103     22     bd6t accuracy  left right 0.2389074
      105     23     bd6t accuracy  left right 0.2820200
      107     24     bd6t accuracy  left right 0.2198767
      109     25     bd6t accuracy  left right 0.4524300
      111     26     bd6t accuracy  left right 0.2300322
      113     27     bd6t accuracy  left right 0.2416458
      115     28     bd6t accuracy  left right 0.2191954
      117     29     bd6t accuracy  left right 0.6286309
      119     30     bd6t accuracy  left right 0.2628880
      141     31     bd6t    speed right  left 0.6651145
      143     32     bd6t    speed right  left 0.2299034
      145     33     bd6t    speed right right 2.4026591
      147     34     bd6t    speed right right 0.3498691
      149     35     bd6t    speed right right 0.3024430
      151     36     bd6t    speed right  left 0.3085880
      153     37     bd6t    speed right  left 0.2907286
      155     38     bd6t    speed right  left 0.2056963
      157     39     bd6t    speed right  left 0.2142630
      159     40     bd6t    speed right  left 0.4062331
      181     41     bd6t  neutral right  left 0.2521459
      183     42     bd6t  neutral right right 0.5737047
      185     43     bd6t  neutral right  left 0.2679879
      187     44     bd6t  neutral right  left 0.6547885
      189     45     bd6t  neutral right  left 0.7484957
      191     46     bd6t  neutral right  left 0.7427598
      193     47     bd6t  neutral right  left 0.2639333
      195     48     bd6t  neutral right  left 0.2970921
      197     49     bd6t  neutral right  left 0.4422466
      199     50     bd6t  neutral right  left 0.9803098
      221     51     bd6t accuracy right  left 0.7469548
      223     52     bd6t accuracy right right 0.2676699
      225     53     bd6t accuracy right right 1.3892407
      227     54     bd6t accuracy right  left 0.2842021
      229     55     bd6t accuracy right  left 0.3571977
      231     56     bd6t accuracy right  left 0.8165039
      233     57     bd6t accuracy right  left 0.2736779
      235     58     bd6t accuracy right  left 0.4250408
      237     59     bd6t accuracy right  left 0.3621805
      239     60     bd6t accuracy right  left 0.3439226

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
          trials subjects        E     S     R        rt
      1        1     as1t    speed  left  left 1.4633807
      3        2     as1t    speed  left  left 1.4760285
      5        3     as1t    speed  left  left 1.1592134
      7        4     as1t    speed  left right 1.0073454
      9        5     as1t    speed  left  left 1.1259135
      11       6     as1t    speed  left  left 1.0885685
      13       7     as1t    speed  left  left 0.9759322
      15       8     as1t    speed  left  left 1.5074960
      17       9     as1t    speed  left  left 1.1934289
      19      10     as1t    speed  left  left 1.0744474
      41      11     as1t  neutral  left right 1.0673134
      43      12     as1t  neutral  left right 1.7074722
      45      13     as1t  neutral  left  left 1.9051035
      47      14     as1t  neutral  left  left 2.0398431
      49      15     as1t  neutral  left  left 2.3613085
      51      16     as1t  neutral  left right 1.9783355
      53      17     as1t  neutral  left  left 1.9133404
      55      18     as1t  neutral  left  left 1.5454823
      57      19     as1t  neutral  left right 1.4862452
      59      20     as1t  neutral  left  left 2.0591050
      81      21     as1t accuracy  left  left 2.0299050
      83      22     as1t accuracy  left  left 1.6495790
      85      23     as1t accuracy  left right 1.6319596
      87      24     as1t accuracy  left right 1.3735664
      89      25     as1t accuracy  left right 1.0601139
      91      26     as1t accuracy  left  left 2.1595326
      93      27     as1t accuracy  left  left 2.0660525
      95      28     as1t accuracy  left  left 1.5216962
      97      29     as1t accuracy  left  left 1.2179661
      99      30     as1t accuracy  left right 1.1660569
      121     31     as1t    speed right right 1.2526878
      123     32     as1t    speed right right 1.3612807
      125     33     as1t    speed right  left 1.1081115
      127     34     as1t    speed right right 0.9069431
      129     35     as1t    speed right right 1.3007431
      131     36     as1t    speed right right 1.1144518
      133     37     as1t    speed right right 1.1973887
      135     38     as1t    speed right right 1.2928114
      137     39     as1t    speed right right 1.0390089
      139     40     as1t    speed right right 1.0903481
      161     41     as1t  neutral right right 2.0526427
      163     42     as1t  neutral right right 1.3163763
      165     43     as1t  neutral right  left 1.4552335
      167     44     as1t  neutral right right 1.3787500
      169     45     as1t  neutral right right 2.5690105
      171     46     as1t  neutral right right 2.3174432
      173     47     as1t  neutral right right 1.9109834
      175     48     as1t  neutral right right 2.1595589
      177     49     as1t  neutral right right 1.7756665
      179     50     as1t  neutral right right 1.3533651
      201     51     as1t accuracy right right 1.4120491
      203     52     as1t accuracy right right 1.2279165
      205     53     as1t accuracy right  left 0.8540025
      207     54     as1t accuracy right right 2.2168005
      209     55     as1t accuracy right right 1.5696998
      211     56     as1t accuracy right  left 2.1831008
      213     57     as1t accuracy right  left 1.6599022
      215     58     as1t accuracy right right 1.4069464
      217     59     as1t accuracy right right 1.1379449
      219     60     as1t accuracy right right 1.4788349
      21       1     bd6t    speed  left  left 1.1431110
      23       2     bd6t    speed  left right 1.1568404
      25       3     bd6t    speed  left  left 1.1525343
      27       4     bd6t    speed  left  left 0.8557349
      29       5     bd6t    speed  left  left 1.2749753
      31       6     bd6t    speed  left  left 0.9191639
      33       7     bd6t    speed  left  left 1.1984658
      35       8     bd6t    speed  left right 1.1423907
      37       9     bd6t    speed  left  left 1.3775807
      39      10     bd6t    speed  left  left 1.8239025
      61      11     bd6t  neutral  left  left 1.7628106
      63      12     bd6t  neutral  left right 1.7751813
      65      13     bd6t  neutral  left  left 1.3634590
      67      14     bd6t  neutral  left  left 1.6609535
      69      15     bd6t  neutral  left  left 1.4420837
      71      16     bd6t  neutral  left  left 2.0303752
      73      17     bd6t  neutral  left  left 2.5350163
      75      18     bd6t  neutral  left  left 2.3413114
      77      19     bd6t  neutral  left  left 1.6085614
      79      20     bd6t  neutral  left  left 1.7790590
      101     21     bd6t accuracy  left  left 1.5055781
      103     22     bd6t accuracy  left  left 1.4714349
      105     23     bd6t accuracy  left  left 1.2633969
      107     24     bd6t accuracy  left  left 1.4022881
      109     25     bd6t accuracy  left right 1.5379162
      111     26     bd6t accuracy  left  left 1.2376553
      113     27     bd6t accuracy  left  left 1.3866129
      115     28     bd6t accuracy  left  left 1.5251978
      117     29     bd6t accuracy  left  left 1.0863382
      119     30     bd6t accuracy  left  left 1.4439774
      141     31     bd6t    speed right  left 1.0662231
      143     32     bd6t    speed right right 1.1241908
      145     33     bd6t    speed right right 0.9894373
      147     34     bd6t    speed right right 1.1363703
      149     35     bd6t    speed right right 1.5094163
      151     36     bd6t    speed right right 1.3607508
      153     37     bd6t    speed right  left 2.0314261
      155     38     bd6t    speed right right 1.2109026
      157     39     bd6t    speed right right 1.5667086
      159     40     bd6t    speed right right 1.4264524
      181     41     bd6t  neutral right  left 1.9479597
      183     42     bd6t  neutral right right 1.5346550
      185     43     bd6t  neutral right right 1.5947515
      187     44     bd6t  neutral right right 2.5420919
      189     45     bd6t  neutral right right 1.5466112
      191     46     bd6t  neutral right  left 1.5146953
      193     47     bd6t  neutral right  left 1.9679805
      195     48     bd6t  neutral right right 1.8031439
      197     49     bd6t  neutral right right 1.5779376
      199     50     bd6t  neutral right  left 1.6736543
      221     51     bd6t accuracy right right 1.6805979
      223     52     bd6t accuracy right right 1.5073002
      225     53     bd6t accuracy right right 1.3659314
      227     54     bd6t accuracy right right 1.6692198
      229     55     bd6t accuracy right right 1.0853028
      231     56     bd6t accuracy right  left 1.1938813
      233     57     bd6t accuracy right right 1.8463224
      235     58     bd6t accuracy right  left 1.1581136
      237     59     bd6t accuracy right right 1.1518188
      239     60     bd6t accuracy right right 1.3311090

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
          trials subjects        E     S     R        rt
      1        1     as1t    speed  left  left 0.9280160
      3        2     as1t    speed  left  left 0.9009727
      5        3     as1t    speed  left  left 1.0128159
      7        4     as1t    speed  left  left 0.9078842
      9        5     as1t    speed  left  left 0.8838565
      11       6     as1t    speed  left  left 0.8349570
      13       7     as1t    speed  left  left 0.8541628
      15       8     as1t    speed  left  left 0.8873348
      17       9     as1t    speed  left  left 1.0874602
      19      10     as1t    speed  left  left 0.9044905
      41      11     as1t  neutral  left right 0.9928104
      43      12     as1t  neutral  left  left 1.4726879
      45      13     as1t  neutral  left  left 1.5136812
      47      14     as1t  neutral  left  left 1.4227448
      49      15     as1t  neutral  left  left 1.3624661
      51      16     as1t  neutral  left  left 1.2007451
      53      17     as1t  neutral  left  left 1.5108485
      55      18     as1t  neutral  left  left 1.1971381
      57      19     as1t  neutral  left  left 1.4237964
      59      20     as1t  neutral  left  left 1.4240631
      81      21     as1t accuracy  left  left 1.1261396
      83      22     as1t accuracy  left  left 1.3806725
      85      23     as1t accuracy  left  left 1.1609113
      87      24     as1t accuracy  left  left 1.3577262
      89      25     as1t accuracy  left right 1.0471066
      91      26     as1t accuracy  left  left 1.1073034
      93      27     as1t accuracy  left  left 1.0989463
      95      28     as1t accuracy  left  left 1.3284815
      97      29     as1t accuracy  left  left 1.0555446
      99      30     as1t accuracy  left  left 1.0741559
      121     31     as1t    speed right right 1.0187293
      123     32     as1t    speed right right 0.9814893
      125     33     as1t    speed right  left 0.8739511
      127     34     as1t    speed right right 0.8109042
      129     35     as1t    speed right right 0.9584358
      131     36     as1t    speed right right 1.1550435
      133     37     as1t    speed right right 0.9387957
      135     38     as1t    speed right right 1.0194077
      137     39     as1t    speed right right 1.2271491
      139     40     as1t    speed right right 0.8929923
      161     41     as1t  neutral right right 1.5256713
      163     42     as1t  neutral right right 2.0518276
      165     43     as1t  neutral right right 1.5871876
      167     44     as1t  neutral right right 1.9652669
      169     45     as1t  neutral right right 1.6067294
      171     46     as1t  neutral right right 1.5536834
      173     47     as1t  neutral right right 1.5169147
      175     48     as1t  neutral right right 1.4231504
      177     49     as1t  neutral right right 1.4704327
      179     50     as1t  neutral right right 1.9840533
      201     51     as1t accuracy right right 1.0845494
      203     52     as1t accuracy right right 1.1063536
      205     53     as1t accuracy right  left 0.7933717
      207     54     as1t accuracy right  left 1.0113856
      209     55     as1t accuracy right right 1.1108572
      211     56     as1t accuracy right right 1.0210001
      213     57     as1t accuracy right right 1.1471407
      215     58     as1t accuracy right right 1.4319051
      217     59     as1t accuracy right right 1.7575233
      219     60     as1t accuracy right right 1.3688643
      21       1     bd6t    speed  left  left 1.0208927
      23       2     bd6t    speed  left  left 0.7744698
      25       3     bd6t    speed  left  left 1.1180789
      27       4     bd6t    speed  left  left 0.7575837
      29       5     bd6t    speed  left  left 1.0709417
      31       6     bd6t    speed  left right 1.1845798
      33       7     bd6t    speed  left right 1.0125895
      35       8     bd6t    speed  left  left 1.0491586
      37       9     bd6t    speed  left  left 0.9465356
      39      10     bd6t    speed  left  left 0.8448592
      61      11     bd6t  neutral  left  left 1.6042024
      63      12     bd6t  neutral  left  left 1.4498685
      65      13     bd6t  neutral  left  left 1.9838021
      67      14     bd6t  neutral  left  left 1.2987800
      69      15     bd6t  neutral  left  left 1.8766217
      71      16     bd6t  neutral  left  left 1.4309855
      73      17     bd6t  neutral  left  left 1.7424045
      75      18     bd6t  neutral  left  left 1.3510569
      77      19     bd6t  neutral  left  left 1.2152129
      79      20     bd6t  neutral  left  left 1.5662376
      101     21     bd6t accuracy  left  left 1.2377503
      103     22     bd6t accuracy  left  left 1.3724097
      105     23     bd6t accuracy  left  left 0.9701447
      107     24     bd6t accuracy  left  left 1.0272444
      109     25     bd6t accuracy  left  left 1.4612821
      111     26     bd6t accuracy  left  left 1.0618854
      113     27     bd6t accuracy  left  left 1.0677462
      115     28     bd6t accuracy  left  left 1.3346883
      117     29     bd6t accuracy  left  left 1.8494961
      119     30     bd6t accuracy  left  left 1.4070013
      141     31     bd6t    speed right right 0.9957356
      143     32     bd6t    speed right right 1.1446996
      145     33     bd6t    speed right right 1.3102483
      147     34     bd6t    speed right  left 0.9437308
      149     35     bd6t    speed right right 1.1165105
      151     36     bd6t    speed right right 0.9776703
      153     37     bd6t    speed right right 0.7328950
      155     38     bd6t    speed right right 1.0700732
      157     39     bd6t    speed right right 1.1634011
      159     40     bd6t    speed right right 0.9587314
      181     41     bd6t  neutral right right 1.2750007
      183     42     bd6t  neutral right right 1.7700724
      185     43     bd6t  neutral right right 1.7163212
      187     44     bd6t  neutral right right 1.6584146
      189     45     bd6t  neutral right right 1.2373955
      191     46     bd6t  neutral right right 1.4939404
      193     47     bd6t  neutral right right 1.6174103
      195     48     bd6t  neutral right right 1.4643984
      197     49     bd6t  neutral right right 1.4014406
      199     50     bd6t  neutral right right 1.7699586
      221     51     bd6t accuracy right right 1.2711660
      223     52     bd6t accuracy right right 1.3422554
      225     53     bd6t accuracy right right 1.2014337
      227     54     bd6t accuracy right right 1.2291637
      229     55     bd6t accuracy right right 1.8717891
      231     56     bd6t accuracy right right 1.2700290
      233     57     bd6t accuracy right right 1.4387550
      235     58     bd6t accuracy right right 1.1099540
      237     59     bd6t accuracy right right 0.9106558
      239     60     bd6t accuracy right  left 1.4617551

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
          trials subjects        E     S     R        rt
      1        1     as1t    speed  left  left 0.6802333
      2        2     as1t    speed  left right 0.3581243
      3        3     as1t    speed  left right 0.8831889
      4        4     as1t    speed  left  left 0.4499327
      5        5     as1t    speed  left right 0.3303770
      6        6     as1t    speed  left  left 0.3646702
      7        7     as1t    speed  left  left 0.6389701
      8        8     as1t    speed  left  left 0.5084518
      9        9     as1t    speed  left  left 0.4733514
      10      10     as1t    speed  left  left 0.2684499
      21      11     as1t  neutral  left  left 1.0005519
      22      12     as1t  neutral  left right 2.2555073
      23      13     as1t  neutral  left  left 0.7989559
      24      14     as1t  neutral  left  left 0.5415385
      25      15     as1t  neutral  left  left 0.3569167
      26      16     as1t  neutral  left  left 0.5587834
      27      17     as1t  neutral  left  left 1.3312508
      28      18     as1t  neutral  left  left 0.3797134
      29      19     as1t  neutral  left right 0.5558977
      30      20     as1t  neutral  left  left 0.5683490
      41      21     as1t accuracy  left  left 0.8586173
      42      22     as1t accuracy  left  left 0.5523000
      43      23     as1t accuracy  left  left 0.7979393
      44      24     as1t accuracy  left  left 0.6868122
      45      25     as1t accuracy  left  left 0.8292358
      46      26     as1t accuracy  left  left 0.4086253
      47      27     as1t accuracy  left right 2.3900265
      48      28     as1t accuracy  left right 0.6865439
      49      29     as1t accuracy  left right 0.9177357
      50      30     as1t accuracy  left  left 0.9141145
      61      31     as1t    speed right right 0.4950596
      62      32     as1t    speed right right 0.4025097
      63      33     as1t    speed right right 0.6047484
      64      34     as1t    speed right right 0.3844647
      65      35     as1t    speed right right 0.3144739
      66      36     as1t    speed right right 0.3884522
      67      37     as1t    speed right  left 0.4383368
      68      38     as1t    speed right  left 0.3667670
      69      39     as1t    speed right  left 0.6913264
      70      40     as1t    speed right right 0.3615293
      81      41     as1t  neutral right right 0.3838494
      82      42     as1t  neutral right right 0.8643106
      83      43     as1t  neutral right right 0.7988730
      84      44     as1t  neutral right right 1.1206340
      85      45     as1t  neutral right right 2.1543090
      86      46     as1t  neutral right right 0.5008352
      87      47     as1t  neutral right right 0.9054355
      88      48     as1t  neutral right  left 1.3175715
      89      49     as1t  neutral right right 0.3347010
      90      50     as1t  neutral right right 0.4547618
      101     51     as1t accuracy right  left 0.5610446
      102     52     as1t accuracy right right 0.5708478
      103     53     as1t accuracy right  left 1.1010207
      104     54     as1t accuracy right  left 2.2889515
      105     55     as1t accuracy right right 0.7317996
      106     56     as1t accuracy right  left 4.9271514
      107     57     as1t accuracy right right 0.2915473
      108     58     as1t accuracy right right 0.4503388
      109     59     as1t accuracy right  left 1.6126562
      110     60     as1t accuracy right  left 0.4261822
      11       1     bd6t    speed  left right 0.3656491
      12       2     bd6t    speed  left  left 0.3693564
      13       3     bd6t    speed  left  left 0.7347282
      14       4     bd6t    speed  left  left 0.4125840
      15       5     bd6t    speed  left right 1.6893529
      16       6     bd6t    speed  left  left 1.0776926
      17       7     bd6t    speed  left  left 0.3417613
      18       8     bd6t    speed  left  left 1.4717633
      19       9     bd6t    speed  left  left 0.3481593
      20      10     bd6t    speed  left  left 0.6959749
      31      11     bd6t  neutral  left  left 0.4481234
      32      12     bd6t  neutral  left  left 2.6101841
      33      13     bd6t  neutral  left  left 0.4612578
      34      14     bd6t  neutral  left  left 0.3420450
      35      15     bd6t  neutral  left  left 0.3711668
      36      16     bd6t  neutral  left right 0.4676017
      37      17     bd6t  neutral  left  left 0.3468602
      38      18     bd6t  neutral  left  left 0.9528720
      39      19     bd6t  neutral  left  left 0.4786007
      40      20     bd6t  neutral  left  left 0.4911437
      51      21     bd6t accuracy  left right 0.6757066
      52      22     bd6t accuracy  left  left 0.3551867
      53      23     bd6t accuracy  left right 1.0540506
      54      24     bd6t accuracy  left  left 0.3110610
      55      25     bd6t accuracy  left  left 0.5502727
      56      26     bd6t accuracy  left  left 2.8111990
      57      27     bd6t accuracy  left right 0.6100624
      58      28     bd6t accuracy  left  left 1.0539004
      59      29     bd6t accuracy  left  left 0.5334544
      60      30     bd6t accuracy  left  left 1.1643393
      71      31     bd6t    speed right  left 0.9775268
      72      32     bd6t    speed right  left 1.2008347
      73      33     bd6t    speed right right 0.6443278
      74      34     bd6t    speed right  left 0.5470503
      75      35     bd6t    speed right right 0.9986496
      76      36     bd6t    speed right right 0.3710058
      77      37     bd6t    speed right  left 1.1326745
      78      38     bd6t    speed right right 0.3531640
      79      39     bd6t    speed right right 0.3785078
      80      40     bd6t    speed right  left 0.8941779
      91      41     bd6t  neutral right right 0.4061499
      92      42     bd6t  neutral right  left 0.4199573
      93      43     bd6t  neutral right  left 0.9842469
      94      44     bd6t  neutral right right 0.7297536
      95      45     bd6t  neutral right  left 0.5216554
      96      46     bd6t  neutral right right 2.4432551
      97      47     bd6t  neutral right right 0.5961414
      98      48     bd6t  neutral right right 0.5885120
      99      49     bd6t  neutral right right 0.5745838
      100     50     bd6t  neutral right right 0.3989805
      111     51     bd6t accuracy right right 1.0897391
      112     52     bd6t accuracy right right 1.6217274
      113     53     bd6t accuracy right right 7.5569366
      114     54     bd6t accuracy right  left 0.7934986
      115     55     bd6t accuracy right right 0.4038146
      116     56     bd6t accuracy right right 0.5446294
      117     57     bd6t accuracy right right 0.8167383
      118     58     bd6t accuracy right right 0.6870093
      119     59     bd6t accuracy right  left 1.5046528
      120     60     bd6t accuracy right  left 0.7933934

