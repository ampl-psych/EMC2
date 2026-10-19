# groups

    Code
      print(group_des)
    Output
      v_S1 ~ age + group 
      a ~ age 

---

    Code
      summary(group_des)
    Output
      
       Sampled Parameters: 
      [1] "v_S1"        "v_S1_age"    "v_S1_group1" "a"           "a_age"      
      
       Design Matrices: 
      $v_S1
                  age group v_S1   v_S1_age v_S1_group1
      as1t -0.4705266     A    1 -0.4705266  -0.7071068
      bd6t -0.4178950     B    1 -0.4178950   0.7071068
      bl1t -0.3652634     A    1 -0.3652634  -0.7071068
      
      $a
                  age a      a_age
      as1t -0.4705266 1 -0.4705266
      bd6t -0.4178950 1 -0.4178950
      bl1t -0.3652634 1 -0.3652634
      

---

    Code
      credint(pri, map = TRUE)
    Output
      $mu
                    2.5%    50%  97.5%
      v_Sleft     -2.658  0.039  2.554
      v_Sright    -2.851 -0.055  2.898
      a_Eaccuracy  0.083  1.395 70.144
      a_Eneutral   0.090  1.395 81.679
      a_Espeed     0.150  1.205 11.473
      t0           0.161  1.200 16.139
      

---

    Code
      credint(pri, selection = "beta")
    Output
      $beta
                    2.5%    50% 97.5%
      v           -1.842 -0.018 1.717
      v_S1        -2.163  0.000 1.951
      v_S1_age    -2.058 -0.013 2.043
      v_S1_group1 -2.027 -0.095 2.104
      a           -1.875  0.041 1.879
      a_age       -1.746 -0.065 1.870
      a_Eneutral  -1.932  0.036 1.941
      a_Eaccuracy -2.093  0.022 2.097
      t0          -2.079  0.003 1.948
      

---

    Code
      round(cbind(model_likelihood(emc, parameters)), 6)
    Output
                  [,1]
      as1t -168.111196
      bd6t   50.724485
      bl1t -111.293013
      hsft  435.177221
      hsgt  431.409756
      kd6t -397.322409
      kd9t  396.859774
      kh6t   10.495294
      kmat  179.417003
      ku4t -232.197437
      na1t  480.291304
      rmbt   -9.524265
      rt2t  335.145715
      rt3t -341.556775
      rt5t -673.368679
      scat   -5.240631
      ta5t  239.577548
      vf1t   34.323981
      zk1t -120.740447

