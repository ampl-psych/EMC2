# make_design

    Code
      design <- make_design(data = data.frame(forstmann, CO = 1:nrow(forstmann)),
      model = LBA, matchfun = function(d) d$S == d$lR, formula = list(v ~ lM, sv ~ lM,
      B ~ E + lR, t0 ~ E2 + CO), contrasts = list(v = list(lM = matrix(c(-1 / 2, 1 /
      2), ncol = 1, dimnames = list(NULL, "d")))), constants = c(sv = log(1)),
      functions = list(E2 = function(d) factor(d$E != "speed", labels = c("speed",
        "nonspeed"))))
    Message <simpleMessage>
      Parameter(s) A not specified in formula and assumed constant.
    Output
      
       Sampled Parameters: 
       [1] "v"             "v_lMd"         "sv_lMTRUE"     "B"            
       [5] "B_Eneutral"    "B_Eaccuracy"   "B_lRright"     "t0"           
       [9] "t0_E2nonspeed" "t0_CO"        
      
       Design Matrices: 
      $v
          lM v v_lMd
        TRUE 1   0.5
       FALSE 1  -0.5
      
      $sv
          lM sv sv_lMTRUE
        TRUE  1         1
       FALSE  1         0
      
      $B
              E    lR B B_Eneutral B_Eaccuracy B_lRright
          speed  left 1          0           0         0
          speed right 1          0           0         1
        neutral  left 1          1           0         0
        neutral right 1          1           0         1
       accuracy  left 1          0           1         0
       accuracy right 1          0           1         1
      
      $t0
             E2 CO t0 t0_E2nonspeed t0_CO
          speed  0  1             0     0
       nonspeed  0  1             1     0
      
      $A
       A
       1
      

