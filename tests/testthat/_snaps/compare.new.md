# compare

    Code
      compare(list(samplers_LNR))
    Output
        MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      1 27   1 -619    1 -603     1         16  -635  -646 -650

# savage-dickey

    Code
      savage_dickey(samplers_LNR, parameter = "m", do_plot = F, H0 = -1)
    Output
      [1] 0.03657688

---

    Code
      savage_dickey(samplers_LNR, fun = function(d) d["m"] - d["m_lMd"], H0 = -0.5,
      do_plot = F)
    Output
      [1] 0.09922073

