# compare

    Code
      compare(list(samplers_LNR), cores_for_props = 1)
    Output
        MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      1 27   1 -619    1 -603     1         16  -635  -646 -650

# savage-dickey

    Code
      round(savage_dickey(samplers_LNR, parameter = "m", do_plot = F, H0 = -1), 3)
    Output
      [1] 0.037

---

    Code
      round(savage_dickey(samplers_LNR, fun = function(d) d["m"] - d["m_lMd"], H0 = -
      0.5, do_plot = F), 3)
    Output
      [1] 0.099

