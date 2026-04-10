# compare

    Code
      compare(list(samples_LNR), cores_for_props = 1)
    Output
          MD wMD  DIC wDIC BPIC wBPIC EffectiveN meanD Dmean minD
      1 -571   1 -607    1 -566     1         41  -648  -656 -690

# savage-dickey

    Code
      round(hypothesis(samples_LNR, parameter = "m", do_plot = F, H0 = -1), 2)
    Output
      [1] 0.05

---

    Code
      round(hypothesis(samples_LNR, fun = function(d) d["m"] - d["m_lMd"], H0 = -0.5,
      do_plot = F), 2)
    Output
      [1] 0.09

