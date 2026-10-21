# exG likelihood

    Code
      round(model_likelihood(emc, parameters), 6)
    Output
                     1
       [1,] -12.660168
       [2,] -12.904050
       [3,] -21.118810
       [4,] -14.216410
       [5,] -13.666829
       [6,] -16.003091
       [7,]  -4.204066
       [8,] -15.062852
       [9,] -16.365170
      [10,] -31.423929
      [11,] -18.316601
      [12,] -18.995033
      [13,]  -9.912352
      [14,] -12.677239
      [15,]  -8.988904
      [16,] -18.232424
      [17,]  -6.269664

# posterior-only stop-signal plot summaries are stable

    Code
      writeLines(summary_lines)
    Output
      sources: posterior
      datasets: posterior
      labels: [0.1,0.2] | (0.2,0.3] | (0.3,0.4] | (0.4,0.5]
      x_plot: 0.25, 0.50, 0.75, 1.00
      median: 0.410, 0.450, 0.475, 0.505

