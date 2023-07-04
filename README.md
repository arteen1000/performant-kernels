# performance-lab

Scored best in class of 250 competitors. Link to relevant [competition](https://polyarch.github.io/cs33/labs/lab4/). Achieved energy-delay improvement of 580.94 (geomean speedup 21.33 and geomean energy-use improvement of 27.24) using purely cache and compiler optimizations for an Intel(R) Xeon(R) CPU E620 @ 2.40 GHz with L1d cache (32K) and L2 cache (256K). Implementation in `kernels.c`.

Here's the actual performance details:
```
T#  Kernel |    Ni    Nj      Nk  S | Time(ms)    CPE #Inst(M) #VecInst(M) |  #L1(M)  #L2(M)  #L3(M) #Mem(M) MemEn(j) | Speedup Energyup
 0 Transp. | 10000 10000       -  - |     71.9  1.438    186.8         0.0 |    82.0     5.5    4.13    0.01  0.00680 |     7.4     13.2
 1 Stencil |   128   128     128  8 |    129.9  0.242    606.3       268.7 |   311.0     6.6    2.17    0.04  0.01772 |    32.7     19.2
 2 Stencil |    64    64      64 20 |    289.8  0.276   1240.6       520.6 |   624.9     9.4    0.44    0.03  0.03270 |    27.6     38.8
 3 Stencil |     4     4 1048576  2 |     26.7  0.398     70.7        32.6 |    29.3     1.9    0.60    0.06  0.00251 |    31.1     55.8
Geomean Speedup: 21.33
Geomean Energyup: 27.24
Energy-Delay Improvement: 580.94
```

Find related histogram optimization project [here](https://github.com/arteen1000/threaded-histograms).
