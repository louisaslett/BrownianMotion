# print() shows the skeleton, layers, jumps and labels

    Code
      print(bm)
    Output
      -- Brownian motion -------------------------- 1 dimension --
      1 point at t = 0, no layers
      
         t  W_t  origin
      1  0    0  fixed

---

    Code
      print(bm)
    Output
      -- Brownian motion ----------------------------------- 1 dimension --
      5 points on [0, 4], 3 layers (75% of time), 1 jump, 1 user label
      
         t      W_t    W_tm  origin  layer                   bounds  labels
      1  0   0.0000          fixed   bessel        [-2.035, 1.0000]
      2  1  -1.0351          sim     bessel        [-2.035, 0.4532]
      3  2  -0.5468          sim     intersection  [-1.547, 1.3164]  peak
      4  3   1.3164  0.3164  sim
      5  4   1.8839          sim
      i Layers and bounds apply from each point to the next.

---

    Code
      print(bm, digits = 2)
    Output
      -- Brownian motion --------------------------- 1 dimension --
      5 points on [0, 4], 3 layers (75% of time), 1 jump, 1 user label
      
         t    W_t  W_tm  origin  layer               bounds  labels
      1  0   0.00        fixed   bessel        [-2.0, 1.00]
      2  1  -1.04        sim     bessel        [-2.0, 0.45]
      3  2  -0.55        sim     intersection  [-1.5, 1.32]  peak
      4  3   1.32  0.32  sim
      5  4   1.88        sim
      i Layers and bounds apply from each point to the next.

# print() shows the first and last points of long paths

    Code
      print(bm)
    Output
      -- Brownian motion -------------------------- 1 dimension --
      31 points on [0, 30], no layers
      
           t      W_t  origin
       1   0   0.0000  fixed
       2   1  -1.2195  sim
       3   2  -0.9162  sim
       4   3  -1.2467  sim
       5   4  -2.6499  sim
       :  21 more points (use print(n = ...) to see more)
      27  26  -1.5562  sim
      28  27  -2.7430  sim
      29  28  -2.6648  sim
      30  29  -2.4821  sim
      31  30  -1.4619  sim

---

    Code
      print(bm, n = 3)
    Output
      -- Brownian motion -------------------------- 1 dimension --
      31 points on [0, 30], no layers
      
           t     W_t  origin
       1   0   0.000  fixed
       2   1  -1.220  sim
       :  28 more points (use print(n = ...) to see more)
      31  30  -1.462  sim

# print() shows localised and bridge layers with their bounds

    Code
      print(bm, n = Inf)
    Output
      -- Brownian motion --------------------------- 1 dimension --
      11 points on [0, 3], 9 layers (66.7% of time)
      
               t      W_t  origin   layer                    bounds
       1  0.0000   0.0000  fixed    localised-bb  [-1.742,  0.5000]
       2  0.3609  -0.7417  passage  localised-bb  [-1.826, -0.2417]
       3  0.5308  -0.8264  passage  localised-bb  [-1.610, -0.3264]
       4  0.6133  -0.6102  passage  localised-bb  [-2.121, -0.1102]
       5  0.9071  -2.1211  passage  bessel-bb     [-2.691, -1.8378]
       6  1.0000  -2.3742  sim
       7  2.0000  -2.3109  sim      localised     [-2.811, -1.8109]
       8  2.2641  -2.8109  passage  localised     [-3.311, -2.3109]
       9  2.4644  -3.3109  passage  localised     [-3.811, -2.8109]
      10  2.8368  -3.8109  passage  bessel        [-4.311, -3.3109]
      11  3.0000  -3.7641  sim
      i Layers and bounds apply from each point to the next.

# print() shows multi-dimensional paths and their transformed view

    Code
      print(bm)
    Output
      -- Brownian motion ------------------------- 2 dimensions --
      3 points on [0, 2], 2 layers (100% of time), 1 user label
      Covariance: unit variances, correlation 0.5
      
         t    W_t[1]   W_t[2]  origin  layer   labels
      1  0   0.00000   0.0000  fixed   bessel
      2  1  -0.01437   0.1555  sim     bessel
      3  2  -0.57807  -0.9333  sim             last
      i Layers are those of the components of bm$Z, from each point to the next.

---

    Code
      print(bm$Z)
    Output
      -- Brownian motion, transformed view Z ----- 2 dimensions --
      3 points on [0, 2], 2 layers (100% of time), 1 user label
      Independent standard components, with W_t = Z_t %*% chol(cov)
      
         t    Z_t[1]   Z_t[2]  origin  layer   labels
      1  0   0.00000   0.0000  fixed   bessel
      2  1  -0.01437   0.1878  sim     bessel
      3  2  -0.57807  -0.7439  sim             last
      i Layers apply from each point to the next.

---

    Code
      print(bm3)
    Output
      -- Brownian motion ------------------------- 3 dimensions --
      2 points on [0, 1], no layers, 1 jump
      Covariance: independent, variances 1, 2, 3
      
         t   W_t[1]  W_t[2]  W_t[3]   jump  origin
      1  0   0.0000  0.0000   0.000         fixed
      2  1  -0.7067  0.7539   2.887  1.732  sim

# summary() of a 1-dimensional path

    Code
      print(s)
    Output
      -- Brownian motion summary ------------------------------ 1 dimension --
      Time      [0, 4.157], length 4.157
      Points    6 (1 fixed, 4 simulated, 1 first passage)
      Spacing   min 0.1574, median 1, max 1
      Jumps     2, at t = 2, 3 (sizes 1, -1)
      Labels    peak (1 point)
      Defaults  refine = TRUE, mult = 1, prefer = "bessel"
      
      Path at the skeleton points
                   t     W_t
        start  0.000   0.000
        end    4.157  -7.652
        min    4.157  -7.652
        max    0.000   0.000
      
      Layers: 4 intervals, covering 75.9% of the time
        type          intervals    time  mean width
        bessel                2  2.0000       3.933
        intersection          1  1.0000       2.915
        localised             1  0.1574       2.000
        Not covered  [2, 3]
        Over the layered time the whole path lies in [-7.652, 1].
        1 first passage requested with first_passage().

# summary() of a path with no layers or a single point

    Code
      print(s)
    Output
      -- Brownian motion summary ------------------------------ 1 dimension --
      Time      t = 0
      Points    1 (1 fixed)
      Spacing   -
      Jumps     none
      Labels    none
      Defaults  refine = TRUE, mult = 1, prefer = "bessel"
      
      Path at the skeleton points
               t  W_t
        start  0    0
        end    0    0
        min    0    0
        max    0    0
      
      Layers
        none

# summary() of a multi-dimensional path and its transformed view

    Code
      print(s)
    Output
      -- Brownian motion summary ----------------------------- 3 dimensions --
      Time      [0, 4], length 4
      Points    5 (1 fixed, 4 simulated)
      Spacing   min 1, median 1, max 1
      Jumps     none
      Labels    mid (1 point)
      Defaults  refine = TRUE, mult = 1, prefer = "bessel"
      
      Covariance
        variances 1, 2, 1, correlations from 0.2 to 0.3536
             [,1] [,2] [,3]
        [1,]  1.0  0.5  0.2
        [2,]  0.5  2.0  0.3
        [3,]  0.2  0.3  1.0
      
      Components at the skeleton points
              start     end     min     max     layer bounds
        W[1]      0   0.012  -1.452  0.8496  [-2.452, 1.850]
        W[2]      0  -1.291  -1.291  2.1395  [-2.549, 4.375]
        W[3]      0  -2.888  -2.888  0.0000  [-3.409, 1.611]
      
      Layers: 3 intervals, covering 75% of the time
        Intervals of each type, in each independent component of Z:
        type          Z[1]  Z[2]  Z[3]
        intersection     3     3     3
        Not covered  [3, 4]
        Over the layered time the whole path lies within the layer bounds above.

---

    Code
      print(z)
    Output
      -- Brownian motion summary, transformed view Z --------- 3 dimensions --
      Time      [0, 4], length 4
      Points    5 (1 fixed, 4 simulated)
      Spacing   min 1, median 1, max 1
      Jumps     none
      Labels    mid (1 point)
      Defaults  refine = TRUE, mult = 1, prefer = "bessel"
      
      Covariance
        independent standard components, with W_t = Z_t %*% chol(cov)
      
      Components at the skeleton points
              start      end      min     max     layer bounds
        Z[1]      0   0.0120  -1.4516  0.8496  [-2.452, 1.850]
        Z[2]      0  -0.9807  -0.9807  1.9295  [-1.000, 2.929]
        Z[3]      0  -2.8328  -2.8328  0.0000  [-2.967, 1.000]
      
      Layers: 3 intervals, covering 75% of the time
        Intervals of each type, in each independent component of Z:
        type          Z[1]  Z[2]  Z[3]
        intersection     3     3     3
        Not covered  [3, 4]
        Over the layered time the whole path lies within the layer bounds above.

