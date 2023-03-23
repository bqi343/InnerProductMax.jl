# Inner Product Maximization in 3 Dimensions

## About

Implementations of [this task](https://cstheory.stackexchange.com/questions/34503/maximizing-inner-product) with $d=3$.

## Dev Instructions

Add these to `startup.jl`: `using Revise, ReTest`. Then run these commands from the repo root:

```
julia --project
include("test/InnerProductMaxTests.jl"); retest()
```

## Table Generation

```
include("test/InnerProductMaxTests.jl"); retest("sphere_perf")
```

Results:

```
         N       Q Method  Preproc Time (s)  Preproc Memory (MB)  Query Time (s)
0     1000    1000  Naive          0.000017             0.088208        0.002308
1     1000    1000   Mine          0.010655             8.234688        0.003923
2     2000    2000  Naive          0.000092             0.176208        0.008973
3     2000    2000   Mine          0.023383            17.890672        0.003160
4     5000    5000  Naive          0.000276             0.440208        0.046815
5     5000    5000   Mine          0.081375            48.884608        0.017122
6    10000   10000  Naive          0.000646             0.880208        0.204886
7    10000   10000   Mine          0.321766           102.224256        0.022775
8    20000   20000  Naive          0.001080             1.760208        0.712584
9    20000   20000   Mine          0.370499           210.870496        0.056952
10   50000   50000  Naive          0.003199             4.400208        4.347602
11   50000   50000   Mine          0.841139           562.539104        0.189053
12  100000  100000  Naive          0.008528             8.800208       29.408380
13  100000  100000   Mine          4.712876          1162.170496        0.430899
```