## Algorithms for Generating Random Simple Orthogonal Polygons

1. Proposed Algorithm

Run the proposed algorithm using:

```bash
python3 genDataStatistic.py

The experimental settings can be configured in genDataStatistic.py:

VERTEX_COUNTS = [10, 20, 50, 100, 200, 500, 1000, 100000]
TRIALS_PER_SIZE = 50

* VERTEX_COUNTS: list of polygon sizes to be tested.
* TRIALS_PER_SIZE: number of trials for each polygon size.

The program produces two output files:

* experiment_results.txt: summary of the experimental results.
* experiment_results_detailed.txt: detailed results for all trials.

2. `inflate_cut_final.py`  
An updated implementation of the Inflate-Cut algorithm presented in the paper “Quadratic-Time Linear-Space Algorithms for Generating Orthogonal Polygons with a Given Number of Vertices” by Ana Paula Tomás and António Leslie Bajuelos (EWCG 2004). Its core routine is implemented in optimized C code and called from Python to improve computational performance.


`inflate_paste_final.py`  
An implementation of the Inflate-Paste algorithm presented in the same paper by Ana Paula Tomás and António Leslie Bajuelos (EWCG 2004).


