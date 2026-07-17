Alorithm for generating randomly simple orthogonal polygons:

1. RUN OUR ALGORITHM OF RANDOM GENERATING ORTHO POLYGONS
python3 genDataStatistic.py 

	VERTEX_COUNTS = [10, 20, 50, 100, 200, 500, 1000, 100000] # a list of trials
    	TRIALS_PER_SIZE = 50 # mumber of trials

there are two output: experiment_results.txt and experiment_results_detailed.txt


2. COMPARISION WITH INFLATE_CUT OF Tomás & Bajuelos (2004)
File inflate_cut
python3 inflate_cut |tee output.txt

3. `inflate_cut_final.py`  
   An updated implementation of the Inflate-Cut algorithm. Its core routine is implemented in optimized C code and called from Python to improve computational performance.

4. `inflate_paste_final.py`  
   An implementation of the Inflate-Paste algorithm presented in the paper *“Quadratic-Time Linear-Space Algorithms for Generating Orthogonal Polygons with a Given Number of Vertices”* by Ana Paula Tomás and António Leslie Bajuelos (EWCG 2004).


