# A memetic algorithm based on reformulation local search for minimum sum-of-squares clustering in networks (EMSSC)
This repository includes the source code of the proposed MA algorithm published in an Information Sciences paper titled with "A memetic algorithm based on reformulation local search for minimum sum-of-squares clustering in networks".

The three sets of 132 benchmark instances consist of 40 instances from OR-LIB, 80 instances from TSP-LIB (20 instances with p \in {10, 50, 100, 200} for each instance, giving a total of 80 instances), and 12 instances from University of Florida Sparse Matrix (UFSM) Collection (3 instanes with p \in {10, 50, 100, 200} for each instance, giving a total of 12 instances). To facilitate the further research, we upload the used 132 instances here.

We made comparisons between MA and some state-of-the-art methods from the following related EMSSC works:

[1] E. Carrizosa, N. Mladenovic´ , R. Todosijevic´ , Variable neighborhood search for minimum sum-of-squares clustering on networks, Eur. J. Oper. Res. 230(2) (2013) 356–363.

[2] [34] A. Nikolaev, N. Mladenovic´ , R. Todosijevic´ , J-means and I-means for minimum sum-of-squares clustering on networks, Optim. Lett. 11 (2) (2017) 359–376.

[3] M.G.C. Resende, R.F. Werneck, A hybrid heuristic for the p-median problem, J. Heuristics 10 (1) (2004) 59–88.

Please cite our work as:

Zhou, Q., Benlic, U., & Wu, Q. (2020). A memetic algorithm based on reformulation local search for minimum sum-of-squares clustering in networks. Information Sciences, 541, 271-296.

** Instructions to use the source code of MA

*** To compile:

q.zhou$ make

q.zhou$

*** To run:

q.zhou$ ./MA_EMSSC ./instance_file ./output_res_file 

(where "instance_file is the instance name, output_res_file is a file used to store the running information)

q.zhou$

*** To clean

q.zhou$ make clean

q.zhou$
