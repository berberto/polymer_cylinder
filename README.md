## Sampling random walks in confied geometries

Scripts for simulations and data analysis for

	Adorisio, M., Pezzotta, A., de Mulatier, C., Micheletti, C., & Celani, A. (2017). ''Exact and Efficient Sampling of Conditioned Walks''. ![*J. Stat. Phys.*, **170** (1), 79–100.](https://doi.org/10.1007/s10955-017-1911-y)

*Authors: Matteo Adorisio and Alberto Pezzotta*


#### To compile and use

1. compile auxiliary programs for random numbers, and go to the main file directory
```bash
	cd devel/random && make && cd ../../main
```
2. compile and link all programs in the tree
```bash
	make
```
3. run on a single processor via `./cylinder`, or `mpirun cylinder` for multiple paralle runs. See `run.py` script for running with different parameter sets.
