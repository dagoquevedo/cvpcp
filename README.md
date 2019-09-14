# Heuristic for the Capacitated Vertex <em>p</em>-Center Problem

This repository contains the computational implementations of the heuristic proposed in ["Improving the quality of heuristic solutions for the capacitated vertex <em>p</em>-center problem through iterated greedy local search with variable
neighborhood descent"](https://doi.org/10.1016/j.cor.2014.12.013), you can find there a formal definition and a experimental approach about this method.

## Prerequisites

This implementation was developed in C++ and run under a OS like Linux or Unix, you must verify that you have install the following tools:

* GNU C++
* Make
* Git

If you don't have this tools, then install the packages `build-essential` and `git`.

## Compilation

Use git to download the latest version of the package to your home `git` folder (or elsewhere), and then make it:

```
$ cd ~/git
$ git clone https://github.com/dagoquevedo/cvpcp
$ make
```

The `make` instruction will compile the code generating an executable with a name `CVPCP`.

## Dataset format

This implementation accept read only five possible dataset formats, which are described below. If you want to run your own set you must transform your dataset to one of these formats.

### Format A

Set proposed in [Beasley, 1990](https://doi.org/10.1057/jors.1990.166). This set instances is available [here](/Datasets/Beasley) and the format is:

```
Format A: Beasley
----------------------------------------
1st line	: set instance n p best_lb
2nd line	: capacities 
Other lines	: node x y demand
----------------------------------------
```
Here `set = 1`. The computational implementation used euclidean distances calculated from these coordinates but rounded to the nearest integer.

### Format B

Set proposed in [Scaparra <em>et al.</em>, 2004](https://doi.org/10.1002/net.20000) and based on the
dataset of [ReVelle <em>et al.</em>, 2005](https://doi.org/10.1016/j.ejor.2003.11.032). This set instances is available [here](/Datasets/GalvaoReVelle) and the format is:
```
Format B: Galvão and ReVelle
----------------------------------------
1st line	: set instance n p best_lb
Other lines	: Capacities for each node 
Other lines	: Demands for each node
Other lines	: Matrix n x n of the distance between any two nodes i and j.
----------------------------------------
```
Here `set = 2`.

### Format C

Set proposed in [Lorena and Senne, 2004](https://doi.org/10.1016/S0305-0548(03)00039-X). This set instances is available [here](/Datasets/Lorena) and the format is:

```
Format C: Lorena and Senne
----------------------------------------
1st line        : set instance n p best
Other lines     : x y capacities demands
----------------------------------------
```

Here `set = 3, 8`. The computational implementation used euclidean distances calculated from these coordinates but rounded to the nearest integer.


### Format D

Set proposed in [Ceselli and Righini, 2005](https://doi.org/10.1002/net.20059). This set instances is available [here](Datasets/OR-Library/) and the format is:

```
Format D: Galvão and ReVelle
----------------------------------------
1st line	: set instance n p best_lb
Other lines	: Matrix n x n of the distance between any two nodes i and j.
Other lines	: Demands for each node
Other lines	: Capacities for each node 
----------------------------------------
```

Here `set = 4, 5, 6, 7`.


## Execution

For execute the heuristic method, run the following:

`$ ./CVPCP {instance} {r_max} {alpha} {q}`

Where,

* `{instance}`: file path of an instance with a valid format, defined [here](/README.md#Dataset-format)
* `{r_max}`: maximum number of iterations
* `{alpha}`: in a perturbation function, represent the percent of nodes to be disconnected from a solution
* `{q}`: in a shake function, generated the <img src="https://latex.codecogs.com/gif.latex?q" /> subsets nearest centers to <img src="https://latex.codecogs.com/gif.latex?i" />, we recommend a value <img src="https://latex.codecogs.com/gif.latex?q=\lceil\ln(p)\rceil+1" />.


### Output

The execution generate a output in a single line with the following information:

`$ [set] [instance] [n] [p] [best_lb] [incumbent] [gap] [time] [memory] [feasible]`

Where,

* `[set]`: number of the set
* `[instance]`: is the identificator of instance
* `[n]`: number of nodes
* `[p]`: number of centers
* `[best_lb]`: best known lower bound
* `[incumbent]`: best value found by the heuristic
* `[gap]`: percent of relative deviation or gap with respect to the best known lower bound
* `[time]`: execution time in seconds
* `[memory]`: maximum memory used
* `[feasible]`: 1 if a feasible solution, 0 in other case

### Example

Execution:

`$ ./CVPCP instance.dat 1000 0.4 3`

Output:

`$ 4 1 50 5 29.00 29.00 0.00 0.3580 917504 1`

## Citation

D. R. Quevedo-Orozco and R. Z. Ríos-Mercado, Improving the quality of heuristic solutions for the capacitated vertex <em>p</em>-center problem through iterated greedy local search and variable neighborhood descent, <em>Computers and Operations Research</em>, 62 (2015), 133-144. [doi.org/10.1016/j.cor.2014.12.013](https://doi.org/10.1016/j.cor.2014.12.013)

## Contact

* dago@yalma.fime.uanl.mx
* roger@yalma.fime.uanl.mx
