# Heuristic for the Capacitated Vertex <em>p</em>-Center Problem

This repository contains the computational implementations of the heuristic proposed in ["Improving the quality of heuristic solutions for the capacitated vertex <em>p</em>-center problem through iterated greedy local search with variable
neighborhood descent"](https://doi.org/10.1016/j.cor.2014.12.013), you can find there a formal definition and a experimental approach about this method.

## Prerequisites

This implementation was developed in C++ and run under a Unix-like OS. You must verify that you have install the following tools:

* GNU C++
* Make
* Git

If you don't have this tools, then install the packages `build-essential` and `git`.

## Compilation

Use git to download the latest version of the package to your home `git` folder (or elsewhere), and then make it:

```
$ cd ~/git
$ git clone https://github.com/dagoquevedo/cvpcp
$ cd cvpcp
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

`$ ./CVPCP {file} {r_max} {alpha} {q}`

Where,

|  Parameter |                                          Description                                          |
|----------|---------------------------------------------------------------------------------------------|
| `{file}` | Instance file path with a valid format, defined here                                    |
| `{r_max}`    | Maximum number of iterations                                                                  |
| `{alpha}`    | Percent of nodes to be disconnected from a solution |
| `{q}`        | Generate <img src="https://latex.codecogs.com/gif.latex?q" /> subsets nearest centers to <img src="https://latex.codecogs.com/gif.latex?i" />, we recommend a value <img src="https://latex.codecogs.com/gif.latex?q=\lceil\ln(p)\rceil+1" />         |
| `{output}`    | Optional. The file path where write the solution|

### Output information

The execution report a output with the following relevant information:

`$ [set] [instance] [n] [p] [best_lb] [incumbent] [gap] [time] [memory] [feasible]`

Where,

|  Output  |                                Description                               |
|-----------|------------------------------------------------------------------------|
| `[set]`       | Set number                                                               |
| `[instance]`  | Instance number                                                          |
| `[n]`         | Number of nodes                                                          |
| `[p]`         | Number of centers                                                        |
| `[best_lb]`   | Best known lower bound                                                   |
| `[incumbent]` | Best value found by the heuristic                                        |
| `[gap]`       | Percent of relative deviation with respect to the best known lower bound |
| `[time]`      | Execution time in seconds                                                |
| `[memory]`    | Maximum memory used                                                      |
| `[feasible]`  | 1 if a feasible solution, 0 in other case                                |

### Output solution

If a path output file was especificated in `{output}` parameter, the application will generate a file with the following format:

```
Format output
----------------------------------------
1st line        : set instance n p
Other lines     : k-center node | nodes assigned to k-center node
----------------------------------------
```

## Citation

D. R. Quevedo-Orozco and R. Z. Ríos-Mercado, Improving the quality of heuristic solutions for the capacitated vertex <em>p</em>-center problem through iterated greedy local search and variable neighborhood descent, <em>Computers and Operations Research</em>, 62 (2015), 133-144. [doi.org/10.1016/j.cor.2014.12.013](https://doi.org/10.1016/j.cor.2014.12.013)

## Contact

* dago@yalma.fime.uanl.mx
* roger@yalma.fime.uanl.mx
