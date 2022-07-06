# MPDS

This repository contains the code used in the experiments of our paper titled _Most Probable Densest Subgraphs_.

## Data

The datasets used in our paper are available (here)[https://drive.google.com/drive/folders/14uzqpv-P0qc9jmxqoZvn14TZrAOAIqss?usp=sharing]. Each dataset is in a text file whose first line contains a single integer `n` denoting the number of nodes, followed by lines of the form `u v p`, each denoting an edge `(u, v)` with probability `p`.

## Code Usage

We show the commands needed to run our various experiments with our provided code.

### Our Proposed Dense Subgraphs

Here we show the commands for returning our proposed dense subgraphs: MPDS and NDS. In each case, the parameter `notion-of-density` can take one of the values `edge`, `clique` and `pattern`, denoting the notion of density. Additional required parameters are `h` for clique density (taking one of the values `3`, `4` and `5`) and `psi` for pattern density (taking one of the values `2-star`, `3-star`, `c3-star` and `diamond`); they are taken as input from the user at runtime. The output written to the specified file consists of subgraphs, one per line, containing the space-separated node set followed by the relevant metric (mentioned below). The running times are displayed on the standard output when the commands are run.

#### MPDS
```
python mpds.py path-to-graph path-to-output number-of-samples number-of-subgraphs notion-of-density
```
Metric in output: (estimated) densest subgraph probability

#### NDS
```
python nds.py path-to-graph path-to-output number-of-samples containment-probability-threshold notion-of-density
```
Metric in output: (estimated) densest subgraph containment probability

#### Exact Algorithms
As mentioned in Section VI-F of our paper, we compare our solution against the exact methods, which are run using the following command:
```
python exact.py path-to-graph path-to-output number-of-subgraphs notion-of-density
```

### Existing Dense Subgraphs

Here we show the commands for returning existing dense subgraphs: EDS, core and truss (Section VI-C of our paper). The output written to the specified file consists of a line containing the space-separated node set of the concerned subgraph.

#### EDS
```
python existing/eds.py path-to-graph path-to-output
```

#### Core
```
python existing/core.py path-to-graph path-to-output eta
```

#### Truss
```
python existing/truss.py path-to-graph path-to-output gamma
```
