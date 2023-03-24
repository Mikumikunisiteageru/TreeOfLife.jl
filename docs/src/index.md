# TreeOfLife.jl

```@docs
TreeOfLife
```

## Types

```@docs
AbstractNode
ChronoNode
CladoNode
AbstractTree
ChronoTree
CladoTree
length
getindex
setindex!
empty
==
```

## Newick format

```@docs
readnewick
writenewick
fromnewick
tonewick
```

#### Internal functions

```@docs
TreeOfLife.ignore_comments
TreeOfLife.from_newick
TreeOfLife.get_tree_type
TreeOfLife.to_newick!
```

## Nexus format

```@docs
readnexus
```

## Methods involving one tree

```@docs
gettips
gettipnames
alldistinct
getage
getages
preorder
postorder
isroot
istip
getname
hassibling
rename
rename!
subtree
getmrca
ismonophyl
phylodiv
cutfromroot
cutfromtips
isbinary
treehash
isomorphic
getdescs
getdescnames
```

#### Internal functions

```@docs
TreeOfLife.calibrate_t_root!
TreeOfLife.calibrate_t_branch!
TreeOfLife.mean_
TreeOfLife.preorder!
TreeOfLife.postorder!
TreeOfLife.get_counts
TreeOfLife.sum_t_branch
```

## Methods involving multiple trees

```@docs
consensus
```

#### Internal functions

```@docs
TreeOfLife.count_clade_ages
TreeOfLife.count_clades
TreeOfLife.construct_tree
```
