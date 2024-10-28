# Hardware

A cluster consists of multiple nodes.
A node is basically a single computer, 
roughly comparable to a powerful desktop machine.
Some nodes are networked together with fast connections (Infiniband)
that enable large jobs to run in parallel on multiple nodes.
Finally, some nodes include GPUs (graphical processing units),
which can accelerate certain compute jobs.

Compute nodes in a cluster 
are connected to a set of "submit nodes"
where users login and submit compute jobs,
and also to a central filesystem that stores files.


## Collab

Nodes on Collab are of seven different types:

- vintage -- older hardware, accessible at no cost via the open queue.
- basic -- CPU nodes, without Infiniband, for jobs that fit on a single node.
- standard -- CPU nodes, with Infiniband, for single-node or multinode jobs.
- high-memory -- CPU nodes with extra memory, for memory-intensive jobs.
- GPU -- standard nodes, with one powerful GPU.
- GPU 2 -- older CPU nodes, with two or more older GPUs.
- interactive -- nodes with graphics cards, that service the Portal.

Each node type consists of different hardware, 
appropriate to its purpose:

| Resource | Cores | Memory (GB) | CPU | GPU | Network | Count |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| Vintage | 24 | 128 | E5-2650v4 | | ethernet | 240 |
| Basic | 64 | 256 | Gold 6430 | | ethernet | 120 |
| Standard | 64 <br> 48 <br> 48 | 384 <br> 384 <br> 512 | EPYC 9354 <br> Gold 6248R <br> Gold 6342 | | Infiniband | 36 <br> 168 <br> 160 |
| Hi-memory | 48 <br> 28 | 1024 <br> 1024 | Gold 6342 <br> E7-4830v4 | | ethernet | 8 <br> 2 |
| GPU | 48 <br> 24 <br> 24 | 384 <br> 768 <br> 512 | Gold 6248R <br> Gold 6132 <br> E5-2680v3 | dual A100 <br> quad V100 <br> V100 | Infiniband | 38 <br> 2 <br> 2 |
| GPU2 | 28 <br> 28 | 256 <br> 512 | E5-2680v4 <br> E5-2680v4 | P100 <br> P100 | Infiniband | 76 <br> 8 |
| Interactive | 36 | 512 | Gold 6354 | A40 | ethernet | 12 |

## Restricted

Roar Restricted consists of two different node types:

| Resource | Cores | Memory (GB) | CPU | GPU | Network | Count |
| ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| Standard | 48 <br> 24 | 384 <br> 256 | Gold 6248R <br> E5-2680v3 | | Infiniband | 12 <br> 48 |
| GPU | 28 | 256 | E5-2680v4 | quad P100 | Infiniband | 3 |