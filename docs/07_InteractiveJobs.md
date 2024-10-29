# Interactive jobs

Sometimes, we want to run a compute-intensive job interactively;
for example, a heavy MATLAB, Mathematica, or Gaussian calculation.  

Submit nodes are not designed for such work.
If you try to run a number-crunching program on a submit node,
it will a) run slow, and b) tie up the submit node for everyone else.

Instead, start an "interactive batch" session,
requested from a submit node like this:

```
salloc -N 1 -n 4 -A <alloc> -t 1:00:00
```

In the above example, option `-N` requests the number of nodes,
`-n` the number of cores per node, and `-t` the total time;
`-A <alloc>` specifies your allocation 
(`-A open` for the open queue).

Once the `salloc` starts, you can run compute-intensive jobs 
interactively on the resources you requested.
Typically, `salloc` jobs are short-duration (because interactive)
and single-node (so a basic node is fine).