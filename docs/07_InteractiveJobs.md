# Interactive jobs

Sometimes, we want to run a compute-intensive job interactively;
for example, a heavy MATLAB, Mathematica, or Gaussian calculation.  

Submit nodes are not designed for such work.
If you try to run a number-crunching program on a submit node,
it will a) run slow, and b) tie up the submit node for everyone else.

Instead, start an "interactive batch" session,
requested from a submit node like this:

```
salloc -N 1 -n 4 -A <alloc> -G 1 -t 1:00:00
```

In the above example,

- `N` = number of nodes
- `n` = cores per node
- `A` = account (`<alloc>` is an allocation, or "open")
- `G` = number of GPUs (omit if not needed)
- `t` = time requested

Once the `salloc` starts, you have the requested resources,
and can run compute-intensive jobs interactively.
