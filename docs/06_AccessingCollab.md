# Accessing Collab

To use the cluster beyond the several applications provided by the Portal,
logging on via a terminal is necessary.
The Portal provide a desktop environment
(the "RHEL 8 Interactive Desktop")
which includes a terminal application,
with a "look and feel" familiar to Windows users.

## Access via SSH

Alternatively, you can access the cluster via [SSH][ssh]
from a terminal application running on your laptop.
On Windows, use puTTY, available [here](https://www.putty.org);
on Mac, Terminal (installed by default), 
or iTerm  (an improved version), available [here](https://iterm2.com).
[ssh]: https://linux.die.net/man/1/ssh

From the terminal, log on to a submit node as
```
ssh <user>@submit.hpc.psu.edu
```

You will be prompted for your password, 
and then for multi-factor authentication (MFA), 
which confirms that you are you.  To set up MFA, 
go [here](https://accounts.psu.edu/2fa).

**Note:  for security reasons, Roar Restricted can only be accessed
via the [RR Portal](16_RoarRestricted.md).**

## X Window apps

To use any application that "opens a window"
(called an  "X Window" or "X11" application), 
you need an additional program on your laptop.
On the Mac, this is XQuartz, available [here](https://www.xquartz.org).
On the PC, you need VcXsrv, available [here](https://sourceforge.net/projects/vcxsrv/).

To use cluster applications that open windows, log on with
option `-X` for "X forwarding":
```
ssh -X <user>@submit.hpc.psu.edu
```
When you log on to Collab from somewhere off campus,
X Window apps can sometimes be slow;  
Portal works better in such circumstances.

## Interactive jobs

Sometimes, we want to run a compute-intensive job interactively;
for example, a heavy MATLAB, Mathematica, or Gaussian calculation.  

Submit nodes are not designed for such work.
If you try to run a number-crunching program on a submit node,
it will a) run slow, and b) tie up the submit node for everyone else.

Instead, from a submit node, start an "interactive batch" session:
```
salloc -N 1 -n 4 -A <alloc> -t 1:00:00
```
Option `-N` requests the number of nodes,
`-n` the number of cores per node, and `-t` the total time;
`-A <alloc>` specifies your allocation 
(`-A open` for the open queue).

Once salloc starts, you can run compute-intensive jobs 
interactively on the resources you requested.
