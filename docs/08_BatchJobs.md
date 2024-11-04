# Batch jobs

For compute jobs that take hours or days to run,
instead of sitting at the terminal waiting for the results,
we submit a "batch job" to the queue manager,
which runs the job when resources are available.

On Roar, the queue manager is called SLURM 
(Simple Linux Utility for Resource Management).  
Besides `salloc` for interactive jobs
(see [Interactive jobs](06_AccessingCollab.md#interactive-jobs)),
here are the basic SLURM commands
(with their equivalents under PBS, the previous queue system):

| Command | Effect| PBS equivalent |
| ---- | ---- | ---- |
|`sbatch <script>` | submit batch job `<script>` | `qsub <script>` |
| `squeue -u <userid>` | check on jobs submitted by `<userid>` | `qstat -u <userid>` |
| `scancel <jobID>` | cancel the job | `qdel <jobID>` |

When you execute `sbatch myJob.sh`, SLURM responds with something like
```
Submitted batch job 25789352
```
To check on your jobs, execute `squeue -u <userID>`; SLURM responds with something like
```
JOBID		PARTITION	NAME		USER	ST	TIME	NODES	NODELIST(REASON)
25789352	open 		myJob.sh	abc123	R	1:18:31	1		p-sc-2008
```
Here ST = status:  PD = pending, R = running, C = completed.  
To cancel the job, execute `scancel 25789352`.

## Batch scripts

Jobs submitted to SLURM are in the form of a "batch script".
A batch script is a shell script that executes commands,
with a preamble of comments
that are SLURM directives `#SBATCH ...` to specify:

- an **allocation** to charge for the job;
- a **queue** (qos) to submit the job to;
- a **partition** (type of nodes) to run on;
- nodes, cores, memory, GPUs, and time;
- and other job-related parameters.

An example is:

```
#!/bin/bash
#SBATCH --account=<alloc>
#SBATCH --qos=regular
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=1gb
#SBATCH --time=4:00:00
#SBATCH --job-name=<name>

# as usual, cd to where submitted
cd $SLURM_SUBMIT_DIR

module load gromacs-2019.6  # need this module

SYSTEM=$SLURM_SUBMIT_DIR/System
gmx grompp -f $SYSTEM/nvt.mdp -c $SYSTEM/min.gro -p $SYSTEM/testJob.top -o nvt.tpr 
gmx mdrun -nt 8 -nb cpu -deffnm nvt
```

Everything after the last `#SBATCH` are commands to be executed;
lines with `#` other than `#SBATCH` are ordinary bash script comments.
Most scripts start with `cd $SLURM_SUBMIT_DIR`,
which is the directory from which the job was submitted.

Use `module load` just as you would interactively
to load any modules the batch script needs.
The first line `#!/bin/bash` executes your `.bashrc` file,
and thus performs any initializations it contains
(most importantly, loading modules).

For examples of common resource requests,
see [Batch examples](09_BatchExamples.md).

**The hardware you reserve by choosing a partition,
specifying numbers of nodes and cores, or requesting a GPU,
affects the cost of your job; see [Prices](05_ChargeAccounts.md/#prices).**

## Queues

The directive `#SBATCH --qos=<queue>` submits batch jobs to a queue, 
or QoS = "Quality of Service" in SLURM-speak.
(Queues are like classes of service on an airline flight:
first, business, economy,...)

Collab has six queues:  open, regular, debug, express, and warp.  
Each serves a different purpose, and has different restrictions.

| queue (QOS) | description | restrictions |
| ---- | ---- | ---- |
| open | no-cost access | Portal and vintage hardware only |
| regular | for "normal" jobs | time < 14 days |
| debug	| for testing, debugging, <br> quick analysis | time < 1 hour |
| express | for rush jobs; <br> 1.5x price | time < 14 days |
| warp | for emergency jobs; <br> 3x price | time < 7 days |

Express and warp queues cost more,
but increase a job’s priority so it runs sooner.

## Job output files

By default, batch job standard output and standard error
are both directed to `slurm-%j.out`, where `%j` is the jobID.
But output and error filenames can be customized:
`#SBATCH -e = <file>` redirects standard error to `<file>`,
and ` #SBATCH -o` likewise redirects standard output.

SLURM variables `%x` (job name) and `%u` (username)
are useful for this purpose.  For example,
```
#SBATCH -eo = %u_%x.out
```
writes both standard output and error to `<username>_<jobname>.out`.

## Resource usage

The SLURM command [`sacct`][sacct]
reports the resources used by a completed batch job,
which helps users learn what resources to request next time.
At the bottom of a batch script, the command
[sacct]: https://slurm.schedmd.com/sacct.html

```
sacct -j $SLURM_JOB_ID --format=JobID,JobName,MaxRSS,Elapsed,TotalCPU,State
```
generates a report in the batch output file of resources used.
(`$SLURM_JOB_ID` is a variable that returns the jobID of the batch job.)
As in the example, sacct takes formatting options to control what it prints;
`sacct --helpformat` lists all the options.

## Timing jobs[](#timing-jobs)

It is good practice to test a new workflow
by running small short jobs before submitting big long jobs.
To help plan your compute usage, 
it is helpful to time such test jobs.

Many well-designed applications display timing information
at the end of the log files they generate.
If this is not the case for your application,
you can find out how long a batch job takes
by sandwiching the commands you execute
between [`date`][date] commands:
[date]: https://man7.org/linux/man-pages/man1/date.1.html
```
date
<commands>
date
```
Your batch standard output file will then contain two "timestamps",
from which you can determine the running time.
To time a single command in a batch file, use [`time`][time]:
[time]: https://www.man7.org/linux/man-pages/man1/time.1.html
```
time <command>
```
which will write timing information to standard output.
