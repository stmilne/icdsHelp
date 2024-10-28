# Batch jobs

For compute jobs that take hours or days to run,
instead of sitting at the terminal waiting for the results,
we submit a "batch job" to the queue manager,
which runs the job when resources are available.

On Roar, the queue manager is called SLURM 
(Simple Linux Utility for Resource Management).
Besides `salloc` for interactive jobs (see [Interactive jobs](07_InteractiveJobs.md)),
the basic SLURM commands are:

- `sbatch <script>` Submits `<script>` to SLURM.
- `squeue -u <userid>` Checks on jobs submitted by ``<userid>.
- `scancel <jobid>` Cancels the job.

## Batch scripts

Jobs submitted to SLURM are in the form of a "batch script".
A batch script is a shell script that executes commands,
with a preamble of comments of the form `#SBATCH ...`
that request resources (nodes, cores, memory, GPUs, and time),
specify an allocation to charge for the job,
and set job-related parameters.

A simple example is:

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --gpus=1
#SBATCH --mem=1gb
#SBATCH --time=4:00:00
#SBATCH --job-name=<name>
#SBATCH --account=<account>
#SBATCH --partition=xxx

module load gromacs-2019.6

cd $SLURM_SUBMIT_DIR
SYSTEM=$SLURM_SUBMIT_DIR/System

gmx grompp -f $SYSTEM/nvt.mdp -c $SYSTEM/min.gro -p $SYSTEM/testJob.top -o nvt.tpr 
gmx mdrun -nt 8 -gpu_id 0 -deffnm nvt
```

Here, everything after the last `#SBATCH` are commands to be executed
(load a module for simulation app Gromacs, `cd` to the submit directory, 
define a variable name,
run the Gromacs preprocessor `grompp`, and run a simulation with `mdrun`).

The first line, `#!/bin/bash`, has the effect of running the script "under bash"
(the default shell on Roar), and loading your `.bashrc` initialization file,
so that whatever initializations you have placed there are executed.

Notes:

- `ntasks` = number of cores
- `partition=sla_prio` if xx
- for open queue, `account=open", and omit `partition`
- `$SLURM_SUBMIT_DIR` = the directory in which you submitted the job

## SLURM directives

SLURM has many directives to request resources
and govern the behavior of batch jobs:

| Resource Directive | Description |
| ---- | ---- |
| `-J` or `--job-name` | Job name |
| `-A` or `--account` | Allocation to charge |
| `-p` or `--partition` | Partition (queue) to use |
| `-N` or `--nodes` | Number of nodes |
| `-n` or `--ntasks` | Number of tasks (cores) |
| `--ntasks-per-node` | Alternately, tasks per node |
| `--mem` | Memory per node |
| `--mem-per-cpu` | Alternately, memory per CPU |
| `-t` or `--time` | Maximum run time |
| `-C` or `--constraint` | Required node features |
| `-e` or `--error` | Redirect standard error to file |
| `-o` or `--output` | Redirect standard output to file |

Slurm defines environment variables that the batch script can access.
The most important of these is `SLURM_SUBMIT_DIR`,
which is the directory from which the job was submitted.
(Typically, the first line of a batch script is `cd $SLURM_SUBMIT_DIR`.)

For more information, see documentation on 
[salloc](https://slurm.schedmd.com/salloc.html) 
and [sbatch](https://slurm.schedmd.com/sbatch.html).


## Job output files

Interactive commands write "standard output" and "standard error" to the terminal.
By default, batch job standard output and standard error
are directed `slurm-%j.out`, where the `%j` is the job ID
(visible when the job is submitted, 
or available from `squeue` when the job is queued or running).
But output and error filenames can be customized, 
using the following symbols:

| Symbol | Description |
| :----: | ---- |
| `%j` | Job ID |
| `%x` | Job name |
| `%u` | Username |

For example,
```
#SBATCH -eo = %u_%x.out
```
writes both standard output and error to `<username>_<jobname>.out`.


## [Queues](#queues)

Batch jobs are submitted to a queue,
indicated by the `#SBATCH` directive:
```
#SBATCH --qos=<queue>
```
Collab has six queues:  open, regular, debug, express, and warp;
each serves a different purpose, and has different restrictions.

| queue (QOS) | description | restrictions |
| ---- | ---- | ---- |
| open | no-cost access | Portal and vintage hardware only |
| regular | for "normal" jobs | time < 14 days |
| debug	| for testing, debugging, <br> quick analysis | time < 1 hour |
| express | for rush jobs; <br> 1.5x price | time < 14 days |
| warp | for emergency jobs; <br> 3x price | time < 7 days |

Express and warp queues cost more,
but increase a job’s priority so it runs sooner.


## Converting from PBS

The Penn State cluster before Collab
used a different queue manager, PBS (Portable Batch System).
SLURM and PBS perform the same functions
(submit jobs, check status, cancel jobs, etc.)
but the commands are different.

| Action | PBS | SLURM |
| ---- | ---- | ---- |
| Request an interactive job | `qsub -I` | `salloc` |
| Submit a batch job | `qsub` | `sbatch` |
| Cancel a job | `qdel` | `scancel` |
| Check job status for user | `qstat -u <user>` | `squeue -u <user>` |
| Hold a job | `qhold` | `scontrol hold` |
| Release a job | `qrls` | `scontrol release` |

Resource requests in PBS batch scripts 
must be rewritten to run under SLURM.
Fortunately, the translation is straightforward.

| Resource Request | PBS | SLURM |
| ---- | ---- | ---- |
| Directive | `#PBS` | `#SBATCH` |
| Number of nodes | `-l nodes` | `-N` or `--nodes` |
| Number of CPUs | `-l ppn` | `-n` or `--ntasks` |
| Amount of memory | `-l mem` | `--mem` or `--mem-per-cpu` |
| Walltime | `-l walltime` | `-t` or `--time` |
| Allocation | `-A` | `-A` or `--account` |

For more information, see the 
[Slurm Rosetta Stone](https://slurm.schedmd.com/rosetta.pdf).


## [Timing jobs](#timing-jobs)

It is good practice to test a new workflow
by running small or short jobs before submitting big or long jobs.
To help plan your compute usage, 
it is helpful to time such test jobs.

Many well-designed applications display timing information
at the end of the log files they generate.
If this is not the case for your application,
you can find out how long a batch job takes
by sandwiching the commands you execute
between `date` commands:
```
date
<commands>
date
```
Your batch standard output file will then contain two "timestamps",
from which you can determine the running time.
