# Batch jobs

For compute jobs that take hours or days to run,
instead of sitting at the terminal waiting for the results,
we submit a "batch job" to the queue manager,
which runs the job when resources are available.

On Roar, the queue manager is called SLURM 
(Simple Linux Utility for Resource Management).  
Besides `salloc` for interactive jobs (see [Interactive jobs](07_InteractiveJobs.md)),
the basic SLURM commands are:

| Command | Effect|
| ---- | ---- |
|`sbatch <script>` | submit `<script>` to SLURM |
| `squeue -u <userid>` | check on jobs submitted by `<userid>` |
| `scancel <jobID>` | cancel the job |

## Batch scripts

Jobs submitted to SLURM are in the form of a "batch script".
A batch script is a shell script that executes commands,
with a preamble of comments of the form `#SBATCH ...` specify:

- an *allocation* to charge for the job;
- a *queue* (qos) to submit the job to;
- a *partition* (type of nodes) to run on;
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

module load gromacs-2019.6

cd $SLURM_SUBMIT_DIR
SYSTEM=$SLURM_SUBMIT_DIR/System

gmx grompp -f $SYSTEM/nvt.mdp -c $SYSTEM/min.gro -p $SYSTEM/testJob.top -o nvt.tpr 
gmx mdrun -nt 8 -nb cpu -deffnm nvt
```

Everything after the last `#SBATCH` are commands to be executed.
Most scripts start with `cd $SLURM_SUBMIT_DIR`,
which is the directory from which the job was submitted.

The first line `#!/bin/bash` executes your `.bashrc` file,
and thus performs any initializations it contains
(e.g., loading modules).

**The hardware you reserve by choosing a partition,
specifying numbers of nodes and cores, or requesting a GPU,
affects the cost of your job; see [Prices](04_Allocations.md/#prices).**


## SLURM directives

Most of the usual SLURM directives appear in the example above.
Here is a table:

| Directive | Default | Description |
| ---- | ---- | ---- |
| `-A` or `--account` | `open` | Allocation to charge |
| `-q` or `--qos` | `open` | Queue to submit to |
| `-p` or `--partition` | `vintage` | Partition (type of nodes) to use |
| `-N` or `--nodes` | 1 | Number of nodes |
| `-G` or `--gpus` | 0 | Number of GPUs |
| `-n` or `--ntasks` | 1 | Number of tasks (cores) |
| `--mem` | ? | Memory per node |
| `-t` or `--time` | 1 hour | Maximum run time |
| `-C` or `--constraint` | none | Required node features |
| `-J` or `--job-name` | `jobID` | Job name |

`#SBATCH` directives have *defaults*, listed above,
so you don't need to specify everything.

`#SBATCH` requests must be *consistent*, or the submission will fail.
Don't ask for a GPU if you requested `partition=basic`;
don't specify `account=open`, and then ask for anything but
`partition=vintage`; and so on.

For more information, see documentation on 
[salloc](https://slurm.schedmd.com/salloc.html) 
and [sbatch](https://slurm.schedmd.com/sbatch.html).

## Queues[](#queues)

The directive `#SBATCH --qos=<queue>` submits batch jobs to a queue, 
or QoS = "Quality of Service" in SLURM-speak.
(Queues are a bit like classes of service on an airline flight:
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

## Timing jobs[](#timing-jobs)

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
