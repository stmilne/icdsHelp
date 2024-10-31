# Batch examples

Here we present examples of `#SBATCH` directives
for common use cases.  
(Information about queues is [here](08_BatchJobs.md/#queues),
and about hardware partitions is [here](03_ComputeHardware.md/#partitions).)

## Open queue

When you run with `account=open`, 
only `partition=vintage` and `qos=open` are possible.

```
#SBATCH --account=open
#SBATCH --partition=vintage
#SBATCH --qos=open
#SBATCH --nodes=1
#SBATCH --tasks=16
```

## Basic nodes

If your CPU job fits within a single 64-core node, 
it is most economical to run on basic nodes.

```
#SBATCH --account=<myAllocation>
#SBATCH --partition=basic
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --tasks=16
```

## Multi-node job

If your CPU job needs multiple nodes
with fast (Infiniband) interconnects,
use standard nodes.

```
#SBATCH --account=<myAllocation>
#SBATCH --partition=standard
#SBATCH --qos=regular
#SBATCH --nodes=4
#SBATCH --tasks=48
```

## GPU job
 
Each node in the GPU partition has two A100 GPUs;
the price is per GPU requested.

```
#SBATCH --account=<myAllocation>
#SBATCH --partition=GPU
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --tasks=8
#SBATCH --gpus=1
```

## Older GPUs

Nodes in the GPU2 partition have one P100 GPU,
which are slower but less expensive than A100s.

```
#SBATCH --account=<myAllocation>
#SBATCH --partition=GPU2
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --tasks=16
#SBATCH --gpus=1
```

## SLURM directives

Most of the usual SLURM directives appear in the examples above.
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

Most `#SBATCH` directives have *defaults*, listed above,
so you don't need to specify everything.

`#SBATCH` requests must be *consistent*, or the submission will fail.
Don't ask for a GPU if you requested `partition=basic`;
don't specify `account=open`, and then ask for anything but
`partition=vintage`; and so on.

For more information, see documentation on 
[salloc](https://slurm.schedmd.com/salloc.html) 
and [sbatch](https://slurm.schedmd.com/sbatch.html).

## PBS to SLURM

Resource requests in old PBS batch scripts 
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


