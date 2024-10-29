# PBS to SLURM

The Penn State cluster before Collab
used a different queue manager, PBS (Portable Batch System).
SLURM and PBS perform the same functions
(submit jobs, check status, cancel jobs, etc.),
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


