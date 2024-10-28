# System overview

Roar consists of two clusters:  Roar Collab, and Roar Restricted.
Roar Collab is for general use; 
Roar Restricted is only for users working with sensitive data
requiring extra security provisions.
For details, see the [Roar Restricted](15_RoarRestricted.md) addendum.

## Two ways to access 

Roar Collab can be accessed in two ways:  via the web-based Roar Collab Portal <br>
<https://rcportal.hpc.psu.edu/pun/sys/dashboard> <br>
and by "secure shell" (`ssh`) from a terminal application running on a laptop.

The Portal is designed mainly for interactive work.
It provides:

- graphical, number-crunching programs, 
such as ANSYS, COMSOL, MATLAB, and RStudio;
- a Windows-like desktop environment;
- a web-based file browser, to upload and download files.

The Portal is easy to use, 
because its programs can be launched and used
without knowing Unix.

But for "batch jobs" -- compute jobs that take a long time,
and run without user intervention --
[access via `ssh`](05_AccessViaSSH.md) is needed to prepare and submit jobs.

## Accounts and allocations

To logon to Roar, you need an [account](03_Accounts.md).
To do work on Roar, you can use the open queue at no cost,
which gives you access to the Portal, 
and to batch jobs on vintage hardware.
But to use any of the newer, more powerful hardware,
you need a paid [allocation](04_Allocations.md).

## Roar uses Unix

The operating system for Roar is Unix.
Unix is text-based; users interact with the system by typing commands.
Compute clusters use Unix,
in part because tasks can be automated with scripts.

This User Guide assumes familiarity with Unix,
which any user who wants to do more than use the Portal needs to learn.
A tutorial website for Unix can be found [here][unix_tutorial].
[unix_tutorial]: https://www.tutorialspoint.com/unix/unix_tutorial.pdf
An introduction to Unix, 
written for research undergraduates and graduate students,
is [here](pdf/unixGuide.pdf).

## Tips for beginners

Newbies make mistakes.
Here are some to avoid:

- _Don’t use submit nodes for heavy computing_.
Submit nodes are for preparing files, submitting jobs, 
examining results, transferring files, and so on.
- _Don’t store files on scratch_.
[Scratch is not backed up](10_FileStorage.md), 
and files more than 30 days old are deleted.
- _Don’t overrun your file storage quota_.
Every account comes with some disk space;
if you fill it, weird errors result.
Keep an eye on your [disk space usage](10_FileStorage.md/#quotas).
- _Don’t waste your compute resources_.
Before a big job, run small test jobs
to make sure everything works.
[Time your test jobs](08_BatchJobs.md/#timing-jobs), and plan ahead.
