# Tips for beginners

Newbies make mistakes.
Here are some to avoid:

!!! warning "Don’t use submit nodes for heavy computing."
     Submit nodes are for preparing files, submitting jobs, 
     examining results, and transferring files.

!!! warning "Don’t store files on scratch."
     [Scratch is not backed up](10_FileStorage.md), 
     and files more than 30 days old are deleted.

!!! warning "Don’t overrun your file storage quota."
     If you fill your allotted disk space, weird errors result.
     Keep an eye on your [disk space usage](10_FileStorage.md/#quotas).

!!! warning "Don’t waste your compute resources."
     Before a big job, run test jobs to make sure your code works.
     [Time the test job](08_BatchJobs.md/#timing-jobs), and plan ahead.