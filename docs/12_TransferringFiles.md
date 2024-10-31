# Transferring files

Computational workflows often require files 
on Roar Collab, OneDrive, or your laptop or desktop machine
to be transferred -- copied from one place to another.

Transfers to and from Roar Restricted follow a special protocol,
described [here](16_RoarRestricted.md).

Multiple tools exist to perform these file transfers.
No single tool is best for all cases;
below, we recommend methods, 
and list approximate transfer rates for large files.

| Transfer | Method | Rate (MB/sec) |
| ---- | ---- | ---- |
| Collab &rarr OneDrive | Globus | 50 | 
| OneDrive &rarr Collab | Globus | 10 |
| Collab &harr laptop | Portal | 25 |
| Collab &harr laptop |sftp | 15 |
| OneDrive &harr laptop | web access |20 |

## Portal

The Collab [Portal][portal] top menu under Files/Home
opens a window that enables convenient file transfer 
between Collab and your laptop.
[portal]: https://rcportal.hpc.psu.edu/pun/sys/dashboard

## sftp

[`sftp`][sftp] (secure file transfer protocol) is a Unix tool
for file transfers.  To launch sftp,
[sftp]: https://man7.org/linux/man-pages/man1/sftp.1.html
```
sftp <username>@<address>
```
where `<address>` is the address of a "remote machine",
and `<username>` is your userid on that machine.
For sftp *to* Collab, the address is `submit.hpc.psu.edu`,
the same as for ssh logon.

Just as for ssh logon, you will be prompted 
for your password on the remote machine,
and multi-factor authentication (MFA)
if the remote machine requires it.

sftp is an interactive program.
Once logged on, you can copy files 
from the local machine to the remote, with
``` 
put <filename> 
```
and from the remote machine to the local, with 
```
get <filename> 
```
where `<filename>` is the name of the file to be copied.

When sftp launches, your location on the local machine
is the folder in which you launched sftp,
and your location on the remote machine is your home directory.
To control where on the remote machine files go to and come from,
you can navigate on the remote machine in sftp with `cd`
and list files with `ls`. 
Likewise, within sftp you can navigate on the local machine
with "local" versions of these commands, `lcd` and `lls`.

!!! tip ""
     "Visual" sftp clients for your laptop
     can be used for file transfer to Collab, 
     as well as to OneDrive or other cloud storage providers.
     Two popular options for both OS X and Windows are
     [Cyberduck][cyberduck] and [FileZilla][filezilla].
[cyberduck]:https://cyberduck.io
[filezilla]:https://filezilla-project.org

## Globus

Globus is a web-based tool for file transfers,
which can move files from Roar Collab
to filesystems at institutions outside Penn State,
including our OneDrive accounts.
To get started, go [here][globus].
[globus]: https://docs.globus.org/how-to/get-started/

Globus moves files between named "endpoints".
The endpoint for Collab is "PennState_ICDS_RC",
and for Penn State OneDrive is "Penn State OACIOR OneDrive Collection 01".

Globus is interactive, 
but has the advantage that time-consuming file transfers 
can be submitted as batch jobs, to be performed later.

## rsync

[`rsync`][rsync] is ...
[rsync]: https://linux.die.net/man/1/rsync

## Packing files

To transfer many files,
or an entire directory of files, possibly with subdirectories,
it is convenient to first make a single large file 
of the entire contents,
using the Unix command `tar`
("Tape ARchive", because in ancient days big files 
would be archived to magnetic tape).  The command
```
tar -cf <tarfilename>.tar <filespec>
```
creates (`-c`) a single .tar file (`-f`) (sometimes called a "tarball"),
containing the entire contents of `<filespec>`.
`tar -cvf` operates "verbosely", 
listing every file added to the tarball.

To make a tarball of everything in the current directory:
```
tar -cf myTar.tar *
```
To make a tarball of everything in `thisFolder`:
```
tar -cf myTar.tar thisFolder/*
```

tar with a `z` option
makes a "zipped" (compressed) tarball.
Compressing text files often makes the resulting tarball much smaller;
compressing binary or image files is usually pointless.

To unpack a tarball (option `-x` for "extract"):
```
tar -xf myTar.gz 
```

