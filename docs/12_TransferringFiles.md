# Transferring files

Computational workflows often require files 
on Roar Collab, Roar archival storage, 
OneDrive, or your laptop or desktop machine
to be transferred -- copied from one place to another.

Transfers to and from Roar Restricted follow a special protocol,
described [here](16_RoarRestricted.md).

Multiple tools exist to perform these file transfers.
No single tool is best for all cases;
below, we recommend methods, 
and list approximate transfer rates for large files.

| Transfer | Method | Rate (MB/sec) |
| ---- | ---- | ---- |
| Collab &harr; Archive | Globus | 50 |
| Collab &rarr; OneDrive | Firefox or Globus | 50 | 
| OneDrive &rarr; Collab | Firefox or Globus | 10 |
| Collab &harr; laptop | Portal Files menu | 25 |
| Collab &harr; laptop | Cyberduck or FileZilla | 15 |
| OneDrive &harr; laptop | web access |20 |

## Firefox

The Firefox browser is available either 
from the [Portal Interactive Desktop][portalID]
or from an [SSH -X session][sshx].
From Firefox, you can access OneDrive
and other similar destinations,
and perform file transfers.
[portalID]: 06_AccessingCollab.md
[sshx]: 06_AccessingCollab.md#access-via-ssh

## Portal Files menu

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
from the local machine to the remote with `put <filename>`,
and from the remote machine to the local with `get <filename>`.

When sftp launches, your location on the local machine
is the folder in which you launched sftp,
and your location on the remote machine is your home directory.
To control where on the remote machine files go to and come from,
within sftp you can navigate on the remote machine with `cd`
and list files with `ls`. 
Likewise, within sftp you can navigate on the local machine
with "local" versions of these commands, `lcd` and `lls`.

!!! tip ""
     "Graphical" sftp clients for your laptop
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

Globus moves files between named "endpoints":

| Filesystem | Endpoint |
| ---- | ---- |
| Collab | PennState_ICDS_RC |
| Archive | PennState_ICDS_Archive | 
| PSU OneDrive | "Penn State OACIOR OneDrive Collection 01" |

Globus is interactive, 
but time-consuming file transfers 
can be submitted as batch jobs.

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
With option `-v` (as in `tar -cvf`), tar operates "verbosely", 
listing every file added.

To make a tarball of everything in the current directory:
```
tar -cf myTar.tar *
```

tar with option `-z` (as in `tar -czf`)
makes a "zipped" (compressed) tarball.
Compressing text files often makes the tarball much smaller;
compressing binary or image files is usually pointless.

To unpack a tarball (option `-x` for "extract"):
```
tar -xf myTar.gz 
```

## rsync

Sometimes, you want to copy a directory of files 
from one place to another,
and then later *update* the copy
with any changes made to the originals,
so that the copy reflects the current version.
If the directory contains many files but only a few are changed,
it would be nice to have a program that automatically updates
only the changed files. [`rsync`][rsync] does this:
[rsync]: https://linux.die.net/man/1/rsync
```
rsync <options> <source-path> <destination-path>
```
copies files from `<source-path>` to `<destination-path>`,
*and deletes files if necessary*,
to make the destination the same as the source.
Later, if you change files on `<source-path>`,
run rsync again to update files on `<destination-path>`.

rsync has several important options:

- `a` – "archive" mode; traverses directories recursively
- `v` – "verbose"; reports which files are copied
- `z` – "zips" (compresses) the files on transfer
- `h` – "human readable" reporting

The source and destination can be on the same filesystem,
or they can be different machines entirely.
From a Unix command line on your laptop 
(running Linux or OS X),
```
rsync /work/newData abc123@submit.hpc.psu.edu:/storage/work/abc123/toAnalyze/
```
which will prompt for your password and MFA.

Note that the *destination* pathnames end with a `/`;
this signifies that the copied files will go into the folder `toAnalyze`.
For more examples, visit [Tecmint][tecmint].
[tecmint]: https://www.tecmint.com/rsync-local-remote-file-synchronization-commands/

!!! warning ""
     With rsync, the source is the original, and the destination is the copy.
     Don't reverse direction, or you will confuse rsync and yourself,
     and wind up clobbering or deleting files.
