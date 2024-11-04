# File storage

Files are stored on Roar in a central filesystem;
Collab and Restricted filesystems are separate.

Users have access to four directories:  home, work, group, and scratch.
These have different purposes:

- **home** – for configuration files, and links to work, group, and scratch.
- **work** – for your own work; 
only you have read-write access to your home and work directories.
- **group** – for collaborative work.  group space is owned by a PI;
by default, group members have read access to all files in group.
Space in group is not free, and is sold in 5 TB increments.
- **scratch**  – for temporary storage of large files.  scratch is *not backed up*, 
and files older than 30 days are *automatically deleted*.

Files in home, work, and group are backed up by a sequence of daily "snapshots". 
Snapshots are not saved beyond xxx days.
If you mistakenly delete a file, you can request it to be restored from a snapshot.
Try to avoid this.

## Archive storage

To store infrequently-used files, low-cost archive storage can be purchased. 
The Collab Archive is only accessible via the [Globus][globus] web interface,
so access is not quick or convenient.
If you store a directory that contains many files, 
you should pack the directory into a single file with `tar`
(see [Packing files](12_TransferringFiles.md#packing-files))
before transferring.
[globus]: 12_TransferringFiles.md#globus

## Quotas

home, work, group, and scratch directories are subject to limits,
both on the total filespace and on the total number of files.
On Collab, the limits are:

| Storage | Path | Space | Files | Backup | Purpose |
| :----: | :----: | :----: | :----: | :----: | :----: |
| Home | /storage/home | 16 GB | 500,000 | Daily  | Configuration files |
| Work | /storage/work | 128 GB | 1 million | Daily  | User data |
| Scratch | /scratch | None | 1 million | None | Temporary files |
| Group | /storage/group | Specific to<br>allocation | 1 million <br>per TB | Daily | Shared data |

On Roar Restricted, there is no scratch or group;
otherwise the limits are the same.


## Checking storage

If you fill or nearly fill your home or work directories,
weird errors will result when you try to run programs or write files.
To avoid this, keep an eye on your file sizes and total usage.

There are two tools to check on your disk usage:

- `check_storage_quotas` reports your total usage;
- [`du`][du] reports the sizes of files and directories.
[du]: https://man7.org/linux/man-pages/man1/du.1.html

`du -sh *` gives "human-readable" sizes (in MB, GB, TB) 
for each item in the current directory.

In managing storage, we often want to know where the big files are.
``
du -sh * | sort -h -r
``
lists directory sizes in order from large to small
(the output of du is "piped" to [sort][sort]).
[sort]: https://man7.org/linux/man-pages/man1/sort.1.html

## File permissions

If you are a PI, you have group space on Collab, located at
`/storage/group/<PIuserID>/default`,
for which the default file permissions and ownership are
```
drwxr-s-- 2 root <PIuserID>
```

The `s` setting in the group permissions means 
every file and folder created within the group folder
will have the same group read `r` permission.
However, files created elsewhere and moved into group 
have the permissions they were created with.
To change them, use [`chmod`][chmod]:
[chmod]: https://man7.org/linux/man-pages/man1/chmod.1p.html
```
chmod g+r <filespec>
```
More generally, chmod can be used to add or remove (+ or -) 
permissions to read (r), write (w), or execute (x),
for the user (u), group (g), or others (o).



