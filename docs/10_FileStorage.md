# File storage

Files are stored on Roar in a central filesystem;
Collab and Restricted filesystems are separate.

Users have access to four directories:  home, work, group, and scratch.
These have different purposes:

- home is for configuration files, and links to work, group, and scratch.
- work is for your own work; 
only you have read-write access to your home and work directories.
- group is for collaborative work.  group space is owned by a PI;
by default, group members have read access to all files in group.
Space in group is not free, and is sold in 5 TB increments.
- scratch is for temporary storage of large files.  scratch is *not backed up*, 
and files older than 30 days are *automatically deleted*.

Files in home, work, and group are backed up by a sequence of daily "snapshots". 
Snapshots are not saved beyond xxx days.
If you mistakenly delete a file, you can request it to be restored from a snapshot.
Try to avoid this.

## [Quotas](#quotas)

home, work, group, and scratch directories are subject to limits,
of the total filespace and the total number of files.
On Roar Collab, the limits are:

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
- `du` reports the sizes of files and directories.

`du -sh *` gives "human-readable" sizes (in MB, GB, TB) 
for each item in the current directory.

In managing storage, we often want to know where the big files are.
``
du -sh * | sort -h -r
``
lists directory sizes in order from large to small
(the output of `du` is "piped" to `sort`).

## File permissions

If you are a PI, you have group space on Collab, located at
`/storage/group/<PIuserID>/default`
for which the default file permissions and ownership are
```
drwxr-s-- 2 root <PIuserID>
```

The `s` setting in the group permissions means 
every file and folder created within group
will by default have the same group read `r` permission.
However, files created elsewhere and moved into group 
have the permissions they were created with.
To change them, use `chmod`:
```
chmod g+r <filespec>
```
More generally, `chmod` can be used to add or remove (+ or -) 
permissions to read (r), write (w), or execute (x),
for the user (u), group (g), or others (o).


## Finding files

On Macs and PCs, powerful tools exist for finding files
that match certain criteria, or contain certain text.
Two Unix commands that provide similar capabilities 
are `find` and `grep`.

`find` looks for files based on their attributes (name, type, size, date),
searches recursively (through all subdirectories of a given directory),
and returns a list of file names.

`grep` looks for files based on their contents,
searching locally by default (but can be made to search recursively),
and returns a list of lines that contain the searched-for pattern.

Like many Unix commands, 
`find` and `grep` have many option that expand their capabilities,
and can be combined with other Unix commands 
using pipes (`|`) and redirection (`>`) to perform complex tasks.

Tutorials exist for [Unix][unix_tutorial], 
[grep][grep_tutorial], and [find][find_tutorial].
Google searches on these commands are often helpful
(e.g., "Unix how to find and delete large files").
Here are a few examples of each.
[grep_tutorial]: https://danielmiessler.com/p/grep/
[find_tutorial]: https://snapshooter.com/learn/linux/find

### find

Find all files with name that matches a pattern:
```
find . -name "*.trr"
```

Find files above a certain size:
```
find . -type f -size +100G
```

Find has a powerful option `-exec`,
which executes a command on files that are found.  
Use this with caution!  List the files first.

Find files, and then compress them:
```
find . -name "*.trr" -exec compress {} ; 
```

Find files of a certain type, get their sizes, and sum them:
```
find . -name "*.trr" -exec du -csh '{}' +
```

Find files of a certain type, get their sizes, and sort by size:
```
find . -name "*.trr" -exec du -sh '{}' + | sort -h -r
```

Find files above a certain size and list them: 
```
find . -type f -size +100G -exec ls -lh {} +
```

Find files of a certain extension and delete them: 
```
find . -name "*.trr" -exec rm -f {} +
```

### grep

In its simplest form, `grep` searches for a "pattern"
in a given file or set of files:
```
grep <pattern> <file>
```

In this form, `grep` finds all lines in `<file>` that contain `<pattern>`.
The pattern can have "wildcards" -- `*` matches anything, for example.

`<file>` can be a single filename, or itself a pattern
that matches multiple files; for example, `*.log` matches all logfiles.

`grep` is often useful to "filter" the output of commands.
`ps -e` lists all running processes (a very long list).
To look for any process that mentions `vmd`:
```
ps -e | grep vmd
```
`grep -R` searches recursively through subdirectories.
Here, `grep` searches the current directory `.` and its subdirectories,
including all files `*.sh`, looking for pattern `sed`, 
and writes the results to `sedExamples`: 
```
grep -R ‘.’ --include *.sh -e ‘sed’ > sedExamples
```
