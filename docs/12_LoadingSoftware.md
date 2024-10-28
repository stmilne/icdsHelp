# Loading software

Collab provides many different applications in its software stack,
which consists of two parts:

- preloaded software, available immediately on login;
- modular software, which must be loaded to be used.

To find out if an application is available, use `which`:
```
which <appName>
```
`which` returns the path to the executable `<appName>`,
if it exists on your `$PATH`.
(The `$PATH` variable contains the set of directories
where Unix looks for executables;
defaults for `$PATH` are initialized at logon.)

If `which` doesn’t find the app you are looking for,
it may need to be loaded via the `module` system.
Here are the various commands to search for,
get info on, list, load, and unload modules:

| Command | Description |
| ---- | ---- |
| `module avail` | List all modules |
| `module show <module_name>` | Show module file contents |
| `module spider <module_name>` | Search for a module |
| `module load <module_name>` | Load a module |
| `module load <module>/<version>` | Load a specific version |
| `module unload <module_name>` | Unload a module |
| `module list` | List all loaded modules |
| `module purge` | Unload all modules |
| `module use <path>` | Add a path to `$MODULEPATH` |


What module files mainly do is add various directories 
to the paths Unix searches when looking for applications.
To see what a given modulefile does, use `module show <moduleName>`.

## Python packages

Python is a widely-used programming language,
preloaded on Collab,
which is greatly extensible in its capabilities
by loading additional packages.

To load Python packages, use pip:
```
pip install --user <package>
```
Python packages are loaded by default into `$HOME/.local`.
If you load a large number of big Python packages,
you may exceed the limited storage space in your home directory.

If this becomes a problem, the `.local` directory can be moved to `~/work`,
and a link placed in your home directory.
To make such a link, in your home directory execute the command
```
ln -s .local work/.local
```
which creates an alias (in Unix-speak, a "soft link") named `.local`
that points to the directory you moved to `work`.
