# Access via ssh

To use the cluster beyond the several applications provided by the Portal,
logging on via a terminal is necessary.
The Portal does provide a desktop environment,
which includes a terminal application;
however, the app takes several minutes to load,
requires an estimate of how long you will be logged on,
and ties up multiple cores on the interactive nodes.

An alternative is to logon to the cluster
from a terminal application running on your laptop.
On Windows, use puTTY, available [here](https://www.putty.org);
on Mac, Terminal (installed by default), 
or iTerm  (an improved version), available [here](https://iterm2.com).

From the terminal, logon as
```
ssh <user>@submit.hpc.psu.edu
```

You will be prompted for your password, 
and then for multi-factor authentication (MFA), 
which confirms that you are you.  To set up MFA, 
go [here](https://accounts.psu.edu/2fa).

**Note:  for security reasons, Roar Restricted can only be accessed
via the [RR Portal](15_RoarRestricted.md/#Roar Restricted/).**

## X Window apps

To use any application that "opens a window"
(called an  "X Window" or "X11" application), 
you need an additional program on your laptop.
On the Mac, this is XQuartz, available [here](https://www.xquartz.org).
On the PC, you need VcXsrv, available [here](https://sourceforge.net/projects/vcxsrv/).

To use cluster applications that open windows, logon with
option `-X` for "X forwarding":
```
ssh -X <user>@submit.hpc.psu.edu
```
When you logon to Collab from somewhere off campus,
X Window apps can sometimes be slow;  
Portal works better in such circumstances.

## User settings

On Macs and PCs, user settings are accessed in application menus,
System Settings panels, and the like.
On Unix machines, settings are stored in text files, 
typically hidden (named starting with `.`),
somewhere in your home directory (`~` or `$HOME`).

Settings files are hidden to keep users from meddling with them.
But one important settings file that you should know about is `~/.bashrc`.

When you log in, commands in `.bashrc` are automatically executed.
Lines may be added to `.bashrc` to define frequently used aliases, 
or load [modules](12_LoadingSoftware.md) so that software you use often
is available automatically.
