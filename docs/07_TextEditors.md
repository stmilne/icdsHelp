# Text editors

Unix is a text-based operating system;
programs, batch scripts, and parameter files are text files.
To work on Unix, you need a good text editor.  There are several options.  

`gedit` is a modern windows-style text editor,
with cursor / mouse / cut-and-paste / menus,
reasonably intuitive for Windows and Mac users.
It is available from the Portal Interactive Desktop
(under Applications/Accessories/Text Editor),
or via `ssh -X` (SSH with X Windows enabled).

Other popular editors on Linux include [`emacs`][emacs]
(under Applications/Accessories/emacs, or via `ssh -X`)
and [`nano`][nano] (from the command line, via `ssh -X`).
[emacs]: https://www.gnu.org/software/emacs/
[nano]: https://www.nano-editor.org

For coding projects,
a popular choice is [Visual Studio Code Server][vscode]
(available as "VS Code Server" from the Portal main page).
Not just an editor, VSCode is an Integrated Development Environment,
that checks for errors as you type,
runs Python code interactively, has a debugger, and so on.
[vscode]: https://code.visualstudio.com/docs/remote/vscode-server

Finally, `vi` is beloved by old-school users.
In vi, you navigate the file not with cursor and mouse, but by key commands.
vi offers powerful search-and-replace functions,
and rapid navigation within the file.
Learning vi is like touch-typing -- 
difficult at first, but remarkably fast.  

Tutorials for vi are available [here][vi1], [here][vi2], and [here][vi3].
Once you know vi, [this summary][vi_summary] is useful.
[vi1]: https://www.cs.colostate.edu/helpdocs/vi.html 
[vi2]: http://heather.cs.ucdavis.edu/~matloff/UnixAndC/Editors/ViIntro.html
[vi3]: http://www.viemu.com/a_vi_vim_graphical_cheat_sheet_tutorial.html.  
[vi_summary]: https://bpb-us-e1.wpmucdn.com/sites.psu.edu/dist/0/79295/files/2020/09/notes-on-vi.pdf

Files can be edited on the laptop and transferred to the cluster.
If you do this, you must use a text editor, not a word processor,
which introduces hidden characters that make a mess.
On the Mac, use BBEdit, available [here](https://www.barebones.com/products/bbedit/);
on the PC, the best alternative is Notepad++, available [here](https://notepad-plus-plus.org).
