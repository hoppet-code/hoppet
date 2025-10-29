# external-code/

This directory contains external code. For anything using git subtrees, you _must_
record the remote and the commit hash, tag or branch used.

## libome

Added 2025-10-29, with the following command, run FROM THE TOP-LEVEL DIRECTORY

```
git subtree add --prefix=external-code/libome git@gitlab.com:hoppet-code/libome-fork.git libome-func_copy_and_log-quietnan.patch --squash -m "Added git@gitlab.com:hoppet-code/libome-fork.git, libome-func_copy_and_log-quietnan.patch"
```