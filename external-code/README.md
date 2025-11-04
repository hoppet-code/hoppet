# external-code/

This directory contains external code. For anything using git subtrees, you _must_
record the remote and the commit hash, tag or branch used.

## libome

### 2025-11-01: updated to 1.0.2 from the main repo

git subtree pull --prefix external-code/libome git@gitlab.com:libome/libome.git v1.0.2 --squash

### 2025-10-29:
Added libome, with the following command, run FROM THE TOP-LEVEL DIRECTORY

```
git subtree add --prefix=external-code/libome git@gitlab.com:hoppet-code/libome-fork.git libome-func_copy_and_log-quietnan.patch --squash -m "Added git@gitlab.com:hoppet-code/libome-fork.git, libome-func_copy_and_log-quietnan.patch"
```