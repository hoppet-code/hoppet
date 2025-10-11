# This is the checklist for making releases

- [ ] run ctest
- [ ] check/update NEWS.md, including release date; commit
- [ ] run `scripts/set-version.sh X.Y.Z`
- [ ] update ChangeLog as per suggestion, adding "Release 2.0.0"
- [ ] commit and push (on dev branch)
- [ ] check that wheels work, by going to actions, "Build Python Wheels"
      and selecting the dev branch (takes c. 15 minutes)
- [ ] from master `git merge origin/dev; git push` OR create pull request
- [ ] once merged, make a release from https://github.com/hoppet-code/hoppet/releases/new
      tagging it as hoppet-X.Y.Z; as text, include the NEWS.md text
- [ ] go to on github actions and press "Build Python Wheels"
- [ ] run `scripts/set-version.sh X.Y.Z+1-dev` and update ChangeLog