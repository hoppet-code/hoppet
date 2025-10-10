# This is the checklist for making releases

- [ ] run ctest
- [ ] check/update NEWS.md, including release date; commit
- [ ] run `scripts/set-version.sh X.Y.Z`
- [ ] update ChangeLog as per suggestion, adding "Release 2.0.0"
- [ ] commit and push (on dev branch)
- [ ] from master `git merge origin/dev; git push` OR create pull request
- [ ] run `scripts/set-version.sh X.Y.Z+1-dev` and update ChangeLog