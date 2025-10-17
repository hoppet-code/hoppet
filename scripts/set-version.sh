#!/bin/bash
#
# This script sets all necessary version-number information across
# Hoppet (including manual). It should be run from the
# top-level Hoppet directory.
#

if [ $# -ne 1 ]
then
 echo "Usage: scripts/set-version.sh version-number"
 exit
fi

version=$1
versionNumber=`echo $version | sed 's/-.*//'`
# if there is a dash in the version number, take everything afterwards
# for versionPreRelease, otherwise set versionPreRelease to an empty string
if echo $version | grep -q '-'
then
  versionPreRelease=$(echo $version | sed 's/.*-/-/')
else
  versionPreRelease=""
fi

changedFiles=""
checkChanged() {
  changedFile=$1
  if [ -n "$changedFile" ]; then
    changedFiles="$changedFiles $changedFile"
  fi
  diff $changedFile.bak $changedFile | colordiff
}

bold=$(tput bold)
normal=$(tput sgr0)

echo ${bold}"------------ Will set Hoppet version to $version -----------"${normal}
echo "versionNumber = $versionNumber"
echo "versionPreRelease = $versionPreRelease"

echo
echo ${bold}"------------ Setting it in CMakeLists.txt -------------------"${normal}
sed -i.bak 's/^project.hoppet VERSION [0-9.]*/project(hoppet VERSION '$versionNumber'/' CMakeLists.txt
sed -i.bak 's/^set.PROJECT_VERSION_PRERELEASE.*/set(PROJECT_VERSION_PRERELEASE "'$versionPreRelease'")/' CMakeLists.txt
checkChanged CMakeLists.txt
#diff CMakeLists.txt.bak CMakeLists.txt

echo ${bold}"------------ Setting it in src/welcome_message.f90 ---------------------"${normal}
sed -i.bak 's/private_version_string = .*/private_version_string = "'$version'"/' src/welcome_message.f90
checkChanged src/welcome_message.f90
#diff src/welcome_message.f90.bak src/welcome_message.f90
#AC_INIT([FastJet],[3.0.2-devel])

echo ${bold}"------------ Setting it in pyproject.toml ---------------------"${normal}
sed -i.bak 's/version = .*/version = "'$version'"/' pyproject.toml
checkChanged pyproject.toml

echo ${bold}"------------ Setting it in pyinterface/hoppet.i ---------------------"${normal}
sed -i.bak 's/__version__ = .*/__version__ = "'$version'"/' pyinterface/hoppet.i
checkChanged pyinterface/hoppet.i


# echo
# echo "------------ Setting it in configure.ac ---------------------"
# #sed -i.bak 's/\(AC_INIT.*\)])/\1-'$extralabel'])/' configure.ac
# sed -i.bak 's/^\(AC_INIT(\[.*\],\[\).*/\1'$version'])/' configure.ac
# diff configure.ac.bak configure.ac 
# #AC_INIT([FastJet],[3.0.2-devel])

# # now make sure the windows config file is consistent
# echo
# echo "------------ Setting it in include/fastjet/config_win.h -----"
# cp -p include/fastjet/config_win.h include/fastjet/config_win.h.bak
# cd src
# ./genconfig.sh ../include/fastjet/config_win.h
# cd ..
# diff include/fastjet/config_win.h.bak include/fastjet/config_win.h 

# echo
# echo "------------ Setting it in Doxyfile -------------------------"
# sed -i.bak 's/^\(PROJECT_NUMBER.*=\).*/\1 '$version'/' Doxyfile
# diff Doxyfile.bak Doxyfile


echo
echo ${bold}"------------ Setting it in doc/HOPPET-doc.tex --------------"${normal}
sed -i.bak 's/^\( *\)[^%]*\(%.*VERSION-NUMBER.*\)/\1version '$version'\2/' doc/HOPPET-doc.tex
checkChanged doc/HOPPET-doc.tex

echo
echo ${bold}"------------ Recommended ChangeLog entry --------------------"${normal}
# NB: -e option of echo ensures that \t translates to a tab character
echo -e "\t* CMakeLists.txt:"
# loop over changed files to print out which ones were changed
for f in $changedFiles; do
  if [ "$f" != "CMakeLists.txt" ]; then
    echo -e "\t* $f:"
  fi
done
echo -e "\tchanged version to $version"
