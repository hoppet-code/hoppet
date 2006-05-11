#!/bin/zsh
# create a tar archive

#version=0.9.0c-20050929-1200
version=1.0
origdir=`pwd | sed 's/.*\///'`
echo "Will make an archive of $origdir/"
dir=$origdir-$version
tarname=$dir.tgz

pushd ..

if [[ -e $tarname ]]
then
  echo "Tarfile $tarname already exists. Not proceeding."
else
  echo "Creating $tarname:"
  if [[ -e $dir ]] 
  then
    echo "Could not create $dir as link to $origdir (former exists already)"
  else
    ln -s $origdir $dir
    tar zcvhf $tarname $dir/**/*.(f90|f|h|alg|sh|c|cc) \
                      $dir/(src|example_f77|testing)/(Makefile|*.pl) \
                      $dir/**/READM*[A-Z] $dir/ChangeLog
    rm $dir
  fi
fi

#tar zcf $tarname
popd


