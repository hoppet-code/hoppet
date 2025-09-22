#!/bin/zsh
# create a tar archive

# deduce version automatically from the appropriate include file
version=`grep HOPPET src/welcome_message.f90 | sed 's/.*HOPPET v. *//' | sed "s/ *' *//"`

origdir=`pwd | sed 's/.*\///'`
echo "Will make an archive of $origdir/"
dirhere=hoppet
dirtar=hoppet-$version
tarname=$dirtar.tgz
tmptarname=tmp-$tarname

# make sure we have Makefile with use CGAL=no

if [[ -e ../$tarname ]]
then
  echo "Tarfile $tarname already exists. Not proceeding."
elif [[ -e /tmp/$dirtar ]]
then
  echo "/tmp/$dirtar already exists, not proceeding"
else
  pushd ..

  echo "Creating tmp-$tarname"
  tar --exclude '.svn*' --exclude '*~' -zcf $tmptarname \
                      $dirhere/(src|example_f77|example_f90|benchmarking|benchmarking)/**/*.(f90|f|h|hh|alg|c|cc|C|eps|cpp|gp|pdf) \
                      $dirhere/doc/HOPPET-v1-doc.tex \
                      $dirhere/doc/*.(eps|sty) \
                      $dirhere/ChangeLog \
                      $dirhere/hoppet-config.in \
                      $dirhere/(src|example_f77|example_f90|benchmarking)/**/Makefile \
                      $dirhere/**/(README|INSTALL|Doxyfile|NEWS|COPYING|mkmk|configure|Makefile) \
                      $dirhere/scripts/*[a-z] \
                      $dirhere/example_f90/*.default_output 

  fulltarloc=`pwd`
  pushd /tmp
  echo "Unpacking it as /tmp/$dirhere"
  tar zxf $fulltarloc/$tmptarname
  mv -v /tmp/$dirhere /tmp/$dirtar
  echo "Repacking it with directory name $dirtar"
  tar zcvf $fulltarloc/$tarname $dirtar
  echo 
  echo "Removing /tmp/$dirhere"
  rm -rf $dirtar
  popd
  rm -v $tmptarname

  echo ""
    # if it's gavin running this then automatically copy the tarfile
    # to the web-space
  # webdir=~salam/www/repository/software/fastjet/
  # if [[ $USER = salam && -e $webdir ]]
  #     then
  #     echo "Copying .tgz file to web-site"
  #     cp -vp $tarname $webdir
  #     echo "************   Remember to edit web page **********"
  # fi

  popd

  # reminders about what to do for svn
  URL=`svn info | grep URL | sed 's/^.*URL: //'`
  tagURL=`echo $URL | sed "s/trunk/tags\/hoppet-$version/"`
  echo "Remember to tag the version:"
  echo "svn copy  -m 'tagged release of release $version' $URL $tagURL"
  echo "Copy it to hepforge:"
  echo "scp -p $fulltarloc/$tarname hepforge:hoppet/public_html/downloads/"
fi

#tar zcf $tarname


