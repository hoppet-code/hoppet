#!/bin/bash
#
# script to generate the html for the documentation
#

inname=HOPPET-v1-doc
outname=${inname}l2h

# remove 12 pt: apparently causes cropping problems?
# lstlisting -> verbatim: (lose colours, but get usable code)
echo '\nonstopmode' > $outname.tex
sed -e 's/\[12pt\]//' -e 's/lstlisting/verbatim/g' < $inname.tex >> $outname.tex


# the following go into ~/.latexhtml-init
# # apparently needed for html 4
# $LOWER_CASE_TAGS=1;
# # make sure we have local copies of icons
# $LOCAL_ICONS = 1;
# # try to get rid of black lines
# $DVIPSOPT = ' -E';


# actually run things
#latex $outname.tex
latex -r $outname.tex
latex -r $outname.tex
latex2html -split '+0'\
           -show_section_numbers  -numbered_footnotes \
           $outname.tex
rm $outname.*
