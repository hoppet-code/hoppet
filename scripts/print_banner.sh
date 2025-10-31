#!/bin/bash
# Print colored banner for hoppet 2.0.0 release

# Some colours for printout
RED='\033[0;31m'
GREEN='\033[0;32m'
PURPLE='\033[1;35m'
NC='\033[0m' # No Color
version=$1 # From CMakeLists.txt


# define a function to print centred text
center_text() {
    local width=$1
    while IFS= read -r line; do
        local padding=$(( (width - ${#line}) / 2 ))
        local postpadding=$(( width - ${#line} - padding ))
        printf "%*s%s%*s\n" $padding "" "$line" $postpadding ""
    done
}
centredRelease="$(center_text 78 <<< "HOPPET ${version} RELEASE")"

echo -e "\n"
echo -e "################################################################################"
#echo -e "#${PURPLE}                             HOPPET ${version} RELEASE                             ${NC}#"
echo -e "#${PURPLE}${centredRelease}${NC}#"
echo -e "#                                                                              #"
echo -e "# Note that since the 2.0.0 release, hoppet’s library name has been renamed    #"
echo -e "# from ${RED}hoppet_v1${NC} to ${GREEN}hoppet${NC}, and similarly for the main module and C++ include  #"
echo -e "# file. Users with an existing installation of hoppet (v1) should make sure    #"
echo -e "# that they link with the new library name.                                    #"
echo -e "#                                                                              #"
echo -e "################################################################################"
