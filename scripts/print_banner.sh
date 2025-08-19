#!/bin/bash
# Print colored banner for hoppet 2.0.0 release

# Some colours for printout
RED='\033[0;31m'
GREEN='\033[0;32m'
PURPLE='\033[1;35m'
NC='\033[0m' # No Color
version=$1 # From CMakeLists.txt

echo -e "\n######################################################################################"
echo -e "#${PURPLE}                                HOPPET ${version} RELEASE                                ${NC}#"
echo -e "#                                                                                    #"
echo -e "# Note that since the 2.0.0 release, hoppet’s library name has been renamed from     #"
echo -e "# ${RED}hoppet_v1${NC} to ${GREEN}hoppet${NC}, and similarly for the main module and C++ include file. Users #"
echo -e "# with an existing installation of hoppet (v1) should make sure that they link with  #"
echo -e "# the new library name.                                                              #"
echo -e "#                                                                                    #"
echo -e "######################################################################################"
