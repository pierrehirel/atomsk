#!/bin/bash

# This script installs Atomsk in system directories
# so that all users can use it. This script must be run
# with root permissions.

BINPATH=/usr/local/bin/
DOCPATH=/usr/local/share/doc/
MPATH=/usr/local/share/man/man1/

if [ -w ${BINPATH} ] ; then

  echo "<?> Atomsk will be installed in ${BINPATH}. Continue? (y/n)"
  read answer

  if [ "${answer}" = "y" ] ; then
    # System configuration file
    cp -rf ./etc/atomsk.conf /etc/

    # atomsk binary
    chmod +x atomsk
    cp atomsk ${BINPATH}

    # atomsk tools
    chmod +x ./tools/*.sh
    cp ./tools/* ${BINPATH}
    echo ">>> Program was successfuly installed in ${BINPATH}"

    # atomsk documentation
    mkdir -p ${DOCPATH}/atomsk
    rm -rf ${DOCPATH}/atomsk/*
    cp -rf ./doc/* ${DOCPATH}/atomsk/
    chmod -R a+r ${DOCPATH}/atomsk/
    echo ">>> Html documentation is available from ${DOCPATH}atomsk/index.html"

    # atomsk man page
    mkdir -p ${MPATH}
    gzip -c ./man/atomsk >${MPATH}/atomsk.1.gz
  
  
  else
  
    echo ">>> Installation aborted."
  
  fi

else

  echo "X!X ERROR: this script must be run as root or with sudo."

fi