#!/bin/bash

# This script installs Atomsk in system directories
# so that all users can use it. This script must be run
# with root permissions.

BINPATH=/usr/local/bin/
DOCPATH=/usr/local/share/doc/
MPATH=/usr/local/share/man/man1/

if [ ! -e 'atomsk' ] ; then

  echo "X!X ERROR: the program 'atomsk' does not exist in current directory."

else

  if [ -w ${BINPATH} ] ; then

    echo "<?> Atomsk will be installed in ${BINPATH}. Continue? (y/n)"
    read answer

    if [ "${answer}" = "y" ] ; then
      # System configuration file
      cp -rf ./etc/atomsk.conf /etc/

      # Atomsk binary
      chmod +x atomsk
      cp atomsk ${BINPATH}

      # Atomsk tools
      chmod +x ./tools/*.sh
      cp ./tools/* ${BINPATH}
      echo ">>> The program was successfuly installed in ${BINPATH}"
      echo "    To run it, enter 'atomsk' in a terminal."

      # Atomsk documentation
      mkdir -p ${DOCPATH}/atomsk
      rm -rf ${DOCPATH}/atomsk/*
      cp -rf ./doc/* ${DOCPATH}/atomsk/
      chmod -R a+r ${DOCPATH}/atomsk/
      echo ">>> The html documentation was installed. You may read it by entering the"
      echo "    following address in your Web browser: ${DOCPATH}atomsk/index.html"

      # Atomsk man page
      mkdir -p ${MPATH}
      gzip -c ./man/atomsk >${MPATH}/atomsk.1.gz
    
    
    else
    
      echo ">>> Installation aborted."
    
    fi

  else

    echo "X!X ERROR: this script must be run as root or with sudo."

  fi

fi