#!/bin/bash

# This script installs Atomsk on the local computer.
# If run with super-user rights (su or sudo), then the
# program and documentation are installed in /usr/local/.
# Otherwise they are installed in the user's home directory.

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
    
      echo ">>> Installation cancelled."
    
    fi

  else
  
    # Install program in user's home directory

    BINPATH=~/bin/atomsk/
    DOCPATH=~/bin/atomsk/doc/

    echo "<!> INFO: run this script with super-user rights (su or sudo)"
    echo "         to install Atomsk system-wide."
    echo ""
    echo "<?> Atomsk will be installed in ${BINPATH}. Continue? (y/n)"
    read answer

    if [ "${answer}" = "y" ] ; then
      # Create the target folder
      mkdir -p ${BINPATH}
      
      # System configuration file
      cp -rf ./etc/atomsk.conf ~/.config/

      # Atomsk binary
      chmod +x atomsk
      cp atomsk ${BINPATH}

      # Atomsk tools
      chmod +x ./tools/*.sh
      cp ./tools/* ${BINPATH}
      echo ">>> The program was successfuly installed in ${BINPATH}"
      echo "    To run it, enter 'atomsk' in a terminal."

      # Atomsk documentation
      mkdir -p ${DOCPATH}
      rm -rf ${DOCPATH}/*
      cp -rf ./doc/* ${DOCPATH}
      chmod -R a+r ${DOCPATH}
      echo ">>> The html documentation was installed. You may read it by entering the"
      echo "    following address in your Web browser: ${DOCPATH}index.html"
      
      # Create alias
      n=$(grep "atomsk" ~/.bashrc | wc -l)
      if [ $n -eq 0 ] ; then
        echo "export PATH=\"\$PATH:${BINPATH}\"" >> ~/.bashrc
        echo ">>> ${BINPATH} was added to your PATH environment variable."
      fi
    
    else
    
      echo ">>> Installation cancelled."

    fi

  fi

fi
