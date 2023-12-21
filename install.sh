#!/bin/bash

# This script installs Atomsk on the local computer.
# If run with super-user rights (su or sudo), then the
# program and documentation are installed in /usr/local/.
# Otherwise they are installed in the user's home directory.

web="https://atomsk.univ-lille.fr"

BINPATH=/usr/local/bin/
DOCPATH=/usr/local/share/doc/
MPATH=/usr/local/share/man/man1/

clear

printf "___________________________________________________________\n"
printf "\e[1m              Atomsk Installation Setup\e[0m\n"
printf "_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _\n"

# Check if executable is present; otherwise, try to compile it
if [ ! -e 'atomsk' ] ; then

  if [ -e "src/atomsk" ] ; then
  
    cp src/atomsk .
    
  else

    echo "<i> INFO: the executable 'atomsk' does not exist in current directory."

    if [ ! -d "src" ] ; then

      printf "<?> Do you wish to download Atomsk from the Web site? (y/n) "
      read answer

      if [ "${answer}" = "y" ] ; then

        wget ${web}/version.txt
        if [ -e "version.txt" ] ; then
          ver=$(more version.txt)
          wget ${web}/code/atomsk_b${ver}.tar.gz
          rm version.txt
          if [ -e "atomsk_b${ver}.tar.gz" ] ; then
            tar -xzf atomsk_b${ver}.tar.gz
            cd atomsk_b${ver}
          else
            echo "X!X ERROR: unable to reach ${web} , aborting."
            exit
          fi
        else
          echo "X!X ERROR: unable to reach ${web} , aborting."
          exit
        fi

      else

        echo ">>> Installation cancelled."
        exit

      fi

    fi
    
    if [ -d "src" ] ; then
    
      printf "<?> Do you wish to compile Atomsk from the source code? (y/n) "
      read answer

      if [ "${answer}" = "y" ] ; then
        
        cd src
        make -j3 atomsk
        cd ..
        
        if [ -e "src/atomsk" ] ; then
        
          cp src/atomsk .
        
        else
        
          echo "X!X ERROR: compilation failed. Please read the documentation before compiling again,"
          echo "          or download a binary version from ${web}."
        
        fi


      else

        echo ">>> Installation cancelled."
        exit
    
      fi
    
    else
    
      echo "X!X ERROR: impossible to install Atomsk. Please download the program from ${web}"
      exit
    
    fi
    
  fi

fi


# Check if we have super-user rights
if [ ! "$UID" -eq 0 ] ; then

  echo "<!> INFO: I must have super-user rights (su) to install Atomsk in your system (${BINPATH})."
  echo "<!>       If you choose no, then Atomsk will be installed in the /bin/ directory in your home (${HOME}/bin)."
  printf "<?> Do you wish to activate super-user rights to install Atomsk system-wide? (y/n) "
  read answer

  if [ "${answer}" = "y" ] ; then
  
    exec sudo "$0" "$@"
  
  fi

fi


# Copy program and files into user's computer
if [ "$UID" -eq 0 ] ; then

  printf "<?> Atomsk will be installed in ${BINPATH}. Continue? (y/n) "
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
    echo ""
    echo ">>> The program was successfuly installed in ${BINPATH}"
    echo "    To run it, enter 'atomsk' in a terminal."

    # Atomsk documentation
    mkdir -p ${DOCPATH}/atomsk
    rm -rf ${DOCPATH}/atomsk/*
    cp -rf ./doc/* ${DOCPATH}/atomsk/
    chmod -R a+r ${DOCPATH}/atomsk/
    echo ""
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

  echo ""
  printf "<?> Atomsk will be installed in ${BINPATH}. Continue? (y/n) "
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
    echo ""
    echo ">>> The program was successfuly installed in ${BINPATH}"
    echo "    To run it, enter 'atomsk' in a terminal."

    # Atomsk documentation
    mkdir -p ${DOCPATH}
    rm -rf ${DOCPATH}/*
    cp -rf ./doc/* ${DOCPATH}
    chmod -R a+r ${DOCPATH}
    echo ""
    echo ">>> The html documentation was installed. You may read it by entering the"
    echo "    following address in your Web browser: ${DOCPATH}index.html"
    
    # Create alias
    n=$(grep "atomsk" ~/.bashrc | wc -l)
    if [ $n -eq 0 ] ; then
      echo "export PATH=\"\$PATH:${BINPATH}\"" >> ~/.bashrc
      source ~/.bashrc
    echo ""
      echo ">>> ${BINPATH} was added to your PATH environment variable."
    fi
  
  else
  
    echo ">>> Installation cancelled."

  fi

fi

printf "___________________________________________________________\n"
