#!/bin/bash

# This script uninstalls Atomsk from system directories.
# Note that Atomsk should have been installed with its script "install.sh".

clear
printf "___________________________________________________________\n"
printf "\e[1m              Atomsk Uninstaller\e[0m\n"
printf "_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _\n"

# Get path to atomsk executable
BINPATH="$(which atomsk)"

if [ "$BINPATH" = "" ] || [ ! -e "${BINPATH}" ] ; then
  
  # Executable not found in PATH
  echo "X!X ERROR: Atomsk is not installed on this system, or not"
  echo "          in one of the folders declared in PATH."

elif [ -w ${BINPATH} ] ; then

  BINPATH="$(dirname ${BINPATH} )"

  # User has rights to access $BINPATH directory
  printf "<?> Atomsk will be REMOVED from ${BINPATH}. Continue? (y/n) "
  read answer

  if [ "${answer}" = "y" ] ; then
  
    # Delete Atomsk binary
    rm -f ${BINPATH}/atomsk

    if [ "$BINPATH" = "/usr/local/bin" ] ; then

      DOCPATH=/usr/local/share/doc/
      MPATH=/usr/local/share/man/man1/
      
      # Delete system configuration file
      rm -f /etc/atomsk.conf

      # Atomsk tools
      rm -f ${BINPATH}/cfg_setA.sh
      rm -f ${BINPATH}/dat_mulvec.sh
      rm -f ${BINPATH}/dat_rm0.sh
      rm -f ${BINPATH}/lmp_atom2charge.sh
      rm -f ${BINPATH}/lmp_charge2atom.sh
      rm -f ${BINPATH}/lmp_ortho2tri.sh
      rm -f ${BINPATH}/qepw_bohr.sh

      # Atomsk documentation
      rm -rf ${DOCPATH}/atomsk/

      # Atomsk man page
      rm -f ${MPATH}/atomsk.1.gz

      echo ">>> Program was successfuly removed."
      
    elif [ "$BINPATH" = "~/bin/atomsk" ] ; then

      # Remove the whole folder  ~/bin/atomsk/
      rm -rf ${BINPATH}

      # Display information for the user
      echo ">>> Program was successfuly removed."
      
      n=$(grep "atomsk" ~/.bashrc | wc -l)
      if [ $n -gt 0 ] ; then
        echo "<i> INFO: your file '~/.bashrc' still contains an alias or the path"
        echo "          related to Atomsk. You may edit this file to remove it."
      fi
      
    else

      # Atomsk was installed elsewhere, like /opt/ or another location
      # User probably put it there intentionally
      # Only delete Atomsk executable (that was already done)
      echo ">>> Program was successfuly removed."
        
    fi
    
    if [ -e "~/.config/atomsk.conf" ] ; then
    
      echo "<i> INFO: a file '~/.config/atomsk.conf' still exists in your personal folder."
      echo "          You may either leave it, or remove it with the following command:"
      echo "            rm -f ~/.config/atomsk.conf"
    
    fi
    
  else

    # User did not answer "y"
    echo ">>> Uninstallation cancelled."
  
  fi
  
else

  # User does not have the right to access $BINPATH  
  echo "X!X ERROR: Atomsk is installed in ${BINPATH}, but you do not have"
  echo "          the rights to un-install it. Run this script with super-user"
  echo "          rights (su or sudo) to un-install it."
  exit

fi
