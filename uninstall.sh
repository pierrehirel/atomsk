#!/bin/bash

# This script uninstalls Atomsk from system directories.
# Note that Atomsk should have been installed with its script "install.sh".

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
      echo ">>> You may also delete the file ~/.config/atomsk.conf"
      echo "    and edit your ~/.bashrc file to remove the folder"
      echo "    ${BINPATH} from your PATH environment variable."
      
    else

      # Atomsk was installed elsewhere, like /opt/ or another location
      # User probably put it there intentionally
      # Only delete Atomsk executable (that was already done)
      echo ">>> Program was successfuly removed."
        
    fi
    
  else

    # User did not answer "y"
    echo ">>> Uninstallation cancelled."
  
  fi
  
else

  # User does not have the right to access £BINPATH  
  echo "X!X ERROR: Atomsk is installed in ${BINPATH}, but you do not have"
  echo "          the rights to un-install it. Run this script with super-user"
  echo "          rights (su or sudo) to un-install it."
  exit

fi
