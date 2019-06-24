#!/bin/bash

# This script uninstalls Atomsk from system directories.
# Note that Atomsk should have been installed with its script "install.sh".
# If run with super-user permissions (su or sudo), then atomsk
# and associated tools and doc are removed from /usr/local/ directory.
# Otherwise, the folder ~/bin/atomsk/ in user's home directory is deleted.

BINPATH=/usr/local/bin/
DOCPATH=/usr/local/share/doc/
MPATH=/usr/local/share/man/man1/

if [ -w ${BINPATH} ] ; then

  echo "<?> Atomsk will be REMOVED from ${BINPATH}. Continue? (y/n)"
  read answer

  if [ "${answer}" = "y" ] ; then
    # System configuration file
    rm -f /etc/atomsk.conf

    # Atomsk binary
    rm -f ${BINPATH}/atomsk

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
  
  else
  
    echo ">>> Uninstallation cancelled."
  
  fi

else

  BINPATH=~/bin/atomsk/
  DOCPATH=~/bin/atomsk/doc/

  echo "<?> Atomsk will be REMOVED from ${BINPATH}. Continue? (y/n)"
  read answer

  if [ "${answer}" = "y" ] ; then

    # Remove the whole folder /atomsk/ in user's ~/bin/
    rm -rf ${BINPATH}

    echo ">>> Program was successfuly removed."
    echo ">>> You may also delete the file ~/.config/atomsk.conf"
    echo "    and edit your ~/.bashrc file to remove the folder"
    echo "    ${BINPATH} from your PATH environment variable."
  
  else
  
    echo ">>> Uninstallation cancelled."
  
  fi

fi
