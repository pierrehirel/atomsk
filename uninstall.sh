#!/bin/bash

# This script uninstalls Atomsk from system directories.
# Note that Atomsk should have been installed with its script "install.sh".
# This script must be run with root permissions.

BINPATH=/usr/local/bin/
DOCPATH=/usr/local/share/doc/
MPATH=/usr/local/share/man/man1/

if [ -w ${BINPATH} ] ; then

  echo "<?> Atomsk will be REMOVED from ${BINPATH}. Continue? (y/n)"
  read answer

  if [ "${answer}" = "y" ] ; then
    # System configuration file
    rm -rf ./etc/atomsk.conf /etc/

    # atomsk binary
    rm -f ${BINPATH}/atomsk

    # atomsk tools
    rm -f ${BINPATH}/cfg_setA.sh
    rm -f ${BINPATH}/dat_mulvec.sh
    rm -f ${BINPATH}/dat_rm0.sh
    rm -f ${BINPATH}/lmp_atom2charge.sh
    rm -f ${BINPATH}/lmp_charge2atom.sh
    rm -f ${BINPATH}/lmp_ortho2tri.sh
    rm -f ${BINPATH}/qepw_bohr.sh

    # atomsk documentation
    rm -rf ${DOCPATH}/atomsk

    # atomsk man page
    rm -f ${MPATH}/atomsk.1.gz

    echo ">>> Program was successfuly removed."
  
  else
  
    echo ">>> Uninstallation aborted."
  
  fi

else

  echo "X!X ERROR: this script must be run as root or with sudo."

fi
