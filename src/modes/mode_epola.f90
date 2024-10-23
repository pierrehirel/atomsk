MODULE mode_epola
!
!**********************************************************************************
!* MODE_EPOLA                                                                     *
!**********************************************************************************
!* This mode is meant to compute the electronic polarization in a core-shell      *
!* system, i.e. the difference in position between an atomic core and its         *
!* associated electronic shell:                                                   *
!*          P = 1/2 * (|Q|+|q|) * (R-r)                                           *
!* where Q and R are the charge and position of the core, and q and r the         *
!* charge and position of the shell.                                              *
!* Note that for rigid ions (or neutral atoms) the polarization will              *
!* always be zero by definition.                                                  *
!**********************************************************************************
!* (C) Dec. 2011 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 19 Sept. 2024                                    *
!**********************************************************************************
!* This program is free software: you can redistribute it and/or modify           *
!* it under the terms of the GNU General Public License as published by           *
!* the Free Software Foundation, either version 3 of the License, or              *
!* (at your option) any later version.                                            *
!*                                                                                *
!* This program is distributed in the hope that it will be useful,                *
!* but WITHOUT ANY WARRANTY; without even the implied warranty of                 *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  *
!* GNU General Public License for more details.                                   *
!*                                                                                *
!* You should have received a copy of the GNU General Public License              *
!* along with this program.  If not, see <http://www.gnu.org/licenses/>.          *
!**********************************************************************************
!
!
!Load modules
USE comv        !global variables
USE constants   !math and physics constants
USE messages
USE files_msg
USE subroutines !subroutines for this program
USE out_cfg     !For writing CFG file
USE out_xsf     !For writing XSF file
!
!
CONTAINS
!
SUBROUTINE E_POLA(H,P,S,AUXNAMES,AUX,comment,outputfile)
!
!Declare variables
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=128):: outputfile, pxsf, pcfg  !Output file names
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: fakeAUXNAMES !fake names for ionic polarizations
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
INTEGER:: i
INTEGER:: q, qs !index of charges for cores and shells in AUX
REAL(dp), DIMENSION(3,3):: H !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: ionpola  !polarization vector+amplitude for each ion
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P, S  !positions of cores, shells
REAL(dp), DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties of atoms
!
!
!Initialize variables
temp = TRIM(outputfile)//'_PE'
pxsf = TRIM(temp)//'.xsf'
pcfg = TRIM(temp)//'.cfg'
q = 0
qs = 0
ALLOCATE(ionpola(SIZE(P,1),4))
ionpola(:,:) = 0.d0
ALLOCATE( fakeAUXNAMES(4) )
ALLOCATE(comment(1))
 comment(1)=''
!
!
msg = 'Entering E_POLA'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(4048,(/""/),(/0.d0/))
!
!
!
100 CONTINUE
!Check that we dispose of all the data for the calculation
!
!Check that shells are present
IF( .NOT.ALLOCATED(S)   .OR. SIZE(S,1)==0 ) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(4817,(/""/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Check that auxiliary properties are present
IF( .NOT.ALLOCATED(AUX)      .OR. SIZE(AUX,1)<=0 .OR.           &
  & .NOT.ALLOCATED(AUXNAMES) .OR. SIZE(AUXNAMES)<=0    ) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(4807,(/""/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Look for charges of cores and shells in AUXNAMES
DO i=1,SIZE(AUXNAMES)
  IF( TRIM(AUXNAMES(i))=="q" ) THEN
    q=i
  ELSEIF( TRIM(AUXNAMES(i))=="qs" ) THEN
    qs=i
  ENDIF
ENDDO
IF( q==0 .OR. qs==0 ) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(4816,(/""/),(/DBLE(q),DBLE(qs)/))
  GOTO 1000
ENDIF
!
!
!
200 CONTINUE
!Compute the polarization vectors
DO i=1,SIZE(P,1)
  IF( NINT(S(i,4)).NE.0 .AND. AUX(i,qs).NE.0.d0 ) THEN
    !This ion has a shell: compute polarization
    ionpola(i,1:3) = 0.5d0*( DABS(AUX(i,q))+DABS(AUX(i,qs)) ) * (P(i,1:3)-S(i,1:3))
    ionpola(i,4) = VECLENGTH(ionpola(i,1:3))
  ELSE
    !If shell is defined but has zero charge, it is weird: output a warning
    IF( AUX(i,qs)==0.d0 ) THEN
      nwarn = nwarn+1
      CALL ATOMSK_MSG(4707,(/""/),(/DBLE(i)/))
    ELSE
      !This ion has no shell ("rigid ion") => polarization is zero
      ionpola(i,:) = 0.d0
    ENDIF
  ENDIF
ENDDO
!
!
!
300 CONTINUE
!Output
 comment(1) = '# Electronic polarization computed with atomsk'
fakeAUXNAMES(1) = "Px"
fakeAUXNAMES(2) = "Py"
fakeAUXNAMES(3) = "Pz"
fakeAUXNAMES(4) = "|P|"
!
!Output to CFG
IF(.NOT.overw) CALL CHECKFILE(pcfg,'writ')
CALL WRITE_CFG(H,P,comment,fakeAUXNAMES,ionpola,pcfg)
CALL ATOMSK_MSG(3002,(/pcfg,"CFG",FILE_SIZE(outputfile)/),(/0.d0/))
!
!Output to XSF (rename aux.properties to trick XSF module)
fakeAUXNAMES(1) = "fx"
fakeAUXNAMES(2) = "fy"
fakeAUXNAMES(3) = "fz"
fakeAUXNAMES(4) = ""
IF(.NOT.overw) CALL CHECKFILE(pxsf,'writ')
CALL WRITE_XSF(H,P,comment,fakeAUXNAMES,ionpola,pxsf)
CALL ATOMSK_MSG(3002,(/pxsf,"XSF",FILE_SIZE(outputfile)/),(/0.d0/))
!
!
CALL ATOMSK_MSG(4049,(/""/),(/0.d0/))
!
!
!
1000 CONTINUE
!
!
!
END SUBROUTINE E_POLA
!
!
END MODULE
