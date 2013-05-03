#!/bin/csh
#  file assignments.
#
#  All binary files should be put on a node's local disk ($SCR directory),
#  for the highest speed access possible.  These .Fxx files are typically
#  not saved for the next run, but they may be big and/or I/O intensive.
#
#  It is convenient to write ASCII output files (PUNCH, RESTART, TRAJECT,
#  and MAKEFP) to the user's permanent disk, on your file server.  They
#  are small, written only by the master process, and are useful outputs
#  for further runs.
#
#                        ASCII input files
#             You must edit a, but will probably skip b+c.
#  Some data files may be read by a run, each is read only once, so
#  that storage of one (1) copy on your file server is appropriate.
#  a) AUXDATA is directory of data sets containing
#     Note: change only AUXDATA, not ERICFMT,MCPPATH,BASPATH,QUANPOL!
#        1. a file of Fm(t) data for ERI computations
#        2. a BASES subdirectory, containing files of some basis sets
#        3. a MCP subdirectory, containing files of MCP bases and potentials
#        4. data sets for the Quantum chemistry Polarizable force field
#        5. a EFP subdirectory, containing standard EFP2 potentials
#  b) The EXTBAS file contains any user-supplied basis sets.
#  c) The NUCBAS or POSBAS files are nuclear or positron basis sets,
#     used by the NEO method.  See NEO's documentation for more details.
#  d) there are 3 places where you might want to uncomment a 'set echo',
#     to see all the file definitions (one is just below).
#
#---quiet---set echo
setenv AUXDATA $GMSPATH/auxdata
setenv  EXTBAS /dev/null
setenv  NUCBAS /dev/null
setenv  POSBAS /dev/null
#
setenv ERICFMT $AUXDATA/ericfmt.dat
setenv MCPPATH $AUXDATA/MCP
setenv BASPATH $AUXDATA/BASES
setenv QUANPOL $AUXDATA/QUANPOL
setenv  MAKEFP $USERSCR/$JOB.efp
setenv   GAMMA $USERSCR/$JOB.gamma
setenv TRAJECT $USERSCR/$JOB.trj
setenv RESTART $USERSCR/$JOB.rst
setenv   INPUT $SCR/$JOB.F05
setenv   PUNCH $USERSCR/$JOB.dat
setenv  AOINTS $SCR/$JOB.F08
setenv  MOINTS $SCR/$JOB.F09
setenv DICTNRY $SCR/$JOB.F10
setenv DRTFILE $SCR/$JOB.F11
setenv CIVECTR $SCR/$JOB.F12
setenv CASINTS $SCR/$JOB.F13
setenv  CIINTS $SCR/$JOB.F14
setenv  WORK15 $SCR/$JOB.F15
setenv  WORK16 $SCR/$JOB.F16
setenv CSFSAVE $SCR/$JOB.F17
setenv FOCKDER $SCR/$JOB.F18
setenv  WORK19 $SCR/$JOB.F19
setenv  DASORT $SCR/$JOB.F20
setenv DIABDAT $SCR/$JOB.F21
setenv DFTINTS $SCR/$JOB.F21
setenv DFTGRID $SCR/$JOB.F22
setenv  JKFILE $SCR/$JOB.F23
setenv  ORDINT $SCR/$JOB.F24
setenv  EFPIND $SCR/$JOB.F25
setenv PCMDATA $SCR/$JOB.F26
setenv PCMINTS $SCR/$JOB.F27
setenv SVPWRK1 $SCR/$JOB.F26
setenv SVPWRK2 $SCR/$JOB.F27
setenv COSCAV  $SCR/$JOB.F26
setenv COSDATA $USERSCR/$JOB.cosmo
setenv COSPOT  $USERSCR/$JOB.pot
setenv  MLTPL  $SCR/$JOB.F28
setenv  MLTPLT $SCR/$JOB.F29
setenv  DAFL30 $SCR/$JOB.F30
setenv  SOINTX $SCR/$JOB.F31
setenv  SOINTY $SCR/$JOB.F32
setenv  SOINTZ $SCR/$JOB.F33
setenv  SORESC $SCR/$JOB.F34
#   35 is used by RESTART, see above
setenv GCILIST $SCR/$JOB.F37
setenv HESSIAN $SCR/$JOB.F38
setenv QMMMTEI $SCR/$JOB.F39
setenv SOCCDAT $SCR/$JOB.F40
setenv  AABB41 $SCR/$JOB.F41
setenv  BBAA42 $SCR/$JOB.F42
setenv  BBBB43 $SCR/$JOB.F43
setenv  REMD   $SCR/$JOB.F44
setenv  MCQD50 $SCR/$JOB.F50
setenv  MCQD51 $SCR/$JOB.F51
setenv  MCQD52 $SCR/$JOB.F52
setenv  MCQD53 $SCR/$JOB.F53
setenv  MCQD54 $SCR/$JOB.F54
setenv  MCQD55 $SCR/$JOB.F55
setenv  MCQD56 $SCR/$JOB.F56
setenv  MCQD57 $SCR/$JOB.F57
setenv  MCQD58 $SCR/$JOB.F58
setenv  MCQD59 $SCR/$JOB.F59
setenv  MCQD60 $SCR/$JOB.F60
setenv  MCQD61 $SCR/$JOB.F61
setenv  MCQD62 $SCR/$JOB.F62
setenv  MCQD63 $SCR/$JOB.F63
setenv  MCQD64 $SCR/$JOB.F64
setenv NMRINT1 $SCR/$JOB.F61
setenv NMRINT2 $SCR/$JOB.F62
setenv NMRINT3 $SCR/$JOB.F63
setenv NMRINT4 $SCR/$JOB.F64
setenv NMRINT5 $SCR/$JOB.F65
setenv NMRINT6 $SCR/$JOB.F66
setenv DCPHFH2 $SCR/$JOB.F67
setenv DCPHF21 $SCR/$JOB.F68
setenv ELNUINT $SCR/$JOB.F67
setenv NUNUINT $SCR/$JOB.F68
setenv   GVVPT $SCR/$JOB.F69
setenv NUMOIN  $SCR/$JOB.F69
setenv NUMOCAS $SCR/$JOB.F70
setenv NUELMO  $SCR/$JOB.F71
setenv NUELCAS $SCR/$JOB.F72

#    next files are for RI-MP2
setenv RIVMAT  $SCR/$JOB.F51
setenv RIT2A   $SCR/$JOB.F52
setenv RIT3A   $SCR/$JOB.F53
setenv RIT2B   $SCR/$JOB.F54
setenv RIT3B   $SCR/$JOB.F55

#    Next files are for GMCQDPT
setenv GMCREF $SCR/$JOB.F70
setenv GMCO2R $SCR/$JOB.F71
setenv GMCROC $SCR/$JOB.F72
setenv GMCOOC $SCR/$JOB.F73
setenv GMCCC0 $SCR/$JOB.F74
setenv GMCHMA $SCR/$JOB.F75
setenv GMCEI1 $SCR/$JOB.F76
setenv GMCEI2 $SCR/$JOB.F77
setenv GMCEOB $SCR/$JOB.F78
setenv GMCEDT $SCR/$JOB.F79
setenv GMCERF $SCR/$JOB.F80
setenv GMCHCR $SCR/$JOB.F81
setenv GMCGJK $SCR/$JOB.F82
setenv GMCGAI $SCR/$JOB.F83
setenv GMCGEO $SCR/$JOB.F84
setenv GMCTE1 $SCR/$JOB.F85
setenv GMCTE2 $SCR/$JOB.F86
setenv GMCHEF $SCR/$JOB.F87
setenv GMCMOL $SCR/$JOB.F88
setenv GMCMOS $SCR/$JOB.F89
setenv GMCWGT $SCR/$JOB.F90
setenv GMCRM2 $SCR/$JOB.F91
setenv GMCRM1 $SCR/$JOB.F92
setenv GMCR00 $SCR/$JOB.F93
setenv GMCRP1 $SCR/$JOB.F94
setenv GMCRP2 $SCR/$JOB.F95
setenv GMCVEF $SCR/$JOB.F96
setenv GMCDIN $SCR/$JOB.F97
setenv GMC2SZ $SCR/$JOB.F98
setenv GMCCCS $SCR/$JOB.F99

#    Next files are used only during closed shell coupled cluster runs.
#    Display the numerous definitions iff they are going to be used.
unset echo
#---quiet---set cctyp=`grep -i 'CCTYP[(=]' $SCR/$JOB.F05 | wc -l`
#---quiet---if ($cctyp > 0) set echo
setenv  CCREST $SCR/$JOB.F70
setenv  CCDIIS $SCR/$JOB.F71
setenv  CCINTS $SCR/$JOB.F72
setenv CCT1AMP $SCR/$JOB.F73
setenv CCT2AMP $SCR/$JOB.F74
setenv CCT3AMP $SCR/$JOB.F75
setenv    CCVM $SCR/$JOB.F76
setenv    CCVE $SCR/$JOB.F77
setenv CCQUADS $SCR/$JOB.F78
setenv QUADSVO $SCR/$JOB.F79
setenv EOMSTAR $SCR/$JOB.F80
setenv EOMVEC1 $SCR/$JOB.F81
setenv EOMVEC2 $SCR/$JOB.F82
setenv  EOMHC1 $SCR/$JOB.F83
setenv  EOMHC2 $SCR/$JOB.F84
setenv EOMHHHH $SCR/$JOB.F85
setenv EOMPPPP $SCR/$JOB.F86
setenv EOMRAMP $SCR/$JOB.F87
setenv EOMRTMP $SCR/$JOB.F88
setenv EOMDG12 $SCR/$JOB.F89
setenv    MMPP $SCR/$JOB.F90
setenv   MMHPP $SCR/$JOB.F91
setenv MMCIVEC $SCR/$JOB.F92
setenv MMCIVC1 $SCR/$JOB.F93
setenv MMCIITR $SCR/$JOB.F94
setenv  EOMVL1 $SCR/$JOB.F95
setenv  EOMVL2 $SCR/$JOB.F96
setenv EOMLVEC $SCR/$JOB.F97
setenv  EOMHL1 $SCR/$JOB.F98
setenv  EOMHL2 $SCR/$JOB.F99
setenv  CCVVVV $SCR/$JOB.F80
#
#    Next files are used only during open shell coupled cluster runs.
#
setenv AMPROCC $SCR/$JOB.F70
setenv ITOPNCC $SCR/$JOB.F71
setenv FOCKMTX $SCR/$JOB.F72
setenv  LAMB23 $SCR/$JOB.F73
setenv   VHHAA $SCR/$JOB.F74
setenv   VHHBB $SCR/$JOB.F75
setenv   VHHAB $SCR/$JOB.F76
setenv    VMAA $SCR/$JOB.F77
setenv    VMBB $SCR/$JOB.F78
setenv    VMAB $SCR/$JOB.F79
setenv    VMBA $SCR/$JOB.F80
setenv  VHPRAA $SCR/$JOB.F81
setenv  VHPRBB $SCR/$JOB.F82
setenv  VHPRAB $SCR/$JOB.F83
setenv  VHPLAA $SCR/$JOB.F84
setenv  VHPLBB $SCR/$JOB.F85
setenv  VHPLAB $SCR/$JOB.F86
setenv  VHPLBA $SCR/$JOB.F87
setenv    VEAA $SCR/$JOB.F88
setenv    VEBB $SCR/$JOB.F89
setenv    VEAB $SCR/$JOB.F90
setenv    VEBA $SCR/$JOB.F91
setenv   VPPPP $SCR/$JOB.F92
setenv INTERM1 $SCR/$JOB.F93
setenv INTERM2 $SCR/$JOB.F94
setenv INTERM3 $SCR/$JOB.F95
setenv ITSPACE $SCR/$JOB.F96
setenv INSTART $SCR/$JOB.F97
setenv  ITSPC3 $SCR/$JOB.F98
#
#    Next files are used only during elongation method runs.
#    Display the numerous definitions iff they are going to be used.
unset echo
set elgtyp = `grep -i NELONG= $SCR/$JOB.F05 | wc -l`
if ($elgtyp > 0) then
    set ELGNAME=$4
    if (null$4 == null) set ELGNAME=ELGFILE
    set echo
    setenv AOINTS   $SCR/$ELGNAME.F08
    setenv ELGDOS   $USERSCR/$JOB.ldos
    setenv ELGDAT   $SCR/$ELGNAME.F71
    setenv ELGPAR   $SCR/$ELGNAME.F72
    setenv ELGCUT   $SCR/$ELGNAME.F74
    setenv ELGVEC   $SCR/$ELGNAME.F75
    setenv EGINTA   $SCR/$ELGNAME.F77
    setenv EGINTB   $SCR/$ELGNAME.F78
    setenv EGTDHF   $SCR/$ELGNAME.F79
    setenv EGTEST   $SCR/$ELGNAME.F80
    unset echo
endif
#
#    Next files are used only during extended TDHF package runs.
#    Display the numerous definitions iff they are going to be used.
unset echo
#---quiet---set txtyp=`grep -i RUNTYP=TDHFX $SCR/$JOB.F05 | wc -l`
#---quiet---if ($txtyp > 0) set echo
setenv  OLI201 $SCR/$JOB.F201
setenv  OLI202 $SCR/$JOB.F202
setenv  OLI203 $SCR/$JOB.F203
setenv  OLI204 $SCR/$JOB.F204
setenv  OLI205 $SCR/$JOB.F205
setenv  OLI206 $SCR/$JOB.F206
setenv  OLI207 $SCR/$JOB.F207
setenv  OLI208 $SCR/$JOB.F208
setenv  OLI209 $SCR/$JOB.F209
setenv  OLI210 $SCR/$JOB.F210
setenv  OLI211 $SCR/$JOB.F211
setenv  OLI212 $SCR/$JOB.F212
setenv  OLI213 $SCR/$JOB.F213
setenv  OLI214 $SCR/$JOB.F214
setenv  OLI215 $SCR/$JOB.F215
setenv  OLI216 $SCR/$JOB.F216
setenv  OLI217 $SCR/$JOB.F217
setenv  OLI218 $SCR/$JOB.F218
setenv  OLI219 $SCR/$JOB.F219
setenv  OLI220 $SCR/$JOB.F220
setenv  OLI221 $SCR/$JOB.F221
setenv  OLI222 $SCR/$JOB.F222
setenv  OLI223 $SCR/$JOB.F223
setenv  OLI224 $SCR/$JOB.F224
setenv  OLI225 $SCR/$JOB.F225
setenv  OLI226 $SCR/$JOB.F226
setenv  OLI227 $SCR/$JOB.F227
setenv  OLI228 $SCR/$JOB.F228
setenv  OLI229 $SCR/$JOB.F229
setenv  OLI230 $SCR/$JOB.F230
setenv  OLI231 $SCR/$JOB.F231
setenv  OLI232 $SCR/$JOB.F232
setenv  OLI233 $SCR/$JOB.F233
setenv  OLI234 $SCR/$JOB.F234
setenv  OLI235 $SCR/$JOB.F235
setenv  OLI236 $SCR/$JOB.F236
setenv  OLI237 $SCR/$JOB.F237
setenv  OLI238 $SCR/$JOB.F238
setenv  OLI239 $SCR/$JOB.F239
unset echo

#    Next files are used only during divide-and-conquer runs
setenv   DCSUB $SCR/$JOB.F250
setenv   DCVEC $SCR/$JOB.F251
setenv   DCEIG $SCR/$JOB.F252
setenv   DCDM  $SCR/$JOB.F253
setenv   DCDMO $SCR/$JOB.F254
setenv   DCQ   $SCR/$JOB.F255
setenv   DCW   $SCR/$JOB.F256
setenv   DCEDM $SCR/$JOB.F257

#    Next files are used only during LMO hyperpolarizability analysis
setenv LHYPWRK $SCR/$JOB.F297
setenv LHYPWK2 $SCR/$JOB.F298
setenv BONDDPF $SCR/$JOB.F299

#    Next two values are used only within the VB2000 add-on code
setenv VB2000PATH $GMSPATH/vb2000
setenv GMSJOBNAME $JOB

# Next files are used during EFMO runs
set efmo = `grep -i IEFMO= $SCR/$JOB.F05 | wc -l`
if ($efmo > 0) then
  setenv EFMOI $SCR/$JOB.F102
  setenv EFMOF $SCR/$JOB.F103
endif

#    Next files are used for explicitly correlated methods 
set pt2r12=`egrep -i '(PT212=.TRUE.|PT2R12=.T.)' $SCR/$JOB.F05 | wc -l`
if ($pt2r12 > 0) then
 set interface=`egrep -i  '(RUNR12=.T.|RUNTYP=.TRUE.|SINGLS=.T.|SINGLES=.TRUE)' $SCR/$JOB.F05 | wc -l`
 if ($interface > 0) then
  set echo
  setenv R12INP $SCR/$JOB
  setenv PT2BAS $SCR/$JOB.F300
  unset echo
 else 
  set echo
  setenv R12INP $USERSCR/$JOB
  setenv PT2BAS $SCR/$JOB.F300
  unset echo
 endif
endif

#    Next files are used only during CIM runs
setenv CIMFILE $USERSCR/$JOB.cim
setenv CIMDMN  $USERSCR/$JOB.dmn
setenv CIMDCT  $SCR/$JOB.Fcdt
setenv CIMAOI  $SCR/$JOB.Fcao
setenv CIMMOI  $SCR/$JOB.Fcmo
