#!/bin/csh
#
set GMSDIR=/u1/mike/gamess
#
#       count the FORTRAN
#
chdir $GMSDIR/source
@ nftotl=0
@ nftotc=0
@ nftotm=0
@ nffile=0
foreach FILE (*.src)
   echo -n $FILE:r..
   set nline=`wc -l $FILE`
   set nline=$nline[1]
   set ncomm=`grep -i ^c $FILE | wc -l`
   set nmach=`grep -i ^\* $FILE | wc -l`
   @ nftotl = $nftotl + $nline
   @ nftotc = $nftotc + $ncomm
   @ nftotm = $nftotm + $nmach
   @ nffile++
end
echo " "
#
#       count the quiche
#
@ nctotl=0
@ ncfile=0
foreach FILE (*.c)
   echo -n $FILE:r..
   set nline=`wc -l $FILE`
   set nline=$nline[1]
   @ nctotl = $nctotl + $nline
   @ ncfile++
end
echo " "
#
#       count the documentation
#
chdir $GMSDIR
@ ndtotl=0
@ ndfile=0
foreach FILE (*.DOC)
   echo -n $FILE:r..
   set nline=`wc -l $FILE`
   set nline=$nline[1]
   @ ndtotl = $ndtotl + $nline
   @ ndfile++
end
echo " "
#
#       count the message passing library
#
chdir $GMSDIR/ddi/src
@ nhdditotl=0
@ nhddifile=0
@ ncdditotl=0
@ ncddifile=0
foreach FILE (*.h)
   echo -n $FILE:r..
   set nline=`wc -l $FILE`
   set nline=$nline[1]
   @ nhdditotl = $nhdditotl + $nline
   @ nhddifile++
end
echo " "
foreach FILE (*.c)
   echo -n $FILE:r..
   set nline=`wc -l $FILE`
   set nline=$nline[1]
   @ ncdditotl = $ncdditotl + $nline
   @ ncddifile++
end
#
chdir $GMSDIR/ddi/kickoff
foreach FILE (*.c)
   echo -n $FILE:r..
   set nline=`wc -l $FILE`
   set nline=$nline[1]
   @ ncdditotl = $ncdditotl + $nline
   @ ncddifile++
end
echo " "
#
#       print the results
#
echo " "
echo "GAMESS has $nftotl lines of FORTRAN, of which $nftotc are comments,"
echo "and $nftotm are machine dependent lines, in $nffile source files."
echo "GAMESS has $nctotl lines in its $ncfile quiche (C language) files."
echo " "
echo "DDI has $nhdditotl lines of include statements, in $nhddifile files."
echo "DDI has $ncdditotl lines of quiche, in $ncddifile C language files."
echo " "
echo "GAMESS has $ndtotl lines in its $ndfile *.DOC files."
echo " "
