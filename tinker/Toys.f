C 21 Apr 10 - NA,KK - allow FMO usage, RESTRAIN-POSITION changes
C 21 Apr 10 - KK  - allow choosing minimizer
C 12 Dec 08 - MWS - add a function to return number of QM/MM atoms
C 16 Oct 06 - MWS - check some commons for sanity with GAMESS
C 29 SEP 98 - CHC - Integrating DLC with TINKER
C 18 MAY 98 - JRS - TOYS: CALL TINKER OPT OF REGION 4, FORM QMMM GRAD
C
C
c  ##################################################################
c  #                                                                #
c  #     TOYS: DRIVER FOR QMMM OPTIMIZATION                         #
c  #                                                                #
c  ##################################################################
      SUBROUTINE TOYS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c
c --From GAMESS
      PARAMETER (MXATM=2000)
C
C --From Tinker
      INCLUDE 'sizes.i'
c
C --For use in this routine
      PARAMETER (ANG2BHR=1.889725989D+00,
     *           BHR2ANG=0.52917724924D+00,
     *           CMA2HB=8.432969D-4,zero=0.0D+00)
C
C -- MAXIMUM SIZES and COMMONS FOR QMMM LINKING
c
      LOGICAL MMONLY,QMMM,IMOMM,SIMOMM
      CHARACTER*6 CTMODE, CTMETH
      DOUBLE PRECISION GRDMIN
c
      PARAMETER (MAXLNK=100, MAXR1=2000)
c
      COMMON /QMMM1/ IMOMM,SIMOMM,NPAIR,NSEQ
      COMMON /QMMM2/ IQMATM(MAXR1),ibasfmo(MAXR1)
      COMMON /QMMM3/ LINKge(3*MAXLNK),blqm(MAXLNK),blmm(MAXLNK)
ckk -minimize- begin
      logical freezmm,frzmm
      double precision MINMET,mintmp
      COMMON /NWTOPT/ OPTPRG,MINMET,GRDMIN,maxhess,frzmm,CTMODE,CTMETH
      character*3 mzmeth
      equivalence (mzmeth,mintmp)
ckk -minimize- end
      COMMON /TINOPT/ mparti,MMONLY,QMMM
c
c -- New COMMON between TOYS and Tinker only
      COMMON /TGRAD/ TEG(3*MAXATM)
c
      DIMENSION TTEG(3*Maxatm), CTEMP(3,MAXLNK)
c
C -- Bits from GAMESS
c
      LOGICAL GOPARR,DSKWRK,MASWRK,LINEAR
c
      COMMON /FUNCT / ENERGY,EG(3*MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,nsymc,LINEAR
      COMMON /DLCFRZ/ FVALUE(50),ITABLE(50),IFTYPE(50),NCONST
      COMMON /FMCOM / XX(1)
c
c -- Tinker COMMONS
c
c     x       current x-coordinate for each atom in the system
c     y       current y-coordinate for each atom in the system
c     z       current z-coordinate for each atom in the system
c     n       total number of atoms in the current system
c     type    atom type number for each atom in the system
c
c
      integer n,type
      real*8 x,y,z
      common /atomst/ 
     #        x(maxatm),y(maxatm),z(maxatm),n,type(maxatm)
c
c ---------------------------------------------
c
ckk -minimize- 
      mintmp=minmet
c
C ------------------------------------------------------------
c  
C      Region 2-3 linking scheme
C      IMOMM by Morokuma et al.
c
C ------------------------------------------------------------
        IF(IMOMM) THEN
C    Region 1
          DO 10 I=1,NSEQ
             X(iqmatm(I))= C(1,I)*BHR2ANG
             Y(iqmatm(I))= C(2,I)*BHR2ANG
             Z(iqmatm(I))= C(3,I)*BHR2ANG 

   10     continue
c
c -- change Region 1-2 bondlengths to Region 1-3 bondlengths
c
          DO 20 I=1,NPAIR
             LR1=LINKge(3*I-2)
             LR2=LINKge(3*I-1)
             LR3=LINKge(3*I  )
             EXTEND=blmm(I)/blqm(I)
C
             DX=(C(1,LR2)-C(1,LR1))*BHR2ANG
             DY=(C(2,LR2)-C(2,LR1))*BHR2ANG
             DZ=(C(3,LR2)-C(3,LR1))*BHR2ANG
C        
             X(LR3)=C(1,LR1)*BHR2ANG +(DX*EXTEND)
             Y(LR3)=C(2,LR1)*BHR2ANG +(DY*EXTEND)
             Z(LR3)=C(3,LR1)*BHR2ANG +(DZ*EXTEND)
  20      CONTINUE
c
c      
C
C -- Call Tinker newton optimization with R1 and R3 INACTIVE
C    Convert Tinker gradient from kcal/mole/ang to Hartrees/Bohr
C
ckk -minimize- begin
ckk          CALL TNEWTX(CTMODE,CTMETH,GRDMIN)
          CALL TNEWTX(OPTPRG,CTMODE,CTMETH,mzmeth,GRDMIN,maxhess)
ckk -minimize- begin
C
          DO 30 I=1,NSEQ
            TTEG(I*3-2)=TEG(iqmatm(I)*3-2)*CMA2HB
            TTEG(I*3-1)=TEG(iqmatm(I)*3-1)*CMA2HB
            TTEG(I*3  )=TEG(iqmatm(I)*3  )*CMA2HB
  30      CONTINUE
          DO 35 I=1,npair
            TTEG(3*(I+nseq)-2)=TEG(linkge(I*3)*3-2)*CMA2HB
            TTEG(3*(I+nseq)-1)=TEG(linkge(I*3)*3-1)*CMA2HB
            TTEG(3*(I+nseq)  )=TEG(linkge(I*3)*3  )*CMA2HB
  35      CONTINUE

C
C -- Convert the GAMESS Cartesian gradient to internals
c
          CALL TRANG(EG,NVAR,3*NAT)
C    CHC
C    There may be other ways to check DLC coordinates.
C
      if (nzvar .gt. nvar) then
C...Dynamic memory allocation.......
       CALL VALFM(LOADFM)
       LCOR   = 1 + LOADFM
       LTCOR  = LCOR + NZVAR
       LS     = LTCOR + NZVAR
       LAST   = LS + NZvar*NZVAR
       NEED   = LAST - LCOR
       CALL GETFM(NEED)
C
       CALL DAREAD(IDAF,IODA,XX(LS),NZVAR*NZVAR,46,0)
       call tfdq(EG,XX(lcor),XX(ls),nvar,nzvar,1)
      endif
c
c -- Swap Region 3 and Region 2 positions in C
c    and recalculate B and B inverse before converting Tinker
c    grad to internals
c
          DO 40 I=1,npair
            CTEMP(1,I)=C(1,linkge(I*3-1))
            CTEMP(2,I)=C(2,linkge(I*3-1))
            CTEMP(3,I)=C(3,linkge(I*3-1))
C
            C(1,linkge(I*3-1))=X(linkge(I*3))*ANG2BHR
            C(2,linkge(I*3-1))=Y(linkge(I*3))*ANG2BHR
            C(3,linkge(I*3-1))=Z(linkge(I*3))*ANG2BHR
  40      CONTINUE
C
          CALL BANDBI
          CALL TRANG(TTEG,NVAR,3*NAT)

C...Transform the DLC steps, EG into Redunadant
C...internal coordinate steps, XX(LCOR) 
       call tfdq(TTEG,XX(lTcor),XX(ls),nvar,nzvar,1)

C
C -- QMMM gradient is the sum of QM and MM gradients in internals
C
cjrs temp diagnostic
c
      if (maswrk) write(iw,970)
        DO 260 I=1,NzVAR
          if (maswrk) WRITE(IW,971) I,xx(LCOR+I-1),XX(LTCOR+I-1),
     *    xx(LCOR+I-1)+XX(LTCOR+i-1)
 260    continue
c
        DO 55 I=1,NzVAR
          xx(LCOR+I-1)=XX(LCOR+I-1)+XX(LTCOR+I-1)
  55    CONTINUE
C
C -- Swap back Region 2 positions for Region 3  in C
C    Recalculate B and B inverse
C    Convert QMMM gradient back to Cartesians
C
          DO 60 I=1,Npair
            C(1,linkge(I*3-1))=CTEMP(1,I)
            C(2,linkge(I*3-1))=CTEMP(2,I)
            C(3,linkge(I*3-1))=CTEMP(3,I)
 60       CONTINUE   
c
          CALL BANDBI
      call tfdq(XX(lcor),EG,XX(ls),nvar,nzvar,2)
          CALL TRANGB(EG,NVAR,3*NAT)
          CALL RETFM(NEED)
          GOTO 600
C
C      End of {IF(IMOMM)} block
C
       END IF   
c ---------------------------------------------------------
C
c        SIMOMM: No explicit link atoms used
C
c ---------------------------------------------------------     
         IF(SIMOMM) THEN
          DO 70 I=1,nseq
          X(iqmatm(I))= C(1,I)*BHR2ANG
          Y(iqmatm(I))= C(2,i)*BHR2ANG
  70      Z(iqmatm(I))= C(3,i)*BHR2ANG
c
          do 75 i=1,3*nat
           tteg(i)=zero
  75      continue
C
C -- Call Tinker newton optimization with R1 and R3 INACTIVE
C    Convert Tinker gradient from kcal/mole/ang to Hartrees/Bohr
C
ckk -minimize- begin
ckk          CALL TNEWTX(CTMODE,CTMETH,GRDMIN)
          CALL TNEWTX(OPTPRG,CTMODE,CTMETH,MINMET,GRDMIN,maxhess)
ckk -minimize- begin
          DO 80 I=1,nseq
            tTEG(I*3-2)= TEG(iqmatm(I)*3-2)*CMA2HB
            tTEG(I*3-1)= TEG(iqmatm(I)*3-1)*CMA2HB
            tTEG(I*3  )= TEG(iqmatm(I)*3  )*CMA2HB
  80      CONTINUE
          CALL TRANG(EG,NVAR,3*NAT)
          CALL TRANG(tTEG,NVAR,3*NAT)
C
      if (maswrk) then
      write(iw,975)
        DO 275 I=1,nvar
           WRITE(IW,971) I,EG(I),tTEG(I), 
     *                     eg(i)+tteg(i)
 275    continue
      endif
C
C -- QMMM gradient is the sum of these two
C
           DO 90 I=1,nvar
              EG(i)=EG(i)+tTEG(I)
  90       CONTINUE
           CALL TRANGB(EG,NVAR,3*NAT)
           GOTO 600
C
C      IF(SIMOMM) ENDIF
C
        END IF  
c
c
c
 600  CONTINUE 
c -- This bit should go at the very end after the gradient conversion
c     and values in C have been shifted back to GAMESS values
c
c
      RETURN
  970 FORMAT(/,/,5X,'INTERNAL',8X,'GAMESS GRAD',7X,'Tinker GRAD',6X,
     #'HYBRID GRAD',/, 4X,'COORDINATE',8X,'(H/B, H/R)',6X,'(H/B, H/R)',
     #9X,'(H/B, H/R)',/,4X,63('-'))
  971 FORMAT(6X,I3,10X,F12.8,5X,F12.8,5X,F12.8)
  975 FORMAT(/,/,5X,'INTERNAL',8X,'GAMESS GRAD',7X,'Tinker GRAD',6X,
     #'HYBRID GRAD',/, 4X,'COORDINATE',8X,'(H/B, H/R)',6X,'(H/B, H/R)',
     #9X,'(H/B, H/R)',/,4X,63('-'))
      END
c
      integer function nqmmmatoms()
      INCLUDE 'sizes.i'
      integer n,type
      real*8 x,y,z
      common /atomst/ 
     #        x(maxatm),y(maxatm),z(maxatm),n,type(maxatm)
c
c       ----- return total number of atoms in QM/MM run -----
c       the first ones are the atoms in the QM region, except
c       capping hydrogens are -not- included in this total.
c
      nqmmmatoms = n
      return
      end
C
C
c  ##################################################################
c  #                                                                #
c  #     TOYSfmo: DRIVER FOR QMMM OPTIMIZATION                         #
c  #                                                                #
c  ##################################################################
      SUBROUTINE TOYSfmo(natfmo,fmoc,egmm,emm)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension fmoc(3,*),egmm(*)
c
casa -restrain- begin
      include 'sizes.i'
      include 'restrn.i'
      include 'potent.i'
casa -restrain- end
c
c --From GAMESS
      PARAMETER (MXATM=2000)
C
C --From Tinker
ckk      INCLUDE 'sizes.i'
casa      PARAMETER (MAXATM=12000)
c
casa -restrain- 
C --For use in this subroutine
      PARAMETER (ANG2BHR=1.889725989D+00,
     *           BHR2ANG=0.52917724924D+00,
     *           CMA2HB=8.432969D-4,zero=0.0D+00)
      parameter (tokcal=627.51D+00)
C
C -- MAXIMUM SIZES and COMMONS FOR QMMM LINKING
c
      LOGICAL MMONLY,QMMM,IMOMM,SIMOMM
      CHARACTER*6 CTMODE, CTMETH
      DOUBLE PRECISION  GRDMIN
c
      PARAMETER (MAXLNK=100, MAXR1=2000)
c
      COMMON /QMMM1/ IMOMM,SIMOMM,NPAIR,NSEQ
      COMMON /QMMM2/ IQMATM(MAXR1),ibasfmo(MAXR1)
ccc   COMMON /QMMM3/ LINKge(3*MAXLNK),blqm(MAXLNK),blmm(MAXLNK)
ckk -minimize- begin
      logical frzmm
      double precision MINMET, mintmp
      COMMON /NWTOPT/ OPTPRG,MINMET,GRDMIN,maxhess,frzmm,CTMODE,CTMETH
      character*3 mzmeth
      equivalence (mzmeth,mintmp)
ckk -minimize- end
      COMMON /TINOPT/ mparti,MMONLY,QMMM
c
c -- New COMMON between TOYS and Tinker only
      COMMON /TGRAD/ TEG(3*MAXATM)
c
ckk -fmoimomm-      DIMENSION TTEG(3*Maxatm), CTEMP(3,MAXLNK)
c
C -- Bits from GAMESS
c
      LOGICAL GOPARR,DSKWRK,MASWRK,LINEAR
      integer ME,MASTER,NPROC,IBTYP,IPTIM
c
      COMMON /FUNCT / ENERGY,EG(3*MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
ckk -fmoimomm-      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,nsymc,LINEAR
ckk -fmoimomm-      COMMON /DLCFRZ/ FVALUE(50),ITABLE(50),IFTYPE(50),NCONST
ckk -fmoimomm-      COMMON /FMCOM / XX(1)
c
c -- Tinker COMMONS
c
c     x       current x-coordinate for each atom in the system
c     y       current y-coordinate for each atom in the system
c     z       current z-coordinate for each atom in the system
c     n       total number of atoms in the current system
c     type    atom type number for each atom in the system
c
c
      integer n,type
      real*8 x,y,z
      common /atomst/ 
     #        x(maxatm),y(maxatm),z(maxatm),n,type(maxatm)
c
C ------------------------------------------------------------
        IF(qmmm .and. IMOMM) THEN
          if(maswrk) write(iw,980)
          call abrt
        endif
c ---------------------------------------------------------
C
c        SIMOMM: No explicit link atoms used
C
c ---------------------------------------------------------     
ckk   -minimize-
      mintmp=minmet
ckk -fmoimomm-         IF(SIMOMM) THEN
      if(qmmm) then
          nattmp=nseq
          use_geom = .true.
          do i=1,nseq
             X(iqmatm(I))= fmoc(1,I)*BHR2ANG
             Y(iqmatm(I))= fmoc(2,i)*BHR2ANG
             Z(iqmatm(I))= fmoc(3,i)*BHR2ANG
casa -restrain- begin
             if( .not. frzmm) then
             npfix = npfix + 1
             ipfix(npfix) = iqmatm(i)
             xpfix(npfix) = X(iqmatm(i))
             ypfix(npfix) = Y(iqmatm(i))
             zpfix(npfix) = Z(iqmatm(i))
             pfix(npfix) =  RESTRAINV
             endif
casa -restrain- end
          enddo
c
          do i=1,3*nseq
             egmm(i)=zero
          enddo
C
C -- Call Tinker newton optimization with R1 and R3 INACTIVE
C    Convert Tinker gradient from kcal/mole/ang to Hartrees/Bohr
C
          if(frzmm) then
             call analyze
          else
ckk -minimize- begin
c             CALL TNEWTX(CTMODE,CTMETH,GRDMIN)
             CALL TNEWTX(OPTPRG,CTMODE,CTMETH,mzmeth,GRDMIN,maxhess)
ckk -minimize- end
          endif 
c
          DO I=1,nseq
            egmm(i*3-2)= TEG(iqmatm(I)*3-2)*CMA2HB
            egmm(i*3-1)= TEG(iqmatm(I)*3-1)*CMA2HB
            egmm(i*3  )= TEG(iqmatm(I)*3  )*CMA2HB
          ENDDO
          emm=ENERGY/tokcal
      endif
C
      if(mmonly) then
         nattmp=natfmo
          do i=1,natfmo
             X(i)= fmoc(1,i)*BHR2ANG
             Y(i)= fmoc(2,i)*BHR2ANG
             Z(i)= fmoc(3,i)*BHR2ANG
          enddo
          call analyze
          emm=ENERGY/tokcal
          do i=1,3*natfmo
            egmm(i)= TEG(i)*CMA2HB
          enddo
      endif
c
      if (maswrk) then
c           write(iw,*) ' MM energy emm ',emm
        write(iw,975)
        DO I=1,nattmp
           WRITE(iw,971) I,egmm(I*3-2),egmm(i*3-1),egmm(i*3) 
        ENDDO
      endif
C
C -- QMMM gradient is the sum of these two
C
ckk -fmoimomm- begin
ckk        DO 90 I=1,nvar
ckk           EG(i)=EG(i)+tTEG(I)
ckk  90       CONTINUE
ckk -fmoimomm- end
ckk -fmoimomm-           CALL TRANGB(EG,NVAR,3*NAT)
C
ckk -fmoimomm-        END IF  
c
c
c
c -- This bit should go at the very end after the gradient conversion
c     and values in C have been shifted back to GAMESS values
c
c
      RETURN
ckk -fmoimomm- begin
ckk  970 FORMAT(/,/,5X,'INTERNAL',8X,'GAMESS GRAD',7X,'Tinker GRAD',6X,
ckk     #'HYBRID GRAD',/, 4X,'COORDINATE',8X,'(H/B, H/R)',6X,'(H/B, H/R)',
ckk     #9X,'(H/B, H/R)',/,4X,63('-'))
ckk -fmoimomm- end
  971 FORMAT(6X,I4,10X,F12.8,5X,F12.8,5X,F12.8)
  975 FORMAT(/,/,5X,'Tinker GRAD (H/B)',
     .                /,4X,63('-'))
  980 FORMAT(6X,' IMOMM is not allowed in fmo. Use SIMOMM')
      END
