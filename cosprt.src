C 23 Jun 10 - RP  - partial COSMO-RS PRINTOUT, w/o full COSMO-RS LICENSE
C
C   must not capitalize in this file so COSMO-RS can read its output!
C
C*MODULE COSPRT  *DECK COSPRT
      SUBROUTINE COSPRT()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C INCOMPLETE PRINTOUT ROUTINE FOR PUBLIC RELEASE, THE FULL COSPRT 
C ROUTINE, WITH .cosmo OUTPUT, IS PROVIDED TO PEOPLE WITH A VALID
C COSMO-RS LICENSE:
C      VISIT http://ocikbws.uzh.ch/gamess FOR MORE INFORMATION.
C
      PARAMETER (NPPA=1082)
      LOGICAL COSBUG,COSWRT,DCOSMO,prfcnd,ioutch
      COMMON /CONV  / DENTOL,EN,ETOT,EHF,EHF0,DIFF,ITER,ICALP,ICBET
      COMMON /COSDAT/ SE2,SECORR,QVCOSMO,ELAST,EMP2COS,EMP2LAST,
     *                COSVOL,COSSAR,EDIEL,EOC1,DEOC_RS,SUMQSC,
     *                SUMQSCOLD,ZSUM,ZSUM2,ZSUM3,FEPSI,RDS,DISEX2,
     *                EPSI,COSRAD,DISEX,OUTCHG,EDIEL_SAVE,
     *                maxnps,ICORR,ITRIP,NQS,MP2TRIP,MP2ITER,
     *                ICFREQ,NSPA,NSPH,NPSD,NPS,NPS2,NDEN,NPSPHER,
     *                COSBUG,COSWRT,DCOSMO,PRFCND,ioutch
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
      DATA DMULTI,DBLCAV/8HDMULTI  ,8HDBLCAV  /
C
      WRITE(IW,*)
      WRITE(IW,*)'====================================================='
      WRITE(IW,*)'  GAMESS/COSab RESULTS'
      WRITE(IW,*)'====================================================='
      WRITE(IW,*)
      IF(PRFCND) THEN
      WRITE(IW,*)'  Dielectric Constant                  =',
     &' INFTY'
      ELSE
      WRITE(IW,'(A,F8.2)')'  Dielectric Constant                  =',
     &EPSI
      END IF
      WRITE(IW,*)'  Basic Grid Size                      =',nppa
      WRITE(IW,*)'  Number of Segments                   =',nspa
      WRITE(IW,*)'  Number of Segments located on spheres=',npspher
      WRITE(IW,*)'  Number of Segments of the cavity     =',nps
      WRITE(IW,'(A,F5.2)')'   Solvent Radius (rsolv)               =',
     &cosrad
      IF(OUTCHG.EQ.DBLCAV) THEN
      WRITE(IW,*)'  Outlying Charge Error Method         =',
     &' double cavity'
      ELSE IF (OUTCHG.EQ.DMULTI) THEN
      WRITE(IW,*)'  Outlying Charge Error Method         =',
     &' distributed multipole'
      ENDIF
      WRITE(IW,*)
      WRITE(IW,*)'====================================================='
      WRITE(IW,*)
      IF(MPLEVL.GT.0) then
      WRITE(IW,'(a,f17.10)')'  Total solvated Energy (a.u.)          =',
     &emp2-eoc1
      ELSE
      WRITE(IW,'(a,f17.10)')'  Total solvated Energy (a.u.)          =',
     &etot-eoc1
      ENDIF
      IF(MPLEVL.GT.0) THEN
      WRITE(IW,'(a,f17.10)')'  Total Energy corrected (a.u.)         =',
     &emp2
      ELSE 
      WRITE(IW,'(a,f17.10)')'  Total Energy corrected (a.u.)         =',
     &etot
      ENDIF
      WRITE(IW,*)
      WRITE(IW,'(a,f17.10)')'  Outlying charge effect (a.u.)         =',
     &eoc1 
      WRITE(IW,'(a,f17.10)')'  Outlying charge effect (kcal/mol)     =',
     &eoc1*627.51D+00
      WRITE(IW,*)
      IF (DCOSMO) THEN
      WRITE(IW,'(a,a,f17.10)')'  Outlying charge corr. (DCOSMO-RS)',
     &' (a.u.)     =',deoc_rs
      WRITE(IW,'(a,a,f17.10)')'  Outlying charge corr. (DCOSMO-RS)',
     &' (kcal/mol) =',deoc_rs*627.51D+00
      WRITE(IW,*)
      ENDIF
      WRITE(IW,'(a,f17.10)')'  Dielectric Energy (a.u.)              =',
     &secorr
      WRITE(IW,'(a,f17.10)')'  Dielectric energy corr. (a.u.)        =',
     &secorr+eoc1
      WRITE(IW,*)
      WRITE(IW,'(a,f17.10)')'  Sum of polarization charges           =',
     &zsum2 
      WRITE(IW,'(a,f17.10)')'  Sum of polarization charges(cor.)     =',
     &zsum
      WRITE(IW,*)
      WRITE(IW,'(a,f8.2)')'  Total surface area of cavity (A**2)   =',
     &cossar
      WRITE(IW,'(a,f8.2)')'  Total volume of cavity (A**3)         =',
     &cosvol
      WRITE(IW,*)
      WRITE(IW,*)'====================================================='
      WRITE(IW,*)
      RETURN
      END
