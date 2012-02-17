c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1992  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  bath.i  --  temperature and pressure bath control values  ##
c     ##                                                            ##
c     ################################################################
c
c
c     kelvin      target value for the system temperature (K)
c     atmsph      target value for the system pressure (atm)
c     tautemp     time constant in psec for temperature bath coupling
c     taupres     time constant in psec for pressure bath coupling
c     compress    isothermal compressibility of medium (atm-1)
c     isothermal  logical flag geverning use of temperature bath
c     isobaric    logical flag governing use of pressure bath
c
c
      real*8 kelvin,atmsph,tautemp,taupres,compress
      logical isothermal,isobaric
      common /bath/ kelvin,atmsph,tautemp,taupres,compress,isothermal,
     &              isobaric
