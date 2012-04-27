C 24 Apr 12 - DGF/MWS - Fix ADD_ARG by adding arg, and fake print
C 12 APR 12 - STUB FILE
C
C STUB FILE FOR THE INTERFACE TO MPQC
C*MODULE EXCSTB     *DECK ADD_ARG
      SUBROUTINE ADD_ARG(arg)
      character*256 arg
      WRITE(6,'(/10x,"**********ERROR**********")')
      WRITE(6,'("[2]_R12 NOT ENABLED, PLEASE RE-COMPILE WITH MPQC"/)')
      CALL ABRT
C        next never executed, but keeps FTNCHEK happy.
      WRITE(6,8000) arg
 8000 format('arg=',A)
      RETURN
      END
C*MODULE EXCSTB     *DECK NOTMAIN
      SUBROUTINE PT2R12_MAIN()
      RETURN
      END
