module asym_cstes

  implicit none
  public
  private :: sisnan, disnan

  integer, parameter :: pr = selected_real_kind( p = 14 )
  !
  !! Flag to known if the code is used as a library or not
  !
  logical :: from_lib = .true.
  !
  !! Constantes
  !
  real (pr), parameter :: hbar_c =197.326968_pr ! [MeVfm]
  real (pr), parameter :: hbarc2 = hbar_c * hbar_c
  real (pr), parameter :: hbarc3 = hbarc2 * hbar_c
  real (pr), parameter :: hbarc4 = hbarc2 * hbarc2
  real (pr), parameter :: hbarc5 = hbarc3 * hbarc2
  real (pr), parameter :: hbarc6 = hbarc3 * hbarc3
  real (pr), parameter :: hbarc7 = hbarc4 * hbarc3
  real (pr), parameter :: hbarc8 = hbarc4 * hbarc4
  real (pr), parameter :: hbarc9 = hbarc5 * hbarc4
  real (pr), parameter :: hbarc10 = hbarc5 * hbarc5

  !
  !!
  !
  real (pr) :: tpi_2, qpi2, six_pi2, tqpi2, tpi13, tpi2_13, zeta3,zeta5
  real (pr) :: tp2hbc2,pi2,pi,tpi23,t23,dpi2,qpi
  real (pr) :: t13,t35,t53
  !
  !!
  !
  real (pr), parameter :: qmax = 4.0_pr   ! fm^-1 Maximu value of q in the plots
  real (pr), parameter :: dqq  = 0.01_pr ! fm^-1 step in q for the plots
  !
  !!
  !
  integer, parameter :: npt=96 ! Maximal number of GH points
  logical, parameter :: flagBnpanalytic=.true.
  real (pr)          :: glp(npt),glw(npt)
  complex (pr)       :: cz


  interface is_nan
     module procedure sisnan
     module procedure disnan
  end interface

contains

  !
  !!
  !
  
  subroutine set_lr_cstes()
    implicit none
    !
    t13=1.0_pr /3.0_pr
    t23=2.0_pr /3.0_pr
    t35=3.0_pr /5.0_pr
    t53=5.0_pr /3.0_pr

    pi = 4.0_pr * atan(1.0_pr)
    zeta3=1.202056903159594285
    zeta5=1.03692775514337
    qpi = 4 * pi
    cz = cmplx( 0.0_pr, 0.0_pr , kind=pr)
    pi2=pi*pi
    tpi_2   =  3 * pi2
    dpi2    =  2 * pi2
    qpi2    =  4 * pi2
    six_pi2 =  6 * pi2
    tqpi2   = 12 * pi2
    tpi13   = ( 3.0_pr * pi2 )**t13
    tpi2_13 = ( 1.5_pr * pi2 )**t13
    tpi23 = ( 1.5_pr * pi2 )**t23
    tp2hbc2 = tpi_2 * hbarc2
    !
    ! we calculate at the beginning the number of G.L.
    ! points and weigths to save time
    call LEGZO(npt,glp,glw)

  end subroutine set_lr_cstes

  !
  !!
  !

  function sisnan(x) result(res)
    implicit none
    real ( selected_real_kind( p = 6 ) ) :: x
    logical :: res
    !
    res = ( x /= x )
    !
  end function sisnan

  !
  !!
  !
  function F0m(m,Itau,Isigma,Istau) result(res)
    implicit none
    real (kind=pr), intent (in) :: m,Itau,Isigma,Istau
    real (kind=pr) :: s1,s2,s3,s4,res

    s1=(1+Itau+Isigma+Istau)**m
    s2=(1+Itau-Isigma-Istau)**m
    s3=(1-Itau+Isigma-Istau)**m
    s4=(1-Itau-Isigma+Istau)**m

    res=(s1+s2+s3+s4)/4.0_pr

    return
  end function F0m

  function F1m(m,Itau,Isigma,Istau) result(res)
    implicit none
    real (kind=pr), intent (in) :: m,Itau,Isigma,Istau
    real (kind=pr) :: s1,s2,s3,s4,res

    s1=(1.0_pr+Itau+Isigma+Istau)**m
    s2=(1.0_pr+Itau-Isigma-Istau)**m
    s3=(1.0_pr-Itau+Isigma-Istau)**m
    s4=(1.0_pr-Itau-Isigma+Istau)**m

    res=(s1+s2-s3-s4)/4.0_pr

    return
  end function F1m
  
  !
  !!
  !
  function G0m(m,Itau,Isigma,Istau) result(res)
    implicit none
    real (kind=pr), intent (in) :: m,Itau,Isigma,Istau
    real (kind=pr) :: s1,s2,s3,s4,res

    s1=(1+Itau+Isigma+Istau)**m
    s2=(1+Itau-Isigma-Istau)**m
    s3=(1-Itau+Isigma-Istau)**m
    s4=(1-Itau-Isigma+Istau)**m

    res=(s1-s2+s3-s4)/4.0_pr

    return
  end function G0m
  
  !
  !!
  !
  function G1m(m,Itau,Isigma,Istau) result(res)
    implicit none
    real (kind=pr), intent (in) :: m,Itau,Isigma,Istau
    real (kind=pr) :: s1,s2,s3,s4,res

    s1=(1+Itau+Isigma+Istau)**m
    s2=(1+Itau-Isigma-Istau)**m
    s3=(1-Itau+Isigma-Istau)**m
    s4=(1-Itau-Isigma+Istau)**m

    res=(s1-s2-s3+s4)/4.0_pr

    return
  end function G1m
  !
  !!
  !

  function disnan(x) result(res)
    implicit none
    real ( selected_real_kind( p = 12 ) ) :: x
    logical :: res
    !
    res = ( x /= x )
    !
  end function disnan

  !
  !!
  !

  function trl(str) result(tstr)
    implicit none
    character ( len = * ) :: str
    character ( len = len(trim(adjustl(str))) ) :: tstr
    !
    tstr = trim(adjustl(str))
    !
  end function trl

  !
  !!
  !

  function trluc(s) result(res)
    implicit none
    character ( len = * ) :: s
    character ( len = len(trim(adjustl(s))) ) :: res
    integer :: i, is
    integer, parameter :: sh = iachar( "a" ) - iachar( "A" )
    integer, parameter :: smin = iachar( "a" )
    integer, parameter :: smax = iachar( "z" )
    !
    res = trim(adjustl(s))
    do i = 1, len(res)
       is = iachar(res(i:i))
       if ( is < smin .or. is > smax  ) cycle
       res(i:i) = char( is - sh )
    end do
    !
  end function trluc


!
!!
!

          SUBROUTINE LEGZO(N,X,W)
!C
!C       =========================================================
!C       Purpose : Compute the zeros of Legendre polynomial Pn(x)
!C                 in the interval [-1,1], and the corresponding
!C                 weighting coefficients for Gauss-Legendre
!C                 integration
!C       Input :   n    --- Order of the Legendre polynomial
!C       Output:   X(n) --- Zeros of the Legendre polynomial
!C                 W(n) --- Corresponding weighting coefficients
!C       =========================================================
!C
        IMPLICIT NONE
        INTEGER          :: N,N0,NR,I,K,J
        DOUBLE PRECISION :: X(N),W(N),Z,Z0,P,F0,F1
        DOUBLE PRECISION :: PF,PD,Q,WP,GD,FD
        N0=(N+1)/2
        DO 45 NR=1,N0
           Z=DCOS(3.1415926D0*(NR-0.25D0)/N)
10         Z0=Z
           P=1.0D0
           DO 15 I=1,NR-1
15            P=P*(Z-X(I))
           F0=1.0D0
           IF (NR.EQ.N0.AND.N.NE.2*INT(N/2)) Z=0.0D0
           F1=Z
           DO 20 K=2,N
              PF=(2.0D0-1.0D0/K)*Z*F1-(1.0D0-1.0D0/K)*F0
              PD=K*(F1-Z*PF)/(1.0D0-Z*Z)
              F0=F1
20            F1=PF
           IF (Z.EQ.0.0) GO TO 40
           FD=PF/P
           Q=0.0D0
           DO 35 I=1,NR
              WP=1.0D0
              DO 30 J=1,NR
                 IF (J.NE.I) WP=WP*(Z-X(J))
30            CONTINUE
35            Q=Q+WP
           GD=(PD-Q*FD)/P
           Z=Z-FD/GD
           IF (DABS(Z-Z0).GT.DABS(Z)*1.0D-15) GO TO 10
40         X(NR)=Z
           X(N+1-NR)=-Z
           W(NR)=2.0D0/((1.0D0-Z*Z)*PD*PD)
45         W(N+1-NR)=W(NR)
        RETURN
        END SUBROUTINE
!---------------------

        Function LEG(N,X)
!C
!C       ===============================================
!C       Purpose: Compute Legendre polynomials Pn(x)
!C                and their derivatives Pn'(x)
!C       Input :  x --- Argument of Pn(x)
!C                n --- Degree of Pn(x) ( n = 0,1,...)
!C       Output:  PN(n) --- Pn(x)
!C                PD(n) --- Pn'(x)
!C       ===============================================
        IMPLICIT NONE

        INTEGER          :: N,k
        DOUBLE PRECISION :: X,P0,P1,PF
        Double precision :: PN(0:max(N,2)),PD(0:max(N,2)),LEG

        if(abs(x).gt.1)then
          write(*,*)'qualcosa di wrong nella Function Leg'
          write(*,*)x
          stop
        endif
        PN(0)=1.0D0
        PN(1)=X
        PD(0)=0.0D0
        PD(1)=1.0D0
        P0=1.0D0
        P1=X
        DO 10 K=2,N
           PF=(2.0D0*K-1.0D0)/K*X*P1-(K-1.0D0)/K*P0
           PN(K)=PF
           IF (DABS(X).EQ.1.0D0) THEN
              PD(K)=0.5D0*X**(K+1)*K*(K+1.0D0)
           ELSE
              PD(K)=K*(P1-X*PF)/(1.0D0-X*X)
           ENDIF
           P0=P1
10         P1=PF
        Leg=PN(N)
        RETURN
        END FUNCTION

!
!!
!


end module asym_cstes

