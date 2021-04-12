module inm
  use asym_cstes ! module constants
  use asym_param ! module c.c. constants of the functional

contains
  !
  !!
  !
  subroutine inm_properties(rho_sat,nom)
    implicit none

    character ( len = * ),intent(in):: nom
    integer            :: id
    real ( kind = pr ) :: rho_sat,kappa,EA
    real ( kind = pr ) :: mn,mp,rap,rap1,rap2
    character(len=50) :: filename

    write( filename, '("INM_properties_",A,".dat")' ) trluc(nom)

    call open_file( id, filename )

    ! saturation density
    call saturation_density(rho_sat)
    ! incompressibility at rho_0
    call incompressibility(rho_sat,kappa)
    ! E/A at rho0
    call EnergyPerParticle(rho_sat/2,rho_sat/2,EA)
    ! effective mass

    call mass_skyrme(rho_sat/2,rho_sat/2,mn,mp)

    rap=1.0_pr/((hbarc2/(2.0_pr*mn))/hb2m)

    write(*,*)'=============================='
    write(*,*)' INM properties               '
    write(*,*)
    write(*,*)' symmetric matter             '
    write(*,*)'rho_0= ',rho_sat
    write(*,*)'K_inf= ',kappa
    write(*,*)'E/A =  ',EA
    write(*,*)'m*/m =',rap
    write(*,*)
    write(*,*)' pure neutron matter'

    call mass_skyrme(rho_sat,0.0_pr,mn,mp)
    rap1=1.0_pr/((hbarc2/(2.0_pr*mn))/hb2m)
    rap2=1.0_pr/((hbarc2/(2.0_pr*mp))/hb2m)

    write(*,*)
    write(*,*)'m*/m(neutr)= ',rap1
    write(*,*)'m*/m(proto)= ',rap2
    write(*,*)'Delta m*/m = ',rap1-rap2
    write(*,*)
    write(*,*)'=============================='


    !
    write(id,*)'=============================='
    write(id,*)' INM properties  - ',trluc(nom)
    write(id,*)
    write(id,*)' symmetric matter '
    write(id,*)'rho_0= ',rho_sat
    write(id,*)'K_inf= ',kappa
    write(id,*)'E/A =  ',EA
    write(id,*)'m*/m =',rap
    write(id,*)
    write(id,*)' pure neutron matter'
    write(id,*)
    write(id,*)'m*/m(neutr)= ',rap1
    write(id,*)'m*/m(proto)= ',rap2
    write(id,*)'Delta m*/m = ',rap1-rap2
    write(id,*)
    write(id,*)'=============================='


    call close_file(id)


    return
  end subroutine inm_properties



  function pres(rho ) result(res)
    implicit none

    real ( kind = pr ) :: rho, res, t364,t25,t85,t83

    t25=2.0_pr/5.0_pr  
    t364=3.0_pr/64.0_pr
    t85 = 8.0_pr / 5
    t83 = 8.0_pr / 3
    !
    res = t25 * hb2m * tpi23 * rho**t23      &
         + c_rho(0,1) * rho                  &
         + tpi23  * c_tau(0) * rho**t53      &
         + c_rho(0,2) * s11 * rho**s11       &
         + c_rho(0,3) * s21 * rho**s21      &
         + 2 * Brho_0 * rho**2               &
         + t85 * tpi23 * Btau_0 * rho**t83  &
         +  3 * t364*v0 * rho**3
    !
    res = res * rho
    !
  end function pres


  subroutine saturation_density(rho_sat)
    !----------------------------------------------------
    ! Abstract: routine that calculates the saturation
    !           density in SNM for a given force
    !
    ! input: through modules
    ! output: rho_sat
    !----------------------------------------------------
    implicit none
    real (pr) :: rho_sat
    real (pr) :: r0, r1, r, p0, p1, p,eps
    real (pr), parameter :: dr = 0.01_pr
    integer, parameter :: max_iter = 250
    integer :: iter

    eps = 10 * spacing(1.0_pr)

    !
    rho_sat = -1.0_pr
    !
    r0 = 0.16_pr - eps
    do
       p0 = pres( r0 )
       if ( p0 < 0.0_pr ) exit
       r0 = r0 - dr
       if ( r0 <= - dr  ) return
       if ( r0 < 0.0_pr ) r0 = 0.0_pr
    end do
    r1 = 0.16_pr + eps
    do
       p1 = pres( r1 )
       if ( p1 > 0.0_pr ) exit
       r1 = r1 + dr
       if ( r1 > 1.6_pr ) return
    end do
    iter = 0
    do
       iter = iter + 1
       r = ( r0 + r1 ) / 2
       p = pres(  r )
       if ( p0 * p < 0.0_pr ) then
          p1 = p
          r1 = r
       else
          p0 = p
          r0 = r
       end if
       if ( abs( r0 - r1 ) < 10 * eps ) exit
       if ( iter > max_iter ) return
    end do
    !
    rho_sat = r
    !
  end subroutine saturation_density


  !
  !!
  !
  subroutine EnergyPerParticle(rhoN,rhoP,EA)
    !
    ! routine that calculates the E/A
    ! in asymmetric non polarized matter
    ! 
    ! input
    ! rhoN,rhoP
    ! output: E/A
    !
    implicit none

    real ( kind = pr ), intent (in) :: rhoN,rhoP
    real ( kind = pr )              :: rho0,Ex,EA
    real ( kind = pr )              :: con,t23,t53,coeff
    real ( kind = pr )              :: rhosig,rhosig2

    rho0=rhoN+rhoP   
    Ex=(rhoN-rhoP)/rho0 
    t23=2.0_pr/3
    t53=5.0_pr/3
    con=3.0_pr/5 * tpi23*rho0**t23

    rhosig =rho0**sigma
    rhosig2=rho0**sigma2

    coeff=c_rho(0,1)+c_rho(0,2)*rhosig+c_rho(0,3)*rhosig2+c_tau(0)*con*F0m(t53,Ex,0.0_pr,0.0_pr) &
         +(c_rho(1,1)+c_rho(1,2)*rhosig+c_rho(1,3)*rhosig2)*Ex**2 &
         +c_tau(1)*con*F1m(t53,Ex,0.0_pr,0.0_pr)*Ex

    EA=con*hb2m*F0m(t53,Ex,0.0_pr,0.0_pr)+rho0*coeff

    return
  end subroutine EnergyPerParticle

  !
  !!
  !

  subroutine incompressibility(rho,Kappa)

    ! calculates the incompressibility 
    ! in the SNM case only
    !
    implicit none


    real ( kind = pr ), intent (in) :: rho
    real ( kind = pr )              :: Kappa


    Kappa=18.0_pr*pres(rho)/rho+9.0_pr*rho**2*der2EA(rho) 

    return
  end subroutine incompressibility


  !
  !!
  !

  function der2EA(rho ) result(res)
    implicit none

    real ( kind = pr ) :: rho,a1,a2,a3,t23,t53,res

    t23=2.0_pr/3
    t53=50_pr/3


    a1=-2.0_pr/15*hb2m*tpi23*rho**t23*F0m(t53,0.0_pr,0.0_pr,0.0_pr)
    a2=sigma*c_rho(0,2)*rho**(sigma-1)+sigma2*c_rho(0,3)*rho**(sigma2-1) &
         +6.0_pr/15*c_tau(0)*tpi23*rho**(-1.0_pr/3)*F0m(t53,0.0_pr,0.0_pr,0.0_pr)
    a2=2.0_pr*a2

    a3=sigma*(sigma-1)*c_rho(0,2)*rho**(sigma-2)    &
         +sigma2*(sigma2-1)*c_rho(0,3)*rho**(sigma2-2) &
         -2.0_pr/15*c_tau(0)*tpi23*rho**(-4.0_pr/3)*F0m(t53,0.0_pr,0.0_pr,0.0_pr)
    a3=a3*rho
    !
    res = a1+a2+a3
    !
  end function der2EA


  !
  !!
  !

  subroutine mass_skyrme(rhoN,rhoP,mstarN,mstarP)
    !=========================================================================
    ! abstrac: routine to calculate the effective mass 
    !
    ! input
    ! rhoN,rhoP= density of neutron and protons
    !
    ! output
    ! mstarN(P)  [Mev^2]
    !========================================================================
    implicit none
    real ( kind = pr ), INTENT(in) :: rhoN,rhoP
    real ( kind = pr ) :: hbmstarN,hbmstarP
    real ( kind = pr ) :: mstarN,mstarP
    real ( kind = pr ) :: Ex,keff(0:1),rho0
    !
    ! 0: neutron
    ! 1: protons
    !
    rho0=rhoN+rhoP
    Ex=(rhoN-rhoP)/rho0
    !
    keff(0)=C_tau(0)*rho0+Ex*C_tau(1)*rho0        & ! 2 body
         +(Btau_0+Ex*Btau_1+Ex**2*Btau_10)*rho0**2  ! 3body
    keff(1)=C_tau(0)*rho0-Ex*C_tau(1)*rho0       & ! 2 body
         +(Btau_0-Ex*Btau_1+Ex**2*Btau_10)*rho0**2  ! 3body
    !
    hbmstarN=hb2m+keff(0)
    hbmstarP=hb2m+keff(1)

    mstarN = hbarc2 / ( 2 * hbmstarN )
    mstarP = hbarc2 / ( 2 * hbmstarP )


    !
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    ! Nota: in Eq.24 of HFBRAD there is an error, the term in t2 multiplys
    !       ALWAYS the terms x2 and NOT x1!!
    !--------------------------------------------------------------------------
    !
    return
  end subroutine mass_skyrme


end module inm
