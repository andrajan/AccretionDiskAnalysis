module asym_param

  use asym_cstes
  use asym_io
  implicit none
  public
  !
  !! Namelist for the input
  !
  !character ( len = 20 ) :: edf = "SLY4"
  !character ( len = 4  ) :: opt = "NONE"
  !character ( len = 3  ) :: mat = "SNM"
  
  logical :: from_rhosat=.false.
  !
  !namelist /lr_input/ edf, opt, mat, rpa, plot, smrul, crit, rho
  !
  !! skyrme_force
  !
  real ( kind = pr ) :: t0, t1, t2, t3, t4, td
  real ( kind = pr ) :: x0, x1, x2, x3, x4, xd
  real ( kind = pr ) :: sigma, sigma2, hb2m
  real ( kind = pr ) :: sigma_n, sigma_d, sigma2_n, sigma2_d
  real ( kind = pr ) :: w, b4, b4p, tt, tu
  !
  !! 3-body
  !
  real ( kind = pr ) :: u0, u1, y1, u2, y21, y22
  !
  !! 4-body
  !
  real ( kind = pr ) :: v0   
  !
  logical :: j2terms, cm1
  integer :: icoul
  !
  !! Coef. from the power of the density
  !
  real (pr) :: s11, s12, s21, s22
  !
  !! C coefficients
  !
  real ( kind = pr ) :: c_rho(-1:1,1:3),c_dr(-1:1),c_tau(-1:1)
  real ( kind = pr ) :: c_s(-1:1,1:3),c_nj(-1:1),c_ds(-1:1)
  real ( kind = pr ) :: c_t(-1:1),c_f(-1:1),c_ns(-1:1)
  real ( kind = pr ) :: a_t(0:1),a_ds(0:1)

  !
  !! 3 body
  !
  real ( kind = pr ) ::  Brho_0,Brho_1,Bs_0,Bs_10,Bs_1
  real ( kind = pr ) ::  Btau_0,Btau_10,Btau_1,BT_0,BT_10,BT_01,BT_1 
  real ( kind = pr ) ::  Btaus_0,Btaus_10,Btaus_1,Bdrho_0,Bdrho_10,Bdrho_1
  real ( kind = pr ) ::  Bds_0,Bds_10,Bds_1,Bdrs_0,Bdrs_01,Bdrs_10,Bdrs_1 
  real ( kind = pr ) ::  Bjp_0,Bjp_10,Bjp_1,BJ_0,BJ_10,BJ_1 
  real ( kind = pr ) ::  BJs_0,BJs_01,BJs_10,BJs_1,BdsJ_0,BdsJ_01,BdsJ_10,BdsJ_1

  !
  !! D wave
  !
  real ( kind = pr ) :: CD_tau(-1:1),CD_T(-1:1)
  logical :: flagToDoDwave=.false.
  !
  !!
  !
  real ( kind = pr ) :: Cr0_in,Cr0g_in,Cr0g2_in,Cr1_in,Cr1g_in,Cr1g2_in
  real ( kind = pr ) :: Cdr0_in,Cdr1_in,Ct0_in,Ct1_in,Cs0_in,Cs0g_in,Cs0g2_in
  real ( kind = pr ) :: Cs1_in,Cs1g_in,Cs1g2_in,CnJ0_in,CnJ1_in,Cds0_in,Cds1_in
  real ( kind = pr ) :: CgT0_in,CgT1_in,CF0_in,CF1_in,Cns0_in,Cns1_in
  real ( kind = pr ) :: CD_tau0_in,CD_tau1_in,CD_T0_in,CD_T1_in

  !
  !!
  !

  interface read_force
     module procedure read_force
  end interface

  interface read_functional
     module procedure read_functional
  end interface

contains

  !
  !!
  !

  subroutine read_force(nom)
    !==================================================================
    ! abstract: routine to read the parameters t_i,x_i
    !           of a Skyrme force
    !
    ! input: name of the force
    ! output: trough module lr_param
    !==================================================================

    implicit none
    character ( len = * ), intent(in) :: nom
    logical :: ok,Dwave
    character ( len = 20 ) :: name
    integer :: id,iwrite
    !
    namelist /skf/ name, t0, t1, t2, t3, t4, td, x0, x1, x2, x3, x4, xd,  &
         w, sigma, sigma_n, sigma_d, sigma2, sigma2_n, sigma2_d, j2terms, &
         hb2m, b4, b4p, tt, tu, cm1, icoul, u0, u1, u2, y1, y21, y22, v0
    !
    iwrite=0 ! to write or not on screen


    inquire( exist = ok, file = "forces.param" )
    if ( .not. ok ) then
       print '(" *** Can''t find the file ""forces.param""...")'
       stop
    end if
    !
    call open_file( id, "forces.param", "old" )
    do
       call lr_reset_force() ! reset force parameters
       read( id, nml = skf )

      !  write(*,*)name
       if ( trluc(name) == trluc(nom) ) exit
       if ( trluc(name) == "FIN" ) then
          print '(" *** ",a,": Unknown force in the ph channel...")', nom
          stop
       end if
    end do
    call close_file(id)
    !
    if ( sigma == 0.0_pr  .and. sigma_d  /= 0.0_pr )  &
         &                        sigma  = sigma_n / sigma_d
    if ( sigma2 == 0.0_pr .and. sigma2_d /= 0.0_pr )  &
         &                        sigma2 = sigma2_n / sigma2_d
    !
    if(iwrite.eq.1)then
       print '("*",19x,"--- Parameters of the force ---",20x,"*")'
       print '("*  ph channel:",57x,"*")'
       print '("*",5x,"t0 = ",f12.6,5x,"x0 = ",f10.6,28x,"*")', t0, x0
       print '("*",5x,"t1 = ",f12.6,5x,"x1 = ",f10.6,28x,"*")', t1, x1
       print '("*",5x,"t2 = ",f12.6,5x,"x2 = ",f10.6,28x,"*")', t2, x2
       print '("*",5x,"t3 = ",f12.6,5x,"x3 = ",f10.6,&
            &5x,"sigma  = ",f9.6,5x,"*")', t3, x3, sigma
       if ( t4 * x4 /= 0.0_pr ) &
            print '("*",5x,"t4 = ",f12.6,5x,"x4 = ",f10.6,&
            &5x,"sigma'' = ",f9.6,5x,"*")', t4, x4, sigma2
       print '("*",5x,"tD = ",f12.6,5x,"xD = ",f10.6,5x,"icoul = ",i1,&
            &14x,"*")', td, xd, icoul
       if ( b4 == 0.0_pr ) then
          print '("*",5x,"W0 = ",f12.6,5x,"J2_terms = ",l1,8x,&
               &"1-body cm = ",l1,10x,"*")',  &
               w, j2terms, cm1
       else
          print '("*",5x,"b4 = ",f12.6,4x,"b4p = ",f12.6,3x,"J2_terms = ",&
               &l1,11x,"*")', b4, b4p, j2terms
       end if
       print '("*",5x,"tt = ",f12.6,5x,"tu = ",f12.6,26x,"*")', tt, tu
       !
       print '("*  ph channel, 3-body part:",44x,"*")' ! 3-body
       print '("*",5x,"u0 = ",f12.6,48x,"*")', u0
       print '("*",5x,"u1 = ",f12.6,5x,"y1 = ",f12.6,26x,"*")', u1, y1
       print '("*",5x,"u2 = ",f12.6,5x,"y2 = ",f12.6,26x,"*")', u2, y21,y22
       !
       write( *, '("*")' )
       print '("*  ph channel, 4-body part:",44x,"*")' ! 4-body
       print '("*",5x,"u0 = ",f15.6,48x,"*")', v0
    endif
    



    Dwave=.false.
    !--- flag to check if we have or not the 3 body term
    if(td.ne.0.0_pr.or.xd.ne.0.0_pr)Dwave=.true.

    if(Dwave)stop 'The code does not allow D wave and tensor'

    !
  end subroutine read_force

  !
  !!
  !

  subroutine lr_reset_force()
    implicit none
    !
    !! Default values
    !
    t0       = 0.0_pr
    t1       = 0.0_pr
    t2       = 0.0_pr
    t3       = 0.0_pr
    t4       = 0.0_pr
    x0       = 0.0_pr
    x1       = 0.0_pr
    x2       = 0.0_pr
    x3       = 0.0_pr
    x4       = 0.0_pr
    w        = 0.0_pr
    sigma    = 0.0_pr
    sigma_n  = 0.0_pr
    sigma_d  = 1.0_pr
    sigma2   = 0.0_pr
    sigma2_n = 0.0_pr
    sigma2_d = 1.0_pr
    b4       = 0.0_pr
    b4p      = 0.0_pr
    tt       = 0.0_pr
    tu       = 0.0_pr
    !
    !! D-wave
    !
    td       = 0.0_pr
    xd       = 0.0_pr
    !
    !! 3-body
    !
    u0       = 0.0_pr
    u1       = 0.0_pr
    u2       = 0.0_pr
    y1       = 0.0_pr
    y21      = 0.0_pr
    y22      = 0.0_pr
    !
    !! 4-body
    !
    v0       = 0.0_pr
    !
    j2terms  = .true.
    cm1      = .true.
    hb2m = 20.73553_pr ! ATTENZIONE (modifica) hbar_c**2 / 2 / nucleon_mass0
    !
    icoul = 1
    !
  end subroutine lr_reset_force

  !
  !!
  !

  subroutine read_functional(nom)
    !==================================================================
    ! abstract: routine to read the parameters C^X_i
    !           of a Skyrme functional (2 and 3 body!!)
    !
    ! input: name of the force
    ! output: trough module lr_param
    !==================================================================

    implicit none
    character ( len = * ), intent(in) :: nom
    logical :: ok,troiscorp,Dwave
    character ( len = 20 ) :: name
    integer :: id,iwrite
    real (kind=pr),parameter :: z=0.0_pr
    !
    namelist /fct/ name,Cr0_in,Cr0g_in,Cr0g2_in,Cr1_in,Cr1g_in,Cr1g2_in,    &
         Cdr0_in,Cdr1_in,Ct0_in,Ct1_in,Cs0_in,Cs0g_in,Cs0g2_in,   &
         Cs1_in,Cs1g_in,Cs1g2_in,CnJ0_in,CnJ1_in,Cds0_in,Cds1_in, &
         CgT0_in,CgT1_in,CF0_in,CF1_in,Cns0_in,Cns1_in,sigma,     &   
         sigma2,J2terms,hb2m,Brho_0,Brho_1,Bs_0,Bs_10,Bs_1,Btau_0,&
         Btau_10,Btau_1,BT_0,BT_10,BT_01,BT_1,Btaus_0,Btaus_10,   &
         Btaus_1,Bdrho_0,Bdrho_10,Bdrho_1,Bds_0,Bds_10,Bds_1,     &
         Bdrs_0,Bdrs_01,Bdrs_10,Bdrs_1,Bjp_0,Bjp_10,Bjp_1,BJ_0,   &
         BJ_10,BJ_1,BJs_0,BJs_01,BJs_10,BJs_1,BdsJ_0,BdsJ_01,     &
         BdsJ_10,BdsJ_1,v0,CD_tau0_in,CD_tau1_in,CD_T0_in,CD_T1_in
    !
    iwrite=0 ! to write or not on screen

    inquire( exist = ok, file = "functional.param" )
    if ( .not. ok ) then
       print '(" *** Can''t find the file ""functional.param""...")'
       stop
    end if
    !
    call open_file( id, "functional.param", "old" )
    do
       call reset_functional() ! reset functional parameters
       read( id, nml = fct )

       if ( trluc(name) == trluc(nom) ) exit
       if ( trluc(name) == "FIN" ) then
          print '(" *** ",a,": Unknown functional in the ph channel...")', nom
          stop
       end if
    end do
    call close_file(id)

    troiscorp=.false.
    !--- flag to check if we have or not the 3 body term
    if(Brho_0.ne.z.or.Brho_1.ne.z.or.Bs_0.ne.z.or.Bs_10.ne.z.or.Bs_1.ne.z.or. &
         Btau_0.ne.z.or.Btau_10.ne.z.or.Btau_1.ne.z.or.BT_0.ne.z.or.BT_10.ne.z   &
         .or.BT_01.ne.z.or.BT_1.ne.z.or.Btaus_0.ne.z.or.Btaus_10.ne.z.or.Btaus_1 &
         .ne.z.or.Bdrho_0.ne.z.or.Bdrho_10.ne.z.or.Bdrho_1.ne.z.or.Bds_0.ne.z.or.&
         Bds_10.ne.z.or.Bds_1.ne.z.or.Bdrs_0.ne.z.or.Bdrs_01.ne.z.or.Bdrs_10.ne.z&
         .or.Bdrs_1.ne.z.or.Bjp_0.ne.z.or.Bjp_10.ne.z.or.Bjp_1.ne.z.or.BJ_0.ne.z & 
         .or.BJ_10.ne.z.or.BJ_1.ne.z.or.BJs_0.ne.z.or.BJs_01.ne.z.or.BJs_10.ne.z &
         .or.BJs_1.ne.z.or.BdsJ_0.ne.z.or.BdsJ_01.ne.z.or.BdsJ_10.ne.z.or.       &
         BdsJ_1.ne.z)troiscorp=.true.

    !--- flag to check if we have or not the D wave
    Dwave=.false.
    if(CD_tau0_in.ne.0.0_pr.or.CD_tau1_in.ne.0.0_pr.or.CD_T0_in.ne.0.0_pr &
         .or.CD_T1_in.ne.0.0_pr)Dwave=.true.
    !
    if(iwrite.eq.1)then
       write(*,*) '--- Parameters of the functional ---'
       write(*,*)'C^{\rho}_0=',Cr0_in,'\rho^\gamma = ',Cr0g_in,'\rho^\gamma2 =',Cr0g2_in
       write(*,*)'C^{\rho}_1=',Cr1_in,'\rho^\gamma =',Cr1g_in,'\rho^\gamma2 =',Cr1g2_in
       write(*,*)'C^{\Delta\rho}_0=',Cdr0_in,'C^{\Delta\rho}_1=',Cdr1_in
       write(*,*)'C^{\tau}_0= ',Ct0_in,'C^{\tau}_1= ',Ct1_in
       write(*,*)'C^s_0= ',Cs0_in,'\rho^\gamma = ',Cs0g_in,'\rho^\gamma2 =',Cs0g2_in
       write(*,*)'C^s_1= ',Cs1_in,'\rho^\gamma = ',Cs1g_in,'\rho^\gamma2 =',Cs1g2_in
       write(*,*)'C^\nablaJ_0= ',CnJ0_in,'C^\nablaJ_1= ',CnJ1_in
       write(*,*)'C^{\Delta s}_0= ',Cds0_in,'C^\{Delta s}_1= ',Cds1_in
       write(*,*)'C^T_0= ',CgT0_in,'C^T_1= ',CgT1_in
       write(*,*)'C^F_0= ',CF0_in,'C^F_1= ',CF1_in
       write(*,*)'C{\nabla s}_0= ',Cns0_in,'C{\nabla s}_1= ',Cns1_in
       write(*,*)'gamma= ',sigma,'gamma2= ',sigma2
       if(troiscorp)then
          write(*,*)'Three body terms',Brho_0,Brho_1,Bs_0,Bs_10,Bs_1,Btau_0,&
               Btau_10,Btau_1,BT_0,BT_10,BT_01,BT_1,Btaus_0,Btaus_10,   &
               Btaus_1,Bdrho_0,Bdrho_10,Bdrho_1,Bds_0,Bds_10,Bds_1,     &
               Bdrs_0,Bdrs_01,Bdrs_10,Bdrs_1,Bjp_0,Bjp_10,Bjp_1,BJ_0,   &
               BJ_10,BJ_1,BJs_0,BJs_01,BJs_10,BJs_1,BdsJ_0,BdsJ_01,     &
               BdsJ_10,BdsJ_1
       endif
       if(v0.ne.0)then
          write(*,*)'4 body = ',v0
       endif
       if(Dwave)then
          write(*,*)'D wave term'
          write(*,*)'CD_tau_0= ',CD_tau0_in,' CD_tau_1= ',CD_tau1_in
          write(*,*)'CD_T_0= ',CD_T0_in,' CD_T_1= ',CD_T1_in
       endif
    endif

    if(Dwave)then
       write(*,*)'The code do not support D wave '
       stop
    endif

    !
  contains

    !
    !!
    !

    subroutine reset_functional()
      implicit none
      real (kind=pr),parameter :: z=0.0_pr
      !
      !! Default values
      !
      Cr0_in=z; Cr0g_in=z; Cr0g2_in=z; Cr1_in=z; Cr1g_in=z; Cr1g2_in=z
      Cdr0_in=z; Cdr1_in=z; Ct0_in=z; Ct1_in=z; Cs0_in=z
      Cs0g_in=z; Cs0g2_in=z; Cs1_in=z; Cs1g_in=z; Cs1g2_in=z;
      CnJ0_in=z; CnJ1_in=z; Cds0_in=z; Cds1_in=z 
      CgT0_in=z; CgT1_in=z; CF0_in=z; CF1_in=z; Cns0_in=z; Cns1_in=z
      sigma=z
      sigma2=z
      !
      !! 3-body
      !
      Brho_0=z; Brho_1=z; Bs_0=z; Bs_10=z; Bs_1=z; Btau_0=z 
      Btau_10=z; Btau_1=z; BT_0=z; BT_10=z; BT_01=z; BT_1=z 
      Btaus_0=z; Btaus_10=z; Btaus_1=z; Bdrho_0=z; Bdrho_10=z
      Bdrho_1=z; Bds_0=z; Bds_10=z; Bds_1=z; Bdrs_0=z; Bdrs_01=z
      Bdrs_10=z; Bdrs_1=z; Bjp_0=z; Bjp_10=z; Bjp_1=z; BJ_0=z
      BJ_10=z; BJ_1=z; BJs_0=z; BJs_01=z; BJs_10=z; BJs_1=z
      BdsJ_0=z; BdsJ_01=z; BdsJ_10=z; BdsJ_1=z
      !
      !! 4 body
      !
      v0=z
      !
      !! D wave
      !
      CD_tau0_in=z; CD_tau1_in=z
      CD_T0_in=z; CD_T1_in=z
      !
      j2terms  = .true.
      cm1      = .true.
      hb2m = 20.73553! ATTENZIONE (modifica) hbar_c**2 / 2 / nucleon_mass0
      !
      icoul = 1
      !
    end subroutine reset_functional

    !
    !!
    !

  end subroutine read_functional

  !
  !!
  !

  subroutine set_coeff_functional(zero,readforce)
    !======================================================
    ! abstract: routine that determines the coefficient C
    !           of the functional (coming from a force)
    !
    !input:
    !zero: character that determines if and which parameters 
    !      are set to zero
    !
    !output
    !trough module
    !=====================================================
    implicit none
    logical, intent(in) :: readforce
    character ( len = 4 ) :: zero
    integer :: iwrite
    ! useful variables


    iwrite=0
    c_rho=0.0_pr ! it's always safe to put to zero the coefficients
    c_dr =0.0_pr
    c_tau=0.0_pr
    c_s  =0.0_pr
    c_nj =0.0_pr
    c_ds =0.0_pr
    c_t  =0.0_pr
    c_f  =0.0_pr
    c_ns =0.0_pr

    CD_tau=0.0_pr
    CD_T=0.0_pr

    !
    !-- Coef. from the power of the density
    !
    s11 = sigma  + 1
    s12 = sigma  + 2
    s21 = sigma2 + 1
    s22 = sigma2 + 2
    !
    !------------------ ATTENTION IMPORTANT!!!!!  -----------------------------
    !  The first index of the vector is the isospin, it can vary from -1 to 1,
    !  when it gets values 0,1 is the real isospin, when it assumes -1 is the
    !  fake isospin used to have a similar structure for pure neutron matter!!
    !  The second index refers to the  density depenendent part that
    !  consistutes the coefficient C^x(\rho)
    !--------------------------------------------------------------------------


    if(readforce)then ! from a force
       call from_force(zero)
    else
       call from_functional(zero)
    endif
    !======= Coefficients for PNM case (2 body)
    ! 
    c_rho(-1,1)=c_rho(0,1)+c_rho(1,1) ! ! C^{\rho}_i
    c_rho(-1,2)=c_rho(0,2)+c_rho(1,2)
    c_rho(-1,3)=c_rho(0,3)+c_rho(1,3)
    !
    c_dr(-1)=c_dr(0)+c_dr(1) ! ! C^{\Delta \rho}_iP
    !
    c_tau(-1)=c_tau(0)+c_tau(1) ! C^{\tau}_i
    !
    c_s(-1,1)=c_s(0,1)+c_s(1,1) !  C^s_i
    c_s(-1,2)=c_s(0,2)+c_s(1,2)
    c_s(-1,3)=c_s(0,3)+c_s(1,3)
    !
    c_nj(-1)=0.5_pr*(c_nj(0)+c_nj(1)) ! C^{\nabla J}_i
    !
    c_ds(-1)=(c_ds(0)+c_ds(1)) ! C^{\delta s}_i
    !
    c_t(-1)=c_t(0)+c_t(1) ! C^T_i
    !
    c_f(-1)=0.5_pr*(c_f(0)+c_f(1)) ! C^F_i
    !
    c_ns(-1)=0.5_pr*(c_ns(0)+c_ns(1)) ! C^{\nabla s}_i

    ! --- Dwave
    CD_tau(-1)=CD_tau(0)+CD_tau(1)
    CD_T(-1)=CD_T(0)+CD_T(1)


    !===========================================================================
    ! writing the coefficients on screen
    !

    if(iwrite.eq.1)then
       write(*,*)
       write(*,*)'************************************************************'
       write(*,*)'Resume of the coefficient fo the functional used in the code'
       write(*,*)'************************************************************'
       write(*,*)
       write(*,*)
       write(*,*)'C^{\rho}_0= ',c_rho(0,1),'+ \rho^\gamma ',c_rho(0,2),'+ \rho^\gamma2 ',c_rho(0,3)
       write(*,*)'C^{\rho}_1= ',c_rho(1,1),'+ \rho^\gamma ',c_rho(1,2),'+ \rho^\gamma2 ',c_rho(1,3)
       write(*,*)'C^{\Delta \rho}_0= ',c_dr(0),' C^{\Delta \rho}_1= ',c_dr(1)
       write(*,*)'C^{\tau}_0= ',c_tau(0),' C^{\tau}_1= ',c_tau(1)
       write(*,*)'C^{s}_0= ',c_s(0,1),'+ \rho^\gamma ',c_s(0,2),'+ \rho^\gamma2 ',c_s(0,3)
       write(*,*)'C^{s}_1= ',c_s(1,1),'+ \rho^\gamma ',c_s(1,2),'+ \rho^\gamma2 ',c_s(1,3)
       write(*,*)'C^{\nabla J}_0= ',c_nj(0),' C^{\nabla J}_1= ',c_nj(1)
       write(*,*)'C^{\delta s}_0= ',c_ds(0),' C^{\delta s}_1= ',c_ds(1)
       write(*,*)'C^{T}_0= ',c_t(0),' C^{T}_1= ',c_t(1)
       write(*,*)'C^{F}_0= ',c_f(0),' C^{F}_1= ',c_f(1)
       write(*,*)'C^{\nabla s}_0= ',c_ns(0),' C^{\nabla s}_1= ',c_ns(1)
       write(*,*)
       write(*,*)' gamma = ',sigma,' gamma2 = ',sigma2
       write(*,*)

    endif
  contains

    subroutine from_force(zero)
      !
      !! building the C^X of functional starting form Skyrme x t parameters
      !
      character ( len = 4 ) :: zero
      real ( kind = pr ),parameter :: a12=1.0_pr/2.0_pr
      real ( kind = pr ),parameter :: a54=5.0_pr/4.0_pr
      real ( kind = pr ),parameter :: a8=8.0_pr
      real ( kind = pr ),parameter :: a16=16.0_pr
      real ( kind = pr ),parameter :: a32=32.0_pr
      real ( kind = pr ),parameter :: a64=64.0_pr
      real ( kind = pr ),parameter :: a128=128.0_pr
      real ( kind = pr ),parameter :: a256=256.0_pr
      real ( kind = pr ) :: b_t(0:1),b_f(0:1),b_ds(0:1),b_ns(0:1)
      real ( kind = pr ) :: t_e,t_o

      a_ds=0.0_pr
      a_t=0.0_pr
      b_t=0.0_pr
      b_f=0.0_pr
      b_ds=0.0_pr
      b_ns=0.0_pr
    
      ! 
      !! C^{\rho}_i
      !
      c_rho(0,1)=3.0_pr/8.0_pr*t0
      c_rho(0,2)=3.0_pr/48.0_pr*t3
      c_rho(0,3)=3.0_pr/48.0_pr*t4
      c_rho(1,1)=-t0/4.0_pr*(a12+x0)
      c_rho(1,2)=-t3/24.0_pr*(a12+x3)
      c_rho(1,3)=-t4/24.0_pr*(a12+x4)
      !
      !! C^{\Delta \rho}_i
      !
      c_dr(0)=-9.0_pr/64.0_pr*t1+t2/16.0_pr*(a54+x2)
      c_dr(1)=3.0_pr/32.0_pr*t1*(a12+x1)+t2/32.0_pr*(a12+x2)
      !
      !! C^{\tau}_i
      !
      c_tau(0)=3.0_pr/16.0_pr*t1+t2/4.0_pr*(a54+x2)
      c_tau(1)=-t1/8.0_pr*(a12+x1)+t2/8.0_pr*(a12+x2)
      !
      !! C^s_i
      !
      c_s(0,1)=-t0/4.0_pr*(a12-x0)
      c_s(0,2)=-t3/24.0_pr*(a12-x3)
      c_s(0,3)=-t4/24.0_pr*(a12-x4)
      c_s(1,1)=-t0/8.0_pr
      c_s(1,2)=-t3/48.0_pr
      c_s(1,3)=-t4/48.0_pr
      !
      !! C^{\nabla J}_i
      !
  
      if(b4.eq.0.0_pr.and.b4p.eq.0.0_pr)then ! classica spin orbit
         c_nj(0)=-3.0_pr/4.0_pr*w
         c_nj(1)=-1.0_pr/4.0_pr*w
      else ! Pauli-violating spin orbit
         ! we use the notation by Reinhard -> b_4(p)^Reinhard=0.5 * b_4(p)^Chabanat!
         c_nj(0)=(-b4-b4p/2.0_pr)
         c_nj(1)=-b4p/2.0_pr
      endif
      
      !
      !! A^{\Delta s}_i
      !
      a_ds(0)=3.0_pr/32.0_pr*t1*(a12-x1)+t2/32.0_pr*(a12+x2)
      a_ds(1)=3.0_pr/64.0_pr*t1+t2/64.0_pr

      if( j2terms ) then
         !
         !! A^T_i
         !
         a_t(0)=-t1/8.0_pr*(a12-x1)+t2/8.0_pr*(a12+x2)
         a_t(1)=-t1/16.0_pr+t2/16.0_pr
         !---- here the tensor (coming from a force)-----
         t_e = tt / 3.0_pr
         t_o = tu / 3.0_pr
         
         !
         !! B^T_i
         !
         b_t(0)= - (  t_e+3.0_pr*t_o)/8.0_pr
         b_t(1)=   (  t_e-       t_o)/8.0_pr
         !
         !! B^F_i
         !
         b_f(0)= 3.0_pr*(t_e+3.0_pr*t_o)/8.0_pr
         b_f(1)=-3.0_pr*(t_e-t_o)       /8.0_pr
         !
         !! B^{\Delta s}_i
         !
         b_ds(0)=3.0_pr*(t_e       -t_o)/32.0_pr
         b_ds(1)= -     (3.0_pr*t_e+t_o)/32.0_pr
         !
         !! B^{\nabla s}_i
         !
         b_ns(0)=9.0_pr*(t_e-t_o)/32.0_pr
         b_ns(1)=-3.0_pr*(3.0_pr*t_e+t_o)/32.0_pr
      else
         a_t(0)= 0.0_pr
         a_t(1)= 0.0_pr
         t_e = 0.0_pr
         t_o = 0.0_pr
         b_t(0)= 0.0_pr
         b_t(1)= 0.0_pr 
         b_f(0)= 0.0_pr
         b_f(1)= 0.0_pr
         b_ds(0)= 0.0_pr
         b_ds(1)= 0.0_pr
         b_ns(0)= 0.0_pr
         b_ns(1)= 0.0_pr
      end if

      !=========== 3 body ============================

      Brho_0 = + 3.0_pr * u0/a16
      Brho_1 = - 3.0_pr * u0/a16
      !
      Bs_0 = - 3.0_pr * u0/a16
      Bs_10= + 3.0_pr * u0/a8
      Bs_1 = - 3.0_pr * u0/a16
      !
      Btau_0 = + 3.0_pr*u1/a32 + 15.0_pr*u2/a64 + 3.0_pr*u2*y21/a16 + 3.0_pr*u2*y22/a32
      Btau_10 = - u1/a32 + u1*y1 /a32 - 5.0_pr*u2/a64 -u2*y21/a16 -7.0_pr*u2*y22/a32
      Btau_1 = - u1/a16  - u1*y1/a32 + u2/a32 + u2*y21/a16 - u2*y22/a16
      !
      BT_0 = -u1/a16 + u1*y1/a32 + u2/a32 + u2*y21/a16 +u2*y22/a8
      BT_10=  u1/a16 - u1*y1/a32 -u2/a32 -u2*y21/a16 -u2*y22/a8
      BT_01=  u1/a16 - u2/a32
      BT_1 = -u1/a16 + u2/a32
      !
      Btaus_0 = -u1/a32 - u1*y1/a32 -5.0_pr*u2/a64 -u2*y21/a16 +5.0_pr*u2*y22/a32
      Btaus_10= -u1/a32 -5.0_pr*u2/a64 -u2*y21/a16 -u2*y22/a32
      Btaus_1 = u1/a16 + u1*y1/a32 - u2/a32 - u2*y21/a16 + u2*y22/a16
      !
      Bdrho_0  = + 15.0_pr*u1/a128 - 15.0_pr*u2/a256 - 3.0_pr*u2*y21/a64 - 3.0_pr* u2 * y22/a128
      Bdrho_10 = - 5.0_pr*u1/a64 + u1*y1/a32 + 5.0_pr*u2/a128 + u2*y21/a32 + 7.0_pr * u2*y22/a64
      Bdrho_1  = - 5.0_pr*u1/a128 - u1*y1/a32 - 7.0_pr*u2/a256- u2*y21/a32 - 5.0_pr *u2*y22/a128
      !
      Bds_0 = - 5.0_pr*u1/a128 + u1*y1/a32 -7.0_pr*u2/a256 -u2*y21/a32 + u2*y22/a128
      Bds_10=   5.0_pr*u1/a64 - u1*y1/a32 +u2/a128 + 3.0_pr*u2*y22/a64
      Bds_1=  - 5.0_pr*u1/a128 - 7.0_pr*u2/a256 - u2*y21/a64 - u2*y22/a128
      !
      Bdrs_0=  -5.0_pr*u1/a64 - u1*y1/a32 + 5.0_pr*u2/a128 + u2*y21/a32 - 5.0_pr*u2*y22/a64
      Bdrs_01= -5.0_pr*u1/a64 + 5.0_pr*u2/a128 +u2*y21/a32 +u2*y22/a64
      Bdrs_10=  5.0_pr*u1/a64 + u2/a128 + u2*y21/a32 + u2*y22/a64
      Bdrs_1 =  5.0_pr*u1/a64 + u1*y1/a32 + u2/a128 -3.0_pr*u2*y22/a64
      !
      Bjp_0 = -3.0_pr*u1/a32  -15.0_pr*u2/a64 -3.0_pr*u2*y21/a16 -3.0_pr*u2*y22/a32
      Bjp_10= u1/a16 - u1*y1/a16 + 5.0_pr*u2/a32 + u2*y21/a8 + 7.0_pr*u2*y22/a16
      Bjp_1 = +u1/a32 +u1*y1/a16 -7.0_pr*u2/a64 -u2*y21/a8 -5.0_pr*u2*y22/a32
      !
      BJ_0  = + u1/a32 - u1*y1/a16 - 7.0_pr*u2/a64 - u2*y21/a8 + u2*y22/a32
      BJ_10 = - u1/a16 + u1*y1/a16 + u2/a32 + 3.0_pr *u2*y22/a16
      BJ_1  = + u1/a32 - 7.0_pr *u2/a64 - u2*y21/a16 - u2*y22/a32
      !
      BJs_0 =  u1/a16 + u1*y1/a16 + 5.0_pr*u2/a32 + u2*y21/a8 - 5.0_pr*u2*y22/a16
      BJs_01=  u1/a16 + 5.0_pr*u2/a32 +u2*y21/a8 + u2*y22/a16
      BJs_10= -u1/a16 + u2/a32 +u2*y21/a8 + u2*y22/a16
      BJs_1=  -u1/a16 - u1*y1/a16 + u2/a32 -3.0_pr*u2*y22/a16
      !
      BdsJ_0 = -3.0_pr*u2/a64 - 3.0_pr*u2*y21/a32 +3.0_pr *u2*y22/a32
      BdsJ_01= u1*y1/a16 -3.0_pr*u2/a64 -u2*y21/a32 + u2*y22/a32
      BdsJ_10= -u1*y1/a32 - 3.0_pr*u2/a64 - u2*y21/a32 + u2*y22/a32
      BdsJ_1 = -u1*y1/a32 -3.0_pr*u2/a64 -u2*y21/a32 + u2*y22/a32 


      !======= D wave

      CD_tau(0)= 3.0_pr/16.0_pr * tD
      CD_tau(1)=-1.0_pr/8.0_pr  * tD * (1.0_pr/2.0_pr+xD)
      CD_T(0)=  -1.0_pr/4.0_pr  * tD * (1.0_pr/2.0_pr-xD)
      CD_T(1)=  -1.0_pr/8.0_pr  * tD

      !
      !==============================================================

      !===== manipulating the coefficients C ========================
      select case(trluc(zero))

      case("CTAU") ! putting C^{\nabla J}_i=0
         c_tau(0)=0.0_pr
         c_tau(1)=0.0_pr
      case("SORB") ! putting C^{\nabla J}_i=0
         c_nj(0)=0.0_pr
         c_nj(1)=0.0_pr
      case("TENS") ! putting all the B^{T,F,\Delta s,\nabla s}_i to zero
         b_t(0)=  0.0_pr
         b_t(1)=  0.0_pr  
         b_f(0)=  0.0_pr
         b_f(1)=  0.0_pr
         b_ds(0)= 0.0_pr
         b_ds(1)= 0.0_pr
         b_ns(0)= 0.0_pr
         b_ns(1)= 0.0_pr

      case("SOTE") ! putting all the C^{\nabla J}_i B^{T,F,\Delta s,\nabla s}_i to zero
         c_nj(0)=0.0_pr
         c_nj(1)=0.0_pr
         c_nj(-1)=0.0_pr
!!$       write(*,*)'Putting the spin-orbit term to zero'
!!$       write(*,*)' C^{\nabla J}_i '
         b_t(0)=  0.0_pr
         b_t(1)=  0.0_pr  
         b_f(0)=  0.0_pr
         b_f(1)=  0.0_pr
         b_ds(0)= 0.0_pr
         b_ds(1)= 0.0_pr
         b_ns(0)= 0.0_pr
         b_ns(1)= 0.0_pr
!!$       write(*,*)' Putting ALL the coefficients B of the tensor to Zero'
!!$       write(*,*)' B^T_i',' B^F_i',' B^{\Delta s}_i',' B^{\nabla s}_i'
      case("3BOD") ! putting all the 3 body terms to zero
         Brho_0 = 0.0_pr ;  Brho_1 = 0.0_pr
         Bs_0 =  0.0_pr  ;  Bs_10=  0.0_pr ; Bs_1 =  0.0_pr 
         Btau_0 = 0.0_pr ;  Btau_10 =  0.0_pr ; Btau_1 =  0.0_pr 
         BT_0 =  0.0_pr  ;  BT_10=   0.0_pr ; BT_01=   0.0_pr ; BT_1 =  0.0_pr 
         Btaus_0 =  0.0_pr ; Btaus_10=  0.0_pr ;  Btaus_1 =  0.0_pr 
         Bdrho_0  =  0.0_pr ; Bdrho_10 =  0.0_pr ; Bdrho_1  =  0.0_pr 
         Bds_0 =  0.0_pr ; Bds_10=  0.0_pr ; Bds_1=  0.0_pr 
         Bdrs_0=   0.0_pr ; Bdrs_01=  0.0_pr ; Bdrs_10=  0.0_pr ; Bdrs_1 =  0.0_pr 
         Bjp_0 =  0.0_pr ;  Bjp_10=  0.0_pr ; Bjp_1 =  0.0_pr 
         BJ_0  =  0.0_pr ;  BJ_10 =  0.0_pr ; BJ_1  =  0.0_pr 
         BJs_0 =  0.0_pr ;  BJs_01= 0.0_pr ;  BJs_10= 0.0_pr ;  BJs_1= 0.0_pr   
         BdsJ_0 =  0.0_pr ; BdsJ_01= 0.0_pr ;BdsJ_10=  0.0_pr ; BdsJ_1 = 0.0_pr  
      case("4BOD") ! putting all the 4 body coefficients to zero
         v0=0.0_pr
      case("34BD") ! putting all the 3 body and 4 body terms to zero
         Brho_0 = 0.0_pr ;  Brho_1 = 0.0_pr
         Bs_0 =  0.0_pr  ;  Bs_10=  0.0_pr ; Bs_1 =  0.0_pr
         Btau_0 = 0.0_pr ;  Btau_10 =  0.0_pr ; Btau_1 =  0.0_pr
         BT_0 =  0.0_pr  ;  BT_10=   0.0_pr ; BT_01=   0.0_pr ; BT_1 =  0.0_pr
         Btaus_0 =  0.0_pr ; Btaus_10=  0.0_pr ;  Btaus_1 =  0.0_pr
         Bdrho_0  =  0.0_pr ; Bdrho_10 =  0.0_pr ; Bdrho_1  =  0.0_pr
         Bds_0 =  0.0_pr ; Bds_10=  0.0_pr ; Bds_1=  0.0_pr
         Bdrs_0=   0.0_pr ; Bdrs_01=  0.0_pr ; Bdrs_10=  0.0_pr ; Bdrs_1 =  0.0_pr
         Bjp_0 =  0.0_pr ;  Bjp_10=  0.0_pr ; Bjp_1 =  0.0_pr
         BJ_0  =  0.0_pr ;  BJ_10 =  0.0_pr ; BJ_1  =  0.0_pr
         BJs_0 =  0.0_pr ;  BJs_01= 0.0_pr ;  BJs_10= 0.0_pr ;  BJs_1= 0.0_pr
         BdsJ_0 =  0.0_pr ; BdsJ_01= 0.0_pr ;BdsJ_10=  0.0_pr ; BdsJ_1 = 0.0_pr
         v0=0.0_pr
      case("DWAV")! putting all the D wave to zero
         CD_tau= 0.0_pr; CD_T= 0.0_pr
      case("NONE")
!!$       write(*,*)'ALL coefficients are taken into account'
      case default
         stop 'I do not understand whic set of coefficients C you want!'

      end select

      !======== summing the coefficents A and B =====
      !
      ! C^{\delta s}_i
      !
      c_ds(0)=a_ds(0)+b_ds(0)
      c_ds(1)=a_ds(1)+b_ds(1)
      !
      ! C^T_i
      !
      c_t(0)=a_t(0)+b_t(0)
      c_t(1)=a_t(1)+b_t(1)
      !
      ! C^F_i
      !
      c_f(0)=b_f(0)
      c_f(1)=b_f(1)
      !
      ! C^{\nabla s}_i
      !
      c_ns(0)=b_ns(0)
      c_ns(1)=b_ns(1)

    end subroutine from_force
    !
    !!
    !
    subroutine from_functional(zero)
      !
      !! building the C^X of functional (not a Force!!!!)
      !
      character ( len = 4 ) :: zero

      ! 
      !! C^{\rho}_i
      !
      c_rho(0,1)=Cr0_in
      c_rho(0,2)=Cr0g_in
      c_rho(0,3)=Cr0g2_in
      c_rho(1,1)=Cr1_in
      c_rho(1,2)=Cr1g_in
      c_rho(1,3)=Cr1g2_in
      !
      !! C^{\Delta \rho}_i
      !
      c_dr(0)=Cdr0_in
      c_dr(1)=Cdr1_in
      !
      !! C^{\tau}_i
      !
      c_tau(0)=Ct0_in
      c_tau(1)=Ct1_in
      !
      !! C^s_i
      !
      c_s(0,1)=Cs0_in
      c_s(0,2)=Cs0g_in
      c_s(0,3)=Cs0g2_in
      c_s(1,1)=Cs1_in
      c_s(1,2)=Cs1g_in
      c_s(1,3)=Cs1g2_in
      !
      !! C^{\nabla J}_i
      !
      c_nj(0)=CnJ0_in
      c_nj(1)=CnJ1_in

      if( j2terms ) then
         !
         !! C^T_i
         !
         c_t(0)= CgT0_in
         c_t(1)= CgT1_in
         !
         !! C^F_i
         !
         c_f(0)= CF0_in
         c_f(1)= CF1_in
         !
         !! C^{\Delta s}_i
         !
         c_ds(0)=Cds0_in
         c_ds(1)=Cds1_in
         !
         !! C^{\nabla s}_i
         !
         c_ns(0)=Cns0_in
         c_ns(1)=Cns1_in
      else
         c_t(0)= 0.0_pr
         c_t(1)= 0.0_pr
         c_f(0)= 0.0_pr
         c_f(1)= 0.0_pr
         c_ds(0)= 0.0_pr
         c_ds(1)= 0.0_pr
         c_ns(0)= 0.0_pr
         c_ns(1)= 0.0_pr
      end if

      !=========== 3 body ===========================================
      ! if present the coefficients are directly read from file
      !==============================================================

      !========== D wave ============================================
    
      CD_tau(0)=CD_tau0_in
      CD_tau(1)=CD_tau1_in
      CD_T(0)=CD_T0_in
      CD_T(1)=CD_T1_in

      !===== manipulating the coefficients C ========================
      select case(zero)

      case("SORB") ! putting C^{\nabla J}_i=0
         c_nj(0)=0.0_pr
         c_nj(1)=0.0_pr
      case("TENS") ! putting all the C^{T,F,\Delta s,\nabla s}_i to zero
         c_t(0)=  0.0_pr
         c_t(1)=  0.0_pr  
         c_f(0)=  0.0_pr
         c_f(1)=  0.0_pr
         c_ds(0)= 0.0_pr
         c_ds(1)= 0.0_pr
         c_ns(0)= 0.0_pr
         c_ns(1)= 0.0_pr
      case("SOTE") ! putting all the C^{\nabla J}_i C^{T,F,\Delta s,\nabla s}_i to zero
         c_nj(0)=0.0_pr
         c_nj(1)=0.0_pr
         c_t(0)=  0.0_pr
         c_t(1)=  0.0_pr  
         c_f(0)=  0.0_pr
         c_f(1)=  0.0_pr
         c_ds(0)= 0.0_pr
         c_ds(1)= 0.0_pr
         c_ns(0)= 0.0_pr
         c_ns(1)= 0.0_pr
      case("3BOD") ! putting all the 3 body terms to zero
         Brho_0 = 0.0_pr ;  Brho_1 = 0.0_pr
         Bs_0 =  0.0_pr  ;  Bs_10=  0.0_pr ; Bs_1 =  0.0_pr 
         Btau_0 = 0.0_pr ;  Btau_10 =  0.0_pr ; Btau_1 =  0.0_pr 
         BT_0 =  0.0_pr  ;  BT_10=   0.0_pr ; BT_01=   0.0_pr ; BT_1 =  0.0_pr 
         Btaus_0 =  0.0_pr ; Btaus_10=  0.0_pr ;  Btaus_1 =  0.0_pr 
         Bdrho_0  =  0.0_pr ; Bdrho_10 =  0.0_pr ; Bdrho_1  =  0.0_pr 
         Bds_0 =  0.0_pr ; Bds_10=  0.0_pr ; Bds_1=  0.0_pr 
         Bdrs_0=   0.0_pr ; Bdrs_01=  0.0_pr ; Bdrs_10=  0.0_pr ; Bdrs_1 =  0.0_pr 
         Bjp_0 =  0.0_pr ;  Bjp_10=  0.0_pr ; Bjp_1 =  0.0_pr 
         BJ_0  =  0.0_pr ;  BJ_10 =  0.0_pr ; BJ_1  =  0.0_pr 
         BJs_0 =  0.0_pr ;  BJs_01= 0.0_pr ;  BJs_10= 0.0_pr ;  BJs_1= 0.0_pr   
         BdsJ_0 =  0.0_pr ; BdsJ_01= 0.0_pr ;BdsJ_10=  0.0_pr ; BdsJ_1 = 0.0_pr
      case("4BOD") ! putting all the 4 body coefficients to zero
         v0=0.0_pr  
      case("DWAV")! putting all the D wave to zero
         CD_tau= 0.0_pr; CD_T= 0.0_pr
      case("NONE")
!!$       write(*,*)'ALL coefficients are taken into account'
      case default
         stop 'I do not understand whic set of coefficients C you want!'
      end select


    end subroutine from_functional


  end subroutine set_coeff_functional




end module asym_param
