module wi

  use asym_cstes ! module constants
  use asym_gauss ! module to integrate
  use asym_param
  use beta

contains
  
  !
  !!
  !

  subroutine w1central(q,rhoN,rhoP,w1) ! [ MeVfm^3]
    !=============================================
    ! abstract: function W_1^(0,1)=W_
    ! input
    ! q: transfer momentum in fm^-1
    ! rhoN,rhoP: densities of neutrons and protons
    ! output
    ! w1: in the four Spin Isospin channels
    !=============================================
    implicit none
    real ( kind = pr ),intent(in) :: q,rhoN,rhoP
    real ( kind = pr )            :: q2,t12,w1(0:1,0:1)
    real ( kind = pr )            :: rho1,rhosig,rhosig2
    real ( kind = pr )            :: rhosigm2,rhosig2m2

    rho1=rhoN-rhoP
    rhosig=(rhoN+rhoP)**(sigma )
    rhosig2=(rhoN+rhoP)**(sigma2)
    rhosigm2=(rhoN+rhoP)**(sigma-2 )
    rhosig2m2=(rhoN+rhoP)**(sigma2-2)
    s11=sigma+1
    s12=sigma+2
    s21=sigma2+1
    s22=sigma2+2
    q2=q*q
    t12=1.0_pr/2.0_pr

    w1(0,0) =2.0_pr*c_rho(0,1) + c_rho(0,2) * s11 * s12 * rhosig  & ! 2 body
         + c_rho(0,3) * s21 * s22 * rhosig2                       & ! 2 body
         + sigma *(sigma -1)*c_rho(1,2)*rhosigm2 *rho1**2         & ! 2 body
         + sigma2*(sigma2-1)*c_rho(1,3)*rhosig2m2*rho1**2         & ! 2 body
         - q2 * ( 2 * c_dr(0) + t12 * c_tau(0) )                    ! 2 body

    w1(0,1) = 2.0_pr* c_rho(1,1)                                      & ! 2 body
         + c_rho(1,2) * 2.0_pr* rhosig + c_rho(1,3) * 2 * rhosig2     & ! 2 body
         - q2 * ( 2 * c_dr(1) + t12 * c_tau(1) )                    ! 2 body

    w1(1,0) = 2 * c_s(0,1)                                        & ! 2 body
         + 2 * c_s(0,2) * rhosig + 2 * c_s(0,3) * rhosig2         & ! 2 body
         - q2 * ( 2 * c_ds(0) + t12*c_t(0) )                        ! 2 body

    w1(1,1) = 2 * c_s(1,1)                                        & ! 2 body
         + 2 * c_s(1,2) * rhosig + 2 * c_s(1,3) * rhosig2         & ! 2 body
         - q2 * ( 2 * c_ds(1) + 0.5_pr * c_t(1) )                   ! 2 body

    w1(0,0) = 4 * w1(0,0) ! See definition of W1 and W2 on Van Giai paper
    w1(0,1) = 4 * w1(0,1)
    w1(1,0) = 4 * w1(1,0)
    w1(1,1) = 4 * w1(1,1)  


    return
  end subroutine w1central
  
  !
  !!
  !

  subroutine w1qq(spin,proj,q,rhoN,rhoP,w1nn,w1pp,w1pn,w1np) ! W1 [ MeVfm^3]
    !=============================================
    ! abstract: function W_1^nn=W_
    ! input
    ! spin,proj: S,M
    ! q: transfer momentum in fm^-1
    ! rhoN,rhoP: densities of neutrons and protons
    ! output
    ! w1nn: function W_1^nn
    ! w1pp: function W_1^pp
    ! w1pn: function W_1^pn
    ! w1np: function W_1^np
    !=============================================

    implicit none
    integer, intent (in)          :: spin,proj
    integer                       :: deltaM
    real ( kind = pr ),intent(in) :: q,rhoN,rhoP
    real ( kind = pr ),intent(out):: w1nn,w1pp,w1pn,w1np
    real ( kind = pr )            :: q2,t12,w1(0:1,0:1)
    real ( kind = pr )            :: rho1

    rho1=rhoN-rhoP
    
    q2=q*q
    t12=1.0_pr/2.0_pr

    call w1central(q,rhoN,rhoP,w1) 
   
   
    if(spin.eq.0)then

       w1nn=t12*(w1(0,0)+w1(0,1)) &
            +8.0_pr*rho1*(sigma*(rhoN+rhoP)**(sigma-1)*c_rho(1,2) &
                         +sigma2*(rhoN+rhoP)**(sigma2-1)*c_rho(1,3))
       w1pp=t12*(w1(0,0)+w1(0,1)) &
            -8.0_pr*rho1*(sigma*(rhoN+rhoP)**(sigma-1)*c_rho(1,2)+sigma2*(rhoN+rhoP)**(sigma2-1)*c_rho(1,3))
       w1pn=t12*(w1(0,0)-w1(0,1))
       w1np=w1pn

    endif
    !
    if(spin.eq.1)then

       deltaM=0
       if(proj.eq.0)deltaM=1

       w1nn=t12*(w1(1,0)+w1(1,1)) &
            +deltaM*q2*(4.0_pr*(c_ns(0)+c_ns(1)) - (c_f(0)+c_f(1)))
       w1pp=w1nn
       w1pn=t12*(w1(1,0)-w1(1,1)) &
            +deltaM*q2*(4.0_pr*(c_ns(0)-c_ns(1)) - (c_f(0)-c_f(1)))
       w1np=w1pn

    endif


    return
  end subroutine w1qq

  !
  !!
  !
  subroutine w2qq(spin,w2nn,w2pp,w2pn,w2np)  ! W2  [ MeVfm^5]
    !=============================================
    ! abstract: function W_2^nn=W_
    ! input
    ! spin: S
    ! output
    ! w2nn: function W_2^nn
    ! w2pp: function W_2^pp
    ! w2pn: function W_2^pn
    ! w2np: function W_2^np
    !=============================================

    implicit none
    integer, intent (in)          :: spin
    real ( kind = pr ),intent(out):: w2nn,w2pp,w2pn,w2np
    real ( kind = pr )            :: t12

    t12=1.0_pr/2.0_pr


    if (spin.eq.0)then

       w2nn=t12*(4.0_pr*C_tau(0)+4.0_pr*C_tau(1)) 
       w2pp=w2nn
       w2pn=t12*(4.0_pr*C_tau(0)-4.0_pr*C_tau(1)) 
       w2np=w2pn

    endif
    if (spin.eq.1)then

       w2nn=t12*(4.0_pr*C_T(0)+4.0_pr*c_T(1)) 
       w2pp=w2nn
       w2pn=t12*(4.0_pr*C_T(0)-4.0_pr*c_T(1)) 
       w2np=w2pn

    endif

  
    return
  end subroutine w2qq


  !
  !!
  !

  subroutine matrix_nochex(Spin,proj,T,rhoN,rhoP,mu,MstarN,MstarP,qMev,omega,Amatrix,Bmatrix,Amatrixp,Bmatrixp,size,asymptotic)
    !---------------------------------------------
    ! abstract: the routine produces the matrices
    ! and vector column used to calculate the
    ! response functions via Cramer's method
    !
    ! input
    ! Spin,proj: S,M
    ! T: temperature
    ! rhoN,rhoP: densitites (neutron/protons)
    ! mstarN,mstarP: effective masses (neutron/proton)
    ! q: transfer momentum
    ! omega: tranfer energy 
    !
    ! output
    ! Amatrix: square matrix containing the interaction
    ! Bmatrix: column matrix containing the source
    ! size: the size of the matrix
    !---------------------------------------------
    implicit none
    integer, intent (in)          :: spin,proj
    logical, intent (in), optional:: asymptotic
    integer                       :: i1,i2,size
    real ( kind = pr ),intent(in) :: T,rhoN,rhoP,qMev,omega,mu(0:1),mstarN,mstarP
    real ( kind = pr )            :: c_fp,c_fm,c_njp,c_njm
    complex ( kind = pr )         :: betann(0:8),betapp(0:8)
    real ( kind = pr )            :: kfn,kfp,q
    complex ( kind = pr )         :: w1tildeNN,w1tildePP,w1tildeNP,w1tildePN
    real ( kind = pr )            :: w1nn(0:1,0:1),w1pp(0:1,0:1),w1pn(0:1,0:1),w1np(0:1,0:1)
    real ( kind = pr )            :: mnn,mpp,mpn,mnp,q2,q4,q3
    real ( kind = pr )            :: nun,nup,kn,kp,nkn,nkp
    real ( kind = pr )            :: w2nn(0:1),w2pp(0:1),w2pn(0:1),w2np(0:1)
    complex ( kind = pr )         :: xnn,xpp,ynp,ypn
    complex ( kind = pr )         :: znn(0:1),zpn(0:1),znp(0:1),zpp(0:1),ZZ(0:1)
    complex ( kind = pr )         :: Apn_p,Apn_m,Anp_p,Anp_m,Bpn_p,Bpn_m,Bnp_p,Bnp_m
    complex ( kind = pr )         :: alpha0n(1:7),alpha1n(1:7),alpha0p(1:7),alpha1p(1:7)
    complex ( kind = pr )         :: auxnn,auxpp
    complex ( kind = pr )         :: spn,smn,spp,smp,SSn,SSp,ZZp
    complex ( kind = pr )         :: Abpn_p,Abpn_m,Abnp_p,Abnp_m,Bbpn_p,Bbpn_m,Bbnp_p,Bbnp_m
    complex ( kind = pr )         :: Amatrix(1:8,1:8),Bmatrix(1:8)
    complex ( kind = pr )         :: Amatrixp(1:8,1:8),Bmatrixp(1:8)


    kfn=tpi13* rhoN**t13  ! [fm^-1]
    kfp=tpi13* rhoP**t13  ! [fm^-1]

    nun= (omega*mstarN)/(qMev*kfN*hbar_c)
    kn= qMev/hbar_c / (2.0_pr*kfN)  
    nup= (omega*mstarP)/(qMev*kfP*hbar_c)
    kp= qMev/hbar_c / (2.0_pr*kfP)  

   ! write(*,*)mu,mstarN,mstarP,rhon,rhop
    nkn=nun/kn-1.0_pr
    nkp=nup/kp-1.0_pr

    Amatrix(:,:)=Cmplx(0.0_pr,0.0_pr, kind =pr)
    Bmatrix(:)=Cmplx(0.0_pr,0.0_pr, kind =pr)
    Amatrixp(:,:)=Cmplx(0.0_pr,0.0_pr, kind =pr)
    Bmatrixp(:)=Cmplx(0.0_pr,0.0_pr, kind =pr)


    alpha0n(:)=Cmplx(0.0_pr,0.0_pr, kind =pr)
    alpha1n(:)=Cmplx(0.0_pr,0.0_pr, kind =pr)
    alpha0p(:)=Cmplx(0.0_pr,0.0_pr, kind =pr)
    alpha1p(:)=Cmplx(0.0_pr,0.0_pr, kind =pr)




        !write(*,*)mu


    if(.not.asymptotic)then 

       call betannall(omega,mstarN,qMev,T,mu(0),kfn,betann)
       call betannall(omega,mstarP,qMev,T,mu(1),kfp,betapp)

    else

       call betannallasymptotic(mstarN,qMev,kfn,betann)
       call betannallasymptotic(mstarP,qMev,kfp,betapp)

    endif
    !---- here i introduce again the unit of measure 

    q=qMev/hbar_c  ! [fm^-1]
    q2=q*q      ! [fm^-2]
    q3=q2*q     ! [fm^-3]
    q4=q2*q2    ! [fm^-4]
    betann(:)=betann(:)/hbarc3 ! [MeV^-1fm^-3]
    betapp(:)=betapp(:)/hbarc3

    c_fp=c_f(0)+c_f(1)   ! [MeV fm^5]
    c_fm=c_f(0)-c_f(1)   ! [MeV fm^5]
    c_njp=c_nj(0)+c_nj(1)   ! [MeV fm^5]
    c_njm=c_nj(0)-c_nj(1)   ! [MeV fm^5]




    do i1=0,1
       do i2=0,i1
          call w1qq(i1,i2,q,rhoN,rhoP,mnn,mpp,mpn,mnp) ! W1 [ MeVfm^3]
          w1nn(i1,i2)=mnn 
          w1np(i1,i2)=mnp
          w1pn(i1,i2)=mpn
          w1pp(i1,i2)=mpp
       enddo
       call w2qq(i1,mnn,mpp,mpn,mnp) ! W2  [ MeVfm^5]
       w2nn(i1)=mnn
       w2np(i1)=mnp
       w2pn(i1)=mpn
       w2pp(i1)=mpp
    enddo

    !
    !  HF response function (no Vph)
    ! 
    if(spin.eq.-1)then

       size=1 ! size of the matrix 

       ! B matrix (source)

       Bmatrix(1)=betann(0)   ! nn and np excitations
       Bmatrixp(1)=betapp(0)  ! pp and pn excitations

       Amatrix(1,1)=1.0_pr
       Amatrixp(1,1)=1.0_pr 

    endif

    !
    !------ Spin=0 M=0 channel
    !
    if(spin.eq.0)then

       size=6 ! size of the matrix 

       xpp=q2*(betapp(2)-betapp(3))*(w2pp(1)-c_fp) ! pure number
       xnn=q2*(betann(2)-betann(3))*(w2nn(1)-c_fp) ! pure number
       ynp=q2*(betann(2)-betann(3))*(w2np(1)-c_fm) ! pure number
       ypn=q2*(betapp(2)-betapp(3))*(w2pn(1)-c_fm) ! pure number

       auxnn=4.0_pr*q4/(1.0_pr+xnn-ynp*ypn/(1.0_pr+xpp))
       auxpp=4.0_pr*q4/(1.0_pr+xpp-ynp*ypn/(1.0_pr+xnn))

       ! the units are [ MeVfm^3] for the \tilde{W1} 
       w1tildeNN=w1nn(0,0)-auxnn*(c_njp*c_njm*(betapp(2)-betapp(3))*ynp/(1.0_pr+xpp) - (betann(2)-betann(3))*c_njp**2 ) &
            -auxpp*(c_njp*c_njm*(betann(2)-betann(3))*ypn/(1.0_pr+xnn) - (betapp(2)-betapp(3))*c_njm**2 )
       w1tildePP=w1pp(0,0)-auxpp*(c_njp*c_njm*(betann(2)-betann(3))*ypn/(1.0_pr+xnn) - (betapp(2)-betapp(3))*c_njp**2 ) &
            -auxnn*(c_njp*c_njm*(betapp(2)-betapp(3))*ynp/(1.0_pr+xpp) - (betann(2)-betann(3))*c_njm**2 )
       w1tildeNP=w1np(0,0)-auxnn*(c_njp**2*(betapp(2)-betapp(3))*ynp/(1.0_pr+xpp) - (betann(2)-betann(3))*c_njp*c_njm ) &
            -auxpp*(c_njm**2*(betann(2)-betann(3))*ypn/(1.0_pr+xnn) - (betapp(2)-betapp(3))*c_njp*c_njm )
       w1tildePN=w1pn(0,0)-auxpp*(c_njp**2*(betann(2)-betann(3))*ypn/(1.0_pr+xnn) - (betapp(2)-betapp(3))*c_njp*c_njm ) &
            -auxnn*(c_njm**2*(betapp(2)-betapp(3))*ynp/(1.0_pr+xpp) - (betann(2)-betann(3))*c_njp*c_njm )

       !    B matrix (source)

       Bmatrix(1)=betann(0)      ! nn and np excitations
       Bmatrix(2)=q2*betann(2)
       Bmatrix(3)=q*betann(1)

       Bmatrixp(1)=betapp(0)      ! pp and pn excitations
       Bmatrixp(2)=q2*betapp(2)
       Bmatrixp(3)=q*betapp(1)

       !     A matrix (interaction)      nn and np excitations

       Amatrix(1,1)=1.0_pr-betann(0)*w1tildeNN-q2*betann(2)*w2nn(0)
       Amatrix(1,2)=-betann(0)*w2nn(0)
       Amatrix(1,3)=2.0_pr*q*betann(1)*w2nn(0)
       Amatrix(1,4)=-betann(0)*w1tildenp-q2*betann(2)*w2np(0)
       Amatrix(1,5)=-betann(0)*w2np(0)
       Amatrix(1,6)=2.0_pr*q*betann(1)*w2np(0)
       !
       Amatrix(2,1)=-q2*betann(2)*w1tildeNN-q4*betann(5)*w2nn(0)
       Amatrix(2,2)=1.0_pr-q2*betann(2)*w2nn(0)
       Amatrix(2,3)=2.0_pr*q3*betann(4)*w2nn(0)
       Amatrix(2,4)=-q2*betann(2)*w1tildenp-q4*betann(5)*w2np(0)
       Amatrix(2,5)=-q2*betann(2)*w2np(0)
       Amatrix(2,6)=2.0_pr*q3*betann(4)*w2np(0)
       !
       Amatrix(3,1)=-q*betann(1)*w1tildeNN-q3*betann(4)*w2nn(0)
       Amatrix(3,2)=-q*betann(1)*w2nn(0)
       Amatrix(3,3)=1.0_pr+2.0_pr*q2*betann(3)*w2nn(0)
       Amatrix(3,4)=-q*betann(1)*w1tildenp-q3*betann(4)*w2np(0)
       Amatrix(3,5)=-q*betann(1)*w2np(0)
       Amatrix(3,6)=2.0_pr*q2*betann(3)*w2np(0)
       !
       Amatrix(4,1)=-betapp(0)*w1tildePN-q2*betapp(2)*w2pn(0)
       Amatrix(4,2)=-betapp(0)*w2pn(0)
       Amatrix(4,3)=2.0_pr*q*betapp(1)*w2pn(0)
       Amatrix(4,4)=1.0_pr-betapp(0)*w1tildePP-q2*betapp(2)*w2pp(0)
       Amatrix(4,5)=-betapp(0)*w2pp(0)
       Amatrix(4,6)=2.0_pr*q*betapp(1)*w2pp(0)
       !
       Amatrix(5,1)=-q2*betapp(2)*w1tildePN-q4*betapp(5)*w2pn(0)
       Amatrix(5,2)=-q2*betapp(2)*w2pn(0)
       Amatrix(5,3)=2.0_pr*q3*betapp(4)*w2pn(0)
       Amatrix(5,4)=-q2*betapp(2)*w1tildePP-q4*betapp(5)*w2pp(0)
       Amatrix(5,5)=1.0_pr-q2*betapp(2)*w2pp(0)
       Amatrix(5,6)=2.0_pr*q3*betapp(4)*w2pp(0)
       !
       Amatrix(6,1)=-q*betapp(1)*w1tildePN-q3*betapp(4)*w2pn(0)
       Amatrix(6,2)=-q*betapp(1)*w2pn(0)
       Amatrix(6,3)=2.0_pr*q2*betapp(3)*w2pn(0)
       Amatrix(6,4)=-q*betapp(1)*w1tildePP-q3*betapp(4)*w2pp(0)
       Amatrix(6,5)=-q*betapp(1)*w2pp(0)
       Amatrix(6,6)=1.0_pr+2.0_pr*q2*betapp(3)*w2pp(0)

       ! pp and pn excitations

       Amatrixp(1,1)=1.0_pr-betapp(0)*w1tildePP-q2*betapp(2)*w2pp(0)
       Amatrixp(1,2)=-betapp(0)*w2pp(0)
       Amatrixp(1,3)=2.0_pr*q*betapp(1)*w2pp(0)
       Amatrixp(1,4)=-betapp(0)*w1tildepn-q2*betapp(2)*w2pn(0)
       Amatrixp(1,5)=-betapp(0)*w2pn(0)
       Amatrixp(1,6)=2.0_pr*q*betapp(1)*w2pn(0)
       !
       Amatrixp(2,1)=-q2*betapp(2)*w1tildePP-q4*betapp(5)*w2pp(0)
       Amatrixp(2,2)=1.0_pr-q2*betapp(2)*w2pp(0)
       Amatrixp(2,3)=2.0_pr*q3*betapp(4)*w2pp(0)
       Amatrixp(2,4)=-q2*betapp(2)*w1tildepn-q4*betapp(5)*w2pn(0)
       Amatrixp(2,5)=-q2*betapp(2)*w2pn(0)
       Amatrixp(2,6)=2.0_pr*q3*betapp(4)*w2pn(0)
       !
       Amatrixp(3,1)=-q*betapp(1)*w1tildePP-q3*betapp(4)*w2pp(0)
       Amatrixp(3,2)=-q*betapp(1)*w2pp(0)
       Amatrixp(3,3)=1.0_pr+2.0_pr*q2*betapp(3)*w2pp(0)
       Amatrixp(3,4)=-q*betapp(1)*w1tildepn-q3*betapp(4)*w2pn(0)
       Amatrixp(3,5)=-q*betapp(1)*w2pn(0)
       Amatrixp(3,6)=2.0_pr*q2*betapp(3)*w2pn(0)
       !
       Amatrixp(4,1)=-betann(0)*w1tildeNP-q2*betann(2)*w2np(0)
       Amatrixp(4,2)=-betann(0)*w2np(0)
       Amatrixp(4,3)=2.0_pr*q*betann(1)*w2np(0)
       Amatrixp(4,4)=1.0_pr-betann(0)*w1tildeNN-q2*betann(2)*w2nn(0)
       Amatrixp(4,5)=-betann(0)*w2nn(0)
       Amatrixp(4,6)=2.0_pr*q*betann(1)*w2nn(0)
       !
       Amatrixp(5,1)=-q2*betann(2)*w1tildeNP-q4*betann(5)*w2np(0)
       Amatrixp(5,2)=-q2*betann(2)*w2np(0)
       Amatrixp(5,3)=2.0_pr*q3*betann(4)*w2np(0)
       Amatrixp(5,4)=-q2*betann(2)*w1tildeNN-q4*betann(5)*w2nn(0)
       Amatrixp(5,5)=1.0_pr-q2*betann(2)*w2nn(0)
       Amatrixp(5,6)=2.0_pr*q3*betann(4)*w2nn(0)
       !
       Amatrixp(6,1)=-q*betann(1)*w1tildeNP-q3*betann(4)*w2np(0)
       Amatrixp(6,2)=-q*betann(1)*w2np(0)
       Amatrixp(6,3)=2.0_pr*q2*betann(3)*w2np(0)
       Amatrixp(6,4)=-q*betann(1)*w1tildeNN-q3*betann(4)*w2nn(0)
       Amatrixp(6,5)=-q*betann(1)*w2nn(0)
       Amatrixp(6,6)=1.0_pr+2.0_pr*q2*betann(3)*w2nn(0)

    endif

    !
    !------ Spin=1 M=1 channel
    !

    if(spin.eq.1.and.proj.eq.1)then

       size=8 ! size of the matrix 

       znn(:)=q2*(betann(2)-betann(3))*w2nn(:)
       znp(:)=q2*(betann(2)-betann(3))*w2np(:)
       zpn(:)=q2*(betapp(2)-betapp(3))*w2pn(:)
       zpp(:)=q2*(betapp(2)-betapp(3))*w2pp(:)

       ZZ(:)=1.0_pr/((1.0_pr+znn(:))*(1.0_pr+zpp(:))-zpn(:)*znp(:))

       Apn_p=q2*ZZ(1)*( (1.0_pr+zpp(1))*c_fp*(betann(2)-betann(3)) - c_fm*znp(1)*(betapp(2)-betapp(3)))
       Apn_m=q2*ZZ(1)*( (1.0_pr+zpp(1))*c_fm*(betann(2)-betann(3)) - c_fp*znp(1)*(betapp(2)-betapp(3)))
       Anp_p=q2*ZZ(1)*( (1.0_pr+znn(1))*c_fp*(betapp(2)-betapp(3)) - c_fm*zpn(1)*(betann(2)-betann(3)))
       Anp_m=q2*ZZ(1)*( (1.0_pr+znn(1))*c_fm*(betapp(2)-betapp(3)) - c_fp*zpn(1)*(betann(2)-betann(3)))

       Bpn_p=q2*ZZ(1)*( (1.0_pr+zpp(1))*c_fp*(betann(4)-betann(6)) - c_fm*znp(1)*(betapp(4)-betapp(6)))
       Bpn_m=q2*ZZ(1)*( (1.0_pr+zpp(1))*c_fm*(betann(4)-betann(6)) - c_fp*znp(1)*(betapp(4)-betapp(6)))
       Bnp_p=q2*ZZ(1)*( (1.0_pr+znn(1))*c_fp*(betapp(4)-betapp(6)) - c_fm*zpn(1)*(betann(4)-betann(6)))
       Bnp_m=q2*ZZ(1)*( (1.0_pr+znn(1))*c_fm*(betapp(4)-betapp(6)) - c_fp*zpn(1)*(betann(4)-betann(6)))


       !

       w1tildenn=w1nn(1,1)-2*q4*ZZ(0)*c_njp*(znp(0)*c_njm*(betapp(2)-betapp(3))-(1.0_pr+zpp(0))*c_njp*(betann(2)-betann(3))) &
            -2*q4*ZZ(0)*c_njm*(zpn(0)*c_njp*(betann(2)-betann(3))-(1.0_pr+znn(0))*c_njm*(betapp(2)-betapp(3))) &
            +c_fp**2*q4*(betann(5)-betann(7))+c_fm**2*q4*(betapp(5)-betapp(7))      &
            -q2*nkn*(c_fp*(znn(1)*Bpn_p+znp(1)*Bnp_m))-q2*nkp*(c_fm*(zpp(1)*Bnp_m+zpn(1)*Bpn_p))
       w1tildenp=w1np(1,1)-2*q4*ZZ(0)*c_njp*(znp(0)*c_njp*(betapp(2)-betapp(3))-(1.0_pr+zpp(0))*c_njm*(betann(2)-betann(3))) &
            -2*q4*ZZ(0)*c_njm*(zpn(0)*c_njm*(betann(2)-betann(3))-(1.0_pr+znn(0))*c_njp*(betapp(2)-betapp(3))) &
            +c_fp*c_fm*q4*((betann(5)-betann(7))+(betapp(5)-betapp(7)))      &
            -q2*nkn*(c_fp*(znn(1)*Bpn_m+znp(1)*Bnp_p))-q2*nkp*(c_fm*(zpp(1)*Bnp_p+zpn(1)*Bpn_m))

       w1tildepp=w1pp(1,1)-2*q4*ZZ(0)*c_njp*(zpn(0)*c_njm*(betann(2)-betann(3))-(1.0_pr+znn(0))*c_njp*(betapp(2)-betapp(3))) &
            -2*q4*ZZ(0)*c_njm*(znp(0)*c_njp*(betapp(2)-betapp(3))-(1.0_pr+zpp(0))*c_njm*(betann(2)-betann(3))) &
            +c_fp**2*q4*(betapp(5)-betapp(7))+c_fm**2*q4*(betann(5)-betann(7))      &
            -q2*nkp*(c_fp*(zpp(1)*Bnp_p+zpn(1)*Bpn_m))-q2*nkn*(c_fm*(znn(1)*Bpn_m+znp(1)*Bnp_p))
       w1tildepn=w1pn(1,1)-2*q4*ZZ(0)*c_njp*(zpn(0)*c_njp*(betann(2)-betann(3))-(1.0_pr+znn(0))*c_njm*(betapp(2)-betapp(3))) &
            -2*q4*ZZ(0)*c_njm*(znp(0)*c_njm*(betapp(2)-betapp(3))-(1.0_pr+zpp(0))*c_njp*(betann(2)-betann(3))) &
            +c_fp*c_fm*q4*((betann(5)-betann(7))+(betapp(5)-betapp(7)))      &
            -q2*nkp*(c_fp*(zpp(1)*Bnp_m+zpn(1)*Bpn_p))-q2*nkn*(c_fm*(znn(1)*Bpn_p+znp(1)*Bnp_m))

       !

       alpha1n(1)=2.0_pr*q2*(c_fp*Bpn_p+c_fm*Bnp_m)
       alpha1p(1)=2.0_pr*q2*(c_fp*Bnp_m+c_fm*Bpn_p)
       alpha1n(3)=2.0_pr*q*(Apn_p*c_fp+c_fm*Anp_m)-2*q*w2nn(1)
       alpha1p(3)=2.0_pr*q*(Apn_p*c_fm+c_fp*Anp_m)-2*q*w2pn(1)
       alpha0n(3)=q*nkn*(c_fp*(znn(1)*Apn_p+znp(1)*Anp_m))+q*nkp*(c_fm*(zpp(1)*Anp_m+zpn(1)*Apn_p)) &
            -2*q3*(c_fp**2*(betann(4)-betann(6))+c_fm**2*(betapp(4)-betapp(6)) )
       alpha0p(3)=q*nkp*(c_fp*(zpp(1)*Anp_m+zpn(1)*Apn_p))+q*nkn*(c_fm*(znn(1)*Apn_p+znp(1)*Anp_m)) &
            -2*q3*(c_fp*c_fm)*((betann(4)-betann(6))+(betapp(4)-betapp(6)))

       alpha1n(5)=2*q2*(c_fp*Bpn_m+c_fm*Bnp_p)
       alpha1p(5)=2*q2*(c_fp*Bnp_p+c_fm*Bpn_m)
       alpha1n(7)=2*q*(Anp_p*c_fm+c_fp*Apn_m)-2*q*w2np(1)
       alpha1p(7)=2*q*(Anp_p*c_fp+c_fm*Apn_m)-2*q*w2pp(1)

       alpha0n(7)=q*nkn*(c_fp*(znn(1)*Apn_m+znp(1)*Anp_p))+q*nkp*(c_fm*(zpp(1)*Anp_p+zpn(1)*Apn_m)) &
            -2*q3*(c_fp*c_fm)*((betann(4)-betann(6))+(betapp(4)-betapp(6)) )
       alpha0p(7)=q*nkp*(c_fp*(zpp(1)*Anp_p+zpn(1)*Apn_m))+q*nkn*(c_fm*(znn(1)*Apn_m+znp(1)*Anp_p)) &
            -2*q3*(c_fp**2*(betapp(4)-betapp(6))+c_fm**2*(betann(4)-betann(6)) )


       !
       !   B matrix (source)
       !
       Bmatrix(1)=betann(0) ! nn and np exicitations
       Bmatrix(2)=q2*betann(2)
       Bmatrix(3)=q*betann(1)
       Bmatrix(4)=q2*(betann(2)-betann(3))

       Bmatrixp(1)=betapp(0) ! nn and np exicitations
       Bmatrixp(2)=q2*betapp(2)
       Bmatrixp(3)=q*betapp(1)
       Bmatrixp(4)=q2*(betapp(2)-betapp(3))

       !
       !  A matrix (interaction)
       !
       !          nn and pn excitations

       Amatrix(1,1)=1.0_pr-betann(0)*w1tildenn-q2*betann(2)*w2nn(1)-c_fp*q2*(betann(2)-betann(3))+betann(1)*alpha1n(1)
       Amatrix(1,2)=-betann(0)*w2nn(1)
       Amatrix(1,3)=-betann(1)*alpha1n(3)-betann(0)*alpha0n(3)
       Amatrix(1,4)=-c_fp*betann(0)
       Amatrix(1,5)=-betann(0)*w1tildenp-q2*betann(2)*w2np(1)-c_fm*q2*(betann(2)-betann(3))+alpha1n(5)*betann(1)
       Amatrix(1,6)=-betann(0)*w2np(1)
       Amatrix(1,7)=-alpha1n(7)*betann(1)-alpha0n(7)*betann(0)
       Amatrix(1,8)=-c_fm*betann(0)
       !
       Amatrix(2,1)=-q2*betann(2)*w1tildenn-q4*betann(5)*w2nn(1)-c_fp*q4*(betann(5)-betann(8))+q2*betann(4)*alpha1n(1)
       Amatrix(2,2)=1.0_pr-q2*betann(2)*w2nn(1)
       Amatrix(2,3)=-q2*betann(4)*alpha1n(3)-q2*betann(2)*alpha0n(3)
       Amatrix(2,4)=-q2*c_fp*betann(2)
       Amatrix(2,5)=-q2*betann(2)*w1tildenp-q4*betann(5)*w2np(1)-c_fm*q4*(betann(5)-betann(8))+q2*alpha1n(5)*betann(4)
       Amatrix(2,6)=-q2*betann(2)*w2np(1)
       Amatrix(2,7)=-q2*alpha1n(7)*betann(4)-q2*alpha0n(7)*betann(2)
       Amatrix(2,8)=-q2*c_fm*betann(2)
       !
       Amatrix(3,1)=-q*betann(1)*w1tildenn-q3*betann(4)*w2nn(1)-c_fp*q3*(betann(4)-betann(6))+q*betann(3)*alpha1n(1)
       Amatrix(3,2)=-q*betann(1)*w2nn(1)
       Amatrix(3,3)=1.0_pr-q*betann(3)*alpha1n(3)-q*betann(1)*alpha0n(3)
       Amatrix(3,4)=-q*c_fp*betann(1)
       Amatrix(3,5)=-q*betann(1)*w1tildenp-q3*betann(4)*w2np(1)-c_fm*q3*(betann(4)-betann(6))+q*alpha1n(5)*betann(3)
       Amatrix(3,6)=-q*betann(1)*w2np(1)
       Amatrix(3,7)=-q*alpha1n(7)*betann(3)-q*alpha0n(7)*betann(1)
       Amatrix(3,8)=-q*c_fm*betann(1)   
       !
       Amatrix(4,1)=-q2*(betann(2)-betann(3))*w1tildenn-q4*(betann(5)-betann(8))*w2nn(1) &
            -c_fp*q4*(betann(5)-2*betann(8)+betann(7))+q2*(betann(4)-betann(6))*alpha1n(1)
       Amatrix(4,2)=-q2*(betann(2)-betann(3))*w2nn(1)
       Amatrix(4,3)=-q2*(betann(4)-betann(6))*alpha1n(3)-q2*(betann(2)-betann(3))*alpha0n(3)
       Amatrix(4,4)=1.0_pr-c_fp*q2*(betann(2)-betann(3))
       Amatrix(4,5)=-q2*(betann(2)-betann(3))*w1tildenp-q4*(betann(5)-betann(8))*w2np(1) &
            -c_fm*q4*(betann(5)-2*betann(8)+betann(7))+alpha1n(5)*q2*(betann(4)-betann(6))
       Amatrix(4,6)=-q2*(betann(2)-betann(3))*w2np(1)
       Amatrix(4,7)=-q2*alpha1n(7)*(betann(4)-betann(6))-alpha0n(7)*q2*(betann(2)-betann(3))
       Amatrix(4,8)=-c_fm*q2*(betann(2)-betann(3))
       !
       Amatrix(5,1)=-betapp(0)*w1tildepn-q2*betapp(2)*w2pn(1)-c_fm*q2*(betapp(2)-betapp(3))+betapp(1)*alpha1p(1)
       Amatrix(5,2)=-betapp(0)*w2pn(1)
       Amatrix(5,3)=-betapp(1)*alpha1p(3)-betapp(0)*alpha0p(3)
       Amatrix(5,4)=-c_fm*betapp(0)
       Amatrix(5,5)=1.0_pr-betapp(0)*w1tildepp-q2*betapp(2)*w2pp(1)-c_fp*q2*(betapp(2)-betapp(3))+alpha1p(5)*betapp(1)
       Amatrix(5,6)=-betapp(0)*w2pp(1)
       Amatrix(5,7)=-alpha1p(7)*betapp(1)-alpha0p(7)*betapp(0)
       Amatrix(5,8)=-c_fp*betapp(0)
       !
       Amatrix(6,1)=-q2*betapp(2)*w1tildepn-q4*betapp(5)*w2pn(1)-c_fm*q4*(betapp(5)-betapp(8))+q2*betapp(4)*alpha1p(1)
       Amatrix(6,2)=-q2*betapp(2)*w2pn(1)
       Amatrix(6,3)=-q2*betapp(4)*alpha1p(3)-q2*betapp(2)*alpha0p(3)
       Amatrix(6,4)=-q2*c_fm*betapp(2)
       Amatrix(6,5)=-q2*betapp(2)*w1tildepp-q4*betapp(5)*w2pp(1)-c_fp*q4*(betapp(5)-betapp(8))+q2*alpha1p(5)*betapp(4)
       Amatrix(6,6)=1.0_pr-q2*betapp(2)*w2pp(1)
       Amatrix(6,7)=-q2*alpha1p(7)*betapp(4)-q2*alpha0p(7)*betapp(2)
       Amatrix(6,8)=-q2*c_fp*betapp(2)
       !
       Amatrix(7,1)=-q*betapp(1)*w1tildepn-q3*betapp(4)*w2pn(1)-c_fm*q3*(betapp(4)-betapp(6))+q*betapp(3)*alpha1p(1)
       Amatrix(7,2)=-q*betapp(1)*w2pn(1)
       Amatrix(7,3)=-q*betapp(3)*alpha1p(3)-q*betapp(1)*alpha0p(3)
       Amatrix(7,4)=-q*c_fm*betapp(1)
       Amatrix(7,5)=-q*betapp(1)*w1tildepp-q3*betapp(4)*w2pp(1)-c_fp*q3*(betapp(4)-betapp(6))+q*alpha1p(5)*betapp(3)
       Amatrix(7,6)=-q*betapp(1)*w2pp(1)
       Amatrix(7,7)=1.0_pr-q*alpha1p(7)*betapp(3)-q*alpha0p(7)*betapp(1)
       Amatrix(7,8)=-q*c_fp*betapp(1)   
       !
       Amatrix(8,1)=-q2*(betapp(2)-betapp(3))*w1tildepn-q4*(betapp(5)-betapp(8))*w2pn(1) &
            -c_fm*q4*(betapp(5)-2*betapp(8)+betapp(7))+q2*(betapp(4)-betapp(6))*alpha1p(1)
       Amatrix(8,2)=-q2*(betapp(2)-betapp(3))*w2pn(1)
       Amatrix(8,3)=-q2*(betapp(4)-betapp(6))*alpha1p(3)-q2*(betapp(2)-betapp(3))*alpha0p(3)
       Amatrix(8,4)=-c_fm*q2*(betapp(2)-betapp(3))
       Amatrix(8,5)=-q2*(betapp(2)-betapp(3))*w1tildepp-q4*(betapp(5)-betapp(8))*w2pp(1) &
            -c_fp*q4*(betapp(5)-2*betapp(8)+betapp(7))+alpha1p(5)*q2*(betapp(4)-betapp(6))
       Amatrix(8,6)=-q2*(betapp(2)-betapp(3))*w2pp(1)
       Amatrix(8,7)=-q2*alpha1p(7)*(betapp(4)-betapp(6))-alpha0p(7)*q2*(betapp(2)-betapp(3))
       Amatrix(8,8)=1.0_pr-c_fp*q2*(betapp(2)-betapp(3))

       !   pp and np excitations

       Amatrixp(1,1)=1.0_pr-betapp(0)*w1tildepp-q2*betapp(2)*w2pp(1)-c_fp*q2*(betapp(2)-betapp(3))+betapp(1)*alpha1p(5)
       Amatrixp(1,2)=-betapp(0)*w2pp(1)
       Amatrixp(1,3)=-betapp(1)*alpha1p(7)-betapp(0)*alpha0p(7)
       Amatrixp(1,4)=-c_fp*betapp(0)
       Amatrixp(1,5)=-betapp(0)*w1tildepn-q2*betapp(2)*w2pn(1)-c_fm*q2*(betapp(2)-betapp(3))+alpha1p(1)*betapp(1)
       Amatrixp(1,6)=-betapp(0)*w2pn(1)
       Amatrixp(1,7)=-alpha1p(3)*betapp(1)-alpha0p(3)*betapp(0)
       Amatrixp(1,8)=-c_fm*betapp(0)
       !
       Amatrixp(2,1)=-q2*betapp(2)*w1tildepp-q4*betapp(5)*w2pp(1)-c_fp*q4*(betapp(5)-betapp(8))+q2*betapp(4)*alpha1p(5)
       Amatrixp(2,2)=1.0_pr-q2*betapp(2)*w2pp(1)
       Amatrixp(2,3)=-q2*betapp(4)*alpha1p(7)-q2*betapp(2)*alpha0p(7)
       Amatrixp(2,4)=-q2*c_fp*betapp(2)
       Amatrixp(2,5)=-q2*betapp(2)*w1tildepn-q4*betapp(5)*w2pn(1)-c_fm*q4*(betapp(5)-betapp(8))+q2*alpha1p(1)*betapp(4)
       Amatrixp(2,6)=-q2*betapp(2)*w2pn(1)
       Amatrixp(2,7)=-q2*alpha1p(3)*betapp(4)-q2*alpha0p(3)*betapp(2)
       Amatrixp(2,8)=-q2*c_fm*betapp(2)
       !
       Amatrixp(3,1)=-q*betapp(1)*w1tildepp-q3*betapp(4)*w2pp(1)-c_fp*q3*(betapp(4)-betapp(6))+q*betapp(3)*alpha1p(5)
       Amatrixp(3,2)=-q*betapp(1)*w2pp(1)
       Amatrixp(3,3)=1.0_pr-q*betapp(3)*alpha1p(7)-q*betapp(1)*alpha0p(7)
       Amatrixp(3,4)=-q*c_fp*betapp(1)
       Amatrixp(3,5)=-q*betapp(1)*w1tildepn-q3*betapp(4)*w2pn(1)-c_fm*q3*(betapp(4)-betapp(6))+q*alpha1p(1)*betapp(3)
       Amatrixp(3,6)=-q*betapp(1)*w2pn(1)
       Amatrixp(3,7)=-q*alpha1p(3)*betapp(3)-q*alpha0p(3)*betapp(1)
       Amatrixp(3,8)=-q*c_fm*betapp(1)   
       !
       Amatrixp(4,1)=-q2*(betapp(2)-betapp(3))*w1tildepp-q4*(betapp(5)-betapp(8))*w2pp(1) &
            -c_fp*q4*(betapp(5)-2*betapp(8)+betapp(7))+q2*(betapp(4)-betapp(6))*alpha1p(5)
       Amatrixp(4,2)=-q2*(betapp(2)-betapp(3))*w2pp(1)
       Amatrixp(4,3)=-q2*(betapp(4)-betapp(6))*alpha1p(7)-q2*(betapp(2)-betapp(3))*alpha0p(7)
       Amatrixp(4,4)=1.0_pr-c_fp*q2*(betapp(2)-betapp(3))
       Amatrixp(4,5)=-q2*(betapp(2)-betapp(3))*w1tildepn-q4*(betapp(5)-betapp(8))*w2pn(1) &
            -c_fm*q4*(betapp(5)-2*betapp(8)+betapp(7))+alpha1p(1)*q2*(betapp(4)-betapp(6))
       Amatrixp(4,6)=-q2*(betapp(2)-betapp(3))*w2pn(1)
       Amatrixp(4,7)=-q2*alpha1p(3)*(betapp(4)-betapp(6))-alpha0p(3)*q2*(betapp(2)-betapp(3))
       Amatrixp(4,8)=-c_fm*q2*(betapp(2)-betapp(3))
       !
       Amatrixp(5,1)=-betann(0)*w1tildenp-q2*betann(2)*w2np(1)-c_fm*q2*(betann(2)-betann(3))+betann(1)*alpha1n(5)
       Amatrixp(5,2)=-betann(0)*w2np(1)
       Amatrixp(5,3)=-betann(1)*alpha1n(7)-betann(0)*alpha0n(7)
       Amatrixp(5,4)=-c_fm*betann(0)
       Amatrixp(5,5)=1.0_pr-betann(0)*w1tildenn-q2*betann(2)*w2nn(1)-c_fp*q2*(betann(2)-betann(3))+alpha1n(1)*betann(1)
       Amatrixp(5,6)=-betann(0)*w2nn(1)
       Amatrixp(5,7)=-alpha1n(3)*betann(1)-alpha0n(3)*betann(0)
       Amatrixp(5,8)=-c_fp*betann(0)
       !
       Amatrixp(6,1)=-q2*betann(2)*w1tildenp-q4*betann(5)*w2np(1)-c_fm*q4*(betann(5)-betann(8))+q2*betann(4)*alpha1n(5)
       Amatrixp(6,2)=-q2*betann(2)*w2np(1)
       Amatrixp(6,3)=-q2*betann(4)*alpha1n(7)-q2*betann(2)*alpha0n(7)
       Amatrixp(6,4)=-q2*c_fm*betann(2)
       Amatrixp(6,5)=-q2*betann(2)*w1tildenn-q4*betann(5)*w2nn(1)-c_fp*q4*(betann(5)-betann(8))+q2*alpha1n(1)*betann(4)
       Amatrixp(6,6)=1.0_pr-q2*betann(2)*w2nn(1)
       Amatrixp(6,7)=-q2*alpha1n(3)*betann(4)-q2*alpha0n(3)*betann(2)
       Amatrixp(6,8)=-q2*c_fp*betann(2)
       !
       Amatrixp(7,1)=-q*betann(1)*w1tildenp-q3*betann(4)*w2np(1)-c_fm*q3*(betann(4)-betann(6))+q*betann(3)*alpha1n(5)
       Amatrixp(7,2)=-q*betann(1)*w2np(1)
       Amatrixp(7,3)=-q*betann(3)*alpha1n(7)-q*betann(1)*alpha0n(7)
       Amatrixp(7,4)=-q*c_fm*betann(1)
       Amatrixp(7,5)=-q*betann(1)*w1tildenn-q3*betann(4)*w2nn(1)-c_fp*q3*(betann(4)-betann(6))+q*alpha1n(1)*betann(3)
       Amatrixp(7,6)=-q*betann(1)*w2nn(1)
       Amatrixp(7,7)=1.0_pr-q*alpha1n(3)*betann(3)-q*alpha0n(3)*betann(1)
       Amatrixp(7,8)=-q*c_fp*betann(1)   
       !
       Amatrixp(8,1)=-q2*(betann(2)-betann(3))*w1tildenp-q4*(betann(5)-betann(8))*w2np(1) &
            -c_fm*q4*(betann(5)-2*betann(8)+betann(7))+q2*(betann(4)-betann(6))*alpha1n(5)
       Amatrixp(8,2)=-q2*(betann(2)-betann(3))*w2np(1)
       Amatrixp(8,3)=-q2*(betann(4)-betann(6))*alpha1n(7)-q2*(betann(2)-betann(3))*alpha0n(7)
       Amatrixp(8,4)=-c_fm*q2*(betann(2)-betann(3))
       Amatrixp(8,5)=-q2*(betann(2)-betann(3))*w1tildenn-q4*(betann(5)-betann(8))*w2nn(1) &
            -c_fp*q4*(betann(5)-2*betann(8)+betann(7))+alpha1n(1)*q2*(betann(4)-betann(6))
       Amatrixp(8,6)=-q2*(betann(2)-betann(3))*w2nn(1)
       Amatrixp(8,7)=-q2*alpha1n(3)*(betann(4)-betann(6))-alpha0n(3)*q2*(betann(2)-betann(3))
       Amatrixp(8,8)=1.0_pr-c_fp*q2*(betann(2)-betann(3))
    endif

    !
    !------ Spin=1 M=0 channel
    ! 

    if(spin.eq.1.and.proj.eq.0)then
       size=8 ! size of the matrix 

       znn(:)=q2*(betann(2)-betann(3))*w2nn(:)
       znp(:)=q2*(betann(2)-betann(3))*w2np(:)
       zpn(:)=q2*(betapp(2)-betapp(3))*w2pn(:)
       zpp(:)=q2*(betapp(2)-betapp(3))*w2pp(:)

       spn=c_fp*q2*(betann(2)-betann(3))    
       smn=c_fm*q2*(betann(2)-betann(3))    
       spp=c_fp*q2*(betapp(2)-betapp(3))    
       smp=c_fm*q2*(betapp(2)-betapp(3))  

       SSn=1.0_pr+znn(1)+3.0_pr*spn-(znp(1)+3.0_pr*smn)*(zpn(1)+3.0_pr*smp)/(1.0_pr+zpp(1)+3.0_pr*spp)    
       SSp=1.0_pr+zpp(1)+3.0_pr*spp-(zpn(1)+3.0_pr*smp)*(znp(1)+3.0_pr*smn)/(1.0_pr+znn(1)+3.0_pr*spn)    

       ZZp=1.0_pr/((1.0_pr+zpp(1)+3.0_pr*spp)*(1.0_pr+znn(1)+3.0_pr*spn)-(znp(1)+3.0_pr*smn)*(zpn(1)+3.0_pr*smp))



       Abpn_p=ZZp*( (1.0_pr+zpp(1)+3.0_pr*spp)*spn - (znp(1)+3.0_pr*smn)*smp)
       Abpn_m=ZZp*( (1.0_pr+zpp(1)+3.0_pr*spp)*smn - (znp(1)+3.0_pr*smn)*spp)
       Abnp_p=ZZp*( (1.0_pr+znn(1)+3.0_pr*spn)*spp - (zpn(1)+3.0_pr*smp)*smn)
       Abnp_m=ZZp*( (1.0_pr+znn(1)+3.0_pr*spn)*smp - (zpn(1)+3.0_pr*smp)*spn)

       Bbpn_p=ZZp*(nkn*(1.0_pr+zpp(1)+3.0_pr*spp)*spn - nkp*(znp(1)+3.0_pr*smn)*smp)
       Bbpn_m=ZZp*(nkn*(1.0_pr+zpp(1)+3.0_pr*spp)*smn - nkp*(znp(1)+3.0_pr*smn)*spp)
       Bbnp_p=ZZp*(nkp*(1.0_pr+znn(1)+3.0_pr*spn)*spp - nkn*(zpn(1)+3.0_pr*smp)*smn) 
       Bbnp_m=ZZp*(nkp*(1.0_pr+znn(1)+3.0_pr*spn)*smp - nkn*(zpn(1)+3.0_pr*smp)*spn) 




       w1tildenn=w1nn(1,0)-q2*c_fp*nkn*(Bbpn_p*(znn(1)+3.0_pr*spn)+Bbnp_m*(znp(1)+3.0_pr*smn)) &
            -q2*c_fm*nkp*(Bbpn_p*(zpn(1)+3.0_pr*smp)+Bbnp_m*(zpp(1)+3.0_pr*spp)) &
            +4.0_pr*q4*(c_fp**2*(betann(8)-betann(7))+c_fm**2*(betapp(8)-betapp(7)))

       w1tildenp=w1np(1,0)-q2*c_fp*nkn*(Bbnp_p*(znp(1)+3.0_pr*smn)+Bbpn_m*(znn(1)+3.0_pr*spn)) &
            -q2*c_fm*nkp*(Bbnp_p*(zpp(1)+3.0_pr*spp)+Bbpn_m*(zpn(1)+3.0_pr*smp)) &
            +4.0_pr*q4*c_fp*c_fm*((betapp(8)-betapp(7))+(betann(8)-betann(7)))

       w1tildepp=w1pp(1,0)-q2*c_fp*nkp*(Bbnp_p*(zpp(1)+3.0_pr*spp)+Bbpn_m*(zpn(1)+3.0_pr*smp)) &
            -q2*c_fm*nkn*(Bbnp_p*(znp(1)+3.0_pr*smn)+Bbpn_m*(znn(1)+3.0_pr*spn)) &
            +4.0_pr*q4*(c_fp**2*(betapp(8)-betapp(7))+c_fm**2*(betann(8)-betann(7)))

       w1tildepn=w1pn(1,0)-q2*c_fp*nkp*(Bbpn_p*(zpn(1)+3.0_pr*smp)+Bbnp_m*(zpp(1)+3.0_pr*spp)) &
            -q2*c_fm*nkn*(Bbpn_p*(znn(1)+3.0_pr*spn)+Bbnp_m*(znp(1)+3.0_pr*smn)) &
            +4.0_pr*q4*c_fp*c_fm*((betapp(8)-betapp(7))+(betann(8)-betann(7)))



       alpha1n(1)=2*q2*(c_fp*Bbpn_p+c_fm*Bbnp_m)
       alpha1p(1)=2*q2*(c_fp*Bbnp_m+c_fm*Bbpn_p)

       alpha1n(3)=2*q*(-w2nn(1)-2.0_pr*c_fp+2.0_pr*Abpn_p*c_fp+2.0_pr*Abnp_m*c_fm)
       alpha1p(3)=2*q*(-w2pn(1)-2.0_pr*c_fm+2.0_pr*Abnp_m*c_fp+2.0_pr*Abpn_p*c_fm)

       alpha0n(3)= 2*q*nkn*c_fp*(Abpn_p*(znn(1)+3.0_pr*spn)+Abnp_m*(znp(1)+3.0_pr*smn)-spn) &
            +2*q*nkp*c_fm*(Abnp_m*(zpp(1)+3.0_pr*spp)+Abpn_p*(zpn(1)+3.0_pr*smp)-smp)
       alpha0p(3)= 2*q*nkp*c_fp*(Abnp_m*(zpp(1)+3.0_pr*spp)+Abpn_p*(zpn(1)+3.0_pr*smp)-smp) &
            +2*q*nkn*c_fm*(Abpn_p*(znn(1)+3.0_pr*spn)+Abnp_m*(znp(1)+3.0_pr*smn)-spn)

       alpha1n(5)=2*q2*(c_fp*Bbpn_m+c_fm*bbnp_p)
       alpha1p(5)=2*q2*(c_fp*Bbnp_p+c_fm*bbpn_m)

       alpha1n(7)=2*q*(-w2np(1)-2.0_pr*c_fm+2.0_pr*Abpn_m*c_fp+2.0_pr*Abnp_p*c_fm)
       alpha1p(7)=2*q*(-w2pp(1)-2.0_pr*c_fp+2.0_pr*Abnp_p*c_fp+2.0_pr*Abpn_m*c_fm)

       alpha0n(7)= 2*q*nkn*c_fp*(Abpn_m*(znn(1)+3.0_pr*spn)+Abnp_p*(znp(1)+3.0_pr*smn)-smn) &
            +2*q*nkp*c_fm*(Abnp_p*(zpp(1)+3.0_pr*spp)+Abpn_m*(zpn(1)+3.0_pr*smp)-spp)
       alpha0p(7)= 2*q*nkp*c_fp*(Abnp_p*(zpp(1)+3.0_pr*spp)+Abpn_m*(zpn(1)+3.0_pr*smp)-spp) &
            +2*q*nkn*c_fm*(Abpn_m*(znn(1)+3.0_pr*spn)+Abnp_p*(znp(1)+3.0_pr*smn)-smn)



       !
       !   B matrix (source)
       !
       Bmatrix(1)=betann(0) ! nn and np exicitations
       Bmatrix(2)=q2*betann(2)
       Bmatrix(3)=q*betann(1)
       Bmatrix(4)=q2*betann(3)

       Bmatrixp(1)=betapp(0) ! nn and np exicitations
       Bmatrixp(2)=q2*betapp(2)
       Bmatrixp(3)=q*betapp(1)
       Bmatrixp(4)=q2*betapp(3)

       !
       !  A matrix (interaction)
       !
       !          nn and pn excitations

       Amatrix(1,1)=1.0_pr-betann(0)*w1tildenn-q2*betann(2)*w2nn(1)-2.0_pr*c_fp*q2*betann(3)+betann(1)*alpha1n(1)
       Amatrix(1,2)=-betann(0)*w2nn(1)
       Amatrix(1,3)=-betann(1)*alpha1n(3)-betann(0)*alpha0n(3)
       Amatrix(1,4)=-c_fp*betann(0)
       Amatrix(1,5)=-betann(0)*w1tildenp-q2*betann(2)*w2np(1)-2.0_pr*c_fm*q2*betann(3)+alpha1n(5)*betann(1)
       Amatrix(1,6)=-betann(0)*w2np(1)
       Amatrix(1,7)=-alpha1n(7)*betann(1)-alpha0n(7)*betann(0)
       Amatrix(1,8)=-c_fm*betann(0)
       !
       Amatrix(2,1)=-q2*betann(2)*w1tildenn-q4*betann(5)*w2nn(1)-2.0_pr*c_fp*q4*betann(8)+q2*betann(4)*alpha1n(1)
       Amatrix(2,2)=1.0_pr-q2*betann(2)*w2nn(1)
       Amatrix(2,3)=-q2*betann(4)*alpha1n(3)-q2*betann(2)*alpha0n(3)
       Amatrix(2,4)=-q2*c_fp*betann(2)
       Amatrix(2,5)=-q2*betann(2)*w1tildenp-q4*betann(5)*w2np(1)-2.0_pr*c_fm*q4*betann(8)+q2*alpha1n(5)*betann(4)
       Amatrix(2,6)=-q2*betann(2)*w2np(1)
       Amatrix(2,7)=-q2*alpha1n(7)*betann(4)-q2*alpha0n(7)*betann(2)
       Amatrix(2,8)=-q2*c_fm*betann(2)
       !
       Amatrix(3,1)=-q*betann(1)*w1tildenn-q3*betann(4)*w2nn(1)-2.0_pr*c_fp*q3*betann(6)+q*betann(3)*alpha1n(1)
       Amatrix(3,2)=-q*betann(1)*w2nn(1)
       Amatrix(3,3)=1.0_pr-q*betann(3)*alpha1n(3)-q*betann(1)*alpha0n(3)
       Amatrix(3,4)=-q*c_fp*betann(1)
       Amatrix(3,5)=-q*betann(1)*w1tildenp-q3*betann(4)*w2np(1)-2.0_pr*c_fm*q3*betann(6)+q*alpha1n(5)*betann(3)
       Amatrix(3,6)=-q*betann(1)*w2np(1)
       Amatrix(3,7)=-q*alpha1n(7)*betann(3)-q*alpha0n(7)*betann(1)
       Amatrix(3,8)=-q*c_fm*betann(1)   
       !
       Amatrix(4,1)=-q2*betann(3)*w1tildenn-q4*betann(8)*w2nn(1) &
            -2.0_pr*c_fp*q4*betann(7)+q2*betann(6)*alpha1n(1)
       Amatrix(4,2)=-q2*betann(3)*w2nn(1)
       Amatrix(4,3)=-q2*betann(6)*alpha1n(3)-q2*betann(3)*alpha0n(3)
       Amatrix(4,4)=1.0_pr/2.0_pr-c_fp*q2*betann(3)
       Amatrix(4,5)=-q2*betann(3)*w1tildenp-q4*betann(8)*w2np(1) &
            -2.0_pr*c_fm*q4*betann(7)+alpha1n(5)*q2*betann(6)
       Amatrix(4,6)=-q2*betann(3)*w2np(1)
       Amatrix(4,7)=-q2*alpha1n(7)*betann(6)-alpha0n(7)*q2*betann(3)
       Amatrix(4,8)=-c_fm*q2*betann(3) 
       !
       Amatrix(5,1)=-betapp(0)*w1tildepn-q2*betapp(2)*w2pn(1)-2.0_pr*c_fm*q2*betapp(3)+betapp(1)*alpha1p(1)
       Amatrix(5,2)=-betapp(0)*w2pn(1)
       Amatrix(5,3)=-betapp(1)*alpha1p(3)-betapp(0)*alpha0p(3)
       Amatrix(5,4)=-c_fm*betapp(0)
       Amatrix(5,5)=1.0_pr-betapp(0)*w1tildepp-q2*betapp(2)*w2pp(1)-2.0_pr*c_fp*q2*betapp(3)+alpha1p(5)*betapp(1)
       Amatrix(5,6)=-betapp(0)*w2pp(1)
       Amatrix(5,7)=-alpha1p(7)*betapp(1)-alpha0p(7)*betapp(0)
       Amatrix(5,8)=-c_fp*betapp(0)
       !
       Amatrix(6,1)=-q2*betapp(2)*w1tildepn-q4*betapp(5)*w2pn(1)-2.0_pr*c_fm*q4*betapp(8)+q2*betapp(4)*alpha1p(1)
       Amatrix(6,2)=-q2*betapp(2)*w2pn(1)
       Amatrix(6,3)=-q2*betapp(4)*alpha1p(3)-q2*betapp(2)*alpha0p(3)
       Amatrix(6,4)=-q2*c_fm*betapp(2)
       Amatrix(6,5)=-q2*betapp(2)*w1tildepp-q4*betapp(5)*w2pp(1)-2.0_pr*c_fp*q4*betapp(8)+q2*alpha1p(5)*betapp(4)
       Amatrix(6,6)=1.0_pr-q2*betapp(2)*w2pp(1)
       Amatrix(6,7)=-q2*alpha1p(7)*betapp(4)-q2*alpha0p(7)*betapp(2)
       Amatrix(6,8)=-q2*c_fp*betapp(2)
       !
       Amatrix(7,1)=-q*betapp(1)*w1tildepn-q3*betapp(4)*w2pn(1)-2.0_pr*c_fm*q3*betapp(6)+q*betapp(3)*alpha1p(1)
       Amatrix(7,2)=-q*betapp(1)*w2pn(1)
       Amatrix(7,3)=-q*betapp(3)*alpha1p(3)-q*betapp(1)*alpha0p(3)
       Amatrix(7,4)=-q*c_fm*betapp(1)
       Amatrix(7,5)=-q*betapp(1)*w1tildepp-q3*betapp(4)*w2pp(1)-2.0_pr*c_fp*q3*betapp(6)+q*alpha1p(5)*betapp(3)
       Amatrix(7,6)=-q*betapp(1)*w2pp(1)
       Amatrix(7,7)=1.0_pr-q*alpha1p(7)*betapp(3)-q*alpha0p(7)*betapp(1)
       Amatrix(7,8)=-q*c_fp*betapp(1)   
       !
       Amatrix(8,1)=-q2*betapp(3)*w1tildepn-q4*betapp(8)*w2pn(1) &
            -2.0_pr*c_fm*q4*betapp(7)+q2*betapp(6)*alpha1p(1)
       Amatrix(8,2)=-q2*betapp(3)*w2pn(1)
       Amatrix(8,3)=-q2*betapp(6)*alpha1p(3)-q2*betapp(3)*alpha0p(3)
       Amatrix(8,4)=-c_fm*q2*betapp(3)
       Amatrix(8,5)=-q2*betapp(3)*w1tildepp-q4*betapp(8)*w2pp(1) &
            -2.0_pr*c_fp*q4*betapp(7)+alpha1p(5)*q2*betapp(6)
       Amatrix(8,6)=-q2*betapp(3)*w2pp(1)
       Amatrix(8,7)=-q2*alpha1p(7)*betapp(6)-alpha0p(7)*q2*betapp(3)
       Amatrix(8,8)=1.0_pr/2.0_pr-c_fp*q2*betapp(3)

       !          pp and np excitations

       Amatrixp(1,1)=1.0_pr-betapp(0)*w1tildepp-q2*betapp(2)*w2pp(1)-2.0_pr*c_fp*q2*betapp(3)+betapp(1)*alpha1p(5)
       Amatrixp(1,2)=-betapp(0)*w2pp(1)
       Amatrixp(1,3)=-betapp(1)*alpha1p(7)-betapp(0)*alpha0p(7)
       Amatrixp(1,4)=-c_fp*betapp(0)
       Amatrixp(1,5)=-betapp(0)*w1tildepn-q2*betapp(2)*w2pn(1)-2.0_pr*c_fm*q2*betapp(3)+alpha1p(1)*betapp(1)
       Amatrixp(1,6)=-betapp(0)*w2pn(1)
       Amatrixp(1,7)=-alpha1p(3)*betapp(1)-alpha0p(3)*betapp(0)
       Amatrixp(1,8)=-c_fm*betapp(0)
       !
       Amatrixp(2,1)=-q2*betapp(2)*w1tildepp-q4*betapp(5)*w2pp(1)-2.0_pr*c_fp*q4*betapp(8)+q2*betapp(4)*alpha1p(5)
       Amatrixp(2,2)=1.0_pr-q2*betapp(2)*w2pp(1)
       Amatrixp(2,3)=-q2*betapp(4)*alpha1p(7)-q2*betapp(2)*alpha0p(7)
       Amatrixp(2,4)=-q2*c_fp*betapp(2)
       Amatrixp(2,5)=-q2*betapp(2)*w1tildepn-q4*betapp(5)*w2pn(1)-2.0_pr*c_fm*q4*betapp(8)+q2*alpha1p(1)*betapp(4)
       Amatrixp(2,6)=-q2*betapp(2)*w2pn(1)
       Amatrixp(2,7)=-q2*alpha1p(3)*betapp(4)-q2*alpha0p(3)*betapp(2)
       Amatrixp(2,8)=-q2*c_fm*betapp(2)
       ! 
       Amatrixp(3,1)=-q*betapp(1)*w1tildepp-q3*betapp(4)*w2pp(1)-2.0_pr*c_fp*q3*betapp(6)+q*betapp(3)*alpha1p(5)
       Amatrixp(3,2)=-q*betapp(1)*w2pp(1)
       Amatrixp(3,3)=1.0_pr-q*betapp(3)*alpha1p(7)-q*betapp(1)*alpha0p(7)
       Amatrixp(3,4)=-q*c_fp*betapp(1)
       Amatrixp(3,5)=-q*betapp(1)*w1tildepn-q3*betapp(4)*w2pn(1)-2.0_pr*c_fm*q3*betapp(6)+q*alpha1p(1)*betapp(3)
       Amatrixp(3,6)=-q*betapp(1)*w2pn(1)
       Amatrixp(3,7)=-q*alpha1p(3)*betapp(3)-q*alpha0p(3)*betapp(1)
       Amatrixp(3,8)=-q*c_fm*betapp(1)   
       !
       Amatrixp(4,1)=-q2*betapp(3)*w1tildepp-q4*betapp(8)*w2pp(1) &
            -2.0_pr*c_fp*q4*betapp(7)+q2*betapp(6)*alpha1p(5)
       Amatrixp(4,2)=-q2*betapp(3)*w2pp(1)
       Amatrixp(4,3)=-q2*betapp(6)*alpha1p(7)-q2*betapp(3)*alpha0p(7)
       Amatrixp(4,4)=1.0_pr/2.0_pr-c_fp*q2*betapp(3)
       Amatrixp(4,5)=-q2*betapp(3)*w1tildepn-q4*betapp(8)*w2pn(1) &
            -2.0_pr*c_fm*q4*betapp(7)+alpha1p(1)*q2*betapp(6)
       Amatrixp(4,6)=-q2*betapp(3)*w2pn(1)
       Amatrixp(4,7)=-q2*alpha1p(3)*betapp(6)-alpha0p(3)*q2*betapp(3)
       Amatrixp(4,8)=-c_fm*q2*betapp(3) 
       !
       Amatrixp(5,1)=-betann(0)*w1tildenp-q2*betann(2)*w2np(1)-2.0_pr*c_fm*q2*betann(3)+betann(1)*alpha1n(5)
       Amatrixp(5,2)=-betann(0)*w2np(1)
       Amatrixp(5,3)=-betann(1)*alpha1n(7)-betann(0)*alpha0n(7)
       Amatrixp(5,4)=-c_fm*betann(0)
       Amatrixp(5,5)=1.0_pr-betann(0)*w1tildenn-q2*betann(2)*w2nn(1)-2.0_pr*c_fp*q2*betann(3)+alpha1n(1)*betann(1)
       Amatrixp(5,6)=-betann(0)*w2nn(1)
       Amatrixp(5,7)=-alpha1n(3)*betann(1)-alpha0n(3)*betann(0)
       Amatrixp(5,8)=-c_fp*betann(0)
       !
       Amatrixp(6,1)=-q2*betann(2)*w1tildenp-q4*betann(5)*w2np(1)-2.0_pr*c_fm*q4*betann(8)+q2*betann(4)*alpha1n(5)
       Amatrixp(6,2)=-q2*betann(2)*w2np(1)
       Amatrixp(6,3)=-q2*betann(4)*alpha1n(7)-q2*betann(2)*alpha0n(7)
       Amatrixp(6,4)=-q2*c_fm*betann(2)
       Amatrixp(6,5)=-q2*betann(2)*w1tildenn-q4*betann(5)*w2nn(1)-2.0_pr*c_fp*q4*betann(8)+q2*alpha1n(1)*betann(4)
       Amatrixp(6,6)=1.0_pr-q2*betann(2)*w2nn(1)
       Amatrixp(6,7)=-q2*alpha1n(3)*betann(4)-q2*alpha0n(3)*betann(2)
       Amatrixp(6,8)=-q2*c_fp*betann(2)
       !
       Amatrixp(7,1)=-q*betann(1)*w1tildenp-q3*betann(4)*w2np(1)-2.0_pr*c_fm*q3*betann(6)+q*betann(3)*alpha1n(5)
       Amatrixp(7,2)=-q*betann(1)*w2np(1)
       Amatrixp(7,3)=-q*betann(3)*alpha1n(7)-q*betann(1)*alpha0n(7)
       Amatrixp(7,4)=-q*c_fm*betann(1)
       Amatrixp(7,5)=-q*betann(1)*w1tildenn-q3*betann(4)*w2nn(1)-2.0_pr*c_fp*q3*betann(6)+q*alpha1n(1)*betann(3)
       Amatrixp(7,6)=-q*betann(1)*w2nn(1)
       Amatrixp(7,7)=1.0_pr-q*alpha1n(3)*betann(3)-q*alpha0n(3)*betann(1)
       Amatrixp(7,8)=-q*c_fp*betann(1)   
       !
       Amatrixp(8,1)=-q2*betann(3)*w1tildenp-q4*betann(8)*w2np(1) &
            -2.0_pr*c_fm*q4*betann(7)+q2*betann(6)*alpha1n(5)
       Amatrixp(8,2)=-q2*betann(3)*w2np(1)
       Amatrixp(8,3)=-q2*betann(6)*alpha1n(7)-q2*betann(3)*alpha0n(7)
       Amatrixp(8,4)=-c_fm*q2*betann(3)
       Amatrixp(8,5)=-q2*betann(3)*w1tildenn-q4*betann(8)*w2nn(1) &
            -2.0_pr*c_fp*q4*betann(7)+alpha1n(1)*q2*betann(6)
       Amatrixp(8,6)=-q2*betann(3)*w2nn(1)
       Amatrixp(8,7)=-q2*alpha1n(3)*betann(6)-alpha0n(3)*q2*betann(3)
       Amatrixp(8,8)=1.0_pr/2.0_pr-c_fp*q2*betann(3)


    endif
    !---------------------------


    return
  end subroutine matrix_nochex

  !
  !!
  !
  subroutine matrix_chex(Spin,proj,T,rhoN,rhoP,Un,Up,mstarN,mstarP,qMev,omega, &
       Amatrix,Bmatrix,Amatrixp,Bmatrixp,size)
    !---------------------------------------------
    ! abstract: the routine produces the matrices
    ! and vector column used to calculate the
    ! response functions via Cramer's method
    !
    ! input
    ! Spin,proj: S,M
    ! T: temperature
    ! rhoN,rhoP: densitites (neutron/protons)
    ! mstarN,mstarP: effective masses (neutron/proton)
    ! q: transfer momentum
    ! omega: tranfer energy 
    !
    ! output
    ! Amatrix: square matrix containing the interaction
    ! Bmatrix: column matrix containing the source
    ! Amatrixp: square matrix containing the interaction
    ! Bmatrixp: column matrix containing the source
    ! size: the size of the matrix
    !---------------------------------------------
    implicit none
    integer, intent (in)          :: spin,proj
    integer                       :: size
    real ( kind = pr ),intent(in) :: T,rhoN,rhoP,mstarN,mstarP,qMev,omega,Un,Up
    complex ( kind = pr )         :: betanp(0:8),betapn(0:8)
    real ( kind = pr )            :: kfn,kfp,mu(0:1),q
    complex ( kind = pr )         :: w1tilde,znp(0:1),zbnp(0:1)
    complex ( kind = pr )         :: w1tildep,zpn(0:1),zbpn(0:1)
    real ( kind = pr )            :: w1(0:1,0:1),w2(0:1,0:1),kf,ef
    real ( kind = pr )            :: q2,q4,q3,q5,q6
    real ( kind = pr )            :: nun,nup,kn,kp,nkn,nkp
    complex ( kind = pr )         :: b2b3,b4b6,alpha11np,alpha30np,alpha31np
    complex ( kind = pr )         :: b2b3p,b4b6p,alpha11pn,alpha30pn,alpha31pn
    complex ( kind = pr )         :: Amatrix(1:8,1:8),Bmatrix(1:8)
    complex ( kind = pr )         :: Amatrixp(1:8,1:8),Bmatrixp(1:8)

    kfn=tpi13* rhoN**t13  ! [fm^-1]
    kfp=tpi13* rhoP**t13  ! [fm^-1]

    nun= (omega*mstarN)/(qMev*kfN*hbar_c)
    kn= (qMev/hbar_c) / (2.0_pr*kfN)  
    nup= (omega*mstarP)/(qMev*kfP*hbar_c)
    kp= (qMev/hbar_c) / (2.0_pr*kfP)  

    kf=tpi13*(( rhoN+rhoP)/2.0_pr)**t13
    ef=kf**2*(hbarc2/((mstarN+mstarP))) ! [MeV]

    nkn=nun/kn-1.0_pr
    nkp=nup/kp-1.0_pr

    Amatrix(:,:)=Cmplx(0.0_pr,0.0_pr, kind =pr)
    Bmatrix(:)=Cmplx(0.0_pr,0.0_pr, kind =pr)
    Amatrixp(:,:)=Cmplx(0.0_pr,0.0_pr, kind =pr)
    Bmatrixp(:)=Cmplx(0.0_pr,0.0_pr, kind =pr)

    q=qMeV/hbar_c
    q2=q*q
    q3=q2*q
    q4=q2*q2
    q5=q3*q2
    q6=q5*q
    call w1central(q,rhoN,rhoP,w1) ! W1 [ MeVfm^3]

    !---- w2 term -----------------------

    W2(0,0)=4.0_pr*C_tau(0)  ! W2 [ MeVfm^5]
    W2(0,1)=4.0_pr*C_tau(1)
    W2(1,0)=4.0_pr*C_T(0)
    W2(1,1)=4.0_pr*C_T(1)

    !------------------------------------

    !call chemical(T,Un,Up,rhoN,rhoP,mstarN,mstarP,mu)
mu=(0d0,0d0)
    !write(*,*)'ciao',T/ef
    if(T/ef.le.0.05_pr)then

       if(flagBnpanalytic)then ! full analytical beta_np
          call betapnTzero(omega,mstarN,mstarP,qMeV,Un,Up,kfn,kfp,betapn) ! np
          call betapnTzero(omega,mstarP,mstarN,qMeV,Up,Un,kfp,kfn,betanp) ! pn
       else! imaginary part analytical / real part integrated
          call betapnallNUMERICTZERO(omega,mstarN,mstarP,qMeV,Un,Up,kfn,kfp,betapn)
          call betapnallNUMERICTZERO(omega,mstarP,mstarN,qMeV,Up,Un,kfp,kfn,betanp)
       endif
       !       write(*,*)betanp(0),'ciccio'
    else
       call betapnall(omega,mstarN,mstarP,qMeV,T,mu(0),mu(1),Un,Up,kfn,kfp,betapn) ! np
       call betapnall(omega,mstarP,mstarN,qMeV,T,mu(1),mu(0),Up,Un,kfp,kfn,betanp) ! pn
    endif


    betanp(:)=betanp(:)/hbarc3! [MeV^-1 fm^-3]
    betapn(:)=betapn(:)/hbarc3


    if(spin.eq.-1)then ! HF response no Vph interaction
       size=1

       Bmatrix(1)=betanp(0)
       Bmatrixp(1)=betapn(0)
       !       write(*,*)'dentro',betanp(0),betapn(0)

       !   np excitations
       Amatrix(1,1)=Cmplx(1.0_pr,0.0_pr, kind =pr)
       Amatrixp(1,1)=Cmplx(1.0_pr,0.0_pr, kind =pr)

    endif

    !
    ! ------   S=0 channel
    !
    if(spin.eq.0)then
       size=3

       W1tilde=W1(0,1)+16.0_pr*q4*(c_nj(1)**2)*(betanp(2)-betanp(3))/ &
            (1.0_pr+q2*(betanp(2)-betanp(3))*(w2(1,1)-2.0_pr*c_f(1)))

       W1tildeP=W1(0,1)+16.0_pr*q4*(c_nj(1)**2)*(betapn(2)-betapn(3))/ &
            (1.0_pr+q2*(betapn(2)-betapn(3))*(w2(1,1)-2.0_pr*c_f(1)))
       !
       !  Bmatrix (source)
       !

       Bmatrix(1)=betanp(0)
       Bmatrix(2)=q2*betanp(2)
       Bmatrix(3)=q*betanp(1)

       Bmatrixp(1)=betapn(0)
       Bmatrixp(2)=q2*betapn(2)
       Bmatrixp(3)=q*betapn(1)

       !
       ! Amatrix (interaction)
       !
       !   np excitations
       Amatrix(1,1)=1.0_pr-betanp(0)*W1tilde-q2*betanp(2)*W2(0,1)
       Amatrix(1,2)=-betanp(0)*W2(0,1)
       Amatrix(1,3)=2.0_pr*q*betanp(1)*W2(0,1)
       !
       Amatrix(2,1)=-q2*betanp(2)*W1tilde-q4*betanp(5)*W2(0,1)
       Amatrix(2,2)=1.0_pr-q2*betanp(2)*W2(0,1)
       Amatrix(2,3)=2.0_pr*q3*betanp(4)*W2(0,1)  
       !
       Amatrix(3,1)=-q*betanp(1)*W1tilde-q3*betanp(4)*W2(0,1)
       Amatrix(3,2)=-q*betanp(1)*W2(0,1)
       Amatrix(3,3)=1.0_pr+2.0_pr*q2*betanp(3)*W2(0,1)

       ! pn excitations
       Amatrixp(1,1)=1.0_pr-betapn(0)*W1tildeP-q2*betapn(2)*W2(0,1)
       Amatrixp(1,2)=-betapn(0)*W2(0,1)
       Amatrixp(1,3)=2.0_pr*q*betapn(1)*W2(0,1)
       !
       Amatrixp(2,1)=-q2*betapn(2)*W1tildeP-q4*betapn(5)*W2(0,1)
       Amatrixp(2,2)=1.0_pr-q2*betapn(2)*W2(0,1)
       Amatrixp(2,3)=2.0_pr*q3*betapn(4)*W2(0,1)  
       !
       Amatrixp(3,1)=-q*betapn(1)*W1tildep-q3*betapn(4)*W2(0,1)
       Amatrixp(3,2)=-q*betapn(1)*W2(0,1)
       Amatrixp(3,3)=1.0_pr+2.0_pr*q2*betapn(3)*W2(0,1)


    endif

    if(spin.eq.1.and.proj.eq.1)then

       size=4
       znp(:)=W2(:,1)*q2*(betanp(2)-betanp(3))
       zpn(:)=W2(:,1)*q2*(betapn(2)-betapn(3))


       W1tilde=W1(1,1)+8.0_pr*q4*(c_nj(1)**2)/(1.0_pr+znp(0))*(betanp(2)-betanp(3)) &
            +q4*(c_f(1)**2)*(4.0_pr*(betanp(5)-betanp(7))                 &
            -8.0_pr*W2(1,1)*q2*(betanp(4)-betanp(6))**2/(1.0_pr+znp(1)))
       W1tildep=W1(1,1)+8.0_pr*q4*(c_nj(1)**2)/(1.0_pr+zpn(0))*(betapn(2)-betapn(3)) &
            +q4*(c_f(1)**2)*(4.0_pr*(betapn(5)-betapn(7))                 &
            -8.0_pr*W2(1,1)*q2*(betapn(4)-betapn(6))**2/(1.0_pr+zpn(1)))


       b4b6=betanp(4)-betanp(6)
       b2b3=betanp(2)-betanp(3)
       b4b6p=betapn(4)-betapn(6)
       b2b3p=betapn(2)-betapn(3)

       !
       !  Bmatrix
       !

       Bmatrix(1)=betanp(0)
       Bmatrix(2)=q2*betanp(2)
       Bmatrix(3)=q*betanp(1)
       Bmatrix(4)=q2*(betanp(2)-betanp(3))

       Bmatrixp(1)=betapn(0)
       Bmatrixp(2)=q2*betapn(2)
       Bmatrixp(3)=q*betapn(1)
       Bmatrixp(4)=q2*(betapn(2)-betapn(3))
       !
       !  Amatrix
       !
       !   np excitations
       Amatrix(1,1)=1.0_pr-betanp(0)*W1tilde-q2*betanp(2)*W2(1,1) &
            -2.0_pr*c_f(1)*q2*(betanp(2)-betanp(3))       &
            +8.0_pr*(c_f(1))**2*q4*betanp(1)*b4b6/(1.0_pr+znp(1))
       Amatrix(1,2)=-betanp(0)*W2(1,1)
       Amatrix(1,3)=2.0_pr*q*betanp(1)*W2(1,1)-8.0_pr*q3*(c_f(1)**2)* &
            (betanp(1)*b2b3-betanp(0)*b4b6)/(1.0_pr+znp(1))
       Amatrix(1,4)=-2.0_pr*c_f(1)*betanp(0)
       !
       Amatrix(2,1)=-q2*betanp(2)*W1tilde-q4*betanp(5)*W2(1,1) &
            -2.0_pr*c_f(1)*q4*(betanp(5)-betanp(8))       &
            +8.0_pr*(c_f(1))**2*q6*betanp(4)*b4b6/(1.0_pr+znp(1))
       Amatrix(2,2)=1.0_pr-q2*betanp(2)*W2(1,1)
       Amatrix(2,3)=2.0_pr*q3*betanp(4)*W2(1,1)+8.0_pr*q5*(c_f(1)**2)*&
            (betanp(4)*betanp(3)-betanp(2)*betanp(6))/(1.0_pr+znp(1))
       Amatrix(2,4)=-2.0_pr*q2*c_f(1)*betanp(2) 
       !
       Amatrix(3,1)=-q*betanp(1)*W1tilde-q3*betanp(4)*W2(1,1) &
            -2.0_pr*c_f(1)*q3*(betanp(4)-betanp(6))       &
            +8.0_pr*(c_f(1))**2*q5*betanp(3)*b4b6/(1.0_pr+znp(1))
       Amatrix(3,2)=-q*betanp(1)*W2(1,1)
       Amatrix(3,3)=1.0_pr+2.0_pr*q2*betanp(3)*W2(1,1)-8.0_pr*q4*(c_f(1)**2)*&
            (betanp(3)*b2b3-betanp(1)*b4b6)/(1.0_pr+znp(1))
       Amatrix(3,4)=-2.0_pr*q*c_f(1)*betanp(1)
       !
       Amatrix(4,1)=-q2*(betanp(2)-betanp(3))*W1tilde-q4*(betanp(5)-betanp(8))*W2(1,1) &
            -2.0_pr*c_f(1)*q4*(betanp(5)-2.0_pr*betanp(8)+betanp(7))       &
            +8.0_pr*(c_f(1))**2*q6*b4b6**2/(1.0_pr+znp(1))
       Amatrix(4,2)=-q2*(betanp(2)-betanp(3))*W2(1,1)
       Amatrix(4,3)=2.0_pr*q3*b4b6*W2(1,1)
       Amatrix(4,4)=1.0_pr-2.0_pr*q2*c_f(1)*(betanp(2)-betanp(3))

       !   pn excitations
       Amatrixp(1,1)=1.0_pr-betapn(0)*W1tildep-q2*betapn(2)*W2(1,1) &
            -2.0_pr*c_f(1)*q2*(betapn(2)-betapn(3))       &
            +8.0_pr*(c_f(1))**2*q4*betapn(1)*b4b6p/(1.0_pr+zpn(1))
       Amatrixp(1,2)=-betapn(0)*W2(1,1)
       Amatrixp(1,3)=2.0_pr*q*betapn(1)*W2(1,1)-8.0_pr*q3*(c_f(1)**2)* &
            (betapn(1)*b2b3p-betapn(0)*b4b6p)/(1.0_pr+zpn(1))
       Amatrixp(1,4)=-2.0_pr*c_f(1)*betapn(0)
       !
       Amatrixp(2,1)=-q2*betapn(2)*W1tildep-q4*betapn(5)*W2(1,1) &
            -2.0_pr*c_f(1)*q4*(betapn(5)-betapn(8))       &
            +8.0_pr*(c_f(1))**2*q6*betapn(4)*b4b6p/(1.0_pr+zpn(1))
       Amatrixp(2,2)=1.0_pr-q2*betapn(2)*W2(1,1)
       Amatrixp(2,3)=2.0_pr*q3*betapn(4)*W2(1,1)+8.0_pr*q5*(c_f(1)**2)*&
            (betapn(4)*betapn(3)-betapn(2)*betapn(6))/(1.0_pr+zpn(1))
       Amatrixp(2,4)=-2.0_pr*q2*c_f(1)*betapn(2) 
       !
       Amatrixp(3,1)=-q*betapn(1)*W1tildep-q3*betapn(4)*W2(1,1) &
            -2.0_pr*c_f(1)*q3*(betapn(4)-betapn(6))       &
            +8.0_pr*(c_f(1))**2*q5*betapn(3)*b4b6p/(1.0_pr+zpn(1))
       Amatrixp(3,2)=-q*betapn(1)*W2(1,1)
       Amatrixp(3,3)=1.0_pr+2.0_pr*q2*betapn(3)*W2(1,1)-8.0_pr*q4*(c_f(1)**2)*&
            (betapn(3)*b2b3p-betapn(1)*b4b6p)/(1.0_pr+zpn(1))
       Amatrixp(3,4)=-2.0_pr*q*c_f(1)*betapn(1)
       !
       Amatrixp(4,1)=-q2*(betapn(2)-betapn(3))*W1tildep-q4*(betapn(5)-betapn(8))*W2(1,1) &
            -2.0_pr*c_f(1)*q4*(betapn(5)-2.0_pr*betapn(8)+betapn(7))       &
            +8.0_pr*(c_f(1))**2*q6*b4b6p**2/(1.0_pr+zpn(1))
       Amatrixp(4,2)=-q2*(betapn(2)-betapn(3))*W2(1,1)
       Amatrixp(4,3)=2.0_pr*q3*(b4b6p*W2(1,1))
       Amatrixp(4,4)=1.0_pr-2.0_pr*q2*c_f(1)*(betapn(2)-betapn(3))

    endif


    if(spin.eq.1.and.proj.eq.0)then
       size=4
       b4b6=betanp(4)-betanp(6)
       b2b3=betanp(2)-betanp(3)

       zbnp(:)=(W2(:,1)+6.0_pr*c_f(1))*q2*b2b3


       W1tilde=W1(1,1)+q2*(8.0_pr*c_ns(1)-2.0_pr*c_f(1)) &
            +16.0_pr*(c_f(1)**2)*q3*(q*(betanp(8)-betanp(7))-q3*(W2(1,1)+6*c_f(1))*b4b6**2/(1.0_pr+zbnp(1)) )

       alpha11np=16.0_pr*(c_f(1)**2)*q3*b4b6/(1.0_pr+zbnp(1))
       alpha30np=alpha11np
       alpha31np=W2(1,1)+4.0_pr*c_f(1) &
            -8.0_pr*(c_f(1)**2)*q2*b2b3/(1.0_pr+zbnp(1))

       !
       b4b6p=betapn(4)-betapn(6)
       b2b3p=betapn(2)-betapn(3)

       zbpn(:)=(W2(:,1)+6.0_pr*c_f(1))*q2*b2b3p


       W1tildep=W1(1,1)+q2*(8.0_pr*c_ns(1)-2.0_pr*c_f(1)) &
            +16.0_pr*(c_f(1)**2)*q3*(q*(betapn(8)-betapn(7))-q3*(W2(1,1)+6*c_f(1))*b4b6p**2/(1.0_pr+zbpn(1)) )


       alpha11pn=16.0_pr*(c_f(1)**2)*q3*b4b6p/(1.0_pr+zbpn(1))
       alpha30pn=alpha11pn
       alpha31pn=W2(1,1)+4.0_pr*c_f(1) &
            -8.0_pr*(c_f(1)**2)*q2*b2b3p/(1.0_pr+zbpn(1))



       !
       !  Bmatrix
       !

       Bmatrix(1)=betanp(0)
       Bmatrix(2)=q2*betanp(2)
       Bmatrix(3)=q*betanp(1)
       Bmatrix(4)= q2*(betanp(2)-betanp(3))

       Bmatrixp(1)=betapn(0)
       Bmatrixp(2)=q2*betapn(2)
       Bmatrixp(3)=q*betapn(1)
       Bmatrixp(4)= q2*(betapn(2)-betapn(3))
       !
       !  Amatrix
       !

       Amatrix(1,1)=1.0_pr-betanp(0)*W1tilde-q2*betanp(2)*W2(1,1) &
            +q*betanp(1)*alpha11np-4.0_pr*c_f(1)*q2*betanp(3)
       Amatrix(1,2)=-betanp(0)*W2(1,1)
       Amatrix(1,3)=2.0_pr*q*betanp(1)*alpha31np+betanp(0)*alpha30np
       Amatrix(1,4)=-4.0_pr*c_f(1)*betanp(0)
       !
       Amatrix(2,1)=-q2*betanp(2)*W1tilde-q4*betanp(5)*W2(1,1) &
            +q3*betanp(4)*alpha11np-4.0_pr*c_f(1)*q4*betanp(8)
       Amatrix(2,2)=1.0_pr-q2*betanp(2)*W2(1,1)
       Amatrix(2,3)=2.0_pr*q3*betanp(4)*alpha31np+q2*betanp(2)*alpha30np
       Amatrix(2,4)=-4.0_pr*c_f(1)*q2*betanp(2)
       !
       Amatrix(3,1)=-q*betanp(1)*W1tilde-q3*betanp(4)*W2(1,1) &
            +q2*betanp(3)*alpha11np-4.0_pr*c_f(1)*q3*betanp(6)
       Amatrix(3,2)=-q*betanp(1)*W2(1,1)
       Amatrix(3,3)=1.0_pr+2.0_pr*q2*betanp(3)*alpha31np+q*betanp(1)*alpha30np
       Amatrix(3,4)=-4.0_pr*c_f(1)*q*betanp(1)
       !
       Amatrix(4,1)=-q2*betanp(3)*W1tilde-q4*betanp(8)*W2(1,1) &
            +q3*betanp(6)*alpha11np-4.0_pr*c_f(1)*q4*betanp(7)
       Amatrix(4,2)=-q2*betanp(3)*W2(1,1)
       Amatrix(4,3)=+2.0_pr*q3*betanp(6)*alpha31np+q2*betanp(3)*alpha30np
       Amatrix(4,4)=1.0_pr-4.0_pr*c_f(1)*q2*betanp(3)

       !
       !  pn excitations
       Amatrixp(1,1)=1.0_pr-betapn(0)*W1tildep-q2*betapn(2)*W2(1,1) &
            +q*betapn(1)*alpha11pn-4.0_pr*c_f(1)*q2*betapn(3)
       Amatrixp(1,2)=-betapn(0)*W2(1,1)
       Amatrixp(1,3)=2.0_pr*q*betapn(1)*alpha31pn+betapn(0)*alpha30pn
       Amatrixp(1,4)=-4.0_pr*c_f(1)*betapn(0)
       !
       Amatrixp(2,1)=-q2*betapn(2)*W1tildep-q4*betapn(5)*W2(1,1) &
            +q3*betapn(4)*alpha11pn-4.0_pr*c_f(1)*q4*betapn(8)
       Amatrixp(2,2)=1.0_pr-q2*betapn(2)*W2(1,1)
       Amatrixp(2,3)=2.0_pr*q3*betapn(4)*alpha31pn+q2*betapn(2)*alpha30pn
       Amatrixp(2,4)=-4.0_pr*c_f(1)*q2*betapn(2)
       !
       Amatrixp(3,1)=-q*betapn(1)*W1tildep-q3*betapn(4)*W2(1,1) &
            +q2*betapn(3)*alpha11pn-4.0_pr*c_f(1)*q3*betapn(6)
       Amatrixp(3,2)=-q*betapn(1)*W2(1,1)
       Amatrixp(3,3)=1.0_pr+2.0_pr*q2*betapn(3)*alpha31pn+q*betapn(1)*alpha30pn
       Amatrixp(3,4)=-4.0_pr*c_f(1)*q*betapn(1)
       !
       Amatrixp(4,1)=-q2*betapn(3)*W1tildep-q4*betapn(8)*W2(1,1) &
            +q3*betapn(6)*alpha11pn-4.0_pr*c_f(1)*q4*betapn(7)
       Amatrixp(4,2)=-q2*betapn(3)*W2(1,1)
       Amatrixp(4,3)=+2.0_pr*q3*betapn(6)*alpha31pn+q2*betapn(3)*alpha30pn
       Amatrixp(4,4)=1.0_pr-4.0_pr*c_f(1)*q2*betapn(3)




    endif



    return
  end subroutine matrix_chex

  !
  !!
end module wi
