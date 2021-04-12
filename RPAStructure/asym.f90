module asym

  use asym_cstes ! module constants
  use asym_gauss ! module to integrate
  use asym_param ! module c.c. constants of the functional
  use beta       ! module for the beta functions nn and np
  use wi         ! module for the A and B matrix
  use determinante
  use inm
  use omp_lib

contains

  subroutine asymmetric(rho,Y,TMeV,cross,energy)
    !-----------------------------------------
    !input
    !
    ! nom: name of the force
    ! rho: density of the system
    ! Y: asymmetry parameter
    ! TMev: temperature [MeV]
    !-----------------------------------------


    implicit none

    real ( kind = pr ),intent(in)     :: rho,y,TMeV,energy
    real ( kind = pr ),intent(out)     :: cross
    integer, parameter                :: itracemeff=0
    integer, parameter                :: imax=1500
    integer                           :: size,Spin,proj,iw,id,id2,idmin,id3,idimnn(0:3),idnnre(0:3),i
    real ( kind = pr )                :: domega,factor,Tsuref
    real ( kind = pr )                :: rhoP,rhoN,mstarN,mstarP,omega,T,q
    real ( kind = pr )                :: UpotN,UpotP,mu(0:1),rho_sat,ef,kf
    real ( kind = pr )                :: kfn,kfp,Emax
    real ( kind = pr )                :: qmin(0:1,0:1,0:1),dmin(0:1,0:1,0:1)
    complex ( kind = pr )             :: Amatrix(1:8,1:8),Bmatrix(1:8)
    complex ( kind = pr )             :: Amatrixp(1:8,1:8),Bmatrixp(1:8)
    complex ( kind = pr )             :: Chinn,Chinp,Chipn,Chipp
    complex ( kind = pr ),allocatable :: Chi(:,:,:,:),ChiHF(:,:)
    character ( len = 80 )            :: filename,filename2,filename3,filenamemin,string
    character ( len = 80 )            :: filenamePNim(0:4),filenamePNre(0:4)
    real (kind=pr) ::start,finish
    ! 
    logical :: isospin_representation
    !

    rhoN=0.5_pr*rho*(1+y)
    rhoP=rho-rhoN

    call potential(rhoN,rhoP,UpotN,UpotP)
   

    call nmfp(rho,TMeV,Y,cross,energy)
       
          return
  end subroutine asymmetric


  !
  !! *************************************************************************
  !

  subroutine potential(rhoN,rhoP,UpotN,UpotP)
    !=========================================================================
    ! abstrac: routine to calculate the potential
    !
    ! input
    ! rhoN,rhoP= density of neutron and protons
    !
    ! output
    ! UpotN,UpotP [Mev ]
    !========================================================================
    implicit none
    real ( kind = pr ), INTENT(in) :: rhoN,rhoP
    integer            :: iso
    real ( kind = pr ) :: Upot(0:1),sig,eps,Ex,tau1
    real ( kind = pr ) :: rho_sig,rho_sig2,sig1_2,sig2_2
    real ( kind = pr ) :: rho0,rho1,tau(0:1),tau0,rho(0:1)
    real ( kind = pr ) :: UpotN,UpotP
    !
    rho(0)=rhoN
    rho(1)=rhoP
    rho0=rhoN+rhoP
    rho1=rhoN-rhoP
    rho_sig = rho0**sigma
    rho_sig2 = rho0**sigma2

    Ex=(rhoN-rhoP)/rho0
    tau0=t35*tpi23*rho0**t53*F0m(t53,Ex,0.0_pr,0.0_pr)
    tau1=t35*tpi23*rho0**t53*F1m(t53,Ex,0.0_pr,0.0_pr)


    tau(0)=(tau0+tau1)/2.0_pr
    tau(1)=(tau0-tau1)/2.0_pr


    !write(*,*)( c_tau(0) - c_tau(1) ) * tau0   
    !write(*,*) 2 * c_tau(1) *tau(0),2 * c_tau(1) *tau(1)


    eps = 10 * spacing(1.0_pr)

    sig1_2 = sigma + 2
    sig2_2 = sigma2 + 2

    !
    ! 0: neutron
    ! 1: protons
    !

    do iso=0,1
       sig=(-1.d0)**iso
       Upot(iso)= 2 * ( c_rho(0,1) - c_rho(1,1) ) * rho0              & ! 2 body
            & + 4 * c_rho(1,1) * rho(iso)                             &
            & + ( c_tau(0) - c_tau(1) ) * tau0                        &
            & + 2 * c_tau(1) * tau(iso)                               &
            & + sig1_2 * ( c_rho(0,2) - c_rho(1,2))*rho_sig*rho0      &
            & + 2 * sigma * c_rho(1,2)                                &
            &     * rho_sig / ( rho0 + eps ) * ( rhon**2 + rhop**2 )  &
            & + 4 * c_rho(1,2) * rho_sig * rho(iso)                   &
            & + sig2_2 * ( c_rho(0,3) -c_rho(1,3))*rho_sig2*rho0      &
            & + 2 * sigma2 * c_rho(1,3)                               &
            &    * rho_sig2 / ( rho0 + eps ) * ( rhon**2 + rhop**2 )  &
            & + 4 * c_rho(1,3) * rho_sig2 * rho(iso)                  & 
            & + 3 * Brho_0 * rho0**2                                  & ! rho rho rho  3 body
            & + Brho_1 * ( rho1**2 +sig* 2 * rho0 * rho1 )            & !
            & + 2 * Btau_0  * rho0 * tau0                             & ! rho rho tau
            & +     Btau_1  * ( rho1 * tau1 +sig*rho0 * tau1 )        &
            & +sig* 2 * Btau_10 * tau0 * rho1    


    enddo
    !
    !--------------------------------------------------------------------------
    !
    UpotN=Upot(0)
    UpotP=Upot(1)

    return
  end subroutine potential

  subroutine dnmfp(E,T,rhoN,rhoP,mu,MstarN,MstarP,q,omega,dcross)

        integer:: spin, proj,size,iw
        real ( kind = pr )     ::T,rhoN,rhoP,mstarN,mstarP,q,omega,E,mu(0:1),kfn
        complex ( kind = pr )             :: Amatrix(1:8,1:8),Bmatrix(1:8)
        complex ( kind = pr )             :: Amatrixp(1:8,1:8),Bmatrixp(1:8),betann(0:8)
        complex ( kind = pr )             :: Chinn,Chinp,Chipn,Chipp,beta0
        REAL ( kind = pr )             :: Chi(0:1,0:1,0:3)

        real ( kind = pr ) :: xRVO,xRAL1,xRAT1,q0,qMev,Enu,aux,ga2, ga=1.255_pr
        
        logical ::asymptotic=.false.
        real ( kind = pr )     dcross,wv,wav,wb,x
       
        ga2=ga*ga
                
        q0=omega
        qMev=q
        Enu=E
        !write(*,*)q0,qMev,Enu
        dcross=0.0
        factor=1.0_pr
        if (T.ne.0d0) then  
        factor=((1.0_pr-exp(-omega/T)))
        end if 
        if(q0.eq.0.0_pr)factor=1.0_pr
        !factor=1.0_pr
        
        if (.true.) then

          do spin=0,1
                do proj=0,spin,1
                call  matrix_nochex(-1,proj,T,rhoN,rhoP,mu,mstarN,mstarP,qMeV,omega,&
                     Amatrix,Bmatrix,Amatrixp,Bmatrixp,size,asymptotic)

                call cramer(Amatrix,Bmatrix,Amatrixp,Bmatrixp,size,Chinn,Chinp, &
                     Chipn,Chipp)
                
                   
                   Chi(spin,proj,0)=aimag(Chinn)/factor ! nn
                   Chi(spin,proj,1)=aimag(Chinp)/factor ! np
                   Chi(spin,proj,2)=aimag(Chipn)/factor ! pn
                   Chi(spin,proj,3)=aimag(Chipp)/factor ! pp
             
                   !write(*,*) spin, proj, aimag(chinn),q,omega
                   enddo
             enddo
                
                wv=chi(0,0,0)+0.08*(chi(0,0,1)+chi(0,0,2))+0.0064*chi(0,0,3)
                wav=chi(1,0,0)-(chi(1,0,1)+chi(1,0,2))*1.0+chi(1,0,3)
                wb=chi(1,1,0)-(chi(1,1,1)+chi(1,1,2))*1.0+chi(1,1,3)
        else

                !!This is the RPA approximation our virial fit was callibrated too


                    kfn=tpi13* rhoN**t13  ! [fm^-1]
                    
                    call betannall(omega,mstarN,qMev,T,mu(0),kfn,betann)
     
                    beta0=betann(0)
                    wav=-2*aimag(beta0/(1-4.5e-5*beta0))/factor/pi/hbarc3
                    wb=wav
                    wv=-2*aimag(beta0/(1-1.76e-5*beta0))/factor/pi/hbarc3
                
        end if
        
    aux=(q0-2*Enu)**2
    xRV0=qMeV*(-qMeV**2+aux)*wv
    xRAL1=ga2*q0**2*(-qMeV**2+aux)*wav/(tiny(1.d0)+qMeV)
    xRAT1=ga2*(qMeV**2-q0**2)*(qMeV**2+aux)*wav/(tiny(1.d0)+qMeV)
    dcross=1.d0/(32.0_pr*Enu**2*pi)*(xRV0+xRAL1+xRAT1)*2
  
 end subroutine dnmfp

  subroutine nmfp(rho,T,Y,avgnmfpout,energy)

        real ( kind = pr )     ::T,rho,rhoN,rhoP,mstarN,mstarP,q,omega,Y,mu(0:1),UpotN,UpotP,energy
        integer                 ::iw,iq,it,npw=npt,npq=npt,nptemp=npt,num_threads
        real(kind=pr)           :: q0i,q0f,qi,qf,xcut
        real ( kind = pr ),dimension(1000) :: wgauss,xgauss, wq0,xq0 ! weights and points of GL integration
        real (kind=pr)      ::dcross,cross,avgnmfp,sflux,avgnmfpout
        real (kind=pr),dimension(:),allocatable  :: llw,llq
        

        xcut=7
        gfactor = 1.0_pr/5.297376431_pr         !
        
    rhoP=rho*Y
    rhoN=rho-rhoP

    
    call mass_skyrme(rhoN,rhoP,mstarN,mstarP)
        
    call chemical(T,0.0_pr,0.0_pr,rhoN,rhoP,mstarN,mstarP,mu)
!$ num_threads=omp_get_max_threads()
!$ call omp_set_num_threads(num_threads)
    

        allocate(llq(npq),llw(npw))

        q0i=-xcut*energy
        q0f=energy
        CALL gset(q0i,q0f,npw,xq0,wq0)

        
        !$omp parallel do default(shared) private(iw,qf,qi,omega,xgauss,wgauss,dcross,llq)        
          do iw=1,npw
             omega=xq0(iw)!  MeV      
             qf=(2.0_pr*energy-omega) ! [MeV]
             qi=abs(omega)  ! [MeV]

             CALL gset(qi,qf,npq,xgauss,wgauss)

             do iq=1,npq
                q=xgauss(iq) ! MeV
                
                call dnmfp(energy,T,rhoN,rhoP,mu,mstarN,mstarP,q,omega,dcross)
                if(isnan(dcross)) dcross=0.0d0
                llq(iq)=dcross*wgauss(iq)
             end do

        llw(iw)=sum(llq)*wq0(iw)

        end do
        !$omp end parallel do
        avgnmfp=gfactor/sum(llw). ![m]
        cross=1.d0/avgnmfp
        


       avgnmfpout=avgnmfp     


     end subroutine nmfp





  subroutine  chemical(T,UpotN,UpotP,rhoN,rhoP,mstarN,mstarP,mu)
    !====================================================
    ! Subroutine to calculate the chemical potential
    !
    ! Input
    ! T: temperature
    ! UpotN(P): Neutron (proton) mean field
    ! rhoN(P): Neutron (proton) density
    ! 
    ! Output
    ! muN(p): Neutron (proton) chemical potential
    !====================================================
    implicit none

    integer                        :: iso,iflag
    real ( kind = pr ), intent(in) :: T,UpotN,UpotP,rhoN,rhoP,mstarN,mstarP
    real ( kind = pr )             :: U(0:1),mu(0:1),muLow(0:1),muUp(0:1)
    real ( kind = pr )             :: kf(0:1),m(0:1),rho(0:1)
    real ( kind = pr )             :: ef(0:1)
    real ( kind = pr )             :: mutest(0:1)

    U(0)=UpotN ; U(1)=UpotP
    m(0)=hbarc2/(2.0_pr*mstarN) ; m(1)=hbarc2/(2.0_pr*mstarP)

    rho(0)=rhoN ; rho(1)=rhoP
    kf(:)=tpi13*rho(:)**t13 
    ef(:)=kf(:)**2*m(:)

    do iso=0,1
       if(T/(ef(iso)).lt.0.02_pr)then
          !write(*,*)'Very low temperature.... ',iso
          mu(iso)=(ef(iso)+U(iso))*(1.d0-pi2/12.0_pr*(T/(ef(iso)+U(iso)))**2.0_pr)
       else
          !-- first ansatz
          muLow(iso)=-Max(40000000000000000.5_pr*T*ef(iso),ef(iso))
          muUp(iso)=Max(1.5_pr*T*ef(iso),ef(iso))
          if(rho(iso).eq.0.0_pr)then
             mutest(iso)=0.0_pr
             go to 100
          endif

          iflag=0
          call  bisezione(ef(iso),muLow(iso),muUp(iso),m(iso),U(iso),T, &
               rho(iso),mutest(iso),iflag)

100       continue

          if(iflag.eq.0)then
             mu(iso)=mutest(iso)
          endif

          if(iflag.eq.1)then
             !write(*,*)'second attemp... for chemical',iso
             iflag=0
             if(mutest(iso).eq.muLow(iso))then
                muUp(iso)=mutest(iso)
                muLow(iso)=50*mutest(iso)
             endif
             if(mutest(iso).eq.muUP(iso))then
                muUp(iso)=50*mutest(iso)
                muLow(iso)=mutest(iso)
             endif

             if(iflag.eq.1)then
                write(*,*)'I do not find chemical potential'
                write(*,*)muUp(iso),muLow(iso),muTest(iso)
                stop
             else
                mu(iso)=mutest(iso)
             endif

          endif
       endif
    enddo
    return
  end subroutine chemical
  !
  !!
  !
  subroutine bisezione(ef,mL,mU,massa,U,T,rho,muX,iflag)
    implicit none
    !====================================================
    ! abstract: bisection method for chemical potential
    ! input:
    ! ef: fermi energy
    ! mL,mU: starting values of chemical potential
    ! massa,U: mean field
    ! T: temperature
    ! rho: desired value of density
    !
    ! output
    ! muX: chemical potential found
    ! iflag: falg when i find (0) or not (1) the
    !        chemical potential
    !====================================================

    real ( kind = pr ),intent (in) :: ef,U,massa,T,rho
    integer                        :: i,ipt,iflag
    integer, parameter             :: iter=1000
    real ( kind = pr )             :: x(npt),w(npt),sum,mL,mU,muX
    real ( kind = pr )             :: integranda,esp,rhotest,k,a,b
    real ( kind = pr), parameter   :: Etol=0.0000001_pr
    real ( kind = pr ),parameter   :: xcut=5.0_pr ! cut off in the integral

    a=0.d0
    b=sqrt((xcut*T+Max(0.0_pr,ef-U))/massa)

    do i=1,iter
       muX=(mL+mU)/2.0_pr

       CALL gset(a,b,npt,x,w)
       sum=0.0_pr
       do ipt=1,npt
          k=x(ipt)
          esp=massa*k**2+U-muX
          integranda=k**2/(1.0_pr+exp(esp/T))
          sum=sum+w(ipt)*integranda/dpi2              
       enddo
       rhotest=2.0_pr*sum ! spin degeneracy
       if(abs(rhotest/rho-1).lt.Etol)then
          go to 100
       else
          if(rhotest.gt.rho)then
             mU=muX
          else
             mL=muX
          endif
       endif
       !write(*,*)rhotest,rho,muX
       if(i.eq.iter)then
          iflag=1
       endif
    enddo
100 continue
        
    


    return
  end subroutine bisezione
  !
  !!

end module asym
