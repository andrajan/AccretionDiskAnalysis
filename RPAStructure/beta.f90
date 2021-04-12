module beta

  use asym_cstes ! module constants
  use asym_gauss ! module to integrate


contains

  subroutine betannall(omega,msn,q,T,mun,kfn,betann)
    !-------------------------------------------
    ! Abstract: subroutine that calculates the
    ! beta_i^nn functions
    !
    ! input:
    ! omega: energy [MeV]
    ! msn: effective mass [MeV]
    ! q: transfer momentum [MeV]
    ! T: temperature [MeV]
    ! mun: chemical potential [MeV]
    !
    ! output
    ! betann(0-8): beta functions for fiexd
    !              omega and q
    !-------------------------------------------

    implicit none
    integer                       :: i,j
    real ( kind = pr ),intent(in) :: omega,msn,q,T,mun,kfn
    real ( kind = pr )            :: emois,eplus,kmois
    real ( kind = pr )            :: aux1,aux2,omegap
    real ( kind = pr )            :: betaIM(0:8),betaR(0:8)
    real ( kind = pr )            :: x(npt),w(npt),sum,e
    real ( kind = pr )            :: x2(npt),w2(npt),sum2,e2
    real ( kind = pr )            :: x3(npt),w3(npt),sum3,e3
    real ( kind = pr ), parameter :: xcut=6.0_pr
    real ( kind = pr ), parameter :: eps=0.10_pr
    real ( kind = pr )            :: xr1(npt),wr1(npt)
    real ( kind = pr )            :: xr2(npt),wr2(npt)
    real ( kind = pr )            :: xr3(npt),wr3(npt)
    real ( kind = pr )            :: msn2,q2,q4,elimP,elimM
    real ( kind = pr )            :: IM(0:8),efn
    complex ( kind = pr )         :: betann(0:8)

    !if(omega.gt.0_pr)then
    emois=(msn*omega/q - q/2.0_pr)**2/(2.0_pr*msn) ![MeV]
    eplus=(msn*omega/q + q/2.0_pr)**2/(2.0_pr*msn) ![MeV]
    !else
    !   emois=(msn*omega/q + q/2.0_pr)**2/(2.0_pr*msn) ![MeV]
    !   eplus=(msn*omega/q - q/2.0_pr)**2/(2.0_pr*msn) ![MeV]
    !endif

    kmois=msn*omega/q - q/2.0_pr                   ![MeV]
    aux1=exp(-(eplus-mun)/T)
    aux2=exp(-(emois-mun)/T)

    efn=((kfn*hbar_c)**2/(2.d0*msn))
    betann(:)=cz

    !write(*,*)emois,eplus,(eplus-mun)/T,(emois-mun)/T

    msn2=msn**2
    q2=q**2
    q4=q**4

    CALL gset(emois,eplus           ,npt,x ,w )
    CALL gset(eplus,eplus+xcut*T+mun,npt,x2,w2)



    !--------------------------------------------
    !           Immginary parts
    !--------------------------------------------

    ! beta0 nn

    betaIM(0)=T*(msn2)/(qpi*q)*Log((1.0_pr+aux1)/(1.0_pr+aux2)) ![MeV^2]


    ! beta1 nn
    betaIM(1)=(kmois/q)*betaIM(0)

    ! beta2 nn
    sum=0.0_pr
    do i=1,npt
       e=x(i)
       sum=sum+Log(1+exp(-(e-mun)/T))*w(i)
    enddo

    betaIM(2)=T*msn2/(qpi*q)*(2.0_pr*msn/q**2)*( &
         emois*Log((1.0_pr+aux1)/(1.0_pr+aux2))-sum ) ![MeV^2]

    ! beta3 nn
    betaIM(3)=((kmois/q)**2)*betaIM(0)

    ! beta4 nn
    betaIM(4)=(kmois/q)*betaIM(2)

    ! beta5 nn
    sum=0.0_pr
    sum2=0.0_pr
    do i=1,npt
       e=x(i)
       sum=sum+e*T*Log(1+exp(-(e-mun)/T))*w(i)
       !
       e2=x2(i)
       sum2=sum2+T*Log(1+exp(-(e2-mun)/T))*w2(i)
    enddo

    betaIM(5)=-4.0_pr*(msn2/q4)*(msn2/(qpi*q))*(             &
         -emois**2*T*Log(1+aux1)+emois**2*T*Log(1+aux2) &
         +2.0_pr*omega*sum2+2.0_pr*sum)

    ! beta6 nn
    betaIM(6)=((kmois/q)**3)*betaIM(0)

    ! beta7 nn
    betaIM(7)=((kmois/q)**4)*betaIM(0)    

    ! beta8 nn
    betaIM(8)=((kmois/q)**2)*betaIM(2)

    !--------------------------------------------
    !            Real parts
    !--------------------------------------------
    betaR(:)=0.0_pr
    !
    ! definition of the limits


    elimP= q**2/(2.d0*msn)+q*kfn*hbar_c/msn ! Maximum energy exchanged on Fermi energy
    elimM=-(q**2/(2.d0*msn)+q*kfn*hbar_c/msn) !

    elimP= (1.0_pr+xcut*T/efn)*elimP
    elimM= (1.0_pr+xcut*T/efn)*elimM


    CALL gset(elimM,    elimP,    npt,xr1,wr1)
    CALL gset(elimM,    omega-eps,npt,xr2,wr2)
    CALL gset(omega+eps,elimP,    npt,xr3,wr3)


    IM(:)=0.0_pr

    do i=1,npt
       omegap=xr1(i)
       emois=(msn*omegap/q - q/2.0_pr)**2/(2.0_pr*msn) ![MeV]
       eplus=(msn*omegap/q + q/2.0_pr)**2/(2.0_pr*msn) ![MeV]
       kmois=msn*omegap/q - q/2.0_pr                   ![MeV]

       aux1=exp(-(eplus-mun)/T)
       aux2=exp(-(emois-mun)/T)

       !---- immaginary parts 
       IM(0)=T*(msn2)/(qpi*q)*Log((1.0_pr+aux1)/(1.0_pr+aux2)) ![MeV]
       IM(1)=(kmois/q)*IM(0)

       CALL gset(emois,eplus,npt,x2,w2)
       CALL gset(eplus,eplus+xcut*T+mun,npt,x3,w3)


       sum=0.0_pr
       sum2=0.0_pr
       sum3=0.0_pr
       do j=1,npt
          e2=x2(j)
          sum=sum+Log(1.0_pr+exp(-(e2-mun)/T))*w2(j)
          sum2=sum2+e2*T*Log(1.0_pr+exp(-(e2-mun)/T))*w2(j)
          !
          e3=x3(j)
          sum3=sum3+T*Log(1.0_pr+exp(-(e3-mun)/T))*w3(j)
       enddo
       IM(2)=T*msn2/(qpi*q)*(2.0_pr*msn/q**2)*( &
            emois*Log((1.0_pr+aux1)/(1.0_pr+aux2))-sum ) ![MeV^2]
       IM(3)=(kmois**2/q2)*IM(0)

       IM(4)=(kmois/q)*IM(2)

       IM(5)=-4.0_pr*(msn2/q4)*(msn2/(qpi*q))*(             &
            -emois**2*T*Log(1.0_pr+aux1)+emois**2*T*Log(1.0_pr+aux2) &
            +2.0_pr*omegap*sum3+2.0_pr*sum2)

       IM(6)=((kmois/q)**3)*IM(0)
       IM(7)=((kmois/q)**4)*IM(0)    
       IM(8)=((kmois/q)**2)*IM(2)


       betaR(:)= betaR(:) +wr1(i)*(IM(:)-betaIM(:))/(omegap-omega) &
            +betaIM(:)/(xr2(i)-omega)*wr2(i)                       &
            +betaIM(:)/(xr3(i)-omega)*wr3(i)
    enddo

    betann(:)=Cmplx(betaR(:)/pi,betaIM(:), kind =pr) ! Mev^2

    return
  end subroutine betannall

  !
  !!
  !

  subroutine betapnall(omega,msn,msp,q,T,mun,mup,Un,Up,kfn,kfp,betapn)
    !-------------------------------------------
    ! Abstract: subroutine that calculates the
    ! beta_i^pn functions
    !
    ! input:
    ! omega: energy [MeV]
    ! msn: effective mass [MeV]
    ! q: transfer momentum [MeV]
    ! T: temperature [MeV]
    ! mun: chemical potential [MeV]
    !
    ! output
    ! betapn(0-8): beta functions for fiexd
    !              omega and q
    !-------------------------------------------

    implicit none
    integer                       :: i,j
    real ( kind = pr ),intent(in) :: omega,msn,msp,q,T,mun,mup,kfn,kfp,Un,Up
    complex ( kind = pr )         :: betapn(0:8)
    real ( kind = pr ), parameter :: limdm=0.1_pr
    real ( kind = pr ), parameter :: xcut=7.0_pr
    real ( kind = pr ), parameter :: eps=0.1_pr
    real ( kind = pr )            :: efn,efp,deltam,w,kmp,kmn
    real ( kind = pr )            :: x1n(npt),w1n(npt),x1p(npt),w1p(npt)
    real ( kind = pr )            :: betaIM(0:8),betaR(0:8),en,ep
    real ( kind = pr )            :: lim1,lim2,e1n,e2n,e1p,e2p
    real ( kind = pr )            :: q2,xf
    real ( kind = pr )            :: q3,q4,q6,q8,sumIMn(0:8),sumIMp(0:8),sumTIM(0:8)
    real ( kind = pr )            :: betaIMn(0:8),betaIMp(0:8),exon,exop
    real ( kind = pr )            :: elimP,elimM,xr1(npt),xr2(npt),xr3(npt)
    real ( kind = pr )            :: wr1(npt),wr2(npt),wr3(npt),omegap
    real ( kind = pr )            :: elim_neu(0:1),elim_pro(0:1)
    real ( kind = pr )            :: valmaxn,valmaxp

    efn=((kfn*hbar_c)**2/(2.d0*msn)) ! [Mev]
    efp=((kfp*hbar_c)**2/(2.d0*msp)) ! [Mev]
    deltam=msp-msn ! [Mev]
    w=omega-(Up-Un) ! [Mev]
    kmp=msp*w/q-q/2.0_pr  ! [Mev] 
    kmn=msn*w/q-q/2.0_pr  ! [Mev]

    betapn(:)=cz
    betaIM(:)=0.0_pr
    betaR(:)=0.0_pr

    q2=q*q   ! MeV^2
    q3=q2*q  ! MeV^3
    q4=q2*q2 ! MeV^4
    q6=q3*q3 ! MeV^6
    q8=q4*q4 ! MeV^8
    xf=msn*msp/(qpi*q) ! [MeV]


    !write(*,*)'Delta M=', deltam
    if(abs(deltam).gt.limdm)then!(remark n and p are inverted expressely here!!!)
       e1n=msn*q**2/(2.0_pr*deltam**2)*(1+dsqrt(1-(2.0_pr*kmp*deltam)/(q*msn)))**2 ! [Mev]
       e2n=msn*q**2/(2.0_pr*deltam**2)*(1-dsqrt(1-(2.0_pr*kmp*deltam)/(q*msn)))**2 ! [Mev]
       e1p=msp*q**2/(2.0_pr*deltam**2)*(1+dsqrt(1-(2.0_pr*(kmn+q)*deltam)/(q*msp)))**2 ! [Mev]
       e2p=msp*q**2/(2.0_pr*deltam**2)*(1-dsqrt(1-(2.0_pr*(kmn+q)*deltam)/(q*msp)))**2 ! [Mev]

       lim1=min(e1n,e2n)
       if(lim1.gt.0.0_pr)lim1=min(lim1,  xcut*T+Max(0.0_pr,mun-Un) )
       if(lim1.lt.0.0_pr)lim1=max(lim1,-(xcut*T+Max(0.0_pr,mun-Un)))

       lim2=max(e1n,e2n)
       if(lim2.gt.0.0_pr)lim2=min(lim2,  xcut*T+Max(0.0_pr,mun-Un) )
       if(lim2.lt.0.0_pr)lim2=max(lim2,-(xcut*T+Max(0.0_pr,mun-Un)))
       CALL gset(lim1,lim2,npt,x1n,w1n)
       !------
       lim1=min(e1p,e2p)
       if(lim1.gt.0.0_pr)lim1=min(lim1,  xcut*T+Max(0.0_pr,mup-Up) )
       if(lim1.lt.0.0_pr)lim1=max(lim1,-(xcut*T+Max(0.0_pr,mup-Up)))

       lim2=max(e1p,e2p)
       if(lim2.gt.0.0_pr)lim2=min(lim2,  xcut*T+Max(0.0_pr,mup-Up) )
       if(lim2.lt.0.0_pr)lim2=max(lim2,-(xcut*T+Max(0.0_pr,mup-Up)))
       CALL gset(lim1,lim2,npt,x1p,w1p)

       if(is_nan(e1n))then
          !write(*,*)'The integral of neutrons does not exist',omega
          e1n=0.0_pr
          e2n=0.0_pr    
       endif

       if(is_nan(e1p))then
          !write(*,*)'The integral of protons does not exist',omega
          e1p=0.0_pr
          e2p=0.0_pr    
       endif
    else 
       !write(*,*)'passo qui?'
       e1n=xcut*T+Max(0.0_pr,(mun-Un))!huge(1.0_pr) ! [Mev]  cut off
       e2n=(msn*omega/q-q/2.0_pr)**2/(2.0_pr*msn) ! [Mev]
       e1p=xcut*T+Max(0.0_pr,(mup-Up))!huge(1.0_pr) ! [Mev] cut off
       e2p=(msp*omega/q+q/2.0_pr)**2/(2.0_pr*msp) ! [Mev]

       lim1=e2n
       lim2=e1n
       if(lim1>lim2)then
          lim1=0.0_pr ; lim2=0.0_pr
       endif
       CALL gset(lim1,lim2,npt,x1n,w1n)
       lim1=e2p
       lim2=e1p
       if(lim1>lim2)then
          lim1=0.0_pr ; lim2=0.0_pr
       endif
       CALL gset(lim1,lim2,npt,x1p,w1p)


    endif

    !--------------------------------------------
    !           Immginary parts
    !--------------------------------------------


    betaIMn(:)=0.0_pr ; betaIMp(:)=0.0_pr
    if((1-(2.0_pr*kmp*deltam)/(q*msn))<0.0_pr) then
       betaIMn(:)=0.0_pr
    else
       if(e1n.ne.e2n)then
          do i=1,npt
             en=x1n(i)
             exon=(1.0_pr+exp((en+Un-mun)/T))

             betaIMn(0)=betaIMn(0)+w1n(i)/exon
             betaIMn(1)=betaIMn(1)+w1n(i)*(kmp/q+deltam*en/q2)/exon
             betaIMn(2)=betaIMn(2)+w1n(i)*(2.0_pr*msn/q2)*en/exon
             betaIMn(3)=betaIMn(3)+w1n(i)*(kmp*q+en*deltam)**2/q4/exon
             betaIMn(4)=betaIMn(4)+w1n(i)*(2.0_pr*msn/q2)*(en*kmp/q+deltam*en**2/q2)/exon
             betaIMn(5)=betaIMn(5)+w1n(i)*(4.0_pr*msn**2/q4)*(en**2)/exon
             betaIMn(6)=betaIMn(6)+w1n(i)*(deltam*en+q*kmp)**3/q6/exon
             betaIMn(7)=betaIMn(7)+w1n(i)*(deltam*en+q*kmp)**4/q8/exon
             betaIMn(8)=betaIMn(8)+w1n(i)*2.0_pr*msn*en*(deltam*en+q*kmp)**2/q6/exon

          enddo
       endif
    endif

    if((1-(2.0_pr*(kmn+q)*deltam)/(q*msp))<0.0_pr) then
       betaIMp(:)=0.0_pr
    else
       if(e1p.ne.e2p)then

          do i=1,npt
             ep=x1p(i)
             exop=(1.0_pr+exp((ep+Up-mup)/T))

             betaIMp(0)=betaIMp(0)+w1p(i)/exop
             betaIMp(1)=betaIMp(1)+w1p(i)*(kmn/q+deltam*ep/q2)/exop
             betaIMp(2)=betaIMp(2)+w1p(i)*(2.0_pr*msn/q2)*(ep-w)/exop
             betaIMp(3)=betaIMp(3)+w1p(i)*(kmn*q+ep*deltam)**2/q4/exop
             betaIMp(4)=betaIMp(4)+w1p(i)*2.0_pr*msn/q2 &
                  *(deltam/q2*ep*(ep-w)+kmn/q*(ep-w))/exop
             betaIMp(5)=betaIMp(5)+w1p(i)*(4.0_pr*msn**2/q4)*((ep-w)**2)/exop
             betaIMp(6)=betaIMp(6)+w1p(i)*(deltam*ep+q*kmn)**3/q6/exop
             betaIMp(7)=betaIMp(7)+w1p(i)*(deltam*ep+q*kmn)**4/q8/exop
             betaIMp(8)=betaIMp(8)+w1p(i)*2.0_pr*msn*(ep-w)*(deltam*ep+q*kmn)**2/q6/exop
          enddo
       endif
    endif


    ! beta np
    betaIM(:)=xf*(-betaIMn(:)+betaIMp(:)) ! MeV^2


    !--------------------------------------------
    !          Real parts 
    !--------------------------------------------

    ! --- neutrons
    elim_neu(0)=q2/(2.d0*msn)+q*kfn*hbar_c/msn
    elim_neu(0)=(1.0d0+xcut*T/efn)*elim_neu(0)

    elim_neu(1)=-(q2/(2.d0*msn)+q*kfn*hbar_c/msn)
    elim_neu(1)=(1.0d0+xcut*T/efn)*elim_neu(1)

    ! --- protons
    elim_pro(0)=q2/(2.d0*msp)+q*kfp*hbar_c/msp
    elim_pro(0)=(1.0d0+xcut*T/efp)*elim_pro(0)

    elim_pro(1)=-(q2/(2.d0*msp)+q*kfp*hbar_c/msp)
    elim_pro(1)=(1.0d0+xcut*T/efp)*elim_pro(1)


    elimP=max(elim_neu(0),elim_pro(0))
    elimM=min(elim_neu(1),elim_pro(1)) !min


    CALL gset(elimM,    elimP,    npt,xr1,wr1)
    CALL gset(elimM,    omega-eps,npt,xr2,wr2)
    CALL gset(omega+eps,elimP,    npt,xr3,wr3)



    do j=1,npt
       omegap=xr1(j)
       w=omegap-(Up-Un) ! [Mev]
       kmp=msp*w/q-q/2.0_pr  ! [Mev] 
       kmn=msn*w/q-q/2.0_pr  ! [Mev]
       valmaxn=xcut*T+Max(0.0_pr,mun-Un)

       if(abs(deltam).gt.limdm)then!(remark km n and p are inverted expressely here!!!)
          e1n=msn*q**2/(2.0_pr*deltam**2)*(1.0_pr+dsqrt(1.0_pr-(2.0_pr*kmp*deltam)/(q*msn)))**2 ! [Mev]
          e2n=msn*q**2/(2.0_pr*deltam**2)*(1.0_pr-dsqrt(1.0_pr-(2.0_pr*kmp*deltam)/(q*msn)))**2 ! [Mev]

          if(is_nan(e1n))e1n=0.0_pr
          if(is_nan(e2n))e2n=0.0_pr
          lim1=min(e1n,e2n)
          lim2=max(e1n,e2n)
          lim2=min(lim2,valmaxn )
          CALL gset(lim1,lim2,npt,x1n,w1n)

       else 
          e1n=xcut*T+Max(0.0_pr,mun-Un) ! [Mev]  cut off
          e2n=(msn*omegap/q-q/2.0_pr)**2/(2.0_pr*msn) ! [Mev]
          if(is_nan(e1n))e1n=0.0_pr
          if(is_nan(e2n))e2n=0.0_pr

          lim1=e2n
          lim2=e1n
          if(lim1>lim2)then
             lim1=0.0_pr ; lim2=0.0_pr
          endif
          CALL gset(lim1,lim2,npt,x1n,w1n)
          !
       endif

       sumIMn(:)=0.0_pr 
       ! neutrons 
       if((1-(2.0_pr*kmp*deltam)/(q*msn))<0.0_pr.and.omegap.le.elimP.and.omegap.ge.elimM) then
          sumIMn(:)=0.0_pr
       else
          if(valmaxn.lt.lim1)then
             sumIMn(:)=0.0_pr
          else
             if(lim1.lt.lim2)then
                do i=1,npt
                   en=x1n(i)
                   exon=(1.0_pr+exp((en+Un-mun)/T))

                   sumIMn(0)=sumIMn(0)+w1n(i)/exon
                   sumIMn(1)=sumIMn(1)+w1n(i)*(kmp/q+deltam*en/q2)/exon
                   sumIMn(2)=sumIMn(2)+w1n(i)*(2.0_pr*msn/q2)*en/exon
                   sumIMn(3)=sumIMn(3)+w1n(i)*(kmp*q+en*deltam)**2/q4/exon
                   sumIMn(4)=sumIMn(4)+w1n(i)*(2.0_pr*msn/q2)*(en*kmp/q+deltam*en**2/q2)/exon
                   sumIMn(5)=sumIMn(5)+w1n(i)*(4.0_pr*msn**2/q4)*(en**2)/exon
                   sumIMn(6)=sumIMn(6)+w1n(i)*(deltam*en+q*kmp)**3/q6/exon
                   sumIMn(7)=sumIMn(7)+w1n(i)*(deltam*en+q*kmp)**4/q8/exon
                   sumIMn(8)=sumIMn(8)+w1n(i)*2.0_pr*msn*en*(deltam*en+q*kmp)**2/q6/exon

                enddo
             endif
          endif
       endif

       !=== protons
       valmaxp=xcut*T+Max(0.0_pr,mup-Up)

       if(abs(deltam).gt.limdm)then!(remark km n and p are inverted expressely here!!!)            
          e1p=msp*q**2/(2.0_pr*deltam**2)*(1.0_pr+dsqrt(1.0_pr-(2.0_pr*(kmn+q)*deltam)/(q*msp)))**2 ! [Mev]
          e2p=msp*q**2/(2.0_pr*deltam**2)*(1.0_pr-dsqrt(1.0_pr-(2.0_pr*(kmn+q)*deltam)/(q*msp)))**2 ! [Mev]
          if(is_nan(e1p))e1p=0.0_pr
          if(is_nan(e2p))e2p=0.0_pr
          !------
          lim1=min(e1p,e2p)
          lim2=max(e1p,e2p)
          lim2=min(lim2, valmaxp )
          CALL gset(lim1,lim2,npt,x1p,w1p)

       else 
          e1p=xcut*T+Max(0.0_pr,mup-Up) ! [Mev] cut off
          e2p=(msp*omegap/q+q/2.0_pr)**2/(2.0_pr*msp) ! [Mev]
          if(is_nan(e1p))e1p=0.0_pr
          if(is_nan(e2p))e2p=0.0_pr
          !
          lim1=e2p
          lim2=e1p
          if(lim1>lim2)then
             lim1=0.0_pr ; lim2=0.0_pr
          endif
          CALL gset(lim1,lim2,npt,x1p,w1p)
       endif

       sumIMp(:)=0.0_pr
       ! protons
       if((1-(2.0_pr*(kmn+q)*deltam)/(q*msp))<0.0_pr.and.omegap.le.elimP.and.omegap.ge.elimM) then
          sumIMp(:)=0.0_pr
       else
          if(valmaxp.lt.lim1)then
             sumIMp(:)=0.0_pr
          else
             if(lim1.lt.lim2)then
                do i=1,npt
                   ep=x1p(i)
                   exop=(1.0_pr+exp((ep+Up-mup)/T))

                   sumIMp(0)=sumIMp(0)+w1p(i)/exop
                   sumIMp(1)=sumIMp(1)+w1p(i)*(kmn/q+deltam*ep/q2)/exop
                   sumIMp(2)=sumIMp(2)+w1p(i)*(2.0_pr*msn/q2)*(ep-w)/exop
                   sumIMp(3)=sumIMp(3)+w1p(i)*(kmn*q+ep*deltam)**2/q4/exop
                   sumIMp(4)=sumIMp(4)+w1p(i)*2.0_pr*msn/q2 &
                        *(deltam/q2*ep*(ep-w)+kmn/q*(ep-w))/exop
                   sumIMp(5)=sumIMp(5)+w1p(i)*(4.0_pr*msn**2/q4)*((ep-w)**2)/exop
                   sumIMp(6)=sumIMp(6)+w1p(i)*(deltam*ep+q*kmn)**3/q6/exop
                   sumIMp(7)=sumIMp(7)+w1p(i)*(deltam*ep+q*kmn)**4/q8/exop
                   sumIMp(8)=sumIMp(8)+w1p(i)*2.0_pr*msn*(ep-w)*(deltam*ep+q*kmn)**2/q6/exop

                enddo
             endif
          endif
       endif

       sumTIM(:)=xf*(-sumIMn(:)+sumIMp(:)) ! MeV

       betaR(:)= betaR(:) +(sumTIM(:)-betaIM(:))/(xr1(j)-omega)*wr1(j) &
            +betaIM(:)/(xr2(j)-omega)*wr2(j)+betaIM(:)/(xr3(j)-omega)*wr3(j)

    enddo

    !-------------------------------------
    !        Result
    !-------------------------------------
    betapn(:)=Cmplx(betaR(:)/pi,betaIM(:), kind =pr) ! MeV^2

    return
  end subroutine betapnall

  !
  !!
  !
  !======= start package Lindhart functions ================================!

  FUNCTION ftheta(x)
    !=======================================================================!
    !  Heavyside theta function
    !  Entry : x (real)
    !  Sortie : ftheta qui vaut 0 ou 1 (reel)
    !=======================================================================!
    IMPLICIT none
    !-----------------------------------------------------------------------!
    real ( kind = pr ), INTENT(in) :: x
    real ( kind = pr ) :: ftheta
    !=======================================================================!
    IF ( x >= 0.0_pr ) THEN
       ftheta = 1.0_pr
    ELSE
       ftheta = 0.0_pr
    END IF
    !
  END FUNCTION ftheta
  !
  !!
  !
  FUNCTION aalind( nu, k ) result(res)
    !=======================================================================
    ! abstract:Fomula C.6 of Garcia-Recio Ann.Phys.214, 293-340 (1992)
    ! Only A+(nu,k) is defined since A-(nu,k) = A+(-nu,k)
    ! input: 
    ! k: momentum
    ! nu=frequency
    ! sgn= +1/-1 
    !=======================================================================!
    IMPLICIT none
    !-----------------------------------------------------------------------!
    REAL( kind = pr), INTENT(in) :: nu, k
    !-----------------------------------------------------------------------!
    COMPLEX( kind = pr ) :: res
    !=======================================================================!
    REAL( kind = pr ) :: fac, xlog, signu, kpn
    !
    if( nu == 0.0_pr ) then
       signu = 0.0_pr
    else
       signu = SIGN( +1.0_pr, nu)
    end if
    !
    kpn = k + nu
    fac = 1 - kpn**2
    xlog = ABS( ( kpn + 1 ) / ( kpn - 1 ) )
    if ( fac > 0.0_pr ) then
       res = CMPLX( log(xlog), - pi * signu, kind = pr )
    else
       res = CMPLX( log(xlog), 0.0_pr, kind = pr )
    end if
    res = res * fac / ( 4 * k )
    !
  END FUNCTION aalind
  !
  !!
  !
  FUNCTION lindpi0(nu, k)
    !=======================================================================
    ! abstract: Modified Lindhardt function Pi0 Eq. C3
    !  of Garcia-Recio Ann.Phys.214, 293-340 (1992)
    !  NOTE: the factor m^*k_F is NOT included in the routine, this is the 
    !  reason of the name "modified", it should be multiplied outsied this 
    !  function
    !
    ! input: 
    ! k: momentum
    ! nu= frequency
    !=======================================================================!
    IMPLICIT none
    !-----------------------------------------------------------------------!
    REAL(pr), INTENT(in) :: nu, k
    !-----------------------------------------------------------------------!
    COMPLEX(pr) :: lindpi0
    !=======================================================================!
    !

    lindpi0 = - ( 1.0_pr + aalind( nu, k ) + aalind( - nu, k) ) / qpi2
    !
  END FUNCTION lindpi0
  !
  !!
  !
  FUNCTION lindpi2(nu, k)
    !=======================================================================
    ! Eq. C.4 see comments on lindpi0
    !=======================================================================!
    IMPLICIT NONE
    !-----------------------------------------------------------------------!
    REAL(pr), INTENT(in) :: nu, k
    !-----------------------------------------------------------------------!
    COMPLEX(pr) :: lindpi2
    !=======================================================================!
    !

    lindpi2 = - ( 3.0_pr + k**2 - nu**2                  &
         + ( 1.0_pr + (k-nu)**2 ) * aalind(   nu, k )  &
         + ( 1.0_pr + (k+nu)**2 ) * aalind( - nu, k ) ) / 2 / qpi2
    !
    return
  END FUNCTION lindpi2
  !
  !!
  !
  FUNCTION lindpi4(nu, k)
    !=======================================================================
    ! Eq. C.5 see comments on lindpi0
    !=======================================================================!
    IMPLICIT none
    !-----------------------------------------------------------------------!
    REAL(pr), INTENT(in) :: nu, k
    !-----------------------------------------------------------------------!
    COMPLEX(pr) :: lindpi4
    real (pr) :: k2, k3, k4, nu2, nu4, kn2, dkn3
    real (pr), parameter :: t49_3 = 49.0_pr / 3.0_pr
    !=======================================================================!
    !
    k2 = k * k
    k3 = k2 * k
    k4 = k2 * k2
    nu2 = nu * nu
    nu4 = nu2 * nu2
    kn2 = k2 * nu2
    dkn3 = 2 * k * nu**3

    lindpi4 = - ( 5 + t49_3 * k2 + k4 - nu2 + 16* kn2 - nu4     &
         &  + ( 1 + k2 + k4 - 4*k*nu - 2*k3*nu + nu2 + 18*kn2 - dkn3 + nu4 ) &
         &    * aalind( nu, k )                          &
         &  + ( 1 + k2 + k4 + 4*k*nu + 2*k3*nu + nu2 + 18*kn2 + dkn3 + nu4 ) &
         &    * aalind( - nu, k ) ) / tqpi2
    !
  END FUNCTION lindpi4
  !================ end of package Lindhart ==================================

  !
  !!
  !


  subroutine betannallasymptotic(msn,qMeV,kfn,betann)
    !-------------------------------------------
    ! Abstract: subroutine that calculates the
    ! beta_i^nn functions in the limit w->0
    !
    ! input:
    ! msn: effective mass [MeV]
    ! q: transfer momentum [MeV]
    !
    ! output
    ! betann(0-8): beta functions for fiexd
    !              omega and q
    !-------------------------------------------

    implicit none
    real ( kind = pr ),intent(in) :: msn,qMeV,kfn
    real ( kind = pr )            :: betaR(0:8),kn,q,arg,ffk
    complex ( kind = pr )         :: betann(0:8)

    q=qMeV/hbar_c+1.0d-50

    betaR(:)=0.0_pr
    kn= q/(2.0_pr*kfn) ! no units

    arg=(kn+1)/(kn-1)  ! compared to Mathematica we need to introduce
    !                  ! an Abs value of
    if(arg.gt.0)then ! arg of Log, moreover we keep trace of
       !             ! the sign of the arg
       ffk=0.5_pr*(1+(1-kn**2)/(2*kn)*log(arg)) ! justified emipirically
    else
       ffk=0.5_pr*(1.0_pr+(1.0_pr-kn**2)/(2.0_pr*kn)*log(-arg))
    endif
    if(kn.eq.1.0_pr)ffk=0.5_pr
    if(is_nan(ffk)) then
            print *, kn,arg,log(-arg)
    end if
    
    if(is_nan(ffk))stop ' NaN in the logarithm!!!'

    ! beta 0

    betaR(0)=-ffk*Kfn*msn/(2*pi2)*hbar_c ! [MeV^2]

    ! beta 1

    betaR(1)=ffk*Kfn*msn/(4*pi2)*hbar_c ! [MeV^2]

    ! beta 2
    betaR(2)=-Kfn*(1.0_pr+ffk+ffk*kn**2)*msn*hbar_c/(16.0_pr*kn**2*pi2)

    ! beta 3
    betaR(3)=-Kfn*(1.0_pr+3.0_pr*ffk*kn**2)*msn*hbar_c/(24.0_pr*kn**2*pi2)

    ! beta 4
    betaR(4)=Kfn*(7.0_pr+3*ffk*(1.0_pr+kn**2))*msn*hbar_c/(96.0_pr*kn**2*pi2)

    ! beta 5
    betaR(5)=-Kfn*(6.0_pr+23.0_pr*kn**2+3*ffk*(1.0_pr+kn**2+kn**4))*msn*hbar_c/(288*kn**4*pi2)

    ! beta 6
    betaR(6)=Kfn*msn*hbar_c*(1.0_pr+ffk*kn**2)/(16.0_pr*kn**2*pi2)

    ! beta 7
    betaR(7)=-Kfn*msn*hbar_c*(3.0_pr+35*kn**2+15*ffk*kn**4)/(480.0_pr*kn**4*pi2)

    ! beta 8
    betaR(8)=-Kfn*msn*hbar_c*(2.0_pr+3.0_pr*(5.0_pr+ffk)*kn**2+3.0_pr*ffk*kn**4)/(192.0_pr*kn**4*pi2)

    betann(:)=Cmplx(betaR(:),0.0_pr, kind = pr)
    
    return
  end subroutine betannallasymptotic

  !
  !!
  !
  subroutine betapnTzero(omega,msn,msp,qMeV,Un,Up,kfn,kfp,betapn)
    !-------------------------------------------
    ! Abstract: subroutine that calculates the
    ! beta_i^pn functions at T=0
    !
    ! input:
    ! omega: energy [MeV]
    ! msn: effective mass [MeV]
    ! q: transfer momentum [MeV]
    !
    ! output
    ! betapn(0-8): beta functions for fiexd
    !              omega and q [MeV^2]
    !-------------------------------------------

    implicit none
    integer                       :: i
    real ( kind = pr ),intent(in) :: omega,msn,msp,qMeV,kfn,kfp,Un,Up
    complex ( kind = pr )         :: betapn(0:8)
    real ( kind = pr )            :: efn,efp,Bp,Bn,w,Abig,q,eps
    real ( kind = pr )            :: Cn,Cp,Deltan,Deltap,rhon,rhop,yy
    real ( kind = pr )            :: betaR(0:8),betaIM(0:8),Delta
    real ( kind = pr )            :: xfn,xfp,limGup,limGdown,limGtup,limGtdown
    real ( kind = pr )            :: alpha1(0:8),alpha2(0:8),kn,kp
    real ( kind = pr )            :: beta1(0:8),beta2(0:8)
    real ( kind = pr )            :: gamma1(0:8),gamma2(0:8)
    real ( kind = pr )            :: delta1(0:8),delta2(0:8)
    real ( kind = pr )            :: arc1,arc2,arc3,arc4
    real ( kind = pr )            :: log1,log2,log3,log4

    efn=((kfn*hbar_c)**2/(2.d0*msn)) ! [Mev]
    efp=((kfp*hbar_c)**2/(2.d0*msp)) ! [Mev]

    rhon=(kfn/tpi13)**3
    rhop=(kfp/tpi13)**3
    yy=(rhon-rhop)/(rhon+rhop)


    q=qMeV/hbar_c! [fm]

    kn=q/(2*kfn)
    kp=q/(2*kfp)

    betapn(:)=cz
    betaIM(:)=0.0_pr
    betaR(:)=0.0_pr

    eps=10.0_pr**(-10.0_pr)

    !---- preliminary parameters

    Bp=1.0_pr/(2*msp) ; Bn=1.0_pr/(2*msn) ! [MeV^-1]
    w=omega-Up+Un ! [MeV]
    Abig=Bn-Bp ! [MeV^-1]
    Cp=w-Bp*qMeV**2 ; Cn=w+Bn*qMeV**2 ![MeV]
    Deltap=Bp**2*qMeV**2-Abig*Cp; Deltan=Bn**2*qMeV**2-Abig*Cn
    Delta=2*Bn*Bp*(0.5_pr*qMeV**2-w*(msp-msn)) ! pure number

    !--- Limit Y->0 we need some manipulation


    !---- Imaginary parts of betapn
    if(Delta.ge.0.0_pr)then
       xfn=1.0_pr/(8*pi*qMeV*Bn) ; xfp=1.0_pr/(8*pi*qMeV*Bp)
       if(abs(yy).le.0.05_pr)then! limit Y=0
          limGup =kfn*hbar_c
          limGdown =kfn*hbar_c*abs(omega*msn/(qMeV*kfn*hbar_c)-kn)
          limGtup=kfn*hbar_c
          limGtdown=kfp*hbar_c*abs(omega*msp/(qMeV*kfp*hbar_c)+kp)
          Abig=eps ! Not 0 for numerical reasons
 
          !stop 'The code has not been tested yet!'
       else! .... normal
          limGup=Min(kfn*hbar_c, (qMeV*Bp+sqrt(Delta))/dabs(Abig)) ! [MeV]
          limGdown=abs((qMeV*Bp-sqrt(Delta))/Abig)     ! [MeV]

          limGtup=Min(kfp*hbar_c, (qMeV*Bn+sqrt(Delta))/dabs(Abig)) ! [MeV]
          limGtdown=abs((qMeV*Bn-sqrt(Delta))/Abig)     ! [MeV]  
       endif
       do i=0,8 ! MeV^2
          BetaIM(i)=-xfp*(Gbig(i,limGup)-Gbig(i,limGdown))*ftheta(limGup-limGdown) &
               +xfn*(Gtildebig(i,limGtup)-Gtildebig(i,limGtdown))*ftheta(limGtup-limGtdown) 
       enddo
    endif

    !--- Real parts of betanp in MeV^2
    arc1=0.0_pr; arc2=0.0_pr; arc3=0.0_pr; arc4=0.0_pr

    if(Delta.le.0.0_pr)then
       arc1=atan((Abig*kfn*hbar_c-Bp*qMeV)/sqrt(-Delta))
       arc2=atan((Abig*kfn*hbar_c+Bp*qMeV)/sqrt(-Delta))
       arc3=atan((Abig*kfp*hbar_c-Bn*qMeV)/sqrt(-Delta))
       arc4=atan((Abig*kfp*hbar_c+Bn*qMeV)/sqrt(-Delta))
    endif

    log1=0.0_pr; log2=0.0_pr
    if(Delta.gt.0.0_pr)then
       log1=abs((sqrt(abs(Delta))-Abig*kfn*hbar_c)**2-(qMeV*Bp)**2)/abs((sqrt(abs(Delta))+Abig*kfn*hbar_c)**2-(qMeV*Bp)**2)
       log1=log(abs(log1))
       !
       log2=((sqrt(abs(Delta))-Abig*kfp*hbar_c)**2-(qMeV*Bn)**2)/((sqrt(abs(Delta))+Abig*kfp*hbar_c)**2-(qMeV*Bn)**2)
       log2=log(abs(log2))

    endif

    log3=dabs(Abig*kfn**2*hbarc2-2*Bp*kfn*hbar_c*qMeV+Cp)/dabs(Abig*kfn**2*hbarc2+2*Bp*kfn*hbar_c*qMeV+Cp)
    log3=log((log3))

    log4=dabs(Abig*kfp**2*hbarc2-2*Bn*kfp*hbar_c*qMeV+Cn)/dabs(Abig*kfp**2*hbarc2+2*Bn*kfp*hbar_c*qMeV+Cn)
    log4=log((log4))


    !==== building the difference pieces of the formula...

    call ReA1(alpha1)
    call ReA2(alpha2)
    call ReB1(beta1) 
    call ReB2(beta2) 
    call ReD1(delta1)
    call ReD2(delta2)
    call ReG1(gamma1)
    call ReG2(gamma2)


    do i=0,8
       BetaR(i)=alpha1(i)+alpha2(i)  &
            +beta1(i)*(arc1+arc2)*ftheta(-Delta)& 
            +beta2(i)*(arc3+arc4)*ftheta(-Delta) &
            +delta1(i)*ftheta(Delta)*log1 &
            +delta2(i)*log2*ftheta(Delta) &
            +gamma1(i)*log3 &
            +gamma2(i)*log4
       !if(i.eq.5) write(*,*)'aa',BetaIM(i),BetaR(i)
      
       betapn(i)=Cmplx( BetaR(i),BetaIM(i), kind =pr)

       if(i.eq.8)then
         write(88,666)omega,alpha1(i)+alpha2(i),  &
            +beta1(i)*(arc1+arc2)*ftheta(-Delta)& 
            +beta2(i)*(arc3+arc4)*ftheta(-Delta), &
            +delta1(i)*ftheta(Delta)*log1 &
            +delta2(i)*log2*ftheta(Delta), &
            +gamma1(i)*log3 &
            +gamma2(i)*log4
666 format(5(1x,E15.6))
       endif
    enddo

    !-- putting togheter




  contains

    function Gbig(i,k) result(res) ! in [MeV^2]
      implicit none
      integer :: i
      real ( kind = pr ) :: res,k! [MeV]
      real ( kind = pr ) :: aux    
      select case(i)

      case (0)
         res=k**2/2
      case(1)
         res=(Abig*k**4+2*k**2*Cp)/(8*qMeV**2*Bp)
      case(2)
         res=k**4/(4*qMeV**2)
      case(3)
         res=(Abig**2*k**6+3*Abig*k**4*Cp+3*k**2*Cp**2)/(24*qMeV**4*Bp**2)
      case(4)
         res=(2*Abig*k**6+3*k**4*Cp)/(24*qMeV**4*Bp)
      case(5)
         res=k**6/(6*qMeV**4)
      case(6)
         res=((Abig*k**2+Cp)**4-Cp**4)/(64*Abig*qMeV**6*Bp**3)
      case(7)
         res=((Abig*k**2+Cp)**5-Cp**5)/(160*Abig*qMeV**8*Bp**4)
      case(8)
         aux=3*Abig**2*k**8+8*Abig*k**6*Cp+6*k**4*Cp**2
         res=aux/(96*qMeV**6*Bp**2)

      case default
         stop 'Error in function Gbig!'
      end select

      return
    end function Gbig
    !
    !!
    !

    function Gtildebig(i,k) result(res) ! in [MeV^2]
      implicit none
      integer :: i
      real ( kind = pr ) :: res,k! [MeV]
      real ( kind = pr ) :: k2,k4,k6,k8
      real ( kind = pr ) :: q2,q4,q6,q8,aux,aux2

      q2=qMeV**2; q4=q2*q2; q6=q4*q2; q8=q6*q2
      k2=k*k; k4=k2*k2; k6=k4*k2; k8=k6*k2
      res=0.0_pr
      select case(i)
      case(0)
         res=k2/2
      case(1)
         res=k2*(Abig*k2-4*q2*Bn+2*Cn)/(8*q2*Bn)
      case(2)
         res=k2*(-Abig*k2+Bn*(k2+2*q2)-2*Cn)/(4*q2*Bn)
      case(3)
         res=(Abig**2*k6+3*Abig*k4*(Cn-2*q2*Bn)+3*k2*(Cn-2*q2*Bn)**2)/(24*q4*Bn**2)
      case(4)
         aux=2*Abig**2*k4-Bn*(Abig*k2*(2*k2+9*q2)+3*Cn*(k2+6*q2))+6*Cn*(Abig*k2+Cn)+6*q2*Bn**2*(k2+2*q2)
         res=-k2*aux/(24*q4*Bn**2)
      case(5)
         aux=3*k2*(Bn-Abig)*(q2*Bn-Cn)+k4*(Abig-Bn)**2+3*(Cn-q2*Bn)**2
         res=k2*aux/(6*q4*Bn**2)
      case(6)
         aux=(Abig*k2-2*q2*Bn+Cn)**4-(Cn-2*q2*Bn)**4
         res=aux/(64*Abig*q6*Bn**3)
      case(7)
         aux=(Abig**4*k8-5*(2*q2*Bn-Cn)*(Abig*k2-2*q2*Bn+Cn)*(Abig**2*k4-(2*q2*Bn-Cn)*(Abig*k2-2*q2*Bn+Cn)))
         res=k2*aux/(160*q8*Bn**4)
      case(8)
         aux=3*Abig**2*k8*(Abig-Bn)+4*Abig*k6*(Cn*(3*Abig-2*Bn)+q2*Bn*(4*Bn-5*Abig))
         aux2=6*k4*(2*q2*Bn-Cn)*(Cn*(Bn-3*Abig)-2*q2*Bn*(Bn-2*Abig)) -12.0_pr*k2*(q2*Bn-Cn)*((Cn-2.0_pr*q2*Bn)**2)
         res=-(aux+aux2)/(96*q6*Bn**3) ! attenzione al segno meno in aux2
      case default
         stop 'Error in function Gtildebig!'
      end select
      return
    end function Gtildebig
    !
    !!
    !

    !
    !! === the different components of Real(betapn)
    !
    subroutine ReA1(alpha1)  ! alpha1(i)  [MeV^2]
      implicit none

      real ( kind = pr ) :: alpha1(0:8)
      real ( kind = pr ) :: q2,q4,q6,q8
      real ( kind = pr ) :: kfn2,kfn4,Bp2,Bp3,Bp4,Bp5,Bp6,Bp7,Bp8
      real ( kind = pr ) :: Abig2,Abig3,Abig4,Abig5

      alpha1(:)=0.0_pr

      q2=q*q; q4=q2*q2; q6=q4*q2; q8=q6*q2
      kfn2=kfn**2; kfn4=kfn2*kfn2
      Abig2=Abig**2; Abig3=Abig2*Abig; Abig4=Abig3*Abig; Abig5=Abig4*Abig
      Bp2=Bp**2; Bp3=Bp2*Bp; Bp4= Bp3*Bp; Bp5= Bp4*Bp; Bp6= Bp5*Bp; Bp7= Bp6*Bp; Bp8= Bp7*Bp

      alpha1(0)=kfn/(4*pi2*Abig)
      !
      alpha1(1)=-kfn*(Abig*(Abig*kfn2*hbarc2+Cp)-4*q2*Bp2*hbarc2)/(16*pi2*Abig2*q2*Bp*hbarc2)
      !
      alpha1(2)=kfn*(Abig*(Abig*kfn2*hbarc2-9*Cp)+12*q2*Bp2*hbarc2)/(24*pi2*Abig3*q2*hbarc2)
      !
      alpha1(3)=-kfn*(3*Abig2*(Abig*kfn2*hbarc2+Cp)**2-4*Abig*q2*Bp2*hbarc2*(Abig*kfn2*hbarc2-6*Cp)-48*q4*Bp4*hbarc4)/ &
           (144*pi2*Abig3*q4*Bp2*hbarc4)
      !
      alpha1(4)=kfn*(-3*Abig2*(2*Abig*kfn2*hbarc2-Cp)*(Abig*kfn2*hbarc2+Cp)+4*Abig*q2*Bp2*hbarc2*(2*Abig*kfn2*hbarc2-21*Cp) &
           +96*q4*Bp4*hbarc4)/(144*pi2*Abig4*q4*Bp*hbarc4)
      !
      alpha1(5)=kfn*(3*Abig4*kfn4*hbarc4+5*(-3*Abig2*Cp*(Abig*kfn2*hbarc2-5*Cp)+4*Abig*q2*Bp2*hbarc2*(Abig*kfn2*hbarc2-15*Cp) &
           +48*q4*Bp4*hbarc4))/(180*pi2*Abig5*q4*hbarc4)
      !
      alpha1(6)=-kfn*(4*Abig3*kfn2*q2*Bp2*hbarc4*(Abig*kfn2*hbarc2+Cp)+3*Abig3*(Abig*kfn2*hbarc2+Cp)**3 &
           -16*Abig*q4*Bp4*hbarc4*(Abig*kfn2*hbarc2-9*Cp)-192*q6*Bp6*hbarc6)/(384*pi2*Abig4*q6*Bp3*hbarc6)
      !
      alpha1(7)=kfn*(-20*Abig4*kfn2*q2*Bp2*hbarc4*(Abig*kfn2*hbarc2+Cp)**2-15*Abig4*(Abig*kfn2*hbarc2+Cp)**4 &
           +16*Abig2*q4*Bp4*hbarc4*(3*Abig2*kfn4*hbarc4-10*Abig*kfn2*Cp*hbarc2+30*Cp**2) &
           +320*Abig*q6*Bp6*hbarc6*(Abig*kfn2*hbarc2-12*Cp)+3840*q8*Bp8*hbarc8)/(4800*pi2*Abig5*q8*Bp4*hbarc8)
      !
      alpha1(8)=kfn*(-15*Abig3*(3*Abig*kfn2*hbarc2-Cp)*(Abig*kfn2*hbarc2+Cp)**2 &
           +4*Abig2*q2*Bp2*hbarc2*(9*Abig2*kfn4*hbarc4+5*Cp*(24*Cp-7*hbarc2*Abig*kfn2)) &
           +240*Abig*q4*Bp4*hbarc4*(Abig*kfn2*hbarc2-13*Cp)+2880*q6*Bp6*hbarc6)/(2880*pi2*Abig5*q6*Bp2*hbarc6)

      alpha1(:)=alpha1(:)*hbar_c
      return
    end subroutine ReA1
!!$    !
    !!  **** alpha 2  *****
    !
    subroutine ReA2(alpha2) ! alpha2(i)    [MeV^2]
      implicit none

      real ( kind = pr ) :: alpha2(0:8)
      real ( kind = pr ) :: kfp2,kfp4,Bn2,Bn3,Bn4,Bn5,Bn6,Bn7,Bn8
      real ( kind = pr ) :: Abig2,Abig3,Abig4,Abig5
      real ( kind = pr ) :: q2,q4,q6,q8

      q2=q*q; q4=q2*q2; q6=q4*q2; q8=q6*q2
      kfp2=kfp*kfp;  kfp4=kfp2*kfp2

      Bn2=Bn*Bn; Bn3=Bn2*Bn; Bn4=Bn3*Bn; Bn5=Bn4*Bn
      Bn6=Bn5*Bn; Bn7=Bn6*Bn; Bn8=Bn7*Bn

      Abig2=Abig*Abig; Abig3=Abig2*Abig; Abig4=Abig3*Abig; Abig5=Abig4*Abig

      alpha2(0)=-kfp/(4*pi2*Abig)
      !
      alpha2(1)=kfp*(Abig2*kfp2*hbarc2+4*q2*Bn*(Abig-Bn)*hbarc2+Abig*Cn)/(16*pi2*Abig2*q2*Bn*hbarc2) 
      !
      alpha2(2)=kfp*(-3*Abig3*kfp2*hbarc2+Bn*hbarc2*(12*q2*Bn*(Abig-Bn)-Abig2*(kfp2+6*q2))-3*Abig*Cn*(Abig-3*Bn))/ &
           (24*pi2*Abig**3*q2*Bn*hbarc2)
      !
      alpha2(3)=kfp*(-18*hbarc2*Abig2*q2*Bn*(Abig*kfp2*hbarc2+Cn)+3*Abig2*(Abig*kfp2*hbarc2+Cn)**2-4*Abig*hbarc2*q2*Bn2*(Abig*(kfp2+9*q2)*hbarc2-6*Cn) &
           +72*Abig*q4*Bn3*hbarc4-48*q4*Bn4*hbarc4)/(144*pi2*Abig3*q4*Bn2*hbarc4)
      !
      alpha2(4)=kfp*(-6*Abig3*(Abig*kfp2*hbarc2+Cn)**2+2*Abig2*hbarc2*q2*Bn2*(7*Abig*kfp2*hbarc2+18*Abig*q2*hbarc2-51*Cn) &
           +3*Abig2*Bn*(Abig*kfp2*hbarc2+Cn)*(2*Abig*kfp2*hbarc2+9*Abig*q2*hbarc2-Cn))/(144*pi2*Abig4*q4*Bn2*hbarc4) &
           +kfp*(-4*Abig*q2*Bn3*hbarc2*(2*Abig*kfp2*hbarc2+27*Abig*q2*hbarc2-21*Cn)+168*Abig*q4*Bn4*hbarc4-96*q4*Bn5*hbarc4)/ &
           (144*pi2*Abig4*q4*Bn2*hbarc4)
      !
      alpha2(5)=kfp*(15*Abig4*(Abig*kfp2*hbarc2+Cn)**2-15*Abig3*Bn*(Abig*kfp2*hbarc2+Cn)*(2*Abig*kfp2*hbarc2+3*Abig*q2*hbarc2-Cn) &
           +20*Abig2*q2*hbarc2*Bn3*(2*Abig*kfp2*hbarc2+9*Abig*q2*hbarc2-21*Cn))/(180*pi2*Abig5*q4*Bn2*hbarc4) &
           +kfp*(-Abig2*Bn2*(Abig2*hbarc4*(3*kfp4+35*kfp2*q2+45*q4)-15*Abig*Cn*hbarc2*(kfp2+17*q2)+75*Cn**2) &
           -20*Abig*q2*hbarc2*Bn4*(Abig*hbarc2*(kfp2+21*q2)-15*Cn)+480*Abig*q4*hbarc4*Bn5-240*q4*Bn6*hbarc4)/ &
           (180*pi2*Abig5*q4*hbarc4*Bn2)
      !
      alpha2(6)=kfp*(-24*Abig3*q2*Bn*hbarc2*(Abig*kfp2*hbarc2+Cn)**2+4*Abig3*hbarc4*q2*Bn2*(kfp2+18*q2)*(Abig*kfp2*hbarc2+Cn) &
           +3*Abig3*(Abig*kfp2*hbarc2+Cn)**3+32*Abig2*q4*hbarc4*Bn3*(Abig*hbarc2*(kfp2+3*q2)-6*Cn))/ &
           (384*pi2*Abig4*q6*Bn3*hbarc6)&
           +kfp*(-16*Abig*q4*Bn4*hbarc4*(Abig*hbarc2*(kfp2+18*q2)-9*Cn)+384*Abig*q6*hbarc6*Bn5-192*q6*hbarc6*Bn6)/ &
           (384*pi2*Abig4*q6*hbarc6*Bn3)
      !
      alpha2(7)=kfp*(-150*Abig4*q2*hbarc2*Bn*(Abig*kfp2*hbarc2+Cn)**3+20*Abig4*q2*hbarc4*Bn2*(kfp2+30*q2)*(Abig*kfp2*hbarc2+Cn)**2 &
           -200*Abig4*q4*hbarc6*Bn3*(kfp2+6*q2)*(Abig*kfp2*hbarc2+Cn)+15*Abig4*(Abig*kfp2*hbarc2+Cn)**4)/ &
           (4800*pi2*Abig5*q8*hbarc8*Bn4) &
           +kfp*(800*Abig2*q6*hbarc6*Bn5*(Abig*(kfp2+6*q2)*hbarc2-9*Cn)-16*Abig2*hbarc4*q4*Bn4*(Abig2*(3*kfp4+50*kfp2*q2+75*q4)*hbarc4 &
           -10*Abig*Cn*(kfp2+30*q2)*hbarc2+30*Cn**2))/(4800*pi2*Abig5*q8*Bn4*hbarc8) &
           +kfp*(-320*Abig*q6*hbarc6*Bn6*(Abig*(kfp2+30*q2)*hbarc2-12*Cn)+9600*Abig*q8*hbarc8*Bn7-3840*q8*Bn8*hbarc8) &
           /(4800*pi2*Abig5*q8*hbarc8*Bn4)
      !
      alpha2(8)=kfp*(-45*Abig4*(Abig*kfp2*hbarc2+Cn)**3-60*Abig3*q2*Bn2*hbarc2*(Abig*kfp2*hbarc2+Cn)*(5*Abig*kfp2*hbarc2+12*hbarc2*Abig*q2-2*Cn) &
           +15*Abig3*Bn*(Abig*kfp2*hbarc2+Cn)**2*(3*Abig*kfp2*hbarc2+20*hbarc2*Abig*q2-Cn) )/(2880*pi2*Abig5*hbarc6*q6*Bn3) &
           +kfp*(80*Abig2*q4*hbarc4*Bn4*(7*Abig*kfp2*hbarc2+36*Abig*q2*hbarc2-69*Cn)-4*Abig2*q2*hbarc2*Bn3*(Abig2*hbarc4* &
           (9*kfp4+130*kfp2*q2+180*q4)+5*Cn*(-7*Abig*hbarc2*kfp2-174*hbarc2*Abig*q2+24*Cn)))/ &
           (2880*pi2*Abig5*hbarc6*q6*Bn3) &
           +kfp*(-240*Abig*q4*hbarc4*Bn5*(Abig*hbarc2*(kfp2+26*q2)-13*Cn)+6720*Abig*q6*hbarc6*Bn6-2880*q6*Bn7*hbarc6)/ &
           (2880*pi2*Abig5*hbarc6*q6*Bn3)   

      alpha2(:)=alpha2(:)*hbar_c
      return
    end subroutine ReA2
    !
    !
    subroutine ReB1(beta1) ! beta1(i)  ! [MeV^2] 

      implicit none

      real ( kind = pr ) :: beta1(0:8)
      real ( kind = pr ) :: q2,q4,q6,q8,SqDelta
      real ( kind = pr ) :: Abig2,Abig3,Abig4,Abig5,Abig6
      real ( kind = pr ), parameter :: pi=2.d0*asin(1.d0)
      real ( kind = pr ), parameter :: qpi2=4.0_pr*pi**2

      q2=q*q; q4=q2*q2; q6=q4*q2; q8=q6*q2

      if(Delta.le.0.0_pr)then
         SqDelta=sqrt(-Delta) ! pure numnber
      else
         SqDelta=0.0_pr
      endif
      Abig2=Abig*Abig; Abig3=Abig2*Abig; Abig4=Abig3*Abig; Abig5=Abig4*Abig; Abig6=Abig5*Abig ! [MeV^(2,3,4,5,6)]

      beta1(0)=-SqDelta/(qpi2*Abig2)
      beta1(1)=-SqDelta*Bp/(qpi2*Abig3)
      beta1(2)=-SqDelta*(Abig*Cp+2*Delta)/(qpi2*Abig4*q2*hbarc2)
      beta1(3)=-SqDelta*(3*Abig*Cp+4*Delta)/(3*qpi2*Abig4*q2*hbarc2)
      beta1(4)=-SqDelta*Bp*(3*Abig*Cp+8*Delta)/(3*qpi2*Abig5*q2*hbarc2)
      beta1(5)=-SqDelta*(Abig*Cp+4*Delta)*(3*Abig*Cp+4*Delta)/(3*qpi2*Abig6*q4*hbarc4)
      beta1(6)=-SqDelta*Bp*(Abig*Cp+2*Delta)/(qpi2*Abig5*q2*hbarc2)
      beta1(7)=-SqDelta*(5*Abig*Cp*(Abig*Cp+4*Delta)+16*Delta**2)/(5*qpi2*Abig6*q4*hbarc4)
      beta1(8)=-SqDelta*(Abig*Cp*(3*Abig*Cp+14*Delta)+12*Delta**2)/(3*qpi2*Abig6*q4*hbarc4)


      !beta1(:)= beta1(:)/hbar_c

      return
    end subroutine ReB1
    !
    subroutine ReB2(beta2) ! beta2(i)  ! [MeV^2]  

      implicit none

      real ( kind = pr ) :: beta2(0:8)
      real ( kind = pr ) :: q2,q4,q6,SqDelta
      real ( kind = pr ) :: Abig2,Abig3,Abig4,Abig5,Abig6
      !real ( kind = pr ), parameter :: pi=2.d0*asin(1.d0)
      !real ( kind = pr ), parameter :: qpi2=4.0_pr*pi**2

      q2=q*q; q4=q2*q2; q6=q4*q2
      if(Delta.le.0.0_pr)then
         SqDelta=sqrt(-Delta) ! pure numnber
         Abig2=Abig*Abig; Abig3=Abig2*Abig; Abig4=Abig3*Abig; Abig5=Abig4*Abig; Abig6=Abig5*Abig ! [MeV^(2,3,4,5,6)]

         beta2(0)=+SqDelta/(qpi2*Abig2)
         beta2(1)=+SqDelta*(Bn-Abig)/(qpi2*Abig3)
         beta2(2)=+SqDelta*(Abig2*q2*hbarc2-2*Abig*hbarc2*q2*Bn+Abig*Cn+2*Delta)/(qpi2*Abig4*q2*hbarc2)
         beta2(3)=+SqDelta*(3*Abig2*q2*hbarc2-6*Abig*q2*hbarc2*Bn+3*Abig*Cn+4*Delta)/(3*qpi2*Abig4*q2*hbarc2)
         beta2(4)=-SqDelta*(Abig-Bn)*(3*Abig2*q2*hbarc2-6*Abig*q2*Bn*hbarc2+3*Abig*Cn+8*Delta)/(3*qpi2*Abig5*q2*hbarc2)
         beta2(5)=+SqDelta*(3*Abig4*q4*hbarc4-4*Abig*q2*hbarc2*Bn*(3*Abig2*q2*hbarc2+3*Abig*Cn+8*Delta) &
              +ABig*Cn*(18*Abig2*q2*hbarc2+3*Abig*Cn+16*Delta)+28*Abig2*Delta*q2*hbarc2      &
              +16*Delta**2)/(3*qpi2*Abig6*q4*hbarc4)
         beta2(6)=-SqDelta*(Abig-Bn)*(Abig2*q2*hbarc2-2*Abig*q2*Bn*hbarc2+Abig*Cn+2*Delta)/(qpi2*Abig5*q2*hbarc2)
         beta2(7)=+SqDelta*(5*Abig4*q4*hbarc4-20*Abig*q2*Bn*hbarc2*(Abig2*q2*hbarc2+Abig*Cn+2*Delta) &
              +5*Abig*Cn*(6*Abig2*q2*hbarc2+Abig*Cn+4*Delta)+40*Abig2*Delta*q2*hbarc2+16*Delta**2)/&
              (5*qpi2*Abig6*q4*hbarc4)
         beta2(8)=(-3*Delta*(Abig4*q4*hbarc4+4*Delta**2)+Abig*(4*Delta*q2*Bn*hbarc2*(3*Abig2*q2*hbarc2+3*Abig*Cn+7*Delta) &
              +Cn*(Abig*Cn*(26*Abig2*q2*hbarc2-22*Abig*Cn-47*Delta)+34*Abig2*Delta*q2*hbarc2-36*Delta**2) &
              +Bn**4*(-26*Abig*q6*hbarc6+22*q4*Cn*hbarc4)))/(3*qpi2*Abig6*q4*hbarc4*SqDelta)
      else
         beta2(:)=0.0_pr
      endif
      !      beta2(:)=beta2(:)/hbar_c

      return
    end subroutine ReB2
    !
    !!
    !
    subroutine ReD1(delta1) ! delta1(i)  [MeV^2] 

      implicit none

      real ( kind = pr ) :: delta1(0:8)
      real ( kind = pr ) :: q2,q4,q6,SqDelta
      real ( kind = pr ) :: Abig2,Abig3,Abig4,Abig5,Abig6
      !real ( kind = pr ), parameter :: pi=2.d0*asin(1.d0)
      !real ( kind = pr ), parameter :: qpi2=4.0_pr*pi**2

      q2=q*q; q4=q2*q2; q6=q4*q2
      if(delta.ge.0.0_pr)then
         SqDelta=sqrt(abs(Delta))

         Abig2=Abig*Abig; Abig3=Abig2*Abig; Abig4=Abig3*Abig; Abig5=Abig4*Abig; Abig6=Abig5*Abig ! [MeV^(2,3,4,5,6)]

         delta1(0)=SqDelta/(2*qpi2*Abig2)
         delta1(1)=SqDelta*Bp/(2*qpi2*Abig3)
         delta1(2)=SqDelta*(Abig*Cp+2*Delta)/(2*qpi2*Abig4*q2*hbarc2)
         delta1(3)=SqDelta*(3*Abig*Cp+4*Delta)/(6*qpi2*Abig4*q2*hbarc2)
         delta1(4)=SqDelta*Bp*(3*Abig*Cp+8*Delta)/(6*qpi2*Abig5*q2*hbarc2)
         delta1(5)=SqDelta*(Abig*Cp+4*Delta)*(3*Abig*Cp+4*Delta)/(6*qpi2*Abig6*q4*hbarc4)
         delta1(6)=SqDelta*Bp*(Abig*Cp+2*Delta)/(2*qpi2*Abig5*q2*hbarc2)
         delta1(7)=SqDelta*(5*Abig*Cp*(Abig*Cp+4*Delta)+16*Delta**2)/(10*qpi2*Abig6*hbarc4*q4)
         delta1(8)=SqDelta*(Abig*Cp*(3*Abig*Cp+14*Delta)+12*Delta**2)/(6*qpi2*Abig6*hbarc4*q4)
      else
         SqDelta=0.0_pr
      endif
      !      delta1(:)=delta1(:)/hbar_c

      return
    end subroutine ReD1
    !
    !!
    !
    subroutine ReD2(delta2) ! delta2(i) [MeV^2]   

      implicit none

      real ( kind = pr ) :: delta2(0:8)
      real ( kind = pr ) :: q2,q4,q6,SqDelta
      real ( kind = pr ) :: Abig2,Abig3,Abig4,Abig5,Abig6
      !real ( kind = pr ), parameter :: pi=2.d0*asin(1.d0)
      !real ( kind = pr ), parameter :: qpi2=4.0_pr*pi**2

      Abig2=Abig*Abig; Abig3=Abig2*Abig; Abig4=Abig3*Abig; Abig5=Abig4*Abig; Abig6=Abig5*Abig ! [MeV^(2,3,4,5,6)]
      q2=q*q; q4=q2*q2; q6=q4*q2
      if(Delta.ge.0.0_pr)then
         SqDelta=sqrt(abs(Delta))

         delta2(0)=-SqDelta/(2*qpi2*Abig2)
         delta2(1)=+SqDelta*(Abig-Bn)/(2*qpi2*Abig3)
         delta2(2)=-SqDelta*(Abig2*q2*hbarc2-2*Abig*q2*Bn*hbarc2+Abig*Cn+2*Delta)/(2*qpi2*Abig4*q2*hbarc2)
         delta2(3)=-SqDelta*(3*Abig2*q2*hbarc2-6*Abig*q2*Bn*hbarc2+3*Abig*Cn+4*Delta)/(6*qpi2*Abig4*q2*hbarc2)
         delta2(4)=+SqDelta*(Abig-Bn)*(3*Abig2*q2*hbarc2-6*Abig*q2*Bn*hbarc2+3*Abig*Cn+8*Delta)/(6*qpi2*Abig**5*q2*hbarc2)
         delta2(5)=-SqDelta*(3*Abig4*q4*hbarc4-4*Abig*q2*Bn*hbarc2*(3*Abig2*q2*hbarc2+3*Abig*Cn+8*Delta) &
              +Abig*Cn*(18*Abig2*q2*hbarc2+3*Abig*Cn+16*Delta)+28*Abig2*Delta*q2*hbarc2+16*Delta**2)/&
              (6*qpi2*Abig6*q4*hbarc4)
         delta2(6)=+(SqDelta*(Abig-Bn))*(Abig2*q2*hbarc2-2.*Abig*q2*Bn*hbarc2+Abig*Cn+2*Delta)/(2*qpi2*Abig5*q2*hbarc2)
         delta2(7)=-SqDelta*(5*Abig4*q4*hbarc4-20*hbarc2*Abig*q2*Bn*(Abig2*q2*hbarc2+Abig*Cn+2*Delta)+5*Abig*Cn*(6*Abig2*q2*hbarc2 &
              +Abig*Cn+4*Delta)+40*Abig2*Delta*q2*hbarc2+16*Delta**2)/(10*qpi2*Abig6*q4*hbarc4)
         delta2(8)=-(3*Delta*(Abig4*q4*hbarc4+4*Delta**2)+Abig*(-4*Delta*q2*hbarc2*Bn*(3*Abig2*q2*hbarc2+3*Abig*Cn+7*Delta) &
              +Cn*(Abig*Cn*(-26*Abig2*q2*hbarc2+22*Abig*Cn+47*Delta)-34*Abig2*Delta*q2*hbarc2+36*Delta**2) &
              +Bn**4*(26*Abig*q6*hbarc6-22*hbarc4*q4*Cn)))/(6*qpi2*Abig6*SqDelta*q4*hbarc4)
      else
         delta2(:)=0.0_pr
      endif
      !      delta2(:)=delta2(:)/hbar_c

    end subroutine ReD2
    !
    !!
    !
    subroutine ReG1(gamma1) ! gamma1(i)   [MeV^2]

      implicit none

      real ( kind = pr ) :: gamma1(0:8)
      real ( kind = pr ) :: q2,q3,q4,q5,q6,q7,q8,q9,q10
      real ( kind = pr ) :: Abig2,Abig3,Abig4,Abig5,Abig6
      real ( kind = pr ) :: kfn2,kfn4,kfn6,Bp2,Bp3,Bp4,Bp5,Bp6,Bp8,Bp10

      kfn2=kfn*kfn; kfn4=kfn2*kfn2; kfn6=kfn4*kfn2
      q2=q*q; q3=q2*q; q4=q2*q2; q5=q4*q; q6=q4*q2
      q7=q6*q; q8=q6*q2; q9=q8*q; q10=q8*q2
      Abig2=Abig*Abig; Abig3=Abig2*Abig; Abig4=Abig3*Abig; Abig5=Abig4*Abig; Abig6=Abig5*Abig ! [MeV^(-2,-3,-4,-5,-6)]
      Bp2=Bp*Bp; Bp3=Bp2*Bp; Bp4=Bp3*Bp; Bp5=Bp4*Bp; Bp6=Bp5*Bp ! [MeV^-2,-3,....]
      Bp8=Bp6*Bp2; Bp10=Bp6*Bp4

      gamma1(0)=(2*q2*Bp**2*hbarc2-Abig*(Abig*kfn2*hbarc2+Cp))/(4*qpi2*Abig2*q*Bp*hbar_c)
      gamma1(1)=-(Abig2*(Abig*kfn2*hbarc2+Cp)**2+4*Abig*q2*hbarc2*Bp2*Cp-8*q4*Bp4*hbarc4)/(16*qpi2*Abig3*q3*Bp2*hbarc3)
      gamma1(2)=(-Abig4*kfn4*hbarc4+Abig2*Cp**2-8*Abig*q2*Bp2*hbarc2*Cp+8*q4*Bp4*hbarc4)/(8*qpi2*Abig4*q3*Bp*hbarc3)
      gamma1(3)=-(Abig3*(Abig*kfn2*hbarc2+Cp)**3+24*Abig*q4*Bp4*hbarc4*Cp-32*q6*Bp6*hbarc6)/(48*qpi2*Abig4*q5*Bp3*hbarc5)
      gamma1(4)=(-Abig3*(2*Abig*kfn2*hbarc2-Cp)*(Abig*kfn2*hbarc2+Cp)**2+12*Abig2*q2*hbarc2*Bp2*Cp**2 &
           -72*Abig*q4*hbarc4*Bp4*Cp+64*q6*hbarc6*Bp6)/(48*qpi2*Abig5*q5*hbarc5*Bp2)
      gamma1(5)=-(Abig6*kfn6*hbarc6+Abig3*Cp**3-18*Abig2*q2*Bp2*hbarc2*Cp**2+48*Abig*q4*hbarc4*Bp4*Cp &
           -32*q6*Bp6*hbarc6)/(12*qpi2*Abig6*q5*hbarc5*Bp)
      gamma1(6)=(-Abig4*(Abig*kfn2*hbarc2+Cp)**4+16*Abig2*q4*hbarc4*Bp4*Cp**2-128*Abig*q6*Bp6*hbarc6*Cp+128*q8*Bp8*hbarc8)/&
           (128*qpi2*Abig5*q7*hbarc7*Bp4)
      gamma1(7)=(-Abig5*(Abig*kfn2*hbarc2+Cp)**5+160*Abig2*hbarc6*q6*Bp6*Cp**2-640*Abig*q8*hbarc8*Bp8*Cp+512*q10*Bp10*hbarc10)/&
           (320*qpi2*Abig6*hbarc9*q9*Bp5)
      gamma1(8)=(-Abig4*(3*Abig*kfn2*hbarc2-Cp)*(Abig*kfn2*hbarc2+Cp)**3+144.*Abig2*q4*hbarc4*Bp4*Cp**2 &
           -512.*Abig*q6*hbarc6*Bp6*Cp+384.*q8*Bp8*hbarc8)/(192.0_pr*qpi2*Abig6*q7*hbarc7*Bp3)


      !      gamma1(:)=gamma1(:)/hbar_c
      return
    end subroutine ReG1
    !
    !!
    !
    subroutine ReG2(gamma2) ! gamma2(i)  [MeV^2] 

      implicit none

      real ( kind = pr ) :: gamma2(0:8)
      real ( kind = pr ) :: q2,q3,q4,q5,q6,q7,q8,q9,q10
      real ( kind = pr ) :: kfp2,kfp4,Bn2,Bn3,Bn4,Bn5,Bn6,Bn7,Bn8,Bn9,Bn10
      !real ( kind = pr ), parameter :: pi=2.d0*asin(1.d0)
      !real ( kind = pr ), parameter :: qpi2=4.0_pr*pi**2
      real ( kind = pr ) :: Abig2,Abig3,Abig4,Abig5,Abig6

      q2=q*q; q3=q2*q; q4=q2*q2; q5=q4*q; q6=q4*q2
      q7=q6*q; q8=q6*q2; q9=q8*q; q10=q8*q2
      kfp2=kfp**2; kfp4=kfp2*kfp2
      Bn2=Bn**2; Bn3=Bn2*Bn; Bn4=Bn3*Bn; Bn5=Bn4*Bn; Bn6=Bn5*Bn
      Bn7=Bn6*Bn; Bn8=Bn7*Bn; Bn9=Bn8*Bn; Bn10=Bn9*Bn
      Abig2=Abig*Abig; Abig3=Abig2*Abig; Abig4=Abig3*Abig; Abig5=Abig4*Abig; Abig6=Abig5*Abig ! [MeV^(2,3,4,5,6)]

      gamma2(0)=(Abig*(Abig*kfp2*hbarc2+Cn)-2*q2*Bn2*hbarc2)/(4*qpi2*Abig2*hbar_c*q*Bn)
      !
      gamma2(1)=(-4*Abig2*q2*hbarc2*Bn*(Abig*kfp2*hbarc2+Cn)+Abig2*(Abig*kfp2*hbarc2+Cn)**2 &
           +4*Abig*q2*Bn2*hbarc2*Cn+8*Abig*q4*Bn3*hbarc4-8*q4*Bn4*hbarc4)/(16*qpi2*Abig3*q3*Bn2*hbarc3)
      !
      gamma2(2)=-(Abig3*(Abig*kfp2*hbarc2+Cn)**2-Abig2*Bn*(Abig2*kfp2*hbarc4*(kfp2+2*q2)+2*hbarc2*Abig*q2*Cn &
           -Cn**2)+4*Abig2*hbarc2*q2*Bn2*Cn+4*Abig*q2*Bn3*hbarc2*(Abig*q2*hbarc2-2*Cn) &
           -8*Abig*q4*Bn4*hbarc4+8*q4*Bn5*hbarc4)/(8*qpi2*Abig4*q3*hbarc3*Bn2)
      !
      gamma2(3)=(12*Abig3*q4*Bn2*hbarc4*(Abig*kfp2*hbarc2+Cn)-6*Abig3*q2*hbarc2*Bn*(Abig*kfp2*hbarc2+Cn)**2 &
           +Abig3*(Abig*kfp2*hbarc2+Cn)**3-24*Abig2*q4*Bn3*Cn*hbarc4-24*Abig*q4*Bn4*hbarc4*(Abig*q2*hbarc2-Cn) &
           +48*Abig*q6*Bn5*hbarc6-32*q6*Bn6*hbarc6)/(48*qpi2*Abig4*q5*Bn3*hbarc5)
      !
      gamma2(4)=(-2*Abig4*(Abig*kfp2*hbarc2+Cn)**3+Abig3*Bn*(Abig*kfp2*hbarc2+Cn)**2*(2*Abig*kfp2*hbarc2+9*hbarc2*Abig*q2-Cn) &
           +12*Abig2*q2*hbarc2*Bn3*Cn*(3*Abig*q2*hbarc2-Cn)+24*Abig2*hbarc4*q4*Bn4*(Abig*q2*hbarc2-4*Cn)) &
           /(48*qpi2*Abig5*q5*hbarc5*Bn3)&
           +(-6*Abig3*q2*hbarc2*Bn2*(Abig2*kfp2*hbarc4*(kfp2+2*q2)+2*Abig*q2*Cn*hbarc2-Cn**2) &
           -72*Abig*q4*hbarc4*Bn5*(Abig*q2*hbarc2-Cn)+112*Abig*q6*Bn6*hbarc6-64*q6*Bn7*hbarc6)&
           /(48*qpi2*Abig5*hbarc5*q5*Bn3)
      !
      gamma2(5)=(Abig5*(Abig*kfp2*hbarc2+Cn)**3-Abig4*Bn*(Abig*kfp2*hbarc2+Cn)**2*(2*Abig*kfp2*hbarc2+3*hbarc2*Abig*q2-Cn) &
           -12*Abig3*hbarc2*q2*Bn3*Cn*(Abig*q2*hbarc2-Cn)+24*Abig2*q4*hbarc4*Bn5*(Abig*q2*hbarc2-3*Cn))/(12*qpi2*Abig6*hbarc5*q5*Bn3) &
           +(-6*Abig2*q2*hbarc2*Bn4*(Abig2*q4*hbarc4-8*Abig*q2*hbarc2*Cn+3*Cn**2)+Abig3*Bn2*(Abig3*kfp2*hbarc6*(kfp4+3*kfp2*q2+3*q4) &
           +3*Abig2*q4*hbarc4*Cn-3*Abig*q2*hbarc2*Cn**2+Cn**3)-8*Abig*q4*hbarc4*Bn6*(7*Abig*q2*hbarc2-6*Cn)) &
           /(12*qpi2*Abig6*q5*hbarc5*Bn3) + (64*Abig*q6*hbarc6*Bn7-32*q6*hbarc6*Bn8)/(12*qpi2*Abig6*q5*hbarc5*Bn3)
      !
      gamma2(6)=(-32*Abig4*q6*hbarc6*Bn3*(Abig*kfp2*hbarc2+Cn)+24*Abig4*q4*hbarc4*Bn2*(Abig*kfp2*hbarc2+Cn)**2 &
           -8*Abig4*q2*Bn*hbarc2*(Abig*kfp2*hbarc2+Cn)**3 &
           +Abig4*(Abig*kfp2*hbarc2+Cn)**4 &
           +64*Abig2*q6*hbarc6*Bn5*(Abig*q2*hbarc2-3.*Cn) &
           +16*Abig2*q4*hbarc4*Bn4*Cn*(6.0_pr*Abig*q2*hbarc2-Cn)-64*Abig*q6*Bn6*hbarc6*(3.0_pr*Abig*q2*hbarc2-2.0_pr*Cn) &
           +256*Abig*q8*hbarc8*Bn7-128*q8*hbarc8*Bn8)/(128.*qpi2*Abig5*q7*hbarc7*Bn4)
      !
      gamma2(7)=(80*Abig5*q8*hbarc8*Bn4*(Abig*kfp2*hbarc2+Cn)-80*abig5*q6*hbarc6*Bn3*(Abig*kfp2*hbarc2+Cn)**2 &
           +40*Abig5*q4*hbarc4*Bn2*(Abig*kfp2*hbarc2+Cn)**3-10*Abig5*q2*hbarc2*Bn*(Abig*kfp2*hbarc2+Cn)**4 &
           +Abig5*(Abig*kfp2*hbarc2+Cn)**5)/(320*qpi2*Abig6*q9*hbarc9*Bn5) &
           +(-160*Abig3*q6*hbarc6*Bn5*Cn*(2*Abig*q2*hbarc2-Cn)+640*Abig2*q8*hbarc8*Bn7*(Abig*q2*hbarc2-2*Cn))/&
           (320*qpi2*Abig6*q9*hbarc9*Bn5) &
           +(-160*Abig2*q6*hbarc6*Bn6*(Abig2*q4*hbarc4-6*Abig*q2*Cn*hbarc2+Cn**2)-640*Abig*q8*hbarc8*Bn8* &
           (2*Abig*q2*hbarc2-Cn)+1280*Abig*q10*hbarc10*Bn9-512*q10*hbarc10*Bn10)/ &
           (320*qpi2*Abig6*q9*hbarc9*Bn5)
      !
      gamma2(8)=(-3*Abig5*(Abig*kfp2*hbarc2+Cn)**4-8*Abig4*q2*hbarc2*Bn2*(Abig*kfp2*hbarc2+Cn)**2*(2*Abig*hbarc2*&
           (kfp2+3*q2)-Cn)+Abig4*Bn*(Abig*kfp2*hbarc2+Cn)**3*(3*Abig*kfp2*hbarc2+20*Abig*q2*hbarc2-Cn) &
           -48*Abig3*q4*hbarc4*Bn4*Cn*(4*Abig*q2*hbarc2-3*Cn) &
           +192*Abig2*q6*hbarc6*Bn6*(2*Abig*q2*hbarc2-5*Cn)-48*Abig2*q4*hbarc4*Bn5*(2*Abig2*q4*hbarc4 &
           -14*Abig*q2*hbarc2*Cn+3*Cn**2)+24*Abig4*q4*hbarc4*Bn3*(Abig2*kfp2*hbarc4*(kfp2+2.*q2) &
           +2.*Abig*q2*hbarc2*Cn-Cn**2) &
           -64*Abig*q6*hbarc6*Bn7*(13.*Abig*q2*hbarc2-8*Cn)+896*Abig*q8*hbarc8*Bn8 &
           -384*q8*hbarc8*Bn9)/(192*qpi2*Abig6*q7*hbarc7*Bn4)           


      !      gamma2(:)=gamma2(:)/hbar_c

      return
    end subroutine ReG2
    !
  end subroutine betapnTzero


  !
  !!
  !

  subroutine betapnallNUMERICTZERO(omega,msn,msp,qMeV,Un,Up,kfn,kfp,betapn)
    !-------------------------------------------
    ! Abstract: subroutine that calculates the
    ! beta_i^pn functions at T=0 (half analytic/half numeric)
    !
    ! input:
    ! omega: energy [MeV]
    ! msn: effective mass [MeV]
    ! q: transfer momentum [MeV]
    ! T: temperature [MeV]
    ! mun: chemical potential [MeV]
    !
    ! output
    ! betapn(0-8): beta functions for fiexd
    !              omega and q
    !-------------------------------------------

    implicit none
    integer                       :: i,j
    real ( kind = pr ),intent(in) :: omega,msn,msp,qMeV,kfn,kfp,Un,Up
    complex ( kind = pr )         :: betapn(0:8)
    real ( kind = pr ), parameter :: limdm=0.1_pr
    real ( kind = pr ), parameter :: T=0.1_pr
    real ( kind = pr ), parameter :: xcut=7.0_pr
    real ( kind = pr ), parameter :: eps=0.1_pr
    real ( kind = pr )            :: efn,efp,deltam,w
    real ( kind = pr )            :: betaIM(0:8),betaR(0:8)
    real ( kind = pr )            :: q2,xf
    real ( kind = pr )            :: q3,q4,q6,q8,sumTIM(0:8)
    real ( kind = pr )            :: betaIMn(0:8),betaIMp(0:8)
    real ( kind = pr )            :: elimP,elimM,xr1(npt),xr2(npt),xr3(npt)
    real ( kind = pr )            :: wr1(npt),wr2(npt),wr3(npt),omegap
    real ( kind = pr )            :: elim_neu(0:1),elim_pro(0:1)
    real ( kind = pr )            :: Abig,Bp,Bn,Cp,Cn,Delta,rhon,rhop,kn,kp,yy
    real ( kind = pr )            :: limGup,limGdown,limGtup,limGtdown
    real ( kind = pr )            :: xfp,xfn,q

    efn=((kfn*hbar_c)**2/(2.d0*msn)) ! [Mev]
    efp=((kfp*hbar_c)**2/(2.d0*msp)) ! [Mev]
    deltam=msp-msn ! [Mev]
    w=omega-(Up-Un) ! [Mev]

    rhon=(kfn/tpi13)**3
    rhop=(kfp/tpi13)**3
    yy=(rhon-rhop)/(rhon+rhop)


    q=qMeV/hbar_c! [fm]
    kn=q/(2*kfn)
    kp=q/(2*kfp)


    betapn(:)=cz
    betaIM(:)=0.0_pr
    betaR(:)=0.0_pr

    q2=q*q   ! MeV^2
    q3=q2*q  ! MeV^3
    q4=q2*q2 ! MeV^4
    q6=q3*q3 ! MeV^6
    q8=q4*q4 ! MeV^8
    xf=msn*msp/(qpi*q) ! [MeV]


    Bp=1.0_pr/(2*msp) ; Bn=1.0_pr/(2*msn) ! [MeV^-1]
    w=omega-Up+Un ! [MeV]
    Abig=Bn-Bp ! [MeV^-1]
    Cp=w-Bp*qMeV**2 ; Cn=w+Bn*qMeV**2 ![MeV]
    Delta=2*Bn*Bp*(0.5_pr*qMeV**2-w*(msp-msn)) ! pure number

    xfn=1.0_pr/(8*pi*qMeV*Bn) ; xfp=1.0_pr/(8*pi*qMeV*Bp)

    !--------------------------------------------
    !           Immginary parts
    !--------------------------------------------

    betaIMn(:)=0.0_pr ; betaIMp(:)=0.0_pr

    if(Delta.ge.0.0_pr)then
       if(abs(yy).le.0.05_pr)then! limit Y=0
          limGup =kfn*hbar_c
          limGdown =kfn*hbar_c*abs(omega*msn/(qMeV*kfn*hbar_c)-kn)
          limGtup=kfn*hbar_c
          limGtdown=kfp*hbar_c*abs(omega*msp/(qMeV*kfp*hbar_c)+kp)
          Abig=eps ! Not 0 for numerical reasons

          !stop 'The code has not been tested yet!'
       else! .... normal
          limGup=Min(kfn*hbar_c, (qMeV*Bp+sqrt(Delta))/dabs(Abig)) ! [MeV]
          limGdown=abs((qMeV*Bp-sqrt(Delta))/Abig)     ! [MeV]

          limGtup=Min(kfp*hbar_c, (qMeV*Bn+sqrt(Delta))/dabs(Abig)) ! [MeV]
          limGtdown=abs((qMeV*Bn-sqrt(Delta))/Abig)     ! [MeV]  
       endif
       do i=0,8 ! MeV^2
          BetaIM(i)=-xfp*(Gbig(i,limGup)-Gbig(i,limGdown))*ftheta(limGup-limGdown) &
               +xfn*(Gtildebig(i,limGtup)-Gtildebig(i,limGtdown))*ftheta(limGtup-limGtdown) 
       enddo
    endif


    !--------------------------------------------
    !          Real parts 
    !--------------------------------------------

    ! --- neutrons
    elim_neu(0)=qMeV**2/(2.d0*msn)+q*kfn*hbarc2/msn
    elim_neu(0)=(1.0d0+xcut*T/efn)*elim_neu(0)
    elim_neu(1)=-(qMeV**2/(2.d0*msn)+q*kfn*hbarc2/msn)
    elim_neu(1)=(1.0d0+xcut*T/efn)*elim_neu(1)

    ! --- protons
    elim_pro(0)=qMeV**2/(2.d0*msp)+q*kfp*hbarc2/msp
    elim_pro(0)=(1.0d0+xcut*T/efp)*elim_pro(0)
    elim_pro(1)=-(qMeV**2/(2.d0*msp)+q*kfp*hbarc2/msp)
    elim_pro(1)=(1.0d0+xcut*T/efp)*elim_pro(1)


    elimP=max(elim_neu(0),elim_pro(0))
    elimM=min(elim_neu(1),elim_pro(1)) !min


    CALL gset(elimM,    elimP,    npt,xr1,wr1)
    CALL gset(elimM,    omega-eps,npt,xr2,wr2)
    CALL gset(omega+eps,elimP,    npt,xr3,wr3)

    betaR(:)= 0.0_pr

    do j=1,npt
       omegap=xr1(j)
       w=omegap-(Up-Un) ! [Mev]
       Cp=w-Bp*qMeV**2 ; Cn=w+Bn*qMeV**2 ![MeV]
       Delta=2*Bn*Bp*(0.5_pr*qMeV**2-w*(msp-msn)) ! pure number

       !           Immginary parts
       sumTIM(:)=0.0_pr
       if(Delta.ge.0.0_pr)then
          if(abs(yy).le.0.05_pr)then! limit Y=0
             limGup =kfn*hbar_c
             limGdown =kfn*hbar_c*abs(omega*msn/(qMeV*kfn*hbar_c)-kn)
             limGtup=kfn*hbar_c
             limGtdown=kfp*hbar_c*abs(omega*msp/(qMeV*kfp*hbar_c)+kp)
             Abig=eps ! Not 0 for numerical reasons

!             stop 'The code has not been tested yet!'
          else! .... normal
             limGup=Min(kfn*hbar_c, (qMeV*Bp+sqrt(Delta))/dabs(Abig)) ! [MeV]
             limGdown=abs((qMeV*Bp-sqrt(Delta))/Abig)     ! [MeV]

             limGtup=Min(kfp*hbar_c, (qMeV*Bn+sqrt(Delta))/dabs(Abig)) ! [MeV]
             limGtdown=abs((qMeV*Bn-sqrt(Delta))/Abig)     ! [MeV]  
          endif
          do i=0,8 ! MeV^2
             sumTIM(i)=-xfp*(Gbig(i,limGup)-Gbig(i,limGdown))*ftheta(limGup-limGdown) &
                  +xfn*(Gtildebig(i,limGtup)-Gtildebig(i,limGtdown))*ftheta(limGtup-limGtdown) 
          enddo
       endif


       betaR(:)= betaR(:) +(sumTIM(:)-betaIM(:))/(xr1(j)-omega)*wr1(j) &
            +betaIM(:)/(xr2(j)-omega)*wr2(j)+betaIM(:)/(xr3(j)-omega)*wr3(j)

    enddo

    !-------------------------------------
    !        Result
    !-------------------------------------
    betapn(:)=Cmplx(betaR(:)/pi,betaIM(:), kind =pr ) ! MeV^2

  contains

    function Gbig(i,k) result(res) ! in [MeV^2]
      implicit none
      integer :: i
      real ( kind = pr ) :: res,k! [MeV]
      real ( kind = pr ) :: aux    
      select case(i)

      case (0)
         res=k**2/2
      case(1)
         res=(Abig*k**4+2*k**2*Cp)/(8*qMeV**2*Bp)
      case(2)
         res=k**4/(4*qMeV**2)
      case(3)
         res=(Abig**2*k**6+3*Abig*k**4*Cp+3*k**2*Cp**2)/(24*qMeV**4*Bp**2)
      case(4)
         res=(2*Abig*k**6+3*k**4*Cp)/(24*qMeV**4*Bp)
      case(5)
         res=k**6/(6*qMeV**4)
      case(6)
         res=((Abig*k**2+Cp)**4-Cp**4)/(64*Abig*qMeV**6*Bp**3)
      case(7)
         res=((Abig*k**2+Cp)**5-Cp**5)/(160*Abig*qMeV**8*Bp**4)
      case(8)
         aux=3*Abig**2*k**8+8*Abig*k**6*Cp+6*k**4*Cp**2
         res=aux/(96*qMeV**6*Bp**2)

      case default
         stop 'Error in function Gbig!'
      end select

      return
    end function Gbig
    !
    !!
    !

    function Gtildebig(i,k) result(res) ! in [MeV^2]
      implicit none
      integer :: i
      real ( kind = pr ) :: res,k! [MeV]
      real ( kind = pr ) :: k2,k4,k6,k8
      real ( kind = pr ) :: q2,q4,q6,q8,aux,aux2

      q2=qMeV**2; q4=q2*q2; q6=q4*q2; q8=q6*q2
      k2=k*k; k4=k2*k2; k6=k4*k2; k8=k6*k2
      res=0.0_pr
      select case(i)
      case(0)
         res=k2/2
      case(1)
         res=k2*(Abig*k2-4*q2*Bn+2*Cn)/(8*q2*Bn)
      case(2)
         res=k2*(-Abig*k2+Bn*(k2+2*q2)-2*Cn)/(4*q2*Bn)
      case(3)
         res=(Abig**2*k6+3*Abig*k4*(Cn-2*q2*Bn)+3*k2*(Cn-2*q2*Bn)**2)/(24*q4*Bn**2)
      case(4)
         aux=2*Abig**2*k4-Bn*(Abig*k2*(2*k2+9*q2)+3*Cn*(k2+6*q2))+6*Cn*(Abig*k2+Cn)+6*q2*Bn**2*(k2+2*q2)
         res=-k2*aux/(24*q4*Bn**2)
      case(5)
         aux=3*k2*(Bn-Abig)*(q2*Bn-Cn)+k4*(Abig-Bn)**2+3*(Cn-q2*Bn)**2
         res=k2*aux/(6*q4*Bn**2)
      case(6)
         aux=(Abig*k2-2*q2*Bn+Cn)**4-(Cn-2*q2*Bn)**4
         res=aux/(64*Abig*q6*Bn**3)
      case(7)
         aux=(Abig**4*k8-5*(2*q2*Bn-Cn)*(Abig*k2-2*q2*Bn+Cn)*(Abig**2*k4-(2*q2*Bn-Cn)*(Abig*k2-2*q2*Bn+Cn)))
         res=k2*aux/(160*q8*Bn**4)
      case(8)
         aux=3*Abig**2*k8*(Abig-Bn)+4*Abig*k6*(Cn*(3*Abig-2*Bn)+q2*Bn*(4*Bn-5*Abig))
         aux2=6*k4*(2*q2*Bn-Cn)*(Cn*(Bn-3*Abig)-2*q2*Bn*(Bn-2*Abig)) -12.0_pr*k2*(q2*Bn-Cn)*((Cn-2.0_pr*q2*Bn)**2)
         res=-(aux+aux2)/(96*q6*Bn**3) ! attenzione al segno meno in aux2
      case default
         stop 'Error in function Gtildebig!'
      end select
      return
    end function Gtildebig
    !
  end subroutine betapnallNUMERICTZERO

  !
  !!
  !=======  END =============

end module beta
