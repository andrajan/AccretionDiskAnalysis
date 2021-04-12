module determinante

  use asym_cstes ! module constants


contains
  subroutine cramer(Amatrix,Bmatrix,Amatrixp,Bmatrixp,size,Chinn,Chinp,Chipn,Chipp) 
    implicit none
    !-----------------------------------------------------------------------
    ! abstract: subroutine that calculates the response functions
    !
    ! input
    ! size: size of the matrix
    ! Amatrix,Amatrixp: square matrix that contains the interaction
    ! Bmatrix,Bmatrixp: column vector that cotains the sources
    !
    ! output
    ! Chinn,Chinp,Chipn,Chipp: response functions
    !-----------------------------------------------------------------------

    integer                         :: size    
    complex ( kind = pr ),intent(in):: Amatrix(1:8,1:8),Bmatrix(1:8)
    complex ( kind = pr ),intent(in):: Amatrixp(1:8,1:8),Bmatrixp(1:8)
    complex ( kind = pr )           :: Chinn,Chinp,Chipn,Chipp

    complex ( kind = pr )           :: Num(1:8,1:8),Deno(1:8,1:8)


    !
    !---- building the matrix
    !

    Num(:,:)=Amatrix(:,:)
    Deno(:,:)=Amatrix(:,:)

    !   calculating chi^nn
    Num(:,1)=Bmatrix(:)

    call det(Num,Deno,size,Chinn)

    Num(:,:)=Amatrix(:,:)
    Deno(:,:)=Amatrix(:,:)

    if(size.eq.6)then  !   calculating chi^pn
       Num(:,4)=Bmatrix(:)
    elseif(size.eq.8)then
       Num(:,5)=Bmatrix(:)
    elseif(size.eq.1)then
       Num(:,1)=1.0_pr
    endif

    call det(Num,Deno,size,Chipn)

    Num(:,:)=Amatrixp(:,:)
    Deno(:,:)=Amatrixp(:,:)

    !   calculating chi^pp
    Num(:,1)=Bmatrixp(:)

    call det(Num,Deno,size,Chipp)

    Num(:,:)=Amatrixp(:,:)
    Deno(:,:)=Amatrixp(:,:)

    if(size.eq.6)then  !   calculating chi^np
       Num(:,4)=Bmatrixp(:)
    elseif(size.eq.8)then
       Num(:,5)=Bmatrixp(:)
    elseif(size.eq.1)then
       Num(:,1)=1.0_pr
    endif

    call det(Num,Deno,size,Chinp)



    return
  end subroutine cramer

  !
  !!
  !

  subroutine cramerNP(Amatrix,Bmatrix,Amatrixp,Bmatrixp,size,Chinp,Chipn) 
    implicit none
    !-----------------------------------------------------------------------
    ! abstract: subroutine that calculates the response functions
    !
    ! input
    ! size: size of the matrix
    ! Amatrix,Amatrixp: square matrix that contains the interaction
    ! Bmatrix,Bmatrixp: column vector that cotains the sources
    !
    ! output
    ! Chinp,Chipn: response functions
    !-----------------------------------------------------------------------

    integer                         :: size    
    complex ( kind = pr ),intent(in):: Amatrix(1:8,1:8),Bmatrix(1:8)
    complex ( kind = pr ),intent(in):: Amatrixp(1:8,1:8),Bmatrixp(1:8)
    complex ( kind = pr )           :: Chinp,Chipn
    complex ( kind = pr )           :: Num(1:8,1:8),Deno(1:8,1:8)

    !
    !---- building the matrix
    !

    Num(:,:)=Amatrix(:,:)
    Deno(:,:)=Amatrix(:,:)

    !   calculating chi^pn
    Num(:,1)=Bmatrix(:)
    !Chipn=(0.0_pr,0.0_pr)
    

    call det(Num,Deno,size,Chipn)


    !   calculating chi^np
    Num(:,:)=Amatrixp(:,:)
    Deno(:,:)=Amatrixp(:,:)


    Num(:,1)=Bmatrixp(:)

    call det(Num,Deno,size,Chinp)



    return
  end subroutine cramerNP
  !
  !!
  !

  subroutine det(Num,Deno,size,Chi)
    implicit none
    integer              :: size,indx(size),i,j
    real ( kind = pr )   :: d,deg
    complex ( kind = pr ):: Num(1:8,1:8),Deno(1:8,1:8)
    complex ( kind = pr ):: Chi,zero
    complex ( kind = pr ):: determinante(0:1),T1(size,size),T2(size,size)
  
    zero=(0.0_pr,0.0_pr)
    deg=2.0_pr

    determinante(:)=CMPLX(1.0_pr,0.0_pr, kind =pr)
    do i=1,size
       do j=1,size
          T1(i,j)=Num(i,j)
          T2(i,j)=Deno(i,j)
       enddo
    enddo

    call cludcmp(T1,size,size,indx,d)


    do i=1,size
       determinante(0)=determinante(0)*T1(i,i)
    enddo
    ! we check the parity of permutations 
    determinante(0)=determinante(0)*Cmplx(d,0.0_pr, kind =pr)

    !
    call cludcmp(T2,size,size,indx,d)

    do i=1,size
       determinante(1)=determinante(1)*T2(i,i)
    enddo
    ! we check the parity of permutations 
    determinante(1)=determinante(1)*Cmplx(d,0.0_pr, kind =pr)

    !write(*,*)determinante(0),determinante(1)
    if(determinante(0).eq.zero.and.determinante(1).eq.zero)then
       Chi=0.0_pr
    else
       Chi=-deg*determinante(0)/determinante(1)/pi
    endif
    return
  end subroutine det


!


    !=======================================================================
  SUBROUTINE cludcmp(a,n,np,indx,d)
    !========================================================================
    !
    !  Numerical Recipes : LU decomposition                          
    !
    ! Given a matrix a(1:n,1:n), with physical dimension np by np, 
    ! this routine replaces it by
    ! the LU decomposition of a rowwise permutation of itself. a and n are 
    ! input. a is output,
    ! arranged as in equation (2.3.14) above; indx(1:n) is an output vector
    ! that records the
    ! row permutation effected by the partial pivoting; d is output as ±1 
    ! depending on whether
    ! the number of row interchanges was even or odd, respectively. 
    ! This routine is used in
    ! combination with lubksb to solve linear equations or invert a matrix.
    !==========================================================================
    !
    INTEGER                     :: n,np,indx(n)
    real ( kind = pr )          :: d
    complex ( kind = pr )       :: a(np,np),sum,dum2
    INTEGER, PARAMETER          :: NMAX=500
    real ( kind = pr ),PARAMETER:: TINY=1.0e-20
    INTEGER                     :: i,imax,j,k
    real ( kind = pr )          :: aamax,dum,vv(NMAX)

    d=1.0_pr
    DO  i=1,n
       aamax=0.
       DO  j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
       ENDDO
       if (aamax.eq.0.) then
          !write(*,*) 'singular matrix in ludcmp'
          a(:,:)=Cmplx(0.0_pr,0.0_pr, kind =pr)
          return
       endif
       vv(i)=1.0_pr/aamax
    ENDDO
    DO  j=1,n
       DO  i=1,j-1
          sum=a(i,j)
          DO  k=1,i-1
             sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
       enddo
       aamax=0.0_pr
       DO  i=j,n
          sum=a(i,j)
          DO  k=1,j-1
             sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
             imax=i
             aamax=dum
          endif
       enddo
       if (j.ne.imax)then
          DO  k=1,n
             dum2=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum2
          enddo
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if(a(j,j).eq.0.)a(j,j)=TINY
       if(j.ne.n)then
          dum2=1./a(j,j)
          DO  i=j+1,n
             a(i,j)=a(i,j)*dum2
          enddo
       endif
    enddo
    RETURN
  END subroutine cludcmp

 end module determinante
