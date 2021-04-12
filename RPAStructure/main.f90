module main

  use asym_cstes ! module constants
  use asym       ! module main code
  use asym_io    ! module io
  use omp_lib
  

implicit none
contains

double precision function anmfp(rho,Y,TMeV,nom,type,readforce,energy)
  double precision,          intent(in) :: rho,Y,TMeV,energy
  character ( len = 10 ),       intent(in) :: nom
  character ( len = 100 ),      intent(in) :: type
  logical,                      intent(in) :: readforce
  double precision                         :: avgnmfp

  integer:: num_threads
  double precision:: start,finish
  

  ! setting the coupling constants
  call set_lr_cstes
    if(readforce) then
       call read_force(nom)
    else
       call read_functional(nom)
    endif
    call set_coeff_functional(trluc(type),readforce)
!  !$ start=omp_get_wtime()
  call asymmetric(rho,y,TMev,avgnmfp,energy)         
!  !$ write(*,*)avgnmfp,energy,int((omp_get_wtime()-start)/60),'minutes'
anmfp=avgnmfp
end function anmfp

end module main


