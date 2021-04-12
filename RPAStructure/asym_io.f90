!------------------------------------------------------------------------
!
! Useful routines to open/close files
!
!------------------------------------------------------------------------
module asym_io
  use asym_cstes

  implicit none

  integer, parameter, private :: fdim = 100, id0 = 10
  integer, private :: current_id = 0, current_index = 0
  integer, dimension(fdim), private :: files_num = 0
  !
  !! Files related stuff
  !
  character ( len = 250 ) :: out, inp, wfo, wfc, ppt
  character ( len =  60 ) :: wrt_pot, rd_pot
  character ( len =   3 ) :: pre
  logical :: io_debug = .false.

  interface open_file
     module procedure open_file
  end interface

  interface close_file
     module procedure close_file
  end interface


  interface error
     module procedure error
  end interface


  contains

 subroutine open_file( id, name, status, form )
    implicit none
    integer, intent(out) :: id
    character ( len = * ), intent(in), optional :: name, status, form
    logical :: ok, yes
    character ( len = 11 ) :: frm
    !
    frm = "formatted"
    if ( present(form) ) then
       if ( len(form) <= 11 ) frm = form
    end if
    !
    if ( current_id == 0 ) then
       current_id = id0 ! look for the first availlable id
    end if
    !
    do
       inquire( unit = current_id, opened = ok )
       if ( .not.ok ) exit
       current_id = current_id + 1
       if ( current_id >= 100 )  &
            call error( "open_file", "Too many files opened..." )
    end do

  !
    current_index = current_index + 1
    files_num(current_index) = current_id
    if ( present(name) ) then
       if ( present(status) ) then
          if ( trluc(status) == "OLD" ) then
             inquire( file = name, exist = yes )
             if ( .not.yes ) then
                print '(" *** File ",a," does not exist !...")', name
                stop
             end if
          end if
          open( unit = current_id, file = name, status = status,  &
               form = trl(frm) )
       else
          open( unit = current_id, file = name, status = "unknown",  &
               form = trl(frm) )
       end if
    else
       if ( present(status) ) then
          open( unit = current_id, status = status )
       else
          open( unit = current_id, status = "unknown",  &
               form = trl(frm) )
       end if
    end if
    id = current_id
   !
    if ( io_debug )  &
         print '(1x,">>>>> file number ",i2," has been opened.")', id
    !
  end subroutine open_file

  !
  !!
  !
  subroutine close_file(id)
    implicit none
    integer, intent(in) :: id
    character ( len = 3 ) :: num
    logical :: ok, ok_ok
    integer :: i
    !
    write( num, '(i3)' ) id
    inquire( unit = id, opened = ok )
    ok_ok = .false.
    do i = 1, fdim
       if ( files_num(i) == id ) then
          ok_ok = .true.
          exit
       end if
    end do
    if ( .not.( ok .and. ok_ok ) )  &
         call error( "close_file", trl(num)//" is not opened..." )
    close( files_num(i) )
    if ( io_debug )  &
         print '(1x,"<<<<< file number ",i2," has been closed.")', files_num(i)
    if ( i == fdim ) then
       files_num(i) = 0
    else
       files_num(i:fdim-1) = files_num(i+1:fdim)
    end if
    !
   if ( current_index > 1 ) then
       current_index = current_index - 1
       current_id = files_num(current_index)
    else
       current_id = 0
       current_index = 0
    end if
    !
  end subroutine close_file

  !
  !!
  !
subroutine error( name, reason )
    implicit none
    character ( len = * ), intent(in), optional :: name, reason
    !
    print '(50("*"))'
    if ( present(name) ) print '(" Error in: ",a)', name
    if ( present(reason) ) print '(1x,a)', reason
    print '(50("*"))'
    stop
    !
  end subroutine error


end module asym_io
