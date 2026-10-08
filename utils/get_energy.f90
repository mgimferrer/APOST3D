!! ********************************************************************** !!
!! program: get_energy                                                    !!
!! purpose: Gaussian 09 version of get_energy_g16 (to be deprecated):    !!
!! reads the reference energy components of a finished calculation from  !!
!! its .log file (route with #P, iop(3/33=3) and Pop=Full) and prints     !!
!! them in .fchk format, to be appended to the .fchk for ENPART: kinetic, !!
!! electron-nuclear and electron-electron energies from the last          !!
!! "KE= ... PE= ... EE=" line, plus the ECP integral matrix               !!
!! when the calculation uses pseudopotentials. Every error goes to stderr !!
!! with exit code 1 and nothing on stdout, so a wrong .log never adds     !!
!! anything to the .fchk.                                                 !!
!! usage: get_energy name.log >> name.fchk                                !!
!! author: MGimf                                                          !!
!! ********************************************************************** !!
program get_energy
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  integer, parameter :: dp=kind(1.0d0)
  character(len=:), allocatable :: logname
  character(len=256) :: line
  character(len=43) :: title
  integer :: iu,ios,lname,nbas,ntri,i,j,imat
  logical :: normal,have_scf,have_ke,have_ecp,ecp_used,ok
  real(dp) :: ekin,epe,eee
  real(dp), allocatable :: xecp(:)

  if(command_argument_count().ne.1) &
    call fail('usage: get_energy name.log >> name.fchk')
  call get_command_argument(1,length=lname)
  allocate(character(len=lname) :: logname)
  call get_command_argument(1,logname)
  open(newunit=iu,file=logname,status='old',action='read',iostat=ios)
  if(ios.ne.0) call fail('cannot open '//logname)

!! one pass: the last KE/PE/EE line, the basis size and the last ECP     !!
!! integral matrix; the .log must end normally after the SCF             !!
  normal=.false.
  have_scf=.false.
  have_ke=.false.
  have_ecp=.false.
  ecp_used=.false.
  nbas=0
  ntri=0
  do
    read(iu,'(a)',iostat=ios) line
    if(ios.ne.0) exit
    if(index(line,'SCF Done:').gt.0) then
      have_scf=.true.
      normal=.false.
    else if(index(line,' KE=').eq.1) then
      call value_after(line,'KE=',ekin,have_ke)
      if(have_ke) call value_after(line,'PE=',epe,have_ke)
      if(have_ke) call value_after(line,'EE=',eee,have_ke)
    else if(nbas.eq.0.and.index(line,'NBasis').gt.0) then
      i=index(line,'NBasis')
      j=index(line(i:),'=')
      ios=1
      if(j.gt.0) read(line(i+j:),*,iostat=ios) nbas
      if(ios.ne.0) nbas=0
    else if(index(line,'Pseudopotential Parameters').gt.0) then
      ecp_used=.true.
    else if(index(line,'ECP Int').gt.0) then
!! the ECP matrix is the block with IMat=1; IMat=2-4 are the x, y, z     !!
!! spin-orbit components, which do not enter the energy                  !!
      i=index(line,'IMat=')
      if(i.gt.0) then
        j=index(line(i:),':')
        if(j.gt.0) line(i+j-1:i+j-1)=' '
        read(line(i+5:),*,iostat=ios) imat
        if(ios.ne.0) call fail('unreadable ECP integral header: '//trim(line))
        if(imat.ne.1) cycle
      end if
      if(nbas.eq.0) call fail('ECP integrals found before NBasis')
      ntri=nbas*(nbas+1)/2
      if(.not.allocated(xecp)) allocate(xecp(ntri))
      call read_lower(iu,nbas,xecp,ok)
      if(.not.ok) call fail('ECP integral matrix incomplete')
      have_ecp=.true.
    else if(index(line,'Normal termination').gt.0) then
      normal=.true.
    else if(index(line,'Error termination').gt.0) then
      normal=.false.
    end if
  end do
  close(iu)

  if(.not.have_scf) call fail('no "SCF Done" line in '//logname)
  if(.not.normal) call fail(logname//' did not end normally after the SCF')
  if(.not.have_ke) &
    call fail('no " KE= ... PE= ... EE=" line (the route needs #P and Pop=Full)')
  if(ecp_used.and..not.have_ecp) &
    call fail('pseudopotentials used but their integrals are not printed (add iop(3/33=3))')

  write(error_unit,'(a)') 'get_energy: Gaussian 09 version, to be deprecated (use get_energy_g16)'
  title='Kinetic Energy'
  write(*,'(A43,A6,ES22.15)') title,'R     ',ekin
  title='Electron-Nuclei Energy'
  write(*,'(A43,A6,ES22.15)') title,'R     ',epe
  title='Electron-Electron Energy'
  write(*,'(A43,A6,ES22.15)') title,'R     ',eee
  if(have_ecp) then
    title='ECP Matrix'
    write(*,'(A43,A6,I12)') title,'R   N=',ntri
    write(*,'(5ES16.8)') (xecp(i),i=1,ntri)
  end if

contains

!! ---- !!
!! value after the first occurrence of key on line (list-directed) !!
  subroutine value_after(line,key,x,ok)
    character(len=*), intent(in) :: line,key
    real(dp), intent(out) :: x
    logical, intent(out) :: ok
    integer :: k,ios
    k=index(line,key)
    ok=.false.
    if(k.eq.0) return
    read(line(k+len(key):),*,iostat=ios) x
    ok=(ios.eq.0)
  end subroutine value_after

!! ---- !!
!! symmetric matrix as Gaussian prints it (lower triangle in blocks of !!
!! five columns, a column-index line per block) into packed row order  !!
  subroutine read_lower(iu,n,x,ok)
    integer, intent(in) :: iu,n
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    character(len=256) :: line
    integer :: d,i,j,irow,ios
    ok=.false.
    do d=1,n,5
      read(iu,'(a)',iostat=ios) line
      if(ios.ne.0) return
      do i=d,n
        read(iu,'(a)',iostat=ios) line
        if(ios.ne.0) return
        read(line,*,iostat=ios) irow,(x(i*(i-1)/2+j),j=d,min(i,d+4))
        if(ios.ne.0.or.irow.ne.i) return
      end do
    end do
    ok=.true.
  end subroutine read_lower

!! ---- !!
  subroutine fail(msg)
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: msg
    write(error_unit,'(a)') 'get_energy: '//msg
    flush(error_unit)
    stop 1
  end subroutine fail

end program get_energy
