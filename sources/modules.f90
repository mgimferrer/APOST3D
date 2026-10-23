!! ********************************************************************** !!
!! F90 MODULES -- basis set, AO/density matrices, integration grid,       !!
!! per-atom EFO/NAO/QTAIM scratch arrays, timers, .inp keyword flags      !!
!!   basis_set         -- basis-set data (primitives, contraction         !!
!!                         coeffs, AO overlap S)                          !!
!!                         CONTAINS build_basis, do_overlap               !!
!!   ao_matrices       -- MO coefficient/density matrices (c/p/cb/ps/pa/  !!
!!                         pb/c_no/occ_no)                                !!
!!                         CONTAINS build_ao_matrices                     !!
!!   integration_grid  -- atom-centered grid (Nrad/Nang/pha/phb/rr00,     !!
!!                         Lebedev th/ph/w, radial wr/xr)                 !!
!!                         CONTAINS build_integration_grid                !!
!!   effao_mod         -- effective-AO population matrices (p0/p0net/     !!
!!                         p0gro/ip0)                                     !!
!!                         CONTAINS allocate_effao                        !!
!!   nao_mod           -- natural-AO matrices (unao/ssnao)                 !!
!!                         CONTAINS allocate_nao                          !!
!!   version_mod       -- version number and date (banner, results file)  !!
!!   elements_mod      -- element symbols for atom labels                 !!
!!   output_mod        -- print_value: 'label (unit) : value' lines       !!
!!   json_mod          -- the results file <job>.apost.json               !!
!!   timing_mod        -- CPU+wall-clock timers                           !!
!!                         CONTAINS get_wall_time, print_timer            !!
!!   input_options_mod -- ~70 .inp keyword flags parsed by read_input()   !!
!!                         data only, no CONTAINS                         !!
!! ********************************************************************** !!

!! ----------------------------------------------------------------------
!! module: output_mod
!! Helpers for the output style guide (docs/source/developer/
!! output-style.md). print_value prints one "label (unit) : value" line,
!! with the label padded to label_width so that the colons of neighbouring
!! lines line up; the unit is optional.
!!   call print_value('Kinetic energy', x, 7, 'au')  -- real, 7 decimals
!!   call print_value('Number of atoms', n)           -- integer
!!   call print_value('Functional', 'B3LYP')          -- text
!! author: MGimf
!! ----------------------------------------------------------------------
   MODULE output_mod
   implicit none
   integer, parameter :: label_width=40
   private
   public :: print_value, label_width

   interface print_value
     module procedure print_value_real, print_value_int, print_value_text
   end interface

   CONTAINS

!! ---- !!
   function padded_label(label,units) result(lab)
   character(len=*), intent(in) :: label
   character(len=*), intent(in), optional :: units
   character(len=:), allocatable :: lab
   lab=trim(label)
   if(present(units)) then
     if(len_trim(units).gt.0) lab=lab//' ('//trim(units)//')'
   end if
   if(len(lab).lt.label_width) lab=lab//repeat(' ',label_width-len(lab))
   end function padded_label

!! ---- !!
   subroutine print_value_real(label,x,ndec,units)
   character(len=*), intent(in) :: label
   real*8, intent(in) :: x
   integer, intent(in) :: ndec
   character(len=*), intent(in), optional :: units
   character(len=32) :: fmt
   write(fmt,'(a,i0,a,i0,a)') '(2x,a,'':'',f',max(15,ndec+9),'.',ndec,')'
   write(*,fmt) padded_label(label,units),x
   end subroutine print_value_real

!! ---- !!
   subroutine print_value_int(label,n,units)
   character(len=*), intent(in) :: label
   integer, intent(in) :: n
   character(len=*), intent(in), optional :: units
   write(*,'(2x,a,'': '',i0)') padded_label(label,units),n
   end subroutine print_value_int

!! ---- !!
   subroutine print_value_text(label,text,units)
   character(len=*), intent(in) :: label,text
   character(len=*), intent(in), optional :: units
   write(*,'(2x,a,'': '',a)') padded_label(label,units),trim(text)
   end subroutine print_value_text

   END MODULE output_mod

!! ***** !!

   MODULE basis_set
   integer :: numprim,nbasis,mmax,natoms
   INTEGER, ALLOCATABLE :: nlm(:,:),iptoat(:),nprimbas(:,:)
   INTEGER, ALLOCATABLE :: ihold(:),llim(:),iulim(:)
   INTEGER, ALLOCATABLE :: iptob_cartesian(:)
   REAL*8, ALLOCATABLE :: expp(:),coefpb(:,:),coord(:,:)
   REAL*8, ALLOCATABLE :: s(:,:),s12p(:,:),s12m(:,:)
   REAL*8 fact(0:30),fact2(0:20)
   DATA fact/1.0D0,1.0D0,2.0D0,6.0D0,24.0D0,120.0D0,720.0D0,5040.0d0,40320.0d0,362880.0d0,&
   3628800.0d0,39916800.0d0, 479001600.0d0, 6227020800.0d0, 87178291200.0d0, 1307674368000.0d0,&
   20922789888000.0d0,355687428096.0d3, 6402373705728.0d3, 121645100408832.0d3, 24329020081766.4d5,&
   510909421717094.d5, 112400072777760.d7, 258520167388849.d8, 620448401733239.d9, 155112100433309.d11,&
   403291461126605.d12, 108888694504183.d14, 304888344611713.d15, 884176199373970.d16, 265252859812191.d18 /
   DATA fact2/1.0d0,1.0d0, 2.0d0, 3.0d0, 8.0d0, 15.0d0, 48.0d0, 105.0d0, 384.0d0, 945.0d0, &
   3840.0d0, 10395.0d0, 46080.0d0, 135135.0d0, 645120.0d0, 2027025.0d0,10321920.0d0, &
   34459425.0d0, 185794560.0d0, 654729075.0d0, 3715891200.0d0/
! numprim -> number of cartesian primitives
! nbasis  -> number of basis functions
! mmax -> max number of primitives per basis function
! natoms -> number of atoms
! nlm(numprim,3) xyz exponents of primitives
! expp(numprim,) gaussian exponent of primitives
! coefpb(numprim,nbasis) primitive contraction coeff of basis functions
! nprimbas(numprim,mmax) list of primitives contributing to each basis function 
! iptoat(numprim) atom to which each primitive belongs
! iptob_cartesian(numprim) -> basis function to which each primitive belongs
! s(nbasis,nbasis) AO overlap matrix
! coord(3,natoms) atomic coordinates 
! ihold(nbasis) atom to which each basis function belongs
! llim,iulim(natoms) first and lad basis function of each atom

   CONTAINS

!! ********************************************************************* !!
!! subroutine: build_basis                                               !!
!! purpose: reads the basis set (shells, primitives, contraction coeffs) !!
!! from the .fchk (unit 15, already open) and builds the primitive/basis !!
!! function maps, pure-to-cartesian coefficients, and normalization used !!
!! by the rest of the code; prints the atom/basis/primitive counts.      !!
!! arguments: none (output via basis_set module's own arrays)            !!
!! author: PSalse                                                        !!
!! ********************************************************************* !!
   SUBROUTINE build_basis()
   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
   PARAMETER(PI=4.0d0*DATAN(1.0d0),TOL=1.0d-8)
   integer, allocatable :: mnsh(:),iatsh(:),mssh(:)
   real*8, allocatable :: expsh(:),c1(:),c2(:),xnorm(:),coefp(:)
   integer :: ncshell,npshell
   integer :: dummy,ilog
   DIMENSION mult(-5:5)
   DATA mult/11,9,7,5,4,1,3,6,10,15,21/

!! reading basis set info !!
   ncshell=int_locate(15,"Number of contract",ilog)
   npshell=int_locate(15,"Number of primi",ilog)
   nbasis=int_locate(15,"Number of basis",ilog)
   natoms=int_locate(15,"Number of atoms",ilog)
   allocate(coord(3,natoms))
   dummy=int_locate(15,"Current cartesian",ilog)
   read(15,*)(coord(1,i),coord(2,i),coord(3,i),i=1,natoms)

   allocate(mnsh(ncshell),iatsh(ncshell),mssh(ncshell))
   allocate(expsh(npshell),c1(npshell),c2(npshell))

   dummy=int_locate(15,"tives per she",ilog)
   read(15,*)(mnsh(i),i=1,ncshell)
   dummy=int_locate(15,"Shell to atom",ilog)
   read(15,*)(iatsh(i),i=1,ncshell)
   dummy=int_locate(15,"Primitive expo",ilog)
   read(15,*)(expsh(i),i=1,npshell)
   dummy=int_locate(15,"Contraction co",ilog)
   read(15,*)(c1(i),i=1,npshell)
   dummy=int_locate(15,"Shell types",ilog)
   read(15,*)(mssh(i),i=1,ncshell)
   dummy=0
   do i=1,ncshell
     if(mssh(i).eq.-1) dummy=1
   end do
   if(dummy.eq.1) then
     dummy=int_locate(15,"P(S=P) Cont",ilog)
     read(15,*)(c2(i),i=1,npshell)
   end if

!! processing basis set !!
   numprim=0
   nbasis=0
   do i=1,ncshell
     if(mssh(i).lt.-1) then
       numprim=numprim+mult(abs(mssh(i)))*mnsh(i)
     else
       numprim=numprim+mult(mssh(i))*mnsh(i)
     end if
     nbasis=nbasis+mult(mssh(i))
   end do

   allocate(ihold(nbasis),llim(natoms),iulim(natoms))
!! basis to atom map !!
   ii=0
   do i=1,ncshell
     do k=1,mult(mssh(i))
       ii=ii+1
       ihold(ii)=iatsh(i)
     end do
   end do
!! setting basis set limits for mulliken !!
   llim(1)=1
   iulim(natoms)=nbasis
   iat=1
   do i=1,nbasis
     if(ihold(i).ne.iat) then
       iulim(iat)=i-1
       llim(iat+1)=i
       iat=iat+1
     end if
   end do

   allocate (nlm(numprim,3),expp(numprim),iptoat(numprim),coefpb(numprim,nbasis),iptob_cartesian(numprim))
   allocate (coefp(numprim),xnorm(numprim))

!! angular momentum of primitives and primitive-to-atom map !!
   nlm=0
   icount=1
   do i=1,ncshell
     if(mssh(i).eq.0) then
       do j=1,mnsh(i)
         iptoat(icount)=iatsh(i)
         icount=icount+1
       end do
     else if (mssh(i).eq.1) then
       do j=1,mnsh(i)
         nlm(icount,1)=1
         nlm(icount+1,2)=1
         nlm(icount+2,3)=1
         do ii=0,2
           iptoat(icount+ii)=iatsh(i)
         end do
         icount=icount+3
       end do
     else if (mssh(i).eq.-1) then
       do j=1,mnsh(i)
         nlm(icount+1,1)=1
         nlm(icount+2,2)=1
         nlm(icount+3,3)=1
         do ii=0,3
           iptoat(icount+ii)=iatsh(i)
         end do
         icount=icount+4
       end do
     else if (abs(mssh(i)).eq.2) then
       do j=1,mnsh(i)
         nlm(icount,1)=2   !dxx
         nlm(icount+1,2)=2 !dyy
         nlm(icount+2,3)=2 !dzz
         nlm(icount+3,1)=1 !dxy
         nlm(icount+3,2)=1
         nlm(icount+4,1)=1 !dxz
         nlm(icount+4,3)=1
         nlm(icount+5,2)=1 !dyz
         nlm(icount+5,3)=1
         do ii=0,5
           iptoat(icount+ii)=iatsh(i)
         end do
         icount=icount+6
       end do
     else if (abs(mssh(i)).eq.3) then
       do j=1,mnsh(i)
         nlm(icount,1)=3    !fxxx
         nlm(icount+1,2)=3  !fyyy
         nlm(icount+2,3)=3  !fzzz
         nlm(icount+3,1)=1  !fxyy
         nlm(icount+3,2)=2
         nlm(icount+4,1)=2  !fxxy
         nlm(icount+4,2)=1
         nlm(icount+5,1)=2  !fxxz
         nlm(icount+5,3)=1
         nlm(icount+6,1)=1  !fxzz
         nlm(icount+6,3)=2
         nlm(icount+7,2)=1  !fyzz
         nlm(icount+7,3)=2
         nlm(icount+8,2)=2  !fyyz
         nlm(icount+8,3)=1
         nlm(icount+9,1)=1  !fxyz
         nlm(icount+9,2)=1
         nlm(icount+9,3)=1
         do ii=0,9
           iptoat(icount+ii)=iatsh(i)
         end do
         icount=icount+10
       end do
     else if (abs(mssh(i)).eq.4) then
       do j=1,mnsh(i)
         nlm(icount,3)=4    !   ZZZZ
         nlm(icount+1,2)=1  !   YZZZ
         nlm(icount+1,3)=3
         nlm(icount+2,2)=2  !   YYZZ
         nlm(icount+2,3)=2
         nlm(icount+3,2)=3  !   YYYZ
         nlm(icount+3,3)=1
         nlm(icount+4,2)=4  !   YYYY
         nlm(icount+5,1)=1  !   XZZZ
         nlm(icount+5,3)=3
         nlm(icount+6,1)=1  !   XYZZ
         nlm(icount+6,2)=1
         nlm(icount+6,3)=2
         nlm(icount+7,1)=1  !   XYYZ
         nlm(icount+7,2)=2
         nlm(icount+7,3)=1
         nlm(icount+8,1)=1  !   XYYY
         nlm(icount+8,2)=3
         nlm(icount+9,1)=2  !   XXZZ
         nlm(icount+9,3)=2
         nlm(icount+10,1)=2 !   XXYZ
         nlm(icount+10,2)=1
         nlm(icount+10,3)=1
         nlm(icount+11,1)=2 !   XXYY
         nlm(icount+11,2)=2
         nlm(icount+12,1)=3 !   XXXZ
         nlm(icount+12,3)=1
         nlm(icount+13,1)=3 !   XXXY
         nlm(icount+13,2)=1
         nlm(icount+14,1)=4 !   XXXX
         do ii=0,14
           iptoat(icount+ii)=iatsh(i)
         end do
         icount=icount+15
       end do
     else
       call apost_stop('angular momentum not implemented')
     end if
   end do

!! list primitive exponents and coefficients !!
   icount=1
   jcount=1
   do i=1,ncshell
     do j=1,mnsh(i)
       kk=mult(abs(mssh(i)))
       if(mssh(i).eq.-1) kk=4
       do k=1,kk
         expp(icount)=expsh(jcount)
         if(mssh(i).eq.-1.and.k.ne.1) then
           coefp(icount)=c2(jcount)
         else
           coefp(icount)=c1(jcount)
         end if
         icount=icount+1
       end do
       jcount=jcount+1
     end do
   end do

!! primitive normalization !!
   do i=1,numprim
     nn=nlm(i,1)
     ll=nlm(i,2)
     mm=nlm(i,3)
     fnn=fact(nn)/fact(2*nn)
     fll=fact(ll)/fact(2*ll)
     fmm=fact(mm)/fact(2*mm)
     xnorm(i)=(2.0d0*expp(i)/PI)**0.75d0*DSQRT((8.0d0*expp(i))**(nn+ll+mm)*fnn*fll*fmm)
   end do

!! generating primitive-to-orbital map, and pure-to-cartesian mapping   !!
!! up to G-type orbitals                                                !!
   numprim=0
   nbasis=0
   coefpb=0.0d0
   do i=1,ncshell
     do j=1,mnsh(i)
       if(mssh(i).ge.-1) then
         do k=1,mult(mssh(i))
           numprim=numprim+1
           coefpb(numprim,nbasis+k)=coefp(numprim)*xnorm(numprim)
           iptob_cartesian(numprim)=nbasis+k ! primitive-to-function map (Cartesian)
         end do
       else if (mssh(i).eq.-2) then ! mapping for pure 5d
         coefpb(numprim+1,nbasis+1)=coefp(numprim+1)*xnorm(numprim+1)*(-0.5d0)
         coefpb(numprim+1,nbasis+4)=coefp(numprim+1)*xnorm(numprim+1)*(sqrt(3.0d0)/2.0d0)
         coefpb(numprim+2,nbasis+1)=coefp(numprim+2)*xnorm(numprim+2)*(-0.5d0)
         coefpb(numprim+2,nbasis+4)=coefp(numprim+2)*xnorm(numprim+2)*(-sqrt(3.0d0)/2.0d0)
         coefpb(numprim+3,nbasis+1)=coefp(numprim+3)*xnorm(numprim+3)
         coefpb(numprim+4,nbasis+5)=coefp(numprim+4)*xnorm(numprim+4)
         coefpb(numprim+5,nbasis+2)=coefp(numprim+5)*xnorm(numprim+5)
         coefpb(numprim+6,nbasis+3)=coefp(numprim+6)*xnorm(numprim+6)
         numprim=numprim+mult(abs(mssh(i)))
       else if (mssh(i).eq.-3) then ! mapping for pure 7f
         coefpb(numprim+1,nbasis+2)=coefp(numprim+1)*xnorm(numprim+1)*(-sqrt(6.0d0)/4.0d0)
         coefpb(numprim+1,nbasis+6)=coefp(numprim+1)*xnorm(numprim+1)*(sqrt(10.0d0)/4.0d0)
         coefpb(numprim+2,nbasis+3)=coefp(numprim+2)*xnorm(numprim+2)*(-sqrt(6.0d0)/4.0d0)
         coefpb(numprim+2,nbasis+7)=coefp(numprim+2)*xnorm(numprim+2)*(-sqrt(10.0d0)/4.0d0)
         coefpb(numprim+3,nbasis+1)=coefp(numprim+3)*xnorm(numprim+3)
         coefpb(numprim+4,nbasis+2)=coefp(numprim+4)*xnorm(numprim+4)*(-sqrt(30.0d0)/20.0d0)
         coefpb(numprim+4,nbasis+6)=coefp(numprim+4)*xnorm(numprim+4)*(-3.0d0*sqrt(2.0d0)/4.0d0)
         coefpb(numprim+5,nbasis+3)=coefp(numprim+5)*xnorm(numprim+5)*(-sqrt(30.0d0)/20.0d0)
         coefpb(numprim+5,nbasis+7)=coefp(numprim+5)*xnorm(numprim+5)*(3.0d0*sqrt(2.0d0)/4.0d0)
         coefpb(numprim+6,nbasis+1)=coefp(numprim+6)*xnorm(numprim+6)*(-3.0d0*sqrt(5.0d0)/10.0d0)
         coefpb(numprim+6,nbasis+4)=coefp(numprim+6)*xnorm(numprim+6)*(sqrt(3.0d0)/2.0d0)
         coefpb(numprim+7,nbasis+2)=coefp(numprim+7)*xnorm(numprim+7)*(sqrt(30.0d0)/5.0d0)
         coefpb(numprim+8,nbasis+3)=coefp(numprim+8)*xnorm(numprim+8)*(sqrt(30.0d0)/5.0d0)
         coefpb(numprim+9,nbasis+1)=coefp(numprim+9)*xnorm(numprim+9)*(-3.0d0*sqrt(5.0d0)/10.0d0)
         coefpb(numprim+9,nbasis+4)=coefp(numprim+9)*xnorm(numprim+9)*(-sqrt(3.0d0)/2.0d0)
         coefpb(numprim+10,nbasis+5)=coefp(numprim+10)*xnorm(numprim+10)
         numprim=numprim+mult(abs(mssh(i)))
       else if (mssh(i).eq.-4) then ! mapping for pure 9f
         coefpb(numprim+1 ,nbasis+1)=coefp(numprim+1 )*xnorm(numprim+1 )
         coefpb(numprim+2 ,nbasis+3)=coefp(numprim+2 )*xnorm(numprim+2 )*(sqrt(70.0d0)/7.0d0)
         coefpb(numprim+3 ,nbasis+1)=coefp(numprim+3 )*xnorm(numprim+3 )*(-3.0d0*sqrt(105.0d0)/35.0d0)
         coefpb(numprim+3 ,nbasis+4)=coefp(numprim+3 )*xnorm(numprim+3 )*(-3.0d0*sqrt(21.0d0)/14.0d0)
         coefpb(numprim+4 ,nbasis+3)=coefp(numprim+4 )*xnorm(numprim+4 )*(-3.0d0*sqrt(70.0d0)/28.0d0)
         coefpb(numprim+4 ,nbasis+7)=coefp(numprim+4 )*xnorm(numprim+4 )*(-sqrt(10.0d0)/4.0d0)
         coefpb(numprim+5 ,nbasis+1)=coefp(numprim+5 )*xnorm(numprim+5 )*(3.0d0/8.0d0)
         coefpb(numprim+5 ,nbasis+4)=coefp(numprim+5 )*xnorm(numprim+5 )*(sqrt(5.0d0)/4.0d0)
         coefpb(numprim+5 ,nbasis+8)=coefp(numprim+5 )*xnorm(numprim+5 )*(sqrt(35.0d0)/8.0d0)
         coefpb(numprim+6 ,nbasis+2)=coefp(numprim+6 )*xnorm(numprim+6 )*(sqrt(70.0d0)/7.0d0)
         coefpb(numprim+7 ,nbasis+5)=coefp(numprim+7 )*xnorm(numprim+7 )*(3.0d0*sqrt(7.0d0)/7.0d0)
         coefpb(numprim+8 ,nbasis+2)=coefp(numprim+8 )*xnorm(numprim+8 )*(-3.0d0*sqrt(14.0d0)/28.0d0)
         coefpb(numprim+8 ,nbasis+6)=coefp(numprim+8 )*xnorm(numprim+8 )*(-3.0d0*sqrt(2.0d0)/4.0d0)
         coefpb(numprim+9 ,nbasis+5)=coefp(numprim+9 )*xnorm(numprim+9 )*(-sqrt(35.0d0)/14.0d0)
         coefpb(numprim+9 ,nbasis+9)=coefp(numprim+9 )*xnorm(numprim+9 )*(-sqrt(5.0d0)/2.0d0)
         coefpb(numprim+10,nbasis+1)=coefp(numprim+10)*xnorm(numprim+10)*(-3.0d0*sqrt(105.0d0)/35.0d0)
         coefpb(numprim+10,nbasis+4)=coefp(numprim+10)*xnorm(numprim+10)*(3.0d0*sqrt(21.0d0)/14.0d0)
         coefpb(numprim+11,nbasis+3)=coefp(numprim+11)*xnorm(numprim+11)*(-3.0d0*sqrt(14.0d0)/28.0d0)
         coefpb(numprim+11,nbasis+7)=coefp(numprim+11)*xnorm(numprim+11)*(3.0d0*sqrt(2.0d0)/4.0d0)
         coefpb(numprim+12,nbasis+1)=coefp(numprim+12)*xnorm(numprim+12)*(3.0d0*sqrt(105.0d0)/140.0d0)
         coefpb(numprim+12,nbasis+8)=coefp(numprim+12)*xnorm(numprim+12)*(-3.0d0*sqrt(3.0d0)/4.0d0)
         coefpb(numprim+13,nbasis+2)=coefp(numprim+13)*xnorm(numprim+13)*(-3.0d0*sqrt(70.0d0)/28.0d0)
         coefpb(numprim+13,nbasis+6)=coefp(numprim+13)*xnorm(numprim+13)*(sqrt(10.0d0)/4.0d0)
         coefpb(numprim+14,nbasis+5)=coefp(numprim+14)*xnorm(numprim+14)*(-sqrt(35.0d0)/14.0d0)
         coefpb(numprim+14,nbasis+9)=coefp(numprim+14)*xnorm(numprim+14)*(sqrt(5.0d0)/2.0d0)
         coefpb(numprim+15,nbasis+1)=coefp(numprim+15)*xnorm(numprim+15)*(3.0d0/8.0d0)
         coefpb(numprim+15,nbasis+4)=coefp(numprim+15)*xnorm(numprim+15)*(-sqrt(5.0d0)/4.0d0)
         coefpb(numprim+15,nbasis+8)=coefp(numprim+15)*xnorm(numprim+15)*(sqrt(35.0d0)/8.0d0)
         numprim=numprim+mult(abs(mssh(i)))
       end if
     end do
     nbasis=nbasis+mult(mssh(i))
   end do

!! max num prim per basis function !!
   mmax=0
   do i=1,ncshell
     ii=mnsh(i)
     if(mssh(i).le.-2) ii=3*ii
     if(mssh(i).le.-4) ii=2*ii
     if(ii.gt.mmax) mmax=ii
   end do
   mmax=mmax+1
   allocate(nprimbas(mmax,nbasis))

   nprimbas=0
   do i=1,nbasis
     npb=0
     do k=1,numprim
       if(abs(coefpb(k,i)).gt.TOL) then
         npb=npb+1
         nprimbas(npb,i)=k
       end if
     end do
   end do

!! calculating overlap matrix !!
   call do_overlap()

!! calculating S^1/2 and S^-1/2 !!
   allocate(s12p(nbasis,nbasis),s12m(nbasis,nbasis))
   call build_Smp(nbasis,s,s12m,s12p,0)

!! deallocating auxiliary arrays !!
   deallocate(coefp,xnorm)
   deallocate(mnsh,iatsh,mssh)
   deallocate(expsh,c1,c2)

   END SUBROUTINE build_basis

!! ***** !!

!! ********************************************************************* !!
!! subroutine: do_overlap                                                !!
!! purpose: primitive-Gaussian overlap matrix sp(numprim,numprim), via   !!
!! the closed-form Gaussian-product/binomial-expansion formula per       !!
!! Cartesian direction (screening primitive pairs that are exactly zero  !!
!! by symmetry along any direction), then contracted through coefpb      !!
!! into the basis-function AO overlap matrix S(nbasis,nbasis). Called    !!
!! once from build_basis().                                              !!
!! arguments: none (numprim/nbasis/coord/nlm/expp/coefpb via module,     !!
!! output S via module)                                                  !!
!! author: PSalse.                                                       !!
!! ********************************************************************* !!
   subroutine do_overlap()
   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
   PARAMETER(PI=4.0d0*DATAN(1.0d0),TOL=1.0d-8)
   dimension AminusB(3)
   real*8, allocatable ::sp(:,:)

!! Not worth OMP... !!
   allocate(sp(numprim,numprim))
   allocate(s(nbasis,nbasis))
   do ia=1,numprim
     do ib=1,ia
       sp(ia,ib)=0.0d0
       do ixyz=1,3 !! screening of primitives: skip pairs exactly zero by symmetry along this direction  !!
         AminusB(ixyz)=coord(ixyz,iptoat(ia))-coord(ixyz,iptoat(ib))
         ii=mod(nlm(ia,ixyz)+nlm(ib,ixyz),2)
         if(abs(AminusB(ixyz)).lt.TOL.and.ii.ne.0) go to 111
       end do
       gamma_p=expp(ia)+expp(ib)
       eta_p=(expp(ia)*expp(ib))/gamma_p
       do_ov=PI**(3.0d0/2.0d0)/gamma_p**(3.0d0/2.0d0)
       do ixyz=1,3 !! accumulate the overlap integral along all 3 directions !!
         do_ov=do_ov*fact(nlm(ia,ixyz))*fact(nlm(ib,ixyz))/(2.0d0**(nlm(ia,ixyz)+nlm(ib,ixyz)))
         if(abs(AminusB(ixyz)).gt.TOL) do_ov=do_ov*exp(-eta_p*AminusB(ixyz)**2.0d0)
         sum_i=0.0d0
         do i1=0,nlm(ia,ixyz)/2
           j1=nlm(ia,ixyz)-2*i1
           do i2=0,nlm(ib,ixyz)/2
             j2=nlm(ib,ixyz)-2*i2
             j=j1+j2
             facij=fact(i1)*fact(j1)*fact(i2)*fact(j2)*expp(ia)**(nlm(ia,ixyz)-i1)*expp(ib)**(nlm(ib,ixyz)-i2)
             sum_r=0.0d0
             if(abs(AminusB(ixyz)).gt.TOL) then !! avoid 0^0 !!
               do ir=0,j/2
                 xfac=eta_p**(j-ir)*(2.0d0*AminusB(ixyz))**(j-2*ir)/(fact(ir)*fact(j-2*ir))
                 if(MOD(ir,2).ne.0) xfac=-xfac
                 sum_r=sum_r+xfac
               end do
             else if(mod(j,2).eq.0) then
               sum_r=eta_p**(j/2)/fact(j/2)
               if(mod(j/2,2).ne.0) sum_r=-sum_r
             end if
             if(MOD(j1,2).ne.0) sum_r=-sum_r
             sum_i=sum_i+sum_r*fact(j)/facij
           end do
         end do
         do_ov=do_ov*sum_i
       end do
       sp(ia,ib)=do_ov
111   if(ia.ne.ib) sp(ib,ia)=sp(ia,ib)
     end do
   end do

!! contract primitive overlaps into the basis-function overlap matrix !!
   do i=1,nbasis
     do j=1,i
       S(i,j)=0.0d0
       k=1
       do while(nprimbas(k,i).ne.0)
         l=1
         do while(nprimbas(l,j).ne.0)
           S(i,j)=S(i,j)+coefpb(nprimbas(k,i),i)*coefpb(nprimbas(l,j),j)*sp(nprimbas(k,i),nprimbas(l,j))
           l=l+1
         end do
         k=k+1
       end do
       if(i.ne.j) s(j,i)=s(i,j)
     end do
   end do
   deallocate(sp)
   end subroutine do_overlap

   END MODULE basis_set 

!! ***** !!

   MODULE ao_matrices
   DOUBLE PRECISION, ALLOCATABLE :: c(:,:), p(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: cb(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: ps(:,:), pa(:,:), pb(:,:)
   DOUBLE PRECISION, ALLOCATABLE :: c_no(:,:), occ_no(:,:)
! c(igr,igr) -> alpha MO coefficients
! p(igr,igr) -> density matrix
! cb(igr,igr) -> beta MO coefficients
! pa(igr,igr) -> alpha density matrix
! pb(igr,igr) -> beta density matrix
! ps(igr,igr) -> spin density matrix (pa-pb)
! c_no(igr,igr) -> natural orbital coefficients
! occ_no(igr,igr) -> natural orbital occupation values
! igr is the number of basis functions

   CONTAINS

!! ********************************************************************* !!
!! subroutine: build_ao_matrices                                         !!
!! purpose: allocates the ao_matrices module's MO-coefficient/density-   !!
!! matrix arrays (c/p/cb/ps/pa/pb/c_no/occ_no) to the actual basis size. !!
!! called from input2.f's input(), right after igr becomes known --      !!
!! same pattern effao_mod/nao_mod's own allocate_* subroutines          !!
!! follow, just predating that naming convention.                        !!
!! arguments: igr (in) -- number of basis functions                      !!
!! author: MMO                                                           !!
!! ********************************************************************* !!
   SUBROUTINE build_ao_matrices(igr)

   ALLOCATE(c(igr,igr),p(igr,igr))
   ALLOCATE(cb(igr,igr))
   ALLOCATE(ps(igr,igr),pa(igr,igr),pb(igr,igr))
   ALLOCATE(c_no(igr,igr),occ_no(igr,igr))

   END SUBROUTINE build_ao_matrices

   END MODULE ao_matrices

!! ***** !!

   MODULE integration_grid
   INTEGER :: Nrad,Nang
   DOUBLE PRECISION :: pha, phb, rr00
   DIMENSION :: leved(32)
   DATA leved/6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974,&
   1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810/
   DOUBLE PRECISION,dimension(1000):: th,ph,w  
   DOUBLE PRECISION,dimension(500):: wr,xr  
!! one-electron grid given in # GRID (RADIAL/ANGULAR), 0 = default; set by !!
!! read_input.f's read_grid, used by build_integration_grid               !!
   INTEGER :: nrad_in=0, nang_in=0
! Nrad -> number of radial points in the atomic grid
! Nang -> number of angular points in the atomic grid
! pha -> rotated grid for zero-error (first roration)
! phb -> rotated grid for zero-error (second rotation)
! rr00 -> nuclear distance at which half of the radial points have been distributed
! th(nang) -> theta: angular coordinate for a given angular plane of points
! ph(nang) -> phi: angular coordinate for a given angular plane of points
! w(nang) -> Lebedev weight of each angular point
! wr(nrad) -> mathematical weight of a given radial surface
! xr(nrad) -> distance from grid center to a given radial surface

   CONTAINS

!! ********************************************************************* !!
!! subroutine: build_integration_grid                                    !!
!! purpose: picks the atom-centered one-electron grid (radial/angular    !!
!! points, rr00): the ENPART/POLAR high-accuracy defaults or the        !!
!! plain ones, each overridden by # GRID's RADIAL/ANGULAR (nrad_in,      !!
!! nang_in); prints the choice, then builds the grid via quad().         !!
!! arguments:                                                            !!
!!   ienpart, ipolar (in) -- select the high-accuracy defaults when any  !!
!!     is set                                                            !!
!! author: MMO, MGimf                                                    !!
!! ********************************************************************* !!
   SUBROUTINE build_integration_grid(ienpart, ipolar)
   use output_mod, only: print_value
   character(len=12) :: cr0
   common /nat/ nat,igr,ifg,nocc,nalf,nb,kop

!! ENPART/POLAR defaults for high-accuracy one-el integrations !!
   if(ienpart.eq.1.or.ipolar.eq.1) then
     nrad=150
     nang=590
     rr00=0.500d0

!! APOST legacy defaults for one-el integrations !!
   else
     nrad=40
     nang=146
     rr00=0.5d0
   end if

!! # GRID's RADIAL/ANGULAR (already checked by read_grid) !!
   if(nrad_in.gt.0) nrad=nrad_in
   if(nang_in.gt.0) nang=nang_in

!! Rotation angles... they are zero for one-el part, just to be consistent !!
!! written as literals: a module subroutine can't include parameter.h  !!
!! (see do_overlap), so ZERO is not defined here                         !!
   pha=0.0d0
   phb=0.0d0

   call print_box('SETTING ATOMIC GRIDS FOR INTEGRATION')
   call print_value('Radial points',nrad)
   call print_value('Angular points',nang)
   write(cr0,'(f7.3)') rr00
   call print_value('r0',trim(adjustl(cr0)),'bohr')
   call print_value('Grid points',nrad*nang*nat)

   call quad(nrad,nang)

   END SUBROUTINE build_integration_grid

!! ---- !!

!! ********************************************************************* !!
!! subroutine: check_grid                                                !!
!! purpose: stops the run unless the grid is one quad() can build: 1-500 !!
!!   radial points (wr/xr hold 500) and a Lebedev-Laikov angular size    !!
!!   (leved(1:19); quad() builds nothing for any other value).           !!
!! arguments:                                                            !!
!!   nrad0, nang0 (in) -- radial and angular points asked for            !!
!!   where        (in) -- where they were given, for the message         !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
   SUBROUTINE check_grid(nrad0,nang0,where)
   integer, intent(in) :: nrad0,nang0
   character(len=*), intent(in) :: where

   if(nrad0.lt.1.or.nrad0.gt.500) then
     write(*,'(2x,a,a,i0,a)') where,': RADIAL ',nrad0, &
       ' must be between 1 and 500'
     call apost_stop(' Invalid radial grid')
   end if
   if(all(leved(1:19).ne.nang0)) then
     write(*,'(2x,a,a,i0,a)') where,': ANGULAR ',nang0, &
       ' is not a Lebedev-Laikov grid; use one of'
     write(*,'(2x,19(1x,i0))') leved(1:19)
     call apost_stop(' Invalid angular grid')
   end if

   END SUBROUTINE check_grid

   END MODULE integration_grid

!! ***** !!

!! ********************************************************************* !!
!! module: effao_mod                                                     !!
!! purpose: replaces the legacy 'common /effao/' block (nmax x nmax,     !!
!!   fixed at compile time regardless of the real system size). Holds    !!
!!   the effective-atomic-orbital population matrices used by effao.f,   !!
!!   print.f and ueos.f. Allocated to the actual (igr,nat) once those    !!
!!   are known from the .fchk, so there is no silent size cap anymore.   !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
   MODULE effao_mod
   real*8, allocatable :: p0(:,:)    !! p0(igr,igr)    -- effective AO population matrix
   real*8, allocatable :: p0net(:,:) !! p0net(igr,nat) -- net population, per basis fn/atom
   real*8, allocatable :: p0gro(:,:) !! p0gro(igr,nat) -- gross population, per basis fn/atom
   real*8, allocatable :: p0fit(:,:) !! p0fit(igr,nat) -- fit (%) of each reprojected EFO (real space)
   integer, allocatable :: ip0(:)    !! ip0(nat)       -- per-atom effective AO bookkeeping
   real*8, allocatable :: p0coef(:,:,:)      !! p0coef(igr,igr,nat) -- per-fragment EFO coefficients (EOS real-space path),
                                             !! persists across ueffao3d_frag's per-fragment loop, unlike scratch p0
   real*8, allocatable :: p0poolcoef(:,:,:)  !! p0poolcoef(igr,igr,2) -- pooled+sorted+igr-truncated EFO coefficients per
                                             !! spin channel, built by eos_analysis, feeds the shared .fchk EFO writer
                                             !! (print.f's rwf_effao_orbprint/uwf_effao_orbprint)

   CONTAINS

   !! ********************************************************************* !!
   !! subroutine: allocate_effao                                            !!
   !! purpose: allocate the effao_mod arrays to the real system size.       !!
   !!   call once, right after igr/nat become known -- see input2.f, next   !!
   !!   to the existing call to build_ao_matrices(igr), which follows the   !!
   !!   same pattern for the ao_matrices module above.                      !!
   !! arguments:                                                            !!
   !!   igr (in) -- number of basis functions                               !!
   !!   nat (in) -- number of atoms                                         !!
   !! author: MGimf                                                         !!
   !! ********************************************************************* !!
   SUBROUTINE allocate_effao(igr,nat)
   integer, intent(in) :: igr,nat

   allocate(p0(igr,igr),p0net(igr,nat),p0gro(igr,nat),ip0(nat))
   allocate(p0fit(igr,nat))
   allocate(p0coef(igr,igr,nat),p0poolcoef(igr,igr,2))

   END SUBROUTINE allocate_effao

   END MODULE effao_mod

!! ***** !!

!! ********************************************************************* !!
!! module: nao_mod                                                       !!
!! purpose: replaces the legacy 'common /nao/' block (same nmax x nmax   !!
!!   cap as effao_mod above). Holds the natural-atomic-orbital matrices  !!
!!   used by effao.f and mulliken.f.                                     !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
   MODULE nao_mod
   real*8, allocatable :: unao(:,:)  !! unao(igr,igr)  -- natural AO transformation matrix
   real*8, allocatable :: ssnao(:,:) !! ssnao(igr,igr) -- overlap matrix in the NAO basis

   CONTAINS

!! ********************************************************************* !!
!! subroutine: allocate_nao                                              !!
!! purpose: allocate the nao_mod arrays to the real system size. call    !!
!!   once, right after igr becomes known -- see allocate_effao above.    !!
!! arguments:                                                            !!
!!   igr (in) -- number of basis functions                               !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
   SUBROUTINE allocate_nao(igr)
   integer, intent(in) :: igr

   allocate(unao(igr,igr),ssnao(igr,igr))

   END SUBROUTINE allocate_nao

   END MODULE nao_mod

!! ***** !!

!! ----------------------------------------------------------------------
!! module: version_mod
!! The release identity, in one place: printed in the banner (kiir) and
!! written in the results file (json_mod). Set version_date (and its ISO
!! form) to the release date when a version is released.
!! author: MGimf
!! ----------------------------------------------------------------------
   MODULE version_mod
   implicit none
   character(len=*), parameter :: version_number='5'
   character(len=*), parameter :: version_full='5.0'
   character(len=*), parameter :: version_date='22 September 2026'
   character(len=*), parameter :: version_date_iso='2026-09-22'
   END MODULE version_mod

!! ***** !!

!! ----------------------------------------------------------------------
!! module: elements_mod
!! Element symbols for the atom labels of the output, right-justified in
!! two characters (' H', 'He'): element_symbol(z), and element_label(z)
!! with a blank on each side ('  H ', ' He '), the four-character form of
!! the matrix and table rows. Atomic numbers outside 1-118 give '??'.
!! author: MGimf
!! ----------------------------------------------------------------------
   MODULE elements_mod
   implicit none
   character(len=2), parameter :: symbols(118) = (/ &
      ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne', &
      'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca', &
      'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
      'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr', &
      'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
      'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd', &
      'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
      'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg', &
      'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
      'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
      'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
      'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' /)

   CONTAINS

!! ---- !!
   function element_symbol(iz) result(sym)
   integer, intent(in) :: iz
   character(len=2) :: sym
   sym='??'
   if(iz.ge.1.and.iz.le.118) sym=symbols(iz)
   end function element_symbol

!! ---- !!
   function element_label(iz) result(lab)
   integer, intent(in) :: iz
   character(len=4) :: lab
   lab=' '//element_symbol(iz)//' '
   end function element_label

   END MODULE elements_mod

!! ***** !!


!! ***** !!

!! ----------------------------------------------------------------------
!! module: json_mod
!! Writes the machine-readable results file <job>.apost.json while the
!! run proceeds (schema: docs/source/developer/json-output.md). The file
!! is written as <job>.apost.json.tmp and renamed when it is closed, at
!! the normal end of the run or by apost_stop, so the final name only
!! ever holds a complete file. Reals are written with 17 significant
!! digits (they read back as the same double), NaN or infinity as null.
!!   json_open(job,inp,fchk)       start the file (header written)
!!   json_begin_object(key)        open an object ('' inside an array)
!!   json_begin_array(key)         open an array
!!   json_end                      close the innermost object or array
!!   json_value(key,value)         real*8, integer, logical or text
!!   json_null(key)                a null value
!!   json_vector(key,v,n)          real*8, integer or text array
!!   json_matrix(key,a,n,m,lda)    real*8 n x m matrix, as a list of rows
!!   json_vector_or_null(key,v,lok,n)  real*8 or integer array, null
!!                                 where lok is false
!!   json_column(key,title,unit,decimals,lsigned)  one column description
!!                                 of a table (inside its columns array)
!!   json_index(key,kind,title,unit,decimals,spin,lsigned,role)
!!                                 adds the value just written under key
!!                                 (in the current object) to the index
!!   json_close(status,message)    "normal" or "stopped"; writes the index
!!                                 and the run block, closes everything
!! author: MGimf
!! ----------------------------------------------------------------------
   MODULE json_mod
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: json_open, json_close, json_begin_object, json_begin_array, &
             json_end, json_value, json_null, json_vector, json_matrix, &
             json_index, json_is_open, json_vector_or_null, json_column, &
             json_run_id, json_atomic_weights

   integer, parameter :: maxdepth=32
   integer :: iu=-1, depth=0
   logical :: is_open=.false.
   logical :: first(0:maxdepth)=.true.
   character(len=1) :: closer(0:maxdepth)=' '
!! key of each open container: the JSON Pointer of the current position !!
   character(len=64) :: keys(0:maxdepth)=' '
   character(len=:), allocatable :: final_name, temp_name
   character(len=:), allocatable :: job_name, inp_name, fchk_name
   character(len=25) :: run_date
!! run id: the start time and 4 hex digits of the system clock, written !!
!! to run.id and to the title line of the orbital .fchk files            !!
   character(len=30) :: run_id=' '

!! index entries, written by json_close (fields as in the index)       !!
   type index_entry
     character(len=128) :: id,analysis,path,title
     character(len=16) :: kind,unit,spin,role
     integer :: decimals
     logical :: signed
   end type index_entry
   type(index_entry), allocatable :: entries(:)
   integer :: nentries=0

!! cell functions of the Becke-type atomic weights (json_atomic_weights), !!
!! written by json_close; nweights=0: no such atomic definition          !!
   character(len=8) :: weights_profile=' '
   integer :: weights_stiffness=0, nweights=0
   real*8, allocatable :: weights_chi(:,:)

   interface json_value
     module procedure value_real, value_int, value_logical, value_text
   end interface
   interface json_vector
     module procedure vector_real, vector_int, vector_text
   end interface
   interface json_vector_or_null
     module procedure vector_real_or_null, vector_int_or_null
   end interface

   CONTAINS

!! ---- !!
!! the run id (blank before json_open)                                  !!
   function json_run_id() result(id)
   character(len=:), allocatable :: id
   id=trim(run_id)
   end function json_run_id

!! ---- !!
   function to_lower(s) result(r)
   character(len=*), intent(in) :: s
   character(len=len(s)) :: r
   integer :: i
   r=s
   do i=1,len(s)
     if(s(i:i).ge.'A'.and.s(i:i).le.'Z') r(i:i)=achar(iachar(s(i:i))+32)
   end do
   end function to_lower

!! ---- !!
   logical function json_is_open()
   json_is_open=is_open
   end function json_is_open

!! ---- !!
!! start the file: remove any results file of an earlier run, then write !!
!! the header into the temporary file                                     !!
   subroutine json_open(job,inp,fchk)
   use version_mod, only: version_full, version_date_iso
   character(len=*), intent(in) :: job,inp,fchk
   integer :: ios,iold
   integer(kind=8) :: iclock
   character(len=4) :: hex4
   character(len=8) :: d
   character(len=10) :: t
   character(len=5) :: z
!! the file goes next to the .fchk (the job name may include a folder); !!
!! the names recorded in it are bare file names                         !!
   final_name=trim(job)//'.apost.json'
   temp_name=final_name//'.tmp'
   job_name=bare(job)
   inp_name=bare(inp)
   fchk_name=bare(fchk)
   call date_and_time(d,t,z)
   run_date=d(1:4)//'-'//d(5:6)//'-'//d(7:8)//'T'//t(1:2)//':'//t(3:4)// &
            ':'//t(5:6)//z(1:3)//':'//z(4:5)
   call system_clock(iclock)
   write(hex4,'(z4.4)') mod(iclock,65536)
   run_id=trim(run_date)//'/'//to_lower(hex4)
   open(newunit=iold,file=final_name,status='old',iostat=ios)
   if(ios.eq.0) close(iold,status='delete')
   open(newunit=iu,file=temp_name,status='replace',action='write', &
        iostat=ios)
   if(ios.ne.0) return
   is_open=.true.
   depth=0
   write(iu,'(a)',advance='no') '{'
   depth=1
   first(1)=.true.
   closer(1)='}'
   call json_value('schema','apost3d-results')
   call json_value('schema_version','0.16')
   call json_begin_object('program')
   call json_value('name','APOST-3D')
   call json_value('version',version_full)
   call json_value('version_date',version_date_iso)
   call json_end()
   end subroutine json_open

!! ---- !!
!! keep the cell functions of the Becke-type atomic weights (known only  !!
!! after the analyses block is open): json_close writes them as the      !!
!! top-level atomic_weights block                                        !!
   subroutine json_atomic_weights(profile,k,chi,n)
   character(len=*), intent(in) :: profile
   integer, intent(in) :: k,n
   real*8, intent(in) :: chi(n,n)
   weights_profile=profile
   weights_stiffness=k
   nweights=n
   if(allocated(weights_chi)) deallocate(weights_chi)
   allocate(weights_chi(n,n))
   weights_chi=chi
   end subroutine json_atomic_weights

!! ---- !!
!! close every open container, write the atomic weights, the index and  !!
!! the run block, rename the file                                       !!
   subroutine json_close(status,message)
   use omp_lib, only: omp_get_max_threads
   character(len=*), intent(in) :: status,message
   if(.not.is_open) return
   do while(depth.gt.1)
     call json_end()
   end do
   if(nweights.gt.0) then
     call json_begin_object('atomic_weights')
     call json_value('profile',trim(weights_profile))
     call json_value('stiffness',weights_stiffness)
     call json_matrix('radius_ratios',weights_chi,nweights,nweights, &
                      nweights)
     call json_end()
   end if
   call write_index()
   call json_begin_object('run')
   call json_value('job_name',job_name)
   call json_value('input_file',inp_name)
   call json_value('wavefunction_file',fchk_name)
   call json_value('date',trim(run_date))
   call json_value('id',trim(run_id))
   call json_value('threads',omp_get_max_threads())
   call json_value('status',trim(status))
   if(len_trim(message).gt.0) then
     call json_value('stop_message',trim(adjustl(message)))
   else
     call json_null('stop_message')
   end if
   call json_end()
   write(iu,'(a)') ''
   write(iu,'(a)') '}'
   close(iu)
   is_open=.false.
   call rename(temp_name,final_name)
   end subroutine json_close

!! ---- !!
!! comma, line break and indentation before a new item; the key if the  !!
!! current container is an object                                     !!
   subroutine item(key)
   character(len=*), intent(in) :: key
   if(.not.first(depth)) write(iu,'(a)',advance='no') ','
   first(depth)=.false.
   write(iu,'(a)',advance='no') new_line('a')//repeat('  ',depth)
   if(closer(depth).eq.'}') write(iu,'(a)',advance='no') quoted(key)//': '
   end subroutine item

!! ---- !!
   subroutine json_begin_object(key)
   character(len=*), intent(in) :: key
   if(.not.is_open) return
   call item(key)
   write(iu,'(a)',advance='no') '{'
   depth=depth+1
   first(depth)=.true.
   closer(depth)='}'
   keys(depth)=key
   end subroutine json_begin_object

!! ---- !!
   subroutine json_begin_array(key)
   character(len=*), intent(in) :: key
   if(.not.is_open) return
   call item(key)
   write(iu,'(a)',advance='no') '['
   depth=depth+1
   first(depth)=.true.
   closer(depth)=']'
   keys(depth)=key
   end subroutine json_begin_array

!! ---- !!
   subroutine json_end()
   if(.not.is_open.or.depth.le.1) return
   if(.not.first(depth)) write(iu,'(a)',advance='no') &
     new_line('a')//repeat('  ',depth-1)
   write(iu,'(a)',advance='no') closer(depth)
   depth=depth-1
   end subroutine json_end

!! ---- !!
   subroutine value_real(key,x)
   character(len=*), intent(in) :: key
   real*8, intent(in) :: x
   if(.not.is_open) return
   call item(key)
   write(iu,'(a)',advance='no') real_text(x)
   end subroutine value_real

!! ---- !!
   subroutine value_int(key,n)
   character(len=*), intent(in) :: key
   integer, intent(in) :: n
   character(len=24) :: buf
   if(.not.is_open) return
   call item(key)
   write(buf,'(i0)') n
   write(iu,'(a)',advance='no') trim(buf)
   end subroutine value_int

!! ---- !!
   subroutine value_logical(key,l)
   character(len=*), intent(in) :: key
   logical, intent(in) :: l
   if(.not.is_open) return
   call item(key)
   if(l) then
     write(iu,'(a)',advance='no') 'true'
   else
     write(iu,'(a)',advance='no') 'false'
   end if
   end subroutine value_logical

!! ---- !!
   subroutine value_text(key,text)
   character(len=*), intent(in) :: key,text
   if(.not.is_open) return
   call item(key)
   write(iu,'(a)',advance='no') quoted(text)
   end subroutine value_text

!! ---- !!
   subroutine json_null(key)
   character(len=*), intent(in) :: key
   if(.not.is_open) return
   call item(key)
   write(iu,'(a)',advance='no') 'null'
   end subroutine json_null

!! ---- !!
   subroutine vector_real(key,v,n)
   character(len=*), intent(in) :: key
   integer, intent(in) :: n
   real*8, intent(in) :: v(n)
   integer :: i
   if(.not.is_open) return
   call item(key)
   write(iu,'(a)',advance='no') '['
   do i=1,n
     if(i.gt.1) write(iu,'(a)',advance='no') ', '
     write(iu,'(a)',advance='no') real_text(v(i))
   end do
   write(iu,'(a)',advance='no') ']'
   end subroutine vector_real

!! ---- !!
   subroutine vector_int(key,v,n)
   character(len=*), intent(in) :: key
   integer, intent(in) :: n
   integer, intent(in) :: v(n)
   integer :: i
   character(len=24) :: buf
   if(.not.is_open) return
   call item(key)
   write(iu,'(a)',advance='no') '['
   do i=1,n
     if(i.gt.1) write(iu,'(a)',advance='no') ', '
     write(buf,'(i0)') v(i)
     write(iu,'(a)',advance='no') trim(buf)
   end do
   write(iu,'(a)',advance='no') ']'
   end subroutine vector_int

!! ---- !!
   subroutine vector_text(key,v,n)
   character(len=*), intent(in) :: key
   integer, intent(in) :: n
   character(len=*), intent(in) :: v(n)
   integer :: i
   if(.not.is_open) return
   call item(key)
   write(iu,'(a)',advance='no') '['
   do i=1,n
     if(i.gt.1) write(iu,'(a)',advance='no') ', '
     write(iu,'(a)',advance='no') quoted(trim(adjustl(v(i))))
   end do
   write(iu,'(a)',advance='no') ']'
   end subroutine vector_text

!! ---- !!
   subroutine vector_real_or_null(key,v,lok,n)
   character(len=*), intent(in) :: key
   integer, intent(in) :: n
   real*8, intent(in) :: v(n)
   logical, intent(in) :: lok(n)
   integer :: i
   if(.not.is_open) return
   call item(key)
   write(iu,'(a)',advance='no') '['
   do i=1,n
     if(i.gt.1) write(iu,'(a)',advance='no') ', '
     if(lok(i)) then
       write(iu,'(a)',advance='no') real_text(v(i))
     else
       write(iu,'(a)',advance='no') 'null'
     end if
   end do
   write(iu,'(a)',advance='no') ']'
   end subroutine vector_real_or_null

!! ---- !!
   subroutine vector_int_or_null(key,v,lok,n)
   character(len=*), intent(in) :: key
   integer, intent(in) :: n
   integer, intent(in) :: v(n)
   logical, intent(in) :: lok(n)
   integer :: i
   character(len=24) :: buf
   if(.not.is_open) return
   call item(key)
   write(iu,'(a)',advance='no') '['
   do i=1,n
     if(i.gt.1) write(iu,'(a)',advance='no') ', '
     if(lok(i)) then
       write(buf,'(i0)') v(i)
       write(iu,'(a)',advance='no') trim(buf)
     else
       write(iu,'(a)',advance='no') 'null'
     end if
   end do
   write(iu,'(a)',advance='no') ']'
   end subroutine vector_int_or_null

!! ---- !!
!! a table column: {"key", "title", "unit" (null for ' '), "decimals",  !!
!! "signed"}                                                            !!
   subroutine json_column(key,title,unit,decimals,lsigned)
   character(len=*), intent(in) :: key,title,unit
   integer, intent(in) :: decimals
   logical, intent(in) :: lsigned
   if(.not.is_open) return
   call json_begin_object('')
   call json_value('key',key)
   call json_value('title',title)
   call text_or_null('unit',unit)
   call json_value('decimals',decimals)
   call json_value('signed',lsigned)
   call json_end()
   end subroutine json_column

!! ---- !!
!! rows 1..n, columns 1..m of a(lda,*)                                   !!
   subroutine json_matrix(key,a,n,m,lda)
   character(len=*), intent(in) :: key
   integer, intent(in) :: n,m,lda
   real*8, intent(in) :: a(lda,*)
   integer :: i
   if(.not.is_open) return
   call json_begin_array(key)
   do i=1,n
     call vector_real('',a(i,1:m),m)
   end do
   call json_end()
   end subroutine json_matrix

!! ---- !!
!! index entry for the value written under key in the current object:  !!
!! path = JSON Pointer to it, id = the path below /analyses joined with !!
!! dots, analysis = the analysis block it is in. kind: scalar,          !!
!! atom_vector, atom_matrix, fragment_vector or fragment_matrix. unit   !!
!! and spin: ' ' is written as null; decimals < 0 too (a table, whose  !!
!! columns carry their own). lsigned: the values can be                 !!
!! negative by nature (a diverging colour scale suits them). role: ' '  !!
!! for a result, 'comparison' for a value printed only to compare with  !!
!! (the key is left out for ' ').                                       !!
   subroutine json_index(key,kind,title,unit,decimals,spin,lsigned,role)
   character(len=*), intent(in) :: key,kind,title,unit,spin,role
   integer, intent(in) :: decimals
   logical, intent(in) :: lsigned
   type(index_entry), allocatable :: tmp(:)
   character(len=128) :: path,id
   integer :: i
   if(.not.is_open) return
   path=' '
   id=' '
   do i=2,depth
     path=trim(path)//'/'//trim(keys(i))
     if(i.ge.3) id=trim(id)//trim(keys(i))//'.'
   end do
   path=trim(path)//'/'//key
   id=trim(id)//key
   if(.not.allocated(entries)) allocate(entries(32))
   if(nentries.eq.size(entries)) then
     allocate(tmp(2*nentries))
     tmp(1:nentries)=entries(1:nentries)
     call move_alloc(tmp,entries)
   end if
   nentries=nentries+1
   entries(nentries)%id=id
   entries(nentries)%analysis=' '
   if(depth.ge.3) entries(nentries)%analysis=keys(3)
   entries(nentries)%path=path
   entries(nentries)%kind=kind
   entries(nentries)%title=title
   entries(nentries)%unit=unit
   entries(nentries)%decimals=decimals
   entries(nentries)%spin=spin
   entries(nentries)%signed=lsigned
   entries(nentries)%role=role
   end subroutine json_index

!! ---- !!
!! the index: one object per entry, in the order they were added       !!
   subroutine write_index()
   integer :: i
   call json_begin_array('index')
   do i=1,nentries
     call json_begin_object('')
     call json_value('id',trim(entries(i)%id))
     call json_value('analysis',trim(entries(i)%analysis))
     call text_or_null('spin',entries(i)%spin)
     if(len_trim(entries(i)%role).gt.0) &
       call json_value('role',trim(entries(i)%role))
     call json_value('path',trim(entries(i)%path))
     call json_value('kind',trim(entries(i)%kind))
     call json_value('title',trim(entries(i)%title))
     call text_or_null('unit',entries(i)%unit)
     if(entries(i)%decimals.ge.0) then
       call json_value('decimals',entries(i)%decimals)
     else
       call json_null('decimals')
     end if
     call json_value('signed',entries(i)%signed)
     call json_end()
   end do
   call json_end()
   end subroutine write_index

!! ---- !!
   subroutine text_or_null(key,text)
   character(len=*), intent(in) :: key,text
   if(len_trim(text).eq.0) then
     call json_null(key)
   else
     call json_value(key,trim(text))
   end if
   end subroutine text_or_null

!! ---- !!
!! a path without its folders                                           !!
   function bare(path) result(b)
   character(len=*), intent(in) :: path
   character(len=:), allocatable :: b
   b=trim(path)
   b=b(index(b,'/',back=.true.)+1:)
   end function bare

!! ---- !!
   function real_text(x) result(t)
   real*8, intent(in) :: x
   character(len=:), allocatable :: t
   character(len=32) :: buf
   if(.not.ieee_is_finite(x)) then
     t='null'
   else
     write(buf,'(es24.16e3)') x
     t=trim(adjustl(buf))
   end if
   end function real_text

!! ---- !!
!! a JSON string: quotes, backslashes and control characters escaped    !!
   function quoted(s) result(q)
   character(len=*), intent(in) :: s
   character(len=:), allocatable :: q
   integer :: i,ic
   character(len=6) :: hex
   q='"'
   do i=1,len(s)
     ic=iachar(s(i:i))
     if(s(i:i).eq.'"') then
       q=q//'\"'
     else if(s(i:i).eq.'\') then
       q=q//'\\'
     else if(ic.eq.10) then
       q=q//'\n'
     else if(ic.eq.9) then
       q=q//'\t'
     else if(ic.lt.32) then
       write(hex,'(a,z4.4)') '\u',ic
       q=q//hex
     else
       q=q//s(i:i)
     end if
   end do
   q=q//'"'
   end function quoted

   END MODULE json_mod

!! ***** !!

!! ********************************************************************* !!
!! module: timing_mod                                                    !!
!! purpose: grep-friendly CPU + wall-clock timers                        !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
   MODULE timing_mod
   CONTAINS

!! ********************************************************************* !!
!! subroutine: get_wall_time                                             !!
!! purpose: wall-clock reading via system_clock (seconds). integer*8     !!
!!   counters avoid the default kind's sub-hour wraparound.              !!
!! arguments:                                                            !!
!!   twall (out) -- current wall-clock reading, seconds                  !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
   SUBROUTINE get_wall_time(twall)
   IMPLICIT REAL*8(A-H,O-Z)
   integer*8 :: icount,icount_rate
   real*8, intent(out) :: twall

   call system_clock(count=icount,count_rate=icount_rate)
   twall=real(icount,8)/real(icount_rate,8)

   END SUBROUTINE get_wall_time

!! ***** !!

!! ********************************************************************* !!
!! subroutine: print_timer                                               !!
!! purpose: fixed-format "TIMING CPU/WALL :: label value" print, so      !!
!!   `grep "TIMING"` finds every timer regardless of caller.             !!
!! arguments:                                                            !!
!!   label (in) -- section name (<=40 chars, keeps the value column      !!
!!                 fixed across call sites)                              !!
!!   tcpu  (in) -- cpu_time() delta, seconds (thread-summed)             !!
!!   twall (in) -- get_wall_time() delta, seconds (true wall-clock)      !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
   SUBROUTINE print_timer(label,tcpu,twall)
   IMPLICIT REAL*8(A-H,O-Z)
   character(len=*), intent(in) :: label
   real*8, intent(in) :: tcpu,twall
   character(len=40) :: label40

!! assignment left-justifies+pads (unlike A40 on write, which right-justifies) !!
   label40=label

   write(*,'(2x,a,a40,f14.2,a2)') 'TIMING CPU  :: ',label40,tcpu,' s'
   write(*,'(2x,a,a40,f14.2,a2)') 'TIMING WALL :: ',label40,twall,' s'

   END SUBROUTINE print_timer

   END MODULE timing_mod

!! ***** !!

!! ********************************************************************* !!
!! module: input_options_mod                                             !!
!! purpose: holds the .inp keyword flags parsed by read_input() (see     !!
!! read_input.f) that main.f's own control flow (validation, iopt(200)   !!
!! population, dispatch) reads directly. Flags already carried by       !!
!! an existing COMMON block (icas, ibcp, aerf, iaccur, nrad22, etc.)     !!
!! are NOT duplicated here -- see main.f's own COMMON declarations.      !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
   MODULE input_options_mod
   IMPLICIT REAL*8(A-H,O-Z)

!! wavefunction source / density choice !!
   integer :: iallpo,ndens0
   character(len=80) :: dmfile1,dmfile2   !! # DM file names, for INPUT SUMMARY

!! atoms in molecules (Hilbert-space + real-space AIM selection) !!
   integer :: imulli,ihirsh,inewbec,istiff,iradmat,itfvc
   integer :: ibecke,itfvcdef

!! QTAIM !!
   integer :: iqtaim,istep,inna,imaxdist,iscreening,ipath

!! miscellaneous / integration control !!
   integer :: iopop,idoint,ilaplacian,ielcount, &
              iatdens,inopop
!! VERBOSE (# METHOD): 0 (default) prints the final results only, 1 also  !!
!! the intermediate matrices (for now only ENPART uses it)                 !!
   integer :: iverbose
   real*8  :: Rmax

!! eff-AOs and EOS !!
   integer :: ieffao,ieffthr,icube,jcubthr,kcubthr,ieos,iueos,ieoscent, &
              iloba
   real*8  :: xthresh
   real*8  :: cubespacing,cubepadding !! # CUBE grid spacing / box margin beyond the atoms (bohr) !!
!! # CUBE NEG_EFOS: cubes of GEOS paired EFOs with net occupation         !!
!! <= -inegcubthr/1000 (flag + optional threshold, default 25)            !!
   integer :: inegefos,inegcubthr

!! local spin and correlated-WF input !!
   integer :: ispin,icorr

!! ENPART -- igrid2e: 1 if # GRID gives RADIAL_2E or ANGULAR_2E; irot2e:  !!
!! where the two-electron rotation angles come from (0 calibrated table,  !!
!! 1 ROTATION_2E, 2 not calibrated for ANGULAR_2E). Set by read_grid.      !!
!! cfuncname: the predefined functional keyword of # ENPART (SVWN,      !!
!! B3LYP, ...), blank when the libxc ids were given with LIBRARY       !!
   character(len=12) :: cfuncname=' '
   integer :: ienpart,ihf,id_xcfunc,id_xfunc,id_cfunc,iecorr, &
              ithrebod,iexact,ihomo,idek,iionic,ianalytical, &
              igrid2e,irot2e
!! TAU_BODEN (# ENPART, exploration, meta-GGA only): the pair tau of the  !!
!! BODEN term -- 1 kernel, 2 von Weizsacker, 3 2*wA*wB*tau (tauboden)     !!
   integer :: itauboden

!! DFT-DM1 (approximate one-particle RDM1) !!
   integer :: idftdm1,id_func_dm1,inatorb_dm1

!! NLOP / static field !!
   integer :: ipolar,ifield

!! atom/fragment restriction !!
   integer :: idoat,idofr

!! correlated-WF DM1/DM2 input !!
   integer :: iorca,ipyscf

!! OSLO; branching iterations per channel: 1 = every channel, 2 = alpha, !!
!! 3 = beta (BRANCH_ITERATION, BRANCH_ALPHA, BRANCH_BETA)                 !!
   integer :: ioslo,ilow2,ifolitol,ioslofchk
   integer, parameter :: maxbranch=20
   integer :: nbranch(3),ibranchit(maxbranch,3)

!! external .fchk sources !!
   integer :: iqchem

!! X-ray scattering !!
   integer :: iscattfact

   END MODULE input_options_mod

!! ***** !!