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
!!                         p0gro/ip0), replaces legacy common /effao/     !!
!!                         CONTAINS allocate_effao                        !!
!!   nao_mod           -- natural-AO matrices (unao/ssnao), replaces      !!
!!                         legacy common /nao/                            !!
!!                         CONTAINS allocate_nao                          !!
!!   stv_mod           -- QTAIM overlap/kinetic matrices (sp/tt),         !!
!!                         replaces legacy common /stv/                   !!
!!                         CONTAINS allocate_stv                          !!
!!   timing_mod        -- CPU+wall-clock timers                           !!
!!                         CONTAINS get_wall_time, print_timer            !!
!!   input_options_mod -- ~70 .inp keyword flags parsed by read_input()   !!
!!                         data only, no CONTAINS                         !!
!!                                                                        !!
!! Worth a look next:                                                     !!
!!   - basis_set/ao_matrices/integration_grid's own top-of-module         !!
!!     variable-doc comments are still old-style single-`!`, not          !!
!!     converted to `!! !!`. Some are openly uncertain about what the     !!
!!     array holds: ao_matrices has "pa/pb(igr,igr) -> density matrix     !!
!!     only for alpha/beta spinorbitals?" and "occ_no(igr,igr) -> ?";     !!
!!     integration_grid has several "probably ..." guesses for th/ph/w.   !!
!!   - a stray `!!MMO- MODULE TESTING STARTS HERE.` marker sits right     !!
!!     before the ao_matrices module statement -- reads like leftover     !!
!!     development scaffolding rather than real documentation.            !!
!! ********************************************************************** !!

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
!MMO- iptob_cartesian(numprim) -> orbital to which each primitive belongs
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
   write(*,'(2x,a,1x,i0)') 'Number of atoms                  :', &
     natoms
   write(*,'(2x,a,1x,i0)') 'Number of basis functions        :', &
     nbasis
   write(*,'(2x,a,1x,i0)') 'Primitive gaussians              :', &
     numprim

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
       stop 'angular momentum not implemented'
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
           iptob_cartesian(numprim)=nbasis+k !MMO- ptob map for cartesian
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
! pa(igr,igr) -> density matrix only for alpha spinorbitals?
! pb(igr,igr) -> density matrix only for beta spinorbitals?
! ps(igr,igr) -> spin density matrix (pa-pb)
! c_no(igr,igr) -> c for Natural Orbitals?
! occ_no(igr,igr) -> natural orbital occupation values
! Note: igr is number is basis functions

   CONTAINS

!! ********************************************************************* !!
!! subroutine: build_ao_matrices                                         !!
!! purpose: allocates the ao_matrices module's MO-coefficient/density-   !!
!! matrix arrays (c/p/cb/ps/pa/pb/c_no/occ_no) to the actual basis size. !!
!! called from input2.f's input(), right after igr becomes known --      !!
!! same pattern effao_mod/nao_mod/stv_mod's own allocate_* subroutines   !!
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
! Nrad -> number of radial points in the atomic grid
! Nang -> number of angular points in the atomic grid
! pha -> rotated grid for zero-error (first roration)
! phb -> rotated grid for zero-error (second rotation)
! rr00 -> nuclear distance at which half of the radial points have been distributed
! th(nang) -> theta: angular coordinate for a given angular plane of points
! ph(nang) -> phi: angular coordinate for a given angular plane of points
! w(nang) -> probably mathematical weight of each angular plane?
! wr(nrad) -> mathematical weight of a given radial surface
! xr(nrad) -> distance from grid center to a given radial surface

   CONTAINS

!! ********************************************************************* !!
!! subroutine: build_integration_grid                                    !!
!! purpose: picks the atom-centered grid size (radial/angular points,    !!
!! rr00) -- from command-line overrides if given, else ENPART/POLAR/     !!
!! EDAIQA high-accuracy defaults, else the plain one-electron defaults;  !!
!! prints the choice, then builds the grid via quad().                   !!
!! arguments:                                                            !!
!!   ienpart, ipolar, ifinegrid (in) -- select the high-accuracy         !!
!!     defaults when any is set. The EDAIQA branch also tests iedaiqa,   !!
!!     which is never assigned anywhere (implicitly-typed local, not a   !!
!!     dummy arg/COMMON/module var) -- undefined behavior, no test       !!
!!     exercises # EDAIQA. Needs sign-off before fixing.                 !!
!! author: MMO, MGimf                                                    !!
!! ********************************************************************* !!
   SUBROUTINE build_integration_grid(ienpart, ipolar,ifinegrid)
   common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
   character*30 integ1,integ2

!! command-line override (argv(2)/argv(3)), rarely used in practice !!
   call getarg(2,integ1)
   call getarg(3,integ2)
   if(integ1.ne.' '.and.integ2.ne.' ') then
     read(integ1,'(i4)') Nrad
     read(integ2,'(i4)') Nang
     if(Nrad.gt.500) stop 'Max number of radial points  is 500 '
     do 111 i=1,18
       npoints=leved(i)
       if(nang.lt.leved(i+1)) goto 211
  111 continue
     npoints=leved(19)
  211 continue
     print *,' Angular points:',npoints
     nang=npoints
     rr00=0.500d0

!! ENPART/POLAR/EDAIQA defaults for high-accuracy one-el integrations !!
   else if(ienpart.eq.1.or.ipolar.eq.1.or.iedaiqa.eq.1) then
     nrad=150
     nang=590
     if(ifinegrid.eq.1) nang=974
     rr00=0.500d0

!! APOST legacy defaults for one-el integrations !!
   else
     nrad=40
     nang=146
     rr00=0.5d0
   end if

!! Rotation angles... they are zero for one-el part, just to be consistent !!
!! written as a literal, not the ZERO symbol -- this module-CONTAINS      !!
!! subroutine can't include parameter.h (see do_overlap above), and ZERO  !!
!! is not declared here, so it was silently an untyped uninitialized      !!
!! variable, not the constant 0.0d0 (found 2026-08-20).                   !!
   pha=0.0d0
   phb=0.0d0

   call print_box('SETTING ATOMIC GRIDS FOR INTEGRATION')
   write(*,'(2x,a,1x,i0)') 'Radial points  :',nrad
   write(*,'(2x,a,1x,i0)') 'Angular points :',nang
   write(*,'(2x,a,1x,f7.3)') 'r0 (radial)    :',rr00
   write(*,'(2x,a,1x,i0)') 'Grid points    :',nrad*nang*nat

   call quad(nrad,nang)

   END SUBROUTINE build_integration_grid

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
   integer, allocatable :: ip0(:)    !! ip0(nat)       -- per-atom effective AO bookkeeping
   real*8, allocatable :: p0coef(:,:,:)     !! p0coef(igr,igr,nat) -- per-fragment EFO coefficients (EOS real-space path),
                                             !! persists across ueffao3d_frag's per-fragment loop, unlike scratch p0
   real*8, allocatable :: p0poolcoef(:,:,:) !! p0poolcoef(igr,igr,2) -- pooled+sorted+igr-truncated EFO coefficients per
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

!! ********************************************************************* !!
!! module: stv_mod                                                       !!
!! purpose: replaces the legacy 'common /stv/' block (same nmax x nmax   !!
!!   cap as effao_mod above). holds the QTAIM overlap/kinetic-energy-    !!
!!   related matrices used by qtaim.f.                                   !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
   MODULE stv_mod
   real*8, allocatable :: sp(:,:) !! sp(igr,igr) -- QTAIM overlap-related matrix
   real*8, allocatable :: tt(:,:) !! tt(igr,igr) -- QTAIM kinetic-energy-related matrix

   CONTAINS

!! ********************************************************************* !!
!! subroutine: allocate_stv                                              !!
!! purpose: allocate the stv_mod arrays to the real system size. call    !!
!!   once, right after igr becomes known -- see allocate_effao above.    !!
!! arguments:                                                            !!
!!   igr (in) -- number of basis functions                               !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
   SUBROUTINE allocate_stv(igr)
   integer, intent(in) :: igr

   allocate(sp(igr,igr),tt(igr,igr))

   END SUBROUTINE allocate_stv

   END MODULE stv_mod

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
!! population, dispatch) reads directly -- replaces what used to be      !!
!! plain implicitly-typed locals in main.f. Flags already carried by     !!
!! an existing COMMON block (icas, ibcp, aerf, iaccur, nrad22, etc.)     !!
!! are NOT duplicated here -- see main.f's own COMMON declarations.      !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
   MODULE input_options_mod
   IMPLICIT REAL*8(A-H,O-Z)

!! wavefunction source / density choice !!
   integer :: iwfn,iallpo,ndens0

!! atoms in molecules (Hilbert-space + real-space AIM selection) !!
   integer :: imulli,ihirsh,inewbec,istiff,iradmat,itfvc

!! QTAIM !!
   integer :: iqtaim,istep,inna,imaxdist,iscreening,ipath

!! miscellaneous / integration control !!
   integer :: iopop,isha,idoint,ipca,ilaplacian,ifinegrid,ielcount, &
              iatdens,inopop
   real*8  :: Rmax

!! eff-AOs and EOS !!
   integer :: ieffao,ieffthr,icube,jcubthr,kcubthr,ieos,iueos,ieoscent, &
              iloba
   real*8  :: xthresh
   real*8  :: cubespacing,cuberadscale !! # CUBE grid spacing (bohr) / atomic-radius padding scale !!

!! local spin and correlated-WF input !!
   integer :: ispin,icorr,idafh

!! ENPART -- iigrid is EDAIQA's own MOD-GRIDTWOEL flag (read_input.f's    !!
!! inline "# EDAIQA" readchar); ENPART's own MOD-GRIDTWOEL request goes    !!
!! through read_gridtwoel("# ENPART", ienpart_gridtwoel) instead, its own !!
!! dedicated flag -- the two used to alias the same variable (whichever   !!
!! section parsed last in read_input.f won), a real bug fixed 2026-08-30. !!
   integer :: ienpart,ihf,id_xcfunc,id_xfunc,id_cfunc,iecorr, &
              ithrebod,iexact,ihomo,idek,iionic,ianalytical,itop,ietop, &
              ipairs,iigrid,ienpart_gridtwoel

!! DFT-DM1 (approximate one-particle RDM1, formerly HIRAO internally) !!
   integer :: idftdm1,id_func_dm1,inatorb_dm1

!! EDAIQA !!
   integer :: iedaiqa,iflip

!! NLOP / static field !!
   integer :: ipolar,ifield

!! atom/fragment restriction !!
   integer :: idoat,idofr

!! correlated-WF DM1/DM2 input !!
   integer :: iorca,ipyscf

!! OSLO !!
   integer :: ioslo,ilow2,ifolitol,ibranch,ioslofchk

!! external .fchk sources !!
   integer :: iqchem,imokit

!! X-ray scattering !!
   integer :: iscattfact

   END MODULE input_options_mod

!! ***** !!