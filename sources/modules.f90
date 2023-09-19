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

   SUBROUTINE build_basis()
   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
   PARAMETER(PI=4.0d0*DATAN(1.0d0),TOL=1.0d-8)
   integer, allocatable :: mnsh(:),iatsh(:),mssh(:)
   real*8, allocatable :: expsh(:),c1(:),c2(:),xnorm(:),coefp(:)
   integer :: ncshell,npshell
   integer :: dummy,ilog 
   DIMENSION mult(-5:5)
   DATA mult/11,9,7,5,4,1,3,6,10,15,21/

!   maxl=int_locate(15,"Highest angular",ilog)
!   maxk=int_locate(15,"Largest degree",ilog)
! reading basis set info
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
!   dummy=int_locate(15,"Coordinatomses of e",ilog)
!   read(15,*)(x(i),y(i),z(i),i=1,ncshell)

! end reading basis set info

! processing basis set
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
        print *,"Number of atoms:",natoms
        print *,"Number of basis functions:",nbasis
        print *,'Total number of primitive gaussians',numprim
       
        allocate(ihold(nbasis),llim(natoms),iulim(natoms))
! basis to atom map
        ii=0
        do i=1,ncshell
         do k=1,mult(mssh(i))
          ii=ii+1
          ihold(ii)=iatsh(i)
         end do
        end do
! setting basis set limits for mulliken
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

! angular momentum of primitives and primtive to atom map 
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

!        do i=1,numprim
!          write(*,'(5i3)') i,(nlm(i,k),k=1,3),iptoat(i)
!        end do

!  list primitive exponents and coefficients
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

!primitive normalization
        do i=1,numprim
          nn=nlm(i,1)
          ll=nlm(i,2)
          mm=nlm(i,3)
          fnn=fact(nn)/fact(2*nn)
          fll=fact(ll)/fact(2*ll)
          fmm=fact(mm)/fact(2*mm)
          xnorm(i)=(2.0d0*expp(i)/PI)**0.75d0*DSQRT((8.0d0*expp(i))**(nn+ll+mm)*fnn*fll*fmm)
        end do

        
!Generating primitive to orbital map.
!And pure to cartesian mapping up to G-type orbitals 

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

! max num prim per basis function
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

!Calculating overlap matrix
        call do_overlap()
! Calculating S^1/2 and S^-1/2
        allocate(s12p(nbasis,nbasis),s12m(nbasis,nbasis))
        call build_Smp(nbasis,s,s12m,s12p,0)

! dealocating auxiliary arrays
        deallocate(coefp,xnorm)
        deallocate(mnsh,iatsh,mssh)
        deallocate(expsh,c1,c2)

   END SUBROUTINE build_basis

   subroutine do_overlap()
   IMPLICIT DOUBLE PRECISION(A-H,O-Z)
   PARAMETER(PI=4.0d0*DATAN(1.0d0),TOL=1.0d-8)
   dimension AminusB(3)
   real*8, allocatable ::sp(:,:)

   allocate(sp(numprim,numprim))
   allocate(s(nbasis,nbasis))
   do ia=1,numprim
     do ib=1,ia
       sp(ia,ib)=0.0d0
       do ixyz=1,3 ! screening of primitives
        AminusB(ixyz)=coord(ixyz,iptoat(ia))-coord(ixyz,iptoat(ib))
        ii=mod(nlm(ia,ixyz)+nlm(ib,ixyz),2)
        if(abs(AminusB(ixyz)).lt.TOL.and.ii.ne.0) go to 111
       end do
       gamma_p=expp(ia)+expp(ib)
       eta_p=(expp(ia)*expp(ib))/gamma_p
       do_ov=PI**(3.0d0/2.0d0)/gamma_p**(3.0d0/2.0d0)
       do ixyz=1,3 !calculating for all 3 directions
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
           if(abs(AminusB(ixyz)).gt.TOL) then !avoid 0^0
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
! basis functions overlap
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


!!MMO- MODULE TESTING STARTS HERE.

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
! occ_no(igr,igr) -> ?
!Note: igr is number is basis functions.

   CONTAINS

   SUBROUTINE build_ao_matrices(igr)

   ALLOCATE(c(igr,igr),p(igr,igr))
   ALLOCATE(cb(igr,igr))
   ALLOCATE(ps(igr,igr),pa(igr,igr),pb(igr,igr))
   ALLOCATE(c_no(igr,igr),occ_no(igr,igr))

   END SUBROUTINE build_ao_matrices

   END MODULE ao_matrices



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
! pha -> rotated grid for zero-error. Change in theta?
! phb -> rotated grid for zero-error. Change in phi?
! rr00 -> nuclear distance at which half of the radial points have been distributed
! th(nang?) -> theta: angular coordinate for a given angular plane of points
! ph(nang?) -> phi: angular coordinate for a given angular plane of points
! w(nang?) -> probably mathematical weight of each angular plane?
! wr(nrad?) -> mathematical weight of a given radial surface
! xr(nrad?) -> distance from grid center to a given radial surface

   CONTAINS

   SUBROUTINE build_integration_grid(ienpart, ipolar,ifinegrid)
   character*30 integ1,integ2
!   INTEGER :: Nrad,Nang
!   DOUBLE PRECISION :: pha, phb, rr00

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
! Enpart defaults for high-accuracyone-el integrations
       else if(ienpart.eq.1.or.ipolar.eq.1.or.iedaiqa.eq.1) then
        nrad=150
        nang=590
        if(ifinegrid.eq.1) nang=974
        rr00=0.500d0
!       else if(ieoscent.eq.1) then
! simple grid to compute centroids of localized orbitals
!        nrad=20
!        nang=50
!        rr00=0.5d0
!! OSLOs testing !!
!       else if(ioslo.eq.1) then
!        nrad=40
!        nang=146
!        rr00=0.5d0
       else
! APOST old defaults for one-el integrations
        nrad=40
        nang=146
        rr00=0.5d0
       end if

       pha=0.0d0
       phb=0.0d0
!       pha=pha*dacos(-1.0d0)/180.0
!       phb=phb*dacos(-1.0d0)/180.0

       print *,' -------------------------------------'
       print *,' SETTING ATOMIC GRIDS FOR INTEGRATION '
       print *,' -------------------------------------'
       print *,' '
       write(*,'(A17,i5)') ' Radial points  :',nrad
       write(*,'(A17,i5)') ' Angular points :',nang
       write(*,'(A17,f7.3 )'),' r0 (radial)    :',rr00

       call quad(nrad,nang)


   END SUBROUTINE build_integration_grid

   END MODULE integration_grid
