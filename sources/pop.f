      subroutine fborder(sat)
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension sat(nbasis,nbasis,natoms)

      dimension rindex(maxat,maxat),tindex(maxat),diag(maxat)
      allocatable tt(:,:,:)
      allocate (tt(nbasis,nbasis,natoms))

C  COMPUTE THE MATRIX PRODUCTS P*S^A
C
      do iat=1,natoms
       do mu=1,nbasis
        do nu=1,nbasis
         x=0.d0
         do itau=1,nbasis
          x=x+P(mu,itau)*sat(itau,nu,iat)
         enddo
         tt(mu,nu,iat)=x
        enddo
       enddo
      enddo

      do iat=1,natoms
       do ibt=iat,natoms
        x=0.d0
        do mu=1,nbasis
         do nu=1,nbasis
           x=x+tt(mu,nu,iat)*tt(nu,mu,ibt)
         enddo
        enddo
        rindex(iat,ibt)=x
        rindex(ibt,iat)=x
       enddo
       diag(iat)=rindex(iat,iat)
      enddo

c Ps contribution
      if(kop.ne.0)then
       do iat=1,natoms
        do mu=1,nbasis
         do nu=1,nbasis
          x=0.d0
          do itau=1,nbasis
           x=x+Ps(mu,itau)*sat(itau,nu,iat)
          enddo
          tt(mu,nu,iat)=x
         enddo
        enddo
       enddo
       do iat=1,natoms
        do ibt=iat,natoms
         x=0.d0
         do mu=1,nbasis
          do nu=1,nbasis
           x=x+tt(mu,nu,iat)*tt(nu,mu,ibt)
          enddo
         enddo
         rindex(iat,ibt)=rindex(iat,ibt)+x
         rindex(ibt,iat)=rindex(iat,ibt)
        enddo
       enddo
      endif
c
c Savinig bo matrix. di will be overwritten in case of correlated calcualtion
      do i=1,natoms
       do j=1,natoms
        if (i.eq.j) then
         bo(i,i)=rindex(i,i)*0.5d0    
        else
         bo(i,j)=rindex(i,j)     
        end if
        di(i,j)=bo(i,j)     
       end do
      end do
c making zeroes for printing purposes
      do i=1,natoms
      rindex(i,i)=0.d0
      enddo

      WRITE(*,6342)
 6342 FORMAT(1x,/21X,'"FUZZY ATOMS" BOND ORDER MATRIX'//)
c      CALL Mprint(rindex,NATOMS,maxat)
      CALL Mprint(bo,NATOMS,maxat)
C
C CALCULATION OF THE VALENCE NUMBERS
C
      DO I=1,NATOMS
       X=0.D0
       DO  J=1,NATOMS
        X=X+rindex(I,J)
       end do
      rindex(I,I)=X
      end do
      do i=1,natoms 
       tindex(i)=rindex(i,i)
       diag(i)=2.d0*qat(i,1)-diag(i)
      enddo
      
      print *,'  '
      print *,'   TOTAL VALENCES '
      print *,'  '
      print *,'    Atom     V_A'
      print *,' -----------------'
      call vprint(diag,nat,maxat,1)
      print *,' ------------------'
      print *,'  '
      print *,'VALENCES USED IN BONDS'
      print *,'(SUM OF BOND ORDERS)'
      print *,'  '
      print *,'    Atom     VB_A'
      print *,' -----------------'
      call vprint(tindex,nat,maxat,1)
      print *,' ------------------'
      print *,'  '
      print *,'   FREE VALENCES'
      print *,' '
      print *,'    Atom     F_A'
      print *,' -----------------'
      do i=1,natoms
      diag(i)=diag(i)-tindex(i)
      enddo
      call vprint(diag,nat,maxat,1)
      print *,' ------------------'
      
      deallocate(tt)
      
      return
      end

      

      subroutine opop(wp,omp,omp2,rho)
      use basis_set, only: natoms
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common/actual/iact,jat,icenter
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(100)
      common /achi/achi(maxat,maxat),ibcp
      dimension wp(nrad*nang*natoms),rho(nrad*nang*natoms)
      dimension omp(nrad*nang*natoms),omp2(nrad*nang*natoms,natoms)


c IOPS
      idono=iopt(4)
      iallpo = iopt(7) 

      iatps=nrad*nang
      itotps=iatps*natoms

c igr= number of basis functions
c iatps = number of grid points per atom
c wp(i) = integration weight of the ith grid point
c chp(i,j) = value of the jth atomic orbital at the ith grid point
c omp(i) =  becke weight of the ith point of the atom to which the point belongs

c Computing  overlap population      
c each atom with its own grid points
       do i=1,natoms
        do j=1,natoms
         op(i,j)=0.0d0
        end do
       end do

       if(iallpo.eq.0) then
        do icenter=1,natoms
         do ifut=iatps*(icenter-1)+1,iatps*icenter
          do jcenter=icenter,natoms
           op(icenter,jcenter)=op(icenter,jcenter)+wp(ifut)*rho(ifut)*omp2(ifut,icenter)*omp2(ifut,jcenter)
          end do
         enddo
        enddo
       else if(iallpo.eq.1) then ! overriding rho
        do kcenter=1,natoms
         do ifut=iatps*(kcenter-1)+1,iatps*kcenter
           rho(ifut)=wp(ifut)*omp(ifut)*rho(ifut)
         end do
        end do
        do ifut=1,itotps 
         do icenter=1,natoms
          do jcenter=icenter,natoms
           op(icenter,jcenter)=op(icenter,jcenter)+rho(ifut)*omp2(ifut,icenter)*omp2(ifut,jcenter)
          end do
         enddo
        enddo
       end if 

       xx0=0.0d0
       do icenter=1,natoms
        xx0=xx0+op(icenter,icenter)
        do jcenter=icenter+1,natoms
         op(jcenter,icenter)=op(icenter,jcenter)
         xx0=xx0+2.0d0*op(icenter,jcenter)
        end do
       end do

      return 
      end

      subroutine fspindec(sat)
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /localspin/xlsa(maxat,maxat),ua(maxat)
      dimension sat(nbasis,nbasis,natoms)
c
      dimension rindex1(maxat,maxat), rindex2(maxat,maxat)
      dimension rindex3(maxat,maxat) ,rindex(maxat,maxat)
      dimension tt(:,:,:)
      dimension tts(:,:)
      allocatable tt
      allocatable tts

c
C
c
C    CALCULATING "FUZZY" BOND-ORDER AND VALENCE INDICES 
C  According to I. MAYER and P. SALVADOR, to be published  
C
c
C  Input parameters: Pa: Total electron density matrix;
c                    Pb: Spin density matrix;
c                    S: Overlap matrix; 
c                    Illim and Iulim: arrays of lower and upper limits of the 
c                                basis orbitals belonging to a given atom;
c                    Natoms: number of the atoms;
c                    Nbasis: number of basis orbitals.
C  
c     Uses also:     qsat(maxat,2): an array, the first column of which
c                    contains "fuzzy atom" populations of individual atoms
C  
c       
      allocate (tt(nbasis,nbasis,natoms))
      allocate (tts(nbasis,nbasis))

      do mu=1,nbasis
       do nu=1,nbasis
        tts(mu,nu)=0.0d0
       enddo 
      enddo

      do iat=1,natoms
       do mu=1,nbasis
        do nu=1,nbasis
         x=0.d0
         do itau=1,nbasis
          x=x+ps(mu,itau)*sat(itau,nu,iat)
         enddo
         tt(mu,nu,iat)=x
         tts(mu,nu)=tts(mu,nu)+x
        enddo
       enddo
      enddo

C Number of efectively unpaired electrons
      sum=0.0d0
      do iat=1,natoms
       x=0.0d0  
       do mu=1,nbasis
        do nu=1,nbasis
         x=x+tt(mu,nu,iat)*tts(nu,mu)
        enddo
       enddo
       ua(iat)=x
       sum=sum+x
      enddo
      print *,'  '
      print *,' EFFECTIVELY UNPAIRED ELECTRONS'
      print *,'  '
      print *,'    Atom     u_A'
      print *,' -----------------'
      call vprint(ua,nat,maxat,1)
      print *,' ------------------'
      write(*,'(a16,f10.5)') ' Sum check N_D = ' ,sum 

CCCCC
C DEPRECATED
CCCCC
      if(1.eq.0) then
C No U decomposition !!!!  
      sum=0.0d0
      do iat=1,natoms
       do ibt=iat,natoms
        x=0.0d0  
        do mu=1,nbasis
         do nu=1,nbasis
          x=x+tt(mu,nu,iat)*tt(nu,mu,ibt)
         enddo
        enddo
        rindex(iat,ibt)=x*0.50d0+(qsat(iat,1)*qsat(ibt,1))*0.25d0
        rindex(ibt,iat)=rindex(iat,ibt)
        sum=sum+rindex(iat,ibt)
        if(iat.ne.ibt) sum=sum+rindex(iat,ibt)
       enddo
      enddo


c     CALL Mprint(rindex,NATOMS,maxat)
c      WRITE(*,6342)
c 6342 FORMAT(1x,/21X,'"FUZZY ATOMS" S^2 DECOMPOSITION (a=0)'//)
c      CALL Mprint(rindex,NATOMS,maxat)
c      write(*,*) ' '
c      write(*,*) '<S^2> = ',sum


C     U Decomposition 
      sum=0.0d0
      do iat=1,natoms
      x=0.0d0
      do mu=1,nbasis
      do nu=1,nbasis
      x=x+tts(mu,nu)*tt(nu,mu,iat)
      enddo
      enddo
      rindex(iat,iat)=x*0.50d0+(qsat(iat,1)*qsat(iat,1))*0.25d0
      sum=sum+rindex(iat,iat)
      do ibt=iat+1,natoms
      rindex(iat,ibt)=(qsat(iat,1)*qsat(ibt,1))*0.25d0
      rindex(ibt,iat)=rindex(iat,ibt)
      sum=sum+2.0d0*rindex(iat,ibt)
      enddo
      enddo

c      WRITE(*,6352)
c 6352 FORMAT(1x,/21X,'"FUZZY ATOMS" S^2 DECOMPOSITION (a=1/2)'//)
c      CALL Mprint(rindex,NATOMS,maxat)
c      write(*,*) ' '
c      write(*,*) '<S^2> = ',sum

C   3/8  U Decomposition 

      sum=0.0d0
      do iat=1,natoms
      do jat=1,natoms
       rindex(iat,jat)=0.0d0
      end do 
      end do 


      do iat=1,natoms
      x=0.0d0
      do mu=1,nbasis
      do nu=1,nbasis
      x=x+tts(mu,nu)*tt(nu,mu,iat)
      enddo
      enddo
      rindex(iat,iat)=x*0.375d0
      do ibt=iat,natoms
      x=0.0d0
      do mu=1,nbasis
      do nu=1,nbasis
      x=x+tt(mu,nu,ibt)*tt(nu,mu,iat)
      enddo
      enddo
      rindex(iat,ibt)=rindex(iat,ibt)+
     &                x*0.125d0+(qsat(iat,1)*qsat(ibt,1))*0.25d0
      rindex(ibt,iat)=rindex(iat,ibt)
      sum=sum+rindex(iat,ibt)
      if(iat.ne.ibt) sum=sum+rindex(iat,ibt)
      enddo
      enddo

c      WRITE(*,6333)
c 6333 FORMAT(1x,/21X,'3/8 U "FUZZY ATOMS" SPIN DECOMPOSITION MATRIX'//)
c      CALL Mprint(rindex,NATOMS,maxat)
c      write(*,*) ' '
c      write(*,*) '<S^2> = ',sum

      end if
CCCCC
C DEPRECATED
CCCCC

C   3/4  U Decomposition 

      sum=0.0d0
      do iat=1,natoms
      do jat=1,natoms
       rindex(iat,jat)=0.0d0
       rindex1(iat,jat)=0.0d0
       rindex2(iat,jat)=0.0d0
      end do 
      end do 


      do iat=1,natoms
       x=0.0d0
       do mu=1,nbasis
        do nu=1,nbasis
         x=x+tts(mu,nu)*tt(nu,mu,iat)
        enddo
       enddo
       rindex(iat,iat)=x*0.750d0
       rindex3(iat,iat)=x*0.750d0
       do ibt=iat,natoms
        x=0.0d0
        do mu=1,nbasis
         do nu=1,nbasis
          x=x+tt(mu,nu,ibt)*tt(nu,mu,iat)
         enddo
        enddo
        rindex(iat,ibt)=rindex(iat,ibt)-x*0.25d0+(qsat(iat,1)*qsat(ibt,1))*0.25d0
        rindex(ibt,iat)=rindex(iat,ibt)
c      rindex1(iat,ibt)=rindex1(iat,ibt)-x*0.25d0
c      rindex2(iat,ibt)=rindex2(iat,ibt)+(qsat(iat,1)*qsat(ibt,1))*0.25d0
c      rindex1(ibt,iat)=rindex1(iat,ibt)
c      rindex2(ibt,iat)=rindex2(iat,ibt)
        sum=sum+rindex(iat,ibt)
        if(iat.ne.ibt) sum=sum+rindex(iat,ibt)
       enddo
      enddo

c deprecated
      if(1.eq.0) then
      print *,' '
      print *,'                  *** RECOMMENDED FORMULATION *** '
      WRITE(*,6090)
 6090 FORMAT(1x,/21X,'"FUZZY ATOMS" S^2 DECOMPOSITION (a=3/4)'//)
      CALL Mprint(rindex3,NATOMS,maxat)
      write(*,*) ' '
      write(*,*) '<S^2> = ',sum


      WRITE(*,6091)
 6091 FORMAT(1x,/21X,'3/4 U2 "FUZZY ATOMS" SPIN DECOMPOSITION MATRIX'//)
      CALL Mprint(rindex1,NATOMS,maxat)
      write(*,*) ' '
      write(*,*) '<S^2> = ',sum


      WRITE(*,6092)
 6092 FORMAT(1x,/21X,'3/4 U3 "FUZZY ATOMS" SPIN DECOMPOSITION MATRIX'//)
      CALL Mprint(rindex2,NATOMS,maxat)
      write(*,*) ' '
      write(*,*) '<S^2> = ',sum
      end if

c LSA
      print *,' '
      WRITE(*,6313)
 6313 FORMAT(1x,/21X,'"FUZZY ATOMS" S^2 DECOMPOSITION (a=3/4)'//)
      CALL Mprint(rindex,NATOMS,maxat)
      write(*,*) ' '
      write(*,'(a20,f10.5)') 'Sum check  <S^2> = ' ,sum
      write(*,*) ' '

      do i=1,natoms
       do j=1,natoms
        xlsa(i,j)=rindex(i,j)
       end do
      end do

C
C NOW DAVIDOSN
C
C
C  COMPUTE THE MATRIX PRODUCTS P*S^A
C
      do iat=1,natoms
      do mu=1,nbasis
      do nu=1,nbasis
      x=0.d0
      do itau=1,nbasis
      x=x+p(mu,itau)*sat(itau,nu,iat)
      enddo
      tt(mu,nu,iat)=x
      enddo
      enddo
      enddo

      do iat=1,natoms
      xx0=0.d0
      do ibt=1,natoms
 
      x=0.d0
      do mu=1,nbasis
      do nu=1,nbasis
      x=x+tt(mu,nu,iat)*tt(nu,mu,ibt)
      enddo
      enddo
      if(iat.ne.ibt) then
       rindex(iat,ibt)=rindex(iat,ibt)-3.0d0/8.0d0*x
       xx0=xx0+x
      end if
      
      enddo
      rindex(iat,iat)=rindex(iat,iat)+3.0d0/8.0d0*xx0
      enddo

      if(kop.ne.0)then

      do iat=1,natoms
      do mu=1,nbasis
      do nu=1,nbasis
      x=0.d0
      do itau=1,nbasis
      x=x+ps(mu,itau)*sat(itau,nu,iat)
      enddo
      tt(mu,nu,iat)=x
      enddo
      enddo
      enddo

      do iat=1,natoms
      xx0=0.d0
      do ibt=1,natoms
 
      x=0.d0
      do mu=1,nbasis
      do nu=1,nbasis
      x=x+tt(mu,nu,iat)*tt(nu,mu,ibt)
      enddo
      enddo

      if(iat.ne.ibt) then 
       rindex(iat,ibt)=rindex(iat,ibt)-3.0d0/8.0d0*x
       xx0=xx0+x
      end if
      
      enddo
      rindex(iat,iat)=rindex(iat,iat)+3.0d0/8.0d0*xx0
      enddo

      end if

      WRITE(*,6343)
 6343 FORMAT(1x,/21X,'"FUZZY ATOMS" DAVIDSON SPIN DEC. MATRIX'//)
      CALL Mprint(rindex,NATOMS,maxat)
      sum=0.0d0
      do iat=1,natoms
      do ibt=1,natoms
       sum=sum+rindex(iat,ibt)
      end do
      end do
      write(*,*) ' '
      write(*,'(a20,f10.5)') 'Sum check  <S^2> = ' ,sum
      write(*,*) ' '

      deallocate(tt,tts)
      
      return
      end

