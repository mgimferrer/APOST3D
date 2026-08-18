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

      call print_box('"FUZZY ATOMS" BOND ORDER MATRIX')
      CALL Mprint2(bo,NATOMS,maxat)
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
      
      call print_valence_table('TOTAL VALENCES',' ','Total valences',
     +  'V_A',diag)
      call print_valence_table('VALENCES USED IN BONDS',
     +  '(SUM OF BOND ORDERS)','Valences used in bonds','VB_A',tindex)
      do i=1,natoms
      diag(i)=diag(i)-tindex(i)
      enddo
      call print_valence_table('FREE VALENCES',' ','Free valences',
     +  'F_A',diag)

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
      common /iops/iopt(200)
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

!! ********************************************************************* !!
!! subroutine: population_density                                        !!
!! purpose: spin- and total-density atomic-population contraction        !!
!! (P*S^A products), writing qat/qsat via COMMON. Split from the print   !!
!! side (population_print_charges/population_print_overlap) so the       !!
!! idoint-guarded print_int call in main.f can sit between them,         !!
!! exactly where it did before this routine existed.                     !!
!! arguments:                                                            !!
!! sat (in) -- per-atom AO overlap matrix (numint_sat/tomull/tolow)      !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine population_density(sat)
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension sat(nbasis,nbasis,natoms)

      call print_box('DOING POPULATION ANALYSIS')
      write(*,'(2x,a)') 'Partial atomic charges'
      write(*,'(2x,a)') 'Atomic spin densities'
      write(*,'(2x,a)') 'Bond orders and Valences'

!! spin density: P^s*S^A contraction, one atom per outer iteration       !!
      if(kop.ne.0.or.nalf.ne.nb) then

!! parallel over kat (natoms) -- each iteration reduces mu,nu into a     !!
!! private x and writes only its own qsat(kat,1) slot; ps/sat shared     !!
!! read-only, no cross-iteration dependency.                             !!
!$OMP PARALLEL DO PRIVATE(mu,nu,x)
        do kat=1,natoms
          x=0.d0
          do mu=1,nbasis
            do nu=1,nbasis
              x=x+ps(mu,nu)*sat(mu,nu,kat)
            enddo
          enddo
          qsat(kat,1)=x
          qsat(kat,2)=0.d0
        enddo
!$OMP END PARALLEL DO

!! second (Ps*S) contribution, O(igr) outer trips only -- left serial,   !!
!! not worth threading, and kat=ihold(mu) can repeat across mu so the    !!
!! qsat(kat,2) accumulation isn't race-free without extra care.          !!
        do mu=1,nbasis
          kat=ihold(mu)
          x=0.d0
          do itau=1,nbasis
            x=x+ps(mu,itau)*s(itau,mu)
          enddo
          qsat(kat,2)=qsat(kat,2)+x
        enddo
      end if

!! total density: P*S^A contraction, same pattern as spin density above  !!
!$OMP PARALLEL DO PRIVATE(mu,nu,x)
      do kat=1,natoms
        x=0.d0
        do mu=1,nbasis
          do nu=1,nbasis
            x=x+p(mu,nu)*sat(mu,nu,kat)
          enddo
        enddo
        qat(kat,1)=x
        qat(kat,2)=0.d0
      enddo
!$OMP END PARALLEL DO

!! second (P*S) contribution -- same O(igr)-only reasoning as above      !!
      do mu=1,nbasis
        kat=ihold(mu)
        x=0.d0
        do itau=1,nbasis
          x=x+p(mu,itau)*s(itau,mu)
        enddo
        qat(kat,2)=qat(kat,2)+x
      enddo

      return
      end

!! ********************************************************************* !!
!! subroutine: print_population_table                                    !!
!! purpose: shared title/column-header/data/sum/fragment-breakdown       !!
!! layout for the three 2-column (apost3d, Mulliken) atomic tables       !!
!! (electron populations, atomic charges, spin populations) -- replaces  !!
!! three near-identical copies of this block, previously inconsistent    !!
!! in small ways (a stray "apost3D" instead of "apost3d" in one of the   !!
!! three, an extra blank line before the fragment breakdown in two of    !!
!! the three but not the first). Header/separator/sum-line widths are    !!
!! sized to match vprint's own fixed data-row format exactly (32 columns !!
!! for the default 6-decimal case, 48 for FULLPRECISION) -- previously a  !!
!! hardcoded 29/35-column layout that didn't cover the actual numbers.    !!
!! arguments:                                                            !!
!! title    (in) -- print_box title, e.g. 'ELECTRON POPULATIONS'         !!
!! fraglabel(in) -- fragment-breakdown label, e.g. 'Electron populations'!!
!! arr      (in) -- the (maxat,2) table to print (qat/dummyvec/qsat)     !!
!! tc1,tc2  (in) -- apost3d/Mulliken column sums, computed by the caller !!
!!                  (the summation itself differs per table)             !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine print_population_table(title,fraglabel,arr,tc1,tc2)
      use input_options_mod, only: idofr
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /printout/iaccur
      character*(*) title,fraglabel
      dimension arr(maxat,2)
      character*80 line

      call print_box(title)
      if(iaccur.eq.0) then
        write(*,'(1x,a7,2a12)') 'Atom','apost3d','Mulliken'
        write(*,'(2x,a)') repeat('-',30)
      else
        write(*,'(1x,a7,2a20)') 'Atom','apost3d','Mulliken'
        write(*,'(2x,a)') repeat('-',46)
      end if
      call vprint(arr,nat,maxat,2)
      if(iaccur.eq.0) then
        write(*,'(2x,a)') repeat('-',30)
        write(*,162) tc1,tc2
      else
        write(*,'(2x,a)') repeat('-',46)
        write(*,172) tc1,tc2
      end if

      if (idofr.eq.1) then
        line ='   FRAGMENT ANALYSIS : '//fraglabel
        call group_by_frag_vec(2,line,arr)
      end if

      return

 162  format(1x,'    Sum',2f12.6)
 172  format(1x,'    Sum',2f20.13)

      end

!! ********************************************************************* !!
!! subroutine: print_valence_table                                       !!
!! purpose: shared title/column-header/data layout for the three         !!
!! single-column atomic valence tables printed by fborder (V_A/VB_A/     !!
!! F_A) -- same consolidation print_population_table already applies to  !!
!! the electron/charge/spin population tables, and for the same reason:  !!
!! the header/separator width now matches vprint's own fixed data-row    !!
!! format (20 columns) instead of a hardcoded value that didn't cover    !!
!! the numbers. Also prints a per-fragment breakdown (DOFRAGS) via        !!
!! group_by_frag_vec, matching the population/bond-order tables --        !!
!! previously the only ones of the six population-analysis tables         !!
!! without one.                                                           !!
!! arguments:                                                            !!
!! title    (in) -- print_box title, e.g. 'TOTAL VALENCES'               !!
!! subtitle (in) -- optional explanatory line under the title, blank     !!
!!                   (' ') if none, e.g. '(SUM OF BOND ORDERS)'          !!
!! fraglabel(in) -- fragment-breakdown label, e.g. 'Total valences'      !!
!! label    (in) -- column header, e.g. 'V_A'                            !!
!! arr      (in) -- the (maxat) vector to print (diag/tindex)            !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine print_valence_table(title,subtitle,fraglabel,label,arr)
      use input_options_mod, only: idofr
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      character*(*) title,subtitle,fraglabel,label
      dimension arr(maxat)
      character*80 line

      call print_box(title)
      if(len_trim(subtitle).gt.0) then
        write(*,'(1x,a)') trim(subtitle)
        write(*,*)
      end if
      write(*,'(1x,a7,a12)') 'Atom',label
      write(*,'(2x,a)') repeat('-',18)
      call vprint(arr,nat,maxat,1)
      write(*,'(2x,a)') repeat('-',18)

      if (idofr.eq.1) then
        line=' FRAGMENT ANALYSIS: '//fraglabel
        call group_by_frag_vec(1,line,arr)
      end if

      return
      end

!! ********************************************************************* !!
!! subroutine: population_print_charges                                  !!
!! purpose: prints electron populations and partial atomic charges       !!
!! (plus per-fragment breakdowns where DOFRAGS is set) from qat, already !!
!! filled by population_density. The ELCOUNT/NCTAIM call that used to    !!
!! sit between these two prints stays in main.f (like idoint/print_int)  !!
!! -- NCTAIM lives in subroutines_mmo.f, which apost3d-eos's link list   !!
!! deliberately excludes, so pop.f itself must not reference it.         !!
!! arguments: none (all via qat/COMMON)                                  !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine population_print_charges()
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension dummyvec(maxat,2)

!! electron populations -- straight sum of qat !!
      tc1=0.d0
      tc2=0.d0
      do i=1,nat
        tc1=tc1+qat(i,1)
        tc2=tc2+qat(i,2)
      enddo
      call print_population_table('ELECTRON POPULATIONS',
     +  'Electron populations',qat,tc1,tc2)

!! partial charges -- nuclear charge minus qat, atom by atom !!
      tc1=0.d0
      tc2=0.d0
      do i=1,nat
        tc1=tc1-qat(i,1)+zn(i)
        tc2=tc2-qat(i,2)+zn(i)
        dummyvec(i,1)=zn(i)-qat(i,1)
        dummyvec(i,2)=zn(i)-qat(i,2)
      enddo
      call print_population_table('TOTAL ATOMIC CHARGES',
     +  'Atomic Charges',dummyvec,tc1,tc2)

      return
      end

!! ********************************************************************* !!
!! subroutine: population_print_overlap                                  !!
!! purpose: prints spin populations (plus per-fragment breakdown) and    !!
!! overlap populations from qsat/op, already filled by                   !!
!! population_density -- op itself is built here via opop/mull_opop.     !!
!! arguments:                                                            !!
!! wp,omp,omp2,rho (in) -- integration grid weights/density, passed      !!
!! straight to opop                                                      !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine population_print_overlap(wp,omp,omp2,rho)
      use basis_set, only: natoms
      use integration_grid
      use input_options_mod, only: idofr,imulli,iopop,icorr
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      dimension wp(nrad*nang*natoms),rho(nrad*nang*natoms)
      dimension omp(nrad*nang*natoms),omp2(nrad*nang*natoms,natoms)
      character*80 line

!! spin populations -- straight sum of qsat, only for open-shell/       !!
!! correlated-with-DM cases                                              !!
      if(kop.ne.0.OR.(icas.eq.1.and.nalf.ne.nb.and.icorr.ne.0)) then
        tc1=0.d0
        tc2=0.d0
        do i=1,nat
          tc1=tc1+qsat(i,1)
          tc2=tc2+qsat(i,2)
        enddo
        call print_population_table('SPIN POPULATIONS',
     +    'Spin Populations',qsat,tc1,tc2)
      end if

!! overlap populations -- op itself comes from opop (real-space) or     !!
!! mull_opop (Mulliken); iopop=0 just takes the diagonal from qat        !!
      if(iopop.eq.0) then
        do i=1,nat
          op(i,i)=qat(i,1)
        end do
      else if(imulli.ne.1) then
        call opop(wp,omp,omp2,rho)
      else
        call mull_opop()
      end if

      if(iopop.eq.1) then
        call print_box('APOST3D OVERLAP POPULATION MATRIX')
        call mprint2(op,nat,maxat)
        if (idofr.eq.1) then
          line ='   FRAGMENT ANALYSIS : Overlap Populations'
          call group_by_frag_mat(0,line,op)
        end if
      end if

      return
      end

!! ********************************************************************* !!
!! subroutine: bond_order_analysis                                       !!
!! purpose: thin wrapper around fborder -- kept separate from            !!
!! population_density/population_print_charges/population_print_overlap  !!
!! on purpose, so a future second bond-order scheme gets its own         !!
!! equally-thin subroutine instead of being wedged into population math. !!
!! arguments:                                                            !!
!! sat (in) -- per-atom AO overlap matrix, passed straight to fborder    !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine bond_order_analysis(sat)
      use basis_set, only: natoms,nbasis
      use input_options_mod, only: idofr
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      dimension sat(nbasis,nbasis,natoms)
      character*80 line

      call fborder(sat)
      if (idofr.eq.1) then
        line ='   FRAGMENT ANALYSIS : Fuzzy Bond Order'
        call group_by_frag_mat(1,line ,bo)
      end if

      return
      end

!! ********************************************************************* !!
!! subroutine: pca_analysis                                              !!
!! purpose: "fuzzy atoms" PCA -- diagonalizes a covariance-like matrix   !!
!! built from the overlap population (op) and delocalization index (di)  !!
!! matrices, printing eigenvectors/eigenvalues and their projection      !!
!! against qat. Needs di already populated (fborder, called earlier via  !!
!! bond_order_analysis).                                                 !!
!! FIXED 2026-08-15 (M. Gimferrer): scr, the eigenvector-output argument !!
!! to diagonalize, was allocated as scr(nat) -- a vector -- but          !!
!! diagonalize needs a full (M,M) matrix there, same as every other call !!
!! site in this codebase. Wrote past scr's allocation whenever nat>1     !!
!! (heap overflow, never caught -- PCA has zero test coverage). Now      !!
!! allocated scr(nat,nat). Confirmed output-preserving: diagonalize's    !!
!! write order fills scr's first nat elements (column 1) before any      !!
!! out-of-bounds write happens, and every print below only ever read     !!
!! scr(1:nat) -- i.e. column 1 -- so the values printed are identical    !!
!! before and after, only the undefined-behavior overflow is gone.       !!
!! STILL OPEN, not touched -- needs M. Gimferrer + P. Salvador to        !!
!! confirm intent first: the two prints right after diagonalize look     !!
!! swapped relative to its actual contract (pca/A0 comes back with       !!
!! eigenvalues on its diagonal, scr/X holds the eigenvectors) -- but the !!
!! code labels pca "PCA EIGENVECTORS" and prints scr as if it held       !!
!! eigenvalues. See CLAUDE.md Known Issues.                              !!
!! arguments: none (all via op/di/qat/COMMON)                            !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine pca_analysis()
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      allocatable pca(:,:),scr(:,:)

      allocate(pca(nat,nat))
      allocate(scr(nat,nat))
      do i=1,nat
        do j=1,nat
          pca(i,j)=op(i,j)-0.50d0*di(i,j)
        end do
      end do

      call print_box('"FUZZY ATOMS" COVARIANCE MATRIX')
      call mprint(pca,nat,nat)
      write(*,*)

      call diagonalize(nat,nat,pca,scr,0)

      call print_box('"FUZZY ATOMS" PCA EIGENVECTORS')
      call mprint(pca,nat,nat)
      write(*,*)
      write(*,'(8f10.4)') (scr(i,1),i=1,nat)
      write(*,*)

      do i=1,nat
        xx=ZERO
        do k=1,nat
          xx=xx+pca(k,i)*qat(k,1)
        end do
        write(*,'(2x,a,i3,a,2f14.6)') 'PC: ',i,' sum: ',xx,xx*scr(i,1)
      end do

      deallocate(pca,scr)

      return
      end

