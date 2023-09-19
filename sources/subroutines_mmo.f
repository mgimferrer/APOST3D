!NOTE: filename (subroutines_mmo.f90) is just a placeholder to easily auto-complete after typing a single letter.
      SUBROUTINE NCTAIM(iatps,wp,omp,omp2,nat,sat,igr,ibaspoint,chp)
      include 'parameter.h'
      DOUBLE PRECISION, ALLOCATABLE :: NCTAIM_qat(:)
      DOUBLE PRECISION :: omp(iatps)
      double precision, dimension(nat*iatps,nat) :: omp2
      double precision, dimension(nat*iatps) :: wp
      double precision, dimension(igr,igr,nat) ::sat
      dimension ibaspoint(iatps*nat)
      double precision, dimension(iatps*nat,igr) ::chp
      common /iops/iopt(100)

      iradmat = iopt(48)
      iunit=32

       ALLOCATE(NCTAIM_qat(nat))
       CALL NCTAIM_reference(NCTAIM_qat,iunit)
       if(iradmat.eq.2) CALL NCTAIM_HirshIt_loop(NCTAIM_qat,omp,omp2,wp,nat,iatps,sat,igr,ibaspoint,chp)
       DEALLOCATE(NCTAIM_qat)

      END SUBROUTINE

      SUBROUTINE NCTAIM_reference(NCTAIM_qat,iunit)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /iops/iopt(100)
      DOUBLE PRECISION, DIMENSION(nat) :: NCTAIM_qat

      iradmat = iopt(48)

       write(*,*) '      +-----------------------------+      '
       write(*,*) '      |     NCTAIM 3D PARTITION     |      '
       write(*,*) '      |  (No-Charge-Transfer AIMs)  |      '
       write(*,*) '      +-----------------------------+      '
       write(*,*)
       write(*,*) 'Partitioning real-space ensuring constant atomic Q (Hirshfeld-Iterative-like version).'
       if(iradmat.eq.1) then
         write(*,*) 'Writing reference file (elcount.inp), containing all atomic populations at F=0.'
         OPEN(UNIT=iunit,FILE="elcount.inp")
         write(iunit,*) (qat(iatom,1),iatom=1,nat)
         CLOSE(iunit)
       else if(iradmat.eq.2) then
         write(*,*) 'Reading reference electron populations of all atoms at F=0 (from reference file elcount.inp).'
         OPEN(UNIT=iunit,FILE="elcount.inp",STATUS="OLD")
         read(iunit,*) (NCTAIM_qat(iatom),iatom=1,nat)
         CLOSE(iunit)
       else if(iradmat.eq.0) then
         write(*,*) 'Incompatible options ELCOUNT=1,RADMAT=0. Stopping.'
         STOP 'Incompatible options. Problem in ELCOUNT user input.'
       end if
      END SUBROUTINE

      SUBROUTINE NCTAIM_HirshIt_loop(NCTAIM_qat,omp,omp2,wp,nat,iatps,sat,igr,ibaspoint,chp)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameter.h'
      DOUBLE PRECISION, DIMENSION(nat) :: NCTAIM_qat
      DOUBLE PRECISION, ALLOCATABLE :: atomic_X(:), prev_atomic_X(:), prev_atomic_Q(:), atomic_err(:)
      DOUBLE PRECISION :: omp(iatps)
      INTEGER, ALLOCATABLE :: int_l(:), int_u(:)
      CHARACTER*3, ALLOCATABLE :: atom_converged(:)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      INTEGER, DIMENSION(0:6) :: PRINT_INFO
      DATA PRINT_INFO/+2,+2,+1,0,-1,-2,-2/
      double precision, dimension(nat*iatps,nat) :: omp2
      double precision, dimension(nat*iatps) :: wp
      dimension sat(igr,igr,nat)
      dimension ibaspoint(iatps*nat),chp(iatps*nat,igr)
       ALLOCATE(atomic_X(nat),prev_atomic_X(nat),prev_atomic_Q(nat),atomic_err(nat),int_l(nat),int_u(nat),atom_converged(nat))
!atomic_X -> the calculation of the weight function is using, for each atom, a combination of two states. For each atom, this is the parameter that determines how much its weight function is using of its upper reference (the rest, i.e., 1-atomic_X is the fraction that corresponds to the atom with one less electron)
!prev_atomic_X -> dummy variable that stores the atomic_X values (explained above) of the previous iteration
!prev_atomic_Q -> dummy variable that stores the atomic charges of the previous iteration
!atomic_err(nat) -> for each atomic charge, difference from reference charge
!int_l -> lower limit to the observed population (e.g. reference with 9 electrons if our population is 9.2)
!int_u -> upper limit to the observed population (e.g. 10 electrons if our population is 9.2)
!atom_converged -> character variable that stores whether an atom has converged or not, in the current iteration

       write(*,*)
       write(*,*) '      +-----------------------------+      '
       write(*,*) '      |  Iterative process begins   |      '
       write(*,*) '      +-----------------------------+      '
       write(*,*)

       ielcount_iteration=0 !iteration counter
       nconv=0 !number of converged atoms, resets every iteration
       xconv_error=0.0d0 !sum over all atoms of difference from convergence
       itotps=iatps*nat

C <NCTAIM> Preparation: convergence criterion.

       xx=0.0d0
       do iatom=1,nat
         xx=xx+qat(iatom,1)
       end do
       xx=ABS(xx-NINT(xx))
!TEST       write(*,*) 'Considering the accuracy for the total electron count, requested convergence for the iterative NCTAIM process will be ',xx
!TEST       xconvergence_criterion=xx
       write(*,*) 'Testing: arbitrary convergence of 1E-5.'
       xconvergence_criterion=1E-5

C <NCTAIM> Preparation: zeroth iteration.

       write(*,*) 'Zeroth iteration: checking for convergence with default atomic regions...'
       write(*,*)
       do iatom=1,nat
         xx=ABS(qat(iatom,1)-NCTAIM_qat(iatom))
         xconv_error=xconv_error+xx
         if(xx.lt.xconvergence_criterion) nconv=nconv+1
       end do
       write(*,'(x,I2,a,I2,a)') nconv,' atoms out of',nat,' match reference populations.'
       write(*,*) 'Sum of individual atomic errors: ',xconv_error
       write(*,*)
       if(nconv.eq.nat) then
         write(*,*) 'Unexpectedly, populations match using the default weight functions.'
       else
         write(*,*) 'Electron populations change with default atomic regions.'
         write(*,*) 'Atom # | Charge           | Lower int. | Upper int. | Molar fraction X'
         do iatom=1,nat
           xcharge=zn(iatom)-qat(iatom,1)
           int_l(iatom)=INT(qat(iatom,1))-INT(zn(iatom))+3
           int_u(iatom)=int_l(iatom)+1
           atomic_X(iatom)=qat(iatom,1)-INT(qat(iatom,1)) !fraction of species of the upper limit (e.g. 0.2 if our population is 9.2)
           write(*,'(3x,I2.2,2x,a,F16.13,a,4x,I2,4x,a,4x,I2,4x,a,F15.10)') iatom,' | ',xcharge,' | ',PRINT_INFO(int_l(iatom)),' | '
     +,PRINT_INFO(int_u(iatom)), ' | ',atomic_X(iatom)
         end do
         write(*,*)
         CALL NCTAIM_HirshIt_wat(iatps,omp2,nat,atomic_X,int_l,int_u) !obtain updated omp2
         ifut=1
         do icenter=1,nat
          do k=1,iatps
           omp(ifut)=omp2(ifut,icenter)
           ifut=ifut+1
          end do
         end do
         CALL numint_sat(igr,itotps,nat,wp,omp,omp2,chp,ibaspoint,sat)
         CALL Population_analysis(sat,nat,igr,qat)
         write(*,*)
         write(*,*) 'Since this is the 0th iteration, making a directed change (of arbitrary magnitude) to the molar fraction (X) 
     +of each atom, to evaluate rate of change of the population.'
         write(*,*) 'Atom # | Initial molar fraction X | Arbitrary guess for 1st iteration'
         do iatom=1,nat
           xguess_direction=qat(iatom,1)-NCTAIM_qat(iatom)
           prev_atomic_X(iatom)=atomic_X(iatom) !store 0th iteration X before updating it
           if(xguess_direction.gt.0) then !population obtained is too large
             atomic_X(iatom)=atomic_X(iatom)-0.01d0
             if(atomic_X(iatom).lt.0) atomic_X(iatom)=prev_atomic_X(iatom)*0.99d0
!alternatively, shift reference integers by +1
!               atomic_X(iatom)=atomic_X(iatom)+1.0d0
!               int_l(iatom)=int_l(iatom)+1
!               int_u(iatom)=int_u(iatom)+1
!             end if
           else if(xguess_direction.lt.0) then !population obtained is too small
             atomic_X(iatom)=atomic_X(iatom)+0.01d0
             if(atomic_X(iatom).gt.1) atomic_X(iatom)=prev_atomic_X(iatom)*1.01d0
!alternatively, shift reference integers by -1
!               atomic_X(iatom)=atomic_X(iatom)-1.0d0
!               int_l(iatom)=int_l(iatom)-1
!               int_u(iatom)=int_u(iatom)-1
!             end if
           end if
           write(*,'(3x,I2.2,2x,a,F19.13,5x,a,F19.13)') iatom,' | ',prev_atomic_X(iatom),' | ',atomic_X(iatom)
           prev_atomic_Q(iatom)=qat(iatom,1) !store 0th iteration Q before exiting subroutine
         end do
         write(*,*)
       end if

C <NCTAIM> Preparation done. Main algorithm.

       do while(nconv.lt.nat)
         write(*,*)
         CALL NCTAIM_HirshIt_wat(iatps,omp2,nat,atomic_X,int_l,int_u) !obtain updated omp2
         ifut=1
         do icenter=1,nat
          do k=1,iatps
           omp(ifut)=omp2(ifut,icenter)
           ifut=ifut+1
          end do
         end do
         CALL numint_sat(igr,itotps,nat,wp,omp,omp2,chp,ibaspoint,sat)
         CALL Population_analysis(sat,nat,igr,qat)
         write(*,*)
         ielcount_iteration=ielcount_iteration+1
         nconv=0
         xconv_error=0.0d0
         write(*,'(x,a,I3.3,a)') ' <<< Iteration ',ielcount_iteration,' >>>'

C STEP 1: CHECK FOR CONVERGENCE.

         do iatom=1,nat
           atom_converged(iatom)="NO"
           atomic_err(iatom)=ABS(qat(iatom,1)-NCTAIM_qat(iatom))
           xconv_error=xconv_error+atomic_err(iatom)
           if(atomic_err(iatom).lt.xconvergence_criterion) then
             nconv=nconv+1
             atom_converged(iatom)="YES"
           end if
         end do
         write(*,'(x,I2,a,I2,a)') nconv,' atoms out of',nat,' match reference populations.'
         write(*,'(x,a,F20.16)') 'Sum of individual atomic errors: ',xconv_error
         write(*,*)
         write(*,*) 'Atom # |    Current X    |     Delta X     |    Delta Q      |      dQ/dX      | Predicted shift || Converged?
     + |    Atom error '

         do iatom=1,nat
C STEP 2: Check change in molar fraction X.
           xdelta_X=atomic_X(iatom)-prev_atomic_X(iatom)

C STEP 3: Check change in population.
           xdelta_Q=qat(iatom,1)-prev_atomic_Q(iatom)

C STEP 4: Calculate slope of the change in Q caused by a change in X.
           xelcount_slope=xdelta_Q/xdelta_X

C STEP 5: Extrapolate shift in X required for atomic charge to match reference.
           xshift=(NCTAIM_qat(iatom)-qat(iatom,1))/xelcount_slope 

C STEP 6: Store current value of X.
           prev_atomic_X(iatom)=atomic_X(iatom)

C STEP 7: Update value of X, so it can be used in the upcoming iteration.
           atomic_X(iatom)=atomic_X(iatom)+xshift

C STEP 8: Store current value of Q.
           prev_atomic_Q(iatom)=qat(iatom,1)

           write(*,'(3x,I2.2,2x,a,F15.12,a,F15.12,a,F15.12,a,F15.09,a,F15.12,a,4x,a,3x,a,F15.12)') iatom,' | ',prev_atomic_X(iatom
     +),' | ',xdelta_X,' | ',xdelta_Q,' | ',xelcount_slope,' | ',xshift,' || ',atom_converged(iatom),' | ',atomic_err(iatom)
         end do

         if(ielcount_iteration.ge.50) then
           write(*,*) 'Convergence problem, maximum number of iterations reached.'
           STOP 'Convergence problem'
         end if
       end do !cycle main algorithm until convergence achieved
      END SUBROUTINE

      SUBROUTINE NCTAIM_HirshIt_wat(iatps,whi,nat,atomic_X,int_l,int_u)
      use integration_grid
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameter.h'
      common /hirsh/y2(50,5,150),radial(50,5,150),xr2(150),ieq(maxat),nradat,nat0,pop(maxat)
      common /coord/coord(3,maxat),zn(maxat),iznuc(maxat)

c      double precision, dimension(maxat) :: pop
      dimension ipollas(0:6)
      double precision, dimension(nat*iatps,nat) :: whi

      DOUBLE PRECISION, DIMENSION(nat) :: xatom_pop
      DIMENSION :: atomic_X(nat)
      INTEGER, DIMENSION(nat) :: int_l, int_u

       xLOWTOLER=1.0D-12
       ipollas(0)=5
       ipollas(1)=5
       ipollas(2)=3
       ipollas(3)=1
       ipollas(4)=2
       ipollas(5)=4
       ipollas(6)=4

       ifut=0
       do icenter=1,nat
         do k=1,nrad
           do nn=1,nang
             ifut=ifut+1
             thx=th(nn)
             fix=ph(nn)
             rr=xr(k)
             xabs0=rr*dsin(thx)*dcos(fix)
             yabs0=rr*dsin(thx)*dsin(fix)
             zabs0=rr*dcos(thx)
             xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
             yabs=-dcos(pha)*dsin(phb)*xabs0+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
             zabs=dsin(pha)*dsin(phb)*xabs0-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0
             xabs=xabs+coord(1,icenter)
             yabs=yabs+coord(2,icenter)
             zabs=zabs+coord(3,icenter)
!WATHIRSH2 FUNCTION TO OBTAIN whi(ifut,icenter), ignoring off-atom terms (i.e. not full omp2)
             xpoint_population=0.0d0
             do iatom=1,nat
               dist=DSQRT((xabs-coord(1,iatom))**2.0d0+(yabs-coord(2,iatom))**2.0d0+(zabs-coord(3,iatom))**2.0d0)
               xatom_pop(iatom)=atomic_X(iatom)*splint(ieq(iatom),ipollas(int_u(iatom)),nradat,dist)+(1-atomic_X(iatom))*
     +splint(ieq(iatom),ipollas(int_l(iatom)),nradat,dist)
               xpoint_population=xpoint_population+xatom_pop(iatom)
             end do
             if(xpoint_population.lt.xLOWTOLER) then
               whi(ifut,icenter)=0.0d0
             else
               whi(ifut,icenter)=xatom_pop(icenter)/xpoint_population
             end if
           end do
         end do
       end do
      END SUBROUTINE

      SUBROUTINE calc_twoel_analytical(itotps,wp,omp2,pcoord,rho,chp2,chp2b,coul,exch_hf)
      USE basis_set
      USE ao_matrices
      USE integration_grid
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot

      dimension :: wp(itotps),pcoord(itotps,3)
      dimension :: omp2(itotps,nat),rho(itotps)
      dimension :: coul(maxat,maxat)
      dimension :: exch_hf(maxat,maxat)
      DIMENSION :: chp2(itotps,500)
      DIMENSION :: chp2b(itotps,500)

      PARAMETER(PII=4.0d0*DATAN(1.0d0))
      DIMENSION E(0:10,0:10,0:10)
      DIMENSION E_coulomb(3,0:10)
      DIMENSION max_ang(3),xnew_coord(3)
      DIMENSION Boys(0:20)
      DIMENSION xR_Boys(0:20,0:10,0:10,0:10)

      iatps = nrad*nang
      idofr = Iopt(40)
      nat   = natoms
      N_primitives = numprim

!!! LOOPING OVER LOWER TRIANGLE TO CALCULATE <Ga|1/rc|Gb> !!!

       coul=0.0d0
       exch_hf=0.0d0
       do ia=1,N_primitives
         iAO_a=iptob_cartesian(ia)
         iatom_a=iptoat(ia)
         prefactor0=coefpb(ia,iAO_a)*2.0d0*PII
         do ib=1,ia
           iAO_b=iptob_cartesian(ib)
           iatom_b=iptoat(ib)
           exp_sum=expp(ia)+expp(ib)
           prefactor1=prefactor0*coefpb(ib,iAO_b)/exp_sum
C -- To account for the full matrix while doing only lower diagonal, counting twice if primitives aren't the same and belong to the same orbital (diagonal term in AO matrix).
           if(iAO_a.eq.iAO_b.and.ia.ne.ib) prefactor1=2.0d0*prefactor1
C -- With the same objective, off-diagonal terms (in the AO matrix) obviously multiplied by 2.
           if(iAO_a.ne.iAO_b) prefactor1=2.0d0*prefactor1
           do idir=1,3
             gauss_dist=coord(idir,iatom_a)-coord(idir,iatom_b) !distance between the 2 gaussians
             if(ABS(gauss_dist).lt.1.0d-15) gauss_dist=0.0d0 !TOL_DIST
             xnew_coord(idir)=(expp(ia)*coord(idir,iatom_a)+expp(ib)*coord(idir,iatom_b))/exp_sum !coordinates of the new gaussian,
!    + sum of the 2
             xdist_a=xnew_coord(idir)-coord(idir,iatom_a) !distance from gaussian A to gaussian sum
             xdist_b=xnew_coord(idir)-coord(idir,iatom_b) !distance from gaussian B to gaussian sum
             imax=nlm(ia,idir)
             jmax=nlm(ib,idir)
             max_ang(idir)=imax+jmax
             xx=(expp(ia)*expp(ib))/exp_sum
             E(0,0,0)=exp(-xx*gauss_dist**2.0d0)
             if(imax.eq.0.and.jmax.eq.0) E_coulomb(idir,0)=E(0,0,0)
             do ii=0,imax-1
               do it=0,ii+1
                 if (((it-1).lt.0).or.((it-1).gt.ii)) then
                   E1=0
                 else
                   E1=E(ii,0,it-1)/(2.0d0*exp_sum)
                 end if
                 if (it.gt.ii) then
                   E2=0
                 else
                   E2=E(ii,0,it)*xdist_a
                 end if
                 if ((it+1).gt.ii) then
                   E3=0
                 else
                   E3=real((it+1))*E(ii,0,it+1)
                 end if
                 E(ii+1,0,it)=E1+E2+E3
                 if(ii.eq.imax-1.and.jmax.eq.0) E_coulomb(idir,it)=E(ii+1,0,it)
               end do
             end do

             do jj=0,jmax-1
               do it=0,imax+jj+1
                 if (((it-1).lt.0).or.((it-1).gt.(jj+imax))) then
                   E1=0
                 else
                   E1=E(imax,jj,it-1)/(2.0d0*exp_sum)
                 end if
                 if (it.gt.(jj+imax)) then
                   E2=0
                 else
                   E2=E(imax,jj,it)*xdist_b
                 end if
                 if ((it+1).gt.(jj+imax)) then
                   E3=0
                 else
                   E3=real((it+1))*E(imax,jj,it+1)
                 end if
                 E(imax,jj+1,it)=E1+E2+E3
                 if(jj+1.eq.jmax) E_coulomb(idir,it)=E(imax,jmax,it)
               end do
             end do
             E_ij=E(imax,jmax,0) !never used
           end do !end of idir loop. it only depends on the 2 gaussians being considered (regardless of position of point charge), so the best point to loop through all points charges (grid points) would probably be here (right after n_max)

!!! R RECURRENCE FOR COULOMB BEGINS, USING THE SAVED E VALUES (ALL VALUES OF T FOR IMAX,JMAX) ACROSS ALL 3 DIRECTIONS !!!
!!! CALCULATING BOYS FUNCTION OF ORDER n_max (lower n has to be obtained through recurrence) !!!

           n_max=max_ang(1)+max_ang(2)+max_ang(3)
           do icenter=1,nat
             do ifut=iatps*(icenter-1)+1,iatps*icenter
               x0=wp(ifut)*omp2(ifut,icenter) !MG
               dx0=pcoord(ifut,1) !MG
               dy0=pcoord(ifut,2) !MG
               dz0=pcoord(ifut,3) !MG
               prefactor_ex=prefactor1*x0
               prefactor_coul=0.50d0*prefactor_ex*rho(ifut) !prefactor1 times point charge (x0 times rho)
               dist_x=xnew_coord(1)-dx0
               dist_y=xnew_coord(2)-dy0
               dist_z=xnew_coord(3)-dz0
               Rmod_CN2=dist_x**2.0d0+dist_y**2.0d0+dist_z**2.0d0 !R, in module SQUARED, between product gaussian (C) and point charge (N)
               if(Rmod_CN2.lt.1.0d-15) then !consider zero distance, TOL_DIST
                 Boys(n_max)=1.0d0/(2.0d0*real(n_max)+1.0d0)
                 xT=0.0d0
               else
                 xT=exp_sum*Rmod_CN2
                 if(n_max.eq.0) then !s-type
                   Boys(n_max)=DSQRT(PII/(4.0d0*xT))*erf(DSQRT(xT))
                 else if(xT.ge.8.0d0.and.xT.le.15.0d0) then !Taylor expansion
                   CALL Boys_expansion(n_max,xT,Boys(n_max))
                 else if(xT.gt.15.0d0) then !asymptotic value
                   istartingfactorial=2*n_max-1 !DATA doublefactorial/_,_,_...
                   idoublefactorial=istartingfactorial
                   ncount=1
                   do while(istartingfactorial-2*ncount.ge.2)
                     idoublefactorial=idoublefactorial*(istartingfactorial-2*ncount)
                     ncount=ncount+1
                   end do
                   xvalue=real(idoublefactorial)/2.0d0**(real(n_max)+1.0d0)
                   Boys(n_max)=xvalue*DSQRT(PII/(xT**(2.0d0*real(n_max)+1.0d0)))
                 else !infinite series
                   xa=real(n_max)+0.50d0
                   gammaincomplete=0.0d0
                   do i=0,30
                     xnewsum=xT**real(i)/((xa+real(i))*FACT(i))
                     if(MOD(i,2).ne.0) xnewsum=-xnewsum
                     gammaincomplete=gammaincomplete+xnewsum
                   end do
                   Boys(n_max)=0.50d0*gammaincomplete
                 end if
               end if
       
!!! DOWNWARD RECURRENCE FOR BOYS FUNCTION OF ALL ORDERS !!!
       
               k=1
               do while(k.le.n_max)
                 Boys(n_max-k)=(2.0d0*xT*Boys(n_max-k+1)+exp(-xT))/(2.0d0*(real(n_max-k+1))-1.0d0)
                 k=k+1
               end do
       
!!! GENERATING R VALUES FOR t=0,u=0,v=0 FOR N=[0,n_max] !!!
               xR_Boys=0.0d0 !probably not needed as it's recursive, checking
               do k=0,n_max
                 xR_Boys(k,0,0,0)=(-2.0d0*exp_sum)**k*Boys(k)
               end do
       
!!! MAIN RECURRENCE LOOP, TO OBTAIN R VALUES FOR ALL VALUES OF n,t,u,v !!!
!!! ALSO INCLUDES FINAL CALCULATIONS, USING ALL THE OBTAINED VALUES OF xR_Boys(0,t,u,v) AND E(imax,jmax,t/u/v) STORED AS E_coulomb(idir,it) !!!
       
               Vab=0.0d0
               do it=0,max_ang(1)
                 do iu=0,max_ang(2)
                   do iv=0,max_ang(3)
                     do k=0,n_max-it-iu-iv
                       if(it.eq.0) then
                         xR_Boys_a=0.0d0
                       else
                         xR_Boys_a=real(it)*xR_Boys(k+1,it-1,iu,iv)
                       end if
                       xR_Boys(k,it+1,iu,iv)=xR_Boys_a+dist_x*xR_Boys(k+1,it,iu,iv)
       
                       if(iu.eq.0) then
                         xR_Boys_a=0.0d0
                       else
                         xR_Boys_a=real(iu)*xR_Boys(k+1,it,iu-1,iv)
                       end if
                       xR_Boys(k,it,iu+1,iv)=xR_Boys_a+dist_y*xR_Boys(k+1,it,iu,iv)
       
                       if(iv.eq.0) then
                         xR_Boys_a=0.0d0
                       else
                         xR_Boys_a=real(iv)*xR_Boys(k+1,it,iu,iv-1)
                       end if
                       xR_Boys(k,it,iu,iv+1)=xR_Boys_a+dist_z*xR_Boys(k+1,it,iu,iv)
                     end do
                     Vab=Vab+E_coulomb(1,it)*E_coulomb(2,iu)*E_coulomb(3,iv)*xR_Boys(0,it,iu,iv)
                   end do
                 end do
               end do
c coulomb
               xresult_coul=prefactor_coul*Vab
               coul(icenter,icenter)=coul(icenter,icenter)+xresult_coul*p(iAO_a,iAO_b)

c exchange
               if(kop.ne.1) then !restricted exchange
                 xresult_ex=-prefactor_ex*Vab
                 do iMO=1,nocc
                   exchange0=xresult_ex*c(iAO_a,iMO) !x2 factor because closed shell calculation using MOs
                   do iMO2=1,nocc
                     exchange1=exchange0*c(iAO_b,iMO2)
                     exchange2=exchange1*chp2(ifut,iMO)*chp2(ifut,iMO2)
                     exch_hf(icenter,icenter)=exch_hf(icenter,icenter)+exchange2
                   end do
                 end do
               else !unrestricted exchange
                 xresult_ex=-prefactor_ex*Vab*0.5d0 !dividing by two because alpha and beta will be counted separately
               
                 do iMO=1,nalf
                   exchange0_alpha=xresult_ex*c(iAO_a,iMO)
                   if(iMO.le.nb) exchange0_beta=xresult_ex*cb(iAO_a,iMO)
                   do iMO2=1,nalf
                     exchange1_alpha=exchange0_alpha*c(iAO_b,iMO2)
                     exchange2_alpha=exchange1_alpha*chp2(ifut,iMO)*chp2(ifut,iMO2)
                     exch_hf(icenter,icenter)=exch_hf(icenter,icenter)+exchange2_alpha
                     if(iMO2.le.nb.and.iMO.le.nb) then
                       exchange1_beta=exchange0_beta*cb(iAO_b,iMO2)
                       exchange2_beta=exchange1_beta*chp2b(ifut,iMO)*chp2b(ifut,iMO2)
                       exch_hf(icenter,icenter)=exch_hf(icenter,icenter)+exchange2_beta
                     end if
                   end do
                 end do
               end if

             end do !closing loop over all integration points of atom icenter
           end do !closing loop over all atoms
         end do
       end do !end of main loop

       write(*,*) " "
       write(*,*) " *** ANALYTICAL COULOMB *** "
 
       coulen=ZERO
       do i=1,nat
         coulen=coulen+coul(i,i)
       end do
 
       write(*,*) " Ecoul TERMS "
       write(*,*) " "
       CALL Mprint(coul,nat,maxat)
       write(*,*) " Coulomb energy : ",coulen
 
       write(*,*) " "
       write(*,*) " *** ANALYTICAL EXCHANGE *** "
 
       exchen_hf=ZERO
       do i=1,nat
         exchen_hf=exchen_hf+exch_hf(i,i)
       end do

       if(kop.ne.1) then
         write(*,*) " HF Ex TERMS  "
       else
         write(*,*) " UHF Ex TERMS  "
       end if
       write(*,*) " "
       CALL Mprint(exch_hf,nat,maxat)
       write(*,*) " Hartree-Fock exchange energy : ",exchen_hf
       write(*,*) " "

      END SUBROUTINE

      SUBROUTINE Boys_expansion(n_max,xT,Boys)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BOYS_FACT(5)
      DATA BOYS_FACT/1.d0,2.d0,6.d0,24.d0,120.d0/
      DIMENSION WA_Boys(11,29)
      DATA (WA_Boys(n,1),n=1,11)/0.0195608d0,0.00364669d0,0.00111862d0,0.000468431d0,0.000242526d0,0.00014577d0,0.000097472d0,
     +0.0000704136d0,0.0000538481d0,0.0000429782d0,0.0000354424d0/,
     $     (WA_Boys(n,2),n=1,11)/0.0186829d0,0.00338106d0,0.00100873d0,0.000412112d0,0.000208955d0,0.000123469d0,0.000081445d0,
     +0.0000582071d0,0.0000441371d0,0.0000349906d0,0.0000286997d0/,
     $     (WA_Boys(n,3),n=1,11)/0.0178681d0,0.00314123d0,0.000911924d0,0.000363529d0,0.000180488d0,0.000104818d0,0.0000681859d0,
     +0.0000481953d0,0.0000362266d0,0.0000285198d0,0.0000232616d0/,
     $     (WA_Boys(n,4),n=1,11)/0.0171104d0,0.00292416d0,0.000826419d0,0.000321513d0,0.000156294d0,0.0000891873d0,0.0000571985d0,
     +0.0000399724d0,0.0000297754d0,0.0000232726d0,0.0000188722d0/,
     $     (WA_Boys(n,5),n=1,11)/0.0164044d0,0.00272721d0,0.000750702d0,0.000285084d0,0.000135686d0,0.0000760629d0,0.0000480782d0,
     +0.0000332091d0,0.000024508d0,0.0000190135d0,0.0000153263d0/,
     $     (WA_Boys(n,6),n=1,11)/0.0157453d0,0.0025481d0,0.000683481d0,0.000253419d0,0.00011809d0,0.0000650203d0,0.0000404947d0,
     +0.0000276383d0,0.0000202022d0,0.0000155529d0,0.0000124595d0/,
     $     (WA_Boys(n,7),n=1,11)/0.015129d0,0.00238485d0,0.000623653d0,0.000225827d0,0.000103031d0,0.0000557101d0,0.0000341779d0,
     +0.000023043d0,0.0000166778d0,0.0000127383d0,0.0000101396d0/,
     $     (WA_Boys(n,8),n=1,11)/0.0145517d0,0.00223574d0,0.000570276d0,0.000201725d0,0.0000901145d0,0.0000478443d0,0.0000289067d0,
     +0.0000192465d0,0.0000137895d0,0.0000104465d0,8.260577971693436E-6/,
     $     (WA_Boys(n,9),n=1,11)/0.0140101d0,0.00209924d0,0.000522541d0,0.000180619d0,0.0000790087d0,0.0000411848d0,0.0000245001d0,
     +0.0000161051d0,0.0000114193d0,8.578378262164428E-6,6.737300687156685E-6/,
     $     (WA_Boys(n,10),n=1,11)/0.0135012d0,0.00197405d0,0.000479752d0,0.000162093d0,0.000069438d0,0.0000355347d0,0.0000208094d0,
     +0.0000135017d0,9.471736878316576E-6,7.053926821346582E-6,5.501217677944795E-6/,
     $     (WA_Boys(n,11),n=1,11)/0.0130222d0,0.00185901d0,0.000441309d0,0.000145792d0,0.000061171d0,0.0000307307d0,0.0000177125d0,
     +0.0000113405d0,7.869163916680592E-6,5.80846024129448E-6,4.497200748449379E-6/,
     $     (WA_Boys(n,12),n=1,11)/0.0125709d0,0.00175308d0,0.000406696d0,0.000131415d0,0.0000540135d0,0.0000266374d0,0.0000151089d0,
     +9.543599724113801E-6,6.5486412555044375E-6,4.789710490139303E-6,3.680861022155174E-6/,
     $     (WA_Boys(n,13),n=1,11)/0.012145d0,0.00165538d0,0.000375463d0,0.000118706d0,0.0000478025d0,0.0000231421d0,0.0000129157d0,
     +8.04699180325652E-6,5.458961812050689E-6,3.955389710843769E-6,3.016431051707373E-6/,
     $     (WA_Boys(n,14),n=1,11)/0.0117426d0,0.0015651d0,0.000347222d0,0.000107447d0,0.0000424005d0,0.000020151d0,0.0000110647d0,
     +6.798376276547872E-6,4.558448846544275E-6,3.2712546857899363E-6,2.475068922112048E-6/,
     $     (WA_Boys(n,15),n=1,11)/0.0113619d0,0.00148155d0,0.000321635d0,0.0000974484d0,0.0000376915d0,0.0000175859d0,
     +9.499443638459046E-6,5.75485047731543E-6,3.81314628329268E-6,2.7095515558230534E-6,2.033499524941453E-6/,
     $     (WA_Boys(n,16),n=1,11)/0.0110013d0,0.00140409d0,0.000298406d0,0.0000885512d0,0.0000335775d0,0.0000153814d0,
     +8.173154150633697E-6,4.8811909545661805E-6,3.1953583574648872E-6,2.247765275090749E-6,1.672925359561938E-6/,
     $     (WA_Boys(n,17),n=1,11)/0.0106594d0,0.00133217d0,0.000277279d0,0.000080617d0,0.0000299754d0,0.0000134827d0,
     +7.04711984415125E-6,4.148441054539189E-6,2.6824702322432485E-6,1.8676134191372295E-6,1.378152893686431E-6/,
     $     (WA_Boys(n,18),n=1,11)/0.0103348d0,0.00126529d0,0.000258027d0,0.0000735268d0,0.0000268145d0,0.0000118439d0,
     +6.089189123641323E-6,3.532764059693503E-6,2.255994760108594E-6,1.5542360428544605E-6,1.1368914084822314E-6/,
     $     (WA_Boys(n,19),n=1,11)/0.0100264d0,0.00120301d0,0.000240454d0,0.000067178d0,0.000024035d0,0.0000104263d0,
     +5.272628628800348E-6,3.014511050397064E-6,1.9008013873868572E-6,1.2955429275308632E-6,9.391899322427788E-7/,
     $     (WA_Boys(n,20),n=1,11)/0.00973295d0,0.00114494d0,0.000224384d0,0.0000614818d0,0.0000215856d0,9.197635833852241E-6,
     +4.5751743306442665E-6,2.577462531412298E-6,1.6044918676611233E-6,1.0816872579180755E-6,7.76984784609612E-7/,
     $     (WA_Boys(n,21),n=1,11)/0.00945357d0,0.00109071d0,0.000209665d0,0.0000563613d0,0.0000194227d0,8.130379517397552E-6,
     +3.9782540122764285E-6,2.208210799121746E-6,1.3568943914649463E-6,9.046409242635738E-7,6.437357693257796E-7/,
     $     (WA_Boys(n,22),n=1,11)/0.0091873d0,0.00104001d0,0.000196161d0,0.0000517497d0,0.0000175089d0,7.201431182817885E-6,
     +3.4663493986594863E-6,1.8956564025560807E-6,1.1496533030678184E-6,7.578515640049953E-7,5.341334540326893E-7/,
     $     (WA_Boys(n,23),n=1,11)/0.0089333d0,0.000992538d0,0.000183753d0,0.0000475888d0,0.0000158121d0,6.391208339399082E-6,
     +3.026472197252E-6,1.6305971804591078E-6,9.758960363489164E-7,6.359653927498272E-7,4.438634874578624E-7/,
     $     (WA_Boys(n,24),n=1,11)/0.00869078d0,0.000948047d0,0.000172333d0,0.0000438278d0,0.0000143048d0,5.683105061880854E-6,
     +2.6477331561605937E-6,1.405392484813605E-6,8.299624811562363E-7,5.346030229794071E-7,3.694167080892861E-7/,
     $     (WA_Boys(n,25),n=1,11)/0.00845904d0,0.000906297d0,0.000161809d0,0.0000404225d0,0.0000129633d0,5.063013022410489E-6,
     +2.320987163294029E-6,1.2136885260823884E-6,7.071848651534656E-7,4.5017798995758116E-7,3.0793603821185376E-7/,
     $     (WA_Boys(n,26),n=1,11)/0.00823742d0,0.000867075d0,0.000152096d0,0.0000373341d0,0.000011767d0,4.5189232096371906E-6,
     +2.0385405792175185E-6,1.0501934551115288E-6,6.037085445372549E-7,3.7974972381135496E-7,2.5709294675278E-7/,
     $     (WA_Boys(n,27),n=1,11)/0.00802531d0,0.000830187d0,0.000143118d0,0.0000345284d0,0.0000106983d0,4.0405942296332715E-6,
     +1.7939095628505054E-6,9.104929579375757E-7,5.163459524955484E-7,3.2090432533957015E-7,2.1498769550217542E-7/,
     $     (WA_Boys(n,28),n=1,11)/0.00782215d0,0.00079546d0,0.00013481d0,0.0000319756d0,9.741947025074953E-6,3.6192756244500433E-6,
     +1.5816202292644905E-6,7.908988842176024E-7,4.424574466830895E-7,2.716578010993421E-7,1.8006872774705342E-7/,
     $     (WA_Boys(n,29),n=1,11)/0.00762742d0,0.000762731d0,0.000127112d0,0.0000296492d0,8.884563981972003E-6,3.24747671603967E-6,
     +1.3970431662671294E-6,6.883248391168378E-7,3.798539981494806E-7,2.303774548112768E-7,1.510674743511661E-7/

       xcheck_min=0.25d0
       do iT=1,29 !find which pretabulated xT is closest to the xT of Boys being calculated
         xpoint=(8.d0+0.25d0*real(iT-1))
         xcheck=ABS(xpoint-xT)
         if(xcheck.lt.xcheck_min) then
           xcheck_min=xcheck
           iT_min=iT
         end if
       end do

       Boys=WA_Boys(n_max,iT_min)
       xdeltaT=xT-(8.d0+0.25d0*real(iT_min-1))
       do i=1,5 !Taylor expand around the found closest value to obtain the Boys function of interest
         Taylor_term=WA_Boys(n_max+i,iT_min)*xdeltaT**real(i)/BOYS_FACT(i)
         if(MOD(i,2).ne.0) Taylor_term=-Taylor_term
         Boys=Boys+Taylor_term
       end do

      END SUBROUTINE

      subroutine numint_dft_analytical(ndim,itotps,wp,omp,omp2,chp,eto,pcoord,sat)
      use xc_f90_types_m
      use xc_f90_lib_m
      USE ao_matrices
      USE integration_grid
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(100)
      common/actual/iact,jat,icenter
      common /achi/achi(maxat,maxat),ibcp
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/exchg/exch(maxat,maxat),xmix
      dimension eto(maxat,maxat)
      dimension wp(itotps),chp(itotps,ndim)
      dimension omp(itotps),omp2(itotps,nat),pcoord(itotps,3)
! (TO DO) Laplacian!!!
      dimension sat(igr,igr,nat)
      character*80 line

      allocatable :: chp2(:,:),scr(:),rho(:),scrx(:,:)
      allocatable :: scr_a(:),rho_a(:),scr2(:)
      allocatable :: sab(:,:,:)

      i5d      =  Iopt(11)
      idofr    =  Iopt(40)
      ithrebod =  Iopt(44)
      itype    =  Iopt(55)
      ianalytical=  Iopt(90)
      iatps=nang*nrad

      ALLOCATE(chp2(itotps,nocc),scr(itotps),rho(itotps))
      ALLOCATE(scr_a(itotps),rho_a(itotps),scr2(itotps))
      ALLOCATE(sab(nocc,nocc,nat))
      ALLOCATE(scrx(igr,nocc))

!! BUILDING RHO!!

      do kk=1,itotps
        xx0=ZERO
        do jj=1,nocc
          xx=ZERO
          do ii=1,igr 
            xx=xx+c(ii,jj)*chp(kk,ii) 
          end do
          chp2(kk,jj)=xx
          xx0=xx0+xx*xx
        end do
        rho(kk)=TWO*xx0
      end do

!! EXCHANGE AND CORRELATION FOR DFT !!

      write(*,*) " *** DFT MOLECULAR ENERGY DECOMPOSITION *** "
      write(*,*) " "

      do iatom=1,nat
        do nu=1,igr
          do ii=1,nocc
            xx=ZERO
            do mu=1,igr
              xx=xx+c(mu,ii)*sat(mu,nu,iatom)
            end do
            scrx(nu,ii)=xx
          end do
         end do
  
         do ii=1,nocc
          do jj=1,nocc
            xx=ZERO
            do nu=1,igr
              xx=xx+c(nu,jj)*scrx(nu,ii)        
            end do
            sab(ii,jj,iatom)=xx
            if(ii.ne.jj) sab(jj,ii,iatom)=xx
          end do
        end do
      end do

!! CHECKING !!

      write(*,*) " CHECKING MO OVERLAPS "
      do ii=1,nocc                 
        do jj=ii,nocc
          x=ZERO
          if(ii.eq.jj) x=-ONE
          do kk=1,nat
            x=x+sab(ii,jj,kk)
          end do
!         if(abs(x).gt.1.0d-2) then
          if(abs(x).gt.1.0d-3) then
            write(*,*) ii,jj,x
            stop " Problem with MOs overlaps "
          end if
        end do
      end do
      write(*,*) " "

!! chp2 : MOs, chpd : MO DERIVATIVES !!

!      write(*,*) 'Checking itype',itype !B3LYP is 15
!      if(itype.gt.1) call grdrho(itotps,nocc,chpd,0) !program works with this line commented... removed allocate

      if(itype.gt.1) then
        ideriv=0
        call sigma(pcoord,chp2,scr) !this line IS NEEDED for analytical as well
      end if

!! CALCULATE "EXACT" (ONCE-CENTER) XC ENERGY !! 
!! Laplacian for the MGGA functionals can be calculated using LAPLACIAN subroutine (qtaim.f). NOT IMPLEMENTED

      call xc(itotps,rho,scr,scr2) 
      xtot=ZERO
      do icenter=1,nat
        x=ZERO
!! ALLPOINTS INTEGRATION ACTIVATED... !!
        do ifut=1,itotps 
          x=x+wp(ifut)*scr2(ifut)*omp(ifut)*omp2(ifut,icenter)
        end do
        exch(icenter,icenter)=x
        xtot=xtot+x
      end do
      write(*,*) " PURE DFT ONE-CENTER CONTRIBUTIONS (EXACT) "
      write(*,*) " "
      CALL Mprint(exch,NAT,maxat)
      write(*,*) " Exact one-center decomposition for DFT xc energy : ",xtot
      write(*,*) " "

!! PRINTING !!

      exchen=0.0d0
      xtot=ZERO
      do ii=1,nat
        exchen=exchen+exch(ii,ii)
        eto(ii,ii)=eto(ii,ii)+exch(ii,ii)
        xtot=xtot+eto(ii,ii)
      end do
      etot=xtot
      write(*,*) " FINAL KS-DFT XC ENERGY COMPONENTS "
      write(*,*) " "
      call Mprint(exch,nat,maxat)
      write(*,*) " Sum of (exact) DFT xc energy components : ",exchen
      write(*,*) " "
      if(xmix.gt.ZERO) then
        write(*,*) " WARNING: HF-exchange part missing "
        write(*,*) " "
      end if
      if (idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Exc Decomposition' 
        call group_by_frag_mat(1,line,exch)
      end if

      DEALLOCATE(chp2,scr,rho,scr_a,rho_a,scr2,sab,scrx)
      end

      SUBROUTINE Population_analysis(aux_sat,nat0,igr0,aux_qat)
      USE basis_set
      USE ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
!      common /ps/ps(nmax,nmax),pa(nmax,nmax),pb(nmax,nmax)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /iops/iopt(100)
      DIMENSION aux_sat(igr0,igr0,nat0)
      DIMENSION aux_qat(maxat,2)
      character*80 line

      idoint = iopt(1)
      ipolar = iopt(46)

      if(ipolar.eq.1) iaccur=1

      if(kop.ne.0.or.nalf.ne.nb) then
        tc1=0.d0
        tc2=0.d0
        do kat=1,nat
          x=0.d0
          do mu=1,igr
            do nu=1,igr
              x=x+ps(mu,nu)*aux_sat(mu,nu,kat)
            end do
          end do
          qsat(kat,1)=x
          qsat(kat,2)=0.d0
        end do
        do mu=1,igr
          kat=ihold(mu)
          x=0.d0
          do itau=1,igr
            x=x+ps(mu,itau)*s(itau,mu)
          end do
          qsat(kat,2)=qsat(kat,2)+x
        end do
      end if

c Total density 
      tc1=0.d0
      tc2=0.d0
      do kat=1,nat
        x=0.d0
        do mu=1,igr
          do nu=1,igr
            x=x+p(mu,nu)*aux_sat(mu,nu,kat)
          end do
        end do
        qat(kat,1)=x
        qat(kat,2)=0.d0
      end do
      do mu=1,igr
        kat=ihold(mu)
        x=0.d0
        do itau=1,igr
          x=x+p(mu,itau)*s(itau,mu)
        end do
        qat(kat,2)=qat(kat,2)+x
      end do

      if(idoint.eq.1)  call print_int(ndim,nat,sat,namepat)

CCCCCCCCCCCCCCC
C ELECTRON POPULATIONS
CCCCCCCCCCCCCCC
      tc1=0.d0
      tc2=0.d0
      do i=1,nat
       tc1=tc1+qat(i,1)
       tc2=tc2+qat(i,2)
       aux_qat(i,1)=qat(i,1)
       aux_qat(i,2)=qat(i,2)
      enddo
      print *,'  '
      print *,'    ELECTRON POPULATIONS'
      print *,'  '
      print *,'  Atom   apost3d      Mulliken'
      print *,' -----------------------------'
      call vprint(qat,nat,maxat,2)
      print *,' -----------------------------'
      if(iaccur.eq.0) then
        write(*,'(1x,a,2(F10.6,2X))') '    Sum  ', tc1, tc2
      else
        write(*,'(1x,a,2(F20.13,2X))') '    Sum  ', tc1, tc2
      end if

      if (idofr.eq.1) then
       line ='   FRAGMENT ANALYSIS : Electron populations'
       call group_by_frag_vec(2,line ,qat)
      end if

      END SUBROUTINE
