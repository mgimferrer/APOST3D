
!! ***** !!

      subroutine numint_one(ndim,itotps,wp,rho,omp,omp2,pcoord,chp,eto)

      use ao_matrices
      use integration_grid

      implicit real*8(a-h,o-z)

      include 'parameter.h'

      common /nat/nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common /actual/iact,jat,icenter
      common /energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common /energ0/ekin0,eelnuc0,evee0,etot0
      common /exchg/exch(maxat,maxat),xmix
      common /efield/field(4),edipole

      integer, intent(in) :: ndim,itotps
      character*80 line

      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3)
      dimension omp(itotps),omp2(itotps,nat),rho(itotps)
! TO DO: xkdens = kinetic energy density
! TO DO: Laplacian!
      dimension xkdens(itotps)
      dimension eto(maxat,maxat)
      dimension epa(maxat,maxat),ekin(maxat,maxat),dip(maxat,3),enuc(maxat,maxat)
      dimension eecp(maxat),ect(maxat,maxat)
      dimension edip(maxat,maxat)

      allocatable :: chp2(:,:),scr(:),chpaux(:,:)

      idofr   = iopt(40)
      ifield  = iopt(47)
      ipoints = itotps
      iatps   = nang*nrad
      xxdip   = ZERO

      write(*,*) " "
      write(*,*) " --------------------------- "
      write(*,*) "  GENERAL ONE-ELECTRON PART  "
      write(*,*) " --------------------------- "
      write(*,*) " "

      ALLOCATE(chp2(itotps,nocc),scr(itotps))

!! TO MOs !!
      do k=1,ipoints
        xx0=ZERO
        do j=1,nocc
          xx=ZERO
          do i=1,igr 
            xx=xx+c(i,j)*chp(k,i) 
          end do
          chp2(k,j)=xx
        end do
      end do

!! ZEROING MATRICES !!
      epa=ZERO
      ekin=ZERO
      enuc=ZERO
      eto=ZERO
      ect=ZERO

!! DOING ENERGY PARTITION !!
!! ONE-ELECTRON PART !!

!! IN CASE OF HAVING PSEUDOPOTENTIAL !!
      call mga_misc(iecp,eecp) !! MAYBE WE SHOULD DO A BETTER VERSION OF THIS... !!
      if(iecp.eq.1) then
        write(*,*) " ADDING ECP ATOMIC ENERGIES TO E-N TERMS "
        do ii=1,nat
          epa(ii,ii)=eecp(ii)
        end do 
      end if

!! IN CASE OF HAVING AN ELECTRIC FIELD !!
      !MMO- minor: is it really necessary to call field_misc again? this has been done on main.f, 
      !MMO- and everything is stored in /common/efield
      call field_misc(ifield)
      !MMO- fixed the following: xxdip should be initialized as 0.0d0 if there's no field... 
      !MMO- as xxdip should be used for the printing of integration error (see below)
      if(ifield.eq.1) then
        write(*,*) " The system is under a static electric field "
        write(*,'(2x,3(x,a3,x,f8.6))') "Fx=",field(2),"Fy=",field(3),"Fz=",field(4)
        write(*,*) " Electron-nuclear and nuclear repulsion terms are affected by E. field "
        xtot=ZERO
        do icenter=1,nat
          x=ZERO
          y=ZERO
          z=ZERO
          zztop=ZERO
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            distx=pcoord(ifut,1)
            disty=pcoord(ifut,2)
            distz=pcoord(ifut,3)
            x=x+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*distx
            y=y+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*disty
            z=z+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*distz
            zztop=zztop+wp(ifut)*omp2(ifut,icenter)*rho(ifut)
          end do
          dip(icenter,1)=-x+zztop*coord(1,icenter)    
          dip(icenter,2)=-y+zztop*coord(2,icenter)    
          dip(icenter,3)=-z+zztop*coord(3,icenter)    
          xtot=(coord(1,icenter)*field(2)+coord(2,icenter)*field(3)+coord(3,icenter)*field(4))    
          ect(icenter,icenter)=-zztop*xtot                
        end do

!! ENERGY CONTRIBUTIONS. ELECTRIC FIELD X,Y,Z COMPONENTS ARE ON field(2-4) !!
        xtot=ZERO
        do i=1,nat
          edip(i,i)=+dip(i,1)*field(2)+dip(i,2)*field(3)+dip(i,3)*field(4)
          xtot=xtot+ect(i,i)+edip(i,i)
        end do
        xxdip=xtot
        write(*,'(2x,a37,x,f14.7)') "Electronic dipole moment energy term:",xxdip
        write(*,*) " "
      else
        !! MG: at the beginning of the subroutine is made zero, necessary again? !!
        xxdip=ZERO !MMO- added this else condition
      end if

!! ELECTRON-NUCLEI ATTRACTION INTEGRALS !!
      xtot=ZERO
      do icenter=1,nat
        do jcenter=1,nat     
          x=ZERO
          zz=zn(jcenter)
          dx0=coord(1,jcenter)
          dy0=coord(2,jcenter)
          dz0=coord(3,jcenter)
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            distx=pcoord(ifut,1)-dx0
            disty=pcoord(ifut,2)-dy0
            distz=pcoord(ifut,3)-dz0
            rr=dsqrt(distx*distx+disty*disty+distz*distz)
            if(rr.gt.1.0d-10) then
              x=x+zz*wp(ifut)*omp2(ifut,icenter)*rho(ifut)/rr
            end if
          end do
          xtot=xtot+x
          epa(icenter,jcenter)=epa(icenter,jcenter)-x
        end do
      end do

      xtot=ZERO
      do i=1,nat
        xtot=xtot+epa(i,i)
        do j=1,i-1  
          epa(i,j)=(epa(i,j)+epa(j,i))
          epa(j,i)=epa(i,j)
          xtot=xtot+epa(i,j)
        end do
      end do
      eelnuc=xtot

!! PRINTING MATRIX !!
      write(*,*) " ----------------------------- "
      write(*,*) "  ELECTRON-NUCLEAR ATTRACTION  "
      write(*,*) " ----------------------------- "
      write(*,*) " "
      call MPRINT2(epa,nat,maxat)
      write(*,'(2x,a23,x,f14.7)') "Electron-nuclei energy:",eelnuc

!! CHECKING INTEGRATION ERROR (IF ENERGIES INCLUDED IN .fchk FILE) !!
      if(eelnuc0.ne.ZERO) then
        !MMO- printing fix: to correctly evaluate error, add xxdip here, because during the decomp process we do not add 
        !MMO- it to any atomic contributions (we avoid atomic terms related to dipole)
        write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",(eelnuc+xxdip-eelnuc0)*tokcal
      else
        eelnuc0=eelnuc
      end if
      write(*,*) " "
      if(idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Electron-nuclei attraction ' 
        call group_by_frag_mat(1,line,epa)
      end if

!! GENERATIONG GRID FOR 2nd DERIVATIVE OVER AOs !!
      ALLOCATE(chpaux(itotps,igr))
      call dpoints(chpaux,pcoord)

!! TO MOs AND MULTIPLY BY MOs AND SUM OVER MOs USING SCRATCH ARRAY !!
!! HENCE WE HAVE DENSITY IN scr EXCEPT FOR THE FACTOR OF TWO !!
      do k=1,ipoints
        x=ZERO
        do j=1,nocc
          xx=ZERO
          do i=1,igr 
            xx=xx+c(i,j)*chpaux(k,i) 
          end do
          x=x+chp2(k,j)*xx
        end do
        scr(k)=x
!! TO DO: Kinetic energy density !!
!       xkdens(k)=scr(k)
!! TO DO: Laplacian !!
      end do

      xtot=ZERO
      do icenter=1,nat
        x=ZERO
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          x=x+wp(ifut)*scr(ifut)*omp2(ifut,icenter)
        end do
        ekin(icenter,icenter)=ekin(icenter,icenter)+x
        xtot=xtot+x
      end do
      ekinen=xtot

!! PRINTING MATRIX !!
      write(*,*) " ---------------- "
      write(*,*) "  KINETIC ENERGY  "
      write(*,*) " ---------------- "
      write(*,*) " "
      call MPRINT2(ekin,nat,maxat)
      write(*,'(2x,a15,x,f14.7)') "Kinetic energy:",ekinen

!! CHECKING INTEGRATION ERROR (IF ENERGIES INCLUDED IN .fchk FILE) !!
      if(ekin0.ne.ZERO) then
        write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",(ekinen-ekin0)*tokcal
      else
        ekin0=ekinen  
      end if
      write(*,*) " "
      if (idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Kinetic energy ' 
        call group_by_frag_mat(1,line,ekin)
      end if
      DEALLOCATE(chpaux)

!! ADDING NUCLEAR REPULSION !!
      erep=ZERO
      do i=1,nat
        eto(i,i)=ekin(i,i)+epa(i,i)
        do j=i+1,nat
          dist=dsqrt((coord(1,i)-coord(1,j))**TWO+(coord(2,i)-coord(2,j))**TWO+(coord(3,i)-coord(3,j))**TWO)
          eto(i,j)=ekin(i,j)+epa(i,j)+zn(i)*zn(j)/dist
          eto(j,i)=eto(i,j)

!! SAVING ALSO FOR PRINTING PURPOSES, I THINK WE SHOULD GIVE IT... COMES FOR FREE !!
          enuc(i,j)=zn(i)*zn(j)/dist
          enuc(j,i)=enuc(i,j)
          erep=erep+zn(i)*zn(j)/dist
        end do
      end do

!! PRINTING MATRIX !!
      write(*,*) " --------------------------- "
      write(*,*) "  NUCLEAR-NUCLEAR REPULSION  "
      write(*,*) " --------------------------- "
      write(*,*) " "
      call MPRINT2(enuc,nat,maxat)
      write(*,'(2x,a25,x,f14.7)') "Nuclear repulsion energy:",erep
      write(*,*) " "

!! IN CASE OF HAVING AN ELECTRIC FIELD !!
      if(ifield.eq.1) then
        xxx=ZERO
        xtot=ZERO
        do i=1,nat
          xxx=zn(i)*(coord(1,i)*field(2)+coord(2,i)*field(3)+coord(3,i)*field(4))
          ect(i,i)=ect(i,i)+xxx
          xtot=xtot+xxx
        end do
        write(*,'(2x,a50,x,f14.7)') "Nuclear repulsion energy including nuclear dipole:",erep+xtot
        write(*,'(2x,a34,x,f14.7)') "Nuclear dipole moment energy term:",xtot
        edipole=xtot+xxdip
        write(*,'(2x,a32,x,f14.7)') "Total dipole moment energy term:",edipole
        write(*,*) " "
        write(*,*) " ----------------------------------------------------------- "
        write(*,*) "  INTRINSIC DIPOLE ENERGY CONTRIBUTION (ORIGIN INDEPENDENT)  "
        write(*,*) " ----------------------------------------------------------- "
        write(*,*) " "
        call MPRINT2(edip,nat,maxat)
        xtot=ZERO
        do i=1,nat
          xtot=xtot+edip(i,i)
        end do
        write(*,'(2x,a28,x,f14.7)') "Total Intrinsic dipole term:",xtot
        write(*,*) " "
        write(*,*) " -------------------------------------------------------- "
        write(*,*) "  CHARGE-TRANSFER ENERGY CONTRIBUTION (ORIGIN DEPENDENT)  "
        write(*,*) " -------------------------------------------------------- "
        write(*,*) " "
        call MPRINT2(ect,nat,maxat)
        xtot=ZERO
        do i=1,nat
          xtot=xtot+ect(i,i)
        end do
        write(*,'(2x,a27,x,f14.7)') "Total Charge-Transfer term:",xtot
        write(*,*) " "
      end if

!! TOTAL ENERGY UP TO HERE !!
      etot=ekinen+eelnuc+erep
      etot0=ekin0+eelnuc0+erep

!! DEALLOCATING !!
      DEALLOCATE(chp2,scr)

      end

!! ***** !!

      subroutine numint_two(ndim,itotps,wp,omp,omp2,pcoord,chp,rho,eto)

      use ao_matrices
      use integration_grid
      use OMP_LIB
      use IFPORT

      implicit real*8(a-h,o-z)

      integer external OMP_GET_NUM_PROCS

      include 'parameter.h'
      parameter (maxpass=5)

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /iops/iopt(100)
      common/actual/iact,jat,icenter
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/energ0/ekin0,eelnuc0,evee0,etot0
      common/exchg/exch(maxat,maxat),xmix
      common/twoel/twoeltoler
      common/efield/field(4),edipole
!! FROM # GRID SECTION !!
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22
      common /modgrid2/thr3

      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3)
      dimension omp(itotps),omp2(itotps,nat),rho(itotps)
      dimension eto(maxat,maxat),Excmp(maxat,maxat)

      dimension coul0(maxat,2*maxpass)
      dimension exch_hf(maxat,maxat)
      dimension coul(maxat,maxat)
      
      character*80 line
      character*100 threadenv

      allocatable :: chp2(:,:), fij(:,:),xocc(:,:)
      allocatable :: chppha(:,:), wppha(:),pcoordpha(:,:),omppha(:)
      allocatable :: chp2pha(:,:),rhopha(:),omp2pha(:,:),ibaspointpha(:)

!! FOR ENPART PARALLEL !!
      allocatable :: chp2s(:,:),chp2phas(:,:) !! SWITCHED ORDER OF COLUMNS AND ROWS !!
      allocatable :: ijpaircount(:,:),istart(:),iend(:) !! ATOM PAIRS INCLUDED, AND FOR SLICING chp MATRICES !!
      allocatable :: exch_hfk(:,:)
      allocatable :: exch_hfij(:,:)
      allocatable :: f3k(:)

      ihf         = Iopt(18) 
      idofr       = Iopt(40)
      ithrebod    = Iopt(44)
      ifield      = Iopt(47)
      ipolar      = Iopt(46)
      ianalytical = Iopt(91)

      !thr2=thr3 !! THRESH FOR NUMERICAL INTEGRATION (DISTANCE) !!
      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if
      idoex=0
      if(xmix.gt.ZERO) idoex=1
      iatps=nang*nrad

!! DOING ENERGY PARTITION !!
      npass=2     
      write(*,*) " "
      write(*,*) " --------------------------- "
      write(*,*) "  GENERAL TWO-ELECTRON PART  "
      write(*,*) " --------------------------- "
      write(*,*) " "
      do i=1,nat
        do j=1,nat
          coul(i,j)=ZERO
          exch_hf(i,j)=ZERO
        end do
        do j=1,2*npass
          coul0(i,j)=ZERO
        end do
      end do

      ALLOCATE(chp2(itotps,nocc))
      ALLOCATE(chp2s(nocc,itotps)) !! TODO: CHANGING THAT NEVER USE CHP2 !!

!! TO MOs !!
      do k=1,itotps
        do j=1,nocc
          xx=ZERO
          do i=1,igr 
            xx=xx+c(i,j)*chp(k,i) 
          end do
          chp2(k,j)=xx
          chp2s(j,k)=xx
        end do
      end do

      if(ianalytical.eq.0) then !MMO- skip if analytical

!! ROTATED GRID, ANGLE CONTROLLED BY # GRID SECTION !!
!! CHECK diatXC PAPER FOR OPTIMIZED VALUES: phb=0.162d0 and later 0.182d0 for 40 146 !!
      phb=phb12
      pha=ZERO 
      write(*,'(2x,a20,x,f10.6,x,f10.6)') "Rotating for angles:",pha,phb
      ALLOCATE(wppha(itotps),omppha(itotps),omp2pha(itotps,nat))
      ALLOCATE(chppha(itotps,ndim),pcoordpha(itotps,3),ibaspointpha(itotps))
      ALLOCATE(chp2pha(itotps,nocc),rhopha(itotps))
      ALLOCATE(chp2phas(nocc,itotps)) !!TODO

      call prenumint(ndim,itotps,nat,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)

!! TO MOs !!
      do k=1,itotps
        do j=1,nocc
          xx=ZERO
          do i=1,igr 
            xx=xx+c(i,j)*chppha(k,i) 
          end do
          chp2pha(k,j)=xx !! MG: TODO
          chp2phas(j,k)=xx
        end do
      end do
      end if !MMO- non-analytical skip ends here

!! COULOMB PART !!
      call cpu_time(xtime)
      if(ianalytical.eq.0) then
        call calc_coul(itotps,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,rho,rhopha,coul)
      else
        call calc_twoel_analytical(itotps,wp,omp2,pcoord,rho,chp2,chp2b,coul,exch_hf)
      end if
      call cpu_time(xtime2)
      write(*,'(a32,f10.1,a2)')'(Elapsed time :: Coulomb energy ',xtime2-xtime,'s)' 
      write(*,*) " "
      xtime=xtime2

!! IN CASE OF HAVING HF-TYPE EXCHANGE !!
      if(ianalytical.eq.0) then !MMO- skip if analytical
      if(idoex.eq.1) then

!! COMPUTING MULTIPOLAR APPROACH HERE !!
        norb2=nocc*(nocc+1)/2
        ALLOCATE(fij(itotps,norb2),xocc(norb2,norb2))
        do ii=1,itotps
          irun=0
          do jj=1,nocc
            do kk=jj,nocc
              irun=irun+1
              fij(ii,irun)=chp2(ii,jj)*chp2(ii,kk)
              if(ii.eq.1) then
                xocc(irun,irun)=-TWO
                if(jj.ne.kk) xocc(irun,irun)=-FOUR
              end if
            end do
          end do
        end do
        call multipolar(nocc,itotps,wp,omp2,pcoord,fij,xocc,Excmp)
        DEALLOCATE(fij,xocc)

!! MG: REORDERING THE LOOPS FOR PARALLELIZATION PURPOSES !!
!! IMPLEMENTATION PERFORMED THANKS TO Dr. R. OSWALD !!
        write(*,*) " ------------------------------------------------- "
        write(*,*) "  EVALUATING HARTREE-FOCK-TYPE EXCHANGE INTEGRALS  "
        write(*,*) " ------------------------------------------------- "
        write(*,*) " "
        write(*,'(2x,a22,x,i6,x,a8)') "Two-el integrations for",nocc*(nocc+1),"MO pairs"
        write(*,'(2x,a36,x,f10.6)') "Threshold for atom pair calculation :",threbod

!! FIRST ONLY SAME CENTER TERMS !!
!! EVALUATING NUMBER OF CORES FOR SPLITTING THE CALCULATION BY THREADS !!
        call getenv('OMP_NUM_THREADS',threadenv)
        if(trim(threadenv)=='') then
          write(*,*) " OMP_NUM_THREADS not set"
          ithreadenv=ZERO
        else
          read(unit=threadenv,FMT='(I4)') ithreadenv
        end if
        !$OMP PARALLEL
        iprocs=OMP_GET_MAX_THREADS()
        ithreads=INT(OMP_GET_NUM_PROCS())
        !$OMP END PARALLEL
        if(ithreadenv.ne.ZERO) then
          write(*,'(2x,a43,x,i3,x,a14,x,i3,x,a24)') "Two-el integration will be distributed over",ithreadenv,"threads out of",
     &    ithreads,"available hardware cores"
          ithreads=ithreadenv
        else
          write(*,'(2x,a43,x,i3,x,a14,x,i3,x,a24)') "Two-el integration will be distributed over",ithreads,"threads out of",
     &    ithreads,"available hardware cores"
        end if

!! TO ENSURE PROPER SLICING BY THREADS !!
        ALLOCATE(f3k(ithreads))
        ALLOCATE(exch_hfk(nat,ithreads))
        ALLOCATE(istart(ithreads),iend(ithreads))
        itilerest=mod(iatps,ithreads)
        ispace=(iatps-itilerest)/ithreads

!! MORE CONVOLUTED LOOP STRUCTURE, AVOIDED PROBLEM OF MAX PARALLEL 8 CORES !!
        exch_hfk=ZERO
        !DIR$ NOPARALLEL
        do icenter=1,nat
          ioffset=(icenter-1)*iatps
          istart=0
          iend=0
          !DIR$ NOPARALLEL
          do ik=1,(ithreads-1)
            istart(ik)=ioffset+((ik-1)*ispace)+1
            iend(ik)=ioffset+(ik*ispace)
          end do
          istart(ithreads)=iend(ithreads-1)+1
          iend(ithreads)=icenter*iatps

!! THIS TWO CALLS ARE CRUCIAL !!
          call omp_set_dynamic(.false.)
          call omp_set_num_threads(ithreads)
          !DIR$ PARALLEL
          do ik=1,ithreads
            do ifut=istart(ik),iend(ik)
              x0=wp(ifut)*omp2(ifut,icenter)
              dx0=pcoord(ifut,1)
              dy0=pcoord(ifut,2)
              dz0=pcoord(ifut,3)
              f3k(ik)=ZERO
              do jfut=iatps*(icenter-1)+1,iatps*icenter
                x1=wppha(jfut)*omp2pha(jfut,icenter)
                dx1=pcoordpha(jfut,1)
                dy1=pcoordpha(jfut,2)
                dz1=pcoordpha(jfut,3)
                dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                if(dist.gt.1.0d-12) then !! MG: Could be controlled using the thr2 variable !!
                  do i=1,nocc-1
                    f2=x0*chp2s(i,ifut)*chp2s(i,ifut)
                    f3k(ik)=f3k(ik)+f2*x1*chp2phas(i,jfut)*chp2phas(i,jfut)/dist
                    do j=i+1,nocc
                      f2=TWO*x0*chp2s(i,ifut)*chp2s(j,ifut)
                      f3k(ik)=f3k(ik)+f2*x1*chp2phas(i,jfut)*chp2phas(j,jfut)/dist !! MG: LC-wPBE programmed in version 3.1 !!
                    end do
                  end do
                  f2=x0*chp2s(nocc,ifut)*chp2s(nocc,ifut)
                  f3k(ik)=f3k(ik)+f2*x1*chp2phas(nocc,jfut)*chp2phas(nocc,jfut)/dist
                end if
              end do
              exch_hfk(icenter,ik)=exch_hfk(icenter,ik)+f3k(ik)
            end do
          end do
        end do

!! ADDING THE TERMS INTO THE ORIGINAL exch_hf MATRIX !!
        do icenter=1,nat
          do isum=1,ithreads
             exch_hf(icenter,icenter)=exch_hf(icenter,icenter)-exch_hfk(icenter,isum)
          end do
        end do

!! NOW PAIRS OF CENTERS !!
        iterms=0
        ipaircounter=0
        ALLOCATE(ijpaircount(nat*nat,2))
        do icenter=1,nat
          do jcenter=icenter+1,nat
            bx0=bo(icenter,jcenter)
            if(bx0.ge.threbod) then
              ipaircounter=ipaircounter+1
              ijpaircount(ipaircounter,1)=icenter
              ijpaircount(ipaircounter,2)=jcenter
            else
              iterms=iterms+1
            end if
          end do
        end do
        write(*,'(2x,a34,x,i5,x,a10)') "Skipping numerical integration for",iterms,"atom pairs"
        write(*,*) " "

!! AGAIN, PREPARING FOR SLICING AND BLOCKING PARALLELIZATION OF SOME LOOPS !!
        ALLOCATE(exch_hfij(ipaircounter,ithreads))
        itilerest=mod(iatps,ithreads)
        ispace=(iatps-itilerest)/ithreads
        exch_hfij=ZERO
        !DIR$ NOPARALLEL
        do numpairnat=1,ipaircounter
          icenter=ijpaircount(numpairnat,1)
          jcenter=ijpaircount(numpairnat,2)
          ioffset=(icenter-1)*iatps
          istart=0
          iend=0
          !DIR$ NOPARALLEL
          do ik=1,(ithreads-1)
            istart(ik)=ioffset+((ik-1)*ispace)+1
            iend(ik)=ioffset+(ik*ispace)
          end do
          istart(ithreads)=iend(ithreads-1)+1
          iend(ithreads)=icenter*iatps
          !DIR$ PARALLEL
          do ik=1,ithreads
            do ifut=istart(ik),iend(ik)
              x0=wp(ifut)*omp2(ifut,icenter)
              dx0=pcoord(ifut,1)
              dy0=pcoord(ifut,2)
              dz0=pcoord(ifut,3)
              chp2snifut=chp2s(nocc,ifut)
              f3k(ik)=ZERO
              do jfut=iatps*(jcenter-1)+1,iatps*jcenter
                x1=wp(jfut)*omp2(jfut,jcenter)
                dx1=pcoord(jfut,1)
                dy1=pcoord(jfut,2)
                dz1=pcoord(jfut,3)
                dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                if(dist.gt.1.0d-12) then !! Could be controlled using the thr2 variable !!
                  chp2snjfut=chp2s(nocc,jfut)
                  do i=1,nocc-1
                    chp2sijfut=chp2s(i,jfut)
                    chp2siifut=chp2s(i,ifut)
                    do j=i+1,nocc
                      f2=TWO*x0*chp2siifut*chp2s(j,ifut)
                      f3k(ik)=f3k(ik)+f2*x1*chp2sijfut*chp2s(j,jfut)/dist !! LC-wPBE programmed in version 3.1 !!
                    end do
                    f2=x0*chp2siifut*chp2siifut
                    f3k(ik)=f3k(ik)+f2*x1*chp2sijfut*chp2sijfut/dist
                  end do
                  f2=x0*chp2snifut*chp2snifut
                  f3k(ik)=f3k(ik)+f2*x1*chp2snjfut*chp2snjfut/dist
                end if
              end do
              exch_hfij(numpairnat,ik)=exch_hfij(numpairnat,ik)+f3k(ik)
            end do
          end do
        end do

!! AGAIN, ADDING THE TERMS INTO THE ORIGINAL exch_hf MATRIX !!
        do numpairnat=1,ipaircounter
          icenter=ijpaircount(numpairnat,1)
          jcenter=ijpaircount(numpairnat,2)
          do isum=1,ithreads
            exch_hf(icenter,jcenter)=exch_hf(icenter,jcenter)-exch_hfij(numpairnat,isum)
          end do
        end do

!! ACCOUNTING THAT AB = BA (FACTOR OF 2), AND IF SOME TERMS COME FROM MULTIPOLAR APPROACH !!
        exchen_hf=ZERO
        do i=1,nat
          exchen_hf=exchen_hf+exch_hf(i,i)
          do j=i+1,nat
            if(bo(i,j).lt.threbod) then 
              exch_hf(i,j)=Excmp(i,j)
            else 
              exch_hf(i,j)=TWO*(exch_hf(i,j))
            end if 
            exch_hf(j,i)=exch_hf(i,j)
            exchen_hf=exchen_hf+exch_hf(i,j)
          end do
        end do

!! PRINTING !!
        write(*,*) " ----------------------------------------- "
        write(*,*) "  HARTREE-FOCK-TYPE EXCHANGE ENERGY TERMS  "
        write(*,*) " ----------------------------------------- "
        write(*,*) " "
        CALL MPRINT2(exch_hf,nat,maxat)
        write(*,'(2x,a29,x,f14.7)') "Hartree-Fock exchange energy:",exchen_hf
        write(*,*) " "
        call cpu_time(xtime2)
        write(*,'(a39,f10.1,a2)')'(Elapsed time :: Exact-exchange energy ',xtime2-xtime,'s)' 
        write(*,*) " "
        xtime=xtime2
      end if
      end if !MMO- non-analytical skip ends here

!! CHECKING ACCURACY OF THE TWO-ELECTRON PART !!
!! HF ONLY !!
      if(ihf.eq.1) then
        evee=coulen+exchen_hf
        write(*,'(2x,a39,x,f14.7)') "Total two-electron part (coul+exch_hf):",evee
        if(evee0.ne.ZERO) then
          twoelerr=(evee-evee0)*tokcal
          write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",twoelerr
          write(*,*) " "
        else
          evee0=evee    
        end if

!! DFT CASE !!
      else
        exchen=ZERO
        do i=1,nat
          exch(i,i)=exch(i,i)+xmix*exch_hf(i,i)
          exchen=exchen+exch(i,i)
          do j=i+1,nat
            exch(i,j)=exch(i,j)+xmix*exch_hf(i,j)
            exch(j,i)=exch(i,j)
            exchen=exchen+exch(i,j)
          end do
        end do
        if(idoex.eq.1) then
          write(*,*) " ------------------------ "
          write(*,*) "  HYBRID KS-DFT XC TERMS  "
          write(*,*) " ------------------------ "
          write(*,*) " "
          call MPRINT2(exch,nat,maxat)
          write(*,'(2x,a34,x,f14.7)') "Total exchange-correlation energy:",exchen
          write(*,*) " "
        end if
        evee=coulen+exchen
        write(*,'(2x,a37,x,f14.7)') "KS-DFT electron-electron energy (au):",evee            
        if(evee0.ne.ZERO) then
          twoelerr=(evee-evee0)*tokcal
          write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",twoelerr
        else
          evee0=evee
        end if
        write(*,*) " "
      end if

!! ZERO ERROR STRATEGY FOR VEE PART !!
      if(ianalytical.eq.1) twoeltoler=232000
      write(*,'(2x,a55,x,f8.2)') "Max error accepted on the two-electron part (kcal/mol):",twoeltoler
      write(*,*) " "

      if(abs(twoelerr).gt.twoeltoler) then

!! NOW ROTATED GRID, ANGLE CONTROLLED BY # GRID SECTION !!
!! CHECK diatXC PAPER FOR OPTIMIZED VALUES: phb=0.162d0 and later 0.182d0 for 40 146 !!
      phb=phb22
      pha=ZERO 
      write(*,'(2x,a20,x,f10.6,x,f10.6)') "Rotating for angles:",pha,phb
      nat0=nat
      call prenumint(ndim,itotps,nat0,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)

!! TO MOs !!
      do k=1,itotps
        do j=1,nocc
          xx=ZERO
          do i=1,igr 
            xx=xx+c(i,j)*chppha(k,i) 
          end do
          chp2phas(j,k)=xx
        end do
      end do

!! DOING INTEGRATIONS FOR ONLY THE ATOMIC TERMS !!
      do icenter=1,nat
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          x0=wp(ifut)*omp2(ifut,icenter)
          dx0=pcoord(ifut,1)
          dy0=pcoord(ifut,2)
          dz0=pcoord(ifut,3)

!! COULOMB PART !!
          f3=ZERO
          do jfut=iatps*(icenter-1)+1,iatps*icenter
            x1=wppha(jfut)*omp2pha(jfut,icenter)
            dx1=pcoordpha(jfut,1)
            dy1=pcoordpha(jfut,2)
            dz1=pcoordpha(jfut,3)
            dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
            if(dist.gt.1.0d-12) f3=f3+rho(ifut)*rhopha(jfut)*x1*x0/dist !! MG: Could be controlled using the thr2 variable !!
          end do
          coul0(icenter,3)=coul0(icenter,3)+f3/TWO
        end do
      end do

!! HF-TYPE EXCHANGE !!
      if(idoex.eq.1) then

!! AS BEFORE, MORE CONVOLUTED LOOP STRUCTURE !!
        exch_hfk=ZERO
        !DIR$ NOPARALLEL
        do icenter=1,nat
          ioffset=(icenter-1)*iatps
          istart=0
          iend=0
          !DIR$ NOPARALLEL
          do ik=1,(ithreads-1)
            istart(ik)=ioffset+((ik-1)*ispace)+1
            iend(ik)=ioffset+(ik*ispace)
          end do
          istart(ithreads)=iend(ithreads-1)+1
          iend(ithreads)=icenter*iatps

!! AGAIN THE TWO CALLS ARE CRUCIAL !!
          call omp_set_dynamic(.false.)
          call omp_set_num_threads(ithreads)
          !DIR$ PARALLEL
          do ik=1,ithreads
            do ifut=istart(ik),iend(ik)
              x0=wp(ifut)*omp2(ifut,icenter)
              dx0=pcoord(ifut,1)
              dy0=pcoord(ifut,2)
              dz0=pcoord(ifut,3)
              f3k(ik)=ZERO
              do jfut=iatps*(icenter-1)+1,iatps*icenter
                x1=wppha(jfut)*omp2pha(jfut,icenter)
                dx1=pcoordpha(jfut,1)
                dy1=pcoordpha(jfut,2)
                dz1=pcoordpha(jfut,3)
                dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                if(dist.gt.1.0d-8) then !! MG: Could be controlled using the thr2 variable. In fact I don't remember why was 10^-8 !!
                  do i=1,nocc-1
                    f2=x0*chp2s(i,ifut)*chp2s(i,ifut)
                    f3k(ik)=f3k(ik)+f2*x1*chp2phas(i,jfut)*chp2phas(i,jfut)/dist
                    do j=i+1,nocc
                      f2=TWO*x0*chp2s(i,ifut)*chp2s(j,ifut)
                      f3k(ik)=f3k(ik)+f2*x1*chp2phas(i,jfut)*chp2phas(j,jfut)/dist !! MG: LC-wPBE programmed in version 3.1 !!
                    end do
                  end do
                  f2=x0*chp2s(nocc,ifut)*chp2s(nocc,ifut)
                  f3k(ik)=f3k(ik)+f2*x1*chp2phas(nocc,jfut)*chp2phas(nocc,jfut)/dist
                end if
              end do
              exch_hfk(icenter,ik)=exch_hfk(icenter,ik)+f3k(ik)
            end do
          end do
        end do

!! ADDING THE TERMS INTO THE ORIGINAL coul0 MATRIX !!
        do icenter=1,nat
          do isum=1,ithreads
            coul0(icenter,4)=coul0(icenter,4)-exch_hfk(icenter,isum)
          end do
        end do
      end if

!! CHECKING NEW VALUES !!
      do i=1,nat
        coul0(i,1)=coul(i,i)
        if(idoex.eq.1) coul0(i,2)=exch_hf(i,i)
      end do
      deltaee=ZERO
      do i=1,nat
        deltaee=deltaee+coul0(i,1)-coul0(i,3)
        if(idoex.eq.1) deltaee=deltaee+xmix*(coul0(i,2)-coul0(i,4))
      end do
      deltaee=deltaee*tokcal
      phabest=ONE-(twoelerr/deltaee)
      write(*,'(2x,a25,x,f8.2)') "New error after rotation:",twoelerr-deltaee
      if(twoelerr-deltaee*twoelerr.gt.ZERO) write(*,*) " WARNING: New error with same sign"
      write(*,'(2x,a18,x,f14.7)') "Damping parameter:",phabest !MMO- not dumping !! ...

!! INTERPOLATING ENERGIES, REPLACING OLD COULOMB AND EXCHANGE TERMS !!
      do i=1,nat
        coul(i,i)=coul0(i,1)*phabest+(ONE-phabest)*coul0(i,3)
        if(idoex.eq.1) exch_hf(i,i)=coul0(i,2)*phabest+(ONE-phabest)*coul0(i,4)
      end do
      write(*,*) " "
      call cpu_time(xtime2)
      write(*,'(a41,f10.1,a2)')'(Elapsed time :: Interpolated zero-error ',xtime2-xtime,'s)' 
      xtime=xtime2

!! PRINTING !!
      write(*,*) " "
      write(*,*) " ------------------------------------------------------- "
      write(*,*) "  INTERPOLATED COULOMB (ELECTRON-ELECTRON) ENERGY TERMS  "
      write(*,*) " ------------------------------------------------------- "
      write(*,*) " "
      call MPRINT2(coul,nat,maxat)
      coulen=ZERO
      do i=1,nat
        do j=i,nat
          coulen=coulen+coul(i,j)
        end do
      end do
      write(*,'(2x,a15,x,f14.7)') "Coulomb energy:",coulen
      write(*,*) " "
      if (idofr.eq.1) then
        line=' FRAGMENT ANALYSIS: Coulomb energy ' 
        call group_by_frag_mat(1,line,coul)
      end if

      if(idoex.eq.1) then
        write(*,*) " ------------------------------------------------------ "
        write(*,*) "  INTERPOLATED HARTREE-FOCK-TYPE EXCHANGE ENERGY TERMS  "
        write(*,*) " ------------------------------------------------------ "
        write(*,*) " "
        call MPRINT2(exch_hf,NAT,maxat)
        exchen_hf=ZERO
        do i=1,nat
          do j=i,nat
            exchen_hf=exchen_hf+exch_hf(i,j)
          end do
        end do
        write(*,'(2x,a29,x,f14.7)') "Hartree-Fock exchange energy:",exchen_hf
        write(*,*) " "

        if(ihf.ne.1) then
          exchen=ZERO
          do i=1,nat
            exch(i,i)=exch(i,i)+xmix*(exch_hf(i,i)-coul0(i,2))
            exchen=exchen+exch(i,i)
            do j=i+1,nat
              exchen=exchen+exch(i,j)
            end do
          end do

          write(*,*) " ------------------------------------- "
          write(*,*) "  INTERPOLATED HYBRID KS-DFT XC TERMS  "
          write(*,*) " ------------------------------------- "
          write(*,*) " "
          call MPRINT2(exch,nat,maxat)
          write(*,'(2x,a34,x,f14.7)') "Total exchange-correlation energy:",exchen
          write(*,*) " "
          if (idofr.eq.1) then
            line=' FRAGMENT ANALYSIS: Final Exc Decomposition ' 
            call group_by_frag_mat(1,line,exch)
          end if
        end if
      else    
        exchen=ZERO
        do i=1,nat
          exchen=exchen+exch(i,i)
          do j=i+1,nat
            exchen=exchen+exch(i,j)
          end do
        end do
      end if
 
!! TOTAL TWO-ELECTRON PART !!
      if(ihf.eq.0) then
        evee=coulen+exchen
        write(*,'(2x,a32,x,f14.7)') "KS-DFT electron-electron energy:",evee             
        twoelerr=(evee-evee0)*tokcal
        write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",twoelerr
      else
        evee=coulen+exchen_hf
        write(*,'(2x,a39,x,f14.7)') "Total two-electron part (coul+exch_HF):",evee            
        twoelerr=(evee-evee0)*tokcal
        write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",twoelerr
      end if
      write(*,*) " "

      end if

!! END ZERO ERROR STRATEGY !!

!! FINAL PRINTING !!
      if(ihf.eq.1) then
        xtot=ZERO
        do i=1,nat
          eto(i,i)=eto(i,i)+coul(i,i)+exch_hf(i,i)
          xtot=xtot+eto(i,i)
          do j=i+1,nat
            eto(i,j)=eto(i,j)+coul(i,j)+exch_hf(i,j)
            eto(j,i)=eto(i,j)
            xtot=xtot+eto(i,j)
          end do
        end do
        etot=xtot

        write(*,*) " -------------------------------------------- "
        write(*,*) "  FUZZY ATOMS Hartree-Fock ENERGY COMPONENTS  "
        write(*,*) " -------------------------------------------- "
        write(*,*) " "
        call MPRINT2(eto,nat,maxat)
        write(*,'(2x,a26,x,f14.7)') "Total integrated energy  :",etot  
        write(*,'(2x,a26,x,f14.7)') "Total energy in Fchk file:",escf  
        if(ifield.eq.1) then
          write(*,'(2x,a20,x,f14.7)') "Total dipole energy:",edipole
          err=(etot-escf+edipole)
        else
          err=(etot-escf)
        end if
        write(*,'(2x,a23,x,f14.7)') "Integration error (au):",err
        write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",err*tokcal
      else
        xtot=ZERO
        do i=1,nat
          eto(i,i)=eto(i,i)+coul(i,i)+xmix*exch_hf(i,i)
          xtot=xtot+eto(i,i)
          do j=i+1,nat
            eto(i,j)=eto(i,j)+coul(i,j)+xmix*exch_hf(i,j)
            eto(j,i)=eto(i,j)
            xtot=xtot+eto(i,j)
          end do
        end do
        etot=xtot

        write(*,*) " -------------------------------------- "
        write(*,*) "  FUZZY ATOMS KS-DFT ENERGY COMPONENTS  "
        write(*,*) " -------------------------------------- "
        write(*,*) " "
        call MPRINT2(eto,nat,maxat)
        write(*,'(2x,a26,x,f14.7)') "Total KS-DFT energy      :",etot  
        write(*,'(2x,a26,x,f14.7)') "Total energy in Fchk file:",escf
        if(ifield.eq.1) then
          write(*,'(2x,a20,x,f14.7)') "Total dipole energy:",edipole
          err=(etot-escf+edipole)
        else
          err=(etot-escf)
        end if
        write(*,'(2x,a23,x,f14.7)') "Integration error (au):",err
        write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",err*tokcal
      end if
      if (idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Energy Decomposition' 
        call group_by_frag_mat(1,line,eto)
      end if

      if(ianalytical.eq.0) then
        DEALLOCATE(chp2,chp2pha,rhopha,chppha,wppha,pcoordpha,omppha,omp2pha)
        DEALLOCATE(chp2s,chp2phas)
!! MG: DEALLOCATING THE PARALLEL PART MISSING (I'VE BEEN LAZY)
      else
!        DEALLOCATE(chp2,rho) !MMO- rho is no longer allocated in the current version?
        DEALLOCATE(chp2)
      end if

      end

!! ***** !!

      subroutine numint_one_uhf(ndim,itotps,wp,rho,omp,omp2,pcoord,chp,eto)
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      character*80 line
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common/actual/iact,jat,icenter
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/energ0/ekin0,eelnuc0,evee0,etot0
      common/exchg/exch(maxat,maxat),xmix
      dimension eto(maxat,maxat)
      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3)
      dimension omp(itotps),omp2(itotps,nat),rho(itotps)
! (TO DO) xkdens = kinetic energy density (MG)
! (TO DO) Laplacian!
      dimension xkdens(itotps)

      dimension Epa(maxat,maxat),ekin(maxat,maxat),enuc(maxat,maxat)
      dimension eecp(maxat)
      allocatable :: chp2(:,:),scr(:),chpaux(:,:), chp3(:,:)

      idofr   =  Iopt(40)
      ipoints =  itotps
      iatps   =  nang*nrad

      write(*,*) " "
      write(*,*) " --------------------------- "
      write(*,*) "  GENERAL ONE-ELECTRON PART  "
      write(*,*) " --------------------------- "
      write(*,*) " "

      ALLOCATE(chp2(itotps,nalf),chp3(itotps,nb))

!! TO MOs !!
      do k=1,ipoints
        do j=1,nalf
          xx=ZERO
          xxb=ZERO
          do i=1,igr 
            xx=xx+c(i,j)*chp(k,i)
            if(j.le.nb) xxb=xxb+cb(i,j)*chp(k,i)
          end do
          chp2(k,j)=xx
          if(j.le.nb) chp3(k,j)=xxb
        end do
      end do

!! ZEROING MATRICES !!
      epa=ZERO
      ekin=ZERO
      eto=ZERO
      enuc=ZERO
      
!! DOING ENERGY PARTITION !!
!! ONE ELECTRON PART !!

!! IN CASE OF HAVING PSEUDOPOTENTIAL !!
      call mga_misc(iecp,eecp) !! MAYBE WE SHOULD DO A BETTER VERSION OF THIS... !!
      if(iecp.eq.1) then
        write(*,*) " ADDING ECP ATOMIC ENERGIES TO E-N TERMS "
        do i=1,nat       
          epa(i,i)=eecp(i)
        end do 
      end if 
      xtot=ZERO
      do icenter=1,nat
        do jcenter=1,nat     
          x=ZERO
          zz=zn(jcenter)
          dx0=coord(1,jcenter)
          dy0=coord(2,jcenter)
          dz0=coord(3,jcenter)
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            distx=pcoord(ifut,1)-dx0
            disty=pcoord(ifut,2)-dy0
            distz=pcoord(ifut,3)-dz0
            rr=dsqrt(distx*distx+disty*disty+distz*distz)
            if(rr.gt.1.0d-12) then
              x=x+zz*wp(ifut)*omp2(ifut,icenter)*rho(ifut)/rr
            end if
          end do
          xtot=xtot+x
          epa(icenter,jcenter)=epa(icenter,jcenter)-x
        end do
      end do
      xtot=ZERO
      do i=1,nat
        xtot=xtot+epa(i,i)
        do j=1,i-1  
          epa(i,j)=(epa(i,j)+epa(j,i))
          epa(j,i)=epa(i,j)
          xtot=xtot+epa(i,j)
        end do
      end do
      eelnuc=xtot

!! PRINTING MATRIX !!
      write(*,*) " ----------------------------- "
      write(*,*) "  ELECTRON-NUCLEAR ATTRACTION  "
      write(*,*) " ----------------------------- "
      write(*,*) " "
      call MPRINT2(epa,nat,maxat)
      write(*,'(2x,a23,x,f14.7)') "Electron-nuclei energy:",eelnuc

!! CHECKING ACCURACY (IF ENERGIES ADDED TO THE fchk FILE)
      if(eelnuc0.ne.ZERO) then
!        write(*,'(a31,f8.2)') 'Integration error (kcal/mol): ',(eelnuc-eelnuc0)*tokcal !MMO- same as restricted case
        write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",(eelnuc+xxdip-eelnuc0)*tokcal
        write(*,*) ' '
      else
        eelnuc0=eelnuc
      end if
      if (idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Electron-nuclei attraction '
        call group_by_frag_mat(1,line ,epa)
      end if

!! GENERATING GRID FOR 2nd DERIVATIVE OVER AOs !!
      ALLOCATE(chpaux(itotps,igr))
      call dpoints(chpaux,pcoord)

!! NOW TO MO AND MULTPLY BY MOs AND SUM OVER MOs USING SCRATCH ARRAY !!
!! HENCE WE HAVE DENSITY IN scr EXCEPT FOR THE FACTOR OF TWO !!
      ALLOCATE(scr(itotps))
      do k=1,ipoints
        x=ZERO
        do j=1,nalf
          xx=ZERO
          xxb=ZERO
          do i=1,igr 
            xx=xx+c(i,j)*chpaux(k,i) 
            if(j.le.nb) xxb=xxb+cb(i,j)*chpaux(k,i) 
          end do
          x=x+chp2(k,j)*xx
          if(j.le.nb) x=x+chp3(k,j)*xxb
        end do
        scr(k)=x/TWO
! (TO DO) Kinetic energy density (MG)
!       xkdens(k)=scr(k)
! (TO DO) Laplacian!
      end do

      xtot=ZERO
      do icenter=1,nat
        x=ZERO
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          x=x+wp(ifut)*scr(ifut)*omp2(ifut,icenter)
        end do
        ekin(icenter,icenter)=ekin(icenter,icenter)+x
        xtot=xtot+x
      end do
      ekinen=xtot

!! PRINTING MATRIX !!
      write(*,*) " ---------------- "
      write(*,*) "  KINETIC ENERGY  "
      write(*,*) " ---------------- "
      write(*,*) " "
      call MPRINT2(ekin,nat,maxat)
      write(*,'(2x,a15,x,f14.7)') "Kinetic energy:",ekinen

!! CHECKING ACCURACY (IF ENERGIES ADDED TO THE fchk FILE ) !!
      if(ekin0.ne.ZERO) then
        write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",(ekinen-ekin0)*tokcal                      
      else
        ekin0=ekinen  
      end if
      write(*,*) " "
      if (idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Kinetic energy ' 
        call group_by_frag_mat(1,line ,ekin)
      end if

!! ADDING NUCLEAR REPULSION !!
      erep=ZERO
      do i=1,nat
        eto(i,i)=ekin(i,i)+epa(i,i)
        do j=i+1,nat
          dist=dsqrt((coord(1,i)-coord(1,j))**2 +(coord(2,i)-coord(2,j))**2+(coord(3,i)-coord(3,j))**2)
          eto(i,j)=ekin(i,j)+epa(i,j)+zn(i)*zn(j)/dist
          eto(j,i)=epa(i,j)

!! SAVING ALSO FOR PRINTING PURPOSES, I THINK WE SHOULD GIVE IT... COMES FOR FREE !!
          enuc(i,j)=zn(i)*zn(j)/dist
          enuc(j,i)=enuc(i,j)
          erep=erep+zn(i)*zn(j)/dist
        end do
      end do

!! PRINTING MATRIX !!
      write(*,*) " --------------------------- "
      write(*,*) "  NUCLEAR-NUCLEAR REPULSION  "
      write(*,*) " --------------------------- "
      write(*,*) " "
      call MPRINT2(enuc,nat,maxat)
      write(*,'(2x,a25,x,f14.7)') "Nuclear repulsion energy:",erep
      write(*,*) " "

!! TOTAL ENERGY UP TO HERE !!
      etot=ekinen+eelnuc+erep
      etot0=ekin0+eelnuc0+erep

      DEALLOCATE(chp2,scr,chpaux,chp3)
      end

!! ***** !!

      subroutine numint_two_uhf(ndim,itotps,wp,omp,omp2,pcoord,chp,rho,eto)

      use ao_matrices
      use integration_grid
      use OMP_LIB
      use IFPORT

      implicit real*8(a-h,o-z)

      integer external OMP_GET_NUM_PROCS

      include 'parameter.h'
      parameter(maxpass=5)

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/twoel/twoeltoler !MMO- adding common
      common /iops/iopt(100)
      common/actual/iact,jat,icenter
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/energ0/ekin0,eelnuc0,evee0,etot0
      common/exchg/exch(maxat,maxat),xmix
      common/efield/field(4),edipole !MMO- adding common
!! FROM # GRID SECTION !!
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22
      common /modgrid2/thr3

      dimension eto(maxat,maxat),Excmp(maxat,maxat)
      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3)
      dimension omp(itotps),omp2(itotps,nat),rho(itotps)

      dimension coul0(maxat,2*maxpass)
      dimension coul(maxat,maxat),exch_hf(maxat,maxat)

      character*80 line
      character*100 threadenv

      allocatable :: chppha(:,:), wppha(:),pcoordpha(:,:),omppha(:)
      allocatable :: chp2pha(:,:),omp2pha(:,:),ibaspointpha(:),rhopha(:)
      allocatable :: chp2(:,:), chp2b(:,:), chp2phab(:,:)
      allocatable :: fij(:,:),fijb(:,:),xocc(:,:),xoccb(:,:),Excmpb(:,:)

!! FOR ENPART PARALLEL !!
      allocatable :: chp2s(:,:),chp2phas(:,:),chp2bs(:,:),chp2phabs(:,:) !! SWITCHED ORDER OF COLUMNS AND ROWS !!
      allocatable :: ijpaircount(:,:),istart(:),iend(:) !! ATOM PAIRS INCLUDED, AND FOR SLICING chp MATRICES !!
      allocatable :: exch_hfk(:,:)
      allocatable :: exch_hfij(:,:)
      allocatable :: f3k(:)

      ihf         = Iopt(18) 
      idofr       = Iopt(40)
      ithrebod    = Iopt(44)
      ipolar      = Iopt(46)
      ifield      = Iopt(47) !MMO- adding iopt
      ianalytical = Iopt(91) !MMO- analytical

      !thr2=thr3 !! THRESH FOR NUMERICAL INTEGRATION (DISTANCE) !!
      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if
      idoex=0
      if(xmix.gt.ZERO) idoex=1
      iatps=nang*nrad

!! DOING ENERGY PARTITION !!
      npass=2
      write(*,*) " "
      write(*,*) " --------------------------- "
      write(*,*) "  GENERAL TWO-ELECTRON PART  "
      write(*,*) " --------------------------- "
      write(*,*) " "
      do i=1,nat
        do j=1,nat
          coul(i,j)=ZERO
          exch_hf(i,j)=ZERO
        end do
        do j=1,2*npass
          coul0(i,j)=ZERO
        end do
      end do

      ALLOCATE(chp2(itotps,nalf),chp2b(itotps,nb))
      ALLOCATE(chp2s(nalf,itotps),chp2bs(nb,itotps)) !! TODO: CHANGING THAT NEVER USE CHP2 !!

!! TO MOs !!
      do k=1,itotps
        do j=1,nalf
          xx=ZERO
          xxb=ZERO
          do i=1,igr 
            xx=xx+c(i,j)*chp(k,i) 
            if(j.le.nb) xxb=xxb+cb(i,j)*chp(k,i) 
          end do
          if(j.le.nb) then
            chp2b(k,j)=xxb
            chp2bs(j,k)=xxb
          end if
          chp2(k,j)=xx
          chp2s(j,k)=xx
        end do
      end do

      if(ianalytical.eq.0) then !MMO- skip if analytical

!! ROTATED GRID, ANGLE CONTROLLED BY # GRID SECTION !!
!! CHECK diatXC PAPER FOR OPTIMIZED VALUES: phb=0.162d0 and later 0.182d0 for 40 146 !!
      phb=phb12
      pha=ZERO 
      write(*,'(2x,a20,x,f10.6,x,f10.6)') "Rotating for angles:",pha,phb
      ALLOCATE(wppha(itotps),omppha(itotps),omp2pha(itotps,nat))
      ALLOCATE(chppha(itotps,ndim),pcoordpha(itotps,3),ibaspointpha(itotps))
      ALLOCATE(chp2pha(itotps,nalf),rhopha(itotps),chp2phab(itotps,nb))
      ALLOCATE(chp2phas(nalf,itotps),chp2phabs(nb,itotps))

      call prenumint(ndim,itotps,nat,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)

!! TO MOs !!
      do k=1,itotps
        xtot=ZERO 
        xtotb=ZERO 
        do j=1,nalf
          xx=ZERO
          xxb=ZERO
          do i=1,igr 
            xx=xx+c(i,j)*chppha(k,i) 
            if(j.le.nb) xxb=xxb+cb(i,j)*chppha(k,i) 
          end do
          chp2pha(k,j)=xx
          chp2phas(j,k)=xx
          if(j.le.nb) then
            chp2phab(k,j)=xxb
            chp2phabs(j,k)=xxb
            xtotb=xtotb+xxb*xxb
          end if 
          xtot=xtot+xx*xx
        end do
        rhopha(k)=xtot+xtotb 
      end do
      end if !MMO- non-analytical skip ends here

!! COULOMB !!
      call cpu_time(xtime)
      if(ianalytical.eq.0) then
        call calc_coul(itotps,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,rho,rhopha,coul)
      else
        CALL calc_twoel_analytical(itotps,wp,omp2,pcoord,rho,chp2,chp2b,coul,exch_hf)
      end if
      call cpu_time(xtime2)
      write(*,'(a32,f10.1,a2)')'(Elapsed time :: Coulomb energy ',xtime2-xtime,'s)'
      write(*,*) " "
      xtime=xtime2

!! EXCHANGE PART !!
      if(ianalytical.eq.0) then !MMO- skip if analytical
      if(idoex.eq.1) then

!! COMPUTING MULTIPOLAR APPROACH HERE !!
        norb2=nalf*(nalf+1)/2
        norb2b=nb*(nb+1)/2
        ALLOCATE(fij(itotps,norb2),xocc(norb2,norb2))
        ALLOCATE(fijb(itotps,norb2b),xoccb(norb2b,norb2b),Excmpb(maxat,maxat))
        do ii=1,itotps
          irun=0
          irunb=0
          do jj=1,nalf
            do kk=jj,nalf
              irun=irun+1
              fij(ii,irun)=chp2(ii,jj)*chp2(ii,kk)
              if(ii.eq.1) then
                xocc(irun,irun)=-ONE
                if(jj.ne.kk) xocc(irun,irun)=-TWO
              end if
              if(jj.le.nb.and.kk.le.nb) then
                irunb=irunb+1
                fijb(ii,irunb)=chp2b(ii,jj)*chp2b(ii,kk)
                if(ii.eq.1) then
                  xoccb(irunb,irunb)=-ONE
                  if(jj.ne.kk) xoccb(irunb,irunb)=-TWO
                end if
              end if 
            end do
          end do
        end do
        call multipolar(nalf,itotps,wp,omp2,pcoord,fij,xocc,Excmp)
        call multipolar(nb,itotps,wp,omp2,pcoord,fijb,xoccb,Excmpb)
        do ii=1,nat
          do jj=ii+1,nat
            Excmp(ii,jj)=Excmp(ii,jj)+Excmpb(ii,jj)
            Excmp(jj,ii)=Excmp(ii,jj)
          end do 
        end do 
        DEALLOCATE(fij,fijb,xocc,xoccb,Excmpb)

!! MG: REORDERING THE LOOPS FOR PARALLELIZATION PURPOSES !!
!! MIMICKING STRATEGY FROM RHF PART !!
        write(*,*) " ------------------------------------------------- "
        write(*,*) "  EVALUATING HARTREE-FOCK-TYPE EXCHANGE INTEGRALS  "
        write(*,*) " ------------------------------------------------- "
        write(*,*) " "
        write(*,'(2x,a22,x,i6,x,a9)') "Two-el integrations for",nalf*(nalf+1),"functions"
        write(*,'(2x,a36,x,f10.6)') "Threshold for atom pair calculation :",threbod

!! FIRST ONLY SAME CENTER TERMS !!
!! EVALUATING NUMBER OF CORES FOR SPLITTING THE CALCULATION BY THREADS !!
        call getenv('OMP_NUM_THREADS',threadenv)
        if(trim(threadenv)=='') then
          write(*,*) " OMP_NUM_THREADS not set"
          ithreadenv=ZERO
        else
          read(unit=threadenv,FMT='(I4)') ithreadenv
        end if
        !$OMP PARALLEL
        iprocs=OMP_GET_MAX_THREADS()
        ithreads=INT(OMP_GET_NUM_PROCS())
        !$OMP END PARALLEL
        if(ithreadenv.ne.ZERO) then
          write(*,'(2x,a43,x,i3,x,a14,x,i3,x,a24)') "Two-el integration will be distributed over",ithreadenv,"threads out of",
     &    ithreads,"available hardware cores"
          ithreads=ithreadenv
        else
          write(*,'(2x,a43,x,i3,x,a14,x,i3,x,a24)') "Two-el integration will be distributed over",ithreads,"threads out of",
     &    ithreads,"available hardware cores"
        end if

!! TO ENSURE PROPER SLICING BY THREADS !!
        ALLOCATE(f3k(ithreads))
        ALLOCATE(exch_hfk(nat,ithreads))
        ALLOCATE(istart(ithreads),iend(ithreads))
        itilerest=mod(iatps,ithreads)
        ispace=(iatps-itilerest)/ithreads

!! MORE CONVOLUTED LOOP STRUCTURE, AVOIDED PROBLEM OF MAX PARALLEL 8 CORES !!
        exch_hfk=ZERO
        !DIR$ NOPARALLEL
        do icenter=1,nat
          ioffset=(icenter-1)*iatps
          istart=0
          iend=0
          !DIR$ NOPARALLEL
          do ik=1,(ithreads-1)
            istart(ik)=ioffset+((ik-1)*ispace)+1
            iend(ik)=ioffset+(ik*ispace)
          end do
          istart(ithreads)=iend(ithreads-1)+1
          iend(ithreads)=icenter*iatps

!! THE TWO CALLS ARE CRUCIAL !!
          call omp_set_dynamic(.false.)
          call omp_set_num_threads(ithreads)
          !DIR$ PARALLEL
          do ik=1,ithreads
            do ifut=istart(ik),iend(ik)
              x0=wp(ifut)*omp2(ifut,icenter)
              dx0=pcoord(ifut,1)
              dy0=pcoord(ifut,2)
              dz0=pcoord(ifut,3)
              f3k(ik)=ZERO
              do jfut=iatps*(icenter-1)+1,iatps*icenter
                x1=wppha(jfut)*omp2pha(jfut,icenter)
                dx1=pcoordpha(jfut,1)
                dy1=pcoordpha(jfut,2)
                dz1=pcoordpha(jfut,3)
                dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                if(dist.gt.1.0d-12) then !! MG: Could be controlled using the thr2 variable !!
                  do i=1,nalf-1
                    f2=x0*chp2s(i,ifut)*chp2s(i,ifut)
                    f3k(ik)=f3k(ik)+f2*x1*chp2phas(i,jfut)*chp2phas(i,jfut)/dist
                    if(i.le.nb) then
                      f2b=x0*chp2bs(i,ifut)*chp2bs(i,ifut)
                      f3k(ik)=f3k(ik)+f2b*x1*chp2phabs(i,jfut)*chp2phabs(i,jfut)/dist
                    end if
                    do j=i+1,nalf
                      f2=TWO*x0*chp2s(i,ifut)*chp2s(j,ifut)
                      f3k(ik)=f3k(ik)+f2*x1*chp2phas(i,jfut)*chp2phas(j,jfut)/dist !! MG: LC-wPBE programmed in version 3.1 !!
                      if(j.le.nb) then
                        f2b=TWO*x0*chp2bs(i,ifut)*chp2bs(j,ifut)
                        f3k(ik)=f3k(ik)+f2b*x1*chp2phabs(i,jfut)*chp2phabs(j,jfut)/dist
                      end if
                    end do
                  end do
                  f2=x0*chp2s(nalf,ifut)*chp2s(nalf,ifut)
                  f3k(ik)=f3k(ik)+f2*x1*chp2phas(nalf,jfut)*chp2phas(nalf,jfut)/dist

!! FOR THE nalf = nb CASE IT IS REQUIRED TO ADD THIS !!
                  if(nalf.eq.nb) then
                    f2b=x0*chp2bs(nb,ifut)*chp2bs(nb,ifut)
                    f3k(ik)=f3k(ik)+f2b*x1*chp2phabs(nb,jfut)*chp2phabs(nb,jfut)/dist
                  end if
                end if
              end do
              exch_hfk(icenter,ik)=exch_hfk(icenter,ik)+f3k(ik)
            end do
          end do
        end do

!! ADDING THE TERMS INTO THE ORIGINAL exch_hf MATRIX !!
        do icenter=1,nat
          do isum=1,ithreads
             exch_hf(icenter,icenter)=exch_hf(icenter,icenter)-exch_hfk(icenter,isum)
          end do
        end do

!! NOW PAIRS OF CENTERS !!
        iterms=0
        ipaircounter=0
        ALLOCATE(ijpaircount(nat*nat,2))
        do icenter=1,nat
          do jcenter=icenter+1,nat
            bx0=bo(icenter,jcenter)
            if(bx0.ge.threbod) then
              ipaircounter=ipaircounter+1
              ijpaircount(ipaircounter,1)=icenter
              ijpaircount(ipaircounter,2)=jcenter
            else
              iterms=iterms+1
            end if
          end do
        end do
        write(*,'(2x,a34,x,i5,x,a10)') "Skipping numerical integration for",iterms,"atom pairs"
        write(*,*) " "

!! AGAIN, PREPARING FOR SLICING AND BLOCKING PARALLELIZATION OF SOME LOOPS !!
        ALLOCATE(exch_hfij(ipaircounter,ithreads))
        itilerest=mod(iatps,ithreads)
        ispace=(iatps-itilerest)/ithreads
        exch_hfij=ZERO
        !DIR$ NOPARALLEL
        do numpairnat=1,ipaircounter
          icenter=ijpaircount(numpairnat,1)
          jcenter=ijpaircount(numpairnat,2)
          ioffset=(icenter-1)*iatps
          istart=0
          iend=0
          !DIR$ NOPARALLEL
          do ik=1,(ithreads-1)
            istart(ik)=ioffset+((ik-1)*ispace)+1
            iend(ik)=ioffset+(ik*ispace)
          end do
          istart(ithreads)=iend(ithreads-1)+1
          iend(ithreads)=icenter*iatps
          !DIR$ PARALLEL
          do ik=1,ithreads
            do ifut=istart(ik),iend(ik)
              x0=wp(ifut)*omp2(ifut,icenter)
              dx0=pcoord(ifut,1)
              dy0=pcoord(ifut,2)
              dz0=pcoord(ifut,3)
              f3k(ik)=ZERO
              do jfut=iatps*(jcenter-1)+1,iatps*jcenter
                x1=wp(jfut)*omp2(jfut,jcenter)
                dx1=pcoord(jfut,1)
                dy1=pcoord(jfut,2)
                dz1=pcoord(jfut,3)
                dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                if(dist.gt.1.0d-12) then !! Could be controlled using the thr2 variable !!
                  do i=1,nalf-1
                    f2=x0*chp2s(i,ifut)*chp2s(i,ifut)
                    f3k(ik)=f3k(ik)+f2*x1*chp2s(i,jfut)*chp2s(i,jfut)/dist
                    if(i.le.nb) then
                      f2b=x0*chp2bs(i,ifut)*chp2bs(i,ifut)
                      f3k(ik)=f3k(ik)+f2b*x1*chp2bs(i,jfut)*chp2bs(i,jfut)/dist
                    end if
                    do j=i+1,nalf
                      f2=TWO*x0*chp2s(i,ifut)*chp2s(j,ifut)
                      f3k(ik)=f3k(ik)+f2*x1*chp2s(i,jfut)*chp2s(j,jfut)/dist !! LC-wPBE programmed in version 3.1 !!
                      if(j.le.nb) then
                        f2b=TWO*x0*chp2bs(i,ifut)*chp2bs(j,ifut)
                        f3k(ik)=f3k(ik)+f2b*x1*chp2bs(i,jfut)*chp2bs(j,jfut)/dist
                      end if
                    end do
                  end do
                  f2=x0*chp2s(nalf,ifut)*chp2s(nalf,ifut)
                  f3k(ik)=f3k(ik)+f2*x1*chp2s(nalf,jfut)*chp2s(nalf,jfut)/dist

!! AGAIN, ONLY WHEN nalf = nb !!
                  if(nalf.eq.nb) then
                    f2b=x0*chp2bs(nb,ifut)*chp2bs(nb,ifut)
                    f3k(ik)=f3k(ik)+f2b*x1*chp2bs(nb,jfut)*chp2bs(nb,jfut)/dist
                  end if
                end if
              end do
              exch_hfij(numpairnat,ik)=exch_hfij(numpairnat,ik)+f3k(ik)
            end do
          end do
        end do

!! AGAIN, ADDING THE TERMS INTO THE ORIGINAL exch_hf MATRIX !!
        do numpairnat=1,ipaircounter
          icenter=ijpaircount(numpairnat,1)
          jcenter=ijpaircount(numpairnat,2)
          do isum=1,ithreads
            exch_hf(icenter,jcenter)=exch_hf(icenter,jcenter)-exch_hfij(numpairnat,isum)
          end do
        end do   

!! ACCOUNTING IF SOME TERMS COME FROM MULTIPOLAR APPROACH... AND FACTORS... !!
        exchen_hf=ZERO
        do i=1,nat
          exch_hf(i,i)=exch_hf(i,i)/TWO
          exchen_hf=exchen_hf+exch_hf(i,i)
          do j=i+1,nat
            if(bo(i,j).lt.threbod) exch_hf(i,j)=Excmp(i,j)
            exch_hf(j,i)=exch_hf(i,j)
            exchen_hf=exchen_hf+exch_hf(i,j)
          end do
        end do
        write(*,*) " ----------------------------------------- "
        write(*,*) "  HARTREE-FOCK-TYPE EXCHANGE ENERGY TERMS  "
        write(*,*) " ----------------------------------------- "
        write(*,*) " "
        call MPRINT2(exch_hf,nat,maxat)
        write(*,'(2x,a29,x,f14.7)') "Hartree-Fock exchange energy:",exchen_hf
        write(*,*) " "
      end if
      end if !MMO- non-analytical skip ends here

!! CHECKING THE ACCURACY OF THE TWO ELECTRON PART !!
!! HF ONLY !!
      if(ihf.eq.1) then
        evee=coulen+exchen_hf
        write(*,'(2x,a39,x,f14.7)') "Total two-electron part (coul+exch_hf):",evee
        if(evee0.ne.ZERO) then
          twoelerr=(evee-evee0)*tokcal
          write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",twoelerr
          write(*,*) " "
        else
          evee0=evee    
        end if

!! DFT CASE !!
      else
        exchen=ZERO
        do i=1,nat
          exch(i,i)=exch(i,i)+xmix*exch_hf(i,i)
          exchen=exchen+exch(i,i)
          do j=i+1,nat
            exch(i,j)=exch(i,j)+xmix*exch_hf(i,j)
            exch(j,i)=exch(i,j)
            exchen=exchen+exch(i,j)
          end do
        end do
        if(idoex.eq.1) then
          write(*,*) " ------------------------ "
          write(*,*) "  HYBRID KS-DFT XC TERMS  "
          write(*,*) " ------------------------ "
          write(*,*) " "
          call MPRINT2(exch,nat,maxat)
          write(*,'(2x,a34,x,f14.7)') "Total exchange-correlation energy:",exchen
          write(*,*) " "
        end if
        evee=coulen+exchen
        write(*,'(2x,a37,x,f14.7)') "KS-DFT electron-electron energy (au):",evee           
        if(evee0.ne.ZERO) then
           twoelerr=(evee-evee0)*tokcal
           write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",twoelerr
        else
          evee0=evee
        end if
        write(*,*) " "
      end if

!! ZERO ERROR STRATEGY FOR THE VEE PART !!
      if(ianalytical.eq.1) twoeltoler=232000
      write(*,'(2x,a55,x,f8.2)') "Max error accepted on the two-electron part (kcal/mol):",twoeltoler
      write(*,*) " "

      if(abs(twoelerr).gt.twoeltoler) then

!! NOW ROTATED GRID, ANGLE CONTROLLED BY # GRID SECTION !!
!! CHECK diatXC PAPER FOR OPTIMIZED VALUES: phb=0.162d0 and later 0.182d0 for 40 146 !!
        phb=phb22
        pha=ZERO 
        write(*,'(2x,a20,x,f10.6,x,f10.6)') "Rotating for angles:",pha,phb
        nat0=nat
        call prenumint(ndim,itotps,nat0,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)

!! TO MOs !!
        do k=1,itotps
          xtot=ZERO 
          xtotb=ZERO 
          do j=1,nalf
            xx=ZERO
            xxb=ZERO
            do i=1,igr 
              xx=xx+c(i,j)*chppha(k,i) 
              if(j.le.nb) xxb=xxb+cb(i,j)*chppha(k,i) 
            end do
            chp2phas(j,k)=xx
            if(j.le.nb) then
              chp2phabs(j,k)=xxb
              xtotb=xtotb+xxb*xxb
            end if 
            xtot=xtot+xx*xx
          end do
          rhopha(k)=xtot+xtotb 
        end do

!! INTERPOLATING ONLY ATOMIC TERMS !!
        do icenter=1,nat
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            x0=wp(ifut)*omp2(ifut,icenter)
            dx0=pcoord(ifut,1)
            dy0=pcoord(ifut,2)
            dz0=pcoord(ifut,3)
            f3=ZERO
            do jfut=iatps*(icenter-1)+1,iatps*icenter
              x1=wppha(jfut)*omp2pha(jfut,icenter)
              dx1=pcoordpha(jfut,1)
              dy1=pcoordpha(jfut,2)
              dz1=pcoordpha(jfut,3)
              dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
              if(dist.gt.1.0d-12) f3=f3+(rho(ifut)*rhopha(jfut))*x1*x0/dist !! MG: Could be controlled using the thr2 variable !!
            end do
            coul0(icenter,3)=coul0(icenter,3)+f3/TWO
          end do
        end do

!! EXCHANGE PART !!
        if(idoex.eq.1) then

!! AS BEFORE, MORE CONVOLUTED LOOP STRUCTURE !!
          exch_hfk=ZERO
          !DIR$ NOPARALLEL
          do icenter=1,nat
            ioffset=(icenter-1)*iatps
            istart=0
            iend=0
            !DIR$ NOPARALLEL
            do ik=1,(ithreads-1)
              istart(ik)=ioffset+((ik-1)*ispace)+1
              iend(ik)=ioffset+(ik*ispace)
            end do
            istart(ithreads)=iend(ithreads-1)+1
            iend(ithreads)=icenter*iatps
  
  !! AGAIN THE TWO CALLS ARE CRUCIAL !!
            call omp_set_dynamic(.false.)
            call omp_set_num_threads(ithreads)
            !DIR$ PARALLEL
            do ik=1,ithreads
              do ifut=istart(ik),iend(ik)
                x0=wp(ifut)*omp2(ifut,icenter)
                dx0=pcoord(ifut,1)
                dy0=pcoord(ifut,2)
                dz0=pcoord(ifut,3)
                f3k(ik)=ZERO
                do jfut=iatps*(icenter-1)+1,iatps*icenter
                  x1=wppha(jfut)*omp2pha(jfut,icenter)
                  dx1=pcoordpha(jfut,1)
                  dy1=pcoordpha(jfut,2)
                  dz1=pcoordpha(jfut,3)
                  dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                  if(dist.gt.1.0d-8) then !! MG: Could be controlled using the thr2 variable. In fact I don't remember why was 10^-8 !!
                    do i=1,nalf-1
                      f2=x0*chp2s(i,ifut)*chp2s(i,ifut)
                      f3k(ik)=f3k(ik)+f2*x1*chp2phas(i,jfut)*chp2phas(i,jfut)/dist
                      if(i.le.nb) then
                        f2b=x0*chp2bs(i,ifut)*chp2bs(i,ifut)
                        f3k(ik)=f3k(ik)+f2b*x1*chp2phabs(i,jfut)*chp2phabs(i,jfut)/dist
                      end if
                      do j=i+1,nalf
                        f2=TWO*x0*chp2s(i,ifut)*chp2s(j,ifut)
                        f3k(ik)=f3k(ik)+f2*x1*chp2phas(i,jfut)*chp2phas(j,jfut)/dist !! MG: LC-wPBE programmed in version 3.1 !!
                        if(j.le.nb) then
                          f2b=TWO*x0*chp2bs(i,ifut)*chp2bs(j,ifut)
                          f3k(ik)=f3k(ik)+f2b*x1*chp2phabs(i,jfut)*chp2phabs(j,jfut)/dist
                        end if
                      end do
                    end do
                    f2=x0*chp2s(nalf,ifut)*chp2s(nalf,ifut)
                    f3k(ik)=f3k(ik)+f2*x1*chp2phas(nalf,jfut)*chp2phas(nalf,jfut)/dist

!! AGAIN, ONLY WHEN nalf = nb !!
                    if(nalf.eq.nb) then
                      f2b=x0*chp2bs(nb,ifut)*chp2bs(nb,ifut)
                      f3k(ik)=f3k(ik)+f2b*x1*chp2phabs(nb,jfut)*chp2phabs(nb,jfut)/dist
                    end if
                  end if
                end do
                exch_hfk(icenter,ik)=exch_hfk(icenter,ik)+f3k(ik)
              end do
            end do
          end do
  
  !! ADDING THE TERMS INTO THE ORIGINAL coul0 MATRIX !!
          do icenter=1,nat
            do isum=1,ithreads
              coul0(icenter,4)=coul0(icenter,4)-exch_hfk(icenter,isum)/TWO
            end do
          end do
        end if
  
!! CHECKING NEW VALUES !!
        do i=1,nat
          coul0(i,1)=coul(i,i)
          if(idoex.eq.1) coul0(i,2)=exch_hf(i,i)
        end do
        deltaee=ZERO
        do i=1,nat
          deltaee=deltaee+coul0(i,1)-coul0(i,3)
          if(idoex.eq.1)deltaee=deltaee+xmix*(coul0(i,2)-coul0(i,4))
        end do
        deltaee=deltaee*tokcal
        phabest=ONE-(twoelerr/deltaee)
        write(*,'(2x,a25,x,f8.2)') "New error after rotation:",twoelerr-deltaee
        if(twoelerr-deltaee*twoelerr.gt.ZERO) write(*,*) " WARNING: New error with same sign"
        write(*,'(2x,a18,x,f14.7)') "Damping parameter:",phabest !MMO- not dumping !! ...

!! INTERPOLATING ENERGIES AND REPLACING OLD COULOMB AND EXCHANGE TERMS !!
        do i=1,nat
          coul(i,i)=coul0(i,1)*phabest+(ONE-phabest)*coul0(i,3)
          if(idoex.eq.1) then
            exch_hf(i,i)=coul0(i,2)*phabest+(ONE-phabest)*coul0(i,4)
          end if
        end do
        write(*,*) " "

!! PRINTING !!
        write(*,*) " ------------------------------------------------------- "
        write(*,*) "  INTERPOLATED COULOMB (ELECTRON-ELECTRON) ENERGY TERMS  "
        write(*,*) " ------------------------------------------------------- "
        write(*,*) " " 
        call MPRINT2(coul,nat,maxat)
        coulen=ZERO
        do i=1,nat
          do j=i,nat
            coulen=coulen+coul(i,j)
          end do
        end do
        write(*,'(2x,a15,x,f14.7)') "Coulomb energy:",coulen
        write(*,*) " "
        if (idofr.eq.1) then
          line='   FRAGMENT ANALYSIS: Coulomb energy ' 
          call group_by_frag_mat(1,line ,coul)
        end if

        if(idoex.eq.1) then
          write(*,*) " ------------------------------------------------------ "
          write(*,*) "  INTERPOLATED HARTREE-FOCK-TYPE EXCHANGE ENERGY TERMS  "
          write(*,*) " ------------------------------------------------------ "
          write(*,*) " "
          call MPRINT2(exch_hf,nat,maxat)
          exchen_hf=ZERO
          do i=1,nat
            do j=i,nat
              exchen_hf=exchen_hf+exch_hf(i,j)
            end do
          end do
          write(*,'(2x,a29,x,f14.7)') "Hartree-Fock exchange energy:",exchen_hf
          write(*,*) " "

          if(ihf.ne.1) then
            exchen=ZERO
            do i=1,nat
              exch(i,i)=exch(i,i)+xmix*(exch_hf(i,i)-coul0(i,2))
              exchen=exchen+exch(i,i)
              do j=i+1,nat
!MMO- commenting the line below this one (instruccions Marti)
!! MG: THIS IS ALREADY DONE IN THE FIRST ROTATION AND HERE ARE NOT TOUCHED !!
!                exch(i,j)=exch(i,j)+xmix*exch_hf(i,j)
                exch(j,i)=exch(i,j)
                exchen=exchen+exch(i,j)
              end do
            end do

            write(*,*) " ------------------------------------- "
            write(*,*) "  INTERPOLATED HYBRID KS-DFT XC TERMS  "
            write(*,*) " ------------------------------------- "
            write(*,*) " "
            call MPRINT2(exch,nat,maxat)
            write(*,'(2x,a34,x,f14.7)') "Total exchange-correlation energy:",exchen
            write(*,*) " "
            if (idofr.eq.1) then
              line='   FRAGMENT ANALYSIS: Final Exc Decomposition ' 
              call group_by_frag_mat(1,line ,exch)
            end if
          end if
        else    
          exchen=ZERO
          do i=1,nat
            exchen=exchen+exch(i,i)
            do j=i+1,nat
              exchen=exchen+exch(i,j)
            end do
          end do
        end if
 
!! TOTAL TWO ELECTRON PART !!
        if(ihf.eq.0) then
          evee=coulen+exchen
          write(*,'(2x,a32,x,f14.7)') "KS-DFT electron-electron energy:",evee              
          twoelerr=(evee-evee0)*tokcal
          write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",twoelerr
        else
          evee=coulen+exchen_hf
          write(*,'(2x,a39,x,f14.7)') "Total two-electron part (coul+exch_HF):",evee            
          twoelerr=(evee-evee0)*tokcal
          write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",twoelerr
        end if
        write(*,*) " "
      end if

!! END ZERO-ERROR STRATEGY !!

!! FINAL PRINTING !!
      if(ihf.eq.1) then
        xtot=ZERO
        do i=1,nat
          eto(i,i)=eto(i,i)+coul(i,i)+exch_hf(i,i)
          xtot=xtot+eto(i,i)
          do j=i+1,nat
            eto(i,j)=eto(i,j)+coul(i,j)+exch_hf(i,j)
            eto(j,i)=eto(i,j)
            xtot=xtot+eto(i,j)
          end do
        end do
        etot=xtot

        write(*,*) " -------------------------------------------- "
        write(*,*) "  FUZZY ATOMS Hartree-Fock ENERGY COMPONENTS  "
        write(*,*) " -------------------------------------------- "
        write(*,*) " "
        call MPRINT2(eto,nat,maxat)
        write(*,'(2x,a26,x,f14.7)') "Total integrated energy  :",etot  
        write(*,'(2x,a26,x,f14.7)') "Total energy in Fchk file:",escf 
  
!MMO - added ifield loop, to match restricted
        if(ifield.eq.1) then
          write(*,'(2x,a20,x,f14.7)') "Total dipole energy:",edipole
          err=(etot-escf+edipole)
        else
          err=(etot-escf)
        end if
        write(*,'(2x,a23,x,f14.7)') "Integration error (au):",err
        write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",err*tokcal
      else
        xtot=ZERO
        do i=1,nat
          eto(i,i)=eto(i,i)+coul(i,i)+xmix*exch_hf(i,i)
          xtot=xtot+eto(i,i)
          do j=i+1,nat
            eto(i,j)=eto(i,j)+coul(i,j)+xmix*exch_hf(i,j)
            eto(j,i)=eto(i,j)
            xtot=xtot+eto(i,j)
          end do
        end do
        etot=xtot

        write(*,*) " -------------------------------------- "
        write(*,*) "  FUZZY ATOMS KS-DFT ENERGY COMPONENTS  "
        write(*,*) " -------------------------------------- "
        write(*,*) " "
        call MPRINT2(eto,nat,maxat)
        write(*,'(2x,a26,x,f14.7)') "Total KS-DFT energy      :",etot  
        write(*,'(2x,a26,x,f14.7)') "Total energy in Fchk file:",escf

!MMO - added ifield loop, to match restricted
        if(ifield.eq.1) then
          write(*,'(2x,a20,x,f14.7)') "Total dipole energy:",edipole
          err=(etot-escf+edipole)
        else
          err=(etot-escf)
        end if
        write(*,'(2x,a23,x,f14.7)') "Integration error (au):",err
        write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",err*tokcal
      end if
      if (idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Energy Decomposition' 
        call group_by_frag_mat(1,line ,eto)
      end if

      if(ianalytic.eq.0) then
        DEALLOCATE(chp2,chp2b,chp2pha,rhopha,chp2phab)
        DEALLOCATE(wppha,omppha,omp2pha,chppha,pcoordpha)
!! MG: DEALLOCATING THE PARALLEL PART MISSING (I'VE BEEN LAZY)
      else
        DEALLOCATE(chp2,chp2b)
      end if
      end

!! ***** !!

      subroutine polar(itotps,nat0,wp,omp,omp2,pcoord,rho)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in) :: itotps,nat0
      include 'parameter.h'
      parameter (maxbonds=300)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /atomrad/atr(maxat),dist(maxat,maxat)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common/actual/iact,jat,icenter
      common/efield/field(4),edipole
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension wp(itotps),pcoord(itotps,3),rho(itotps)
      dimension omp(itotps),omp2(itotps,nat0)
      dimension dipint(maxat,3)
      dimension eint(maxat),dip(maxat,3),dipct(3),diptot(3)         
c MM
      dimension itracker(maxbonds),jtracker(maxbonds)
      real*8 MAT
      dimension MAT(maxbonds,maxbonds+1),ichi(maxat,maxat)
      dimension csolv(maxbonds),qqnorm(maxat),ct(maxat,3)
      dimension scr(3),iringatoms(20)
      integer unique_qij,B(maxat,maxat)
      logical ilog
c

      idofr   =   iopt(40)
      ifield  =   iopt(47)

      ipoints=itotps
      iatps=nang*nrad

      do i=1,3
        dipct(i)=0.0d0
      end do


       print *,' '
       print *,' --------------------------------------'
       print *,'       DIPOLE MOMENT DECOMPOSITION'
       print *,' --------------------------------------'
       print *,' '
       xtot=0.0d0
       do icenter=1,nat
        x=0.0d0
        y=0.0d0
        z=0.0d0
        zztop=0.0d0
        do ifut=iatps*(icenter-1)+1,iatps*icenter
         distx=pcoord(ifut,1)
         disty=pcoord(ifut,2)
         distz=pcoord(ifut,3)
         x=x+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*distx
         y=y+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*disty
         z=z+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*distz
         zztop=zztop+wp(ifut)*omp2(ifut,icenter)*rho(ifut)
        end do
        dip(icenter,1)=-x+zn(icenter)*coord(1,icenter)    
        dip(icenter,2)=-y+zn(icenter)*coord(2,icenter)    
        dip(icenter,3)=-z+zn(icenter)*coord(3,icenter)    
        dipint(icenter,1)=-x+zztop*coord(1,icenter)    
        dipint(icenter,2)=-y+zztop*coord(2,icenter)    
        dipint(icenter,3)=-z+zztop*coord(3,icenter)    
        do j=1,3
         dipct(j)=dipct(j)+(zn(icenter)-zztop)*coord(j,icenter)    
        end do
       enddo
       write(*,*) ' '
       write(*,*) 'Origin-dependent atomic dipole moments'
       call mprint_nlop(dip,nat,maxat)
       write(*,*) ' '
       write(*,*) 'Total dipole moment: '
       do j=1,3  
        diptot(j)=0.0d0
        do i=1,nat
         diptot(j)=diptot(j)+dip(i,j)
        end do
       end do
       write(*,*) '               X                   Y                      Z' 
       write(*,'(3f24.15)') (diptot(i),i=1,3)
       write(*,*) ' '
C
       write(*,*) ' '
       write(*,*) 'Intrinsic atomic dipole moments'
       call mprint_nlop(dipint,nat,maxat)
       write(*,*) ' '
       write(*,*) 'Total intrinsic dipole moment: '
       do j=1,3  
        diptot(j)=0.0d0
        do i=1,nat
         diptot(j)=diptot(j)+dipint(i,j)
        end do
       end do
       write(*,*) '               X                   Y                      Z' 
       write(*,'(3f24.15)') (diptot(i),i=1,3)
       write(*,*) ' '
C
C Generalized Keith algorithm to decompose the charge-transfer term
C Code by MMontilla

C Recalculating geometrical connectivities 
C in case TFVC is not used

       do iatom=1,nat-1
         do jatom=iatom+1,nat
           ichi(iatom,jatom)=1
           rijx2=(coord(1,jatom)+coord(1,iatom))/2.0d0
           rijy2=(coord(2,jatom)+coord(2,iatom))/2.0d0
           rijz2=(coord(3,jatom)+coord(3,iatom))/2.0d0
           disth=dist(iatom,jatom)/2.0d0
           do ii=1,nat
             if(ii.ne.iatom.and.ii.ne.jatom) then
               xdist=sqrt((rijx2-coord(1,ii))**2+(rijy2-coord(2,ii))**2+(rijz2-coord(3,ii))**2)
               if(xdist.lt.disth) then
                 ichi(iatom,jatom)=0
                 go to 123
               end if
             end if
           end do
123        continue
         end do
       end do

       do i=1,nat
         do j=i,nat
           B(i,j)=ichi(i,j)
           B(j,i)=-B(i,j)
         end do
       end do

C Finding unique q_ij contributions, aka Q(A|B) contributions.
       unique_qij=0
       do i=1,Nat
         do j=i+1,Nat
           unique_qij=unique_qij+B(i,j)
         end do
       end do
       write(*,*) 'Unique Q(A|B) contributions found:',unique_qij

C Initializing coeff. matrix.
       do i=1,unique_qij
         do j=1,unique_qij
           MAT(i,j)=0.0d0
         end do
       end do

C Checking for lack of equations (unique_qij=number of unknowns;nat-1=number of independent equations
C ;excess unknowns require additional equations)
       iexcess=unique_qij-nat+1
       if(iexcess.ne.0) then 
         write(*,*) 'WARNING, CYCLIC CONNECTIVITY.'
         write(*,*) 'Will need ',iexcess,'additional equations.'
         write(*,*) ' '

         open(45,FILE="ichi.inp",status="old",iostat=ii)
         if(ii.ne.0) then
           write(*,*) 'Printing connectivity matrix:'
           do i=1,nat
             write(*,'(40i2)') (b(i,j),j=1,nat)
           end do
           stop 'error in polar, ichi.inp file is needed'
         end if
         read(45,'(a80)') 
         read(45,*) N

C Introducing (-N) modifications on the connectiviy matrix.
         if(N.lt.0) then
           write(*,*) 'Removing ',abs(N),' bonds.'
           do k=1,abs(N)
             read(45,*) i,j
             B(i,j)=0
             B(j,i)=0
             unique_qij=unique_qij-1
             write(*,*) 'Removed bond between atoms',i,j            
           end do
         end if
       end if

C Labeling the unique q_ij, or Q(A|B), contributions.
       l=1
       do i=1,Nat
         do j=i+1,Nat
           if(B(i,j).eq.1.0) then
             itracker(l)=i
             jtracker(l)=j
             write(*,*) 'Tracker pair ',l,':',i,j
             l=l+1
           end if
         end do
       end do
       if(l-1.ne.unique_qij) stop 'problem in polar'
C End of modifications. Final connectiviy matrix.

C Additional equations for rings if needed. The same iexcess check is performed here,
C  because the loop is never needed if iexcess is 0.
       iextra=0
       if(iexcess.ne.0) then 
         read(45,*) N
         if(N.ne.0) then
           write(*,*) 'Reading ',N,'additional ring equations.'
           do ii=1,N
             read(45,*) iringsize
             read(45,*) (iringatoms(k),k=1,iringsize)
             iextra=iextra+1
             do i=1,iringsize-1
               j=i+1
               do m=1,unique_qij
                 if(iringatoms(i).eq.itracker(m).and.iringatoms(j).eq.jtracker(m)) then
                   Mat(iextra,m)=1.0d0
                 else if(iringatoms(j).eq.itracker(m).and.iringatoms(i).eq.jtracker(m)) then
                   Mat(iextra,m)=-1.0d0
                 end if
               end do
             end do
C Particular case (connectivity between last and first atoms in the ring). 
             i=iringsize
             j=1
             do m=1,unique_qij
               if(iringatoms(i).eq.itracker(m).and.iringatoms(j).eq.jtracker(m)) then
                 Mat(iextra,m)=1.0d0
               else if(iringatoms(j).eq.itracker(m).and.iringatoms(i).eq.jtracker(m)) then
                 Mat(iextra,m)=-1.0d0
               end if
             end do
C Condition for cycles. 
             Mat(iextra,unique_qij+1)=0.0d0
           end do
         end if
       end if
C End reading.
121    if(iexcess.ne.0) close(45)

C Building the matrix of coefficients.
       do i=1,Nat
         do j=1,Nat
           if(B(i,j).eq.1.or.B(i,j).eq.-1) then
             do m=1,unique_qij
               if(i.eq.itracker(m).and.j.eq.jtracker(m).or.j.eq.itracker(m).and.i.eq.jtracker(m)) then
                 Mat(i+iextra,m)=B(i,j)
               end if
             end do
           end if
         end do
       end do

       if(iextra+nat-1.ne.unique_qij) stop 'problem building coeff matrix in polar'

c renormalizing atomic populations
        xx=0.0d0
        do i=1,nat
         xx=xx+qat(i,1)
        end do
        do i=1,nat
         qqnorm(i)=zn(i)-qat(i,1)*(nalf+nb)/xx
        end do
        
        do i=1,nat       
          Mat(i+iextra,unique_qij+1)=qqnorm(i)
        end do

        write(*,*) 'Echoing matrix of coefficients.'
        do i=1,Nat
            write(*,'(8f10.6)') (Mat(i,j),j=1,unique_qij+1)
        end do

        call SOLVESYSTEM(unique_qij,maxbonds,MAT,csolv,ilog)
        if(ilog.eqv..false.) stop 'error in polar, CT decomposition'


C Echoing the solutions of the system of equations, for testing purposes.
c        do i=1,unique_qij
c            write(*,*) 'Unkown number',i,'is',csolv(i)
c        end do

C Final calculations.
       do i=1,nat
         do j=1,3
           CT(i,j)=0.0d0
           scr(j)=0.0d0
         end do
       end do

       do i=1,Nat
         do j=1,Nat
           if(B(i,j).NE.0) then
             Rabx=(coord(1,i)-coord(1,j))/2.0d0
             Raby=(coord(2,i)-coord(2,j))/2.0d0
             Rabz=(coord(3,i)-coord(3,j))/2.0d0
             do m=1,unique_qij
               if(i.eq.itracker(m).and.j.eq.jtracker(m).or.j.eq.itracker(m).and.i.eq.jtracker(m)) l=m
             end do
             CT(i,1)=CT(i,1)+B(i,j)*csolv(l)*Rabx
             CT(i,2)=CT(i,2)+B(i,j)*csolv(l)*Raby
             CT(i,3)=CT(i,3)+B(i,j)*csolv(l)*Rabz
           end if
         end do
       end do

C Echoing final Keith results.
       write(*,*) ' '
       write(*,*) 'Charge-transfer atomic dipole moments (renormalized charges):'
       call mprint_nlop(ct,nat,maxat)
       write(*,*) ' '
       write(*,*) 'Total CT dipole moment (renormalized charges):'
        do i=1,Nat
         do j=1,3
          scr(j)=scr(j)+CT(i,j)
         end do
        end do
        write(*,'(3f24.15)') (scr(i),i=1,3)
       write(*,*) '  ' 
       write(*,*) 'Total CT dipole moment (unscaled charges):'
       write(*,*) '               X                   Y                      Z' 
       write(*,'(3f24.15)') (dipct(i),i=1,3)

C PSS
       write(*,*) ' '
       write(*,*) 'Total dipole moment:'
       write(*,*) '               X                   Y                      Z' 
       write(*,'(3f24.15)') (dipct(i)+diptot(i),i=1,3)
       write(*,*) ' '

c  energetics 
       if(ifield.eq.1) then
       dipen=0.0d0
       do icenter=1,nat
        eint(icenter)=dip(icenter,1)*field(2)+dip(icenter,2)*field(3)+dip(icenter,3)*field(4)  
        dipen=dipen+eint(icenter)
       end do
       cten=cten+dipct(1)*field(2)+dipct(2)*field(3)+dipct(3)*field(4)

       write(*,*) 'Intrinsic atomic dipole energy terms :'
       do i=1,nat
        write(*,'(a6,i4,f22.15)') 'Atom: ',i,eint(i)
       end do
       write(*,*) ' '
       write(*,*) 'Total Intrinsic dipole energy term :',dipen

       write(*,*) 'Total Charge-transfer  energy term :',cten
       edipole=cten+dipen
       write(*,*) 'Total dipole energy term :',edipole

       end if

       end

!! ***** !!

      subroutine multipolar(norb,itotps,wp,omp2,pcoord,fij,xocc,Excmp)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 muamub,muar,mubr,muaqbr,mubqar
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(100)

      dimension :: Excmp(maxat,maxat),wp(itotps),pcoord(itotps,3)
      dimension :: fij(itotps,norb*(norb+1)/2),omp2(itotps,nat)
      dimension :: xocc(norb*(norb+1)/2,norb*(norb+1)/2)

      allocatable :: dip(:,:),quadp(:,:,:),rvect(:),sij(:),atdist(:,:) 
      allocatable :: Excmp1(:,:),Excmp2(:,:),Excmp3(:,:),Excmp4(:,:)
      allocatable :: Excmp5(:,:),Excmp6(:,:)

      iphf  = Iopt(65)
      iatps = nang*nrad
      ALLOCATE(dip(nat,3),quadp(nat,3,3),rvect(3),sij(nat),atdist(nat,nat))
      ALLOCATE(Excmp1(nat,nat),Excmp2(nat,nat),Excmp3(nat,nat))
      ALLOCATE(Excmp4(nat,nat),Excmp5(nat,nat),Excmp6(nat,nat))

      do iat=1,nat
        atdist(iat,iat)=ZERO
        do jat=iat+1,nat
          xx=ZERO
          do ii=1,3
            rvect(ii)=coord(ii,jat)-coord(ii,iat)
            xx=xx+rvect(ii)*rvect(ii)
          end do
          atdist(iat,jat)=dsqrt(xx)
          atdist(jat,iat)=atdist(iat,jat)
          Excmp(iat,jat)=ZERO
          Excmp1(iat,jat)=ZERO
          Excmp2(iat,jat)=ZERO
          Excmp3(iat,jat)=ZERO
          Excmp4(iat,jat)=ZERO
          Excmp5(iat,jat)=ZERO
          Excmp6(iat,jat)=ZERO
        end do
      end do

!! LOOP OVER MO PAIRS !!

      irun=0
      do i=1,norb
      do j=i,norb

      irun=irun+1
      ffact=xocc(irun,irun) 

!! LOOP OVER ATOMS !!

      do icenter=1,nat

!! QUADRUPOLE TERMS !!

        xx=ZERO
        yy=ZERO
        zz=ZERO
        xy=ZERO
        xz=ZERO
        yz=ZERO

!! DIPOLE TERMS !!

        x=ZERO
        y=ZERO
        z=ZERO

!! CHARGE TERM !!

        sij(icenter)=ZERO
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          distx=pcoord(ifut,1)-coord(1,icenter)
          disty=pcoord(ifut,2)-coord(2,icenter)
          distz=pcoord(ifut,3)-coord(3,icenter)
          wccij=wp(ifut)*omp2(ifut,icenter)*fij(ifut,irun)
          xx=xx+wccij*distx*distx
          yy=yy+wccij*disty*disty
          zz=zz+wccij*distz*distz
          xy=xy+wccij*distx*disty
          xz=xz+wccij*distx*distz
          yz=yz+wccij*disty*distz
          x=x+wccij*distx
          y=y+wccij*disty
          z=z+wccij*distz
          sij(icenter)=sij(icenter)+wccij
        end do
        dip(icenter,1)=x
        dip(icenter,2)=y
        dip(icenter,3)=z
        quadp(icenter,1,1)=xx-(yy+zz)/TWO
        quadp(icenter,1,2)=(THREE/TWO)*xy
        quadp(icenter,2,1)= quadp(icenter,1,2)
        quadp(icenter,1,3)=(THREE/TWO)*xz
        quadp(icenter,3,1)= quadp(icenter,1,3)
        quadp(icenter,2,2)=yy-(xx+zz)/TWO
        quadp(icenter,2,3)=(THREE/TWO)*yz
        quadp(icenter,3,2)= quadp(icenter,2,3)
        quadp(icenter,3,3)=zz-(yy+xx)/TWO
      end do

!! iat == A, jat == B !!

      do iat=1,nat
        do jat=iat+1,nat
          xx=atdist(iat,jat)
          do ii=1,3
            rvect(ii)=coord(ii,jat)-coord(ii,iat)
          end do

!! CHARGE-CHARGE !!

          Excmp1(iat,jat)=Excmp1(iat,jat)+ffact*sij(iat)*sij(jat)/xx

!! CHARGE-DIPOLE !!

          muamub=ZERO
          muar=ZERO
          mubr=ZERO
          do ii=1,3
            muamub=muamub+dip(iat,ii)*dip(jat,ii)
            muar=muar+dip(iat,ii)*rvect(ii)
            mubr=mubr+dip(jat,ii)*rvect(ii)
          end do
          Excmp2(iat,jat)=Excmp2(iat,jat)+ffact*(muar*sij(jat)-mubr*sij(iat))/(xx**THREE)

!! DIPOLE-DIPOLE !!

          Excmp3(iat,jat)=Excmp3(iat,jat)-ffact*(THREE*muar*mubr/(xx**FIVE)-muamub/(xx**THREE))

!! CHARGE-QUADRUPOLE !!

          rqar=ZERO
          rqbr=ZERO
          do ii=1,3
            do jj=1,3
              rqar=rqar+rvect(ii)*quadp(iat,ii,jj)*rvect(jj)
              rqbr=rqbr+rvect(ii)*quadp(jat,ii,jj)*rvect(jj)
            end do
          end do
          Excmp4(iat,jat)=Excmp4(iat,jat)+ffact*(rqar*sij(jat)+rqbr*sij(iat))/(xx**FIVE)

!! DIPOLE-QUADRUPOLE !!

          muaqbr=ZERO
          mubqar=ZERO
          do ii=1,3
            do jj=1,3
              muaqbr=muaqbr+rvect(ii)*quadp(jat,ii,jj)*dip(iat,jj)
              mubqar=mubqar+rvect(ii)*quadp(iat,ii,jj)*dip(jat,jj)
            end do
          end do
          xm=-FIVE*(mubr*rqar-muar*rqbr)+TWO*xx*xx*(mubqar-muaqbr)
          Excmp5(iat,jat)=Excmp5(iat,jat)+ffact*xm/(xx**7.0d0)

!! QUADRUPOLE-QUADRUPOLE !!

          qaqb=ZERO
          rqaqbr=ZERO
          do ii=1,3
            do jj=1,3
              qaqb=qaqb+quadp(jat,ii,jj)*quadp(iat,ii,jj)
              do kk=1,3
                rqaqbr=rqaqbr+rvect(ii)*quadp(jat,ii,kk)*quadp(iat,kk,jj)*rvect(jj)
              end do
            end do
          end do
          xm=(35.0d0/THREE)*(rqar*rqbr)/(xx**9.0d0)+(TWO/THREE)*qaqb/(xx**FIVE)-(60.0/9.0d0)*rqaqbr/(xx**7.0d0)
          Excmp6(iat,jat)=Excmp6(iat,jat)+ffact*xm
        end do
      end do

!! END LOOP OVER MO PAIRS !!

      end do
      end do

!! COMPLETING MATRICES !!

      do i=1,nat
        do j=i+1,nat
          Excmp1(j,i)=Excmp1(i,j)
          Excmp2(j,i)=Excmp2(i,j)
          Excmp3(j,i)=Excmp3(i,j)
          Excmp4(j,i)=Excmp4(i,j)
          Excmp5(j,i)=Excmp5(i,j)
          Excmp6(j,i)=Excmp6(i,j)
          Excmp(i,j)=Excmp1(i,j)+Excmp2(i,j)+Excmp3(i,j)+Excmp4(i,j)+Excmp5(i,j)+Excmp6(i,j)
          Excmp(j,i)=Excmp(i,j)
        end do
      end do
      DEALLOCATE(dip,quadp,rvect,sij,atdist,Excmp1,Excmp2,Excmp3,Excmp4,Excmp5,Excmp6)

      end

!! ***** !!

      subroutine calc_coul(itotps,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,rho,rhopha,coul)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      character*80 line
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common /modgrid2/thr3 !! We can make that the distance threshold is controlled using this !!

      dimension :: wp(itotps),wppha(itotps),pcoord(itotps,3),pcoordpha(itotps,3)
      dimension :: omp2(itotps,nat),omp2pha(itotps,nat),rho(itotps),rhopha(itotps)
      dimension :: coul(maxat,maxat)

      iatps = nrad*nang
      idofr = Iopt(40)

      f2=ZERO
      do icenter=1,nat
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          x0=wp(ifut)*omp2(ifut,icenter)
          dx0=pcoord(ifut,1)
          dy0=pcoord(ifut,2)
          dz0=pcoord(ifut,3)

!! SAME CENTER !!
          f3=ZERO
          do jfut=iatps*(icenter-1)+1,iatps*icenter
            x1=wppha(jfut)*omp2pha(jfut,icenter)
            dx1=pcoordpha(jfut,1)
            dy1=pcoordpha(jfut,2)
            dz1=pcoordpha(jfut,3)
            dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
            if(dist.gt.1.0d-12) f3=f3+rho(ifut)*rhopha(jfut)*x1*x0/dist !! MG: Could be controlled using the thr2 variable !!
          end do
          coul(icenter,icenter)=coul(icenter,icenter)+f3

!! PAIRS OF CENTERS !!
          do jcenter=icenter+1,nat
            f3=ZERO
            do jfut=iatps*(jcenter-1)+1,iatps*jcenter
              x1=wp(jfut)*omp2(jfut,jcenter)
              dx1=pcoord(jfut,1)
              dy1=pcoord(jfut,2)
              dz1=pcoord(jfut,3)
              dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
              if(dist.gt.1.0d-12 ) f3=f3+rho(ifut)*rho(jfut)*x1*x0/dist !! MG: Could be controlled using the thr2 variable !!
            end do
            coul(icenter,jcenter)=coul(icenter,jcenter)+f3
          end do
        end do
      end do

      coulen=ZERO
      do i=1,nat
        coul(i,i)=coul(i,i)/TWO
        coulen=coulen+coul(i,i)
        do j=i+1,nat
          coul(j,i)=coul(i,j)
          coulen=coulen+coul(i,j)
        end do
      end do

!! PRINTING !!
      write(*,*) " "
      write(*,*) " ------------------------------------------ "
      write(*,*) "  COULOMB (ELECTRON-ELECTRON) ENERGY TERMS  "
      write(*,*) " ------------------------------------------ "
      write(*,*) " "
      call MPRINT2(coul,nat,maxat)
      write(*,'(2x,a15,x,f14.7)') "Coulomb energy:",coulen
      write(*,*) " "
      if(idofr.eq.1) then
        line=' FRAGMENT ANALYSIS: Coulomb energy '
        call group_by_frag_mat(1,line,coul)
      end if

      end

!! ***** !!

