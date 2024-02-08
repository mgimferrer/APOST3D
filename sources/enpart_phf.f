      
!! ***** !!

      subroutine numint_one_rphf(ndim,itotps,wp,omp2,pcoord,chp,rho,eto)

      use ao_matrices
      use integration_grid

      implicit real*8(a-h,o-z)

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /coord/coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common /energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common /energ0/ekin0,eelnuc0,evee0,etot0
      common /efield/field(4),edipole
      
      integer, intent(in) :: ndim,itotps
      character*80 line
      
      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3)
      dimension omp2(itotps,nat),rho(itotps)
      dimension epa(maxat,maxat),ekin(maxat,maxat),dip(maxat,3),enuc(maxat,maxat)
      dimension eto(maxat,maxat)
      dimension eecp(maxat),ect(maxat,maxat)
      dimension edip(maxat,maxat)

      allocatable :: chp2(:,:),scr(:),chpaux(:,:)

      idofr  = iopt(40)
      ifield = iopt(47)
      iatps  = nang*nrad

      write(*,*) " "
      write(*,*) " ------------------------------------- "
      write(*,*) "  POST-HARTREE-FOCK ONE-ELECTRON PART  "
      write(*,*) " ------------------------------------- "
      write(*,*) " "

!! CALCULATING NUMBER OF CORE+ACTIVE ORBITALS !!
      imaxorb=0
      do ii=1,igr
        if(occ_no(ii,ii).gt.1.0d-8) imaxorb=imaxorb+1
      end do 
      if(2*imaxorb.ne.nspinorb) stop 'Inconsistency in orbital space'

      ALLOCATE(chp2(itotps,norb),scr(itotps))

!! CALCULATING THE DENSITY FROM NATURAL ORBITALS !!
! in principle no need. rho() has already been computed from P-matrix
! but just in case
      do ii=1,itotps
        xx=ZERO
        do jj=1,norb
          xx1=ZERO
          do kk=1,igr
            xx1=xx1+c_no(kk,jj)*chp(ii,kk)
          end do 
          chp2(ii,jj)=xx1
          xx=xx+xx1*xx1*occ_no(jj,jj)
        end do 
        rho(ii)=xx
      end do 

!! ZEROING MATRICES !!
      epa=ZERO
      ekin=ZERO
      enuc=ZERO
      eto=ZERO
      ect=ZERO

!! IN CASE OF HAVING PSEUDOPOTENTIAL !!
      call mga_misc(iecp,eecp)
      if(iecp.eq.1) then
        write(*,*) " ADDING ECP ATOMIC ENERGIES TO E-N TERMS "
        do ii=1,nat
          epa(ii,ii)=eecp(ii)
        end do 
      end if

!! IN CASE OF HAVING AN EXTERNAL ELECTRIC FIELD !!
!! MG: ADAPTED TO THE VERSION IN enpart.f !!
!! MG: NOT CHECKED BY ME, NO IDEA IF WORKS PROPERLY (will ever be used?) !!
      call field_misc(ifield)
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
        do ii=1,nat
          edip(ii,ii)=+dip(ii,1)*field(2)+dip(ii,2)*field(3)+dip(ii,3)*field(4)
          xtot=xtot+ect(ii,ii)+edip(ii,ii)
        end do
        xxdip=xtot
        write(*,'(2x,a37,x,f14.7)') "Electronic dipole moment energy term:",xxdip
        write(*,*) " "
      else
        xxdip=ZERO
      end if

!! ELECTRON-NUCLEI ATTRACTION INTEGRALS !!
      do icenter=1,nat
        do jcenter=1,nat     
          xx=ZERO
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
              xx=xx+zz*wp(ifut)*omp2(ifut,icenter)*rho(ifut)/rr
            end if
          end do
          epa(icenter,jcenter)=epa(icenter,jcenter)-xx
        end do
      end do

      xtot=ZERO
      do ii=1,nat
        xtot=xtot+epa(ii,ii)
        do jj=1,ii-1  
          epa(ii,jj)=(epa(ii,jj)+epa(jj,ii))
          epa(jj,ii)=epa(ii,jj)
          xtot=xtot+epa(ii,jj)
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

!! KINETIC ENERGY !!
      do ii=1,itotps
        xx=ZERO
        do jj=1,norb
          xx1=ZERO
          do kk=1,igr 
            xx1=xx1+c_no(kk,jj)*chpaux(ii,kk) 
          end do
          xx=xx+chp2(ii,jj)*xx1*occ_no(jj,jj)
        end do
        scr(ii)=xx
      end do
      xtot=ZERO
      do icenter=1,nat
        xx=ZERO
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          xx=xx+wp(ifut)*scr(ifut)*omp2(ifut,icenter)
        end do
        ekin(icenter,icenter)=ekin(icenter,icenter)+xx/TWO
        xtot=xtot+xx/TWO
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
      DEALLOCATE(chp2,scr)
      
      end

!! ***** !!

      subroutine numint_two_rphf(ndim,itotps,wp,omp2,pcoord,chp,rho,eto,dm1,dm2)

      use ao_matrices
      use integration_grid

      implicit real*8(a-h,o-z)

      include 'parameter.h'
      parameter (maxpass=5)

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(100)
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/energ0/ekin0,eelnuc0,evee0,etot0
      common/twoel/twoeltoler
      common/efield/field(4),edipole
!! FROM # GRID SECTION !!
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22
      common /modgrid2/thr3

      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3),rho(itotps)
      dimension  omp2(itotps,nat),eto(maxat,maxat),Excmp(maxat,maxat)
      dimension  coul0(maxat,2*maxpass),coul(maxat,maxat)
      dimension  Ex(maxat,maxat),Ec(maxat,maxat),Exc(maxat,maxat)
      dimension  scr(maxat,maxat),xvee(maxat,maxat)
      dimension  dm1(nspinorb,nspinorb),dm2(norb,norb,norb,norb)

      character*80 line

      allocatable :: chp2(:,:),chp3(:,:)
      allocatable :: chp2pha(:,:),rhopha(:),chp3pha(:,:)
      allocatable :: chppha(:,:),wppha(:),pcoordpha(:,:),omppha(:)
      allocatable :: omp2pha(:,:),ibaspointpha(:)
      allocatable :: cumulab(:,:)
      allocatable :: ffa1(:,:),ffb2(:,:)

      idofr    = iopt(40)
      ithrebod = iopt(44)
      ifield   = iopt(47)
      ipolar   = iopt(46)
      itop     = iopt(58)
      iecorr   = iopt(59)
      iorca    = iopt(43)
      ipyscf   = iopt(87)
      ispinsep = iopt(88)
      iatps    = nang*nrad
      npass    = 2

      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if
      
!! ZEROING MATRICES !!

      coul=ZERO
      Exc=ZERO
      Ex=ZERO
      Ec=ZERO
      scr=ZERO
      xvee=ZERO
      coul0=ZERO

!! DOING ENERGY PARTITION !!  
      write(*,*) " "
      write(*,*) " ------------------------------------- "
      write(*,*) "  POST-HARTREE-FOCK TWO-ELECTRON PART  "
      write(*,*) " ------------------------------------- "
      write(*,*) " "

!! ROTATED GRID, ANGLE CONTROLLED BY # GRID SECTION !!
!! CHECK diatXC PAPER FOR OPTIMIZED VALUES: phb=0.162d0 and later 0.182d0 for 40 146 !!
      phb=phb12
      pha=ZERO 
      write(*,'(2x,a20,x,f10.6,x,f10.6)') "Rotating for angles:",pha,phb
      ALLOCATE(wppha(itotps),omppha(itotps),omp2pha(itotps,nat))
      ALLOCATE(chppha(itotps,ndim),pcoordpha(itotps,3),ibaspointpha(itotps))
      ALLOCATE(chp2pha(itotps,norb),rhopha(itotps))

      call prenumint(ndim,itotps,nat,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)
      
      DEALLOCATE(omppha,ibaspointpha)
      ALLOCATE(chp2(itotps,norb))
      ALLOCATE(chp3(itotps,norb),chp3pha(itotps,norb))

!! CALCULATING THE DENSITY FROM NATURAL ORBITALS !!
!! HERE CHP2 AND CHP2PHA ARE THE NOs, CHP3 AND CHP3PHA ARE THE MOs !!
!! MOs CONSTRUCTION REQUIRED FOR CORRELATION ENERGY CALCULATION !!
      do ii=1,itotps
        xx=ZERO
        xxpha=ZERO
        do jj=1,norb
          xx1=ZERO
          xx1pha=ZERO
          xx3=ZERO
          xx3pha=ZERO
          do kk=1,igr
            xx1=xx1+c_no(kk,jj)*chp(ii,kk)
            xx1pha=xx1pha+c_no(kk,jj)*chppha(ii,kk)
            xx3=xx3+c(kk,jj)*chp(ii,kk)
            xx3pha=xx3pha+c(kk,jj)*chppha(ii,kk)
          end do
          chp2(ii,jj)=xx1
          chp2pha(ii,jj)=xx1pha
          chp3(ii,jj)=xx3
          chp3pha(ii,jj)=xx3pha
          xx=xx+xx1*xx1*occ_no(jj,jj)
          xxpha=xxpha+xx1pha*xx1pha*occ_no(jj,jj)
        end do
        rho(ii)=xx
        rhopha(ii)=xxpha
      end do

!! COULOMB PART !!
      call cpu_time(xtime)
      call calc_coul(itotps,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,rho,rhopha,coul)
      call cpu_time(xtime2)
      write(*,'(a32,f10.1,a2)')'(Elapsed time :: Coulomb energy ',xtime2-xtime,'s)'
      write(*,*) " "
      xtime=xtime2

!! EXCHANGE-CORRELATION PART !!
      write(*,*) " ------------------------------------------------------------------ "
      write(*,*) "  EVALUATING POST-HARTREE-FOCK-TYPE EXCHANGE-CORRELATION INTEGRALS  "
      write(*,*) " ------------------------------------------------------------------ "
      write(*,*) " "

!! MONADIC DIAGONALIZATION !!
      norb2=norb*(norb+1)/2
      ALLOCATE(cumulab(norb2,norb2),ffa1(itotps,norb2),ffb2(itotps,norb2))
      call monadic_diag(0,itotps,chp3,chp3pha,dm1,dm2,cumulab,ffa1,ffb2)

!! COMPUTING Exc, MULTIPOLAR APPROACH AND COMPLETING THE MATRIX !!
      call numint_two_pmon(itotps,norb2,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,cumulab,ffa1,ffb2,Exc)
      call multipolar(norb,itotps,wp,omp2,pcoord,ffa1,cumulab,Excmp)

      exch_corr=ZERO
      do ii=1,nat
        Exc(ii,ii)=Exc(ii,ii)/TWO
        exch_corr=exch_corr+Exc(ii,ii)
        do jj=ii+1,nat
          if(bo(ii,jj).lt.threbod) Exc(ii,jj)=Excmp(ii,jj)
          Exc(jj,ii)=Exc(ii,jj)
          exch_corr=exch_corr+Exc(ii,jj)
        end do
      end do

!! PRINTING !!
      write(*,*) " ---------------------------------------------------------- "
      write(*,*) "  POST-HARTREE-FOCK-TYPE EXCHANGE-CORRELATION ENERGY TERMS  "
      write(*,*) " ---------------------------------------------------------- "
      write(*,*) " "
      call MPRINT2(Exc,nat,maxat)
      write(*,'(2x,a46,x,f14.7)') "Post-Hartree-Fock exchange-correlation energy:",exch_corr
      write(*,*) " "
      call cpu_time(xtime2)
      write(*,'(a40,f10.1,a2)')'(Elapsed time :: PHF Exc energy ',xtime2-xtime,'s)' 
      write(*,*) " "
      xtime=xtime2
      DEALLOCATE(cumulab,ffa1,ffb2)

!! EXCHANGE ENERGY PART ONLY (IF ASKED IN THE INPUT FILE) !!
!! REPEATING THE PROCESS THAN FROM Exc !!
      if(iecorr.eq.1) then
        write(*,*) " --------------------------------------------------------- "
        write(*,*) "  EVALUATING POST-HARTREE-FOCK-TYPE CORRELATION INTEGRALS  "
        write(*,*) " --------------------------------------------------------- "
        write(*,*) " "
        ALLOCATE(cumulab(norb2,norb2),ffa1(itotps,norb2),ffb2(itotps,norb2))
        call monadic_diag(1,itotps,chp3,chp3pha,dm1,dm2,cumulab,ffa1,ffb2)
        call numint_two_pmon(itotps,norb2,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,cumulab,ffa1,ffb2,Ex)
        call multipolar(norb,itotps,wp,omp2,pcoord,ffa1,cumulab,Excmp)

        xx=ZERO
        xxc=ZERO
        do ii=1,nat
          Ex(ii,ii)=Ex(ii,ii)/TWO
          Ec(ii,ii)=Exc(ii,ii)-Ex(ii,ii)
          xx=xx+Ex(ii,ii)
          xxc=xxc+Ec(ii,ii)
          do jj=ii+1,nat
            if(bo(ii,jj).lt.threbod) Ex(ii,jj)=Excmp(ii,jj)
            Ex(jj,ii)=Ex(ii,jj)
            Ec(ii,jj)=Exc(ii,jj)-Ex(ii,jj)
            Ec(jj,ii)=Ec(ii,jj)
            xx=xx+Ex(ii,jj)
            xxc=xxc+Ec(ii,jj)
          end do
        end do

!! PRINTING !!
!! FIRST EXCHANGE !!
        write(*,*) " ---------------------------------------------- "
        write(*,*) "  POST-HARTREE-FOCK-TYPE EXCHANGE ENERGY TERMS  "
        write(*,*) " ---------------------------------------------- "
        write(*,*) " "
        call MPRINT2(Ex,nat,maxat)
        write(*,'(2x,a34,x,f14.7)') "Post-Hartree-Fock exchange energy:",xx
        write(*,*) " "

!! NOW CORRELATION !!
        write(*,*) " ------------------------------------------------- "
        write(*,*) "  POST-HARTREE-FOCK-TYPE CORRELATION ENERGY TERMS  "
        write(*,*) " ------------------------------------------------- "
        write(*,*) " "
        call MPRINT2(Ec,nat,maxat)
        write(*,'(2x,a37,x,f14.7)') "Post-Hartree-Fock correlation energy:",xxc
        write(*,*) " "

!! TIMINGS... !!
        call cpu_time(xtime2)
        write(*,'(a40,f10.1,a2)')'(Elapsed time :: PHF Ex & Ec energy ',xtime2-xtime,'s)'
        write(*,*) " "
        xtime=xtime2
        DEALLOCATE(cumulab,ffa1,ffb2)
      end if

!! ZERO ERROR STRATEGY START !!
      DEALLOCATE(wppha,omp2pha,chppha,chp2pha,pcoordpha,rhopha,chp3pha)

!! CHECKING ACCURACY OF THE TWO-EL PART !!
!! EXCHANGE-CORRELATION TOGETHER CASE !!
      if(icorr.ne.1) then
        evee=coulen+exch_corr
        write(*,'(2x,a35,x,f14.7)') "Total two-electron part (coul+exc):",evee
        if(evee0.ne.ZERO) then
          twoelerr=(evee-evee0)*tokcal
          write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",twoelerr
          write(*,*) " "
        else
          evee0=evee    
        end if

!! ZERO ERROR STRATEGY FOR VEE PART !!
        write(*,'(2x,a55,x,f8.2)') "Max error accepted on the two-electron part (kcal/mol):",twoeltoler
        write(*,*) " "

        if(abs(twoelerr).gt.twoeltoler) then

!! NEW ROTATED GRID, ANGLE CONTROLLED BY # GRID SECTION !!
!! CHECK diatXC PAPER FOR OPTIMIZED VALUES: phb=0.162d0 and later 0.182d0 for 40 146 !!
        phb=phb22
        pha=ZERO
        write(*,'(2x,a20,x,f10.6,x,f10.6)') "Rotating for angles:",pha,phb
        ALLOCATE(wppha(itotps),omppha(itotps),omp2pha(itotps,nat))
        ALLOCATE(chppha(itotps,ndim),pcoordpha(itotps,3),ibaspointpha(itotps))
        ALLOCATE(chp3pha(itotps,norb),rhopha(itotps))
        call prenumint(ndim,itotps,nat,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)

!! COMPUTING MOs AND DENSITY !!
        do ii=1,itotps
          xxpha=ZERO
          do jj=1,norb
            xx1pha=ZERO
            xx3pha=ZERO
            do kk=1,igr
              xx1pha=xx1pha+c_no(kk,jj)*chppha(ii,kk)
              xx3pha=xx3pha+c(kk,jj)*chppha(ii,kk)
            end do
            chp3pha(ii,jj)=xx3pha
            xxpha=xxpha+xx1pha*xx1pha*occ_no(jj,jj)
          end do
          rhopha(ii)=xxpha
        end do

!! DOING INTEGRATIONS FOR ONLY THE ATOMIC TERMS !!
!! COULOMB PART !!
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
              if(dist.gt.1.0d-12) f3=f3+rho(ifut)*rhopha(jfut)*x1*x0/dist
            end do
            coul0(icenter,3)=coul0(icenter,3)+f3/TWO
          end do
        end do

!! EXCHANGE-CORRELATION PART !!
!! NASTY TRICK: CONTROLLED BY ITHREBOD, BUT RECOMPUTING ffb2 IS REQUIRED !!
        ALLOCATE(cumulab(norb2,norb2),ffa1(itotps,norb2),ffb2(itotps,norb2))
        call monadic_diag(0,itotps,chp3,chp3pha,dm1,dm2,cumulab,ffa1,ffb2)
        ithr=iopt(44)
        iopt(44)=1000000
        call numint_two_pmon(itotps,norb2,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,cumulab,ffa1,ffb2,scr)
        iopt(44)=ithr
        do ii=1,nat
          coul0(ii,4)=scr(ii,ii)/TWO
        end do
        DEALLOCATE(ffa1,ffb2,cumulab)
        DEALLOCATE(wppha,omp2pha,chppha,pcoordpha,rhopha,chp3pha)

!! CHECKING NEW VALUES !!
        do ii=1,nat
          coul0(ii,1)=coul(ii,ii)
          coul0(ii,2)=Exc(ii,ii)
        end do
        deltaee=ZERO
        do ii=1,nat
          deltaee=deltaee+coul0(ii,1)-coul0(ii,3)+coul0(ii,2)-coul0(ii,4)
        end do
        deltaee=deltaee*tokcal
        phabest=ONE-(twoelerr/deltaee)
        write(*,'(2x,a25,x,f8.2)') "New error after rotation:",twoelerr-deltaee
        if(twoelerr-deltaee*twoelerr.gt.ZERO) write(*,*) " WARNING: New error with same sign"
        write(*,'(2x,a18,x,f14.7)') "Damping parameter:",phabest  

!! INTERPOLATING ENERGIES, REPLACING OLD COULOMB AND EXCHANGE-CORRELATION TERMS !!
        do ii=1,nat
          coul(ii,ii)=coul0(ii,1)*phabest+(ONE-phabest)*coul0(ii,3)
          Exc(ii,ii)=coul0(ii,2)*phabest+(ONE-phabest)*coul0(ii,4)
        end do
        write(*,*) " "

!! PRINTING !!
!! FIRST COULOMB !!
        write(*,*) " "
        write(*,*) " ------------------------------------------------------- "
        write(*,*) "  INTERPOLATED COULOMB (ELECTRON-ELECTRON) ENERGY TERMS  "
        write(*,*) " ------------------------------------------------------- "
        write(*,*) " "
        call MPRINT2(coul,nat,maxat)
        coulen=ZERO
        do ii=1,nat
          do jj=ii,nat
            coulen=coulen+coul(ii,jj)
          end do
        end do
        write(*,'(2x,a15,x,f14.7)') "Coulomb energy:",coulen
        write(*,*) " "
        if (idofr.eq.1) then
          line=' FRAGMENT ANALYSIS: Coulomb energy ' 
          call group_by_frag_mat(1,line,coul)
        end if

!! NOW EXCHANGE-CORRELATION !!
        write(*,*) " ----------------------------------------------------------------------- "
        write(*,*) "  INTERPOLATED POST-HARTREE-FOCK-TYPE EXCHANGE-CORRELATION ENERGY TERMS  "
        write(*,*) " ----------------------------------------------------------------------- "
        write(*,*) " "
        call MPRINT2(Exc,nat,maxat)
        exch_corr=ZERO
        do ii=1,nat
          do jj=ii,nat
            exch_corr=exch_corr+Exc(ii,jj)
          end do
        end do
        write(*,'(2x,a46,x,f14.7)') "Post-Hartree-Fock exchange-correlation energy:",exch_corr
        write(*,*) " "
        if(idofr.eq.1) then
          line=' FRAGMENT ANALYSIS: Final Exchange-Correlation Decomposition ' 
          call group_by_frag_mat(1,line,Exc)
        end if

        evee=coulen+exch_corr
        write(*,'(2x,a35,x,f14.7)') "Total two-electron part (coul+exc):",evee
        twoelerr=(evee-evee0)*tokcal
        write(*,'(2x,a29,x,f8.2)') "Integration error (kcal/mol):",twoelerr

!! END ZES !!
        end if
        
!! EXCHANGE AND CORRELATION SEPARATED !!
      else 
        write(*,*) " INFO: ZES IS NOT PREPARED FOR Ex AND Ec SEPARATED, YET "
      end if

!! FINAL PRINTING !!
      xtot=ZERO
      do ii=1,nat
        eto(ii,ii)=eto(ii,ii)+coul(ii,ii)+Exc(ii,ii)
        xtot=xtot+eto(ii,ii)
        do jj=ii+1,nat
          eto(ii,jj)=eto(ii,jj)+coul(ii,jj)+Exc(ii,jj)
          eto(jj,ii)=eto(ii,jj)
          xtot=xtot+eto(ii,jj)
        end do
      end do
      etot=xtot

      write(*,*) " ------------------------------------------------- "
      write(*,*) "  FUZZY ATOMS Post-Hartree-Fock ENERGY COMPONENTS  "
      write(*,*) " ------------------------------------------------- "
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
      if (idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Energy Decomposition' 
        call group_by_frag_mat(1,line,eto)
      end if

!! DEALLOCATING !!
      DEALLOCATE(chp2,chp3)

      end

!! ***** !!

      subroutine monadic_diag(ic,itotps,chp2,chp2pha,dm1,dm2,cumulab,fij,fijb)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass

      dimension :: chp2(itotps,norb),chp2pha(itotps,norb)
      dimension :: dm1(nspinorb,nspinorb),dm2(norb,norb,norb,norb)
      dimension :: cumulab(norb*(norb+1)/2,norb*(norb+1)/2)
      dimension :: fij(itotps,norb*(norb+1)/2),fijb(itotps,norb*(norb+1)/2)

      allocatable :: cumul(:,:,:,:),eigenmat(:,:),fa1(:,:),fb2(:,:)
      allocatable :: sfdm1(:,:),psdm1(:,:) 

!! TRANSFORMATION OF THE DM2 INTO THE CUMMULANT (XC IF ic = 0, X IF ic = 1 AND C IF ic = 2) !!

      ALLOCATE(cumul(norb,norb,norb,norb))
      ALLOCATE(sfdm1(norb,norb),psdm1(norb,norb))

C P AND Ps in MO basis
      do i=1,norb
       do k=1,norb
        sfdm1(i,k)=dm1((i-1)*2+1,(k-1)*2+1)+dm1((i-1)*2+2,(k-1)*2+2)
        psdm1(i,k)=dm1((i-1)*2+1,(k-1)*2+1)-dm1((i-1)*2+2,(k-1)*2+2)
       end do
      end do

      do ii=1,norb
        do jj=1,norb
          do kk=1,norb
            do ll=1,norb
              xx1=sfdm1(ii,jj)*sfdm1(kk,ll)
              xx2=sfdm1(ii,kk)*sfdm1(jj,ll)/TWO
              xx3=psdm1(ii,kk)*psdm1(jj,ll)/TWO
              if(ic.eq.0) then 
                cumul(ii,jj,kk,ll)=dm2(ii,jj,kk,ll)-xx1
              else if(ic.eq.1) then
                cumul(ii,jj,kk,ll)=-xx2-xx3
              else if(ic.eq.2) then
                cumul(ii,jj,kk,ll)=dm2(ii,jj,kk,ll)-xx1+xx2+xx3
              end if 
            end do 
          end do 
        end do 
      end do 

      DEALLOCATE(sfdm1,psdm1)
!! GROUPING TERMS FOR J >/ I AND L >/ K, NECESSARY FOR COMPACTING LATER !!
!! COMPACTING INDEXES FOR THE CUMULANT (I,J --> IJ AND K,L --> KL) !!

      ij=1
      do ii=1,norb
        do jj=ii,norb
          kl=1
          do kk=1,norb
            do ll=kk,norb
              if(jj.ne.ii) cumul(ii,jj,kk,ll)=cumul(ii,jj,kk,ll)+cumul(jj,ii,kk,ll)
              if(ll.ne.kk) cumul(ii,jj,kk,ll)=cumul(ii,jj,kk,ll)+cumul(ii,jj,ll,kk)
              if(jj.ne.ii.and.ll.ne.kk) cumul(ii,jj,kk,ll)=cumul(ii,jj,kk,ll)+cumul(jj,ii,ll,kk)
              cumulab(ij,kl)=cumul(ii,jj,kk,ll)
              kl=kl+1
            end do
          end do
          ij=ij+1
        end do
      end do
      DEALLOCATE(cumul)

!! DIAGONALIZING THE CORRESPONDING CUMULANT !!

      norb2=norb*(norb+1)/2
      ALLOCATE(eigenmat(norb2,norb2))
      call diagonalize(norb2,norb2,cumulab,eigenmat,0)

!! FINISHING INDEX COMPACTATION !!
!! PRODUCT OF ORBITALS FROM FIRST GRID IN fa1, FROM SECOND (ROTATED) IN fb2 !!

      ALLOCATE(fa1(itotps,norb2),fb2(itotps,norb2))
      do ifut=1,itotps
        ij=1
        do ii=1,norb
          do jj=ii,norb
            fa1(ifut,ij)=chp2(ifut,ii)*chp2(ifut,jj)
            fb2(ifut,ij)=chp2pha(ifut,ii)*chp2pha(ifut,jj)
            ij=ij+1
          end do
        end do
      end do  

!! TRANSFORMATION OF fa1 AND fb2 FROM DIAGONALIZATION OF COMPACTED EXCHANGE-CORRELATION CUMULANT !!

      do ifut=1,itotps
        do ii=1,norb2
          xx=ZERO
          xxb=ZERO
          do kk=1,norb2
            xx=xx+eigenmat(kk,ii)*fa1(ifut,kk)
            xxb=xxb+eigenmat(kk,ii)*fb2(ifut,kk)
          end do 
          fij(ifut,ii)=xx
          fijb(ifut,ii)=xxb
        end do
      end do
      DEALLOCATE(fa1,fb2,eigenmat)

      end 

! *****

      subroutine numint_two_pmon(itotps,iorb2,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,cumulab,ffa1,ffb2,ef)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(100)

      dimension :: wp(itotps),wppha(itotps),pcoord(itotps,3),pcoordpha(itotps,3)
      dimension :: omp2(itotps,nat),omp2pha(itotps,nat),cumulab(iorb2,iorb2)
      dimension :: ffa1(itotps,iorb2),ffb2(itotps,iorb2)

      dimension :: ef(maxat,maxat)

      iatps    = nrad*nang
      iterms   = 0
      ithrebod = Iopt(44)

      !thr2=thr3 !! THRESH FOR NUMERICAL INTEGRATION (DISTANCE) !!
      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if
      
      write(*,*) " Two-el integrations for ",iorb2," MO pairs "
      write(*,'(a38,f10.6)') " Threshold for atom pair calculation : ",threbod

!! TWO-EL NUMERICAL INTEGRATION FOR PHF and HF Exc, Ex and Ec !!

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
            if(dist.gt.1.0d-12) then
              xxw=x0*x1/dist
              do ij=1,iorb2
                f3=f3+cumulab(ij,ij)*xxw*ffa1(ifut,ij)*ffb2(jfut,ij)
              end do
            end if
          end do
          ef(icenter,icenter)=ef(icenter,icenter)+f3

!! PAIRS OF CENTERS !!

          do jcenter=icenter+1,nat
            bx0=bo(icenter,jcenter)
            if(bx0.ge.threbod) then
              f3=ZERO
              do jfut=iatps*(jcenter-1)+1,iatps*jcenter
                x1=wp(jfut)*omp2(jfut,jcenter)
                dx1=pcoord(jfut,1)
                dy1=pcoord(jfut,2)
                dz1=pcoord(jfut,3)
                dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                if(dist.gt.1.0d-12 ) then
                  xxw=x0*x1/dist
                  do ij=1,iorb2
                    f3=f3+cumulab(ij,ij)*xxw*ffa1(ifut,ij)*ffa1(jfut,ij)
                  end do
                end if
              end do
              ef(icenter,jcenter)=ef(icenter,jcenter)+f3
            else
              iterms=iterms+1
            end if
          end do
        end do
      end do
      write(*,*) " Skipping diatomic exchange integrals for ",iterms
      write(*,*) " "

      end 

! *****



