      
!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: numint_one_rphf
!! One-electron IQA terms of a correlated (CASSCF/CISD) wavefunction:
!! electron-nuclear attraction, kinetic energy (from the natural
!! orbitals) and nuclear repulsion, printed and added to eto.
!! arguments:
!!   ndim   (in)  number of basis functions
!!   itotps (in)  number of grid points
!!   wp     (in)  integration weights
!!   omp2   (in)  atomic weight of each point for every atom
!!   pcoord (in)  coordinates of the points
!!   chp    (in)  basis functions on the grid
!!   rho    (inout) density on the grid (recomputed from the NOs)
!!   eto    (out) one-electron IQA matrix
!! author: PSalse, MGimf
!! ----------------------------------------------------------------------
      subroutine numint_one_rphf(ndim,itotps,wp,omp2,pcoord,chp,rho,eto)

      use ao_matrices
      use integration_grid
      use output_mod, only: print_value

      implicit real*8(a-h,o-z)

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /coord/coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(200)
      common /energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common /energ0/ekin0,eelnuc0,evee0,etot0
      
      integer, intent(in) :: ndim,itotps
      logical lref
      common /printout/iaccur
      
      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3)
      dimension omp2(itotps,nat),rho(itotps)
      dimension epa(maxat,maxat),ekin(maxat,maxat),enuc(maxat,maxat)
      dimension eto(maxat,maxat)
      dimension eecp(maxat),tkin(maxat)

      allocatable :: chp2(:,:),scr(:),chpaux(:,:)

      idofr  = iopt(40)
      iatps  = nang*nrad

      ndec=7
      if(iaccur.ne.0) ndec=13

!! number of occupied (core + active) natural orbitals !!
      imaxorb=0
      do ii=1,igr
        if(occ_no(ii,ii).gt.1.0d-8) imaxorb=imaxorb+1
      end do 
      if(2*imaxorb.ne.nspinorb) call apost_stop('Inconsistency in orbital space')

      ALLOCATE(chp2(itotps,norb),scr(itotps))

!! density from the natural orbitals (rho is already built from P; it is !!
!! recomputed here from the same orbitals as the energy terms)           !!
!! parallel over grid points: each ii writes only chp2(ii,:) and rho(ii) !!
!$OMP PARALLEL DO PRIVATE(ii,xx,jj,xx1,kk)
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
!$OMP END PARALLEL DO

!! energy matrices !!
      epa=ZERO
      ekin=ZERO
      enuc=ZERO
      eto=ZERO

!! pseudopotential (ECP) energies, added to the E-N one-center terms !!
      call mga_misc(iecp,eecp)
      if(iecp.eq.1) then
        do ii=1,nat
          epa(ii,ii)=eecp(ii)
        end do 
      end if


!! electron-nuclear attraction !!
!! parallel over the (icenter,jcenter) pairs (COLLAPSE(2), constant bounds): !!
!! each pair writes only its own epa(icenter,jcenter)                        !!
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,jcenter,xx,zz,dx0,dy0,dz0,
!$OMP&  ifut,distx,disty,distz,rr)
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
!$OMP END PARALLEL DO

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

!! the matrix, its total and (with reference energies) its error !!
      if(iecp.eq.1) then
        call enpart_matrix('ELECTRON-NUCLEAR ATTRACTION','Includes the'//
     +    ' pseudopotential (ECP) energy of each atom on the diagonal.',
     +    ' ',epa,'Electron-nuclear attraction',eelnuc)
      else
        call enpart_matrix('ELECTRON-NUCLEAR ATTRACTION',' ',' ',epa,
     +    'Electron-nuclear attraction',eelnuc)
      end if
      lref=eelnuc0.ne.ZERO
      errk=ZERO
      if(lref) then
        errk=(eelnuc-eelnuc0)*tokcal
        call print_value('Integration error',errk,2,'kcal/mol')
      else
        eelnuc0=eelnuc
      end if
      call enpart_frag('ELECTRON-NUCLEAR ATTRACTION',epa)
      call json_enpart_term('electron_nuclear','Electron-nuclear '//
     +  'attraction',epa,eelnuc,lref,errk)

!! second derivatives of the AOs on the grid, for the kinetic energy !!
      ALLOCATE(chpaux(itotps,igr))
      call dpoints(chpaux,pcoord)

!! kinetic energy density from the natural orbitals !!
!! parallel over grid points: each ii writes only scr(ii) !!
!$OMP PARALLEL DO PRIVATE(ii,xx,jj,xx1,kk)
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
!$OMP END PARALLEL DO
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

!! one-center terms only: printed as a column per atom !!
      do i=1,nat
        tkin(i)=ekin(i,i)
      end do
      call enpart_atom_table('KINETIC ENERGY','One-center terms only:'//
     +  ' the kinetic energy has no pair terms.',' ',1,tkin,'Kinetic',
     +  ' ',.false.)
      call print_value('Kinetic energy',ekinen,ndec,'au')
      lref=ekin0.ne.ZERO
      errk=ZERO
      if(lref) then
        errk=(ekinen-ekin0)*tokcal
        call print_value('Integration error',errk,2,'kcal/mol')
      else
        ekin0=ekinen
      end if
      if(idofr.eq.1) call group_by_frag_vec(1,'FRAGMENT KINETIC ENERGY',
     +  tkin,'Kinetic',' ')
!! the results file keeps the kinetic energy as a (diagonal) matrix !!
      call json_enpart_term('kinetic','Kinetic energy',ekin,ekinen,lref,
     +  errk)
      DEALLOCATE(chpaux)

!! nuclear repulsion, and the one-electron IQA matrix !!

      erep=ZERO
      do i=1,nat
        eto(i,i)=ekin(i,i)+epa(i,i)
        do j=i+1,nat
          dist=dsqrt((coord(1,i)-coord(1,j))**TWO+(coord(2,i)-coord(2,j))**TWO+(coord(3,i)-coord(3,j))**TWO)
          eto(i,j)=ekin(i,j)+epa(i,j)+zn(i)*zn(j)/dist
          eto(j,i)=eto(i,j)

!! the nuclear-repulsion pair matrix, for the output !!
          enuc(i,j)=zn(i)*zn(j)/dist
          enuc(j,i)=enuc(i,j)
          erep=erep+zn(i)*zn(j)/dist
        end do
      end do

      call enpart_matrix('NUCLEAR-NUCLEAR REPULSION',' ',' ',enuc,
     +  'Nuclear-nuclear repulsion',erep)
      call enpart_frag('NUCLEAR-NUCLEAR REPULSION',enuc)
      call json_enpart_term('nuclear_repulsion','Nuclear-nuclear '//
     +  'repulsion',enuc,erep,.false.,0.0d0)


!! one-electron total, and its reference counterpart !!
      etot=ekinen+eelnuc+erep
      etot0=ekin0+eelnuc0+erep
      DEALLOCATE(chp2,scr)
      
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: numint_two_rphf
!! Two-electron IQA terms of a correlated wavefunction: Coulomb, and the
!! exchange-correlation from the cumulant of the 2-RDM (with CORRELATION
!! also exchange and correlation separately), the zero-error correction,
!! the printing of the final matrices and of the IQA matrix.
!! arguments:
!!   ndim,itotps,wp,omp2,pcoord,chp,rho (in)  grid and density, as in
!!                                       numint_one_rphf
!!   eto    (inout) IQA matrix: one-electron terms in, total out
!!   dm1    (in)  spin-orbital 1-RDM
!!   dm2    (in)  spinless 2-RDM
!! author: PSalse, MGimf
!! ----------------------------------------------------------------------
      subroutine numint_two_rphf(ndim,itotps,wp,omp2,pcoord,chp,rho,eto,dm1,dm2)

      use ao_matrices
      use integration_grid
      use timing_mod
      use output_mod, only: print_value
      use input_options_mod, only: iverbose

      implicit real*8(a-h,o-z)

      include 'parameter.h'
      parameter (maxpass=5)

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(200)
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/energ0/ekin0,eelnuc0,evee0,etot0
      common/twoel/twoeltoler
!! rotation angles of the two-electron grid (# GRID) !!
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22

      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3),rho(itotps)
      dimension  omp2(itotps,nat),eto(maxat,maxat),Excmp(maxat,maxat)
      dimension  coul0(maxat,2*maxpass),coul(maxat,maxat)
      dimension  Ex(maxat,maxat),Ec(maxat,maxat),Exc(maxat,maxat)
      dimension  scr(maxat,maxat)
      dimension  dm1(nspinorb,nspinorb),dm2(norb,norb,norb,norb)

      allocatable :: coulb(:,:),excb(:,:)

      allocatable :: chp2(:,:),chp3(:,:)
      allocatable :: chp2pha(:,:),rhopha(:),chp3pha(:,:)
      allocatable :: chppha(:,:),wppha(:),pcoordpha(:,:),omppha(:)
      allocatable :: omp2pha(:,:),ibaspointpha(:)
      allocatable :: cumulab(:,:)
      allocatable :: ffa1(:,:),ffb2(:,:)

      idofr    = iopt(40)
      ithrebod = iopt(44)
      iecorr   = iopt(59)
      iatps    = nang*nrad

      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if
      
!! energy matrices !!

      coul=ZERO
      Exc=ZERO
      Ex=ZERO
      Ec=ZERO
      scr=ZERO
      coul0=ZERO

!! the atom pairs integrated and those from the multipolar expansion !!
      npairs=0
      nskip=0
      do i=1,nat
        do j=i+1,nat
          if(bo(i,j).ge.threbod) then
            npairs=npairs+1
          else
            nskip=nskip+1
          end if
        end do
      end do
      call print_value('Pairs integrated (exchange-correlation)',npairs)
      call print_value('Pairs from the multipolar expansion',nskip)

!! second, rotated grid for the one-center terms (angles: # GRID) !!
      phb=phb12
      pha=ZERO 
      ALLOCATE(wppha(itotps),omppha(itotps),omp2pha(itotps,nat))
      ALLOCATE(chppha(itotps,ndim),pcoordpha(itotps,3),ibaspointpha(itotps))
      ALLOCATE(chp2pha(itotps,norb),rhopha(itotps))

      call prenumint(ndim,itotps,nat,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)
      
      DEALLOCATE(omppha,ibaspointpha)
      ALLOCATE(chp2(itotps,norb))
      ALLOCATE(chp3(itotps,norb),chp3pha(itotps,norb))

!! on both grids: natural orbitals (chp2, chp2pha) for the density, MOs !!
!! (chp3, chp3pha) for the cumulant, whose RDMs are in the MO basis      !!
!! parallel over grid points: each ii writes only its own row of chp2,  !!
!! chp2pha, chp3, chp3pha and its rho(ii), rhopha(ii)                   !!
!$OMP PARALLEL DO PRIVATE(ii,xx,xxpha,jj,xx1,xx1pha,xx3,xx3pha,kk)
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
!$OMP END PARALLEL DO

!! Coulomb !!
      call cpu_time(xtime)
      call get_wall_time(wxtime)
      call calc_coul(itotps,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,rho,rhopha,coul)
      call cpu_time(xtime2)
      call get_wall_time(wxtime2)
      write(*,'(a)') ''
      call print_timer('ENPART Coulomb',xtime2-xtime,wxtime2-wxtime)
      xtime=xtime2
      wxtime=wxtime2

!! exchange-correlation: diagonalize the cumulant (monadic_diag), !!
!! integrate its pair functions, multipolar estimate for the pairs  !!
!! below THREBOD                                                    !!
      norb2=norb*(norb+1)/2
      ALLOCATE(cumulab(norb2,norb2),ffa1(itotps,norb2),ffb2(itotps,norb2))
      call monadic_diag(0,itotps,chp3,chp3pha,dm1,dm2,cumulab,ffa1,ffb2)

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

      call cpu_time(xtime2)
      call get_wall_time(wxtime2)
      call print_timer('ENPART exchange-correlation',xtime2-xtime,
     +  wxtime2-wxtime)
      xtime=xtime2
      wxtime=wxtime2
      DEALLOCATE(cumulab,ffa1,ffb2)

!! CORRELATION: the exchange part from the exchange cumulant, the same !!
!! way; correlation = exchange-correlation - exchange                  !!
      if(iecorr.eq.1) then
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

        exenx=xx
        exenc=xxc

        call cpu_time(xtime2)
        call get_wall_time(wxtime2)
        call print_timer('ENPART exchange and correlation',xtime2-xtime,
     +    wxtime2-wxtime)
        xtime=xtime2
        wxtime=wxtime2
        DEALLOCATE(cumulab,ffa1,ffb2)
      end if

      DEALLOCATE(wppha,omp2pha,chppha,chp2pha,pcoordpha,rhopha,chp3pha)

!! two-electron accuracy check and zero-error correction, on the total !!
!! exchange-correlation (with CORRELATION, exchange and correlation are !!
!! not corrected: see enpart_phf_mats below)                            !!
      evee=coulen+exch_corr
!! the error stays zero (no correction) without a reference energy !!
      iref=0
      if(evee0.ne.ZERO) iref=1
      twoelerr=ZERO
      if(iref.eq.1) then
        twoelerr=(evee-evee0)*tokcal
      else
        evee0=evee
      end if
      call enpart_ee_head(evee,twoelerr,iref,twoeltoler)
!! first-grid values, kept for the results file !!
      eveefg=evee
      errfg=twoelerr
      err2=ZERO
      phabest=ZERO

!! the correction changes only the one-center Coulomb and xc terms; the !!
!! first-grid matrices are kept for VERBOSE 1                           !!
      izerr=0
      if(iref.eq.1.and.abs(twoelerr).gt.twoeltoler) izerr=1
      if(izerr.eq.1.and.iverbose.eq.1) then
        ALLOCATE(coulb(maxat,maxat),excb(maxat,maxat))
        do ii=1,nat
          do jj=1,nat
            coulb(ii,jj)=coul(ii,jj)
            excb(ii,jj)=Exc(ii,jj)
          end do
        end do
        xcoulb=coulen
        xexcb=exch_corr
      end if

      if(izerr.eq.1) then

!! zero-error correction: the one-center terms again on a grid rotated !!
!! by the second angle of # GRID, then interpolated                       !!
      phb=phb22
      pha=ZERO
      ALLOCATE(wppha(itotps),omppha(itotps),omp2pha(itotps,nat))
      ALLOCATE(chppha(itotps,ndim),pcoordpha(itotps,3),ibaspointpha(itotps))
      ALLOCATE(chp3pha(itotps,norb),rhopha(itotps))
      call prenumint(ndim,itotps,nat,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)

!! MOs and density on the rotated grid !!
!! parallel over grid points: each ii writes only chp3pha(ii,:), rhopha(ii) !!
!$OMP PARALLEL DO PRIVATE(ii,xxpha,jj,xx1pha,xx3pha,kk)
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
!$OMP END PARALLEL DO

!! one-center terms only: Coulomb !!
!! parallel over the points ifut of each atom: f3 is PRIVATE, xacc a  !!
!! reduction added to coul0(icenter,3) after the loop (on one thread   !!
!! the summation order is the serial one)                               !!
      do icenter=1,nat
        xacc=ZERO
!$OMP PARALLEL DO PRIVATE(ifut,x0,dx0,dy0,dz0,f3,jfut,x1,dx1,dy1,dz1,
!$OMP&  dist) REDUCTION(+:xacc)
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
          xacc=xacc+f3/TWO
        end do
!$OMP END PARALLEL DO
        coul0(icenter,3)=coul0(icenter,3)+xacc
      end do

!! one-center terms only: exchange-correlation, with ffb2 recomputed for !!
!! the rotated grid; THREBOD is raised for the call so that every pair   !!
!! is skipped                                                            !!
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

!! the change of the one-center terms between the grids, and the damping !!
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
      err1=twoelerr
      err2=twoelerr-deltaee

!! interpolated one-center Coulomb and exchange-correlation terms !!
      do ii=1,nat
        coul(ii,ii)=coul0(ii,1)*phabest+(ONE-phabest)*coul0(ii,3)
        Exc(ii,ii)=coul0(ii,2)*phabest+(ONE-phabest)*coul0(ii,4)
      end do

!! corrected totals and error !!
      coulen=ZERO
      exch_corr=ZERO
      do ii=1,nat
        do jj=ii,nat
          coulen=coulen+coul(ii,jj)
          exch_corr=exch_corr+Exc(ii,jj)
        end do
      end do
      evee=coulen+exch_corr
      twoelerr=(evee-evee0)*tokcal
      call enpart_ee_tail(err1,err2,phabest,evee,twoelerr,.true.)
      write(*,'(a)') ''
      call cpu_time(xtime2)
      call get_wall_time(wxtime2)
      call print_timer('ENPART zero-error correction',xtime2-xtime,
     +  wxtime2-wxtime)
      xtime=xtime2
      wxtime=wxtime2

      end if

      call json_enpart_ee(eveefg,errfg,iref,twoeltoler,izerr,err2,
     +  phabest,evee,twoelerr)

!! final two-electron matrices; with VERBOSE 1 first the first-grid ones !!
!! (with CORRELATION, exchange and correlation separately: only where   !!
!! no zero-error correction applies, which corrects their sum)          !!
      if(izerr.eq.1.and.iverbose.eq.1) then
        call enpart_phf_mats(coulb,xcoulb,excb,xexcb,Ex,exenx,Ec,exenc,
     +    iecorr.eq.1,', BEFORE THE ZERO-ERROR CORRECTION')
        DEALLOCATE(coulb,excb)
      end if
      call enpart_phf_mats(coul,coulen,Exc,exch_corr,Ex,exenx,Ec,exenc,
     +  iecorr.eq.1.and.izerr.eq.0,' ')
      if(iecorr.eq.1.and.izerr.eq.1) then
        write(*,'(a)') ''
        write(*,'(2x,a)') 'CORRELATION: the zero-error correction applies'//
     +    ' to the exchange-correlation'
        write(*,'(2x,a)') 'as a whole, so exchange and correlation are'//
     +    ' given separately only before'
        write(*,'(2x,a)') 'it (VERBOSE 1).'
      end if

!! IQA matrix: one-electron terms (in eto) + Coulomb + xc !!
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

      call enpart_iqa(eto,etot,escf,0,ZERO,Exc,.true.)

      DEALLOCATE(chp2,chp3)

      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: monadic_diag
!! Builds the cumulant of the 2-RDM (exchange-correlation, exchange or
!! correlation), compacts and diagonalizes it, and evaluates its pair
!! functions (eigenvectors on the orbital products) on both grids.
!! arguments:
!!   ic      (in)  0 exchange-correlation, 1 exchange, 2 correlation
!!   itotps  (in)  number of grid points
!!   chp2    (in)  MOs on the first grid
!!   chp2pha (in)  MOs on the rotated grid
!!   dm1,dm2 (in)  the RDMs
!!   cumulab (out) diagonalized compacted cumulant (eigenvalues on the
!!                 diagonal)
!!   fij     (out) pair functions on the first grid
!!   fijb    (out) pair functions on the rotated grid
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine monadic_diag(ic,itotps,chp2,chp2pha,dm1,dm2,cumulab,fij,fijb)

      IMPLICIT REAL*8(A-H,O-Z)

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/ icas,ncasel,ncasorb,nspinorb,norb,icisd,icass

      dimension :: chp2(itotps,norb),chp2pha(itotps,norb)
      dimension :: dm1(nspinorb,nspinorb),dm2(norb,norb,norb,norb)
      dimension :: cumulab(norb*(norb+1)/2,norb*(norb+1)/2)
      dimension :: fij(itotps,norb*(norb+1)/2),fijb(itotps,norb*(norb+1)/2)

      allocatable :: cumul(:,:,:,:),eigenmat(:,:),fa1(:,:),fb2(:,:)
      allocatable :: sfdm1(:,:),psdm1(:,:) 

!! cumulant of the 2-RDM: exchange-correlation (ic=0), exchange (ic=1) !!
!! or correlation (ic=2)                                                  !!
      ALLOCATE(cumul(norb,norb,norb,norb))
      ALLOCATE(sfdm1(norb,norb),psdm1(norb,norb))

!! total and spin 1-RDMs in the MO basis !!
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

!! fold the cumulant onto the pairs i<=j and k<=l, then compact them !!
!! into one index each (ij, kl)                                      !!
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

!! diagonalize the compacted cumulant !!
      norb2=norb*(norb+1)/2
      ALLOCATE(eigenmat(norb2,norb2))
      call diagonalize(norb2,norb2,cumulab,eigenmat,0)

!! orbital products on the first grid (fa1) and the rotated one (fb2) !!
      ALLOCATE(fa1(itotps,norb2),fb2(itotps,norb2))
!! parallel over grid points: each ifut writes only fa1(ifut,:), fb2(ifut,:) !!
!$OMP PARALLEL DO PRIVATE(ifut,ij,ii,jj)
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
!$OMP END PARALLEL DO

!! pair functions: the products transformed by the eigenvectors !!
!! parallel over grid points: each ifut writes only fij(ifut,:), fijb(ifut,:) !!
!$OMP PARALLEL DO PRIVATE(ifut,ii,xx,xxb,kk)
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
!$OMP END PARALLEL DO
      DEALLOCATE(fa1,fb2,eigenmat)

      end 

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: numint_two_pmon
!! Integrates the cumulant pair functions over two electrons: one-center
!! terms with the rotated grid for the second electron, pair terms for
!! the atom pairs at or above THREBOD (the others are left to the
!! multipolar expansion). Adds them to ef.
!! arguments:
!!   itotps        (in)  number of grid points
!!   iorb2         (in)  number of pair functions
!!   wp,wppha      (in)  weights on the first and the rotated grid
!!   omp2,omp2pha  (in)  atomic weights on both grids
!!   pcoord,pcoordpha (in) coordinates on both grids
!!   cumulab       (in)  cumulant eigenvalues (diagonal)
!!   ffa1,ffb2     (in)  pair functions on the first and the rotated grid
!!   ef            (inout) exchange(-correlation) matrix (one-center terms
!!                 counted twice, pair terms once per pair)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine numint_two_pmon(itotps,iorb2,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,cumulab,ffa1,ffb2,ef)

      use integration_grid

      IMPLICIT REAL*8(A-H,O-Z)

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(200)

      dimension :: wp(itotps),wppha(itotps),pcoord(itotps,3),pcoordpha(itotps,3)
      dimension :: omp2(itotps,nat),omp2pha(itotps,nat),cumulab(iorb2,iorb2)
      dimension :: ffa1(itotps,iorb2),ffb2(itotps,iorb2)
      dimension :: ef(maxat,maxat),efrow(maxat)

      ithrebod = iopt(44)
      iatps    = nrad*nang

      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if
      

!! two-electron integration of the pair functions, weighted by the !!
!! cumulant eigenvalues                                           !!

!! parallel over the points ifut of each atom: f3 and the point loop  !!
!! are PRIVATE; the row efrow(icenter:nat) of the atom is a reduction,  !!
!! added to ef after the loop (on one thread the summation order is the !!
!! serial one)                                                          !!
      do icenter=1,nat
        do jcenter=1,nat
          efrow(jcenter)=ZERO
        end do
!$OMP PARALLEL DO PRIVATE(ifut,x0,dx0,dy0,dz0,f3,jfut,x1,dx1,dy1,dz1,
!$OMP&  dist,xxw,ij,jcenter,bx0) REDUCTION(+:efrow)
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
          efrow(icenter)=efrow(icenter)+f3

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
                if(dist.gt.1.0d-12) then
                  xxw=x0*x1/dist
                  do ij=1,iorb2
                    f3=f3+cumulab(ij,ij)*xxw*ffa1(ifut,ij)*ffa1(jfut,ij)
                  end do
                end if
              end do
              efrow(jcenter)=efrow(jcenter)+f3
            end if
          end do
        end do
!$OMP END PARALLEL DO
        do jcenter=icenter,nat
          ef(icenter,jcenter)=ef(icenter,jcenter)+efrow(jcenter)
        end do
      end do

      end 

!! ***** !!

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: enpart_phf_mats
!! Prints the two-electron matrices of a correlated wavefunction:
!! Coulomb, exchange-correlation and, with CORRELATION, exchange and
!! correlation separately.
!! arguments:
!!   coul,xcoul (in)  Coulomb matrix and total
!!   exc,xexc   (in)  exchange-correlation matrix and total
!!   ex,xex     (in)  exchange matrix and total (lxc)
!!   ec,xec     (in)  correlation matrix and total (lxc)
!!   lxc        (in)  logical, print exchange and correlation separately
!!   cstage     (in)  appended to the titles (' ' for the final values)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine enpart_phf_mats(coul,xcoul,exc,xexc,ex,xex,ec,xec,lxc,
     +  cstage)
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      dimension coul(maxat,maxat),exc(maxat,maxat),ex(maxat,maxat)
      dimension ec(maxat,maxat)
      logical lxc
      character*(*) cstage
      character*80 ct
      logical lfin

!! the results file gets the final matrices only !!
      lfin=len_trim(cstage).eq.0
      ct='COULOMB ENERGY'//cstage
      call enpart_matrix(ct,' ',' ',coul,'Coulomb energy',xcoul)
      call enpart_frag(ct,coul)
      if(lfin) call json_enpart_term('coulomb','Coulomb energy',coul,
     +  xcoul,.false.,0.0d0)
      ct='EXCHANGE-CORRELATION ENERGY'//cstage
      call enpart_matrix(ct,'From the cumulant of the 2-RDM.',' ',exc,
     +  'Exchange-correlation energy',xexc)
      call enpart_frag(ct,exc)
      if(lfin) call json_enpart_term('xc','Exchange-correlation energy',
     +  exc,xexc,.false.,0.0d0)
      if(lxc) then
        ct='EXCHANGE ENERGY'//cstage
        call enpart_matrix(ct,'The exchange part of the exchange-'//
     +    'correlation.',' ',ex,'Exchange energy',xex)
        call enpart_frag(ct,ex)
        if(lfin) call json_enpart_term('exchange','Exchange energy',ex,
     +    xex,.false.,0.0d0)
        ct='CORRELATION ENERGY'//cstage
        call enpart_matrix(ct,'Exchange-correlation minus exchange.',
     +    ' ',ec,'Correlation energy',xec)
        call enpart_frag(ct,ec)
        if(lfin) call json_enpart_term('correlation','Correlation '//
     +    'energy',ec,xec,.false.,0.0d0)
      end if
      end
