!! ************************************************************************** !!
!! KS-DFT EXCHANGE-CORRELATION ENERGY PARTITIONING SUBROUTINES                !!
!! RKS (closed-shell):                                                        !!
!!   numint_dft      -- BODEN approximation + "exact" one-center XC term      !!
!!   xc              -- XC energy density at each grid point, via libxc       !!
!!   sigma           -- density-gradient contraction, GGA/hybrid-GGA only     !!
!!   grdrho          -- MO gradient at each grid point, GGA/hybrid-GGA only   !!
!!   grdboden        -- BODEN gradient for one atom pair, GGA/hybrid-GGA only !!
!! UKS (open-shell):                                                          !!
!!   numint_dft_uks  -- UKS twin of numint_dft                                !!
!!   xc_uks          -- UKS twin of xc                                        !!
!!   sigma_uks       -- UKS twin of sigma                                     !!
!!   ugrdrho         -- UKS twin of grdrho                                    !!
!!   grdboden_uks    -- UKS twin of grdboden                                  !!
!! shared, authored PSalse, MGimf.:                                           !!
!!   func_info_print -- functional identification/classification, called once !!
!!                      by main.f before dispatching to either half above     !!
!! ************************************************************************** !!

!! ***** !!

!! *********************************************************************** !!
!! subroutine: xc                                                          !!
!! purpose: restricted XC energy density at each grid point, via libxc.    !!
!!   Exchange and correlation are queried separately when the input        !!
!!   functional specifies them as distinct libxc ids (id_xfunc/id_cfunc),  !!
!!   otherwise as one combined exchange-correlation id (id_xcfunc). MGGA   !!
!!   Laplacian/kinetic-energy-density terms not implemented. See xc_uks    !!
!!   for open-shell.                                                       !!
!! arguments:                                                              !!
!!   npt    (in)  -- number of grid points                                 !!
!!   scr_a  (in)  -- electron density at each point                        !!
!!   scr    (in)  -- density gradient contraction (sigma) at each point,   !!
!!                   only read for GGA/hybrid-GGA functionals              !!
!!   scr2   (out) -- XC energy density at each point (already multiplied   !!
!!                   by the density)                                       !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine xc(npt,scr_a,scr,scr2)
      use xc_f90_types_m
      use xc_f90_lib_m
      implicit real*8(a-h,o-z)
      TYPE(xc_f90_pointer_t) :: xc_func
      TYPE(xc_f90_pointer_t) :: xc_info
      include 'parameter.h'
      integer npt
      real*8 scr_a(npt),scr(npt),scr2(npt)
      common /iops/iopt(200)
      allocatable :: scr2c(:)

      id_xcfunc = iopt(60)
      id_cfunc  = iopt(62)
      id_xfunc  = iopt(61)

      ALLOCATE(scr2c(npt))
      do ii=1,npt
        scr2c(ii)=ZERO
      end do

!! exchange-correlation, from libxc. !!
      if(id_xcfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_xcfunc,XC_UNPOLARIZED)
        select case (xc_f90_info_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f90_lda_exc(xc_func,npt,scr_a(1),scr2(1))
        case(XC_FAMILY_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2(1))
        case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2(1))
        case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! correlation, when specified as a separate libxc id. !!
      if(id_cfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_cfunc,XC_UNPOLARIZED)
        select case (xc_f90_info_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f90_lda_exc(xc_func,npt,scr_a(1),scr2c(1))
        case(XC_FAMILY_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2c(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2c(1))
        case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2c(1))
        case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2c(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! exchange, when specified as a separate libxc id. !!
      if(id_xfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_xfunc,XC_UNPOLARIZED)
        select case (xc_f90_info_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f90_lda_exc(xc_func,npt,scr_a(1),scr2(1))
        case(XC_FAMILY_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2(1))
        case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2(1))
        case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! combine exchange and correlation, multiply by the density. !!
      do ii=1,npt
        scr2(ii)=(scr2c(ii)+scr2(ii))*scr_a(ii)
      end do
      DEALLOCATE(scr2c)

      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: numint_dft                                                  !!
!! purpose: KS-DFT exchange-correlation energy partition, RKS -- BODEN     !!
!!   (bond-order density) approximation for atom-pair XC contributions,    !!
!!   plus an "exact" one-center XC term via direct integration, combined   !!
!!   into the final atomic/diatomic XC decomposition. See numint_dft_uks   !!
!!   for open-shell.                                                       !!
!! arguments:                                                              !!
!!   ndim   (in)    -- number of basis functions (leading dim of chp)      !!
!!   itotps (in)    -- total number of grid points                        !!
!!   wp     (in)    -- integration weight of each grid point               !!
!!   rho    (in)    -- electron density at each grid point                 !!
!!   omp    (in)    -- becke/tfvc weight of each grid point for its own atom !!
!!   omp2   (in)    -- becke/tfvc (or hirshfeld) weight of each point for  !!
!!                      every atom                                         !!
!!   chp    (in)    -- basis-function values at each grid point            !!
!!   eto    (inout) -- total energy matrix, accumulated on top of the      !!
!!                      one- and two-electron parts already in it          !!
!!   pcoord (in)    -- xyz coordinates of each grid point                  !!
!!   sat    (in)    -- per-atom AO overlap matrix (from numint_sat)        !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine numint_dft(ndim,itotps,wp,rho,omp,omp2,chp,eto,pcoord,sat)
      use xc_f90_types_m
      use xc_f90_lib_m
      use ao_matrices
      use integration_grid
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(200)
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/exchg/exch(maxat,maxat),xmix
      dimension eto(maxat,maxat)
      dimension wp(itotps),chp(itotps,ndim)
      dimension omp(itotps),omp2(itotps,nat),rho(itotps),pcoord(itotps,3)
      dimension exch2(maxat,maxat)
      dimension sat(igr,igr,nat)
      character*80 line

      allocatable :: chp2(:,:),scr(:),scrx(:,:)
      allocatable :: scr_a(:),rho_a(:),scr2(:)
      allocatable :: sab(:,:,:),chpd(:,:,:)

      idofr    =  Iopt(40)
      ithrebod =  Iopt(44)
      itype    =  Iopt(55)
      iatps    =  nang*nrad

      ALLOCATE(chp2(itotps,nocc),scr(itotps))
      ALLOCATE(scr_a(itotps),rho_a(itotps),scr2(itotps))
      ALLOCATE(sab(nocc,nocc,nat),chpd(itotps,nocc,3))
      ALLOCATE(scrx(igr,nocc))

      exch2=ZERO

!! threshold for BODEN calculation, from the bond order between a pair of atoms. !!
      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if

      call print_box('RKS-DFT ENERGY DECOMPOSITION')
      write(*,*) " USING BOND ORDER DENSITY APPROACH "
      write(*,'(2x,a36,1x,f10.6)') "Threshold for atom pair calculation:",threbod

!! per-atom AO overlap projected onto MOs (sab), used by the BODEN kernel !!
!! below. parallel over iatom: each iteration only writes its own         !!
!! scrx(:,:)/sab(:,:,iatom), independent across atoms.                    !!
!$OMP PARALLEL DO PRIVATE(iatom,nu,ii,mu,xx,jj,scrx)
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
!$OMP END PARALLEL DO

!! sanity check: sab summed over atoms must reproduce the MO overlap     !!
!! (identity for ii=jj, zero otherwise). !!
      do ii=1,nocc
        do jj=ii,nocc
          x=ZERO
          if(ii.eq.jj) x=-ONE
          do kk=1,nat
            x=x+sab(ii,jj,kk)
          end do
          if(abs(x).gt.1.0d-3) then
            write(*,*) ii,jj,xx
            stop " PROBLEM WITH MOs OVERLAPS "
          end if
        end do
      end do

      call ao_to_mo_grid(itotps,igr,nocc,c,chp,chp2)

      if(itype.gt.1) call grdrho(pcoord,chpd)

      call print_box('BOND ORDER DENSITY FOR ALL ATOM PAIRS')
      write(*,*) " --------------------------- "
      write(*,*) "  Atom   Atom   BODEN value  "
      write(*,*) " --------------------------- "
      do iatom=1,nat
        do jatom=iatom+1,nat

!! BODEN for this atom pair, skipped below the THREBOD bond-order threshold. !!
          bx0=bo(iatom,jatom)
          if(bx0.ge.threbod) then
            x1=ZERO
!! ff2 is PRIVATE, scr_a(jfut) written at a unique index per iteration --  !!
!! no false-sharing risk. x1 is diagnostic-only, REDUCTION is safe. jfut  !!
!! is computed in the body since COLLAPSE needs loop-invariant bounds.    !!
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,jloc,jfut,x2,wa,wb,ff2,i,j,ff) REDUCTION(+:x1)
            do icenter=1,nat
              do jloc=1,iatps
                jfut=iatps*(icenter-1)+jloc
                x2=wp(jfut)*omp(jfut)
                wa=omp2(jfut,iatom)
                wb=omp2(jfut,jatom)
                ff2=ZERO
                do i=1,nocc
                  do j=1,nocc
                    ff=sab(i,j,jatom)*wa+sab(j,i,iatom)*wb
                    ff2=ff2+ff*chp2(jfut,i)*chp2(jfut,j)
                  end do
                end do

!! DEFINING A-B BODEN AS A-B + B-A, HENCE FACTOR OF 2
                scr_a(jfut)=ff2
                if(iatom.ne.jatom) scr_a(jfut)=TWO*ff2
                x1=x1+x2*scr_a(jfut)
              end do
            end do
!$OMP END PARALLEL DO
            write(*,'(4x,i3,4x,i3,4x,f10.7)') iatom,jatom,x1

!! BODEN gradient, for GGA functionals -- scr_a: BODEN for this atom pair, !!
!! scr: sigma BODEN, scr2: XC functional value. MGGA Laplacian term not   !!
!! implemented. !!
            if(itype.gt.1) call grdboden(itotps,omp2,chp2,chpd,scr,iatom,jatom,sab)
            call xc(itotps,scr_a,scr,scr2)

            x1=ZERO
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,iloc,ifut) REDUCTION(+:x1)
            do icenter=1,nat
              do iloc=1,iatps
                ifut=iatps*(icenter-1)+iloc
                x1=x1+wp(ifut)*scr2(ifut)*omp2(ifut,icenter)*omp(ifut)
              end do
            end do
!$OMP END PARALLEL DO
            exch2(iatom,jatom)=x1
          end if
        end do
      end do
      write(*,*) " --------------------------- "

      call print_box('DIATOMIC PURE KS-DFT XC TERMS (BODEN)')
      exchen=ZERO
      do ii=1,nat
        do jj=ii,nat
          if(ii.ne.jj) exch2(jj,ii)=exch2(ii,jj)
          exchen=exchen+exch2(ii,jj)
        end do
      end do
      call MPRINT2(exch2,nat,maxat)

      if(itype.gt.1) call sigma(pcoord,chp2,scr)

!! "exact" one-center XC energy, by direct integration (no BODEN approximation). !!
!! MGGA Laplacian term not implemented (would need qtaim.f's laplacian routine). !!
      call xc(itotps,rho,scr,scr2)
      xtot=ZERO
!! parallel over icenter: each iteration writes only its own exch(icenter,icenter), !!
!! xtot is a genuine running total. !!
!$OMP PARALLEL DO PRIVATE(icenter,x,ifut) REDUCTION(+:xtot)
      do icenter=1,nat
        x=ZERO
        do ifut=1,itotps
          x=x+wp(ifut)*scr2(ifut)*omp(ifut)*omp2(ifut,icenter)
        end do
        exch(icenter,icenter)=x
        xtot=xtot+x
      end do
!$OMP END PARALLEL DO

      call print_box('PURE KS-DFT XC ONE-CENTER TERMS (EXACT)')
      call MPRINT2(exch,nat,maxat)
      write(*,'(2x,a35,x,f14.7)') "KS-DFT exchange-correlation energy:",xtot
      write(*,*) " "

      write(*,*) " REARRANGING ATOMIC COMPONENTS "
      do ii=1,nat
        x0=ZERO
        do jj=1,nat
          if(ii.ne.jj) then
            exch(ii,jj)=exch2(ii,jj)
            x0=x0+exch(ii,jj)
          end if
        end do
        exch(ii,ii)=exch(ii,ii)-x0/TWO
      end do

      call print_box('FINAL PURE KS-DFT EXCHANGE-CORRELATION ENERGY COMPONENTS')
      exchen=ZERO
      do ii=1,nat
        do jj=ii,nat
          exchen=exchen+exch(ii,jj)
        end do
      end do
      call MPRINT2(exch,nat,maxat)
      write(*,'(2x,a47,x,f14.7)') "Sum of pure KS-DFT exchange-correlation energy:",exchen
      if(xmix.gt.ZERO) then
        write(*,*) " "
        write(*,*) " WARNING: HF-exchange part missing "
      end if
      if (idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Exc Decomposition'
        call group_by_frag_mat(1,line,exch)
      end if
      write(*,*) " "

      xtot=ZERO
      do ii=1,nat
        eto(ii,ii)=eto(ii,ii)+exch(ii,ii)
        xtot=xtot+eto(ii,ii)
        do jj=ii+1,nat
          eto(ii,jj)=eto(ii,jj)+exch(ii,jj)
          eto(jj,ii)=eto(ii,jj)
          xtot=xtot+eto(ii,jj)
        end do
      end do
      etot=xtot
      DEALLOCATE(chp2,scr,scr_a,rho_a,scr2,sab,chpd,scrx)
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: sigma                                                       !!
!! purpose: density-gradient contraction (|grad(rho)|^2, "sigma" in        !!
!!   libxc's convention) at each grid point, for GGA/hybrid-GGA XC         !!
!!   functionals. Builds the AO gradient (x/y/z in turn) from primitives,  !!
!!   transforms to MOs, and contracts with the already-transformed chp2.   !!
!!   See sigma_uks for open-shell.                                         !!
!! arguments:                                                              !!
!!   pcoord (in)  -- xyz coordinates of each grid point                    !!
!!   chp2   (in)  -- MO values at each grid point (from the caller)        !!
!!   scr    (out) -- sigma at each grid point                              !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine sigma(pcoord,chp2,scr)
      use basis_set
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chp2(natoms*nrad*nang,nocc),scr(natoms*nang*nrad)
      dimension pcoord(natoms*nang*nrad,3)
      allocatable:: chpd(:,:),chp(:,:)

      ipoints=natoms*nrad*nang
      ALLOCATE(chpd(ipoints,igr),chp(ipoints,igr))

      do kk=1,ipoints
        scr(kk)=ZERO
      end do

      do ixyz=1,3

!! AO gradient (ixyz component) from primitives. parallel over irun: each !!
!! iteration writes only its own chp(irun,:), independent across points.  !!
!$OMP PARALLEL DO PRIVATE(irun,iact,iactat,x,y,z,rr,f,k,ipr,nn,ll,mm,alpha,dx)
        do irun=1,ipoints
          do iact=1,nbasis
            iactat=ihold(iact)
            x=pcoord(irun,1)-coord(1,iactat)
            y=pcoord(irun,2)-coord(2,iactat)
            z=pcoord(irun,3)-coord(3,iactat)
            rr=dsqrt(x**2.0d0+y**2.0d0+z**2.0d0)
            f=0.d0
            k=1
            do while(nprimbas(k,iact).ne.0)
              ipr=nprimbas(k,iact)
              nn=nlm(ipr,1)
              ll=nlm(ipr,2)
              mm=nlm(ipr,3)
              alpha=expp(ipr)
              if(ixyz.eq.1) then
                dx=-2.0d0*alpha*(x**(nn+1))
                if(nn.ge.1) dx=dx+nn*x**(nn-1)
                dx=dx*(y**ll)*(z**mm)*dexp(-expp(ipr)*(rr**2))
              else if(ixyz.eq.2) then
                dx=-2.0d0*alpha*(y**(ll+1))
                if(ll.ge.1) dx=dx+ll*y**(ll-1)
                dx=dx*(x**nn)*(z**mm)*dexp(-expp(ipr)*(rr**2))
              else
                dx=-2.0d0*alpha*(z**(mm+1))
                if(mm.ge.1) dx=dx+mm*z**(mm-1)
                dx=dx*(x**nn)*(y**ll)*dexp(-expp(ipr)*(rr**2))
              end if
              f=f+dx*coefpb(ipr,iact)
              k=k+1
            enddo
            chp(irun,iact)=f
          enddo
        enddo
!$OMP END PARALLEL DO

        call ao_to_mo_grid(ipoints,igr,nocc,c,chp,chpd)

!! contract with chp2 (already MO-transformed by the caller). parallel   !!
!! over k: each iteration only reads chp2(k,:)/chpd(k,:) and writes its  !!
!! own scr(k) -- safe across the three ixyz passes since they run        !!
!! serially, only the point loop within each pass is threaded.          !!
!$OMP PARALLEL DO PRIVATE(k,i,j)
        do k=1,ipoints
          do i=1,nocc
            do j=1,nocc
              scr(k)=scr(k)+chp2(k,i)*chp2(k,j)*chpd(k,i)*chpd(k,j)
            end do
          end do
        end do
!$OMP END PARALLEL DO
      end do

      do kk=1,ipoints
        scr(kk)=16.0d0*scr(kk)
      end do
      deallocate(chpd,chp)
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: grdrho                                                      !!
!! purpose: MO gradient (x/y/z in turn) at each grid point, for GGA/       !!
!!   hybrid-GGA XC functionals. Builds the AO gradient from primitives,    !!
!!   then transforms to MOs. See ugrdrho for open-shell.                   !!
!! arguments:                                                              !!
!!   pcoord (in)  -- xyz coordinates of each grid point                    !!
!!   chpd   (out) -- MO gradient at each grid point, one slice per         !!
!!                   Cartesian component                                   !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine grdrho(pcoord,chpd)
      use basis_set
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chpd(natoms*nrad*nang,nocc,3),pcoord(natoms*nrad*nang,3)
      allocatable:: chp(:,:)

      itotps=natoms*nrad*nang
      ALLOCATE(chp(itotps,igr))

      do ixyz=1,3

!! AO gradient (ixyz component) from primitives. parallel over irun: each !!
!! iteration writes only its own chp(irun,:), independent across points.  !!
!$OMP PARALLEL DO PRIVATE(irun,iact,iactat,x,y,z,rr,f,k,ipr,nn,ll,mm,alpha,dx)
        do irun=1,itotps
          do iact=1,nbasis
            iactat=ihold(iact)
            x=pcoord(irun,1)-coord(1,iactat)
            y=pcoord(irun,2)-coord(2,iactat)
            z=pcoord(irun,3)-coord(3,iactat)
            rr=dsqrt(x**2.0d0+y**2.0d0+z**2.0d0)
            f=0.d0
            k=1
            do while(nprimbas(k,iact).ne.0)
              ipr=nprimbas(k,iact)
              nn=nlm(ipr,1)
              ll=nlm(ipr,2)
              mm=nlm(ipr,3)
              alpha=expp(ipr)
              if(ixyz.eq.1) then
                dx=-2.0d0*alpha*(x**(nn+1))
                if(nn.ge.1) dx=dx+nn*x**(nn-1)
                dx=dx*(y**ll)*(z**mm)*dexp(-expp(ipr)*(rr**2))
              else if(ixyz.eq.2) then
                dx=-2.0d0*alpha*(y**(ll+1))
                if(ll.ge.1) dx=dx+ll*y**(ll-1)
                dx=dx*(x**nn)*(z**mm)*dexp(-expp(ipr)*(rr**2))
              else
                dx=-2.0d0*alpha*(z**(mm+1))
                if(mm.ge.1) dx=dx+mm*z**(mm-1)
                dx=dx*(x**nn)*(y**ll)*dexp(-expp(ipr)*(rr**2))
              end if
              f=f+dx*coefpb(ipr,iact)
              k=k+1
            enddo
            chp(irun,iact)=f
          enddo
        enddo
!$OMP END PARALLEL DO

        call ao_to_mo_grid(itotps,igr,nocc,c,chp,chpd(1,1,ixyz))

      end do
      DEALLOCATE(chp)
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: ugrdrho                                                     !!
!! purpose: MO gradient (x/y/z in turn), alpha and beta, at each grid      !!
!!   point, for GGA/hybrid-GGA XC functionals. Builds the AO gradient      !!
!!   from primitives, then transforms to MOs. See grdrho for               !!
!!   restricted/closed-shell.                                              !!
!! arguments:                                                              !!
!!   pcoord (in)  -- xyz coordinates of each grid point                    !!
!!   chpd   (out) -- alpha MO gradient at each grid point, one slice per   !!
!!                   Cartesian component                                   !!
!!   chpbd  (out) -- beta MO gradient at each grid point, one slice per    !!
!!                   Cartesian component                                   !!
!! author: MGimf.                                                          !!
!! *********************************************************************** !!
      subroutine ugrdrho(pcoord,chpd,chpbd)
      use basis_set
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chpd(natoms*nrad*nang,nalf,3),pcoord(natoms*nrad*nang,3)
      dimension chpbd(natoms*nrad*nang,nb,3)
      allocatable:: chp(:,:)

      itotps=natoms*nrad*nang
      ALLOCATE(chp(itotps,igr))

      do ixyz=1,3

!! AO gradient (ixyz component) from primitives. parallel over irun: each !!
!! iteration writes only its own chp(irun,:), independent across points.  !!
!$OMP PARALLEL DO PRIVATE(irun,iact,iactat,x,y,z,rr,f,k,ipr,nn,ll,mm,alpha,dx)
        do irun=1,itotps
          do iact=1,nbasis
            iactat=ihold(iact)
            x=pcoord(irun,1)-coord(1,iactat)
            y=pcoord(irun,2)-coord(2,iactat)
            z=pcoord(irun,3)-coord(3,iactat)
            rr=dsqrt(x**2.0d0+y**2.0d0+z**2.0d0)
            f=0.d0
            k=1
            do while(nprimbas(k,iact).ne.0)
              ipr=nprimbas(k,iact)
              nn=nlm(ipr,1)
              ll=nlm(ipr,2)
              mm=nlm(ipr,3)
              alpha=expp(ipr)
              if(ixyz.eq.1) then
                dx=-2.0d0*alpha*(x**(nn+1))
                if(nn.ge.1) dx=dx+nn*x**(nn-1)
                dx=dx*(y**ll)*(z**mm)*dexp(-expp(ipr)*(rr**2))
              else if(ixyz.eq.2) then
                dx=-2.0d0*alpha*(y**(ll+1))
                if(ll.ge.1) dx=dx+ll*y**(ll-1)
                dx=dx*(x**nn)*(z**mm)*dexp(-expp(ipr)*(rr**2))
              else
                dx=-2.0d0*alpha*(z**(mm+1))
                if(mm.ge.1) dx=dx+mm*z**(mm-1)
                dx=dx*(x**nn)*(y**ll)*dexp(-expp(ipr)*(rr**2))
              end if
              f=f+dx*coefpb(ipr,iact)
              k=k+1
            enddo
            chp(irun,iact)=f
          enddo
        enddo
!$OMP END PARALLEL DO

        call ao_to_mo_grid(itotps,igr,nalf,c,chp,chpd(1,1,ixyz))
        call ao_to_mo_grid(itotps,igr,nb,cb,chp,chpbd(1,1,ixyz))

      end do
      DEALLOCATE(chp)
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: grdboden                                                    !!
!! purpose: gradient of the bond-order density (BODEN) for one atom pair,  !!
!!   at each grid point -- sigma input for the GGA/hybrid-GGA XC           !!
!!   evaluation in numint_dft's BODEN loop. See grdboden_uks for           !!
!!   open-shell.                                                           !!
!! arguments:                                                              !!
!!   itotps (in)  -- total number of grid points                          !!
!!   omp2   (in)  -- becke/tfvc (or hirshfeld) weight of each point for    !!
!!                   every atom                                            !!
!!   chp2   (in)  -- MO values at each grid point                          !!
!!   chpd   (in)  -- MO gradient at each grid point (from grdrho)          !!
!!   scr    (out) -- BODEN gradient at each grid point                     !!
!!   iat,kat (in) -- the atom pair                                         !!
!!   sab    (in)  -- per-atom MO overlap matrix (from numint_dft)          !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine grdboden(itotps,omp2,chp2,chpd,scr,iat,kat,sab)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chp2(itotps,nocc),scr(itotps)
      dimension chpd(itotps,nocc,3),omp2(itotps,nat)
      dimension sab(nocc,nocc,nat)

      do kk=1,itotps
        scr(kk)=ZERO
      end do
      do ixyz=1,3

!! parallel over irun: each iteration only reads its own omp2(irun,:)/   !!
!! chpd(irun,:,ixyz)/chp2(irun,:) and writes its own scr(irun) -- safe   !!
!! across the three ixyz passes since they run serially, only the point  !!
!! loop within each pass is threaded.                                    !!
!$OMP PARALLEL DO PRIVATE(irun,w1,w2,xx,ii,jj,gab)
        do irun=1,itotps
          w1=omp2(irun,iat)
          w2=omp2(irun,kat)
          xx=ZERO
          do ii=1,nocc
            do jj=1,nocc
              gab=sab(ii,jj,iat)*w2+sab(ii,jj,kat)*w1
              xx=xx+gab*chpd(irun,jj,ixyz)*chp2(irun,ii)
            end do
          end do
          scr(irun)=scr(irun)+xx*xx
        end do
!$OMP END PARALLEL DO
      end do

      xfact=FOUR
      if(iat.ne.kat) xfact=xfact*FOUR
      do kk=1,itotps
        scr(kk)=xfact*scr(kk)
      end do
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: xc_uks                                                      !!
!! purpose: unrestricted XC energy density at each grid point, via libxc.  !!
!!   Exchange and correlation are queried separately when the input        !!
!!   functional specifies them as distinct libxc ids (id_xfunc/id_cfunc),  !!
!!   otherwise as one combined exchange-correlation id (id_xcfunc). MGGA   !!
!!   Laplacian/kinetic-energy-density terms not implemented. See xc for    !!
!!   restricted/closed-shell.                                              !!
!! arguments:                                                              !!
!!   npt    (in)  -- number of grid points                                 !!
!!   scr_ab (in)  -- electron density at each point, order alpha then beta !!
!!   scr    (in)  -- density gradient contraction (sigma) at each point,   !!
!!                   order up-up/up-down/down-down, only read for GGA/     !!
!!                   hybrid-GGA functionals                                 !!
!!   scr2   (out) -- XC energy density at each point (already multiplied   !!
!!                   by the density)                                       !!
!! author: MGimf.                                                          !!
!! *********************************************************************** !!
      subroutine xc_uks(npt,scr_ab,scr,scr2)
      use xc_f90_types_m
      use xc_f90_lib_m
      implicit real*8(a-h,o-z)
      TYPE(xc_f90_pointer_t) :: xc_func
      TYPE(xc_f90_pointer_t) :: xc_info
      include 'parameter.h'
      integer npt
      common /iops/iopt(200)
      dimension scr_ab(2,npt),scr(3,npt),scr2(npt)
      allocatable :: scr2c(:)

      id_xcfunc = iopt(60)
      id_cfunc  = iopt(62)
      id_xfunc  = iopt(61)

      ALLOCATE(scr2c(npt))
      do ii=1,npt
        scr2c(ii)=ZERO
      end do

!! exchange-correlation, from libxc. !!
      if(id_xcfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_xcfunc,XC_POLARIZED)
        select case (xc_f90_info_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f90_lda_exc(xc_func,npt,scr_ab(1,1),scr2(1))
        case(XC_FAMILY_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2(1))
        case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2(1))
        case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! correlation, when specified as a separate libxc id. !!
      if(id_cfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_cfunc,XC_POLARIZED)
        select case (xc_f90_info_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f90_lda_exc(xc_func,npt,scr_ab(1,1),scr2c(1))
        case(XC_FAMILY_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2c(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2c(1))
        case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2c(1))
        case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2c(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! exchange, when specified as a separate libxc id. !!
      if(id_xfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_xfunc,XC_POLARIZED)
        select case (xc_f90_info_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f90_lda_exc(xc_func,npt,scr_ab(1,1),scr2(1))
        case(XC_FAMILY_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2(1))
        case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2(1))
        case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! combine exchange and correlation, multiply by the density. !!
      do ii=1,npt
        scr2(ii)=(scr2c(ii)+scr2(ii))*(scr_ab(1,ii)+scr_ab(2,ii))
      end do
      DEALLOCATE(scr2c)
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: numint_dft_uks                                              !!
!! purpose: KS-DFT exchange-correlation energy partition, UKS twin of      !!
!!   numint_dft -- BODEN (bond-order density) approximation for atom-pair  !!
!!   XC contributions (alpha+beta), plus an "exact" one-center XC term     !!
!!   via direct integration, combined into the final atomic/diatomic XC    !!
!!   decomposition. See numint_dft for restricted/closed-shell.            !!
!! arguments:                                                              !!
!!   ndim   (in)    -- number of basis functions (leading dim of chp)      !!
!!   itotps (in)    -- total number of grid points                        !!
!!   wp     (in)    -- integration weight of each grid point               !!
!!   omp    (in)    -- becke/tfvc weight of each grid point for its own atom !!
!!   omp2   (in)    -- becke/tfvc (or hirshfeld) weight of each point for  !!
!!                      every atom                                         !!
!!   chp    (in)    -- basis-function values at each grid point            !!
!!   pcoord (in)    -- xyz coordinates of each grid point                  !!
!!   eto    (inout) -- total energy matrix, accumulated on top of the      !!
!!                      one- and two-electron parts already in it          !!
!! author: MGimf.                                                          !!
!! *********************************************************************** !!
      subroutine numint_dft_uks(ndim,itotps,wp,omp,omp2,chp,pcoord,eto)
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(200)
      common/actual/iact,jat,icenter
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/exchg/exch(maxat,maxat),xmix
      character*80 line
      dimension eto(maxat,maxat)
      dimension wp(itotps),chp(itotps,ndim)
      dimension omp(itotps),omp2(itotps,nat),pcoord(itotps,3)
      dimension exch2(maxat,maxat)
      allocatable :: chp2(:,:),chp3(:,:),rho(:,:),scrall(:,:),scr2(:)
      allocatable :: chpd(:,:,:),sab(:,:,:),sab2(:,:,:),scr_bod(:,:)
      allocatable :: chpbd(:,:,:)

      idofr    = Iopt(40)
      ithrebod = Iopt(44)
      itype    = Iopt(55)
      iatps    = nang*nrad

      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if

      ALLOCATE(chp2(itotps,nalf),rho(2,itotps),scr_bod(2,itotps))
      ALLOCATE(chpd(itotps,nalf,3),chpbd(itotps,nb,3),sab2(nb,nb,nat))
      ALLOCATE(scrall(3,itotps),chp3(itotps,nb),sab(nalf,nalf,nat))
      ALLOCATE(scr2(itotps))

      exch2=ZERO

      call ao_to_mo_grid(itotps,igr,nalf,c,chp,chp2)
      call ao_to_mo_grid(itotps,igr,nb,cb,chp,chp3)

!! rho (order: alpha, beta), reduced from the MO amplitudes just built.  !!
!! parallel over grid points: each k only reads its own chp2(k,:)/       !!
!! chp3(k,:) and writes only its own rho(:,k).                           !!
!$OMP PARALLEL DO PRIVATE(kk,jj,xx0,xx0b)
      do kk=1,itotps
        xx0=ZERO
        xx0b=ZERO
        do jj=1,nalf
          xx0=xx0+chp2(kk,jj)*chp2(kk,jj)
        end do
        do jj=1,nb
          xx0b=xx0b+chp3(kk,jj)*chp3(kk,jj)
        end do
        rho(1,kk)=xx0
        rho(2,kk)=xx0b
      end do
!$OMP END PARALLEL DO

      call print_box('UKS-DFT ENERGY DECOMPOSITION')
      write(*,*) " USING BOND ORDER DENSITY APPROACH "
      write(*,'(2x,a36,1x,f10.6)') "Threshold for atom pair calculation:",threbod

!! per-atom MO overlap, alpha (sab) and beta (sab2). parallel over       !!
!! iatom: each iteration only writes its own sab(:,:,iatom)/             !!
!! sab2(:,:,iatom), independent across atoms.                            !!
!$OMP PARALLEL DO PRIVATE(iatom,ii,jj,xx,xxb,ifut)
      do iatom=1,nat
        do ii=1,nalf
          do jj=ii,nalf
            xx=ZERO
            xxb=ZERO
            do ifut=1,itotps
              xx=xx+chp2(ifut,ii)*chp2(ifut,jj)*wp(ifut)*omp(ifut)*omp2(ifut,iatom)
              if(jj.le.nb.and.ii.le.nb) xxb=xxb+chp3(ifut,ii)*chp3(ifut,jj)*wp(ifut)*omp(ifut)*omp2(ifut,iatom)
            end do
            sab(ii,jj,iatom)=xx
            if(ii.ne.jj) sab(jj,ii,iatom)=xx
            if(jj.le.nb.and.ii.le.nb) then
              sab2(ii,jj,iatom)=xxb
              if(ii.ne.jj) sab2(jj,ii,iatom)=xxb
            end if
          end do
        end do
      end do
!$OMP END PARALLEL DO

!! sanity check: sab/sab2 summed over atoms must reproduce the MO        !!
!! overlap (identity for ii=jj, zero otherwise), alpha then beta.        !!
      do ii=1,nalf
        do jj=ii,nalf
          xa=ZERO
          if(ii.eq.jj) xa=-ONE
          do kk=1,nat
            xa=xa+sab(ii,jj,kk)
          end do
          if(abs(xa).gt.1.0d-3) then
            write(*,*) ii,jj,xa
            stop " PROBLEM WITH ALPHA MOs OVERLAPS "
          end if
        end do
      end do
      do ii=1,nb
        do jj=ii,nb
          xb=ZERO
          if(ii.eq.jj) xb=-ONE
          do kk=1,nat
            xb=xb+sab2(ii,jj,kk)
          end do
          if(abs(xb).gt.1.0d-3) then
            write(*,*) ii,jj,xb
            stop " PROBLEM WITH BETA MOs OVERLAPS "
          end if
        end do
      end do

      if(itype.gt.1) call ugrdrho(pcoord,chpd,chpbd)

      call print_box('BOND ORDER DENSITY FOR ALL ATOM PAIRS')
      write(*,*) " --------------------------- "
      write(*,*) "  Atom   Atom   BODEN value  "
      write(*,*) " --------------------------- "
      do iatom=1,nat
        do jatom=iatom+1,nat

!! BODEN for this atom pair, skipped below the THREBOD bond-order threshold. !!
          bx0=bo(iatom,jatom)
          if(bx0.ge.threbod) then
            xx=ZERO
            xxb=ZERO
!! MG: same false-sharing-free parallelization as the RHF twin
!! (numint_dft) above -- ff2/ff2b are genuine per-grid-point local scalars,
!! scr_bod(:,jfut) is written at a unique index per iteration, xx/xxb are
!! plain reductions. jfut computed inside the loop body, same reason as
!! the RHF twin above (gfortran's COLLAPSE needs loop-invariant bounds). !!
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,jloc,jfut,x2,wa,wb,ff2,ff2b,ii,jj,ff,ffb) REDUCTION(+:xx,xxb)
            do icenter=1,nat
              do jloc=1,iatps
                jfut=iatps*(icenter-1)+jloc
                x2=wp(jfut)*omp(jfut)
                wa=omp2(jfut,iatom)
                wb=omp2(jfut,jatom)
                ff2=ZERO
                ff2b=ZERO
                do ii=1,nalf
                  do jj=1,nalf
                    ff=sab(ii,jj,jatom)*wa+sab(jj,ii,iatom)*wb
                    ff2=ff2+ff*chp2(jfut,ii)*chp2(jfut,jj)
                    if(jj.le.nb.and.ii.le.nb) then
                      ffb=sab2(ii,jj,jatom)*wa+sab2(jj,ii,iatom)*wb
                      ff2b=ff2b+ffb*chp3(jfut,ii)*chp3(jfut,jj)
                    end if
                  end do
                end do

!! A-B BODEN defined as A-B + B-A, hence the factor of 2. scr_bod holds !!
!! the BODEN matrix, order alpha then beta. !!
                scr_bod(1,jfut)=ff2/TWO
                scr_bod(2,jfut)=ff2b/TWO
                if(iatom.ne.jatom) then
                  scr_bod(1,jfut)=ff2
                  scr_bod(2,jfut)=ff2b
                end if
                xx=xx+x2*scr_bod(1,jfut)
                xxb=xxb+x2*scr_bod(2,jfut)
              end do
            end do
!$OMP END PARALLEL DO
            write(*,'(4x,i3,4x,i3,4x,f10.7)') iatom,jatom,xx+xxb

!! BODEN gradient, for unrestricted GGA functionals. !!
            if(itype.gt.1) call grdboden_uks(itotps,chp2,chp3,omp2,chpd,chpbd,scrall,iatom,jatom,sab,sab2)
            call xc_uks(itotps,scr_bod,scrall,scr2)
            x1=ZERO
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,iloc,ifut) REDUCTION(+:x1)
            do icenter=1,nat
              do iloc=1,iatps
                ifut=iatps*(icenter-1)+iloc
                x1=x1+wp(ifut)*scr2(ifut)*omp2(ifut,icenter)*omp(ifut)
              end do
            end do
!$OMP END PARALLEL DO
            exch2(iatom,jatom)=x1
          end if
        end do
      end do
      write(*,*) " --------------------------- "

      call print_box('DIATOMIC PURE KS-DFT XC TERMS (BODEN)')
      exchen=ZERO
      do ii=1,nat
        exchen=exchen+exch2(ii,ii)
        do jj=ii,nat
          if(ii.ne.jj) then
            exch2(jj,ii)=exch2(ii,jj)
            exchen=exchen+exch2(ii,jj)
          end if
        end do
      end do
      CALL MPRINT2(exch2,nat,maxat)

      if(itype.gt.1) call sigma_uks(pcoord,chp2,chp3,scrall)

!! "exact" one-center XC energy, by direct integration (no BODEN approximation). !!
      call xc_uks(itotps,rho,scrall,scr2)
      xtot=ZERO
!! parallel over icenter: each iteration writes only its own exch(icenter,icenter), !!
!! xtot is a genuine running total. !!
!$OMP PARALLEL DO PRIVATE(icenter,x,ifut) REDUCTION(+:xtot)
      do icenter=1,nat
        x=ZERO
        do ifut=1,itotps
          x=x+wp(ifut)*scr2(ifut)*omp(ifut)*omp2(ifut,icenter)
        end do
        exch(icenter,icenter)=x
        xtot=xtot+x
      end do
!$OMP END PARALLEL DO

      call print_box('PURE KS-DFT XC ONE-CENTER TERMS (EXACT)')
      call MPRINT2(exch,nat,maxat)
      write(*,'(2x,a35,x,f14.7)') "KS-DFT exchange-correlation energy:",xtot
      write(*,*) " "

      write(*,*) " REARRANGING ATOMIC COMPONENTS "
      do ii=1,nat
        x0=ZERO
        do jj=1,nat
          if(ii.ne.jj) then
            exch(ii,jj)=exch2(ii,jj)
            x0=x0+exch(ii,jj)
          end if
        end do
        exch(ii,ii)=exch(ii,ii)-x0/TWO
      end do

      call print_box('FINAL PURE KS-DFT EXCHANGE-CORRELATION ENERGY COMPONENTS')
      exchen=ZERO
      do ii=1,nat
        do jj=ii,nat
          exchen=exchen+exch(ii,jj)
        end do
      end do
      call MPRINT2(exch,nat,maxat)
      write(*,'(2x,a47,x,f14.7)') "Sum of pure KS-DFT exchange-correlation energy:",exchen
      if(xmix.gt.ZERO) then
        write(*,*) " "
        write(*,*) " WARNING: HF-exchange part missing "
      end if
      if (idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Exc Decomposition'
        call group_by_frag_mat(1,line,exch)
      end if
      write(*,*) " "

      xtot=ZERO
      do ii=1,nat
        eto(ii,ii)=eto(ii,ii)+exch(ii,ii)
        xtot=xtot+eto(ii,ii)
        do jj=ii+1,nat
          eto(ii,jj)=eto(ii,jj)+exch(ii,jj)
          eto(jj,ii)=eto(ii,jj)
          xtot=xtot+eto(ii,jj)
        end do
      end do
      DEALLOCATE(chp2,rho,scr_bod,chpd,chpbd,sab2,scrall,chp3,sab,scr2)
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: grdboden_uks                                                !!
!! purpose: gradient of the bond-order density (BODEN) for one atom pair,  !!
!!   alpha and beta, at each grid point -- sigma input for the GGA/        !!
!!   hybrid-GGA XC evaluation in numint_dft_uks's BODEN loop. See          !!
!!   grdboden for restricted/closed-shell.                                 !!
!! arguments:                                                              !!
!!   itotps (in)  -- total number of grid points                          !!
!!   chp2   (in)  -- alpha MO values at each grid point                    !!
!!   chp3   (in)  -- beta MO values at each grid point                     !!
!!   omp2   (in)  -- becke/tfvc (or hirshfeld) weight of each point for    !!
!!                   every atom                                            !!
!!   chpd   (in)  -- alpha MO gradient at each grid point (from ugrdrho)   !!
!!   chpbd  (in)  -- beta MO gradient at each grid point (from ugrdrho)    !!
!!   scrall (out) -- sigma BODEN at each grid point, order                 !!
!!                   up-up/up-down/down-down                               !!
!!   iat,kat (in) -- the atom pair                                         !!
!!   sab,sab2 (in) -- per-atom MO overlap matrix, alpha/beta (from         !!
!!                   numint_dft_uks)                                       !!
!! author: MGimf.                                                          !!
!! *********************************************************************** !!
      subroutine grdboden_uks(itotps,chp2,chp3,omp2,chpd,chpbd,scrall,iat,kat,sab,sab2)
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chp2(itotps,nalf),chp3(itotps,nb),scrall(3,itotps)
      dimension chpd(itotps,nalf,3),chpbd(itotps,nb,3)
      dimension sab(nalf,nalf,nat),sab2(nb,nb,nat),omp2(itotps,nat)

      do kk=1,itotps
        scrall(1,kk)=ZERO
        scrall(2,kk)=ZERO
        scrall(3,kk)=ZERO
      end do

      do ixyz=1,3

!! sigma BODEN, order up-up/up-down/down-down. parallel over irun: each  !!
!! iteration only reads its own omp2(irun,:)/chpd(irun,:,ixyz)/          !!
!! chpbd(irun,:,ixyz)/chp2(irun,:)/chp3(irun,:) and writes its own       !!
!! scrall(:,irun) -- safe across the three ixyz passes since they run    !!
!! serially, only the point loop within each pass is threaded.           !!
!$OMP PARALLEL DO PRIVATE(irun,w1,w2,xxa,xxb,ii,jj,gab,gab2)
        do irun=1,itotps
          w1=omp2(irun,iat)
          w2=omp2(irun,kat)
          xxa=ZERO
          xxb=ZERO
          do ii=1,nalf
            do jj=1,nalf
              gab=sab(ii,jj,iat)*w2+sab(ii,jj,kat)*w1
              xxa=xxa+gab*chpd(irun,jj,ixyz)*chp2(irun,ii)
              if(jj.le.nb.and.ii.le.nb) then
                gab2=sab2(ii,jj,iat)*w2+sab2(ii,jj,kat)*w1
                xxb=xxb+gab2*chpbd(irun,jj,ixyz)*chp3(irun,ii)
              end if
            end do
          end do
          scrall(1,irun)=scrall(1,irun)+xxa*xxa
          scrall(2,irun)=scrall(2,irun)+xxa*xxb
          scrall(3,irun)=scrall(3,irun)+xxb*xxb
        end do
!$OMP END PARALLEL DO
      end do
      if(iat.ne.kat) then
        do kk=1,itotps
          scrall(1,kk)=FOUR*scrall(1,kk)
          scrall(2,kk)=FOUR*scrall(2,kk)
          scrall(3,kk)=FOUR*scrall(3,kk)
        end do
      end if
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: func_info_print                                             !!
!! purpose: prints the density-functional identification block (name,     !!
!!   references, exchange/correlation type, family) via libxc's own info   !!
!!   query functions, and classifies the functional into itype (1 LDA,    !!
!!   2 GGA, 3 meta-GGA) for the caller's gradient-calculation dispatch.    !!
!! arguments:                                                              !!
!!   id_func (in)  -- libxc functional id (from iopt, read by the caller)  !!
!!   itype   (out) -- functional family classification, see above          !!
!!   ifirst  (in)  -- 1 if this is the first of up to 3 consecutive calls  !!
!!                   from the caller (subbox, no leading blank -- the      !!
!!                   caller's own preceding box already supplies one), 0   !!
!!                   otherwise (print_box, own leading blank needed since  !!
!!                   the previous call's content isn't blank-terminated)   !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine func_info_print(id_func,itype,ifirst)

      use xc_f90_types_m
      use xc_f90_lib_m

      implicit real*8(a-h,o-z)

      TYPE(xc_f90_pointer_t) :: xc_func
      TYPE(xc_f90_pointer_t) :: xc_info

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common/exchg/exch(maxat,maxat),xmix

      character*120 name_ref
      character*80 name_func

      if(kop.ne.1) call xc_f90_func_init(xc_func,xc_info,id_func,XC_UNPOLARIZED)
      if(kop.eq.1) call xc_f90_func_init(xc_func,xc_info,id_func,XC_POLARIZED)
      call xc_f90_hyb_exx_coef(xc_func,xmix)
      call xc_f90_info_name(xc_info,name_func)

      if(ifirst.eq.1) then
        call print_subbox('DENSITY FUNCTIONAL INFORMATION')
      else
        call print_box('DENSITY FUNCTIONAL INFORMATION')
      end if
      write(*,*) " Functional name --> ",trim(name_func)
      ii=0
      call xc_f90_info_refs(xc_info,ii,name_ref)
      do while(ii.ge.0)
        write(*,'(2x,a15,i2,a1,x,a120)') "Reference --> [",ii,"]",name_ref
        call xc_f90_info_refs(xc_info,ii,name_ref)
      end do
      select case(xc_f90_info_kind(xc_info))
      case(XC_EXCHANGE)
        write(*,*) " Functional type --> exchange"
      case(XC_CORRELATION)
        write(*,*) " Functional type --> correlation"
      case(XC_EXCHANGE_CORRELATION)
        write(*,*) " Functional type --> exchange-correlation"
      case(XC_KINETIC)
        write(*,*) " Functional type --> kinetic energy"
      case default
        write(*,*) " Functional type --> unknown"
      end select

!! itype: 1 LDA, 2 GGA, 3 meta-GGA. Hybrid functionals get xmix>0 as well. !!
!! Caller stores the returned value in iopt(55) for later use. !!
      itype=0
      select case(xc_f90_info_family(xc_info))
      case(XC_FAMILY_UNKNOWN)
        write(*,*) " Family of the functional unknown"
        stop
      case(XC_FAMILY_NONE)
        write(*,*) " Non-Family identified functional"
        stop
      case(XC_FAMILY_LDA)
        write(*,*) " LDA functional selected"
        itype=1
      case(XC_FAMILY_GGA)
        write(*,*) " GGA functional selected"
        itype=2
      case(XC_FAMILY_HYB_GGA)
        write(*,*) " HYBRID-GGA functional selected"
        write(*,'(2x,a26,x,f5.3)') "HF-type exchange coeff -->",xmix
        itype=2
      case(XC_FAMILY_MGGA)
        write(*,*) " META-GGA functional selected"
        itype=3
        write(*,*) " META-GGA still in development!!!" !! MG: to-do !!
        stop
      case(XC_FAMILY_HYB_MGGA)
        write(*,*) " HYBRID-META-GGA functional selected"
        write(*,'(2x,a26,x,f5.3)') "HF-type exchange coeff -->",xmix
        itype=3
        write(*,*) " META-GGA still in development!!!" !! MG: to-do !!
        stop
      end select
      call xc_f90_func_end(xc_func)

      end 

!! ***** !!

!! *********************************************************************** !!
!! subroutine: sigma_uks                                                   !!
!! purpose: density-gradient contraction (sigma, order up-up/up-down/      !!
!!   down-down) at each grid point, for GGA/hybrid-GGA XC functionals.     !!
!!   Builds the AO gradient (x/y/z in turn) from primitives, transforms    !!
!!   to MOs (alpha and beta), and contracts with the already-transformed   !!
!!   chp2/chp3. See sigma for restricted/closed-shell.                     !!
!! arguments:                                                              !!
!!   pcoord (in)  -- xyz coordinates of each grid point                    !!
!!   chp2   (in)  -- alpha MO values at each grid point (from the caller)  !!
!!   chp3   (in)  -- beta MO values at each grid point (from the caller)   !!
!!   scr    (out) -- sigma at each grid point, order up-up/up-down/down-down !!
!! author: MGimf.                                                          !!
!! *********************************************************************** !!
      subroutine sigma_uks(pcoord,chp2,chp3,scr)
      use basis_set
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chp2(natoms*nrad*nang,nalf),scr(3,natoms*nang*nrad)
      dimension chp3(natoms*nrad*nang,nb)
      dimension pcoord(natoms*nang*nrad,3)
      allocatable:: chpd(:,:),chp(:,:),chpbd(:,:)

      ipoints=natoms*nrad*nang
      ALLOCATE(chpd(ipoints,nalf),chp(ipoints,igr),chpbd(ipoints,nb))

      do kk=1,ipoints
        do jj=1,3
          scr(jj,kk)=ZERO
        end do
      end do

      do ixyz=1,3

!! AO gradient (ixyz component) from primitives. parallel over irun: each !!
!! iteration writes only its own chp(irun,:), independent across points.  !!
!$OMP PARALLEL DO PRIVATE(irun,iact,iactat,x,y,z,rr,f,k,ipr,nn,ll,mm,alpha,dx)
        do irun=1,ipoints
          do iact=1,nbasis
            iactat=ihold(iact)
            x=pcoord(irun,1)-coord(1,iactat)
            y=pcoord(irun,2)-coord(2,iactat)
            z=pcoord(irun,3)-coord(3,iactat)
            rr=dsqrt(x**2.0d0+y**2.0d0+z**2.0d0)
            f=0.d0
            k=1
            do while(nprimbas(k,iact).ne.0)
              ipr=nprimbas(k,iact)
              nn=nlm(ipr,1)
              ll=nlm(ipr,2)
              mm=nlm(ipr,3)
              alpha=expp(ipr)
              if(ixyz.eq.1) then
                dx=-2.0d0*alpha*(x**(nn+1))
                if(nn.ge.1) dx=dx+nn*x**(nn-1)
                dx=dx*(y**ll)*(z**mm)*dexp(-expp(ipr)*(rr**2))
              else if(ixyz.eq.2) then
                dx=-2.0d0*alpha*(y**(ll+1))
                if(ll.ge.1) dx=dx+ll*y**(ll-1)
                dx=dx*(x**nn)*(z**mm)*dexp(-expp(ipr)*(rr**2))
              else
                dx=-2.0d0*alpha*(z**(mm+1))
                if(mm.ge.1) dx=dx+mm*z**(mm-1)
                dx=dx*(x**nn)*(y**ll)*dexp(-expp(ipr)*(rr**2))
              end if
              f=f+dx*coefpb(ipr,iact)
              k=k+1
            enddo
            chp(irun,iact)=f
          enddo
        enddo
!$OMP END PARALLEL DO

        call ao_to_mo_grid(ipoints,igr,nalf,c,chp,chpd)
        call ao_to_mo_grid(ipoints,igr,nb,cb,chp,chpbd)

!! sigma for libxc, order up-up/up-down/down-down. parallel over k: each !!
!! iteration only reads chp2(k,:)/chp3(k,:)/chpd(k,:)/chpbd(k,:) and     !!
!! writes its own scr(:,k) -- safe across the three ixyz passes since    !!
!! they run serially, only the point loop within each pass is threaded. !!
!$OMP PARALLEL DO PRIVATE(k,i,j)
        do k=1,ipoints
          do i=1,nalf
            do j=1,nalf
              scr(1,k)=scr(1,k)+chp2(k,i)*chp2(k,j)*chpd(k,i)*chpd(k,j)
              if(j.le.nb) scr(2,k)=scr(2,k)+chp2(k,i)*chp3(k,j)*chpd(k,i)*chpbd(k,j)
              if(i.le.nb.and.j.le.nb) scr(3,k)=scr(3,k)+chp3(k,i)*chp3(k,j)*chpbd(k,i)*chpbd(k,j)
            end do
          end do
        end do
!$OMP END PARALLEL DO
      end do

      do k=1,ipoints
        do j=1,3
          scr(j,k)=4.0d0*scr(j,k)
        end do
      end do
      DEALLOCATE(chpd,chpbd,chp)
      end

!! ***** !!
