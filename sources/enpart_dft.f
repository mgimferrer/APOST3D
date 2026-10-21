!! ************************************************************************** !!
!! KS-DFT exchange-correlation energy partitioning (ENPART)                   !!
!! RKS (closed-shell):                                                        !!
!!   numint_dft      -- BODEN approximation + "exact" one-center XC term      !!
!!   xc              -- XC energy density at each grid point, via libxc       !!
!!   sigma           -- density-gradient contraction, GGA and meta-GGA only   !!
!!   grdrho          -- MO gradient at each grid point, GGA and meta-GGA only !!
!!   grdboden        -- BODEN gradient for one atom pair, GGA and meta-GGA    !!
!!   tauboden        -- BODEN tau for one atom pair, meta-GGA only            !!
!! UKS (open-shell):                                                          !!
!!   numint_dft_uks  -- UKS twin of numint_dft                                !!
!!   xc_uks          -- UKS twin of xc                                        !!
!!   sigma_uks       -- UKS twin of sigma                                     !!
!!   ugrdrho         -- UKS twin of grdrho                                    !!
!!   grdboden_uks    -- UKS twin of grdboden                                  !!
!!   tauboden_uks    -- UKS twin of tauboden                                  !!
!! shared, authored PSalse, MGimf.:                                           !!
!!   func_info_print -- functional identification/classification, called once !!
!!                      by main.f before dispatching to either half above     !!
!! exploration of the meta-GGA pair tau (TAU_BODEN), to be reduced to one:    !!
!!   tau_diagnostics(_uks), print_tau_diagnostics                             !!
!! ************************************************************************** !!

!! ***** !!

!! *********************************************************************** !!
!! subroutine: xc                                                          !!
!! purpose: restricted XC energy density at each grid point, via libxc.    !!
!!   Exchange and correlation are queried separately when the input        !!
!!   functional specifies them as distinct libxc ids (id_xfunc/id_cfunc),  !!
!!   otherwise as one combined exchange-correlation id (id_xcfunc). The    !!
!!   density-Laplacian ingredient libxc's MGGA call always requires is     !!
!!   passed as a local zero array -- confirmed against libxc 7.1.2's own   !!
!!   functional flags that M06-2X/TPSS (the only MGGA family wired in so   !!
!!   far) declare XC_FLAGS_NEEDS_TAU only, never XC_FLAGS_NEEDS_LAPLACIAN. !!
!!   See xc_uks for open-shell.                                            !!
!! arguments:                                                              !!
!!   npt    (in)  -- number of grid points                                 !!
!!   scr_a  (in)  -- electron density at each point                        !!
!!   scr    (in)  -- density gradient contraction (sigma) at each point,   !!
!!                   only read for GGA/hybrid-GGA/meta-GGA functionals     !!
!!   tau    (in)  -- kinetic energy density at each point, only read for   !!
!!                   meta-GGA functionals                                  !!
!!   scr2   (out) -- XC energy density at each point (already multiplied   !!
!!                   by the density)                                       !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine xc(npt,scr_a,scr,tau,scr2)
      use xc_f03_lib_m
      implicit real*8(a-h,o-z)
      TYPE(xc_f03_func_t) :: xc_func
      TYPE(xc_f03_func_info_t) :: xc_info
      include 'parameter.h'
      integer npt
      real*8 scr_a(npt),scr(npt),tau(npt),scr2(npt)
      common /iops/iopt(200)
      allocatable :: scr2c(:),xlapl(:)

      id_xcfunc = iopt(60)
      id_cfunc  = iopt(62)
      id_xfunc  = iopt(61)

      ALLOCATE(scr2c(npt),xlapl(npt))
      do ii=1,npt
        scr2c(ii)=ZERO
        xlapl(ii)=ZERO
      end do

!! exchange-correlation, from libxc. !!
      if(id_xcfunc.ne.0) then
        call xc_f03_func_init(xc_func,id_xcfunc,XC_UNPOLARIZED)
        xc_info = xc_f03_func_get_info(xc_func)
        select case (xc_f03_func_info_get_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f03_lda_exc(xc_func,int(npt,8),scr_a(1),scr2(1))
        case(XC_FAMILY_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_a(1),scr(1),scr2(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_a(1),scr(1),scr2(1))
        case(XC_FAMILY_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_a(1),scr(1),xlapl(1),tau(1),scr2(1))
        case(XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_a(1),scr(1),xlapl(1),tau(1),scr2(1))
        end select
        call xc_f03_func_end(xc_func)
      end if

!! correlation, when specified as a separate libxc id. !!
      if(id_cfunc.ne.0) then
        call xc_f03_func_init(xc_func,id_cfunc,XC_UNPOLARIZED)
        xc_info = xc_f03_func_get_info(xc_func)
        select case (xc_f03_func_info_get_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f03_lda_exc(xc_func,int(npt,8),scr_a(1),scr2c(1))
        case(XC_FAMILY_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_a(1),scr(1),scr2c(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_a(1),scr(1),scr2c(1))
        case(XC_FAMILY_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_a(1),scr(1),xlapl(1),tau(1),scr2c(1))
        case(XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_a(1),scr(1),xlapl(1),tau(1),scr2c(1))
        end select
        call xc_f03_func_end(xc_func)
      end if

!! exchange, when specified as a separate libxc id. !!
      if(id_xfunc.ne.0) then
        call xc_f03_func_init(xc_func,id_xfunc,XC_UNPOLARIZED)
        xc_info = xc_f03_func_get_info(xc_func)
        select case (xc_f03_func_info_get_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f03_lda_exc(xc_func,int(npt,8),scr_a(1),scr2(1))
        case(XC_FAMILY_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_a(1),scr(1),scr2(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_a(1),scr(1),scr2(1))
        case(XC_FAMILY_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_a(1),scr(1),xlapl(1),tau(1),scr2(1))
        case(XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_a(1),scr(1),xlapl(1),tau(1),scr2(1))
        end select
        call xc_f03_func_end(xc_func)
      end if

!! combine exchange and correlation, multiply by the density. !!
      do ii=1,npt
        scr2(ii)=(scr2c(ii)+scr2(ii))*scr_a(ii)
      end do
      DEALLOCATE(scr2c,xlapl)

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
      use ao_matrices
      use integration_grid
      use output_mod, only: print_value
      use input_options_mod, only: iverbose,itauboden
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
      dimension exch2(maxat,maxat),Excmp(maxat,maxat)
      dimension sat(igr,igr,nat)
      character*80 c1,c2
      dimension xoc(maxat,2)
      allocatable :: ocm(:,:)

      allocatable :: chp2(:,:),scr(:),scrx(:,:)
      allocatable :: scr_a(:),rho_a(:),scr2(:),scr_tau(:),tau_ab(:)
      allocatable :: sab(:,:,:),chpd(:,:,:)
      dimension clipab(maxat,maxat),zab(maxat,maxat)

      ithrebod =  Iopt(44)
      itype    =  Iopt(55)
      iatps    =  nang*nrad

      ALLOCATE(chp2(itotps,nocc),scr(itotps))
      ALLOCATE(scr_a(itotps),rho_a(itotps),scr2(itotps),scr_tau(itotps))
      ALLOCATE(tau_ab(itotps))
      clipab=ZERO
      zab=ZERO
      ALLOCATE(sab(nocc,nocc,nat),chpd(itotps,nocc,3))
      ALLOCATE(scrx(igr,nocc))

      exch2=ZERO

!! threshold for BODEN calculation, from the bond order between a pair of atoms. !!
      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if


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
            write(*,'(2x,a,i0,a,i0,a,es10.2)') 'Atomic overlaps of the MOs ',
     +        ii,' and ',jj,' do not add up to the MO overlap: ',x
            call apost_stop(" PROBLEM WITH MOs OVERLAPS ")
          end if
        end do
      end do

      call ao_to_mo_grid(itotps,igr,nocc,c,chp,chp2)

      if(itype.gt.1) call grdrho(pcoord,chpd)

!! kinetic energy density (tau), meta-GGA only -- reduction over chpd,     !!
!! already computed by grdrho above, no new AO/derivative evaluation.      !!
!! Total tau for doubly occupied MOs: 2*(1/2)*sum|grad phi|^2, same total- !!
!! density convention as sigma. parallel over k: each iteration only reads !!
!! its own chpd(k,:,:) and writes its own scr_tau(k).                      !!
      if(itype.eq.3) then
!$OMP PARALLEL DO PRIVATE(k,ii,ixyz,x)
        do k=1,itotps
          x=ZERO
          do ii=1,nocc
            do ixyz=1,3
              x=x+chpd(k,ii,ixyz)*chpd(k,ii,ixyz)
            end do
          end do
          scr_tau(k)=x
        end do
!$OMP END PARALLEL DO
      end if

!! pairs below the THREBOD bond-order threshold get no BODEN: their whole !!
!! xc term comes from the multipolar expansion of the KS determinant's    !!
!! exchange density (closed shell: -2 phi_i^2, -4 phi_i*phi_j), as for HF. !!
      nskip=0
      do iatom=1,nat
        do jatom=iatom+1,nat
          if(bo(iatom,jatom).lt.threbod) nskip=nskip+1
        end do
      end do
      if(nskip.gt.0) call multipolar_mo(nocc,itotps,wp,omp2,pcoord,chp2,-TWO,-FOUR,Excmp)

      npairs=nat*(nat-1)/2-nskip
      call print_box('DFT EXCHANGE-CORRELATION (BODEN)')
      call print_value('Pairs from the bond order density',npairs)
      call print_value('Pairs from the multipolar expansion',nskip)
      if(nskip.gt.0) call print_skipped_pairs(xmix)
      if(iverbose.eq.1) then
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Integral of the bond order density (BODEN) of'//
     +    ' each atom pair:'
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Atom pair       BODEN'
        write(*,'(2x,a)') repeat('-',21)
      end if
      do iatom=1,nat
        do jatom=iatom+1,nat

!! BODEN for this atom pair, skipped below the THREBOD bond-order threshold. !!
          bx0=bo(iatom,jatom)
          if(bx0.lt.threbod) then
!! a hybrid's xmix share is added with the HF-type exchange (numint_two). !!
            exch2(iatom,jatom)=(ONE-xmix)*Excmp(iatom,jatom)
            if(iverbose.eq.1) write(*,'(2x,i3,a2,i3,a13)') iatom,' -',
     +        jatom,'multipolar'
          else
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

!! the A-B BODEN counts A-B and B-A, hence the factor of 2 !!
                scr_a(jfut)=ff2
                if(iatom.ne.jatom) scr_a(jfut)=TWO*ff2
                x1=x1+x2*scr_a(jfut)
              end do
            end do
!$OMP END PARALLEL DO
            if(iverbose.eq.1) write(*,'(2x,i3,a2,i3,f13.7)') iatom,' -',
     +        jatom,x1

!! BODEN gradient, for GGA functionals -- scr_a: BODEN for this atom pair, !!
!! scr: sigma BODEN, tau_ab: tau BODEN (meta-GGA), scr2: XC functional   !!
!! value. tau_ab is only read by xc for meta-GGA functionals.            !!
            if(itype.gt.1) call grdboden(itotps,omp2,chp2,chpd,scr,iatom,jatom,sab)
            if(itype.eq.3) then
              call tauboden(itotps,omp2,chpd,scr_a,scr,scr_tau,tau_ab,
     +          iatom,jatom,sab)
              call tau_diagnostics(itotps,wp,omp,scr_a,scr,tau_ab,
     +          clipab(iatom,jatom),zab(iatom,jatom))
            end if
            call xc(itotps,scr_a,scr,tau_ab,scr2)

            x1=ZERO
!! xc energy of the pair on the grid: x1 is a reduction, nothing else !!
!! is written                                                         !!
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
      if(iverbose.eq.1) write(*,'(2x,a)') repeat('-',21)
      if(itype.eq.3) call print_tau_diagnostics(threbod,clipab,zab)

      exchen=ZERO
      do ii=1,nat
        do jj=ii,nat
          if(ii.ne.jj) exch2(jj,ii)=exch2(ii,jj)
          exchen=exchen+exch2(ii,jj)
        end do
      end do
      if(iverbose.eq.1) then
        call enpart_matrix('DFT EXCHANGE-CORRELATION, PAIR TERMS',
     +    'Off-diagonal: the pair terms (BODEN, or the multipolar'//
     +    ' expansion below',
     +    'the THREBOD threshold). No one-center terms here.',exch2,
     +    'DFT pair terms',exchen)
        call enpart_frag('DFT EXCHANGE-CORRELATION, PAIR TERMS',exch2)
      end if

      if(itype.gt.1) call sigma(pcoord,chp2,scr)

!! "exact" one-center XC energy, by direct integration (no BODEN approximation). !!
      call xc(itotps,rho,scr,scr_tau,scr2)
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

!! one-center terms by direct integration, kept for the tables below !!
      ALLOCATE(ocm(maxat,maxat))
      ocm=ZERO
      do ii=1,nat
        xoc(ii,1)=exch(ii,ii)
        ocm(ii,ii)=exch(ii,ii)
      end do
      if(iverbose.eq.1) then
        call enpart_matrix('DFT EXCHANGE-CORRELATION, ONE-CENTER TERMS',
     +    'Diagonal: the functional integrated within each atom, before'//
     +    ' half of',
     +    'each pair term is moved out of it.',ocm,
     +    'DFT one-center terms',xtot)
        call enpart_frag('DFT EXCHANGE-CORRELATION, ONE-CENTER TERMS',
     +    ocm)
      end if
      DEALLOCATE(ocm)

!! self term: the one-center term minus half of each pair term of the atom !!
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

      do ii=1,nat
        xoc(ii,2)=exch(ii,ii)
      end do
      call enpart_atom_table('DFT EXCHANGE-CORRELATION BY ATOM',
     +  'In the atom: the functional integrated within the atom. Self'//
     +  ' term: the same',
     +  'minus half of each pair term of the atom (the diagonal below).',
     +  2,xoc,'In the atom','Self term',' ',.true.)

      exchen=ZERO
      do ii=1,nat
        do jj=ii,nat
          exchen=exchen+exch(ii,jj)
        end do
      end do
      if(xmix.gt.ZERO) then
        write(c1,'(a,f5.3,a)') 'DFT part only, without the ',xmix,
     +    ' x HF-type exchange: the complete'
        c2='terms are in EXCHANGE-CORRELATION ENERGY, further below.'
      else
        c1='The whole exchange-correlation energy (pure functional).'
        c2=' '
      end if
      call enpart_matrix('DFT EXCHANGE-CORRELATION ENERGY',c1,c2,exch,
     +  'DFT exchange-correlation',exchen)
      call enpart_frag('DFT EXCHANGE-CORRELATION ENERGY',exch)
      call json_enpart_term('xc_dft','DFT exchange-correlation',exch,
     +  exchen,.false.,0.0d0)

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
      DEALLOCATE(chp2,scr,scr_a,rho_a,scr2,scr_tau,sab,chpd,scrx,tau_ab)
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
!!   otherwise as one combined exchange-correlation id (id_xcfunc). The    !!
!!   density-Laplacian ingredient libxc's MGGA call always requires is     !!
!!   passed as a local zero array -- confirmed against libxc 7.1.2's own   !!
!!   functional flags that M06-2X/TPSS (the only MGGA family wired in so   !!
!!   far) declare XC_FLAGS_NEEDS_TAU only, never XC_FLAGS_NEEDS_LAPLACIAN. !!
!!   See xc for restricted/closed-shell.                                   !!
!! arguments:                                                              !!
!!   npt    (in)  -- number of grid points                                 !!
!!   scr_ab (in)  -- electron density at each point, order alpha then beta !!
!!   scr    (in)  -- density gradient contraction (sigma) at each point,   !!
!!                   order up-up/up-down/down-down, only read for GGA/     !!
!!                   hybrid-GGA/meta-GGA functionals                       !!
!!   tau    (in)  -- kinetic energy density at each point, order alpha     !!
!!                   then beta, only read for meta-GGA functionals         !!
!!   scr2   (out) -- XC energy density at each point (already multiplied   !!
!!                   by the density)                                       !!
!! author: MGimf.                                                          !!
!! *********************************************************************** !!
      subroutine xc_uks(npt,scr_ab,scr,tau,scr2)
      use xc_f03_lib_m
      implicit real*8(a-h,o-z)
      TYPE(xc_f03_func_t) :: xc_func
      TYPE(xc_f03_func_info_t) :: xc_info
      include 'parameter.h'
      integer npt
      common /iops/iopt(200)
      dimension scr_ab(2,npt),scr(3,npt),tau(2,npt),scr2(npt)
      allocatable :: scr2c(:),xlapl(:,:)

      id_xcfunc = iopt(60)
      id_cfunc  = iopt(62)
      id_xfunc  = iopt(61)

      ALLOCATE(scr2c(npt),xlapl(2,npt))
      do ii=1,npt
        scr2c(ii)=ZERO
        xlapl(1,ii)=ZERO
        xlapl(2,ii)=ZERO
      end do

!! exchange-correlation, from libxc. !!
      if(id_xcfunc.ne.0) then
        call xc_f03_func_init(xc_func,id_xcfunc,XC_POLARIZED)
        xc_info = xc_f03_func_get_info(xc_func)
        select case (xc_f03_func_info_get_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f03_lda_exc(xc_func,int(npt,8),scr_ab(1,1),scr2(1))
        case(XC_FAMILY_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),scr2(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),scr2(1))
        case(XC_FAMILY_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),xlapl(1,1),tau(1,1),scr2(1))
        case(XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),xlapl(1,1),tau(1,1),scr2(1))
        end select
        call xc_f03_func_end(xc_func)
      end if

!! correlation, when specified as a separate libxc id. !!
      if(id_cfunc.ne.0) then
        call xc_f03_func_init(xc_func,id_cfunc,XC_POLARIZED)
        xc_info = xc_f03_func_get_info(xc_func)
        select case (xc_f03_func_info_get_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f03_lda_exc(xc_func,int(npt,8),scr_ab(1,1),scr2c(1))
        case(XC_FAMILY_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),scr2c(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),scr2c(1))
        case(XC_FAMILY_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),xlapl(1,1),tau(1,1),scr2c(1))
        case(XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),xlapl(1,1),tau(1,1),scr2c(1))
        end select
        call xc_f03_func_end(xc_func)
      end if

!! exchange, when specified as a separate libxc id. !!
      if(id_xfunc.ne.0) then
        call xc_f03_func_init(xc_func,id_xfunc,XC_POLARIZED)
        xc_info = xc_f03_func_get_info(xc_func)
        select case (xc_f03_func_info_get_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f03_lda_exc(xc_func,int(npt,8),scr_ab(1,1),scr2(1))
        case(XC_FAMILY_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),scr2(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f03_gga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),scr2(1))
        case(XC_FAMILY_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),xlapl(1,1),tau(1,1),scr2(1))
        case(XC_FAMILY_HYB_MGGA)
          call xc_f03_mgga_exc(xc_func,int(npt,8),scr_ab(1,1),scr(1,1),xlapl(1,1),tau(1,1),scr2(1))
        end select
        call xc_f03_func_end(xc_func)
      end if

!! combine exchange and correlation, multiply by the density. !!
      do ii=1,npt
        scr2(ii)=(scr2c(ii)+scr2(ii))*(scr_ab(1,ii)+scr_ab(2,ii))
      end do
      DEALLOCATE(scr2c,xlapl)
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
      use output_mod, only: print_value
      use input_options_mod, only: iverbose,itauboden
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(200)
      common/actual/iact,jat,icenter
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/exchg/exch(maxat,maxat),xmix
      character*80 c1,c2
      dimension xoc(maxat,2)
      allocatable :: ocm(:,:)
      dimension eto(maxat,maxat)
      dimension wp(itotps),chp(itotps,ndim)
      dimension omp(itotps),omp2(itotps,nat),pcoord(itotps,3)
      dimension exch2(maxat,maxat),Excmp(maxat,maxat),Excmpb(maxat,maxat)
      allocatable :: chp2(:,:),chp3(:,:),rho(:,:),scrall(:,:),scr2(:)
      allocatable :: chpd(:,:,:),sab(:,:,:),sab2(:,:,:),scr_bod(:,:)
      allocatable :: chpbd(:,:,:),tau(:,:),tau_ab(:,:)
      dimension clipab(maxat,maxat),zab(maxat,maxat)

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
      ALLOCATE(scr2(itotps),tau(2,itotps),tau_ab(2,itotps))
      clipab=ZERO
      zab=ZERO

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
            write(*,'(2x,a,i0,a,i0,a,es10.2)') 'Atomic overlaps of the alpha MOs ',
     +        ii,' and ',jj,' do not add up to the MO overlap: ',xa
            call apost_stop(" PROBLEM WITH ALPHA MOs OVERLAPS ")
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
            write(*,'(2x,a,i0,a,i0,a,es10.2)') 'Atomic overlaps of the beta MOs ',
     +        ii,' and ',jj,' do not add up to the MO overlap: ',xb
            call apost_stop(" PROBLEM WITH BETA MOs OVERLAPS ")
          end if
        end do
      end do

      if(itype.gt.1) call ugrdrho(pcoord,chpd,chpbd)

!! kinetic energy density (tau), alpha/beta, meta-GGA only -- reduction    !!
!! over chpd/chpbd, already computed by ugrdrho above, no new AO/          !!
!! derivative evaluation. parallel over k: each iteration only reads its   !!
!! own chpd(k,:,:)/chpbd(k,:,:) and writes its own tau(:,k).                !!
      if(itype.eq.3) then
!$OMP PARALLEL DO PRIVATE(k,ii,ixyz,xa,xb)
        do k=1,itotps
          xa=ZERO
          xb=ZERO
          do ii=1,nalf
            do ixyz=1,3
              xa=xa+chpd(k,ii,ixyz)*chpd(k,ii,ixyz)
            end do
          end do
          do ii=1,nb
            do ixyz=1,3
              xb=xb+chpbd(k,ii,ixyz)*chpbd(k,ii,ixyz)
            end do
          end do
          tau(1,k)=HALF*xa
          tau(2,k)=HALF*xb
        end do
!$OMP END PARALLEL DO
      end if

!! pairs below the THREBOD bond-order threshold get no BODEN: their whole !!
!! xc term comes from the multipolar expansion of the KS determinant's    !!
!! exchange density (per spin: -1 phi_i^2, -2 phi_i*phi_j), as for UHF.   !!
      nskip=0
      do iatom=1,nat
        do jatom=iatom+1,nat
          if(bo(iatom,jatom).lt.threbod) nskip=nskip+1
        end do
      end do
      if(nskip.gt.0) then
        call multipolar_mo(nalf,itotps,wp,omp2,pcoord,chp2,-ONE,-TWO,Excmp)
        call multipolar_mo(nb,itotps,wp,omp2,pcoord,chp3,-ONE,-TWO,Excmpb)
      end if

      npairs=nat*(nat-1)/2-nskip
      call print_box('DFT EXCHANGE-CORRELATION (BODEN)')
      call print_value('Pairs from the bond order density',npairs)
      call print_value('Pairs from the multipolar expansion',nskip)
      if(nskip.gt.0) call print_skipped_pairs(xmix)
      if(iverbose.eq.1) then
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Integral of the bond order density (BODEN) of'//
     +    ' each atom pair:'
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Atom pair       BODEN'
        write(*,'(2x,a)') repeat('-',21)
      end if
      do iatom=1,nat
        do jatom=iatom+1,nat

!! BODEN for this atom pair, skipped below the THREBOD bond-order threshold. !!
          bx0=bo(iatom,jatom)
          if(bx0.lt.threbod) then
!! a hybrid's xmix share is added with the HF-type exchange (numint_two_uhf). !!
            exch2(iatom,jatom)=(ONE-xmix)*(Excmp(iatom,jatom)+Excmpb(iatom,jatom))
            if(iverbose.eq.1) write(*,'(2x,i3,a2,i3,a13)') iatom,' -',
     +        jatom,'multipolar'
          else
            xx=ZERO
            xxb=ZERO
!! parallel as in the RHF twin (numint_dft): ff2/ff2b are PRIVATE       !!
!! per-point scalars, scr_bod(:,jfut) is written at a unique index per  !!
!! iteration, xx/xxb are reductions; jfut is computed in the body       !!
!! because COLLAPSE needs loop-invariant bounds.                        !!
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
            if(iverbose.eq.1) write(*,'(2x,i3,a2,i3,f13.7)') iatom,' -',
     +        jatom,xx+xxb

!! BODEN gradient, for unrestricted GGA functionals; tau BODEN for      !!
!! meta-GGA functionals (tau_ab, only read by xc_uks for meta-GGAs).     !!
            if(itype.gt.1) call grdboden_uks(itotps,chp2,chp3,omp2,chpd,chpbd,scrall,iatom,jatom,sab,sab2)
            if(itype.eq.3) then
              call tauboden_uks(itotps,omp2,chpd,chpbd,scr_bod,scrall,
     +          tau,tau_ab,iatom,jatom,sab,sab2)
              call tau_diagnostics_uks(itotps,wp,omp,scr_bod,scrall,
     +          tau_ab,clipab(iatom,jatom),zab(iatom,jatom))
            end if
            call xc_uks(itotps,scr_bod,scrall,tau_ab,scr2)
            x1=ZERO
!! xc energy of the pair on the grid: x1 is a reduction, nothing else !!
!! is written                                                         !!
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
      if(iverbose.eq.1) write(*,'(2x,a)') repeat('-',21)
      if(itype.eq.3) call print_tau_diagnostics(threbod,clipab,zab)

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
      if(iverbose.eq.1) then
        call enpart_matrix('DFT EXCHANGE-CORRELATION, PAIR TERMS',
     +    'Off-diagonal: the pair terms (BODEN, or the multipolar'//
     +    ' expansion below',
     +    'the THREBOD threshold). No one-center terms here.',exch2,
     +    'DFT pair terms',exchen)
        call enpart_frag('DFT EXCHANGE-CORRELATION, PAIR TERMS',exch2)
      end if

      if(itype.gt.1) call sigma_uks(pcoord,chp2,chp3,scrall)

!! "exact" one-center XC energy, by direct integration (no BODEN approximation). !!
      call xc_uks(itotps,rho,scrall,tau,scr2)
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

!! one-center terms by direct integration, kept for the tables below !!
      ALLOCATE(ocm(maxat,maxat))
      ocm=ZERO
      do ii=1,nat
        xoc(ii,1)=exch(ii,ii)
        ocm(ii,ii)=exch(ii,ii)
      end do
      if(iverbose.eq.1) then
        call enpart_matrix('DFT EXCHANGE-CORRELATION, ONE-CENTER TERMS',
     +    'Diagonal: the functional integrated within each atom, before'//
     +    ' half of',
     +    'each pair term is moved out of it.',ocm,
     +    'DFT one-center terms',xtot)
        call enpart_frag('DFT EXCHANGE-CORRELATION, ONE-CENTER TERMS',
     +    ocm)
      end if
      DEALLOCATE(ocm)

!! self term: the one-center term minus half of each pair term of the atom !!
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

      do ii=1,nat
        xoc(ii,2)=exch(ii,ii)
      end do
      call enpart_atom_table('DFT EXCHANGE-CORRELATION BY ATOM',
     +  'In the atom: the functional integrated within the atom. Self'//
     +  ' term: the same',
     +  'minus half of each pair term of the atom (the diagonal below).',
     +  2,xoc,'In the atom','Self term',' ',.true.)

      exchen=ZERO
      do ii=1,nat
        do jj=ii,nat
          exchen=exchen+exch(ii,jj)
        end do
      end do
      if(xmix.gt.ZERO) then
        write(c1,'(a,f5.3,a)') 'DFT part only, without the ',xmix,
     +    ' x HF-type exchange: the complete'
        c2='terms are in EXCHANGE-CORRELATION ENERGY, further below.'
      else
        c1='The whole exchange-correlation energy (pure functional).'
        c2=' '
      end if
      call enpart_matrix('DFT EXCHANGE-CORRELATION ENERGY',c1,c2,exch,
     +  'DFT exchange-correlation',exchen)
      call enpart_frag('DFT EXCHANGE-CORRELATION ENERGY',exch)
      call json_enpart_term('xc_dft','DFT exchange-correlation',exch,
     +  exchen,.false.,0.0d0)

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
      DEALLOCATE(chp2,rho,scr_bod,chpd,chpbd,sab2,scrall,chp3,sab,scr2,tau)
      DEALLOCATE(tau_ab)
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
!! subroutine: tauboden                                                    !!
!! purpose: kinetic energy density (tau) of the bond-order density (BODEN) !!
!!   for one atom pair, at each grid point -- tau input for the meta-GGA   !!
!!   XC evaluation in numint_dft's BODEN loop. EXPLORATION: TAU_BODEN      !!
!!   (# ENPART) chooses the definition, to be reduced to one:              !!
!!   1 kernel: tau_AB = sum_ij K_ij grad(phi_i).grad(phi_j), the tau of    !!
!!     the same frozen-weight kernel K_ij = S^B_ij w_A + S^A_ij w_B that   !!
!!     gives the BODEN and grdboden's sigma (closed shell: 2 x 1/2);       !!
!!   2 von Weizsacker: tau_AB = sigma_AB/(8 rho_AB);                        !!
!!   3 weights: tau_AB = 2 w_A w_B tau (A-B plus B-A).                      !!
!!   See tauboden_uks for open-shell.                                      !!
!! arguments:                                                              !!
!!   itotps (in)  -- total number of grid points                          !!
!!   omp2   (in)  -- atomic weight of each point for every atom            !!
!!   chpd   (in)  -- MO gradient at each grid point (from grdrho)          !!
!!   rhoab  (in)  -- BODEN of the pair at each grid point                  !!
!!   sigab  (in)  -- sigma BODEN of the pair (from grdboden)               !!
!!   taumol (in)  -- molecular tau at each grid point                      !!
!!   tauab  (out) -- tau BODEN of the pair at each grid point              !!
!!   iat,kat (in) -- the atom pair                                         !!
!!   sab    (in)  -- per-atom MO overlap matrix (from numint_dft)          !!
!! author: MGimf.                                                          !!
!! *********************************************************************** !!
      subroutine tauboden(itotps,omp2,chpd,rhoab,sigab,taumol,tauab,
     +  iat,kat,sab)
      use input_options_mod, only: itauboden
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chpd(itotps,nocc,3),omp2(itotps,nat)
      dimension rhoab(itotps),sigab(itotps),taumol(itotps),tauab(itotps)
      dimension sab(nocc,nocc,nat)

      select case(itauboden)
      case(1)
!! parallel over irun: each iteration only reads its own omp2(irun,:)/  !!
!! chpd(irun,:,:) and writes its own tauab(irun). K is symmetric, so    !!
!! the off-diagonal ij terms are counted twice.                          !!
!$OMP PARALLEL DO PRIVATE(irun,w1,w2,xx,ii,jj,gab,dd)
        do irun=1,itotps
          w1=omp2(irun,iat)
          w2=omp2(irun,kat)
          xx=ZERO
          do ii=1,nocc
            do jj=ii,nocc
              gab=sab(ii,jj,iat)*w2+sab(ii,jj,kat)*w1
              dd=chpd(irun,ii,1)*chpd(irun,jj,1)
     +          +chpd(irun,ii,2)*chpd(irun,jj,2)
     +          +chpd(irun,ii,3)*chpd(irun,jj,3)
              if(ii.ne.jj) dd=TWO*dd
              xx=xx+gab*dd
            end do
          end do
          tauab(irun)=xx
        end do
!$OMP END PARALLEL DO
      case(2)
        do irun=1,itotps
          tauab(irun)=ZERO
          if(rhoab(irun).gt.ZERO) tauab(irun)=sigab(irun)/(8.0d0*rhoab(irun))
        end do
      case(3)
        do irun=1,itotps
          tauab(irun)=TWO*omp2(irun,iat)*omp2(irun,kat)*taumol(irun)
        end do
      end select
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: tauboden_uks                                                !!
!! purpose: UKS twin of tauboden, tau BODEN per spin (order alpha, beta):  !!
!!   1 kernel: tau^s_AB = 1/2 sum_ij K^s_ij grad(phi^s_i).grad(phi^s_j);   !!
!!   2 von Weizsacker: tau^s_AB = sigma^ss_AB/(8 rho^s_AB);                 !!
!!   3 weights: tau^s_AB = 2 w_A w_B tau^s.                                 !!
!! arguments:                                                              !!
!!   itotps (in)  -- total number of grid points                          !!
!!   omp2   (in)  -- atomic weight of each point for every atom            !!
!!   chpd   (in)  -- alpha MO gradient at each grid point (from ugrdrho)   !!
!!   chpbd  (in)  -- beta MO gradient at each grid point (from ugrdrho)    !!
!!   rhoab  (in)  -- BODEN of the pair, alpha then beta                    !!
!!   sigab  (in)  -- sigma BODEN of the pair (from grdboden_uks)           !!
!!   taumol (in)  -- molecular tau, alpha then beta                        !!
!!   tauab  (out) -- tau BODEN of the pair, alpha then beta                !!
!!   iat,kat (in) -- the atom pair                                         !!
!!   sab    (in)  -- per-atom alpha MO overlap matrix                      !!
!!   sab2   (in)  -- per-atom beta MO overlap matrix                       !!
!! author: MGimf.                                                          !!
!! *********************************************************************** !!
      subroutine tauboden_uks(itotps,omp2,chpd,chpbd,rhoab,sigab,taumol,
     +  tauab,iat,kat,sab,sab2)
      use input_options_mod, only: itauboden
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chpd(itotps,nalf,3),chpbd(itotps,nb,3),omp2(itotps,nat)
      dimension rhoab(2,itotps),sigab(3,itotps),taumol(2,itotps)
      dimension tauab(2,itotps)
      dimension sab(nalf,nalf,nat),sab2(nb,nb,nat)

      select case(itauboden)
      case(1)
!! parallel over irun: each iteration only reads its own omp2(irun,:)/  !!
!! chpd(irun,:,:)/chpbd(irun,:,:) and writes its own tauab(:,irun).     !!
!$OMP PARALLEL DO PRIVATE(irun,w1,w2,xa,xb,ii,jj,gab,dd)
        do irun=1,itotps
          w1=omp2(irun,iat)
          w2=omp2(irun,kat)
          xa=ZERO
          do ii=1,nalf
            do jj=ii,nalf
              gab=sab(ii,jj,iat)*w2+sab(ii,jj,kat)*w1
              dd=chpd(irun,ii,1)*chpd(irun,jj,1)
     +          +chpd(irun,ii,2)*chpd(irun,jj,2)
     +          +chpd(irun,ii,3)*chpd(irun,jj,3)
              if(ii.ne.jj) dd=TWO*dd
              xa=xa+gab*dd
            end do
          end do
          xb=ZERO
          do ii=1,nb
            do jj=ii,nb
              gab=sab2(ii,jj,iat)*w2+sab2(ii,jj,kat)*w1
              dd=chpbd(irun,ii,1)*chpbd(irun,jj,1)
     +          +chpbd(irun,ii,2)*chpbd(irun,jj,2)
     +          +chpbd(irun,ii,3)*chpbd(irun,jj,3)
              if(ii.ne.jj) dd=TWO*dd
              xb=xb+gab*dd
            end do
          end do
          tauab(1,irun)=HALF*xa
          tauab(2,irun)=HALF*xb
        end do
!$OMP END PARALLEL DO
      case(2)
        do irun=1,itotps
          tauab(1,irun)=ZERO
          tauab(2,irun)=ZERO
!! no density cutoff: tau=0 at a point libxc still evaluates gives NaN !!
!! in its polarized TPSS correlation                                    !!
          if(rhoab(1,irun).gt.ZERO)
     +      tauab(1,irun)=sigab(1,irun)/(8.0d0*rhoab(1,irun))
          if(rhoab(2,irun).gt.ZERO)
     +      tauab(2,irun)=sigab(3,irun)/(8.0d0*rhoab(2,irun))
        end do
      case(3)
        do irun=1,itotps
          tauab(1,irun)=TWO*omp2(irun,iat)*omp2(irun,kat)*taumol(1,irun)
          tauab(2,irun)=TWO*omp2(irun,iat)*omp2(irun,kat)*taumol(2,irun)
        end do
      end select
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: tau_diagnostics                                             !!
!! purpose: EXPLORATION. For one atom pair: the share of the BODEN (in %)  !!
!!   at points where sigma_AB > 8 rho_AB tau_AB, which libxc clips to the  !!
!!   bound (Fermi-hole curvature), and the BODEN-weighted mean of          !!
!!   z = tau_W/tau as libxc sees it (min(z,1)). See tau_diagnostics_uks.   !!
!! arguments:                                                              !!
!!   itotps (in)  -- total number of grid points                          !!
!!   wp     (in)  -- integration weight of each grid point                 !!
!!   omp    (in)  -- weight of each grid point for its own atom            !!
!!   rhoab  (in)  -- BODEN of the pair                                     !!
!!   sigab  (in)  -- sigma BODEN of the pair                               !!
!!   tauab  (in)  -- tau BODEN of the pair                                 !!
!!   clip   (out) -- clipped share of the BODEN, %                         !!
!!   zmean  (out) -- BODEN-weighted mean of min(tau_W/tau,1)               !!
!! author: MGimf.                                                          !!
!! *********************************************************************** !!
      subroutine tau_diagnostics(itotps,wp,omp,rhoab,sigab,tauab,clip,
     +  zmean)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      dimension wp(itotps),omp(itotps)
      dimension rhoab(itotps),sigab(itotps),tauab(itotps)

      den=ZERO
      xc1=ZERO
      xz=ZERO
!! parallel over irun: reductions only !!
!$OMP PARALLEL DO PRIVATE(irun,ww,z) REDUCTION(+:den,xc1,xz)
      do irun=1,itotps
        if(rhoab(irun).gt.1.0d-12) then
          ww=wp(irun)*omp(irun)*rhoab(irun)
          den=den+ww
          if(sigab(irun).gt.8.0d0*rhoab(irun)*tauab(irun)*(ONE+1.0d-10)) then
            xc1=xc1+ww
            z=ONE
          else
            z=sigab(irun)/(8.0d0*rhoab(irun)*tauab(irun))
          end if
          xz=xz+ww*z
        end if
      end do
!$OMP END PARALLEL DO
      clip=ZERO
      zmean=ZERO
      if(den.gt.ZERO) then
        clip=100.0d0*xc1/den
        zmean=xz/den
      end if
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: tau_diagnostics_uks                                         !!
!! purpose: UKS twin of tau_diagnostics, both spins pooled (each point     !!
!!   and spin weighted by its BODEN).                                      !!
!! arguments: as tau_diagnostics, with rhoab/tauab (alpha, beta) and       !!
!!   sigab (up-up, up-down, down-down)                                     !!
!! author: MGimf.                                                          !!
!! *********************************************************************** !!
      subroutine tau_diagnostics_uks(itotps,wp,omp,rhoab,sigab,tauab,
     +  clip,zmean)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      dimension wp(itotps),omp(itotps)
      dimension rhoab(2,itotps),sigab(3,itotps),tauab(2,itotps)

      den=ZERO
      xc1=ZERO
      xz=ZERO
!! parallel over irun: reductions only !!
!$OMP PARALLEL DO PRIVATE(irun,is,ww,ss,z) REDUCTION(+:den,xc1,xz)
      do irun=1,itotps
        do is=1,2
          if(rhoab(is,irun).gt.1.0d-12) then
            ww=wp(irun)*omp(irun)*rhoab(is,irun)
            ss=sigab(2*is-1,irun)
            den=den+ww
            if(ss.gt.8.0d0*rhoab(is,irun)*tauab(is,irun)*(ONE+1.0d-10)) then
              xc1=xc1+ww
              z=ONE
            else
              z=ss/(8.0d0*rhoab(is,irun)*tauab(is,irun))
            end if
            xz=xz+ww*z
          end if
        end do
      end do
!$OMP END PARALLEL DO
      clip=ZERO
      zmean=ZERO
      if(den.gt.ZERO) then
        clip=100.0d0*xc1/den
        zmean=xz/den
      end if
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: print_tau_diagnostics                                       !!
!! purpose: EXPLORATION. Prints, for every BODEN pair, the clipped share   !!
!!   and the mean tau_W/tau from tau_diagnostics(_uks).                    !!
!! arguments:                                                              !!
!!   threbod (in) -- bond-order threshold of the BODEN pairs               !!
!!   clipab  (in) -- clipped share of each pair's BODEN, %                 !!
!!   zab     (in) -- BODEN-weighted mean of min(tau_W/tau,1) of each pair  !!
!! author: MGimf.                                                          !!
!! *********************************************************************** !!
      subroutine print_tau_diagnostics(threbod,clipab,zab)
      use input_options_mod, only: itauboden
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      dimension clipab(maxat,maxat),zab(maxat,maxat)
      character*40 cname(3)
      data cname /'kernel (sum K_ij grad phi_i.grad phi_j)',
     +  'von Weizsacker (sigma/8rho)','weights (2 w_A w_B tau)'/

      write(*,'(a)') ''
      write(*,'(2x,a,i0,a,a)') 'Pair tau (TAU_BODEN ',itauboden,'): ',
     +  trim(cname(itauboden))
      write(*,'(2x,a)') 'Clipped: share of the BODEN where sigma > 8 rho tau'//
     +  ' (libxc clips sigma).'
      write(*,'(2x,a)') '<z>: BODEN-weighted mean of min(tau_W/tau,1).'
      write(*,'(a)') ''
      write(*,'(2x,a)') 'Atom pair   Clipped (%)       <z>'
      write(*,'(2x,a)') repeat('-',34)
      do iatom=1,nat
        do jatom=iatom+1,nat
          if(bo(iatom,jatom).ge.threbod) write(*,'(2x,i3,a2,i3,f14.4,f10.4)')
     +      iatom,' -',jatom,clipab(iatom,jatom),zab(iatom,jatom)
        end do
      end do
      write(*,'(2x,a)') repeat('-',34)
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: print_skipped_pairs                                         !!
!! purpose: explains, under the pair counts of the DFT part, where the xc  !!
!!   term of the atom pairs below the THREBOD threshold comes from.        !!
!! arguments:                                                              !!
!!   xmix  (in) -- exact-exchange fraction of the functional               !!
!! author: MGimf                                                           !!
!! *********************************************************************** !!
      subroutine print_skipped_pairs(xmix)
      implicit real*8(a-h,o-z)
      include 'parameter.h'

      write(*,'(a)') ''
      write(*,'(2x,a)') 'Pairs below the THREBOD threshold take their whole'//
     +  ' exchange-correlation'
      write(*,'(2x,a)') 'from the multipolar expansion, not from BODEN.'
      if(xmix.gt.ZERO) then
        write(*,'(2x,a,f5.3,a,f5.3,a)') 'The DFT part holds its (1 - ',
     +    xmix,') share; the ',xmix,' share is in'
        write(*,'(2x,a)') 'EXACT-EXCHANGE SHARE.'
      end if
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: func_info_print                                             !!
!! purpose: prints one component of the DENSITY FUNCTIONAL section (libxc  !!
!!   name and id, kind, family, exact-exchange fraction, references) via   !!
!!   libxc's own info query functions, and classifies the functional into  !!
!!   itype (1 LDA, 2 GGA, 3 meta-GGA) for the caller's dispatch. Stops    !!
!!   cleanly on anything else (invalid id, range-separated, VV10...).      !!
!! arguments:                                                              !!
!!   id_func (in)  -- libxc functional id (from iopt, read by the caller)  !!
!!   itype   (out) -- functional family classification, see above          !!
!!   ifirst  (in)  -- 1 for the first of up to 3 consecutive calls (opens  !!
!!                   the section), 0 otherwise (a blank line before)       !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine func_info_print(id_func,itype,ifirst)

      use xc_f03_lib_m
      use output_mod, only: print_value
      use json_mod

      implicit real*8(a-h,o-z)

      TYPE(xc_f03_func_t) :: xc_func
      TYPE(xc_f03_func_info_t) :: xc_info
      TYPE(xc_f03_func_reference_t) :: xc_ref

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common/exchg/exch(maxat,maxat),xmix

      character*1024 name_ref
      character*128 name_func
      character*256 cbuf
      character*20 ckind,cfam,clab

      if(kop.ne.1) call xc_f03_func_init(xc_func,id_func,XC_UNPOLARIZED,ierr)
      if(kop.eq.1) call xc_f03_func_init(xc_func,id_func,XC_POLARIZED,ierr)
!! an unknown id leaves xc_func without info -- querying it segfaults. !!
      if(ierr.ne.0) then
        write(*,'(a)') ''
        write(*,'(2x,a,1x,i0,1x,a)') 'ERROR: functional id',id_func,
     +    'is not a valid libxc functional id.'
        write(*,'(2x,a)') 'See https://libxc.gitlab.io/functionals/ for the list of ids.'
        call apost_stop('INVALID LIBXC FUNCTIONAL ID. REVISE inp')
      end if
      xc_info = xc_f03_func_get_info(xc_func)
      xmix = xc_f03_hyb_exx_coef(xc_func)
      name_func = xc_f03_func_info_get_name(xc_info)

      select case(xc_f03_func_info_get_kind(xc_info))
      case(XC_EXCHANGE)
        ckind='exchange'
      case(XC_CORRELATION)
        ckind='correlation'
      case(XC_EXCHANGE_CORRELATION)
        ckind='exchange-correlation'
      case(XC_KINETIC)
        ckind='kinetic energy'
      case default
        ckind='unknown'
      end select

!! itype: 1 LDA, 2 GGA, 3 meta-GGA. Hybrid functionals get xmix>0 as well. !!
!! Caller stores the returned value in iopt(55) for later use. !!
      iflags=xc_f03_func_info_get_flags(xc_info)
      irsh=iand(iflags,XC_FLAGS_HYB_CAM+XC_FLAGS_HYB_CAMY+XC_FLAGS_HYB_LC
     +  +XC_FLAGS_HYB_LCY)
      itype=0
      select case(xc_f03_func_info_get_family(xc_info))
      case(XC_FAMILY_LDA)
        cfam='LDA'
        itype=1
      case(XC_FAMILY_GGA)
        cfam='GGA'
        itype=2
      case(XC_FAMILY_HYB_GGA)
        cfam='hybrid GGA'
        itype=2
      case(XC_FAMILY_MGGA)
        cfam='meta-GGA'
        itype=3
      case(XC_FAMILY_HYB_MGGA)
        cfam='hybrid meta-GGA'
        itype=3
      case default
        cfam='other'
        itype=-1
      end select

      if(ifirst.eq.1) then
        call print_box('DENSITY FUNCTIONAL')
      else
        write(*,'(a)') ''
      end if
      write(cbuf,'(a,a,i0,a)') trim(name_func),' (libxc id ',id_func,')'
      call print_value('libxc functional',trim(cbuf))
      call print_value('Kind',trim(ckind))
      call print_value('Family',trim(cfam))
!! for a range-separated hybrid this would be the long-range fraction only. !!
      if(itype.ge.2.and.xmix.gt.ZERO.and.irsh.eq.0)
     +  call print_value('Exact-exchange fraction',xmix,3)

!! xc_f03_func_info_get_references uses lookahead semantics: the call that !!
!! returns the LAST reference also sets ii=-1 on that same call. irefref   !!
!! numbers the printed references, ii is the next-index cursor.            !!
      ii=0
      irefref=0
      do while(ii.ge.0)
        xc_ref = xc_f03_func_info_get_references(xc_info,ii)
        irefref=irefref+1
        name_ref = xc_f03_func_reference_get_ref(xc_ref)
        write(clab,'(a,i0,a)') 'Reference [',irefref,']'
        call print_value(trim(clab),trim(adjustl(name_ref)))
      end do

!! results file: one object per libxc component (main.f opened the array) !!
      call json_begin_object(' ')
      call json_value('libxc_id',id_func)
      call json_value('name',trim(name_func))
      call json_value('kind',trim(ckind))
      call json_value('family',trim(cfam))
      call json_end()

!! everything past this point is a clean stop for functionals APOST-3D !!
!! can't decompose yet -- before any energy is computed, never silently. !!
      ikind=xc_f03_func_info_get_kind(xc_info)
      call xc_f03_func_end(xc_func)

      if(itype.lt.0) then
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Only LDA, GGA, meta-GGA and their global hybrids are supported.'
        call apost_stop('UNSUPPORTED FUNCTIONAL FAMILY. REVISE inp')
      end if
      if(ikind.eq.XC_KINETIC) then
        write(*,'(a)') ''
        write(*,'(2x,a)') 'This is a kinetic-energy functional, not an exchange or'
        write(*,'(2x,a)') 'correlation one.'
        call apost_stop('KINETIC-ENERGY FUNCTIONAL GIVEN. REVISE inp')
      end if
!! range-separated hybrids: hyb_exx_coef is only the long-range/full-     !!
!! range fraction, the short-range exact exchange is never integrated;    !!
!! VV10: libxc's exc omits the nonlocal correlation part.                 !!
      if(irsh.ne.0) then
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Range-separated hybrids (CAM-B3LYP, wB97X, HSE, LC-...) are'
        write(*,'(2x,a)') 'not supported: only global hybrids are.'
        call apost_stop('RANGE-SEPARATED HYBRID NOT SUPPORTED. REVISE inp')
      end if
      if(iand(iflags,XC_FLAGS_VV10).ne.0) then
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Functionals with VV10 nonlocal correlation are not supported.'
        call apost_stop('VV10 FUNCTIONAL NOT SUPPORTED. REVISE inp')
      end if

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
