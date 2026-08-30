!! ********************************************************************* !!
!! subroutine: dft_dm1                                                   !!
!! purpose: DFT-DM1 (`# METHOD`/`DFT-DM1`, formerly referred to           !!
!!   internally as HIRAO) -- builds a Hirao-style approximate one-       !!
!!   particle RDM1 for UHF/UKS-DFT from the local exchange-energy        !!
!!   density, then reports: the RDM1's exchange energy and bond-order    !!
!!   (delocalization-index) matrices, the exact HF-type exchange from    !!
!!   the real KS orbitals for comparison, and (opt-in, `# DFT-DM1`'s     !!
!!   `NATORB`) the RDM1 projected onto the AO basis and diagonalized for !!
!!   natural-orbital occupations. Builds its own rotated second grid      !!
!!   internally (see `# DFT-DM1`'s `MOD-GRIDTWOEL`/`# GRID`); the first   !!
!!   grid (`wp`/`omp2`/`pcoord`/`chp`) is the caller's (`main.f`).        !!
!! arguments:                                                            !!
!!   itotps (in) -- total grid points on the first grid                  !!
!!   wp, omp2, pcoord, chp (in) -- first-grid quadrature weight, fuzzy-   !!
!!     atom partition weight, coordinates, and AO values                 !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine dft_dm1(itotps,wp,omp2,pcoord,chp)
      use ao_matrices
      use integration_grid
      use basis_set, only: s
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /filename/name0
      common /iops/iopt(200)
      common /dm1opt/densthresh_dm1
      character*60 name0
      integer*8 :: npairtot,npairskip

      dimension :: wp(itotps),omp2(itotps,nat),pcoord(itotps,3),chp(itotps,igr)
      dimension :: exch_hf(maxat,maxat),bondorder(maxat,maxat),exact_exch(maxat,maxat)

!! automatic (stack) per-pair scratch for the double loops below --      !!
!! OMP-PRIVATE needs this, not allocatable (see enpart.f's rvect).       !!
      dimension :: eval_ao(igr),gx_ao(igr),gy_ao(igr),gz_ao(igr)
      dimension :: chp2v(nalf),chp2bv(nb)
      dimension :: rhoab(2,1),scrpt(3,1),excpt(1),excbpt(1)

      allocatable :: wppha(:),omp2pha(:,:),pcoordpha(:,:),chppha(:,:),omp(:),ibaspoint(:)
      allocatable :: chp2(:,:),chp2pha(:,:),chp2b(:,:),chp2phab(:,:)
      allocatable :: rhoscr(:)

!! NATORB block only (RDM1 projected onto the AO basis, alpha/beta       !!
!! separately) -- kept out of the always-used arrays above since this    !!
!! whole block is opt-in and unrelated to anything else in the file.     !!
      allocatable :: dm1_ao(:,:),dm1b_ao(:,:),s0no(:,:),smno(:,:),spno(:,:)
      allocatable :: cno_a(:,:)

      ifunc  = Iopt(57)
      inatorb= Iopt(63)
      iatps  = nrad*nang

      if(kop.ne.1) stop " Only implemented for unrestricted WF "

!! GENERAL INFORMATION !!

      call print_subbox('GENERAL INFORMATION')
      write(*,'(2x,a,1x,i0)')   'Number of basis functions           :',igr
      write(*,'(2x,a,1x,i0)')   'Number of occupied alpha MOs        :',nalf
      write(*,'(2x,a,1x,i0)')   'Number of occupied beta MOs         :',nb
      write(*,'(2x,a,1x,i0)')   'Number of radial points             :',nrad
      write(*,'(2x,a,1x,i0)')   'Number of angular points per radial :',nang
      write(*,'(2x,a,1x,f6.3)') 'Gauss-Legendre R0 parameter         :',rr00
      call func_info_print(ifunc,itype,0)

!! GENERATING SECOND GRID (ROTATED) FOR NUMERICAL INTEGRATION !!
!! chp: value of AO jj at point ii (first grid); chppha: same, second grid !!

      ndim  = igr
      iatps = nang*nrad
      pha   = ZERO
      phb   = 0.162d0
      call quad(nrad,nang)
      ALLOCATE(wppha(itotps),omp(itotps),omp2pha(itotps,nat))
      ALLOCATE(pcoordpha(itotps,3),ibaspoint(itotps),chppha(itotps,igr))
      ALLOCATE(rhoscr(itotps))
      call prenumint(ndim,itotps,nat,wppha,omp,omp2pha,chppha,rhoscr,pcoordpha,ibaspoint,0)
      DEALLOCATE(ibaspoint,omp,rhoscr)

!! TRANSFORMATION TO MOs, BOTH GRIDS -- needed later for the exact       !!
!! HF-type exchange comparison against the real KS orbitals.             !!

      ALLOCATE(chp2(itotps,nalf),chp2pha(itotps,nalf))
      ALLOCATE(chp2b(itotps,nb),chp2phab(itotps,nb))
      do kk=1,itotps
        do imo=1,nalf
          xx=ZERO
          xxb=ZERO
          xxpha=ZERO
          xxphab=ZERO
          do ibf=1,igr
            xx=xx+c(ibf,imo)*chp(kk,ibf)
            xxpha=xxpha+c(ibf,imo)*chppha(kk,ibf)
            if(imo.le.nb) then
              xxb=xxb+cb(ibf,imo)*chp(kk,ibf)
              xxphab=xxphab+cb(ibf,imo)*chppha(kk,ibf)
            end if
          end do
          chp2(kk,imo)=xx
          chp2pha(kk,imo)=xxpha
          if(imo.le.nb) then
            chp2b(kk,imo)=xxb
            chp2phab(kk,imo)=xxphab
          end if
        end do
      end do

!! GENERATING KS-DFT RDM1, ONE PAIR OF GRID POINTS AT A TIME !!

      call print_box('GENERATING KS-DFT RDM1')
        xexch=ZERO
        xexchb=ZERO
        npairtot=0
        npairskip=0
        do icenter=1,nat
          do jcenter=1,nat
            f3=ZERO
            f3b=ZERO
            fbo=ZERO
            fbob=ZERO

!! parallel over ifut: PRIVATE scratch is all automatic (no per-thread   !!
!! allocation needed); gpoints/drho_xyz/sigma_uks_xyz's "current basis   !!
!! function" state is THREADPRIVATE'd at its own declaration; c/cb are   !!
!! shared but read-only; everything else accumulates via REDUCTION.      !!
!$OMP PARALLEL DO PRIVATE(jfut,Rx,Ry,Rz,rhoa,rhob,rhoab,imo,ibf,xx,xxb,
!$OMP&  scraa,scrab,scrbb,scrpt,xfact,excpt,excbpt,r12,x1,xx1,xx1b,
!$OMP&  xksigaa,xksigbb,xx0,xx0b,xkagga,xkbgga,xbf,xbfb,xxrdm1,xxrdm1b,
!$OMP&  x0,eval_ao,gx_ao,gy_ao,gz_ao,chp2v,chp2bv)
!$OMP&  REDUCTION(+:f3,f3b,fbo,fbob,xexch,xexchb,npairtot,npairskip)
            do ifut=iatps*(icenter-1)+1,iatps*icenter
              x0=wp(ifut)*omp2(ifut,icenter)
              do jfut=iatps*(jcenter-1)+1,iatps*jcenter
                Rx=(pcoord(ifut,1)+pcoordpha(jfut,1))/TWO
                Ry=(pcoord(ifut,2)+pcoordpha(jfut,2))/TWO
                Rz=(pcoord(ifut,3)+pcoordpha(jfut,3))/TWO
!! gpoints (qtaim.f) already evaluates AOs via basis_set's coefpb, which  !!
!! bakes the pure-d/f transform into the contraction coefficients -- no  !!
!! separate 5d/6d reorder step needed (see sigma_uks's own AO-gradient   !!
!! loop for the same pattern), so the old gordermat call is dropped.     !!
                call gpoints(Rx,Ry,Rz,gx_ao,gy_ao,gz_ao,eval_ao)
                call calc_uhf_dens(eval_ao,rhoa,rhob)
                npairtot=npairtot+1

!! prune on the density AT THE MIDPOINT R -- not at ifut/jfut themselves !!
!! (a point being in a low-density tail on its own grid doesn't mean R,  !!
!! the actual point the RDM1 kernel below is evaluated at, is negligible !!
!! too). Controlled by # DFT-DM1's DENSTHRESH (default 1e-8). Skips the  !!
!! rest of this pair's cost (sigma_uks_xyz, xc_uks_for_dm1, the Bessel-  !!
!! kernel exchange accumulation) but not gpoints/calc_uhf_dens itself,   !!
!! since R's density isn't known until after that call.                 !!
                if(abs(rhoa+rhob).lt.densthresh_dm1) then
                  npairskip=npairskip+1
                  cycle
                end if

                rhoab(1,1)=rhoa
                rhoab(2,1)=rhob

!! COMPUTING MOs FOR SIGMA CALCULATION !!

                if(itype.gt.1) then
                  do imo=1,nalf
                    xx=ZERO
                    xxb=ZERO
                    do ibf=1,igr
                      xx=xx+c(ibf,imo)*eval_ao(ibf)
                      if(imo.le.nb) xxb=xxb+cb(ibf,imo)*eval_ao(ibf)
                    end do
                    chp2v(imo)=xx
                    if(imo.le.nb) chp2bv(imo)=xxb
                  end do
                  call sigma_uks_xyz(Rx,Ry,Rz,chp2v,chp2bv,scraa,scrab,scrbb)
                  scrpt(1,1)=scraa
                  scrpt(2,1)=scrab
                  scrpt(3,1)=scrbb
                end if

!! COMPUTING BOTH RDM1 AND EXCHANGE HERE !!

                call xc_uks_for_dm1(1,1,ifunc,rhoab,scrpt,excpt)
                call xc_uks_for_dm1(2,1,ifunc,rhoab,scrpt,excbpt)
                xfact=FOUR/THREE
                r12=(pcoord(ifut,1)-pcoordpha(jfut,1))**TWO
                r12=r12+((pcoord(ifut,2)-pcoordpha(jfut,2))**TWO)
                r12=r12+((pcoord(ifut,3)-pcoordpha(jfut,3))**TWO)
                r12=dsqrt(r12)
                x1=wppha(jfut)*omp2pha(jfut,jcenter)

!! COMPUTING KsGGA (1 = ALPHA, 2 = BETA), GENERAL FOR ALL FUNCTIONALS FROM THE LIBRARY !!

                xx1=(rhoa**xfact)
                xx1b=(rhob**xfact)
                if(ABS(xx1).gt.1.0d-14) then
                  xksigaa=-TWO*excpt(1)/xx1
                  xx0=(9.0d0*pi/xksigaa)**HALF
                else
                  xksigaa=ZERO
                  xx0=ZERO
                end if
                if(ABS(xx1b).gt.1.0d-14) then
                  xksigbb=-TWO*excbpt(1)/xx1b
                  xx0b=(9.0d0*pi/xksigbb)**HALF
                else
                  xksigbb=ZERO
                  xx0b=ZERO
                end if
                xkagga=xx0*(rhoa**(ONE/THREE))
                xkbgga=xx0b*(rhob**(ONE/THREE))
                xx0=xkagga*r12
                xx0b=xkbgga*r12

!! xbf = J_1(X)/X !!
!! FIRST ALPHA !!

                if(ABS(xx0).lt.thresh) then
                  xbf=ONE/THREE
                else
                  xbf=(dsin(xx0)-xx0*dcos(xx0))/(xx0**THREE)
                end if

!! NOW BETA !!
                if(ABS(xx0b).lt.thresh) then
                  xbfb=ONE/THREE
                else
                  xbfb=(dsin(xx0b)-xx0b*dcos(xx0b))/(xx0b**THREE)
                end if

!! RDM1 CONTAINING BOTH ALPHA AND BETA !!

                xxrdm1=THREE*xbf*rhoa
                xxrdm1b=THREE*xbfb*rhob

!! COMPUTING ONLY ONCE FROM RDM1 !!

                if(r12.gt.thresh) then
                  xexch=xexch-(xxrdm1*xxrdm1*x0*x1/r12)
                  xexchb=xexchb-(xxrdm1b*xxrdm1b*x0*x1/r12)
                  f3=f3+(xxrdm1*xxrdm1*x0*x1/r12)
                  f3b=f3b+(xxrdm1b*xxrdm1b*x0*x1/r12)
                end if

!! BOND ORDER (DELOCALIZATION INDEX) FROM THE SAME RDM1 -- same double   !!
!! integral as the exchange energy above but without the 1/r12 weight,   !!
!! so no r12>thresh guard is needed here (there's no singularity to      !!
!! avoid once r12 isn't a denominator).                                  !!
                fbo=fbo+(xxrdm1*xxrdm1*x0*x1)
                fbob=fbob+(xxrdm1b*xxrdm1b*x0*x1)
              end do
            end do
!$OMP END PARALLEL DO
            if(icenter.ne.jcenter) then
              f3=TWO*f3
              f3b=TWO*f3b
              fbo=TWO*fbo
              fbob=TWO*fbob
            end if
            exch_hf(icenter,jcenter)=-(f3+f3b)/TWO
            bondorder(icenter,jcenter)=fbo+fbob
          end do
        end do

        write(*,'(2x,a,1x,i14)') "Grid-point pairs evaluated:",npairtot
        write(*,'(2x,a,1x,i14,1x,a,1x,f6.2,1x,a)') "Pairs pruned (R below DENSTHRESH):",
     $npairskip,"(",100.0d0*dble(npairskip)/dble(npairtot),"%)"
        write(*,'(2x,a,1x,i14,1x,a,1x,f6.2,1x,a)') "Pairs kept (full cost paid):",
     $npairtot-npairskip,"(",100.0d0*dble(npairtot-npairskip)/dble(npairtot),"%)"

      call flush

      call print_box('DFT-DM1 EXCHANGE ENERGY (RDM1)')
      call MPRINT2(exch_hf,nat,maxat)
      write(*,'(2x,a,1x,f14.7)') 'Exchange energy (RDM1) :',(xexch+xexchb)/TWO

      call print_box('DFT-DM1 BOND ORDER (RDM1)')
      call MPRINT2(bondorder,nat,maxat)

!! EXACT HF-TYPE EXCHANGE FROM THE REAL KS ORBITALS, FOR COMPARISON --    !!
!! same real-space same-spin exchange integral as enpart.f's numint_two, !!
!! reusing dft_dm1's own grids. Cheap: chp2/chp2pha/chp2b/chp2phab (real !!
!! MO values, both grids) are already computed, nothing here is approximated.!!
      call print_box('EXACT HF-TYPE EXCHANGE ENERGY (KS ORBITALS)')
      xexact=ZERO
      xexactb=ZERO
      do icenter=1,nat
        do jcenter=1,nat
          fex=ZERO
          fexb=ZERO
!$OMP PARALLEL DO PRIVATE(jfut,x1,r12,i,pab,pabb)
!$OMP&  REDUCTION(+:fex,fexb,xexact,xexactb)
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            x0=wp(ifut)*omp2(ifut,icenter)
            do jfut=iatps*(jcenter-1)+1,iatps*jcenter
              x1=wppha(jfut)*omp2pha(jfut,jcenter)
              r12=(pcoord(ifut,1)-pcoordpha(jfut,1))**TWO
              r12=r12+((pcoord(ifut,2)-pcoordpha(jfut,2))**TWO)
              r12=r12+((pcoord(ifut,3)-pcoordpha(jfut,3))**TWO)
              r12=dsqrt(r12)
              if(r12.gt.thresh) then
                pab=ZERO
                do i=1,nalf
                  pab=pab+chp2(ifut,i)*chp2pha(jfut,i)
                end do
                pabb=ZERO
                do i=1,nb
                  pabb=pabb+chp2b(ifut,i)*chp2phab(jfut,i)
                end do
                xexact=xexact-(pab*pab*x0*x1/r12)
                xexactb=xexactb-(pabb*pabb*x0*x1/r12)
                fex=fex+(pab*pab*x0*x1/r12)
                fexb=fexb+(pabb*pabb*x0*x1/r12)
              end if
            end do
          end do
!$OMP END PARALLEL DO
          if(icenter.ne.jcenter) then
            fex=TWO*fex
            fexb=TWO*fexb
          end if
          exact_exch(icenter,jcenter)=-(fex+fexb)/TWO
        end do
      end do
      call MPRINT2(exact_exch,nat,maxat)
      write(*,'(2x,a,1x,f14.7)') 'Exchange energy (exact) :',(xexact+xexactb)/TWO

!! RDM1 -> AO BASIS -> DIAGONALIZE -> NATURAL ORBITAL OCCUPATIONS !!

!! opt-in (# DFT-DM1's NATORB) -- redoes the main loop's per-pair work,  !!
!! then additionally accumulates an igr x igr matrix element per         !!
!! surviving pair, so it's costly even parallelized; kept separate and   !!
!! opt-in rather than folded into the main loop above.                   !!
      if(inatorb.eq.1) then
        call print_box('PROJECTING THE RDM1 ONTO THE AO BASIS')

        ALLOCATE(dm1_ao(igr,igr),dm1b_ao(igr,igr))
        dm1_ao=ZERO
        dm1b_ao=ZERO

        do icenter=1,nat
          do jcenter=1,nat
!! same PRIVATE/THREADPRIVATE reasoning as the main double loop above;   !!
!! dm1_ao/dm1b_ao add an array REDUCTION (gfortran/OpenMP 4.5+ supports  !!
!! this for already-allocated arrays, which they are before this loop). !!
!$OMP PARALLEL DO PRIVATE(jfut,Rx,Ry,Rz,rhoa,rhob,rhoab,imo,ibf,xx,xxb,
!$OMP&  scraa,scrab,scrbb,scrpt,xfact,excpt,excbpt,r12,x1,xx1,xx1b,
!$OMP&  xksigaa,xksigbb,xx0,xx0b,xkagga,xkbgga,xbf,xbfb,xxrdm1,xxrdm1b,
!$OMP&  x0,eval_ao,gx_ao,gy_ao,gz_ao,chp2v,chp2bv,mu,nu)
!$OMP&  REDUCTION(+:dm1_ao,dm1b_ao)
            do ifut=iatps*(icenter-1)+1,iatps*icenter
              x0=wp(ifut)*omp2(ifut,icenter)
              do jfut=iatps*(jcenter-1)+1,iatps*jcenter
                Rx=(pcoord(ifut,1)+pcoordpha(jfut,1))/TWO
                Ry=(pcoord(ifut,2)+pcoordpha(jfut,2))/TWO
                Rz=(pcoord(ifut,3)+pcoordpha(jfut,3))/TWO
                call gpoints(Rx,Ry,Rz,gx_ao,gy_ao,gz_ao,eval_ao)
                call calc_uhf_dens(eval_ao,rhoa,rhob)
                if(abs(rhoa+rhob).lt.densthresh_dm1) cycle

                rhoab(1,1)=rhoa
                rhoab(2,1)=rhob
                if(itype.gt.1) then
                  do imo=1,nalf
                    xx=ZERO
                    xxb=ZERO
                    do ibf=1,igr
                      xx=xx+c(ibf,imo)*eval_ao(ibf)
                      if(imo.le.nb) xxb=xxb+cb(ibf,imo)*eval_ao(ibf)
                    end do
                    chp2v(imo)=xx
                    if(imo.le.nb) chp2bv(imo)=xxb
                  end do
                  call sigma_uks_xyz(Rx,Ry,Rz,chp2v,chp2bv,scraa,scrab,scrbb)
                  scrpt(1,1)=scraa
                  scrpt(2,1)=scrab
                  scrpt(3,1)=scrbb
                end if

                call xc_uks_for_dm1(1,1,ifunc,rhoab,scrpt,excpt)
                call xc_uks_for_dm1(2,1,ifunc,rhoab,scrpt,excbpt)
                xfact=FOUR/THREE
                xx1=(rhoa**xfact)
                xx1b=(rhob**xfact)
                if(ABS(xx1).gt.1.0d-14) then
                  xksigaa=-TWO*excpt(1)/xx1
                  xx0=(9.0d0*pi/xksigaa)**HALF
                else
                  xksigaa=ZERO
                  xx0=ZERO
                end if
                if(ABS(xx1b).gt.1.0d-14) then
                  xksigbb=-TWO*excbpt(1)/xx1b
                  xx0b=(9.0d0*pi/xksigbb)**HALF
                else
                  xksigbb=ZERO
                  xx0b=ZERO
                end if
                xkagga=xx0*(rhoa**(ONE/THREE))
                xkbgga=xx0b*(rhob**(ONE/THREE))
                r12=(pcoord(ifut,1)-pcoordpha(jfut,1))**TWO
                r12=r12+((pcoord(ifut,2)-pcoordpha(jfut,2))**TWO)
                r12=r12+((pcoord(ifut,3)-pcoordpha(jfut,3))**TWO)
                r12=dsqrt(r12)
                xx0=xkagga*r12
                xx0b=xkbgga*r12
                if(ABS(xx0).lt.thresh) then
                  xbf=ONE/THREE
                else
                  xbf=(dsin(xx0)-xx0*dcos(xx0))/(xx0**THREE)
                end if
                if(ABS(xx0b).lt.thresh) then
                  xbfb=ONE/THREE
                else
                  xbfb=(dsin(xx0b)-xx0b*dcos(xx0b))/(xx0b**THREE)
                end if
                xxrdm1=THREE*xbf*rhoa
                xxrdm1b=THREE*xbfb*rhob

                x1=wppha(jfut)*omp2pha(jfut,jcenter)
                do mu=1,igr
                  do nu=1,igr
                    dm1_ao(mu,nu)=dm1_ao(mu,nu)+chp(ifut,mu)*xxrdm1*chppha(jfut,nu)*x0*x1
                    dm1b_ao(mu,nu)=dm1b_ao(mu,nu)+chp(ifut,mu)*xxrdm1b*chppha(jfut,nu)*x0*x1
                  end do
                end do
              end do
            end do
!$OMP END PARALLEL DO
          end do
        end do

!! symmetrize -- gamma(r1,r2)=gamma(r2,r1) exactly for a real single-     !!
!! determinant RDM1, but the two independently-rotated grids don't       !!
!! enforce mu/nu symmetry numerically pair by pair.                      !!
        do mu=1,igr
          do nu=mu+1,igr
            xx=(dm1_ao(mu,nu)+dm1_ao(nu,mu))/TWO
            dm1_ao(mu,nu)=xx
            dm1_ao(nu,mu)=xx
            xx=(dm1b_ao(mu,nu)+dm1b_ao(nu,mu))/TWO
            dm1b_ao(mu,nu)=xx
            dm1b_ao(nu,mu)=xx
          end do
        end do

!! total RDM1 = alpha + beta, matching gennatural's convention (util.f): !!
!! APOST-3D natural orbitals always come from one combined density, not !!
!! separate alpha/beta sets, even for open-shell/unrestricted cases.     !!
        do mu=1,igr
          do nu=1,igr
            dm1_ao(mu,nu)=dm1_ao(mu,nu)+dm1b_ao(mu,nu)
          end do
        end do

!! dm1_ao (double AO integrals, chp*xxrdm1*chppha) is D = S*P*S relative !!
!! to the coefficient-built P (trace(PS)=N) that gennatural's S^1/2      !!
!! transform expects -- needs S^-1/2 here, not S^1/2, or every           !!
!! occupation inflates (found via H2's alpha trace: ~16 instead of ~1).  !!
        ALLOCATE(s0no(igr,igr),smno(igr,igr),spno(igr,igr))
        ALLOCATE(cno_a(igr,igr))
        s0no=s
        call build_Smp(igr,s0no,smno,spno,0)

        call to_lowdin_basis(igr,smno,dm1_ao)
        call diagonalize(igr,igr,dm1_ao,cno_a,0)
        call to_AO_basis(igr,igr,smno,cno_a)

!! this approximate RDM1 isn't guaranteed positive-semidefinite, so a    !!
!! negative tail can appear -- printed in full, not cut off, to show it. !!
        call print_subbox('DFT-DM1 NATURAL ORBITALS')
        write(*,'(2x,a)') 'Occupation numbers (all, including negative-tail artifacts):'
        write(*,'(2x,8f10.5)') (dm1_ao(ii,ii),ii=1,igr)

        DEALLOCATE(dm1_ao,dm1b_ao,s0no,smno,spno,cno_a)
      end if

      end

! *****

!! ********************************************************************* !!
!! subroutine: xc_uks_for_dm1                                            !!
!! purpose: evaluates one spin channel's LDA/GGA exchange energy density !!
!!   (per electron, pre-multiplied by that spin's density) via libxc,    !!
!!   for an arbitrary set of npt points.                                 !!
!! arguments:                                                            !!
!!   isigma (in) -- 1=alpha, 2=beta                                      !!
!!   npt (in) -- number of points                                        !!
!!   id_xfunc (in) -- libxc functional ID                                !!
!!   scr_ab (in) -- density, (alpha,beta)                                !!
!!   scr (in) -- sigma (grad-rho.grad-rho), (up-up,up-down,down-down)    !!
!!   scr2 (out) -- exchange energy density for spin isigma               !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine xc_uks_for_dm1(isigma,npt,id_xfunc,scr_ab,scr,scr2)
      use xc_f90_types_m
      use xc_f90_lib_m
      implicit real*8(a-h,o-z)
      TYPE(xc_f90_pointer_t) :: xc_func
      TYPE(xc_f90_pointer_t) :: xc_info
      include 'parameter.h'

      dimension :: scr_ab(2,npt),scr(3,npt),scr2(npt)

      call xc_f90_func_init(xc_func,xc_info,id_xfunc,XC_POLARIZED)
      select case (xc_f90_info_family(xc_info))
        case(XC_FAMILY_LDA)
          call xc_f90_lda_exc(xc_func,npt,scr_ab(1,1),scr2(1))
        case(XC_FAMILY_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2(1))
        case(XC_FAMILY_HYB_GGA)
          call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2(1))
        case(XC_FAMILY_MGGA)
!          call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2(1))
        case(XC_FAMILY_HYB_MGGA)
!          call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2(1))
      end select
      call xc_f90_func_end(xc_func)

!! MULTIPLYING ONLY ALPHA RHO !!

      do ii=1,npt
        scr2(ii)=scr2(ii)*scr_ab(isigma,ii)
      end do
      end

! *****

!! calc_rhf_dens/calc_uhf_dens/calc_phf_dens: density at a point from    !!
!! its AO values (eval_ao), for RHF/UHF/post-HF wavefunctions. Moved     !!
!! here from tools.f -- feature-specific helpers live with their driver, !!
!! matching enpart_dft.f/oslo.f. Only calc_uhf_dens is called today      !!
!! (dft_dm1 is UHF-only); the other two are kept for a possible future   !!
!! restricted/CASSCF DM1 variant.                                        !!

!! ********************************************************************* !!
!! subroutine: calc_rhf_dens                                             !!
!! purpose: RHF density at a point, from its AO values.                  !!
!! arguments: eval_ao (in), rho (out)                                    !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine calc_rhf_dens(eval_ao,rho)
      use ao_matrices
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop

      dimension :: eval_ao(igr)

      xx0=ZERO
      do imo=1,nocc
        xx=ZERO
        do ibf=1,igr
          xx=xx+c(ibf,imo)*eval_ao(ibf)
        end do
        xx0=xx0+xx*xx
      end do
      rho=TWO*xx0

      end

! *****

!! ********************************************************************* !!
!! subroutine: calc_uhf_dens                                             !!
!! purpose: UHF alpha/beta density at a point, from its AO values.       !!
!! arguments: eval_ao (in), rhoa, rhob (out)                             !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine calc_uhf_dens(eval_ao,rhoa,rhob)
      use ao_matrices
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop

      dimension :: eval_ao(igr)

      xx0=ZERO
      xx0b=ZERO
      do imo=1,nalf
        xx=ZERO
        xxb=ZERO
        do ibf=1,igr
          xx=xx+c(ibf,imo)*eval_ao(ibf)
          if(imo.le.nb) xxb=xxb+cb(ibf,imo)*eval_ao(ibf)
        end do
        if(imo.le.nb) xx0b=xx0b+xxb*xxb
        xx0=xx0+xx*xx
      end do
      rhoa=xx0
      rhob=xx0b

      end

! *****

!! ********************************************************************* !!
!! subroutine: calc_phf_dens                                             !!
!! purpose: post-HF density at a point (natural orbitals/occupations),   !!
!!   from its AO values.                                                 !!
!! arguments: eval_ao (in), rho (out)                                    !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine calc_phf_dens(eval_ao,rho)
      use ao_matrices
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/  nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/  icas,ncasel,ncasorb,nspinorb,norb,icisd,icass

      dimension :: eval_ao(igr)

      xx=ZERO
      do ino=1,norb
        xx1=ZERO
        do ibf=1,igr
          xx1=xx1+c_no(ibf,ino)*eval_ao(ibf)
        end do
        xx=xx+xx1*xx1*occ_no(ino,ino)
      end do
      rho=xx

      end

! *****

!! ********************************************************************* !!
!! subroutine: sigma_uks_xyz                                             !!
!! purpose: UKS sigma (grad-rho.grad-rho contractions) at an arbitrary   !!
!!   point -- distinct from enpart_dft.f's sigma_uks, which only ever    !!
!!   runs on a precomputed grid, not the (r1+r2)/2 midpoints dft_dm1's   !!
!!   double loop needs.                                                  !!
!! arguments:                                                            !!
!!   xabs, yabs, zabs (in) -- point to evaluate at                       !!
!!   chp2, chp3 (in) -- alpha/beta MO values at that point               !!
!!   scraa, scrab, scrbb (out) -- sigma (up-up, up-down, down-down)      !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine sigma_uks_xyz(xabs,yabs,zabs,chp2,chp3,scraa,scrab,scrbb)
      use ao_matrices
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /actual_bf/ iact,jat,icenter
!! DFT-DM1's double loop (dft_dm1.f) calls this under OMP -- each thread  !!
!! needs its own iact, not one shared across all of them.                !!
!$OMP THREADPRIVATE(/actual_bf/)
      common /nat/    nat,igr,ifg,nocc,nalf,nb,kop

      dimension :: chp2(nalf),chp3(nb)

!! automatic (stack), not allocatable -- called concurrently from an OMP !!
!! loop; an allocatable here would need re-ALLOCATEing per thread, an    !!
!! automatic array just works (same fix as dft_dm1's own scratch).       !!
      dimension :: chpd(igr),chp(igr),chpbd(igr)

!! GENERATING GRID FOR 2nd DERIVATIVE OVER AOs !!

      scraa=ZERO
      scrab=ZERO
      scrbb=ZERO
      do ixyz=1,3
        do ii=1,igr
          iact=ii
          chp(iact)=drho_xyz(xabs,yabs,zabs,ixyz)
        end do

!! TO MOs !!

        do j=1,nalf
          xx=ZERO
          xxb=ZERO
          do i=1,igr
            xx=xx+c(i,j)*chp(i)
            if(j.le.nb) xxb=xxb+cb(i,j)*chp(i)
          end do
          chpd(j)=xx
          if(j.le.nb) chpbd(j)=xxb
        end do

!! CALCULATION OF SIGMA IN A GIVEN POINT !!

        do i=1,nalf
          do j=1,nalf
            scraa=scraa+chp2(i)*chp2(j)*chpd(i)*chpd(j)
            if(j.le.nb) scrab=scrab+chp2(i)*chp3(j)*chpd(i)*chpbd(j)
            if(i.le.nb.and.j.le.nb) scrbb=scrbb+chp3(i)*chp3(j)*chpbd(i)*chpbd(j)
          end do
        end do
      end do
      scraa=FOUR*scraa
      scrab=FOUR*scrab
      scrbb=FOUR*scrbb

      end

! *****

!! ********************************************************************* !!
!! function: drho_xyz                                                    !!
!! purpose: AO gradient (ixyz component) at an arbitrary point, for one  !!
!!   basis function (iact, via common/actual_bf/ -- same convention      !!
!!   qtaim.f's gpoints/gxfunct use). Same primitive-loop math as         !!
!!   enpart_dft.f's sigma_uks, for a single point instead of a grid.     !!
!! arguments:                                                            !!
!!   xabs, yabs, zabs (in) -- point to evaluate at                       !!
!!   ixyz (in) -- gradient component, 1=x, 2=y, 3=z                      !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      function drho_xyz(xabs,yabs,zabs,ixyz)
      use basis_set
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /actual_bf/ iact,jat,icenter
!! called (via gxfunct/gyfunct/gzfunct and sigma_uks_xyz) from DFT-DM1's !!
!! OMP-parallelized loop -- each thread needs its own iact.              !!
!$OMP THREADPRIVATE(/actual_bf/)

      iactat=ihold(iact)
      fx=ZERO
      x=xabs-coord(1,iactat)
      y=yabs-coord(2,iactat)
      z=zabs-coord(3,iactat)
      rr=dsqrt(x**TWO+y**TWO+z**TWO)
      k=1
      do while(nprimbas(k,iact).ne.0)
        ipr=nprimbas(k,iact)
        nn=nlm(ipr,1)
        ll=nlm(ipr,2)
        mm=nlm(ipr,3)
        alpha=expp(ipr)
        if(ixyz.eq.1) then
          dx=-TWO*alpha*(x**(nn+1))
          if(nn.ge.1) dx=dx+nn*x**(nn-1)
          dx=dx*(y**ll)*(z**mm)*dexp(-expp(ipr)*(rr**2))
        else if(ixyz.eq.2) then
          dx=-TWO*alpha*(y**(ll+1))
          if(ll.ge.1) dx=dx+ll*y**(ll-1)
          dx=dx*(x**nn)*(z**mm)*dexp(-expp(ipr)*(rr**2))
        else
          dx=-TWO*alpha*(z**(mm+1))
          if(mm.ge.1) dx=dx+mm*z**(mm-1)
          dx=dx*(x**nn)*(y**ll)*dexp(-expp(ipr)*(rr**2))
        end if
        fx=fx+dx*coefpb(ipr,iact)
        k=k+1
      enddo
      drho_xyz=fx

      end

! *****

