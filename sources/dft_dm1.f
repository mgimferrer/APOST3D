      subroutine dft_dm1(itotps,wp,omp2,pcoord,chp)
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /filename/name0
      common /iops/iopt(200)
      common /dm1opt/densthresh_dm1
      character*80 ofile,ofile2
      character*60 name0
      integer*8 :: npairtot,npairskip

      dimension :: wp(itotps),omp2(itotps,nat),pcoord(itotps,3),chp(itotps,igr)
      dimension :: exch_hf(maxat,maxat)

!! automatic (stack, non-allocatable) per-pair scratch for the main RDM1 !!
!! double loop below -- sized by igr/nalf/nb, known at subroutine entry. !!
!! Kept out of the allocatable lists on purpose: an OMP-PRIVATE          !!
!! allocatable array gets each thread an unassociated copy that would    !!
!! need its own per-thread ALLOCATE; a plain automatic array just works, !!
!! same fix already used for enpart.f's multipolar (see its rvect).      !!
      dimension :: eval_ao(igr),gx_ao(igr),gy_ao(igr),gz_ao(igr)
      dimension :: chp2v(nalf),chp2bv(nb)
      dimension :: rhoab(2,1),scrpt(3,1),excpt(1),excbpt(1)

      allocatable :: wppha(:),omp2pha(:,:),pcoordpha(:,:),chppha(:,:),omp(:),ibaspoint(:)
      allocatable :: scr(:,:),scrpha(:,:)
      allocatable :: chp2(:,:),chp2pha(:,:),rho(:,:),rhopha(:,:),exc(:),excpha(:)
      allocatable :: chp2b(:,:),chp2phab(:,:),excb(:),excbpha(:)
      allocatable :: rdm1(:,:),rdm1b(:,:),dm1_mo(:,:),dm1_norm(:,:)
      allocatable :: rhoscr(:)

      ofile  = trim(name0)//".dm1"
      ofile2 = trim(name0)//".dm1norm"
      ifunc  = Iopt(57)
      iatps  = nrad*nang

      if(kop.ne.1) stop " Only implemented for unrestricted WF "

!! INITIAL INFORMATION PRINTING !!

      write(*,*) " "
      write(*,*) " COMPUTING DM1 APPROXIMATION FOR RKS-DFT FUNCTIONALS "
      write(*,*) " "
      write(*,*) " GENERAL INFORMATION "
      write(*,*) " "
!     write(*,*) " OUTPUT FILE FOR DM1 : ",trim(ofile)
!     write(*,*) " OUTPUT FILE FOR NORMALIZED DM1 : ",trim(ofile2)
      write(*,*) " NUMBER OF BASIS FUNCTIONS : ",igr
      write(*,*) " NUMBER OF OCCUPIED ALPHA MOs : ",nalf
      write(*,*) " NUMBER OF OCCUPIED BETA MOs : ",nb
      write(*,*) " NUMBER OF RADIAL POINTS : ",nrad
      write(*,*) " NUMBER OF ANGULAR (PER RADIAL) POINTS : ",nang
      write(*,*) " R0 PARAMETER FOR THE GAUSS-LEGENDRE QUADRATURE : ",rr00
      if(ifunc.eq.999) then 
        write(*,*) " HF FUNCTIONAL SELECTED "
        write(*,*) " "
      else
        call func_info_print(ifunc,itype,1)
      end if

!! GENERATING SECOND GRID (ROTATED) FOR NUMERICAL INTEGRATION !!
!! chp: VALUE OF jj AO IN THE ii POINT (FIRST GRID) !!
!! chppha: VALUE OF jj AO IN THE ii POINT (SECOND GRID) !!

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

!! TRANSFORMATION TO MOs !!

      ALLOCATE(chp2(itotps,nalf),chp2pha(itotps,nalf),exc(itotps),excpha(itotps))
      ALLOCATE(chp2b(itotps,nb),chp2phab(itotps,nb),excb(itotps),excbpha(itotps))
      ALLOCATE(rho(2,itotps),rhopha(2,itotps))
      do kk=1,itotps
        xx0=ZERO
        xx0b=ZERO
        xx0pha=ZERO
        xx0phab=ZERO
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
          xx0=xx0+xx*xx
          xx0pha=xx0pha+xxpha*xxpha
          if(imo.le.nb) then
            chp2b(kk,imo)=xxb
            chp2phab(kk,imo)=xxphab
            xx0b=xx0b+xxb*xxb
            xx0phab=xx0phab+xxphab*xxphab
          end if 
        end do
        rho(1,kk)=xx0
        rho(2,kk)=xx0b
        rhopha(1,kk)=xx0pha
        rhopha(2,kk)=xx0phab
      end do

!! CHECKING ONE- AND TWO-ELECTRON NUMERICAL INTEGRATION ACCURACY !!

      if(ifunc.ne.999) then
        if(itype.gt.1) then
          ALLOCATE(scr(3,itotps),scrpha(3,itotps))
          call sigma_uks(pcoord,chp2,chp2b,scr)
          call sigma_uks(pcoordpha,chp2pha,chp2phab,scrpha)
        end if
        call xc_uks_for_dm1(1,itotps,ifunc,rho,scr,exc)
        call xc_uks_for_dm1(2,itotps,ifunc,rho,scr,excb)
!! debug check requested by Marti -- same LDA/GGA exchange-energy-       !!
!! density evaluation, applied to the ROTATED grid's own density instead !!
!! of the first grid's.                                                 !!
        call xc_uks_for_dm1(1,itotps,ifunc,rhopha,scrpha,excpha)
        call xc_uks_for_dm1(2,itotps,ifunc,rhopha,scrpha,excbpha)
        if(itype.gt.1) DEALLOCATE(scr,scrpha)
      end if

      xx=ZERO
      xxpha=ZERO
      xlsda=ZERO
      xlsdapha=ZERO
      do icenter=1,nat
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          xw=wp(ifut)*omp2(ifut,icenter)
          xwpha=wppha(ifut)*omp2pha(ifut,icenter)
          xx=xx+xw*(rho(1,ifut)+rho(2,ifut))
          xxpha=xxpha+xwpha*(rhopha(1,ifut)+rhopha(2,ifut))
          if(ifunc.ne.999) then
            xlsda=xlsda+xw*(exc(ifut)+excb(ifut))
            xlsdapha=xlsdapha+xwpha*(excpha(ifut)+excbpha(ifut))
          end if
        end do
      end do
      write(*,*) " Integrated Density from First Grid (Alpha+Beta) : ",xx
!! debug checks requested by Marti -- confirm the ROTATED grid alone     !!
!! (same mechanism used as "r2" throughout the double loop below) is a   !!
!! valid, correctly-weighted representation of the density/exchange-     !!
!! energy on its own, independent of any r1/r2 pairing question.         !!
      write(*,*) " Integrated Density from Rotated Grid (Alpha+Beta) : ",xxpha
      if(ifunc.ne.999) then
        write(*,*) " One-el KS-Exchange (Alpha+Beta) : ",xlsda
        write(*,*) " One-el KS-Exchange from Rotated Grid (Alpha+Beta) : ",xlsdapha
      end if

      write(*,*) " "
      write(*,*) " CHECKING TWO-ELECTRON INTEGRALS "
      write(*,*) " "

!! ALWAYS USING THE ROTATED GRID FOR SECOND ELECTRON !!

      f2=ZERO
      do icenter=1,nat
        do jcenter=1,nat
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            x0=wp(ifut)*omp2(ifut,icenter)
            f3=ZERO
            do jfut=iatps*(jcenter-1)+1,iatps*jcenter
              x1=wppha(jfut)*omp2pha(jfut,jcenter)
              f3=f3+rho(1,ifut)*rhopha(1,jfut)*x1*x0
            end do
            f2=f2+f3
          end do
        end do
      end do
      write(*,*) " Calculated square number of electrons (N^2) : ",f2
      write(*,*) " "

      DEALLOCATE(rho,rhopha,exc,excb,excpha,excbpha)

!! CALCULATION OF THE HF RDM1 !!

!     write(*,*) " large allocating matrices "
!     ALLOCATE(rdm1(itotps,itotps),rdm1b(itotps,itotps))
      ALLOCATE(rdm1(1,1),rdm1b(1,1))
!     write(*,*) " done "
      if(ifunc.eq.999) then
        write(*,*) " GENERATING HF RDM1 "
        write(*,*) " "
!! can be removed... if not, need to be adapted to unrestricted !!
!       do icenter=1,nat
!         do jcenter=1,nat
!           do ifut=iatps*(icenter-1)+1,iatps*icenter
!             x0=wp(ifut)*omp2(ifut,icenter)
!             do jfut=iatps*(jcenter-1)+1,iatps*jcenter
!               x1=wppha(jfut)*omp2pha(jfut,jcenter)
!               do imo=1,nocc
!                 do imo2=imo,nocc

!! ONLY ALPHA PART !!

!                   xf=chp2(ifut,imo)*chp2pha(jfut,imo2)
!                   if(imo2.ne.imo) xf=TWO*xf
!                   rdm1(ifut,jfut)=rdm1(ifut,jfut)+xf
!                 end do
!               end do
!             end do
!           end do
!         end do
!       end do
!! ------- !!

!! NOW FOR KS-DFT FUNCTIONALS !!

      else

!! CALCULATION OF THE DENSITY AND ITS GRADIENTS AT THE R ((r1+r2)/2) POINTS !!

        write(*,*) " GENERATING KS-DFT RDM1 "
        write(*,*) " "
        xexch=ZERO
        xexchb=ZERO
        npairtot=0
        npairskip=0
        do icenter=1,nat
          do jcenter=1,nat
            f3=ZERO
            f3b=ZERO

!! parallel over ifut: gx_ao/gy_ao/gz_ao/eval_ao/chp2v/chp2bv/rhoab/scrpt/ !!
!! excpt/excbpt are all automatic (not allocatable, see declaration       !!
!! above), so PRIVATE gives each thread its own real copy, no per-thread  !!
!! (re)allocation needed. gpoints/drho_xyz/sigma_uks_xyz pass which basis  !!
!! function is "current" through common/actual/ (qtaim.f, dft_dm1.f) --   !!
!! THREADPRIVATE'd at their own declarations so each thread gets its own. !!
!! c/cb (ao_matrices) are shared but read-only here. jfut/f3/f3b/xexch/   !!
!! xexchb/npairtot/npairskip are the only cross-iteration accumulators,   !!
!! all via REDUCTION; everything else is written fresh every ifut.        !!
!$OMP PARALLEL DO PRIVATE(jfut,Rx,Ry,Rz,rhoa,rhob,rhoab,imo,ibf,xx,xxb,
!$OMP&  scraa,scrab,scrbb,scrpt,xfact,excpt,excbpt,r12,x1,xx1,xx1b,
!$OMP&  xksigaa,xksigbb,xx0,xx0b,xkagga,xkbgga,xbf,xbfb,xxrdm1,xxrdm1b,
!$OMP&  x0,eval_ao,gx_ao,gy_ao,gz_ao,chp2v,chp2bv)
!$OMP&  REDUCTION(+:f3,f3b,xexch,xexchb,npairtot,npairskip)
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
              end do
            end do
!$OMP END PARALLEL DO
            if(icenter.ne.jcenter) then
              f3=TWO*f3
              f3b=TWO*f3b
            end if
            exch_hf(icenter,jcenter)=-(f3+f3b)/TWO
          end do
        end do

        write(*,*) " "
        write(*,'(2x,a,1x,i14)') "Grid-point pairs evaluated:",npairtot
        write(*,'(2x,a,1x,i14,1x,a,1x,f6.2,1x,a)') "Pairs pruned (R below DENSTHRESH):",
     $npairskip,"(",100.0d0*dble(npairskip)/dble(npairtot),"%)"
        write(*,'(2x,a,1x,i14,1x,a,1x,f6.2,1x,a)') "Pairs kept (full cost paid):",
     $npairtot-npairskip,"(",100.0d0*dble(npairtot-npairskip)/dble(npairtot),"%)"
        write(*,*) " "

      end if
      call flush 

!! PROJECTION OF THE RDM1 INTO THE MOs FORMING DM1 !!

!     write(*,*) " COMPUTING DM1 (in MO's) DIRECTLY "
!     write(*,*) " "
!     ALLOCATE(dm1_mo(igr,igr))
!     do imo=1,igr
!       do imo2=imo,igr
!         xdm1=ZERO
!         do icenter=1,nat
!           do jcenter=1,nat
!             do ifut=iatps*(icenter-1)+1,iatps*icenter
!               x0=wp(ifut)*omp2(ifut,icenter)
!               do jfut=iatps*(jcenter-1)+1,iatps*jcenter
!                 x1=wppha(jfut)*omp2pha(jfut,jcenter)
!                 xf=chp2(ifut,imo)*chp2pha(jfut,imo2)
!                 xdm1=xdm1+rdm1(ifut,jfut)*xf*x0*x1
!               end do
!             end do
!           end do
!         end do
!         dm1_mo(imo,imo2)=xdm1
!         dm1_mo(imo2,imo)=xdm1
!       end do
!     end do

!! EVALUATING THE TRACE OF THE DM1 OBTAINED !!

!     xtr=ZERO
!     do ii=1,igr
!       do jj=1,igr
!         if(ABS(dm1_mo(ii,jj)).gt.1.0d-4) write(*,*) " i,j,dm1(i,j) : ",ii,jj,dm1_mo(ii,jj)
!       end do
!       xtr=xtr+dm1_mo(ii,ii)
!     end do
!     write(*,*) " "
!     write(*,*) " TRACE OF THE DM1 MATRIX : ",xtr
!     write(*,*) " "
!     call flush

!! NORMALIZING THE DM1 !!

!     ALLOCATE(dm1_norm(igr,igr))
!     xnocc=REAL(nocc)*TWO
!     xx0=ZERO
!     do ii=1,igr
!       do jj=1,igr
!         dm1_norm(ii,jj)=dm1_mo(ii,jj)*xnocc/xtr
!       end do
!       xx0=xx0+dm1_norm(ii,ii)
!     end do
!     write(*,*) " TRACE OF THE NORMALIZED DM1 : ",xx0
!     write(*,*) " "
!     call flush

!! COMPUTING THE EXCHANGE ENERGY !!

!     write(*,*) " EXCHANGE ENERGY CALCULATION "
!     write(*,*) " "
!     xexch=ZERO
!     xexchb=ZERO
!     xexch2=ZERO
!     xexch3=ZERO
!     do icenter=1,nat
!       do jcenter=1,nat
!         f3=ZERO
!         f3b=ZERO
!         do ifut=iatps*(icenter-1)+1,iatps*icenter
!           x0=wp(ifut)*omp2(ifut,icenter)
!           do jfut=iatps*(jcenter-1)+1,iatps*jcenter
!             r12=(pcoord(ifut,1)-pcoordpha(jfut,1))**TWO
!             r12=r12+((pcoord(ifut,2)-pcoordpha(jfut,2))**TWO)
!             r12=r12+((pcoord(ifut,3)-pcoordpha(jfut,3))**TWO)
!             r12=dsqrt(r12)
!             x1=wppha(jfut)*omp2pha(jfut,jcenter)

!! COMPUTING ONLY ONCE FROM RDM1 !!

!             if(r12.gt.thresh) then
!               xexch=xexch-(rdm1(ifut,jfut)*rdm1(ifut,jfut)*x0*x1/r12)
!               xexchb=xexchb-(rdm1b(ifut,jfut)*rdm1b(ifut,jfut)*x0*x1/r12)
!               f3=f3+(rdm1(ifut,jfut)*rdm1(ifut,jfut)*x0*x1/r12)
!               f3b=f3b+(rdm1b(ifut,jfut)*rdm1b(ifut,jfut)*x0*x1/r12)
!             end if

!! NOW FROM DM1 AND NORMALIZED DM1 !!

!             if(r12.gt.thresh) then
!               do imo=1,igr
!                 do imo2=imo,igr
!                   xx0=dm1_mo(imo,imo2)*chp2(ifut,imo)*chp2pha(jfut,imo2)
!                   xx1=dm1_norm(imo,imo2)*chp2(ifut,imo)*chp2pha(jfut,imo2)
!                   if(imo2.ne.imo) then
!                     xx0=TWO*xx0
!                     xx1=TWO*xx1
!                   end if
!                   xexch2=xexch2-xx0*xx0*x0*x1/r12
!                   xexch3=xexch3-xx1*xx1*x0*x1/r12
!                  end do 
!                end do 
!             end if
!           end do
!         end do
!         if(icenter.ne.jcenter) then
!           f3=TWO*f3
!           f3b=TWO*f3b
!         end if 
!         exch_hf(icenter,jcenter)=-(f3+f3b)/TWO

!         if(ifunc.eq.999) exch_hf(icenter,jcenter)=exch_hf(icenter,jcenter)/TWO
!         exch_hf(jcenter,icenter)=exch_hf(icenter,jcenter)
!       end do
!     end do

!! checkings that can be removed !!
      if(ifunc.eq.999) then
        write(*,*) " EXCHANGE ENERGY (RDM1) : ",-xexch/FOUR
!       write(*,*) " EXCHANGE ENERGY (DM1) : ",xexch2/FOUR
!       write(*,*) " EXCHANGE ENERGY (NORMALIZED DM1) : ",xexch3/FOUR
        write(*,*) " "

      else
        write(*,*) " EXCHANGE ENERGY (RDM1) : ",(xexch+xexchb)/TWO
!       write(*,*) " EXCHANGE ENERGY (DM1) : ",xexch2/TWO
!       write(*,*) " EXCHANGE ENERGY (NORMALIZED DM1) : ",xexch3/TWO
        write(*,*) " "
      end if

!! printing terms !!
      write(*,*) " EXCHANGE ENERGY CONTRIBUTIONS FROM HIRAO'S RDM1 "
      write(*,*) " "
      call Mprint(exch_hf,nat,maxat)
      write(*,*) " "

!! PRINTING !!
! Unformatted printing of the DM1_MO

!     open(unit=2,file=ofile,status='unknown',form='unformatted')
!     do ii=1,igr
!      do jj=1,igr
!       write(2) 2*ii-1,2*jj-1,DM1_MO(ii,jj)
!       write(2) 2*ii,2*jj,DM1_MO(ii,jj)
!      end do
!     end do
!     write(2) 0,0,0.0E0
!     close(2)

! Unformatted printing of the DM1_NORM

!     open(unit=3,file=ofile2,status='unknown',form='unformatted')
!     do ii=1,igr
!      do jj=1,igr
!       write(3) 2*ii-1,2*jj-1,DM1_NORM(ii,jj)
!       write(3) 2*ii,2*jj,DM1_NORM(ii,jj)
!      end do
!     end do
!     write(3) 0,0,0.0E0
!     close(3)

! End of DM1_MO calculation in restricted case

!     DEALLOCATE(rho,rho1,rho2,c,chp,chp2,chp3,dchp)
!     DEALLOCATE(RDM1,DM1_MO,DM1_NORM,xdist)
!     DEALLOCATE(nquant,xI)
!     DEALLOCATE(rgrad)
!     DEALLOCATE(wp,pcoord)
!     DEALLOCATE(wp2,pcoord2)

      end

! *****

      subroutine gen_ksgga(ifunc,rho0,xksgga)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'

      if(ifunc.eq.1) then
        xkfunc=(THREE/(FOUR*pi))**(ONE/THREE)
        xkfunc=THREE*xkfunc
      else if(ifunc.eq.2) then
        write(*,*) " B88 CALCULATION"
      else if(ifunc.eq.3) then
        write(*,*) " PBE CALCULATION"
      else if(ifunc.eq.4) then
        write(*,*) " PKZB CALCULATION"
      else if(ifunc.eq.5) then
        write(*,*) " TPSS CALCULATION"
      end if

!! COMPUTING THE ksGGA !!

      xx0=(9.0d0*pi/xkfunc)**HALF
      xksgga=xx0*(rho0**(ONE/THREE))

      end

! *****

      subroutine xc_uks_for_dm1(isigma,npt,id_xfunc,scr_ab,scr,scr2)
      use xc_f90_types_m
      use xc_f90_lib_m
      implicit real*8(a-h,o-z)
      TYPE(xc_f90_pointer_t) :: xc_func
      TYPE(xc_f90_pointer_t) :: xc_info
      include 'parameter.h'

      dimension :: scr_ab(2,npt),scr(3,npt),scr2(npt)

!! scr_ab = RHO (ORDER: ALPHA, BETA) !!
!! scr    = SIGMA RHO (ORDER: UP-UP, UP-DOWN, DOWN-DOWN) !!
!! scr2   = EX (PER ELECTRON, MULTIPLICATION BY RHO REMAINING) !!

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

!     subroutine gen_ksigmagga_uks(itype,itotps,chp2,chp2b,rho,xksiggga)
!     IMPLICIT REAL*8(A-H,O-Z)
!     include 'parameter.h'
!     common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
!     common /iops/iopt(100)

!     dimension :: chp2(itotps,nalf),chp2b(itotps,nb),rho(2,itotps)
!     dimension :: scr(3,itotps),exc(itotps),excb(itotps)
!     dimension :: xksiggga(2,itotps)

!     ifunc  = Iopt(57)

!     if(itype.gt.1) call sigma_uks(itotps,chp2,chp2b,scr)
!     call xc_uks_for_dm1(1,itotps,ifunc,rho,scr,exc)
!     call xc_uks_for_dm1(2,itotps,ifunc,rho,scr,excb)

!     xfact=FOUR/THREE
!     do ii=1,itotps
!       xksiggga(1,ii)=-TWO*exc(ii)/(rho(1,ii)**xfact)
!       xksiggga(2,ii)=-TWO*excb(ii)/(rho(2,ii)**xfact)
!     end do

!     end

! *****

!! moved here from tools.f -- DM1/HIRAO-specific helpers belong with their !!
!! driver, matching how enpart_dft.f/oslo.f keep each feature's private   !!
!! subroutines in the same file. Not called by dft_dm1 today (it's UHF-   !!
!! only) but kept for a possible future restricted/CASSCF DM1 variant.   !!

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

!! single-point UKS sigma (grad-rho.grad-rho contractions), needed at the !!
!! (r1+r2)/2 midpoints dft_dm1's double loop evaluates -- distinct from   !!
!! enpart_dft.f's sigma_uks, which only ever runs on a precomputed grid.  !!

      subroutine sigma_uks_xyz(xabs,yabs,zabs,chp2,chp3,scraa,scrab,scrbb)
      use ao_matrices
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /actual/ iact,jat,icenter
!! DFT-DM1's double loop (dft_dm1.f) calls this under OMP -- each thread  !!
!! needs its own iact, not one shared across all of them.                !!
!$OMP THREADPRIVATE(/actual/)
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

!! single-point AO gradient (ixyz component), one basis function (iact,   !!
!! via common/actual/, same convention qtaim.f's gpoints/gxfunct use) per !!
!! call. Rewritten against basis_set's nlm/coefpb/nprimbas/ihold/coord -- !!
!! same primitive-loop math as enpart_dft.f's sigma_uks, just for a       !!
!! single arbitrary point instead of the whole precomputed grid. Replaces !!
!! the old common/data//coeff//lim//hold/-based version, which targeted   !!
!! commons that no longer exist in this codebase.                        !!

      function drho_xyz(xabs,yabs,zabs,ixyz)
      use basis_set
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /actual/ iact,jat,icenter
!! called (via gxfunct/gyfunct/gzfunct and sigma_uks_xyz) from DFT-DM1's !!
!! OMP-parallelized loop -- each thread needs its own iact.              !!
!$OMP THREADPRIVATE(/actual/)

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

