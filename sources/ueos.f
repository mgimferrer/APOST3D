!! *********************************************************************** !!
!! GEOS ("Generalized Effective Oxidation State", `.inp` keyword GEOS,     !!
!! formerly EOS-U) -- effective atomic/fragment orbitals from the paired   !!
!! and unpaired densities separately (open-shell systems), then            !!
!! electron/oxidation-state assignment across both channels together.      !!
!!   effao3d_u    -- real-space (3D grid) EFOs from the paired/unpaired    !!
!!                   densities (Takatsuka's definition), one fragment at   !!
!!                   a time                                                !!
!!   ueos_analysis -- pools every fragment's paired/unpaired EFOs, assigns !!
!!                   electrons to minimize the RMSD against ideal (2 for   !!
!!                   paired, 1 for unpaired) occupations, then derives     !!
!!                   fragment oxidation states                             !!
!! *********************************************************************** !!

!! ***** !!

!! ********************************************************************* !!
!! subroutine: effao3d_u                                                 !!
!! purpose: computes real-space (3D grid) effective fragment orbitals    !!
!!   (EFOs) separately from the paired and unpaired one-particle          !!
!!   densities, using Takatsuka's definition of the unpaired density      !!
!!   (u = n(2-n) per natural-orbital occupation n; paired = total - u).   !!
!!   Same per-fragment scheme as ueffao3d_frag (effao.f), run twice       !!
!!   (icase=1 paired, icase=2 unpaired). Results are stored into          !!
!!   effao_mod (p0/p0net/p0gro/ip0) and, if iueos=1, into the local       !!
!!   up0net/up0gro/iup0 arrays consumed by ueos_analysis below.           !!
!! arguments:                                                             !!
!!   itotps (in) -- total number of grid points (nat*iatps)               !!
!!   ndim   (in) -- number of basis functions (leading dim of chp/sat)    !!
!!   omp    (in) -- becke/tfvc weight of each grid point for its own atom !!
!!   chp    (in) -- basis-function values at each grid point              !!
!!   sat    (in) -- per-atom AO overlap matrix (from numint_sat)          !!
!!   wp     (in) -- integration weight of each grid point                 !!
!!   omp2   (in) -- becke/tfvc (or hirshfeld) weight of each point for    !!
!!                  every atom                                            !!
!!   iueos  (in) -- 1 to also run the electron/oxidation-state assignment !!
!!                  (ueos_analysis) once both densities are done          !!
!! author:                                                                 !!
!! ********************************************************************* !!
      subroutine effao3d_u(itotps,ndim,omp,chp,sat,wp,omp2,iueos)

      use basis_set
      use ao_matrices
      use integration_grid
      use effao_mod, only: p0,p0net,p0gro,ip0

      implicit real*8(a-h,o-z)
      include 'parameter.h'

      integer,intent(in) :: itotps,ndim

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /iops/iopt(200)

      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat),wp(itotps)
      dimension sat(ndim,ndim,nat)

      allocatable :: Pno(:,:),Uno(:,:)
      allocatable :: S0(:,:),Sm(:,:),Splus(:,:),c0(:,:),pp0(:,:)
      allocatable :: scr(:),s0all(:)
      allocatable :: iup0(:,:),up0net(:,:,:),up0gro(:,:,:),up0coef(:,:,:,:)
      allocatable :: poolcoef(:,:,:)
      character(len=30) :: lbl30
      character(len=20) :: ctype

      icube = iopt(13)
      iatps = nang*nrad

!! EFO occupation cutoff for the net-population sum below !!
      xminocc=1.0d-4

      call print_box('DOING EFFAO-3D FROM U FUNCTION')
      write(*,'(2x,a)') 'EFFAO-U: paired and unpaired densities treated separately'

!! build the total (Pno) and Takatsuka unpaired (Uno) density matrices,  !!
!! in the AO basis, from the natural orbitals already in ao_matrices.    !!
      nnorb=0
      do ii=1,igr
        if(occ_no(ii,ii).ge.thresh) nnorb=nnorb+1
      end do
      write(*,'(2x,a,1x,i0)') 'Number of natural orbitals (NOs):',nnorb

      ALLOCATE(Pno(igr,igr),Uno(igr,igr))
      do mu=1,igr
        do nu=1,igr
          xx=ZERO
          xx2=ZERO
          do ii=1,nnorb
            xx=xx+occ_no(ii,ii)*c_no(mu,ii)*c_no(nu,ii)
            xx2=xx2+occ_no(ii,ii)*(TWO-occ_no(ii,ii))*c_no(mu,ii)*c_no(nu,ii)
          end do
          Pno(mu,nu)=xx
          Uno(mu,nu)=xx2
        end do
      end do

      ALLOCATE(scr(itotps))
      ALLOCATE(S0(igr,igr),s0all(igr),Sm(igr,igr),Splus(igr,igr))
      ALLOCATE(c0(igr,igr),pp0(igr,igr))
      if(iueos.eq.1) ALLOCATE(iup0(2,icufr),up0net(2,igr,icufr),up0gro(2,igr,icufr),
     +  up0coef(igr,igr,2,icufr),poolcoef(igr,igr,2))

!! icase=1: paired density (total-unpaired). icase=2: unpaired density. !!
      do icase=1,2
        if(icase.eq.1) then
          call print_box('EFFAOs FROM THE PAIRED DENSITY')
        else if(icase.eq.2) then
          call print_box('EFFAOs FROM THE UNPAIRED DENSITY')
        end if

!! per-fragment loop kept serial on purpose: icufr can be small on a     !!
!! large system -- the O(igr^2)/O(igr^3) work inside each iteration is   !!
!! threaded instead, same lesson already applied throughout effao.f.     !!
        do iicenter=1,icufr

!! W_A (fragment iicenter's total becke/tfvc weight at each grid point) !!
          scr=ZERO
          do ifut=1,itotps
            do icenter=1,nfrlist(iicenter)
              scr(ifut)=scr(ifut)+omp2(ifut,ifrlist(icenter,iicenter))
            end do
          end do

!! net AO overlap block (S0), ALLPOINTS integration. parallel over mu:  !!
!! each mu writes only its own S0(mu,*)/S0(*,mu), no two mu iterations  !!
!! collide. wp/chp/scr/omp shared read-only, xx private.                !!
!$OMP PARALLEL DO PRIVATE(mu,nu,jcenter,ifut,xx)
          do mu=1,igr
            do nu=1,mu
              xx=ZERO
              do jcenter=1,nat
                do ifut=iatps*(jcenter-1)+1,iatps*jcenter
                  xx=xx+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*scr(ifut)*scr(ifut)*omp(ifut)
                end do
              end do
              S0(mu,nu)=xx
              S0(nu,mu)=xx
            end do
          end do
!$OMP END PARALLEL DO
          call build_Smp(igr,S0,Sm,Splus,0)

!! select the paired or unpaired AO density for this icase !!
          do ii=1,igr
            do jj=1,igr
              if(icase.eq.1) pp0(ii,jj)=Pno(ii,jj)-Uno(ii,jj)
              if(icase.eq.2) pp0(ii,jj)=Uno(ii,jj)
            end do
          end do

          call to_lowdin_basis(igr,Splus,pp0)
          call diagonalize(igr,igr,pp0,c0,0)
          call to_AO_basis(igr,igr,Sm,c0)

!! keep EFOs above xminocc; pp0's diagonal is already sorted decreasing.  !!
!! imaxo starts at 0 so a channel with no EFO above threshold (e.g. the  !!
!! unpaired channel on a restricted wavefunction, where Uno is ~0 for    !!
!! integer NO occupations) correctly ends up with imaxo=0 instead of     !!
!! carrying over a stale value from the previous fragment/icase.         !!
          imaxo=0
          ii=1
          do while(pp0(ii,ii).ge.xminocc.and.ii.le.igr)
            imaxo=ii
            ii=ii+1
          end do
          xmaxo=ZERO
          do ii=1,igr
            xmaxo=xmaxo+pp0(ii,ii)
          end do

          write(*,'(2x,a11,x,i3,x,a2)') "** FRAGMENT",iicenter,"**"
          write(*,*) " "
          lbl30="Net occupation for fragment"
          write(*,'(2x,a30,i4,f11.5)') lbl30,iicenter,xmaxo
          write(*,'(2x,a22,x,f10.5)') "Net occupation using >",xminocc
          write(*,60) (pp0(mu,mu),mu=1,imaxo)
          write(*,*) " "

!! gross occupation of each EFO via sat, scaled by its net occupation.  !!
!! parallel over ii: independent per EFO, writes only s0all(ii); xx0 is !!
!! a genuine REDUCTION. c0/sat/pp0 shared, read-only.                   !!
          xx0=ZERO
!$OMP PARALLEL DO PRIVATE(ii,icenter,jcenter,jj,kk,xx,xxx) REDUCTION(+:xx0)
          do ii=1,imaxo
            xxx=ZERO
            do icenter=1,nfrlist(iicenter)
              jcenter=ifrlist(icenter,iicenter)
              xx=ZERO
              do jj=1,igr
                do kk=1,igr
                  xx=xx+c0(kk,ii)*sat(kk,jj,jcenter)*c0(jj,ii)
                end do
              end do
              xxx=xxx+xx
            end do
            xxx=xxx*pp0(ii,ii)
            s0all(ii)=xxx
            xx0=xx0+xxx
          end do
!$OMP END PARALLEL DO
          lbl30="Gross occupation for fragment"
          write(*,'(2x,a30,i4,f11.5)') lbl30,iicenter,xx0
          write(*,60) (s0all(mu),mu=1,imaxo)

!! skip this trailing blank on the very last fragment of icase=1 -- the  !!
!! next thing printed is icase=2's own print_box, which already opens    !!
!! with a leading blank of its own (see Code Style: no double blanks     !!
!! around print_box).                                                    !!
          if(.not.(icase.eq.1.and.iicenter.eq.icufr)) write(*,*) " "

!! store for cube generation, and for ueos_analysis if requested !!
          do kk=1,imaxo
            do mu=1,igr
              p0(mu,kk)=c0(mu,kk)
            end do
            p0net(kk,iicenter)=pp0(kk,kk)
            p0gro(kk,iicenter)=s0all(kk)
            if(iueos.eq.1) then
              up0net(icase,kk,iicenter)=pp0(kk,kk)
              up0gro(icase,kk,iicenter)=s0all(kk)
              do mu=1,igr
                up0coef(mu,kk,icase,iicenter)=c0(mu,kk)
              end do
            end if
          end do
          ip0(iicenter)=imaxo
          if(iueos.eq.1) iup0(icase,iicenter)=imaxo

!! optional: write a cube file for this fragment's paired/unpaired EFOs !!
          if(icase.eq.1) iicase=3
          if(icase.eq.2) iicase=4
          if(icube.eq.1) call cubegen_new(iicenter,iicase)
        end do
      end do

!! electron/oxidation-state assignment across both channels together,   !!
!! using gross populations                                              !!
      if(iueos.eq.1) then
        call ueos_analysis(iup0,up0gro,up0coef,poolcoef)

!! pooled paired/unpaired EFOs as fake Alpha/Beta MOs in a .fchk, for    !!
!! visualization in any standard viewer -- printed by default, same as  !!
!! OSLO's own .fchk output, no separate keyword needed. Restricted      !!
!! wavefunctions have no Beta blocks to splice into at all, hence the   !!
!! kop branch (the unpaired channel is then trivially ~empty, same as   !!
!! for a plain restricted-wavefunction EFFAO run).                      !!
        ctype="-GEOS-EFOs"
        if(kop.eq.0) then
          call rwf_geos_orbprint(poolcoef(:,:,1),ctype)
        else
          call uwf_geos_orbprint(poolcoef(:,:,1),poolcoef(:,:,2),ctype)
        end if
      end if

      DEALLOCATE(S0,Splus,Sm)
      DEALLOCATE(scr,c0,pp0,s0all)
      DEALLOCATE(Pno,Uno)
      DEALLOCATE(iup0,up0net,up0gro,up0coef,poolcoef)

60    FORMAT("  OCCUP.",8f9.4)

      end

!! ****** !!

!! ********************************************************************* !!
!! subroutine: ueos_analysis                                             !!
!! purpose: pools every fragment's paired and unpaired EFO gross         !!
!!   occupations (from effao3d_u above), assigns integer electrons (2    !!
!!   per paired EFO, 1 per unpaired) to whichever pooled slots minimize  !!
!!   the RMSD against these occupations, then derives fragment           !!
!!   oxidation states from the resulting electron counts. Also pools     !!
!!   and sorts each EFO's coefficient vector the same way, truncated to  !!
!!   igr and zero-padded, feeding the GEOS .fchk splicer (rwf_geos_      !!
!!   orbprint/uwf_geos_orbprint) below.                                  !!
!! arguments:                                                            !!
!!   iup0    (in)  -- number of EFOs kept per (channel, fragment)        !!
!!   up0gro  (in)  -- gross occupation of each EFO, per (channel, index, !!
!!                    fragment)                                          !!
!!   up0coef (in)  -- AO coefficients of each EFO, per (basis fn, index, !!
!!                    channel, fragment)                                 !!
!!   poolcoef (out) -- pooled+sorted EFO coefficients, per (basis fn,    !!
!!                    pooled index truncated to igr, channel)            !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine ueos_analysis(iup0,up0gro,up0coef,poolcoef)

      implicit real*8(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /loba/ oxi(maxat),errsav(maxat),elec(maxat),effpop(maxat)
      common /loba2/occup(nmax,2),iorbat(nmax,2),lorb(2),confi0

      dimension iup0(2,icufr),up0gro(2,igr,icufr)
      dimension up0coef(igr,igr,2,icufr),poolcoef(igr,igr,2)
      dimension iorbslot(nmax,2)

      allocatable :: tmp_occup(:),tmp_iorbat(:),tmp_iorbslot(:)
      allocatable :: elec2(:,:),elec_id(:,:)
      allocatable :: tmp_elec_id(:,:)
      allocatable :: elec_frg_count(:,:)

!! pool every fragment's EFOs into one list per channel (1=paired,      !!
!! 2=unpaired); iorbslot keeps each pooled slot's local (within-        !!
!! fragment) EFO index, needed below to fetch its coefficient vector.   !!
      do icase=1,2
        iorb=0
        do ii=1,icufr
          do kk=1,iup0(icase,ii)
            iorb=iorb+1
            occup(iorb,icase)=up0gro(icase,kk,ii)
            iorbat(iorb,icase)=ii
            iorbslot(iorb,icase)=kk
          end do
        end do
        lorb(icase)=iorb
      end do
      write(*,'(2x,a38,x,i4)') "Total number of eff-AO-s for analysis:",lorb(1)+lorb(2)

!! sort each channel's pooled EFOs by decreasing occupation, keeping    !!
!! iorbat/iorbslot aligned so each slot still knows its source fragment !!
!! and local EFO index.                                                 !!
      do icase=1,2
        ALLOCATE(tmp_occup(lorb(icase)))
        ALLOCATE(tmp_iorbat(lorb(icase)))
        ALLOCATE(tmp_iorbslot(lorb(icase)))
        do ii=1,lorb(icase)
          tmp_occup(ii)=occup(ii,icase)
          tmp_iorbat(ii)=iorbat(ii,icase)
          tmp_iorbslot(ii)=iorbslot(ii,icase)
        end do

        do ii=1,lorb(icase)-1
          do jj=ii+1,lorb(icase)
            if(tmp_occup(ii).lt.tmp_occup(jj)) then
              tmp_swap=tmp_occup(ii)
              tmp_occup(ii)=tmp_occup(jj)
              tmp_occup(jj)=tmp_swap

              tmp_swap=tmp_iorbat(ii)
              tmp_iorbat(ii)=tmp_iorbat(jj)
              tmp_iorbat(jj)=tmp_swap

              tmp_swap=tmp_iorbslot(ii)
              tmp_iorbslot(ii)=tmp_iorbslot(jj)
              tmp_iorbslot(jj)=tmp_swap
            end if
          end do
        end do

        do ii=1,lorb(icase)
          occup(ii,icase)=tmp_occup(ii)
          iorbat(ii,icase)=tmp_iorbat(ii)
          iorbslot(ii,icase)=tmp_iorbslot(ii)
        end do
        DEALLOCATE(tmp_occup,tmp_iorbat,tmp_iorbslot)
      end do

!! build the final, igr-wide pooled coefficient matrix per channel from !!
!! the sorted order above -- zero-initialized, so any channel with      !!
!! fewer than igr pooled EFOs above threshold is correctly zero-padded  !!
!! (always the opposite in practice: more pooled EFOs than igr).        !!
      poolcoef=ZERO
      do icase=1,2
        do jj=1,MIN(lorb(icase),igr)
          ii=iorbat(jj,icase)
          kk=iorbslot(jj,icase)
          do mu=1,igr
            poolcoef(mu,jj,icase)=up0coef(mu,kk,icase,ii)
          end do
        end do
      end do

!! ideal-occupation matrix, sized to the larger of the two EFO counts   !!
!! (previously always lorb(1)/paired, which ran elec_id out of bounds   !!
!! whenever there were more unpaired than paired EFOs -- fixed 2026-08-21) !!
      ilorb=MAX(lorb(1),lorb(2))
      ALLOCATE(elec_id(ilorb,2))
      ALLOCATE(elec2(ilorb,2))
      elec_id=ZERO
!! section, not whole-array -- keeps elec2 at ilorb (a bare "=occup"    !!
!! would auto-reallocate to occup's own (nmax,2) shape, F2003 semantics) !!
      elec2=occup(1:ilorb,1:2)

!! ideal assignment: nunp electrons forced unpaired (|nalf-nb|), the    !!
!! rest as electron pairs                                               !!
      nnn=nalf+nb
      nunp=nalf-nb
      if(nunp.ne.0) then
        do ii=1,nunp
          elec_id(ii,2)=ONE
          nnn=nnn-1
        end do
      end if

      npair=nnn/2
      do ii=1,npair
        elec_id(ii,1)=TWO
        nnn=nnn-2
      end do
      write(*,'(2x,a,i0,a,i0,a)') 'Ideal assignment: ',nunp,
     +  ' forced-unpaired electron(s), ',npair,' electron pair(s)'

!! RMSD between the actual pooled occupations and this ideal assignment !!
      xrmsd=ZERO
      do icase=1,2
        do ii=1,lorb(icase)
          xx1=(elec2(ii,icase)-elec_id(ii,icase))
          xx1=xx1*xx1
          xrmsd=xrmsd+xx1
        end do
      end do
      xrmsd=xrmsd/REAL(nalf+nb)
      xrmsd=dsqrt(xrmsd)
      write(*,'(2x,a,f8.5)') 'Initial RMSD:',xrmsd
      write(*,*)

!! iteratively move the least-occupied paired electron pair to the two  !!
!! most-occupied unpaired slots, keeping the move only while it lowers  !!
!! the RMSD -- stops at the first move that doesn't help.               !!
      ALLOCATE(tmp_elec_id(ilorb,2))
      tmp_elec_id=elec_id
      do while(npair.gt.0)
        tmp_elec_id(npair,1)=ZERO
        tmp_elec_id(nunp+1,2)=ONE
        tmp_elec_id(nunp+2,2)=ONE

        xrmsd2=ZERO
        do icase=1,2
          do ii=1,lorb(icase)
            xx1=(elec2(ii,icase)-tmp_elec_id(ii,icase))
            xx1=xx1*xx1
            xrmsd2=xrmsd2+xx1
          end do
        end do
        xrmsd2=xrmsd2/REAL(nalf+nb)
        xrmsd2=dsqrt(xrmsd2)

        write(*,'(2x,a,i0,a,f7.4,a,i0,a,f7.4,a,i0,a,f7.4,a)')
     +    'Trying: paired EFO (frag ',iorbat(npair,1),', occ ',elec2(npair,1),
     +    ') -> unpaired (frag ',iorbat(nunp+1,2),', occ ',elec2(nunp+1,2),
     +    ' / frag ',iorbat(nunp+2,2),', occ ',elec2(nunp+2,2),')'
        if((xrmsd2-xrmsd).lt.ZERO) then
          write(*,'(4x,a,f8.5,a)') 'RMSD = ',xrmsd2,' -- accepted'
          elec_id=tmp_elec_id
          xrmsd=xrmsd2
        else
          write(*,'(4x,a,f8.5,a,f8.5)') 'RMSD = ',xrmsd2,
     +      ' -- rejected, keeping RMSD = ',xrmsd
          go to 69
        end if

        nunp=nunp+2
        npair=npair-1
      end do
69    continue
      write(*,*)
      DEALLOCATE(tmp_elec_id)

!! electrons assigned to each fragment, split by paired/unpaired for    !!
!! the per-channel EOS analysis below                                   !!
      ALLOCATE(elec_frg_count(icufr,2))
      elec=ZERO
      elec_frg_count=ZERO
      do icase=1,2
        do ii=1,lorb(icase)
          if(elec_id(ii,icase).gt.ZERO) then
            ifrg=iorbat(ii,icase)
            elec(ifrg)=elec(ifrg)+elec_id(ii,icase)
            elec_frg_count(ifrg,icase)=elec_frg_count(ifrg,icase)+elec_id(ii,icase)
          end if
        end do
      end do

      write(*,'(2x,a)') 'Electrons assigned per fragment (paired / unpaired):'
      do ifrg=1,icufr
        write(*,'(4x,a,i0,a,f5.2,a,f5.2)') 'Fragment ',ifrg,': ',
     +    elec_frg_count(ifrg,1),' / ',elec_frg_count(ifrg,2)
      end do

!! no trailing blank here -- the next thing printed is icase=1's own     !!
!! print_box below, which already opens with a leading blank.            !!
      do icase=1,2
        if(icase.eq.1) then
          call print_box('EOS ANALYSIS FOR PAIRED ELECTRONS')
        else if(icase.eq.2) then
          call print_box('EOS ANALYSIS FOR UNPAIRED ELECTRONS')
        end if
        write(*,*) "  Frag.  Elect.  Last occ.  First unocc.  "
        write(*,*) " ---------------------------------------- "

!! last occupied / first unoccupied EFO gross occupation per fragment,  !!
!! same convention as effao.f's eos_analysis                            !!
        if(icase.eq.1) xlast=TWO
        if(icase.eq.2) xlast=ONE
        ilast=0
        do ifrg=1,icufr
          if(icase.eq.1) then
            nn=INT(elec_frg_count(ifrg,icase))/2
            if(elec_frg_count(ifrg,icase)/TWO-nn.gt.thresh) nn=nn+1
          end if
          if(icase.eq.2) then
            nn=INT(elec_frg_count(ifrg,icase))
            if(elec_frg_count(ifrg,icase)-nn.gt.thresh) nn=nn+1
          end if

          if(nn+1.gt.iup0(icase,ifrg)) then
            write(*,10) ifrg,elec_frg_count(ifrg,icase),up0gro(icase,nn,ifrg)
          else
            if(nn.ne.0) then
              write(*,15) ifrg,elec_frg_count(ifrg,icase),up0gro(icase,nn,ifrg),up0gro(icase,nn+1,ifrg)
            else
              write(*,15) ifrg,elec_frg_count(ifrg,icase),ZERO,up0gro(icase,nn+1,ifrg)
            end if
          end if
          if(nn.ne.0) then
            if(up0gro(icase,nn,ifrg).lt.xlast) then
              xlast=up0gro(icase,nn,ifrg)
              ilast=ifrg
            end if
          end if
        end do

        xfirst=ZERO
        do ifrg=1,icufr
          if(ifrg.ne.ilast) then
            if(icase.eq.1) then
              nn=INT(elec_frg_count(ifrg,icase))/2
              if(elec_frg_count(ifrg,icase)/TWO-nn.gt.thresh) nn=nn+1
            end if
            if(icase.eq.2) then
              nn=INT(elec_frg_count(ifrg,icase))
              if(elec_frg_count(ifrg,icase)-nn.gt.thresh) nn=nn+1
            end if
            if(nn+1.ne.0.and.up0gro(icase,nn+1,ifrg).gt.xfirst) xfirst=up0gro(icase,nn+1,ifrg)
          end if
        end do
        write(*,*) " ---------------------------------------- "

!! reliability index: paired occupations run 0-2 instead of 0-1, so the !!
!! gap is halved before the same +0.5 scaling used for unpaired/EOS.    !!
        if(icase.eq.1) confi=100.0*dmin1(1.0d0,(xlast-xfirst)/TWO+0.5d0)
        if(icase.eq.2) confi=100.0*dmin1(1.0d0,xlast-xfirst+0.5d0)
        write(*,'(3x,a24,x,f7.3)') "RELIABILITY INDEX R(%) =",confi
        if(icase.eq.1) then
          write(*,'(2x,a)') 'INFO: (paired) occ. values halved in R(%) calculation'
          confi0=confi
        end if
      end do

      zztot=ZERO
      do ifrg=1,icufr
        zzn=ZERO
        do jfrg=1,nfrlist(ifrg)
          zzn=zzn+zn(ifrlist(jfrg,ifrg))
        end do
        oxi(ifrg)=zzn-elec(ifrg)
        zztot=zztot+oxi(ifrg)
      end do

      call print_box('FRAGMENT OXIDATION STATES')
      write(*,*) "  Frag.  Oxidation State  "
      write(*,*) " ------------------------ "
      do ifrg=1,icufr
        write(*,20) ifrg,oxi(ifrg)
      end do
      write(*,*) " ------------------------ "
      write(*,'(3x,a4,x,f4.1)') "Sum:",zztot
      write(*,*) " "

      confi2=dmin1(confi,confi0)
      write(*,'(2x,a32,x,f7.3)') "OVERALL RELIABILITY INDEX R(%) =",confi2
      write(*,'(2x,a34,x,f7.3)') "ELECTRONIC ASSIGNMENT RMSD VALUE =",xrmsd

      DEALLOCATE(elec_id,elec2,elec_frg_count)

10    FORMAT(3x,i3,3x,f6.2,4x,f12.3,'    < thresh',2f12.3)
15    FORMAT(3x,i3,3x,f6.2,4x,f6.3,4x,f6.3)
20    FORMAT(3x,i3,6x,f8.2)

      end

!! ****** !!

!! ********************************************************************* !!
!! subroutine: rwf_geos_orbprint                                         !!
!! purpose: writes the restricted-wavefunction GEOS .fchk -- splices the !!
!!   original .fchk's structure, replacing "Alpha Orbital Energies" and  !!
!!   "Alpha MO coefficients" with the pooled paired-channel EFOs (from   !!
!!   ueos_analysis), sorted by decreasing gross occupation and used in   !!
!!   place of a real orbital energy, copying everything else through     !!
!!   unchanged (including Total SCF Density -- these are visualization   !!
!!   orbitals, not a real wavefunction). Reuses the same splice pattern  !!
!!   as OSLO's rwf_orbprint (oslo.f), extended to also cover the Orbital !!
!!   Energies block. Called whenever GEOS+DOFRAGS runs on a restricted   !!
!!   wavefunction -- see uwf_geos_orbprint for the unrestricted twin.    !!
!! arguments:                                                             !!
!!   pcoef (in) -- (igr,igr) pooled+sorted paired-channel EFO            !!
!!                 coefficients, zero-padded beyond the actual count     !!
!!   ctype (in) -- filename suffix, e.g. "-GEOS-EFOs"                    !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine rwf_geos_orbprint(pcoef,ctype)

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /filename/name0
      common /loba2/occup(nmax,2),iorbat(nmax,2),lorb(2),confi0

      character*80 line
      character*60 name0,name1
      character*20 ctype

      dimension pcoef(igr,igr)
      allocatable :: energ(:)

      iqchem   = iopt(95)
      imokit   = iopt(79)
      indepigr = int_locate(15,"Number of independ",ilog)
      norb     = igr*indepigr

!! fake orbital energies: gross occupation of each pooled EFO, already  !!
!! sorted decreasing by ueos_analysis, zero beyond the actual count.    !!
      ALLOCATE(energ(indepigr))
      energ=ZERO
      do ii=1,MIN(lorb(1),indepigr)
        energ(ii)=occup(ii,1)
      end do

!! Name of the .fchk file !!
      name1=trim(name0)//trim(ctype)//".fchk"
      open(unit=69,file=name1)
      rewind(69)
      rewind(15)

      read(15,'(a80)') line

!! standard/MOKIT .fchk layout: Alpha Orbital Energies precedes Alpha   !!
!! MO coefficients. Q-Chem's is the other way around -- see uwf_geos_   !!
!! orbprint's header for the same iqchem branch used by OSLO's own      !!
!! printers.                                                             !!
      if(iqchem.eq.0) then
        do while(index(line,"Alpha Orbital").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do
        write(69,11) "Alpha Orbital Energies","R","N= ",indepigr
        write(69,13) (energ(ii),ii=1,indepigr)

        do while(index(line,"Alpha MO co").eq.0)
          read(15,'(a80)') line
        end do
        write(69,12) "Alpha MO coefficients","R","N= ",norb
        write(69,13) ((pcoef(ii,jj),ii=1,igr),jj=1,indepigr)

        if(imokit.eq.0) then
          do while(index(line,"Orthonormal basis").eq.0)
            read(15,'(a80)') line
          end do
        end if
      else
        do while(index(line,"Alpha MO co").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do
        write(69,12) "Alpha MO coefficients","R","N= ",norb
        write(69,13) ((pcoef(ii,jj),ii=1,igr),jj=1,indepigr)

        do while(index(line,"Alpha Orbital").eq.0)
          read(15,'(a80)') line
        end do
        write(69,11) "Alpha Orbital Energies","R","N= ",indepigr
        write(69,13) (energ(ii),ii=1,indepigr)
      end if

!! copy everything else through unchanged, Total SCF Density included !!
      do while(.true.)
        write(69,'(a80)') line
        read(15,'(a80)',end=99) line
      end do
99    continue
      close(69)
      DEALLOCATE(energ)

!! Printing formats !!
11    FORMAT(a23,20x,a1,3x,a3,i11)
12    FORMAT(a21,22x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))

      end

!! ****** !!

!! ********************************************************************* !!
!! subroutine: uwf_geos_orbprint                                         !!
!! purpose: writes the unrestricted-wavefunction GEOS .fchk -- splices   !!
!!   the original .fchk's structure, replacing "Alpha/Beta Orbital       !!
!!   Energies" and "Alpha/Beta MO coefficients" with the pooled paired   !!
!!   (-> Alpha) and unpaired (-> Beta) channel EFOs (from ueos_analysis),!!
!!   sorted by decreasing gross occupation and used in place of a real   !!
!!   orbital energy, copying everything else through unchanged           !!
!!   (including Total SCF/Spin SCF Density -- these are visualization    !!
!!   orbitals, not a real wavefunction). Reuses the same splice pattern  !!
!!   as OSLO's uwf_orbprint (oslo.f), extended to also cover the Orbital !!
!!   Energies blocks. Q-Chem-format branch is unverified -- no active    !!
!!   test exercises an unrestricted Q-Chem source, and the codebase's    !!
!!   own unrestricted OSLO printer only ever resumes after a single      !!
!!   "Alpha Orbital" marker for that format, never distinguishing an     !!
!!   Alpha/Beta split there.                                             !!
!! arguments:                                                             !!
!!   pcoef_a (in) -- (igr,igr) pooled+sorted paired-channel (-> Alpha)   !!
!!                   EFO coefficients, zero-padded beyond the actual     !!
!!                   count                                                !!
!!   pcoef_b (in) -- (igr,igr) pooled+sorted unpaired-channel (-> Beta)  !!
!!                   EFO coefficients, zero-padded beyond the actual     !!
!!                   count                                                !!
!!   ctype   (in) -- filename suffix, e.g. "-GEOS-EFOs"                  !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine uwf_geos_orbprint(pcoef_a,pcoef_b,ctype)

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /filename/name0
      common /loba2/occup(nmax,2),iorbat(nmax,2),lorb(2),confi0

      character*80 line
      character*60 name0,name1
      character*20 ctype

      dimension pcoef_a(igr,igr),pcoef_b(igr,igr)
      allocatable :: energ_a(:),energ_b(:)

      iqchem   = iopt(95)
      imokit   = iopt(79)
      indepigr = int_locate(15,"Number of independ",ilog)
      norb     = igr*indepigr

!! fake orbital energies: gross occupation of each pooled EFO, already  !!
!! sorted decreasing by ueos_analysis, zero beyond the actual count.    !!
      ALLOCATE(energ_a(indepigr),energ_b(indepigr))
      energ_a=ZERO
      energ_b=ZERO
      do ii=1,MIN(lorb(1),indepigr)
        energ_a(ii)=occup(ii,1)
      end do
      do ii=1,MIN(lorb(2),indepigr)
        energ_b(ii)=occup(ii,2)
      end do

!! Name of the .fchk file !!
      name1=trim(name0)//trim(ctype)//".fchk"
      open(unit=69,file=name1)
      rewind(69)
      rewind(15)

      read(15,'(a80)') line

      if(iqchem.eq.0) then

!! standard/MOKIT layout: both Orbital Energies blocks precede both MO  !!
!! coefficient blocks.                                                  !!
        do while(index(line,"Alpha Orbital").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do
        write(69,11) "Alpha Orbital Energies","R","N= ",indepigr
        write(69,13) (energ_a(ii),ii=1,indepigr)

        do while(index(line,"Beta Orbital").eq.0)
          read(15,'(a80)') line
        end do
        write(69,11) "Beta Orbital Energies ","R","N= ",indepigr
        write(69,13) (energ_b(ii),ii=1,indepigr)

        do while(index(line,"Alpha MO co").eq.0)
          read(15,'(a80)') line
        end do
        write(69,12) "Alpha MO coefficients","R","N= ",norb
        write(69,13) ((pcoef_a(ii,jj),ii=1,igr),jj=1,indepigr)

        do while(index(line,"Beta MO coef").eq.0)
          read(15,'(a80)') line
        end do
        write(69,12) "Beta MO coefficients ","R","N= ",norb
        write(69,13) ((pcoef_b(ii,jj),ii=1,igr),jj=1,indepigr)

        if(imokit.eq.0) then
          do while(index(line,"Orthonormal basis").eq.0)
            read(15,'(a80)') line
          end do
        end if
      else

!! Q-Chem layout: both MO coefficient blocks precede both Orbital       !!
!! Energies blocks -- unverified, see header.                           !!
        do while(index(line,"Alpha MO co").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do
        write(69,12) "Alpha MO coefficients","R","N= ",norb
        write(69,13) ((pcoef_a(ii,jj),ii=1,igr),jj=1,indepigr)

        do while(index(line,"Beta MO coef").eq.0)
          read(15,'(a80)') line
        end do
        write(69,12) "Beta MO coefficients ","R","N= ",norb
        write(69,13) ((pcoef_b(ii,jj),ii=1,igr),jj=1,indepigr)

        do while(index(line,"Alpha Orbital").eq.0)
          read(15,'(a80)') line
        end do
        write(69,11) "Alpha Orbital Energies","R","N= ",indepigr
        write(69,13) (energ_a(ii),ii=1,indepigr)

        do while(index(line,"Beta Orbital").eq.0)
          read(15,'(a80)') line
        end do
        write(69,11) "Beta Orbital Energies ","R","N= ",indepigr
        write(69,13) (energ_b(ii),ii=1,indepigr)
      end if

!! copy everything else through unchanged, Total/Spin SCF Density       !!
!! included                                                              !!
      do while(.true.)
        write(69,'(a80)') line
        read(15,'(a80)',end=99) line
      end do
99    continue
      close(69)
      DEALLOCATE(energ_a,energ_b)

!! Printing formats !!
11    FORMAT(a23,20x,a1,3x,a3,i11)
12    FORMAT(a21,22x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))

      end

!! ****** !!
