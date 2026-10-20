!! ----------------------------------------------------------------------
!! GEOS (generalized effective oxidation states) and EFFAO-U: effective
!! fragment orbitals of the paired and unpaired densities, and the
!! electron assignment across both.
!!   effao3d_u     -- EFOs of the paired and unpaired densities
!!   ueos_analysis -- electron assignment, oxidation states, .fchk order
!! ----------------------------------------------------------------------

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: effao3d_u
!! Real-space EFOs of the paired and unpaired densities (unpaired: u =
!! n(2-n) per natural orbital, Takatsuka; paired: total minus unpaired),
!! one fragment at a time, as ueffao3d_frag does for alpha and beta.
!! Paired EFOs with a clearly negative net occupation (below
!! xminocc_neg; the paired density is a difference of densities) are
!! kept apart: reported, written to the .fchk file, left out of the
!! assignment. With GEOS (iueos=1) ueos_analysis then assigns the
!! electrons; with EFFAO-U only the EFOs are printed.
!! arguments:
!!   itotps (in)  number of grid points
!!   ndim   (in)  number of basis functions
!!   omp    (in)  (itotps) weight of each point for its own atom
!!   chp    (in)  (itotps,ndim) basis functions at the points
!!   sat    (in)  (ndim,ndim,nat) atomic overlap matrices
!!   wp     (in)  (itotps) integration weights
!!   omp2   (in)  (itotps,nat) atomic weights of each point
!!   iueos  (in)  1: GEOS (assignment and .fchk); 0: EFFAO-U
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine effao3d_u(itotps,ndim,omp,chp,sat,wp,omp2,iueos)

      use basis_set
      use ao_matrices
      use integration_grid
      use effao_mod, only: p0,p0net,p0gro,ip0
      use input_options_mod, only: inegefos,inegcubthr
      use output_mod, only: print_value
      use json_mod

      implicit real*8(a-h,o-z)
      include 'parameter.h'

      integer,intent(in) :: itotps,ndim

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /iops/iopt(200)
      common /loba2/occup(nmax,2),iorbat(nmax,2),lorb(2),confi0

      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat),wp(itotps)
      dimension sat(ndim,ndim,nat)

      allocatable :: Pno(:,:),Uno(:,:)
      allocatable :: S0(:,:),Sm(:,:),Splus(:,:),c0(:,:),pp0(:,:)
      allocatable :: scr(:),s0all(:)
      allocatable :: iup0(:,:),up0net(:,:,:),up0gro(:,:,:),up0coef(:,:,:,:)
      allocatable :: up0fit(:,:,:)
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort
      allocatable :: poolcoef(:,:,:)
      allocatable :: ineg_frg(:),xneg_net(:,:),xneg_gro(:,:),cneg(:,:,:)
      allocatable :: xneg_fit(:,:)
      allocatable :: Sinv(:,:),cprojfrag(:,:),xfitfrag(:)
      allocatable :: p0sav(:,:),occsav(:,:)
      allocatable :: occn(:),fitp(:)
      character(len=20) :: ctype

      icube = iopt(13)
      iatps = nang*nrad

!! EFO occupation cutoff for the net-population sum below !!
      xminocc=1.0d-4

!! independent threshold for significant negative net occupation, paired !!
!! channel only (icase=1): Pno-Uno is a difference of two PSD matrices,   !!
!! not itself PSD, so genuine negative eigenvalues are possible after the !!
!! congruence transform below. Uno alone (icase=2) is PSD and never       !!
!! triggers this. Hand-tuned for now, expect per-system adjustment.       !!
      xminocc_neg=-0.025d0

      call efo_section_header(0,xminocc)

!! build the total (Pno) and Takatsuka unpaired (Uno) density matrices,  !!
!! in the AO basis, from the natural orbitals already in ao_matrices.    !!
      nnorb=0
      do ii=1,igr
        if(occ_no(ii,ii).ge.thresh) nnorb=nnorb+1
      end do
      call print_value('Natural orbitals',nnorb)

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
      ALLOCATE(c0(igr,igr),pp0(igr,igr),occn(igr),fitp(igr))
      if(iueos.eq.1) ALLOCATE(iup0(2,icufr),up0net(2,igr,icufr),up0gro(2,igr,icufr),
     +  up0coef(igr,igr,2,icufr),poolcoef(igr,igr,2),up0fit(2,igr,icufr))
      if(iueos.eq.1) ALLOCATE(ineg_frg(icufr),xneg_net(igr,icufr),xneg_gro(igr,icufr),
     +  cneg(igr,igr,icufr),xneg_fit(igr,icufr))

!! S^-1 for reproject_efos_ao, built once for all fragments/channels !!
      if(iueos.eq.1) then
        ALLOCATE(Sinv(igr,igr),cprojfrag(igr,igr),xfitfrag(igr))
        call build_Sinv(igr,Sinv)
      end if

!! EFFAO-U writes its EFOs to an "effao" block (GEOS: ueos_analysis) !!
      if(iueos.eq.0) then
        call atomic_definition_names(ckey,cname,cshort,cspace)
        call json_begin_object('effao')
        call json_value('scheme',trim(ckey))
      end if

!! icase=1: paired density (total-unpaired). icase=2: unpaired density. !!
      do icase=1,2
        if(icase.eq.1) then
          call print_box('PAIRED DENSITY')
        else if(icase.eq.2) then
          call print_box('UNPAIRED DENSITY')
        end if

!! per-fragment loop kept serial on purpose: icufr can be small on a    !!
!! large system -- the O(igr^2)/O(igr^3) work inside each iteration is  !!
!! threaded instead, same lesson already applied throughout effao.f.    !!
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

!! keep EFOs above xminocc; pp0's diagonal is already sorted decreasing  !!
!! imaxo starts at 0 so a channel with no EFO above threshold (e.g. the  !!
!! unpaired channel on a restricted wavefunction, where Uno is ~0 for    !!
!! integer NO occupations) correctly ends up with imaxo=0 instead of     !!
!! carrying over a stale value from the previous fragment/icase.         !!
          imaxo=0
          do ii=1,igr
            if(pp0(ii,ii).lt.xminocc) exit
            imaxo=ii
          end do

!! significant-negative-occupation EFOs, paired channel only: pp0's       !!
!! diagonal is sorted decreasing, so these sit as a contiguous block at   !!
!! the tail, past index imaxo.                                           !!
          ineg=0
          if(iueos.eq.1.and.icase.eq.1) then
            jj=igr
            do while(jj.gt.imaxo.and.pp0(jj,jj).le.xminocc_neg)
              ineg=ineg+1
              jj=jj-1
            end do
          end if

!! net occupation of the kept EFOs (xmaxo) and of all of them         !!
!! (xnetall, negative paired EFOs included); Tr(D S^F) of this channel's !!
!! density (xgroall): the differences are what the cutoff leaves out    !!
          xmaxo=ZERO
          do ii=1,imaxo
            xmaxo=xmaxo+pp0(ii,ii)
          end do
          xnetall=ZERO
          do ii=1,igr
            xnetall=xnetall+pp0(ii,ii)
          end do
          xgroall=ZERO
          do icenter=1,nfrlist(iicenter)
            jcenter=ifrlist(icenter,iicenter)
            do jj=1,igr
              do kk=1,igr
                if(icase.eq.1) xx=Pno(kk,jj)-Uno(kk,jj)
                if(icase.eq.2) xx=Uno(kk,jj)
                xgroall=xgroall+xx*sat(jj,kk,jcenter)
              end do
            end do
          end do

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

!! gross occupation of the significant-negative-occupation EFOs found     !!
!! above (paired channel only) -- same sat contraction as the positive    !!
!! EFOs, kept in its own loop since these live at the opposite end of     !!
!! pp0's sorted diagonal and are never mixed into imaxo/s0all/xx0/the     !!
!! oxidation-state assignment (see ueos_analysis).                        !!
          if(ineg.gt.0) then
!! gross occupations stay on the raw c0; only the .fchk copy (cneg) is   !!
!! projected                                                             !!
            call reproject_efos_ao(igr,itotps,chp,wp,omp,scr,Sinv,
     +        c0,igr-ineg+1,igr,cprojfrag,xfitfrag)
!$OMP PARALLEL DO PRIVATE(ii,iorb,icenter,jcenter,jj,kk,xx,xxx)
            do ii=1,ineg
              iorb=igr-ineg+ii
              xxx=ZERO
              do icenter=1,nfrlist(iicenter)
                jcenter=ifrlist(icenter,iicenter)
                xx=ZERO
                do jj=1,igr
                  do kk=1,igr
                    xx=xx+c0(kk,iorb)*sat(kk,jj,jcenter)*c0(jj,iorb)
                  end do
                end do
                xxx=xxx+xx
              end do
              xxx=xxx*pp0(iorb,iorb)
              xneg_gro(ii,iicenter)=xxx
              xneg_net(ii,iicenter)=pp0(iorb,iorb)
              xneg_fit(ii,iicenter)=xfitfrag(ii)
              do kk=1,igr
                cneg(kk,ii,iicenter)=cprojfrag(kk,iorb)
              end do
            end do
!$OMP END PARALLEL DO
          end if
          if(iueos.eq.1.and.icase.eq.1) ineg_frg(iicenter)=ineg

!! store for cube generation (raw c0 -- cubegen_new applies the fragment !!
!! weight itself) and, with GEOS, for ueos_analysis (projected, for the  !!
!! .fchk); Fit (%) = 100*(1-xfit), per EFO                               !!
          if(iueos.eq.1.and.imaxo.gt.0) then
            call reproject_efos_ao(igr,itotps,chp,wp,omp,scr,Sinv,
     +        c0,1,imaxo,cprojfrag,xfitfrag)
            do mu=1,imaxo
              fitp(mu)=100.0d0*(1.0d0-xfitfrag(mu))
            end do
          end if
          do mu=1,imaxo
            occn(mu)=pp0(mu,mu)
          end do
          if(iicenter.gt.1) write(*,'(a)') ''
          call print_efo_block('Fragment',iicenter,imaxo,xmaxo,occn,
     +      .true.,xx0,s0all,iueos.eq.1,fitp,xnetall-xmaxo,xgroall-xx0)
          do kk=1,imaxo
            do mu=1,igr
              p0(mu,kk)=c0(mu,kk)
            end do
            p0net(kk,iicenter)=pp0(kk,kk)
            p0gro(kk,iicenter)=s0all(kk)
            if(iueos.eq.1) then
              up0net(icase,kk,iicenter)=pp0(kk,kk)
              up0gro(icase,kk,iicenter)=s0all(kk)
              up0fit(icase,kk,iicenter)=fitp(kk)
              do mu=1,igr
                up0coef(mu,kk,icase,iicenter)=cprojfrag(mu,kk)
              end do
            end if
          end do
          ip0(iicenter)=imaxo
          if(iueos.eq.1) iup0(icase,iicenter)=imaxo

!! optional: write a cube file for this fragment's paired/unpaired EFOs !!
          if(icase.eq.1) iicase=3
          if(icase.eq.2) iicase=4
          if(icube.eq.1) call cubegen_new(iicenter,iicase)

!! # CUBE NEG_EFOS: paired EFOs with net occ. <= -inegcubthr/1000, own   !!
!! selection independent of xminocc_neg (report/.fchk set), numbered     !!
!! from the most negative. Staged into p0's leading columns from the raw !!
!! c0 (cubegen_new applies the fragment weight itself), borrowed slots   !!
!! restored after.                                                       !!
          inegc=0
          if(icube.eq.1.and.inegefos.eq.1.and.icase.eq.1) then
            xnegthr=-REAL(inegcubthr)/1000.0d0
            jj=igr
            do while(jj.gt.imaxo.and.pp0(jj,jj).le.xnegthr)
              inegc=inegc+1
              jj=jj-1
            end do
          end if
          if(inegc.gt.0) then
            ALLOCATE(p0sav(igr,inegc),occsav(inegc,2))
            do ii=1,inegc
              iorb=igr-ii+1
              do mu=1,igr
                p0sav(mu,ii)=p0(mu,ii)
                p0(mu,ii)=c0(mu,iorb)
              end do
              occsav(ii,1)=p0net(ii,iicenter)
              occsav(ii,2)=p0gro(ii,iicenter)
!! gross occupation: same sat contraction as the report's xneg_gro       !!
              xxx=ZERO
              do icenter=1,nfrlist(iicenter)
                jcenter=ifrlist(icenter,iicenter)
                do jj=1,igr
                  do kk=1,igr
                    xxx=xxx+c0(kk,iorb)*sat(kk,jj,jcenter)*c0(jj,iorb)
                  end do
                end do
              end do
              p0net(ii,iicenter)=pp0(iorb,iorb)
              p0gro(ii,iicenter)=xxx*pp0(iorb,iorb)
            end do
            ip0(iicenter)=inegc
            call cubegen_new(iicenter,5)
            do ii=1,inegc
              do mu=1,igr
                p0(mu,ii)=p0sav(mu,ii)
              end do
              p0net(ii,iicenter)=occsav(ii,1)
              p0gro(ii,iicenter)=occsav(ii,2)
            end do
            ip0(iicenter)=imaxo
            DEALLOCATE(p0sav,occsav)
          end if
        end do
        if(iueos.eq.0) call json_efo_channel(icase+2,.true.,.false.)
      end do
      if(iueos.eq.0) call json_end()

!! electron/oxidation-state assignment across both channels together,   !!
!! using gross populations                                              !!
      if(iueos.eq.1) then
        call ueos_analysis(iup0,up0gro,up0coef,poolcoef,
     +    ineg_frg,xneg_net,xneg_gro,cneg,xminocc_neg,xneg_fit,
     +    up0net,up0fit)

!! pooled paired/unpaired EFOs as fake Alpha/Beta MOs in a .fchk, for    !!
!! visualization in any standard viewer -- printed by default, same as  !!
!! OSLO's own .fchk output, no separate keyword needed. A restricted    !!
!! wavefunction (kop=0) has no Beta blocks to splice into, but its      !!
!! unpaired channel isn't always trivially empty -- e.g. a restricted-  !!
!! orbital post-HF wavefunction (CASSCF/FCI/CISD, kop=0 since pySCF's   !!
!! writer never emits a Beta MO set for these regardless of spin state) !!
!! can have real fractional natural-orbital occupations and hence a     !!
!! genuinely populated unpaired channel. When that happens (lorb(2)>0), !!
!! rwfu_effao_orbprint inserts a synthetic Beta MO set (unpaired) next  !!
!! to the real Alpha one (paired) so both channels still end up in one  !!
!! file -- otherwise fall back to the plain restricted, single-channel  !!
!! writer, same as for an ordinary restricted-wavefunction EFFAO run.   !!
        ctype="-GEOS-EFOs"
        if(kop.eq.0) then
          if(lorb(2).gt.0) then
            call rwfu_effao_orbprint(poolcoef(:,:,1),poolcoef(:,:,2),ctype)
          else
            call rwf_effao_orbprint(poolcoef(:,:,1),ctype)
          end if
        else
          call uwf_effao_orbprint(poolcoef(:,:,1),poolcoef(:,:,2),ctype)
        end if
      end if

      DEALLOCATE(S0,Splus,Sm)
      DEALLOCATE(scr,c0,pp0,s0all,occn,fitp)
      DEALLOCATE(Pno,Uno)
      if(iueos.eq.1) then
        DEALLOCATE(iup0,up0net,up0gro,up0coef,poolcoef,up0fit)
        DEALLOCATE(ineg_frg,xneg_net,xneg_gro,cneg,xneg_fit)
        DEALLOCATE(Sinv,cprojfrag,xfitfrag)
      end if

      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: ueos_analysis
!! GEOS electron assignment: pools the paired and unpaired EFOs of all
!! fragments, sorted by gross occupation; starts from every electron
!! paired (except the |Na-Nb| unpaired ones) and moves one pair at a
!! time to the two most occupied unpaired EFOs while that lowers the
!! RMSD between occupations and assigned electrons; prints both
!! channels, the oxidation states and R, and orders the EFOs (negative
!! paired ones last) for the .fchk file.
!! arguments:
!!   iup0     (in)   (2,icufr) EFOs kept per channel and fragment
!!   up0gro   (in)   (2,igr,icufr) their gross occupations
!!   up0coef  (in)   (igr,igr,2,icufr) their coefficients (reprojected)
!!   poolcoef (out)  (igr,igr,2) pooled coefficients, .fchk order
!!   ineg_frg (in)   (icufr) negative paired EFOs of each fragment
!!   xneg_net, xneg_gro (in)  (igr,icufr) their net, gross occupations
!!   cneg     (in)   (igr,igr,icufr) their coefficients
!!   xminocc_neg (in)  net occupation below which a paired EFO is negative
!!   xneg_fit (in)   (igr,icufr) fit fraction of each negative EFO
!!   up0net   (in)   (2,igr,icufr) net occupations of the kept EFOs
!!   up0fit   (in)   (2,igr,icufr) fit (%) of the kept EFOs
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine ueos_analysis(iup0,up0gro,up0coef,poolcoef,
     +  ineg_frg,xneg_net,xneg_gro,cneg,xminocc_neg,xneg_fit,
     +  up0net,up0fit)

      use output_mod, only: print_value
      use json_mod
      implicit real*8(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /filename/name0
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /loba/ oxi(maxat),errsav(maxat),elec(maxat),effpop(maxat)
      common /loba2/occup(nmax,2),iorbat(nmax,2),lorb(2),confi0

      dimension iup0(2,icufr),up0gro(2,igr,icufr)
      dimension up0coef(igr,igr,2,icufr),poolcoef(igr,igr,2)
      dimension iorbslot(nmax,2)
      dimension ineg_frg(icufr),xneg_net(igr,icufr),xneg_gro(igr,icufr)
      dimension cneg(igr,igr,icufr),xneg_fit(igr,icufr)
      dimension up0net(2,igr,icufr),up0fit(2,igr,icufr)
      character*60 name0
      character*80 cfile
      character*8 cspin
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort
      integer, allocatable :: ifr(:),icol(:)
      real*8, allocatable :: xn(:),xf(:),xg(:),xe(:)
      logical, allocatable :: lcol(:),lok(:)

      allocatable :: tmp_occup(:)
      integer, allocatable :: tmp_iorbat(:),tmp_iorbslot(:)
      allocatable :: elec2(:,:),elec_id(:,:)
      allocatable :: tmp_elec_id(:,:)
      allocatable :: elec_frg_count(:,:)
      dimension xlo(maxfrag),xfu(maxfrag),ilo(maxfrag),ifu(maxfrag)
      dimension nlo(maxfrag),ncount(maxfrag)
      dimension npool(2)
      allocatable :: xnegnet_pool(:),xneggro_pool(:),neg_frg(:),neg_slot(:)
      character*60 cval

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
      call atomic_definition_names(ckey,cname,cshort,cspace)
      call json_begin_object('geos')
      call json_value('scheme',trim(ckey))
      call print_box('GEOS ELECTRON ASSIGNMENT')
      call print_value('EFOs in the assignment',lorb(1)+lorb(2))

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
!! occupations equal within 1.0d-8 (symmetry-equivalent fragments, whose !!
!! grid integrals differ in the last digits) are a tie and are not       !!
!! swapped, so the order is the same on every machine                    !!
            if(tmp_occup(ii).lt.tmp_occup(jj)-1.0d-8) then
              tmp_swap=tmp_occup(ii)
              tmp_occup(ii)=tmp_occup(jj)
              tmp_occup(jj)=tmp_swap

              itmp_swap=tmp_iorbat(ii)
              tmp_iorbat(ii)=tmp_iorbat(jj)
              tmp_iorbat(jj)=itmp_swap

              itmp_swap=tmp_iorbslot(ii)
              tmp_iorbslot(ii)=tmp_iorbslot(jj)
              tmp_iorbslot(jj)=itmp_swap
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

!! pool every fragment's significant-negative-occupation paired EFOs      !!
!! (effao3d_u) into one flat list, then sort by decreasing net occupation !!
!! (least-negative first, most-negative last) -- independent of the       !!
!! positive-channel pooling above, and never fed into occup/iorbat/lorb   !!
!! or the electron-assignment logic below.                                !!
      nneg=0
      do ifrg=1,icufr
        nneg=nneg+ineg_frg(ifrg)
      end do

      if(nneg.gt.0) then
        ALLOCATE(xnegnet_pool(nneg),xneggro_pool(nneg),neg_frg(nneg),neg_slot(nneg))
        iorb=0
        do ifrg=1,icufr
          do kk=1,ineg_frg(ifrg)
            iorb=iorb+1
            xnegnet_pool(iorb)=xneg_net(kk,ifrg)
            xneggro_pool(iorb)=xneg_gro(kk,ifrg)
            neg_frg(iorb)=ifrg
            neg_slot(iorb)=kk
          end do
        end do

        do ii=1,nneg-1
          do jj=ii+1,nneg
            if(xnegnet_pool(ii).lt.xnegnet_pool(jj)-1.0d-8) then
              tmp_swap=xnegnet_pool(ii)
              xnegnet_pool(ii)=xnegnet_pool(jj)
              xnegnet_pool(jj)=tmp_swap

              tmp_swap=xneggro_pool(ii)
              xneggro_pool(ii)=xneggro_pool(jj)
              xneggro_pool(jj)=tmp_swap

              itmp_swap=neg_frg(ii)
              neg_frg(ii)=neg_frg(jj)
              neg_frg(jj)=itmp_swap

              itmp_swap=neg_slot(ii)
              neg_slot(ii)=neg_slot(jj)
              neg_slot(jj)=itmp_swap
            end if
          end do
        end do

        write(*,'(2x,a,i0,a,f7.4,a)') 'WARNING: ',nneg,' paired '//
     +    'EFO(s) with net occupation below ',xminocc_neg,
     +    ' left out of the assignment'
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Frag.  Net occ.  Gross occ.        Fit (%)'
        write(*,'(2x,a)') repeat('-',42)
        do ii=1,nneg
          write(*,'(2x,i5,f10.4,f12.4,f15.2)') neg_frg(ii),
     +      xnegnet_pool(ii),xneggro_pool(ii),
     +      100.0d0*(1.0d0-xneg_fit(neg_slot(ii),neg_frg(ii)))
        end do
        write(*,'(2x,a)') repeat('-',42)
      end if

!! pooled coefficient matrix per channel, in the sorted order above,    !!
!! zero-padded. The .fchk file holds one orbital per independent basis  !!
!! function (nindep, fewer than igr with linear dependencies), so the   !!
!! paired channel keeps the last nneg of those slots for the negative   !!
!! EFOs pooled above, after the positive ones that fit before them.     !!
      nindep=nindep_fchk(igr)
      poolcoef=ZERO
      do icase=1,2
        nwin=nindep
        if(icase.eq.1) nwin=nindep-nneg
        do jj=1,MIN(lorb(icase),nwin)
          ii=iorbat(jj,icase)
          kk=iorbslot(jj,icase)
          do mu=1,igr
            poolcoef(mu,jj,icase)=up0coef(mu,kk,icase,ii)
          end do
        end do
      end do

      if(nneg.gt.0) then
        do jj=1,nneg
          ii=neg_frg(jj)
          kk=neg_slot(jj)
          do mu=1,igr
            poolcoef(mu,nindep-nneg+jj,1)=cneg(mu,kk,ii)
          end do
        end do
      end if

!! ideal-occupation matrix, sized to the larger of the two EFO counts !!
      ilorb=MAX(lorb(1),lorb(2))
      ALLOCATE(elec_id(ilorb,2))
      ALLOCATE(elec2(ilorb,2))
      elec_id=ZERO

!! section, not whole-array -- keeps elec2 at ilorb (a bare "=occup"    !!
!! would auto-reallocate to occup's own (nmax,2) shape, F2003 semantics) !!
      elec2=occup(1:ilorb,1:2)

!! patch the paired channel's tail slots with the negative EFOs' gross    !!
!! occupation, matching where poolcoef placed their coefficients above -- !!
!! purely cosmetic, feeds the .fchk writers' fake orbital energies        !!
!! (print.f, occup(:,1)) with zero changes needed there. Done only after  !!
!! elec2 has already captured the real (unpatched) values above, so the   !!
!! electron-assignment logic below never sees these EFOs.                 !!
      if(nneg.gt.0) then
        do jj=lorb(1)+1,nindep-nneg
          occup(jj,1)=ZERO
        end do
        do jj=1,nneg
          occup(nindep-nneg+jj,1)=xneggro_pool(jj)
        end do
      end if

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
      write(*,'(a)') ''
      write(cval,'(i0,a,i0,a)') npair,' electron pairs, ',nunp,
     +  ' unpaired electrons'
      call print_value('Ideal assignment',trim(cval))

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
      call print_value('RMSD of the ideal assignment',xrmsd,5)

!! iteratively move the least-occupied paired electron pair to the two  !!
!! most-occupied unpaired slots, keeping the move only while it lowers  !!
!! the RMSD. code stops at the first move that doesn't help.            !!
      ALLOCATE(tmp_elec_id(ilorb,2))
      tmp_elec_id=elec_id
      if(npair.gt.0) then
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Pairs moved to the unpaired density '//
     +    '(kept while the RMSD decreases):'
        write(*,'(2x,a)') 'Try  Frag.  Paired  Frag.  Unpaired  '//
     +    'Frag.  Unpaired     RMSD  Result'
        write(*,'(2x,a)') repeat('-',71)
      end if
      itry=0
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

        itry=itry+1
        if((xrmsd2-xrmsd).lt.ZERO) then
          write(*,100) itry,iorbat(npair,1),elec2(npair,1),
     +      iorbat(nunp+1,2),elec2(nunp+1,2),iorbat(nunp+2,2),
     +      elec2(nunp+2,2),xrmsd2,'accepted'
          elec_id=tmp_elec_id
          xrmsd=xrmsd2
        else
          write(*,100) itry,iorbat(npair,1),elec2(npair,1),
     +      iorbat(nunp+1,2),elec2(nunp+1,2),iorbat(nunp+2,2),
     +      elec2(nunp+2,2),xrmsd2,'rejected'
          go to 69
        end if

        nunp=nunp+2
        npair=npair-1
      end do
69    continue
      if(itry.gt.0) write(*,'(2x,a)') repeat('-',71)
      call print_value('RMSD of the electron assignment',xrmsd,5)
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

!! the paired EFO count goes up to nindep so that the .fchk writers    !!
!! (print.f, called after this returns) reach the negative EFOs placed  !!
!! at the end; every use of lorb(1) in the assignment is done by now    !!
      npool(1)=lorb(1)
      npool(2)=lorb(2)
      if(nneg.gt.0) lorb(1)=MAX(lorb(1),nindep)

      do icase=1,2
        if(icase.eq.1) then
          call print_box('GEOS, PAIRED ELECTRONS')
        else if(icase.eq.2) then
          call print_box('GEOS, UNPAIRED ELECTRONS')
        end if

!! last occupied / first unoccupied EFO gross occupation per fragment, !!
!! as in effao.f's eos_analysis: with nn its occupied EFOs (pairs in    !!
!! the paired density) rounded up, the nn-th (xlo) and (nn+1)-th (xfu)  !!
!! of its EFOs in gross order, as met in the sorted pooled list         !!
        do ifrg=1,icufr
          if(icase.eq.1) then
            nn=INT(elec_frg_count(ifrg,icase))/2
            if(elec_frg_count(ifrg,icase)/TWO-nn.gt.thresh) nn=nn+1
          end if
          if(icase.eq.2) then
            nn=INT(elec_frg_count(ifrg,icase))
            if(elec_frg_count(ifrg,icase)-nn.gt.thresh) nn=nn+1
          end if
          nlo(ifrg)=nn
          ncount(ifrg)=0
          ilo(ifrg)=0
          ifu(ifrg)=0
          xlo(ifrg)=ZERO
          xfu(ifrg)=ZERO
        end do
        do ii=1,npool(icase)
          ifrg=iorbat(ii,icase)
          ncount(ifrg)=ncount(ifrg)+1
          if(ncount(ifrg).eq.nlo(ifrg)) then
            xlo(ifrg)=elec2(ii,icase)
            ilo(ifrg)=1
          else if(ncount(ifrg).eq.nlo(ifrg)+1) then
            xfu(ifrg)=elec2(ii,icase)
            ifu(ifrg)=1
          end if
        end do

        if(icase.eq.1) xlast=TWO
        if(icase.eq.2) xlast=ONE
        ilast=0
        do ifrg=1,icufr
          if(ilo(ifrg).eq.1) then
            if(xlo(ifrg).lt.xlast) then
              xlast=xlo(ifrg)
              ilast=ifrg
            end if
          end if
        end do

        xfirst=ZERO
        do ifrg=1,icufr
          if(ifrg.ne.ilast.and.ifu(ifrg).eq.1) then
            if(xfu(ifrg).gt.xfirst) xfirst=xfu(ifrg)
          end if
        end do

!! reliability index: paired occupations run 0-2 instead of 0-1, so the !!
!! gap is halved before the same +0.5 scaling used for unpaired/EOS.    !!
        if(icase.eq.1) confi=100.0*dmin1(1.0d0,(xlast-xfirst)/TWO+0.5d0)
        if(icase.eq.2) confi=100.0*dmin1(1.0d0,xlast-xfirst+0.5d0)
        call print_eos_frontier(icufr,elec_frg_count(1,icase),xlo,ilo,
     +    xfu,ifu,confi)
        if(icase.eq.1) then
          write(*,'(2x,a)') 'Paired occupations (0 to 2) are halved in R'
          confi0=confi
        end if

!! results file: this channel's R, frontier table and pooled EFOs, in  !!
!! the order of the .fchk orbitals of the channel (paired: Alpha,      !!
!! unpaired: Beta), with the electrons each one gets                   !!
        cspin='paired'
        if(icase.eq.2) cspin='unpaired'
        call json_begin_object(trim(cspin))
        call json_value('reliability_percent',confi)
        call json_index('reliability_percent','scalar','Reliability '//
     +    'index R, '//trim(cspin)//' ('//trim(cshort)//')','%',3,
     +    trim(cspin),.false.,' ')
        allocate(lok(icufr))
        call json_begin_object('frontier')
        call json_value('rows','fragment')
        call json_value('n_rows',icufr)
        call json_begin_array('columns')
        call json_column('last_occupied','Last occupied EFO',' ',3,
     +    .false.)
        call json_column('first_unoccupied','First unoccupied EFO',' ',
     +    3,.false.)
        call json_end()
        call json_begin_object('data')
        do ifrg=1,icufr
          lok(ifrg)=ilo(ifrg).eq.1
        end do
        call json_vector_or_null('last_occupied',xlo(1:icufr),lok,icufr)
        do ifrg=1,icufr
          lok(ifrg)=ifu(ifrg).eq.1
        end do
        call json_vector_or_null('first_unoccupied',xfu(1:icufr),lok,
     +    icufr)
        call json_end()
        call json_end()
        call json_index('frontier','table','GEOS frontier EFOs, '//
     +    trim(cspin)//' ('//trim(cshort)//')',' ',-1,trim(cspin),
     +    .false.,' ')
        deallocate(lok)
        n=npool(icase)
        nwin=nindep
        if(icase.eq.1) nwin=nindep-nneg
        allocate(ifr(n),icol(n),xn(n),xf(n),xg(n),xe(n),lcol(n))
        do ii=1,n
          ifr(ii)=iorbat(ii,icase)
          xn(ii)=up0net(icase,iorbslot(ii,icase),ifr(ii))
          xf(ii)=up0fit(icase,iorbslot(ii,icase),ifr(ii))
          xg(ii)=elec2(ii,icase)
          xe(ii)=elec_id(ii,icase)
          icol(ii)=ii
          lcol(ii)=ii.le.nwin
        end do
        if(icase.eq.1) then
          call json_efo_table(n,ifr,xn,.true.,xg,.true.,xf,.true.,xe,
     +      icol,lcol,trim(cspin),'alpha')
        else
          call json_efo_table(n,ifr,xn,.true.,xg,.true.,xf,.true.,xe,
     +      icol,lcol,trim(cspin),'beta')
        end if
        deallocate(ifr,icol,xn,xf,xg,xe,lcol)
        call json_end()
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

      confi2=dmin1(confi,confi0)
      call print_oxidation_states(icufr,oxi,zztot,.true.,confi2)

!! results file: the negative paired EFOs (last orbitals of the Alpha  !!
!! block of the .fchk), per-fragment results, overall R, RMSD, file    !!
      allocate(ifr(nneg),icol(nneg),xn(nneg),xg(nneg),xf(nneg))
      do jj=1,nneg
        ifr(jj)=neg_frg(jj)
        xn(jj)=xnegnet_pool(jj)
        xg(jj)=xneggro_pool(jj)
        xf(jj)=100.0d0*(1.0d0-xneg_fit(neg_slot(jj),neg_frg(jj)))
        icol(jj)=nindep-nneg+jj
      end do
      call json_begin_object('negative_efos')
      call json_value('rows','orbital')
      call json_value('n_rows',nneg)
      call json_value('fchk_orbitals','alpha')
      call json_begin_array('columns')
      call json_column('fragments','Fragments',' ',0,.false.)
      call json_column('shares','Fragment shares',' ',5,.false.)
      call json_column('net','Net occupation',' ',4,.true.)
      call json_column('gross','Gross occupation',' ',4,.true.)
      call json_column('fit_percent','Fit','%',2,.false.)
      call json_column('fchk_column','Column in the .fchk file',' ',0,
     +  .false.)
      call json_end()
      call json_begin_object('data')
      call json_begin_array('fragments')
      do jj=1,nneg
        call json_vector('',ifr(jj:jj),1)
      end do
      call json_end()
      call json_begin_array('shares')
      do jj=1,nneg
        call json_vector('',(/ONE/),1)
      end do
      call json_end()
      call json_vector('net',xn,nneg)
      call json_vector('gross',xg,nneg)
      call json_vector('fit_percent',xf,nneg)
      call json_vector('fchk_column',icol,nneg)
      call json_end()
      call json_end()
      call json_index('negative_efos','table','GEOS negative paired '//
     +  'EFOs ('//trim(cshort)//')',' ',-1,'paired',.false.,' ')
      deallocate(ifr,icol,xn,xg,xf)
      call json_begin_object('fragments')
      call json_vector('oxidation_states',oxi(1:icufr),icufr)
      call json_index('oxidation_states','fragment_vector',
     +  'GEOS oxidation states ('//trim(cshort)//')',' ',2,'total',
     +  .true.,' ')
      call json_vector('electrons_paired',elec_frg_count(1:icufr,1),
     +  icufr)
      call json_index('electrons_paired','fragment_vector',
     +  'GEOS electrons, paired ('//trim(cshort)//')',' ',2,'paired',
     +  .false.,' ')
      call json_vector('electrons_unpaired',elec_frg_count(1:icufr,2),
     +  icufr)
      call json_index('electrons_unpaired','fragment_vector',
     +  'GEOS electrons, unpaired ('//trim(cshort)//')',' ',2,
     +  'unpaired',.false.,' ')
      call json_end()
      call json_value('reliability_percent',confi2)
      call json_index('reliability_percent','scalar','Reliability '//
     +  'index R ('//trim(cshort)//')','%',3,'total',.false.,' ')
      call json_value('rmsd',xrmsd)
      call json_index('rmsd','scalar','RMSD of the GEOS electron '//
     +  'assignment ('//trim(cshort)//')',' ',5,' ',.false.,' ')
      cfile=trim(name0)//'-GEOS-EFOs.fchk'
      cfile=cfile(index(cfile,'/',back=.true.)+1:)
      call json_value('orbital_file',trim(cfile))
      call json_value('orbital_file_id',json_run_id())
      call json_end()

      DEALLOCATE(elec_id,elec2,elec_frg_count)
      if(nneg.gt.0) DEALLOCATE(xnegnet_pool,xneggro_pool,neg_frg,neg_slot)

100   FORMAT(2x,i3,i7,f8.4,i7,f10.4,i7,f10.4,f9.5,2x,a)

      end

!! ***** !!