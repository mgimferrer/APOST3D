!! ----------------------------------------------------------------------
!! Effective fragment orbitals (EFOs) and effective oxidation states
!! (EOS), in real space and in Hilbert space:
!!   efo_section_header -- opens the section (definition, densities,
!!                         occupation cutoff)
!!   ueffao3d_frag      -- real-space EFOs of each fragment
!!   ueffaolow_frag     -- Lowdin, Davidson-Lowdin or NAO EFOs
!!   ueffaomull_frag    -- Mulliken EFOs
!!   ueffao3d           -- real-space EFOs of single atoms (DOATOMS)
!!   eos_analysis       -- electron assignment and oxidation states
!!   efo_json_open, json_efo_table, json_efo_channel -- the results file
!! ----------------------------------------------------------------------

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: efo_section_header
!! Opens the EFFAO/EOS/GEOS section: the atomic definition, the
!! densities the EFOs are computed from, and the occupation cutoff.
!! arguments:
!!   idobeta (in)  1 when the beta EFOs are computed separately
!!   xcut    (in)  EFO occupation cutoff
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine efo_section_header(idobeta,xcut)
      use input_options_mod, only: ieffao,idoat
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort

      if(idoat.ne.0) then
        call print_box('EFFECTIVE ATOMIC ORBITALS (EFAOs)')
      else
        call print_box('EFFECTIVE FRAGMENT ORBITALS (EFOs)')
      end if
      call atomic_definition_names(ckey,cname,cshort,cspace)
      if(cspace.eq.'hilbert') then
        call print_value('Atomic definition',
     +    trim(ckey)//' (Hilbert space)')
      else
        call print_value('Atomic definition',trim(ckey)//' (real space)')
      end if
      if(ieffao.eq.1) then
        call print_value('Densities','total')
      else if(ieffao.eq.3) then
        call print_value('Densities','paired and unpaired')
      else if(idobeta.eq.1.or.idoat.ne.0) then
        call print_value('Densities','alpha and beta separately')
      else
        call print_value('Densities','alpha (beta equal to alpha)')
      end if
      if(idoat.ne.0) then
        call print_value('EFAO occupation cutoff',xcut,5)
      else
        call print_value('EFO occupation cutoff',xcut,5)
      end if
      if(ieffao.eq.1) write(*,'(a)') ''
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: efo_json_open
!! Opens the block of the EFOs in the results file: "eos" with EOS,
!! "effao" otherwise, with the atomic definition. main.f closes it once
!! every channel is written.
!! arguments: none
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine efo_json_open()
      use input_options_mod, only: ieos
      use json_mod
      implicit real*8(a-h,o-z)
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort

      call atomic_definition_names(ckey,cname,cshort,cspace)
      if(ieos.eq.1) then
        call json_begin_object('eos')
      else
        call json_begin_object('effao')
      end if
      call json_value('scheme',trim(ckey))
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: json_efo_table
!! Writes an "efos" table, one row per EFO, and its index entry: its
!! fragment (as the one-element lists fragments and shares, the form
!! every orbital table uses), net occupation, and when given the gross
!! occupation, the fit of the reprojected EFO, the electrons the EOS
!! assignment gives it and its orbital in the EFO .fchk file (null when
!! the file has none), whose block fchk_orbitals names.
!! arguments:
!!   n      (in)  number of EFOs
!!   ifr    (in)  (n) fragment of each EFO
!!   xn     (in)  (n) net occupations
!!   lgross (in)  logical: gross occupations given in xg
!!   xg     (in)  (n) gross occupations
!!   lfit   (in)  logical: fits given in xf
!!   xf     (in)  (n) fits, in %
!!   lelec  (in)  logical: assigned electrons in xe and .fchk columns
!!   xe     (in)  (n) assigned electrons
!!   icol   (in)  (n) .fchk column of each EFO
!!   lcol   (in)  (n) logical: the EFO is in the .fchk file
!!   cspin  (in)  channel, for the index
!!   cfchk  (in)  'alpha' or 'beta': the block of the .fchk file that
!!                holds these EFOs; ' ' when there is no file
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_efo_table(n,ifr,xn,lgross,xg,lfit,xf,lelec,xe,
     +  icol,lcol,cspin,cfchk)
      use json_mod
      implicit real*8(a-h,o-z)
      dimension ifr(*),xn(*),xg(*),xf(*),xe(*),icol(*)
      logical lgross,lfit,lelec,lcol(*)
      character*(*) cspin,cfchk
      dimension xunit(1)
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort

      call atomic_definition_names(ckey,cname,cshort,cspace)
      xunit(1)=1.0d0
      call json_begin_object('efos')
      call json_value('rows','orbital')
      call json_value('n_rows',n)
      if(len_trim(cfchk).gt.0) call json_value('fchk_orbitals',
     +  trim(cfchk))
      call json_begin_array('columns')
      call json_column('fragments','Fragments',' ',0,.false.)
      call json_column('shares','Fragment shares',' ',5,.false.)
      call json_column('net','Net occupation',' ',4,.false.)
      if(lgross) call json_column('gross','Gross occupation',' ',4,
     +  .false.)
      if(lfit) call json_column('fit_percent','Fit','%',2,.false.)
      if(lelec) then
        call json_column('electrons','Assigned electrons',' ',2,.false.)
        call json_column('fchk_column','Column in the .fchk file',' ',
     +    0,.false.)
      end if
      call json_end()
      call json_begin_object('data')
      call json_begin_array('fragments')
      do k=1,n
        call json_vector('',ifr(k:k),1)
      end do
      call json_end()
      call json_begin_array('shares')
      do k=1,n
        call json_vector('',xunit,1)
      end do
      call json_end()
      call json_vector('net',xn(1:n),n)
      if(lgross) call json_vector('gross',xg(1:n),n)
      if(lfit) call json_vector('fit_percent',xf(1:n),n)
      if(lelec) then
        call json_vector('electrons',xe(1:n),n)
        call json_vector_or_null('fchk_column',icol(1:n),lcol(1:n),n)
      end if
      call json_end()
      call json_end()
      call json_index('efos','table','EFOs, '//trim(cspin)//' ('//
     +  trim(cshort)//')',' ',-1,cspin,.false.,' ')
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: json_efo_channel
!! Writes one channel of the "effao" block (EFFAO without EOS, or
!! EFFAO-U): the EFOs of every fragment, in fragment order, from
!! effao_mod.
!! arguments:
!!   icase  (in)  0 total density, 1 alpha, 2 beta, 3 paired, 4 unpaired
!!   lgross (in)  logical: the scheme has gross occupations
!!   lfit   (in)  logical: the scheme has fits (real space)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_efo_channel(icase,lgross,lfit)
      use effao_mod, only: p0net,p0gro,p0fit,ip0
      use json_mod
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,
     +  jfrlist(maxat)
      logical lgross,lfit
      character*8 cspin
      integer, allocatable :: ifr(:),icol(:)
      real*8, allocatable :: xn(:),xg(:),xf(:)
      logical, allocatable :: lcol(:)

      n=0
      do i=1,icufr
        n=n+ip0(i)
      end do
      allocate(ifr(n),xn(n),xg(n),xf(n),icol(n),lcol(n))
      k=0
      do i=1,icufr
        do j=1,ip0(i)
          k=k+1
          ifr(k)=i
          xn(k)=p0net(j,i)
          xg(k)=p0gro(j,i)
          xf(k)=p0fit(j,i)
        end do
      end do
      if(icase.eq.0) then
        cspin='total'
      else if(icase.eq.1) then
        cspin='alpha'
      else if(icase.eq.2) then
        cspin='beta'
      else if(icase.eq.3) then
        cspin='paired'
      else
        cspin='unpaired'
      end if
      call json_begin_object(trim(cspin))
      call json_efo_table(n,ifr,xn,lgross,xg,lfit,xf,.false.,xf,icol,
     +  lcol,trim(cspin),' ')
      call json_end()
      deallocate(ifr,xn,xg,xf,icol,lcol)
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: ueffao3d_frag
!! Real-space effective fragment orbitals: for each fragment, the
!! density pk in the basis orthogonalized with the fragment-weighted
!! overlap, diagonalized; prints the net and gross occupations and the
!! fit of each EFO and stores them in effao_mod (p0, p0net, p0gro, ip0,
!! and p0coef, the EFOs reprojected for the .fchk file).
!! arguments:
!!   itotps (in)  number of grid points
!!   ndim   (in)  number of basis functions
!!   omp    (in)  (itotps) weight of each point for its own atom
!!   chp    (in)  (itotps,ndim) basis functions at the points
!!   sat    (in)  (ndim,ndim,nat) atomic overlap matrices
!!   wp     (in)  (itotps) integration weights
!!   omp2   (in)  (itotps,nat) atomic weights of each point
!!   pk     (in)  (ndim,ndim) density: P, Pa or Pb
!!   icase  (in)  0 total density, 1 alpha, 2 beta
!! author: PSalse, ERaco, MGimf
!! ----------------------------------------------------------------------
      subroutine ueffao3d_frag(itotps,ndim,omp,chp,sat,wp,omp2,pk,icase)

      use integration_grid
      use effao_mod, only: p0,p0net,p0gro,ip0,p0coef,p0fit
      use input_options_mod, only: ieos

      implicit real*8(a-h,o-z)

      include 'parameter.h'

      integer,intent(in) :: itotps,ndim

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(200)
      common /atlist/iatlist(maxat),icuat
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat)
      dimension wp(itotps)
      dimension sat(ndim,ndim,nat)
      dimension pk(ndim,ndim)

      allocatable :: s0(:,:),sm(:,:),c0(:,:),splus(:,:),pp0(:,:)
      allocatable :: s0all(:)
      allocatable :: scr(:)
      allocatable :: Sinv(:,:),cprojfrag(:,:),xfitfrag(:)
      allocatable :: occn(:),fitp(:)

      icube   = Iopt(13)
      ieffthr = Iopt(24)
      iatps   = nang*nrad

      xminocc=REAL(ieffthr)/1000.0d0

      if(icase.eq.1) then
        call print_box('ALPHA DENSITY')
      else if(icase.eq.2) then
        call print_box('BETA DENSITY')
      end if

      ALLOCATE(scr(iatps*nat),occn(ndim),fitp(ndim))
      ALLOCATE(s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      ALLOCATE(c0(ndim,ndim),pp0(ndim,ndim))

!! S^-1 for reproject_efos_ao, built once for all fragments !!
      ALLOCATE(Sinv(igr,igr),cprojfrag(igr,igr),xfitfrag(igr))
      call build_Sinv(igr,Sinv)

!! per-fragment loop kept serial on purpose: icufr can be small on a      !!
!! large system, so the loops INSIDE each iteration are threaded instead. !!
      do iicenter=1,icufr

!! scr(ifut): fragment's total becke/tfvc weight at each grid point.    !!
        do ifut=1,iatps*nat
          scr(ifut)=ZERO
          do icenter=1,nfrlist(iicenter)
            scr(ifut)=scr(ifut)+ omp2(ifut,ifrlist(icenter,iicenter))
          end do
        end do

!! net AO overlap block (s0) of the fragment, integrated over every   !!
!! grid point of the molecule                                        !!

!! parallel over mu: each mu writes only its own s0(mu,*)/s0(*,mu), no  !!
!! two mu iterations collide. nu nested serial (bounds depend on mu).   !!
!! wp/chp/scr/omp shared read-only, x private. dominant cost here.      !!
!$OMP PARALLEL DO PRIVATE(mu,nu,jcenter,ifut,x)
        do mu=1,ndim
          do nu=1,mu
            x=ZERO
            do jcenter=1,nat
              do ifut=iatps*(jcenter-1)+1,iatps*jcenter
                x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*scr(ifut)*scr(ifut)*omp(ifut)
              end do
            end do
            s0(mu,nu)=x
            s0(nu,mu)=x
          end do
        end do
!$OMP END PARALLEL DO

!! S^-1/2 and S^+1/2 of the fragment overlap block !!
        call build_Smp(igr,S0,Sm,Splus,0)

!! transform pk into the orthogonalized fragment basis, diagonalize to  !!
!! get the EFOs (c0) and net occupations (pp0 diagonal), back to AO.    !!
        do i=1,igr
          do j=1,igr
            pp0(i,j)=pk(i,j)
          end do
        end do
        call to_lowdin_basis(igr,Splus,pp0)
        call diagonalize(igr,igr,pp0,C0,0)
        call to_AO_basis(igr,igr,Sm,C0)

!! keep EFOs above the EFF_THRESH cutoff (xminocc); pp0's diagonal is   !!
!! already sorted decreasing, so this is a simple prefix. imaxo starts   !!
!! at 0: a fragment with no EFO above the cutoff keeps none.             !!
        imaxo=0
        do i=1,igr
          if(pp0(i,i).lt.xminocc) exit
          imaxo=i
        end do
!! net occupation of the kept EFOs (xmaxo) and of all of them (xnetall, !!
!! the fragment's net population): the difference is what the cutoff     !!
!! leaves out                                                           !!
        xmaxo=ZERO
        do i=1,imaxo
          xmaxo=xmaxo+pp0(i,i)
        end do
        xnetall=ZERO
        do i=1,igr
          xnetall=xnetall+pp0(i,i)
        end do

!! the fragment's gross population in this density, Tr(P S^F), for the  !!
!! gross part of what the cutoff leaves out                             !!
        xgroall=ZERO
        do icenter=1,nfrlist(iicenter)
          jcenter=ifrlist(icenter,iicenter)
          do j=1,igr
            do k=1,igr
              xgroall=xgroall+pk(k,j)*sat(j,k,jcenter)
            end do
          end do
        end do
!! gross occupation of each EFO via sat (all atoms), scaled by its net   !!
!! occupation -- this is what eos_analysis now sorts/selects on.         !!
!! parallel over i: independent per EFO, writes only s0all(i); xx0 is a  !!
!! genuine REDUCTION. c0/sat/pp0/nfrlist/ifrlist shared, read-only.      !!
        xx0=ZERO
!$OMP PARALLEL DO PRIVATE(i,icenter,jcenter,j,k,xx,xxx) REDUCTION(+:xx0)
        do i=1,imaxo
          xxx=ZERO
          do icenter=1,nfrlist(iicenter)
            jcenter=ifrlist(icenter,iicenter)
            xx=ZERO
            do j=1,igr
              do k=1,igr
                xx=xx+c0(k,i)*sat(k,j,jcenter)*c0(j,i)
              end do
            end do
            xxx=xxx+xx
          end do
          xxx=xxx*pp0(i,i)
          s0all(i)=xxx
          xx0=xx0+xxx
        end do
!$OMP END PARALLEL DO

!! store this fragment's EFOs/occupations into effao_mod (see header):   !!
!! p0 keeps the raw c0 (cubegen_new applies the fragment weight itself), !!
!! p0coef (persists across fragments, pooled for the .fchk) gets the     !!
!! projected coefficients. Fit (%) = 100*(1-xfit), per EFO.              !!
        if(imaxo.gt.0) then
          call reproject_efos_ao(igr,itotps,chp,wp,omp,scr,Sinv,
     +      c0,1,imaxo,cprojfrag,xfitfrag)
        end if
        do mu=1,imaxo
          occn(mu)=pp0(mu,mu)
          fitp(mu)=100.0d0*(1.0d0-xfitfrag(mu))
        end do
        if(iicenter.gt.1) write(*,'(a)') ''
        call print_efo_block('Fragment',iicenter,imaxo,xmaxo,occn,
     +    .true.,xx0,s0all,.true.,fitp,xnetall-xmaxo,xgroall-xx0)
        do k=1,imaxo
          do mu=1,igr
            p0(mu,k)=c0(mu,k)
            p0coef(mu,k,iicenter)=cprojfrag(mu,k)
          end do
          p0net(k,iicenter)=pp0(k,k)
          p0gro(k,iicenter)=s0all(k)
          p0fit(k,iicenter)=fitp(k)
        end do
        ip0(iicenter)=imaxo

!! optional: write a cube file for this fragment's EFOs (# CUBE section) !!
        if(icube.eq.1) call cubegen_new(iicenter,icase)

      end do !! end of the per-fragment loop !!
      if(ieos.eq.0) call json_efo_channel(icase,.true.,.true.)

      DEALLOCATE(scr,s0,sm,c0,splus,pp0,s0all,occn,fitp)
      DEALLOCATE(Sinv,cprojfrag,xfitfrag)

      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: eos_analysis
!! Effective oxidation states, one spin channel per call: pools the EFOs
!! of all fragments, sorted by gross occupation, gives the electrons of
!! the channel to the most occupied ones (sharing them evenly among
!! EFOs at the frontier with occupations within thres), prints the
!! electrons of each fragment and the reliability index R, and after
!! the beta channel the oxidation states and the EFO .fchk file.
!! arguments:
!!   idobeta (in)  0: the beta channel is not computed and is assigned
!!                 as the alpha one
!!   icase   (in)  1 alpha, 2 beta
!!   thres   (in)  occupation difference below which frontier EFOs are
!!                 treated as degenerate (EOS_THRESH)
!! author: PSalse, ERaco, MGimf
!! ----------------------------------------------------------------------
      subroutine eos_analysis(idobeta,icase,thres)

      use effao_mod, only: p0gro,ip0,p0coef,p0poolcoef,p0net,p0fit
      use output_mod, only: print_value
      use json_mod
      use input_options_mod, only: imulli

      implicit real*8(a-h,o-z)

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /loba/ oxi(maxat),errsav(maxat),elec(maxat),effpop(maxat)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /loba2/occup(nmax,2),iorbat(nmax,2),lorb(2),confi0
      common /filename/name0

      character(len=20) :: ctype
      character*60 name0
      character*80 cfile
      character*5 cspin
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort
      dimension iorbslot(nmax,2)
      dimension xlo(maxfrag),xfu(maxfrag),ilo(maxfrag),ifu(maxfrag)
      dimension nlo(maxfrag),ncount(maxfrag)
      allocatable :: occup2(:)
!! alpha electrons of each fragment, kept for the results file !!
      dimension elec1(maxfrag)
      save elec1
      integer, allocatable :: ifr(:),icol(:)
      real*8, allocatable :: xn(:),xf(:)
      logical, allocatable :: lcol(:),lok(:)
      logical lreal,lfound


      if(idobeta.eq.0.and.icase.eq.2) then
        call print_box('BETA DENSITY')
        write(*,'(2x,a)') 'The beta density equals the alpha one: '//
     +    'the beta electrons are assigned as the alpha ones'
        lorb(2)=lorb(1)
        do i=1,lorb(2)
          occup(i,2)=occup(i,1)
          iorbat(i,2)= iorbat(i,1)
        end do
        confi=confi0
        call json_begin_object('beta')
        call json_value('same_as_alpha',.true.)
        call json_end()
        go to 99
      else if(icase.eq.2.and.nb.eq.0) then
        call print_box('EOS, BETA ELECTRONS')
        write(*,'(2x,a)') 'No beta electrons to assign'
        do i=1,icufr
          elec(i)=ZERO
        end do
        lorb(2)=0
        p0poolcoef(:,:,2)=ZERO
        call json_begin_object('beta')
        call json_value('no_electrons',.true.)
        call json_end()
!! no beta occupied/unoccupied boundary: the overall R is the alpha one !!
        confi=confi0
        go to 99
      end if

!! the pooled list holds nmax EFOs (EFF_THRESH 0 or below keeps all of  !!
!! every fragment's)                                                    !!
      iorb=0
      do i=1,icufr
        iorb=iorb+ip0(i)
      end do
      if(iorb.gt.nmax) then
        write(*,'(2x,a,i0,a,i0,a)') 'EOS: ',iorb,' EFOs above the '//
     +    'EFF_THRESH cutoff, at most ',nmax,': raise EFF_THRESH'
        call apost_stop(' EOS: too many EFOs (EFF_THRESH)')
      end if

      iorb=0
      do i=1,icufr
        elec(i)=ZERO
        do k=1,ip0(i)
          iorb=iorb+1
          occup(iorb,icase)=p0gro(k,i)
          iorbat(iorb,icase)=i
          iorbslot(iorb,icase)=k
        end do
      end do
      if(icase.eq.1) then
        call print_box('EOS, ALPHA ELECTRONS')
      else
        call print_box('EOS, BETA ELECTRONS')
      end if
      call print_value('EFOs in the assignment',iorb)

      lorb(icase)=iorb
      ALLOCATE(occup2(iorb))

!! pool every fragment's EFOs into one list, sorted by occupation        !!
!! (descending); iorbat/iorbslot track which fragment and local EFO      !!
!! index each slot came from.                                            !!
      do i=1,iorb-1
        do j=i+1,iorb
          if (occup(j,icase).gt.occup(i,icase))then
            xkk=occup(j,icase)
            occup(j,icase)=occup(i,icase)
            occup(i,icase)=xkk
            ikk=iorbat(j,icase)
            iorbat(j,icase)=iorbat(i,icase)
            iorbat(i,icase)=ikk
            ikk=iorbslot(j,icase)
            iorbslot(j,icase)=iorbslot(i,icase)
            iorbslot(i,icase)=ikk
          end if
        end do
      end do

!! pooled+sorted coefficient columns, truncated to igr and zero-padded  !!
!! -- feeds the shared .fchk EFO writer (print.f), same convention as   !!
!! GEOS's ueos_analysis. Left serial: at most igr*igr copies, once per  !!
!! icase, negligible next to the diagonalizations done above.           !!
      p0poolcoef(:,:,icase)=ZERO
      do jj=1,MIN(iorb,igr)
        ii2=iorbat(jj,icase)
        kk2=iorbslot(jj,icase)
        do mu=1,igr
          p0poolcoef(mu,jj,icase)=p0coef(mu,kk2,ii2)
        end do
      end do

!! nnn is the LO EFO index in the pooled, sorted list (alpha/beta          !!
!! electron count); k/kk scan outward from it for EFOs within thres of     !!
!! occup(nnn) -- a quasi-degenerate block spanning the occupied/unoccupied !!
!! boundary, which occup2 below turns into fractional occupations.         !!
      k=0
      nnn=nalf
      if(icase.eq.2) nnn=nb
!! a large EFF_THRESH can leave fewer pooled EFOs than electrons !!
      if(nnn.gt.iorb) then
        if(icase.eq.2) then
          write(*,'(2x,a,i0,a,i0,a)') 'EOS: only ',iorb,
     +      ' EFOs above the EFF_THRESH cutoff for ',nnn,
     +      ' beta electrons: lower EFF_THRESH'
        else
          write(*,'(2x,a,i0,a,i0,a)') 'EOS: only ',iorb,
     +      ' EFOs above the EFF_THRESH cutoff for ',nnn,
     +      ' alpha electrons: lower EFF_THRESH'
        end if
        call apost_stop(' EOS: fewer EFOs than electrons (EFF_THRESH)')
      end if
333   k=k+1
      if(nnn+k.le.iorb) then
        if(dabs(occup(nnn,icase)-occup(nnn+k,icase)).lt.thres) go to 333
      end if
      k=k-1
      if(k.eq.0) then
        call print_value('Electron assignment','integer')
        ifrac=0
        do i=1,iorb
          if(i.le.nnn) then
            occup2(i)=ONE
          else
            occup2(i)=ZERO
          end if
        end do
      else
        call print_value('Electron assignment','fractional')
        ifrac=1
        kk=0
334     kk=kk+1
        if(nnn-kk.ge.1) then
          if(dabs(occup(nnn,icase)-occup(nnn-kk,icase)).lt.thres) go to 334
        end if
        kk=kk-1
!! the kk+1 electrons of the block's occupied EFOs, spread evenly over !!
!! its kk+k+1 quasi-degenerate EFOs                                    !!
        frac=float(kk+1)/float(kk+k+1)
        if(kk.eq.0) then
          write(*,'(2x,a,i0,a,f7.5,a)') 'WARNING: 1 electron shared '//
     +      'by ',kk+k+1,' EFOs with occupations within ',thres,
     +      ' of each other'
        else
          write(*,'(2x,a,i0,a,i0,a,f7.5,a)') 'WARNING: ',kk+1,
     +      ' electrons shared by ',kk+k+1,' EFOs with occupations '//
     +      'within ',thres,' of each other'
        end if
        do i=1,iorb
          if(i.lt.nnn-kk) then
            occup2(i)=ONE
          else if(i.gt.nnn+k) then
            occup2(i)=ZERO
          else
            occup2(i)=frac
          end if
        end do
      end if

      do i=1,iorb
        elec(iorbat(i,icase))=elec(iorbat(i,icase))+occup2(i)
      end do


!! frontier EFOs of each fragment: with nn its electrons rounded up, the  !!
!! nn-th (last occupied, xlo) and (nn+1)-th (first unoccupied, xfu) of   !!
!! its EFOs in gross order, i.e. as met in the pooled list, which is      !!
!! sorted by gross occupation; ilo/ifu: 1 if that EFO exists              !!
      do i=1,icufr
        nn=int(elec(i))
        if(elec(i)-nn.gt.thresh) nn=nn+1
        nlo(i)=nn
        ncount(i)=0
        ilo(i)=0
        ifu(i)=0
        xlo(i)=ZERO
        xfu(i)=ZERO
      end do
      do ii=1,iorb
        i=iorbat(ii,icase)
        ncount(i)=ncount(i)+1
        if(ncount(i).eq.nlo(i)) then
          xlo(i)=occup(ii,icase)
          ilo(i)=1
        else if(ncount(i).eq.nlo(i)+1) then
          xfu(i)=occup(ii,icase)
          ifu(i)=1
        end if
      end do

!! xlast: worst (smallest) "last occupied EFO" gross occupation across    !!
!! fragments -- the weakest-defined occupied/unoccupied boundary.         !!
      xlast=ONE
      ilast=0
      do i=1,icufr
        if(ilo(i).eq.1) then
          if(xlo(i).lt.xlast) then
            xlast=xlo(i)
            ilast=i
          end if
        end if
      end do

!! xfirst: best (largest) "first unoccupied EFO" occupation among the     !!
!! other fragments -- how close the runner-up comes to being occupied.    !!
      xfirst=ZERO
      do i=1,icufr
        if(i.ne.ilast.and.ifu(i).eq.1) then
          if(xfu(i).gt.xfirst) xfirst=xfu(i)
        end if
      end do

!! reliability index: 100% if the occupied/unoccupied gap is >=0.5,      !!
!! scaled down as xlast and xfirst get closer together.                  !!
      confi=100.0*min(1.0d0,xlast-xfirst+0.5d0)
      write(*,'(a)') ''
      call print_eos_frontier(icufr,elec,xlo,ilo,xfu,ifu,confi)
      if(icase.eq.1) then
        do i=1,icufr
          elec1(i)=elec(i)
        end do
      end if

!! results file: this channel's R, frontier table and pooled EFOs (in  !!
!! the order of the .fchk columns, with the electrons each one gets)   !!
      call atomic_definition_names(ckey,cname,cshort,cspace)
      cspin='alpha'
      if(icase.eq.2) cspin='beta'
      lreal=imulli.eq.0
      call json_begin_object(trim(cspin))
      call json_value('reliability_percent',confi)
      call json_index('reliability_percent','scalar','Reliability '//
     +  'index R, '//trim(cspin)//' ('//trim(cshort)//')','%',3,
     +  trim(cspin),.false.,' ')
      if(ifrac.eq.0) then
        call json_value('assignment','integer')
      else
        call json_value('assignment','fractional')
      end if
      allocate(lok(icufr))
      call json_begin_object('frontier')
      call json_value('rows','fragment')
      call json_value('n_rows',icufr)
      call json_begin_array('columns')
      call json_column('last_occupied','Last occupied EFO',' ',3,
     +  .false.)
      call json_column('first_unoccupied','First unoccupied EFO',' ',3,
     +  .false.)
      call json_end()
      call json_begin_object('data')
      do i=1,icufr
        lok(i)=ilo(i).eq.1
      end do
      call json_vector_or_null('last_occupied',xlo(1:icufr),lok,icufr)
      do i=1,icufr
        lok(i)=ifu(i).eq.1
      end do
      call json_vector_or_null('first_unoccupied',xfu(1:icufr),lok,
     +  icufr)
      call json_end()
      call json_end()
      call json_index('frontier','table','EOS frontier EFOs, '//
     +  trim(cspin)//' ('//trim(cshort)//')',' ',-1,trim(cspin),
     +  .false.,' ')
      deallocate(lok)
!! the .fchk file holds the first min(igr, independent functions) EFOs !!
      ncolf=min(igr,int_locate(15,"Number of independ",lfound))
      allocate(ifr(iorb),icol(iorb),xn(iorb),xf(iorb),lcol(iorb))
      do jj=1,iorb
        ifr(jj)=iorbat(jj,icase)
        xn(jj)=p0net(iorbslot(jj,icase),ifr(jj))
        xf(jj)=ZERO
        if(lreal) xf(jj)=p0fit(iorbslot(jj,icase),ifr(jj))
        icol(jj)=jj
        lcol(jj)=jj.le.ncolf
      end do
      call json_efo_table(iorb,ifr,xn,lreal,occup(1:iorb,icase),lreal,
     +  xf,.true.,occup2,icol,lcol,trim(cspin),trim(cspin))
      deallocate(ifr,icol,xn,xf,lcol)
      call json_end()
      DEALLOCATE(occup2)
      if(icase.eq.1) confi0=confi

99    continue

      zztot=ZERO
      do i=1,icufr
        if(icase.eq.1) then
          zzn=ZERO
          do j=1,nfrlist(i)
            zzn=zzn+zn(ifrlist(j,i))
          end do
          oxi(i)=zzn-elec(i)
        else
          oxi(i)=oxi(i)-elec(i)
          zztot=zztot+oxi(i)
        end if
      end do

      if(icase.eq.2) then
        confi2=dmin1(confi,confi0)
        call print_oxidation_states(icufr,oxi,zztot,.true.,confi2)

!! results file: per-fragment results, overall R, the EFO file !!
        call atomic_definition_names(ckey,cname,cshort,cspace)
        call json_begin_object('fragments')
        call json_vector('oxidation_states',oxi(1:icufr),icufr)
        call json_index('oxidation_states','fragment_vector',
     +    'EOS oxidation states ('//trim(cshort)//')',' ',2,'total',
     +    .true.,' ')
        call json_vector('electrons_alpha',elec1(1:icufr),icufr)
        call json_index('electrons_alpha','fragment_vector',
     +    'EOS electrons, alpha ('//trim(cshort)//')',' ',2,'alpha',
     +    .false.,' ')
        call json_vector('electrons_beta',elec(1:icufr),icufr)
        call json_index('electrons_beta','fragment_vector',
     +    'EOS electrons, beta ('//trim(cshort)//')',' ',2,'beta',
     +    .false.,' ')
        call json_end()
        call json_value('reliability_percent',confi2)
        call json_index('reliability_percent','scalar','Reliability '//
     +    'index R ('//trim(cshort)//')','%',3,'total',.false.,' ')
        cfile=trim(name0)//"-EOS-EFOs.fchk"
        cfile=cfile(index(cfile,'/',back=.true.)+1:)
        call json_value('orbital_file',trim(cfile))
        call json_value('orbital_file_id',json_run_id())

!! pooled EOS EFOs as fake Alpha/Beta MOs in a .fchk, for visualization !!
!! in any standard viewer -- printed by default, same as OSLO/GEOS's    !!
!! own .fchk output, for every population scheme. A                    !!
!! restricted wavefunction (kop=0) has no Beta blocks to splice into,   !!
!! but its beta channel isn't always empty -- e.g. a restricted-orbital !!
!! open-shell CASSCF source (idobeta=1 despite kop=0, main.f's own      !!
!! dispatch condition) genuinely populates it. Gated on idobeta itself, !!
!! not lorb(2): the idobeta=0 branch above sets lorb(2)=lorb(1) purely  !!
!! to duplicate occup/iorbat for this print, without ever touching      !!
!! p0poolcoef(:,:,2) -- lorb(2)>0 there does NOT mean real coefficients !!
!! exist. rwfu_effao_orbprint inserts a synthetic Beta MO set next to   !!
!! the real Alpha one when idobeta=1, same convention as GEOS's writer; !!
!! otherwise fall back to the plain restricted, single-channel writer.  !!
        ctype="-EOS-EFOs"
        if(kop.eq.0) then
          if(idobeta.eq.1) then
            call rwfu_effao_orbprint(p0poolcoef(:,:,1),p0poolcoef(:,:,2),ctype)
          else
            call rwf_effao_orbprint(p0poolcoef(:,:,1),ctype)
          end if
        else
          call uwf_effao_orbprint(p0poolcoef(:,:,1),p0poolcoef(:,:,2),ctype)
        end if
      end if

      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: ueffaolow_frag
!! Hilbert-space EFOs with Lowdin, Davidson-Lowdin or NAO atoms: for
!! each fragment, the block of the density on its orthogonal basis
!! functions, diagonalized; prints the occupations and stores them in
!! effao_mod (EFOs already in the AO basis, net = gross).
!! arguments:
!!   icase (in)  0 total density, 1 alpha, 2 beta
!! author: PSalse, ERaco, MGimf
!! ----------------------------------------------------------------------
      subroutine ueffaolow_frag(icase)

      use basis_set
      use ao_matrices
      use effao_mod, only: p0,p0net,p0gro,ip0,p0coef
      use nao_mod, only: unao,ssnao
      use input_options_mod, only: ieos !! replaces common /nao/ -- see modules.f90 !!

      implicit real*8(A-H,O-Z)

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      dimension iao_frag(nmax)

      allocatable :: s0(:,:), sm(:,:), c0(:,:), splus(:,:), pp0(:,:)
      allocatable :: occn(:)

      icube   = Iopt(13)
      ieffthr = Iopt(24)
      imulli  = Iopt(5)

      xminocc=REAL(ieffthr)/1000.0d0

      ALLOCATE(s0(igr,igr),sm(igr,igr),splus(igr,igr))
      ALLOCATE(c0(igr,igr),pp0(igr,igr),occn(igr))

!! AO-to-fragment map. !!
      iao_frag=0
      do ifrag=1,icufr
        do icenter=1,nfrlist(ifrag)
          jcenter=ifrlist(icenter,ifrag)
          do mu=llim(jcenter),iulim(jcenter)
            iao_frag(mu)=ifrag
          end do
        end do
      end do

      if(icase.eq.0) then
        s0=p
      else if (icase.eq.1) then
        s0=pa
      else if (icase.eq.2) then
        s0=pb
      end if

      if(imulli.eq.4) then
!! NAO-to-AO matrix -- note the transpose, ssnao is nonsymmetric. !!
        do i=1,igr
          do j=1,igr
            splus(i,j)=ssnao(j,i)
            sm(i,j)=unao(i,j)
          end do
        end do
      else if(imulli.eq.3) then
!! Davidson-Lowdin orthogonal basis, same as its populations (tolow) !!
        call davidson_lowdin(igr,splus,sm)
      else
        splus=s12p
        sm=s12m
      end if
      if(icase.eq.1) then
        call print_box('ALPHA DENSITY')
      else if(icase.eq.2) then
        call print_box('BETA DENSITY')
      end if

!! transform P with Splus into the orthogonalized basis -- EFOs are      !!
!! back-transformed to the AO basis after diagonalization, below.        !!
      call to_lowdin_basis(igr,splus,s0)

!! loop over fragments -- one diagonalization per fragment avoids having !!
!! to reorder eigenvalues afterward.                                     !!
      do ifrag=1,icufr

!! block-diagonal P for this fragment. parallel over (mu,nu): each      !!
!! iteration writes only its own pp0(mu,nu), independent across pairs.  !!
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(mu,nu)
        do mu=1,igr
          do nu=1,igr
            pp0(mu,nu)=ZERO
            if(iao_frag(mu).eq.ifrag.and.iao_frag(nu).eq.ifrag) pp0(mu,nu)=s0(mu,nu)
          end do
        end do
!$OMP END PARALLEL DO

        call diagonalize(igr,igr,pp0,C0,0)
!! back-transform: EFOs expanded over the full basis set. !!
        call to_AO_basis(igr,igr,sm,C0)

!! number of basis functions of the fragment: at most that many EFOs !!
        imaxeff=0
        do i=1,igr
          if(iao_frag(i).eq.ifrag) imaxeff=imaxeff+1
        end do

!! kept EFOs: up to all imaxeff of the fragment, above the cutoff;     !!
!! xmaxo sums their net occupations, xnetall all of them (the matrix is !!
!! zero outside the fragment's block, so this is its net population)    !!
        imaxo=0
        do i=1,imaxeff
          if(pp0(i,i).lt.xminocc) exit
          imaxo=i
        end do
        xmaxo=ZERO
        do i=1,imaxo
          xmaxo=xmaxo+pp0(i,i)
        end do
        xnetall=ZERO
        do i=1,igr
          xnetall=xnetall+pp0(i,i)
        end do

        do mu=1,imaxo
          occn(mu)=pp0(mu,mu)
        end do
        if(ifrag.gt.1) write(*,'(a)') ''
        call print_efo_block('Fragment',ifrag,imaxo,xmaxo,occn,.false.,
     +    ZERO,occn,.false.,occn,xnetall-xmaxo,ZERO)

!! save this fragment's EFOs (coefficients and occupations) into the    !!
!! shared effao_mod arrays; the coefficients are already in the AO      !!
!! basis, so p0coef (pooled for the .fchk) takes them as they are       !!
        do k=1,imaxo
          do mu=1,igr
            p0(mu,k)=c0(mu,k)
            p0coef(mu,k,ifrag)=c0(mu,k)
          end do
          p0net(k,ifrag)=pp0(k,k)
          p0gro(k,ifrag)=pp0(k,k)
        end do
        ip0(ifrag)=imaxo

!! cube file for visualization, if requested. !!
        if(icube.eq.1) call cubegen_new(ifrag,icase)

      end do

      if(ieos.eq.0) call json_efo_channel(icase,.false.,.false.)

      DEALLOCATE(s0,sm,c0,splus,pp0,occn)

      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: ueffaomull_frag
!! Hilbert-space EFOs with Mulliken atoms: for each fragment, the block
!! of the density on its basis functions, orthogonalized within the
!! fragment and diagonalized; prints the occupations and stores them in
!! effao_mod (net = gross).
!! arguments:
!!   icase (in)  0 total density, 1 alpha, 2 beta
!! author: PSalse, ERaco, MGimf
!! ----------------------------------------------------------------------
      subroutine ueffaomull_frag(icase)

      use basis_set
      use ao_matrices
      use effao_mod, only: p0,p0net,p0gro,ip0,p0coef
      use input_options_mod, only: ieos

      implicit real*8(A-H,O-Z)

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      dimension iao_frag(nmax)

      allocatable :: s0(:,:),sm(:,:),c0(:,:),splus(:,:),pp0(:,:)
      allocatable :: pk(:,:),occn(:)

      icube   = Iopt(13)
      ieffthr = Iopt(24)

      xminocc=real(ieffthr)/1000.0d0

      ALLOCATE(s0(igr,igr),sm(igr,igr),splus(igr,igr))
      ALLOCATE(c0(igr,igr),pp0(igr,igr),pk(igr,igr),occn(igr))

      if(icase.eq.0) then
        pk=p
      else if(icase.eq.1) then
        call print_box('ALPHA DENSITY')
        pk=pa
      else if(icase.eq.2) then
        call print_box('BETA DENSITY')
        pk=pb
      end if

!! AO-to-fragment map. !!
      iao_frag=0
      do ifrag=1,icufr
        do icenter=1,nfrlist(ifrag)
          jcenter=ifrlist(icenter,ifrag)
          do mu=llim(jcenter),iulim(jcenter)
            iao_frag(mu)=ifrag
          end do
        end do
      end do

!! loop over fragments -- one diagonalization per fragment avoids having !!
!! to reorder eigenvalues afterward.                                     !!
      do ifrag=1,icufr

!! block-diagonal S and P (AO basis) for this fragment. parallel over    !!
!! (mu,nu): each iteration writes only its own pp0(mu,nu)/s0(mu,nu),     !!
!! independent across pairs.                                             !!
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(mu,nu)
        do mu=1,igr
          do nu=1,igr
            pp0(mu,nu)=ZERO
            s0(mu,nu)=ZERO
            if(iao_frag(mu).eq.ifrag.and.iao_frag(nu).eq.ifrag) then
              pp0(mu,nu)=pk(mu,nu)
              s0(mu,nu)=s(mu,nu)
            end if
          end do
        end do
!$OMP END PARALLEL DO

!! S0^-1/2, then transform the block P0 with it. !!
        call build_Smp(igr,s0,Sm,Splus,0)
        call to_lowdin_basis(igr,Splus,pp0)
        call diagonalize(igr,igr,pp0,C0,0)
!! back-transform: EFOs expanded over the full basis set. !!
        call to_AO_basis(igr,igr,Sm,C0)

        imaxeff=0
        do i=1,igr
          if(iao_frag(i).eq.ifrag) imaxeff=imaxeff+1
        end do
!! kept EFOs: up to all imaxeff of the fragment, above the cutoff;     !!
!! xmaxo sums their net occupations, xnetall all of them (the matrix is !!
!! zero outside the fragment's block, so this is its net population)    !!
        imaxo=0
        do i=1,imaxeff
          if(pp0(i,i).lt.xminocc) exit
          imaxo=i
        end do
        xmaxo=ZERO
        do i=1,imaxo
          xmaxo=xmaxo+pp0(i,i)
        end do
        xnetall=ZERO
        do i=1,igr
          xnetall=xnetall+pp0(i,i)
        end do

        do mu=1,imaxo
          occn(mu)=pp0(mu,mu)
        end do
        if(ifrag.gt.1) write(*,'(a)') ''
        call print_efo_block('Fragment',ifrag,imaxo,xmaxo,occn,.false.,
     +    ZERO,occn,.false.,occn,xnetall-xmaxo,ZERO)

!! save this fragment's EFOs (coefficients and occupations) into the    !!
!! shared effao_mod arrays; the coefficients are already in the AO      !!
!! basis, so p0coef (pooled for the .fchk) takes them as they are       !!
        do k=1,imaxo
          do mu=1,igr
            p0(mu,k)=c0(mu,k)
            p0coef(mu,k,ifrag)=c0(mu,k)
          end do
          p0net(k,ifrag)=pp0(k,k)
          p0gro(k,ifrag)=pp0(k,k)
        end do
        ip0(ifrag)=imaxo

!! cube file for visualization, if requested. !!
        if(icube.eq.1) call cubegen_new(ifrag,icase)

      end do

      if(ieos.eq.0) call json_efo_channel(icase,.false.,.false.)

      DEALLOCATE(s0,sm,c0,splus,pp0,pk,occn)

      end

!! ----------------------------------------------------------------------
!! subroutine: ueffao3d
!! Real-space effective atomic orbitals (EFAOs) of single atoms
!! (DOATOMS, the atoms of # ATOMS): as ueffao3d_frag, with the atomic
!! weight in place of the fragment weight and without the .fchk file.
!! arguments: as ueffao3d_frag
!! author:
!! ----------------------------------------------------------------------
      subroutine ueffao3d(itotps,ndim,omp,chp,sat,wp,omp2,pk,icase)

      use integration_grid
      use effao_mod, only: p0,p0net,p0gro,ip0

      implicit real*8(a-h,o-z)

      include 'parameter.h'

      integer,intent(in) :: itotps,ndim

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(200)
      common /atlist/iatlist(maxat),icuat

      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat)
      dimension wp(itotps)
      dimension sat(ndim,ndim,nat)
      dimension pk(ndim,ndim)

      allocatable :: s0(:,:),sm(:,:),c0(:,:),splus(:,:),pp0(:,:)
      allocatable :: s0all(:),occn(:)

      icube   = Iopt(13)
      ieffthr = Iopt(24)
      iatps   = nang*nrad

      xmaxocc=REAL(ieffthr)/1000.0d0

      if(icase.eq.1) then
        call print_box('ALPHA DENSITY')
      else if(icase.eq.2) then
        call print_box('BETA DENSITY')
      end if

      nat0=nat
      if(icuat.ne.nat) nat0=icuat

      ALLOCATE(s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      ALLOCATE(c0(ndim,ndim),pp0(ndim,ndim),occn(ndim))

!! per-atom loop kept serial on purpose: icuat/nat0 can be small on a    !!
!! large system (or a DOATOMS-restricted handful) -- the O(igr^2)/O(igr^3)!!
!! work inside each iteration is threaded instead, same lesson already   !!
!! applied to ueffao3d_frag's per-fragment loop above.                   !!
      do iicenter=1,nat0
        icenter=iatlist(iicenter)

!! net AO overlap block (s0) at this atom, ALLPOINTS integration.        !!
!! parallel over mu: each mu writes only its own s0(mu,*)/s0(*,mu), no   !!
!! two mu iterations collide. nu nested serial (bounds depend on mu).    !!
!! wp/chp/omp/omp2 shared read-only, x private. dominant cost here.      !!
!$OMP PARALLEL DO PRIVATE(mu,nu,jcenter,ifut,x)
        do mu=1,ndim
          do nu=1,mu
            x=ZERO
            do jcenter=1,nat
              do ifut=iatps*(jcenter-1)+1,iatps*jcenter
                x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*omp2(ifut,icenter)*omp2(ifut,icenter)*omp(ifut)
              end do
            end do
            s0(mu,nu)=x
            s0(nu,mu)=x
          end do
        end do
!$OMP END PARALLEL DO

!! S^-1/2 and S^+1/2 of the atomic overlap block !!
        call build_Smp(igr,S0,Sm,Splus,0)

!! transform pk into the orthogonalized atomic basis, diagonalize to    !!
!! get the EFAOs (c0) and net occupations (pp0 diagonal), back to AO.   !!
        do i=1,igr
          do j=1,igr
            pp0(i,j)=pk(i,j)
          end do
        end do
        call to_lowdin_basis(igr,Splus,pp0)
        call diagonalize(igr,igr,pp0,C0,0)
        call to_AO_basis(igr,igr,Sm,C0)

!! keep EFAOs above the EFF_THRESH cutoff (xmaxocc); pp0's diagonal is  !!
!! already sorted decreasing, so this is a simple prefix. imaxo starts   !!
!! at 0: an atom with no EFAO above the cutoff keeps none.               !!
        imaxo=0
        do i=1,igr
          if(pp0(i,i).lt.xmaxocc) exit
          imaxo=i
        end do
!! kept EFAOs (xmaxo) and all of them (xnetall); Tr(P S^A) (xgroall) !!
!! for the gross part of what the cutoff leaves out                    !!
        xmaxo=ZERO
        do i=1,imaxo
          xmaxo=xmaxo+pp0(i,i)
        end do
        xnetall=ZERO
        do i=1,igr
          xnetall=xnetall+pp0(i,i)
        end do
        xgroall=ZERO
        do j=1,igr
          do k=1,igr
            xgroall=xgroall+pk(k,j)*sat(j,k,icenter)
          end do
        end do

!! gross occupation of each EFAO via sat, scaled by its net occupation. !!
!! parallel over i: independent per EFAO, writes only s0all(i); xx0 is  !!
!! a genuine REDUCTION. c0/sat/pp0 shared, read-only.                   !!
        xx0=ZERO
!$OMP PARALLEL DO PRIVATE(i,j,k,xx) REDUCTION(+:xx0)
        do i=1,imaxo
          xx=ZERO
          do j=1,igr
            do k=1,igr
              xx=xx+c0(k,i)*sat(k,j,icenter)*c0(j,i)
            end do
          end do
          xx=xx*pp0(i,i)
          s0all(i)=xx
          xx0=xx0+xx
        end do
!$OMP END PARALLEL DO
        do mu=1,imaxo
          occn(mu)=pp0(mu,mu)
        end do
        if(iicenter.gt.1) write(*,'(a)') ''
        call print_efo_block('Atom',icenter,imaxo,xmaxo,occn,.true.,
     +    xx0,s0all,.false.,occn,xnetall-xmaxo,xgroall-xx0)

!! store this atom's EFAOs/occupations into effao_mod (see header). !!
        do k=1,imaxo
          do mu=1,igr
            p0(mu,k)=c0(mu,k)
          end do
          p0net(k,icenter)=pp0(k,k)
          p0gro(k,icenter)=s0all(k)
        end do
        ip0(icenter)=imaxo

!! optional: write a cube file for this atom's EFAOs (# CUBE section) !!
        if(icube.eq.1) call cubegen_new(icenter,icase)

      end do !! end of the per-atom loop !!

      DEALLOCATE(s0,sm,c0,splus,pp0,s0all,occn)

      end

