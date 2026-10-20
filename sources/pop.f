!! ----------------------------------------------------------------------
!! subroutine: fborder
!! Bond orders of every atom pair, Tr(P S^A P S^B) + Tr(Ps S^A Ps S^B),
!! with the localization index on the diagonal, and the total, bonded
!! and free valences of each atom (I. Mayer, P. Salvador, Chem. Phys.
!! Lett. 2004, 383, 368). Works for any atomic definition: sat holds it.
!! Prints the matrix and the three valence tables.
!! arguments:
!!   sat (in)   (nbasis,nbasis,natoms) atomic overlap matrices S^A
!!   val (out)  (maxat,3) total, bonded and free valences
!! author:
!! ----------------------------------------------------------------------
      subroutine fborder(sat,val)
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/ op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/ qat(maxat,2),qsat(maxat,2)
      dimension sat(nbasis,nbasis,natoms),val(maxat,3)

      dimension rindex(maxat,maxat),tindex(maxat),diag(maxat)
      allocatable tt(:,:,:)
      allocate(tt(nbasis,nbasis,natoms))

!! P*S^A matrix products, contracted into the R-index trace below. !!
      do iat=1,natoms
        do mu=1,nbasis
          do nu=1,nbasis
            x=0.d0
            do itau=1,nbasis
              x=x+p(mu,itau)*sat(itau,nu,iat)
            enddo
            tt(mu,nu,iat)=x
          enddo
        enddo
      enddo

      do iat=1,natoms
        do ibt=iat,natoms
          x=0.d0
          do mu=1,nbasis
            do nu=1,nbasis
              x=x+tt(mu,nu,iat)*tt(nu,mu,ibt)
            enddo
          enddo
          rindex(iat,ibt)=x
          rindex(ibt,iat)=x
        enddo
!! diag(iat) is taken before the Ps term below: the total valence      !!
!! follows Mayer's open-shell definition (from P only), while the bond  !!
!! orders get the Ps term, so the one-center spin part shows up as the  !!
!! free valence (total valence minus the sum of bond orders).           !!
        diag(iat)=rindex(iat,iat)
      enddo

!! Ps contribution (open-shell only) !!
      if(kop.ne.0) then
        do iat=1,natoms
          do mu=1,nbasis
            do nu=1,nbasis
              x=0.d0
              do itau=1,nbasis
                x=x+ps(mu,itau)*sat(itau,nu,iat)
              enddo
              tt(mu,nu,iat)=x
            enddo
          enddo
        enddo
        do iat=1,natoms
          do ibt=iat,natoms
            x=0.d0
            do mu=1,nbasis
              do nu=1,nbasis
                x=x+tt(mu,nu,iat)*tt(nu,mu,ibt)
              enddo
            enddo
            rindex(iat,ibt)=rindex(iat,ibt)+x
            rindex(ibt,iat)=rindex(iat,ibt)
          enddo
        enddo
      endif

!! bo keeps the bond orders; di starts as a copy and is replaced by the !!
!! 2-RDM delocalization indices in a correlated local spin analysis     !!
      do i=1,natoms
        do j=1,natoms
          if(i.eq.j) then
            bo(i,i)=rindex(i,i)*0.5d0
          else
            bo(i,j)=rindex(i,j)
          endif
          di(i,j)=bo(i,j)
        enddo
      enddo

!! zero the diagonal so it can be reused as the valence sum below !!
      do i=1,natoms
        rindex(i,i)=0.d0
      enddo

      call print_box('BOND ORDER MATRIX')
      call mprint2(bo,natoms,maxat)

!! valence numbers !!
      do i=1,natoms
        x=0.d0
        do j=1,natoms
          x=x+rindex(i,j)
        enddo
        rindex(i,i)=x
      enddo
      do i=1,natoms
        tindex(i)=rindex(i,i)
        diag(i)=2.d0*qat(i,1)-diag(i)
      enddo

      call print_atom_table('TOTAL VALENCES',' ',
     +  'FRAGMENT TOTAL VALENCES',diag)
      call print_atom_table('VALENCES USED IN BONDS',
     +  'Sum of the bond orders of each atom',
     +  'FRAGMENT VALENCES USED IN BONDS',tindex)
      do i=1,natoms
        val(i,1)=diag(i)
        val(i,2)=tindex(i)
        diag(i)=diag(i)-tindex(i)
        val(i,3)=diag(i)
      enddo
      call print_atom_table('FREE VALENCES',' ',
     +  'FRAGMENT FREE VALENCES',diag)

      deallocate(tt)

      return
      end

      

!! ----------------------------------------------------------------------
!! subroutine: opop
!! Overlap populations of every atom pair, the integral of w_A w_B rho
!! over the molecular grid (I. Mayer, P. Salvador, Chem. Phys. Lett.
!! 2004, 383, 368). By default each atom's own grid points are used;
!! with ALLPOINTS rho is weighted once and every pair is summed over the
!! whole grid.
!! arguments:
!!   wp   (in)     (itotps) integration weights
!!   omp  (in)     (itotps) Becke weight of each point (ALLPOINTS only)
!!   omp2 (in)     (itotps,natoms) atomic weight w_A of each point
!!   rho  (inout)  (itotps) density; ALLPOINTS rescales it in place
!! author:
!! ----------------------------------------------------------------------
      subroutine opop(wp,omp,omp2,rho)
      use basis_set, only: natoms
      use ao_matrices
      use integration_grid
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /actual/ iact,jat,icenter
      common /ovpop/ op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/ iopt(200)
      common /achi/ achi(maxat,maxat),ibcp
      dimension wp(nrad*nang*natoms),rho(nrad*nang*natoms)
      dimension omp(nrad*nang*natoms),omp2(nrad*nang*natoms,natoms)

      iallpo=iopt(7)

      iatps=nrad*nang
      itotps=iatps*natoms

!! not parallelized: called once per calculation on a small (natoms^2) !!
!! matrix, negligible next to the numerical-integration hot paths.     !!
      do i=1,natoms
        do j=1,natoms
          op(i,j)=0.0d0
        enddo
      enddo

      if(iallpo.eq.0) then
!! default scheme: each atom integrated over its own grid points !!
        do icenter=1,natoms
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            do jcenter=icenter,natoms
              op(icenter,jcenter)=op(icenter,jcenter)+wp(ifut)*rho(ifut)*omp2(ifut,icenter)*omp2(ifut,jcenter)
            enddo
          enddo
        enddo
      else if(iallpo.eq.1) then
!! ALLPOINTS scheme: pre-weight rho once, then re-sum it over the full !!
!! grid for every atom pair.                                          !!
        do kcenter=1,natoms
          do ifut=iatps*(kcenter-1)+1,iatps*kcenter
            rho(ifut)=wp(ifut)*omp(ifut)*rho(ifut)
          enddo
        enddo
        do ifut=1,itotps
          do icenter=1,natoms
            do jcenter=icenter,natoms
              op(icenter,jcenter)=op(icenter,jcenter)+rho(ifut)*omp2(ifut,icenter)*omp2(ifut,jcenter)
            enddo
          enddo
        enddo
      endif

      do icenter=1,natoms
        do jcenter=icenter+1,natoms
          op(jcenter,icenter)=op(icenter,jcenter)
        enddo
      enddo

      return
      end

!! ----------------------------------------------------------------------
!! subroutine: fspindec
!! Local spin analysis of a single-determinant wavefunction: effectively
!! unpaired electrons of each atom and the a=3/4 decomposition of <S^2>
!! into atomic and diatomic terms (E. Ramos-Cordoba et al., J. Chem.
!! Theory Comput. 2012, 8, 1270). spincorr (corr.f) does the same for
!! correlated wavefunctions.
!! arguments:
!!   sat (in)  (nbasis,nbasis,natoms) atomic overlap matrices S^A
!! author:
!! ----------------------------------------------------------------------
      subroutine fspindec(sat)
      use basis_set
      use ao_matrices
      use input_options_mod, only: idofr
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /printout/iaccur
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /qat/ qat(maxat,2),qsat(maxat,2)
      common /localspin/ xlsa(maxat,maxat),ua(maxat)
      dimension sat(nbasis,nbasis,natoms)
      dimension rindex(maxat,maxat)
      allocatable tt(:,:,:),tts(:,:)

      allocate(tt(nbasis,nbasis,natoms))
      allocate(tts(nbasis,nbasis))

      do mu=1,nbasis
        do nu=1,nbasis
          tts(mu,nu)=0.0d0
        enddo
      enddo

!! serial: the O(nbasis^2*natoms^2) loops below are cheap next to the !!
!! numerical integration                                              !!

!! Ps*S^A products, tts is their sum over atoms !!
      do iat=1,natoms
        do mu=1,nbasis
          do nu=1,nbasis
            x=0.d0
            do itau=1,nbasis
              x=x+ps(mu,itau)*sat(itau,nu,iat)
            enddo
            tt(mu,nu,iat)=x
            tts(mu,nu)=tts(mu,nu)+x
          enddo
        enddo
      enddo

!! effectively unpaired electrons !!
      sum=0.0d0
      do iat=1,natoms
        x=0.0d0
        do mu=1,nbasis
          do nu=1,nbasis
            x=x+tt(mu,nu,iat)*tts(nu,mu)
          enddo
        enddo
        ua(iat)=x
        sum=sum+x
      enddo

      call atomic_definition_names(ckey,cname,cshort,cspace)
      ndec=6
      if(iaccur.ne.0) ndec=13
      call print_atom_table('EFFECTIVELY UNPAIRED ELECTRONS',' ',' ',ua)
      call print_value('Total unpaired electrons N_D',sum,ndec)
      xnd=sum
      if(idofr.eq.1) call group_by_frag_vec(1,
     +  'FRAGMENT UNPAIRED ELECTRONS',ua,cshort,' ')

!! a=3/4 U decomposition (I. Mayer, P. Salvador local-spin formula) !!
      sum=0.0d0
      do iat=1,natoms
        do jat=1,natoms
          rindex(iat,jat)=0.0d0
        enddo
      enddo

      do iat=1,natoms
        x=0.0d0
        do mu=1,nbasis
          do nu=1,nbasis
            x=x+tts(mu,nu)*tt(nu,mu,iat)
          enddo
        enddo
        rindex(iat,iat)=x*0.750d0
        do ibt=iat,natoms
          x=0.0d0
          do mu=1,nbasis
            do nu=1,nbasis
              x=x+tt(mu,nu,ibt)*tt(nu,mu,iat)
            enddo
          enddo
          rindex(iat,ibt)=rindex(iat,ibt)-x*0.25d0+(qsat(iat,1)*qsat(ibt,1))*0.25d0
          rindex(ibt,iat)=rindex(iat,ibt)
          sum=sum+rindex(iat,ibt)
          if(iat.ne.ibt) sum=sum+rindex(iat,ibt)
        enddo
      enddo

      call print_box('S^2 DECOMPOSITION (a=3/4)')
      call mprint2(rindex,natoms,maxat)
      call print_value('Total <S^2>',sum,ndec)

      do i=1,natoms
        do j=1,natoms
          xlsa(i,j)=rindex(i,j)
        enddo
      enddo
      if(idofr.eq.1) call group_by_frag_mat(0,
     +  'FRAGMENT S^2 DECOMPOSITION (a=3/4)',xlsa)
      call json_local_spin(.false.,xnd,sum,xlsa,xlsa,0.d0,0.d0)

      deallocate(tt,tts)

      return
      end

!! ----------------------------------------------------------------------
!! subroutine: population_density
!! Electron and spin populations of each atom, Tr(P S^A) and Tr(Ps S^A)
!! in /qat/ column 1, and their Mulliken counterparts in column 2. Opens
!! the population analysis section with the atomic definition used.
!! arguments:
!!   sat (in)  (nbasis,nbasis,natoms) atomic overlap matrices S^A
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine population_density(sat)
      use basis_set
      use ao_matrices
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension sat(nbasis,nbasis,natoms)
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort

      call print_box('DOING POPULATION ANALYSIS')
      call atomic_definition_names(ckey,cname,cshort,cspace)
      if(cspace.eq.'hilbert') then
        call print_value('Atomic definition',
     +    trim(ckey)//' (Hilbert space)')
      else
        call print_value('Atomic definition',trim(ckey)//' (real space)')
      end if

!! spin density: P^s*S^A contraction, one atom per outer iteration       !!
      if(kop.ne.0.or.nalf.ne.nb) then

!! parallel over kat (natoms) -- each iteration reduces mu,nu into a     !!
!! private x and writes only its own qsat(kat,1) slot; ps/sat shared     !!
!! read-only, no cross-iteration dependency.                             !!
!$OMP PARALLEL DO PRIVATE(mu,nu,x)
        do kat=1,natoms
          x=0.d0
          do mu=1,nbasis
            do nu=1,nbasis
              x=x+ps(mu,nu)*sat(mu,nu,kat)
            enddo
          enddo
          qsat(kat,1)=x
          qsat(kat,2)=0.d0
        enddo
!$OMP END PARALLEL DO

!! second (Ps*S) contribution, O(igr) outer trips only -- left serial,   !!
!! not worth threading, and kat=ihold(mu) can repeat across mu so the    !!
!! qsat(kat,2) accumulation isn't race-free without extra care.          !!
        do mu=1,nbasis
          kat=ihold(mu)
          x=0.d0
          do itau=1,nbasis
            x=x+ps(mu,itau)*s(itau,mu)
          enddo
          qsat(kat,2)=qsat(kat,2)+x
        enddo
      end if

!! total density: P*S^A contraction, same pattern as spin density above  !!
!$OMP PARALLEL DO PRIVATE(mu,nu,x)
      do kat=1,natoms
        x=0.d0
        do mu=1,nbasis
          do nu=1,nbasis
            x=x+p(mu,nu)*sat(mu,nu,kat)
          enddo
        enddo
        qat(kat,1)=x
        qat(kat,2)=0.d0
      enddo
!$OMP END PARALLEL DO

!! second (P*S) contribution -- same O(igr)-only reasoning as above      !!
      do mu=1,nbasis
        kat=ihold(mu)
        x=0.d0
        do itau=1,nbasis
          x=x+p(mu,itau)*s(itau,mu)
        enddo
        qat(kat,2)=qat(kat,2)+x
      enddo

      return
      end

!! ----------------------------------------------------------------------
!! subroutine: print_population_table
!! Prints a per-atom population table (electrons, charges or spin): one
!! column for the atomic definition used and one for Mulliken (only the
!! Mulliken column when that is the definition), a Total row, and the
!! same per fragment with DOFRAGS.
!! arguments:
!!   title     (in)  section title, e.g. 'ELECTRON POPULATIONS'
!!   fragtitle (in)  title of the fragment table
!!   arr       (in)  (maxat,2) values: atomic definition, Mulliken
!!   tc1,tc2   (in)  column totals
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine print_population_table(title,fragtitle,arr,tc1,tc2)
      use input_options_mod, only: idofr,imulli
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /printout/iaccur
      character*(*) title,fragtitle
      dimension arr(maxat,2)
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort
      character*30 fmthdr,fmttot

      call atomic_definition_names(ckey,cname,cshort,cspace)
      ncol=2
      if(imulli.eq.1) ncol=1

!! atom label 7 characters wide, as in vprint; value columns 12 wide    !!
!! (6 decimals), 20 wide (13) with FULLPRECISION                        !!
      if(iaccur.eq.0) then
        iw=12
        ndec=6
      else
        iw=20
        ndec=13
      end if
      write(fmthdr,'(a,i0,a,i0,a)') '(2x,a,',ncol,'a',iw,')'
      write(fmttot,'(a,i0,a,i0,a,i0,a)') '(2x,a,',ncol,'f',iw,'.',ndec,
     +  ')'

      call print_box(title)
      if(ncol.eq.2) then
        write(*,fmthdr) 'Atom   ',trim(cshort),'Mulliken'
      else
        write(*,fmthdr) 'Atom   ',trim(cshort)
      end if
      write(*,'(2x,a)') repeat('-',7+iw*ncol)
      call vprint(arr,nat,maxat,ncol)
      write(*,'(2x,a)') repeat('-',7+iw*ncol)
      if(ncol.eq.2) then
        write(*,fmttot) 'Total  ',tc1,tc2
      else
        write(*,fmttot) 'Total  ',tc1
      end if

      if(idofr.eq.1) call group_by_frag_vec(ncol,fragtitle,arr,cshort,
     +  'Mulliken')
      end

!! ----------------------------------------------------------------------
!! subroutine: print_atom_table
!! Prints a one-column per-atom table (valences, unpaired electrons),
!! headed by the name of the atomic definition, and the same per
!! fragment with DOFRAGS.
!! arguments:
!!   title     (in)  section title, e.g. 'TOTAL VALENCES'
!!   caption   (in)  line printed under the title, ' ' for none
!!   fragtitle (in)  title of the fragment table, ' ' for none (the
!!                   caller prints it)
!!   arr       (in)  (maxat) values
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine print_atom_table(title,caption,fragtitle,arr)
      use input_options_mod, only: idofr
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /printout/iaccur
      character*(*) title,caption,fragtitle
      dimension arr(maxat)
      character*20 ckey,cspace,fmthdr
      character*60 cname
      character*12 cshort

      call atomic_definition_names(ckey,cname,cshort,cspace)
      call print_box(title)
      if(len_trim(caption).gt.0) then
        write(*,'(2x,a)') trim(caption)
        write(*,'(a)') ''
      end if

!! atom label 7 characters wide, as in vprint; value column 12 wide    !!
!! (6 decimals), 20 wide (13) with FULLPRECISION                       !!
      iw=12
      if(iaccur.ne.0) iw=20
      write(fmthdr,'(a,i0,a)') '(2x,a,a',iw,')'
      write(*,fmthdr) 'Atom   ',trim(cshort)
      write(*,'(2x,a)') repeat('-',7+iw)
      call vprint(arr,nat,maxat,1)
      write(*,'(2x,a)') repeat('-',7+iw)

      if(idofr.eq.1.and.len_trim(fragtitle).gt.0) then
        call group_by_frag_vec(1,fragtitle,arr,cshort,' ')
      end if
      end

!! ----------------------------------------------------------------------
!! subroutine: population_print_charges
!! Prints the electron populations and the atomic charges (nuclear
!! charge minus population) from /qat/, per atom and per fragment.
!! arguments: none (/qat/, /coord/)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine population_print_charges()
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension chg(maxat,2)

      tc1=0.d0
      tc2=0.d0
      do i=1,nat
        tc1=tc1+qat(i,1)
        tc2=tc2+qat(i,2)
      enddo
      call print_population_table('ELECTRON POPULATIONS',
     +  'FRAGMENT ELECTRON POPULATIONS',qat,tc1,tc2)

      tc1=0.d0
      tc2=0.d0
      do i=1,nat
        tc1=tc1-qat(i,1)+zn(i)
        tc2=tc2-qat(i,2)+zn(i)
        chg(i,1)=zn(i)-qat(i,1)
        chg(i,2)=zn(i)-qat(i,2)
      enddo
      call print_population_table('TOTAL ATOMIC CHARGES',
     +  'FRAGMENT CHARGES',chg,tc1,tc2)
      end

!! ----------------------------------------------------------------------
!! subroutine: population_print_overlap
!! Prints the spin populations (open-shell wavefunctions, and correlated
!! ones with a spin density) and, with OPOP, the overlap population
!! matrix, built here (opop in real space, mull_opop for Mulliken).
!! Without OPOP only the diagonal of op is set, to the populations.
!! arguments:
!!   wp,omp,omp2,rho (in)  integration grid and density, for opop
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine population_print_overlap(wp,omp,omp2,rho)
      use basis_set, only: natoms
      use integration_grid
      use input_options_mod, only: idofr,imulli,iopop,icorr
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      dimension wp(nrad*nang*natoms),rho(nrad*nang*natoms)
      dimension omp(nrad*nang*natoms),omp2(nrad*nang*natoms,natoms)

      if(kop.ne.0.OR.(icas.eq.1.and.nalf.ne.nb.and.icorr.ne.0)) then
        tc1=0.d0
        tc2=0.d0
        do i=1,nat
          tc1=tc1+qsat(i,1)
          tc2=tc2+qsat(i,2)
        enddo
        call print_population_table('SPIN POPULATIONS',
     +    'FRAGMENT SPIN POPULATIONS',qsat,tc1,tc2)
      end if

      if(iopop.eq.0) then
        do i=1,nat
          op(i,i)=qat(i,1)
        end do
      else if(imulli.ne.1) then
        call opop(wp,omp,omp2,rho)
      else
        call mull_opop()
      end if

      if(iopop.eq.1) then
        call print_box('OVERLAP POPULATION MATRIX')
        call mprint2(op,nat,maxat)
        if(idofr.eq.1) call group_by_frag_mat(0,
     +    'FRAGMENT OVERLAP POPULATIONS',op)
      end if
      end

!! ----------------------------------------------------------------------
!! subroutine: json_populations
!! Writes the populations block of the results file: electron
!! populations, charges and (when printed) spin populations of the atoms
!! and fragments, for the atomic definition used and, as a comparison,
!! for Mulliken.
!! arguments: none (/qat/, /coord/)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_populations()
      use input_options_mod, only: idofr,imulli,icorr
      use json_mod
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /printout/iaccur
      dimension chg(maxat,2)
      character*20 ckey,cspace,crole
      character*60 cname
      character*12 cshort,cs
      logical lspin

      call atomic_definition_names(ckey,cname,cshort,cspace)
      ndec=6
      if(iaccur.ne.0) ndec=13
!! spin populations exactly when population_print_overlap prints them !!
      lspin=kop.ne.0.or.(icas.eq.1.and.nalf.ne.nb.and.icorr.ne.0)
      do i=1,nat
        chg(i,1)=zn(i)-qat(i,1)
        chg(i,2)=zn(i)-qat(i,2)
      end do

!! set 1: the atomic definition used; set 2: Mulliken, in its own     !!
!! object and marked as a comparison, unless Mulliken is the atomic    !!
!! definition                                                          !!
      nset=2
      if(imulli.eq.1) nset=1
      call json_begin_object('populations')
      call json_value('scheme',trim(ckey))
      do iset=1,nset
        cs=cshort
        crole=' '
        if(iset.eq.2) then
          cs='Mulliken'
          crole='comparison'
          call json_begin_object('mulliken')
        end if
        call json_begin_object('atoms')
        call json_vector('electrons',qat(1:nat,iset),nat)
        call json_index('electrons','atom_vector',
     +    'Electron populations ('//trim(cs)//')',' ',ndec,' ',.false.,
     +    crole)
        call json_vector('charges',chg(1:nat,iset),nat)
        call json_index('charges','atom_vector',
     +    'Atomic charges ('//trim(cs)//')',' ',ndec,' ',.true.,crole)
        if(lspin) then
          call json_vector('spin',qsat(1:nat,iset),nat)
          call json_index('spin','atom_vector',
     +      'Spin populations ('//trim(cs)//')',' ',ndec,' ',.true.,
     +      crole)
        end if
        call json_end()
        if(idofr.eq.1) then
          call json_begin_object('fragments')
          call json_frag_vector('electrons',qat(1,iset))
          call json_index('electrons','fragment_vector',
     +      'Fragment electron populations ('//trim(cs)//')',' ',ndec,
     +      ' ',.false.,crole)
          call json_frag_vector('charges',chg(1,iset))
          call json_index('charges','fragment_vector',
     +      'Fragment charges ('//trim(cs)//')',' ',ndec,' ',.true.,
     +      crole)
          if(lspin) then
            call json_frag_vector('spin',qsat(1,iset))
            call json_index('spin','fragment_vector',
     +        'Fragment spin populations ('//trim(cs)//')',' ',ndec,' ',
     +        .true.,crole)
          end if
          call json_end()
        end if
        if(iset.eq.2) call json_end()
      end do
      call json_end()
      end

!! ----------------------------------------------------------------------
!! subroutine: bond_order_analysis
!! Bond orders and valences (fborder), the fragment bond orders with
!! DOFRAGS, and their block in the results file.
!! arguments:
!!   sat (in)  (nbasis,nbasis,natoms) atomic overlap matrices S^A
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine bond_order_analysis(sat)
      use basis_set, only: natoms,nbasis
      use input_options_mod, only: idofr
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      dimension sat(nbasis,nbasis,natoms),val(maxat,3)

      call fborder(sat,val)
      if(idofr.eq.1) call group_by_frag_mat(1,'FRAGMENT BOND ORDERS',bo)
      call json_bond_orders(val)
      end

!! ----------------------------------------------------------------------
!! subroutine: json_bond_orders
!! Writes the bond_orders block of the results file: the bond order
!! matrix and the valences of the atoms and fragments, and the overlap
!! populations with OPOP.
!! arguments:
!!   val (in)  (maxat,3) total, bonded and free valences (fborder)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_bond_orders(val)
      use input_options_mod, only: idofr,iopop,imulli
      use json_mod
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /printout/iaccur
      dimension val(maxat,3)
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort
      character*20 cval(3)
      character*40 tval(3),tfrag(3)
      logical lsignop

      call atomic_definition_names(ckey,cname,cshort,cspace)
      ndec=6
      if(iaccur.ne.0) ndec=13
      cval(1)='valence_total'
      cval(2)='valence_bonded'
      cval(3)='valence_free'
      tval(1)='Total valences'
      tval(2)='Valences used in bonds'
      tval(3)='Free valences'
      tfrag(1)='Fragment total valences'
      tfrag(2)='Fragment valences used in bonds'
      tfrag(3)='Fragment free valences'
!! Mulliken bond orders and overlap populations can be negative !!
      lsignop=imulli.eq.1

      call json_begin_object('bond_orders')
      call json_value('scheme',trim(ckey))
      call json_begin_object('atoms')
      call json_atom_matrix('matrix',bo,.true.,'localization index',
     +  'bond order',-1,0.d0)
      call json_index('matrix','atom_matrix',
     +  'Bond orders ('//trim(cshort)//')',' ',ndec,' ',lsignop,' ')
      do k=1,3
        call json_vector(trim(cval(k)),val(1:nat,k),nat)
        call json_index(trim(cval(k)),'atom_vector',
     +    trim(tval(k))//' ('//trim(cshort)//')',' ',ndec,' ',.false.,
     +    ' ')
      end do
      call json_end()

      if(idofr.eq.1) then
        call json_begin_object('fragments')
        call json_frag_matrix('matrix',1,bo,.true.,
     +    'localization index','bond order')
        call json_index('matrix','fragment_matrix',
     +    'Fragment bond orders ('//trim(cshort)//')',' ',ndec,' ',
     +    lsignop,' ')
        do k=1,3
          call json_frag_vector(trim(cval(k)),val(1,k))
          call json_index(trim(cval(k)),'fragment_vector',
     +      trim(tfrag(k))//' ('//trim(cshort)//')',' ',ndec,' ',
     +      .false.,' ')
        end do
        call json_end()
      end if

      if(iopop.eq.1) then
        call json_begin_object('overlap_populations')
        call json_begin_object('atoms')
        call json_atom_matrix('matrix',op,.true.,
     +    'self overlap population','overlap population',-1,0.d0)
        call json_index('matrix','atom_matrix',
     +    'Overlap populations ('//trim(cshort)//')',' ',ndec,' ',
     +    lsignop,' ')
        call json_end()
        if(idofr.eq.1) then
          call json_begin_object('fragments')
          call json_frag_matrix('matrix',0,op,.true.,
     +      'self overlap population','overlap population')
          call json_index('matrix','fragment_matrix',
     +      'Fragment overlap populations ('//trim(cshort)//')',' ',ndec,
     +      ' ',lsignop,' ')
          call json_end()
        end if
        call json_end()
      end if
      call json_end()
      end

!! ----------------------------------------------------------------------
!! subroutine: json_local_spin
!! Writes the local_spin block of the results file: effectively unpaired
!! electrons and the <S^2> decomposition of the atoms and fragments
!! (/localspin/), with, for a correlated wavefunction, Mayer's bond
!! orders from the natural orbitals and the LI/DI matrix from the 2-RDM.
!! arguments:
!!   lcorr (in)  logical, correlated wavefunction (spincorr)
!!   xnd   (in)  total number of effectively unpaired electrons
!!   xs2   (in)  total <S^2> of the decomposition
!!   bonew (in)  (maxat,maxat) Mayer's bond orders (lcorr only)
!!   di2   (in)  (maxat,maxat) LI/DI matrix (lcorr only)
!!   xbonew,xdi (in)  their totals, diagonal and one triangle
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_local_spin(lcorr,xnd,xs2,bonew,di2,xbonew,xdi)
      use input_options_mod, only: idofr
      use json_mod
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /localspin/ xlsa(maxat,maxat),ua(maxat)
      common /printout/iaccur
      dimension bonew(maxat,maxat),di2(maxat,maxat)
      logical lcorr
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort

      call atomic_definition_names(ckey,cname,cshort,cspace)
      ndec=6
      if(iaccur.ne.0) ndec=13

      call json_begin_object('local_spin')
      call json_value('scheme',trim(ckey))
      if(lcorr) then
        call json_value('wavefunction','correlated')
      else
        call json_value('wavefunction','single determinant')
      end if
      call json_value('total_unpaired_electrons',xnd)
      call json_index('total_unpaired_electrons','scalar',
     +  'Total unpaired electrons N_D ('//trim(cshort)//')',' ',ndec,
     +  ' ',.false.,' ')
      call json_value('total_s2',xs2)
      call json_index('total_s2','scalar',
     +  'Total <S^2> ('//trim(cshort)//')',' ',ndec,' ',.false.,' ')

      call json_begin_object('atoms')
      call json_vector('unpaired_electrons',ua(1:nat),nat)
      call json_index('unpaired_electrons','atom_vector',
     +  'Effectively unpaired electrons ('//trim(cshort)//')',' ',ndec,
     +  ' ',.false.,' ')
      call json_atom_matrix('s2',xlsa,.true.,'local spin',
     +  'spin coupling',-1,0.d0)
      call json_index('s2','atom_matrix',
     +  'S^2 decomposition ('//trim(cshort)//')',' ',ndec,' ',.true.,
     +  ' ')
      if(lcorr) then
        call json_atom_matrix('new_bond_orders',bonew,.true.,
     +    'localization index','bond order',1,xbonew)
        call json_index('new_bond_orders','atom_matrix',
     +    'New bond orders ('//trim(cshort)//')',' ',ndec,' ',.false.,
     +    ' ')
        call json_atom_matrix('li_di',di2,.true.,'localization index',
     +    'delocalization index',1,xdi)
        call json_index('li_di','atom_matrix',
     +    'LI/DI from the 2-RDM ('//trim(cshort)//')',' ',ndec,' ',
     +    .false.,' ')
      end if
      call json_end()

      if(idofr.eq.1) then
        call json_begin_object('fragments')
        call json_frag_vector('unpaired_electrons',ua)
        call json_index('unpaired_electrons','fragment_vector',
     +    'Fragment unpaired electrons ('//trim(cshort)//')',' ',ndec,
     +    ' ',.false.,' ')
        call json_frag_matrix('s2',0,xlsa,.true.,'local spin',
     +    'spin coupling')
        call json_index('s2','fragment_matrix',
     +    'Fragment S^2 decomposition ('//trim(cshort)//')',' ',ndec,
     +    ' ',.true.,' ')
        if(lcorr) then
          call json_frag_matrix('li_di',1,di2,.true.,
     +      'localization index','delocalization index')
          call json_index('li_di','fragment_matrix',
     +      'Fragment LI/DI from the 2-RDM ('//trim(cshort)//')',' ',
     +      ndec,' ',.false.,' ')
        end if
        call json_end()
      end if
      call json_end()
      end
