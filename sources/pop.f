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
!! Theory Comput. 2012, 8, 1270), and the "DAVIDSON SPIN DEC." matrix,
!! meant as the Clark-Davidson decomposition of the same paper but whose
!! formula does not match it yet (under review). spincorr (corr.f) does
!! the same for correlated wavefunctions.
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

!! "DAVIDSON SPIN DEC." matrix (meant as Clark-Davidson, formula under   !!
!! review, see the header) -- rindex is safely reused here, xlsa already !!
!! holds the a=3/4 result saved above. P*S^A products first, then the Ps !!
!! correction (open-shell only).                                         !!
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
        xx0=0.d0
        do ibt=1,natoms
          x=0.d0
          do mu=1,nbasis
            do nu=1,nbasis
              x=x+tt(mu,nu,iat)*tt(nu,mu,ibt)
            enddo
          enddo
          if(iat.ne.ibt) then
            rindex(iat,ibt)=rindex(iat,ibt)-3.0d0/8.0d0*x
            xx0=xx0+x
          endif
        enddo
        rindex(iat,iat)=rindex(iat,iat)+3.0d0/8.0d0*xx0
      enddo

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
          xx0=0.d0
          do ibt=1,natoms
            x=0.d0
            do mu=1,nbasis
              do nu=1,nbasis
                x=x+tt(mu,nu,iat)*tt(nu,mu,ibt)
              enddo
            enddo
            if(iat.ne.ibt) then
              rindex(iat,ibt)=rindex(iat,ibt)-3.0d0/8.0d0*x
              xx0=xx0+x
            endif
          enddo
          rindex(iat,iat)=rindex(iat,iat)+3.0d0/8.0d0*xx0
        enddo
      endif

      call print_box('"FUZZY ATOMS" DAVIDSON SPIN DEC. MATRIX')
      call mprint(rindex,natoms,maxat)
      sum=0.0d0
      do iat=1,natoms
        do ibt=1,natoms
          sum=sum+rindex(iat,ibt)
        enddo
      enddo
      write(*,*)
      write(*,'(2x,a,f10.5)') 'Sum check <S^2> = ',sum

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
!! subroutine: bond_order_analysis
!! Bond orders and valences (fborder), and the fragment bond orders with
!! DOFRAGS.
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
      end

!! ----------------------------------------------------------------------
!! subroutine: pca_analysis
!! Principal component analysis of the covariance matrix op - di/2
!! (overlap populations and bond orders): prints the matrix, its
!! eigenvectors (columns) and eigenvalues, and the projection of the
!! atomic populations on each component. Needs op and di (fborder).
!! Under discussion (what it should print): not documented, not in the
!! results file.
!! arguments: none (/ovpop/, /qat/)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine pca_analysis()
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      allocatable pca(:,:),scr(:,:)

      allocate(pca(nat,nat))
      allocate(scr(nat,nat))
      do i=1,nat
        do j=1,nat
          pca(i,j)=op(i,j)-0.50d0*di(i,j)
        end do
      end do

      call print_box('"FUZZY ATOMS" COVARIANCE MATRIX')
      call mprint2(pca,nat,nat)

      call diagonalize(nat,nat,pca,scr,0)

      call print_box('"FUZZY ATOMS" PCA EIGENVECTORS')
      call mprint2(scr,nat,nat)
      write(*,*)
      write(*,'(2x,a)') 'Eigenvalues:'
      write(*,'(2x,8f10.4)') (pca(i,i),i=1,nat)

!! the second column multiplies by scr(i,1), component i of the first !!
!! eigenvector (not eigenvalue i): kept until its intent is settled    !!
      call print_box('PCA SUM CHECKS')
      write(*,'(2x,a)') '  PC   Charge projection    x eigenvector 1'
      write(*,'(2x,a)') repeat('-',44)
      do i=1,nat
        xx=ZERO
        do k=1,nat
          xx=xx+pca(k,i)*qat(k,1)
        end do
        write(*,'(3x,i3,2f18.6)') i,xx,xx*scr(i,1)
      end do
      write(*,'(2x,a)') repeat('-',44)

      deallocate(pca,scr)

      return
      end

