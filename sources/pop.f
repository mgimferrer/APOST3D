!! ********************************************************************* !!
!! subroutine: fborder                                                   !!
!! purpose: builds the "fuzzy atoms" bond order matrix (bo/di, via the   !!
!!   R-index trace of P*S^A products, Ps-corrected for open-shell) and   !!
!!   the three derived valence tables (total/used-in-bonds/free).        !!
!!   called from bond_order_analysis for every AIM scheme -- sat already !!
!!   encodes which one (numint_sat/tomull/tolow/...).                    !!
!! arguments:                                                            !!
!!   sat (in) -- per-atom AO overlap matrix                              !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine fborder(sat)
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/ op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/ qat(maxat,2),qsat(maxat,2)
      dimension sat(nbasis,nbasis,natoms)

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
!! diag(iat) is captured here from the P-only R-index, before the Ps    !!
!! correction below is folded into rindex -- TOTAL VALENCES (built from !!
!! diag, see the tail of this routine) is therefore P-only even for     !!
!! open-shell (kop!=0), while the bond order matrix/valences-in-bonds   !!
!! (built from rindex) do get the Ps correction. Possibly inconsistent  !!
!! for open-shell systems; not changed pending confirmation.            !!
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

!! save the bond order matrix -- di gets overwritten later for a !!
!! correlated (CAS/CISD) calculation. !!
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

      call print_box('"FUZZY ATOMS" BOND ORDER MATRIX')
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

      call print_valence_table('TOTAL VALENCES',' ','Total valences',
     +  'V_A',diag)
      call print_valence_table('VALENCES USED IN BONDS',
     +  '(SUM OF BOND ORDERS)','Valences used in bonds','VB_A',tindex)
      do i=1,natoms
        diag(i)=diag(i)-tindex(i)
      enddo
      call print_valence_table('FREE VALENCES',' ','Free valences',
     +  'F_A',diag)

      deallocate(tt)

      return
      end

      

!! ********************************************************************* !!
!! subroutine: opop                                                      !!
!! purpose: builds the atom-pair overlap population matrix (op) by       !!
!!   numerical integration of the density weighted by each atom's        !!
!!   fuzzy-partition weight -- the default scheme (iallpo=0) integrates  !!
!!   each atom's own grid; ALLPOINTS (iallpo=1) instead pre-weights rho   !!
!!   once and re-sums it over the full grid for every atom pair.         !!
!! arguments:                                                            !!
!!   wp        (in) -- integration weight of each grid point             !!
!!   omp       (in) -- becke weight of each grid point (ALLPOINTS only)  !!
!!   omp2      (in) -- per-atom fuzzy weight of each grid point          !!
!!   rho       (in) -- electron density at each grid point (ALLPOINTS    !!
!!                      rescales this array in place)                    !!
!! author:                                                                !!
!! ********************************************************************* !!
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

!! ********************************************************************* !!
!! subroutine: fspindec                                                  !!
!! purpose: "fuzzy atoms" local-spin decomposition for a single-         !!
!!   determinant density (the SPIN keyword's default path -- spincorr    !!
!!   in corr.f is the sibling correlated-wavefunction, icas/icisd        !!
!!   branch). Prints effectively unpaired electrons (u_A), the a=3/4     !!
!!   <S^2> decomposition (I. Mayer, P. Salvador), and its Davidson-      !!
!!   Lowdin-basis twin.                                                  !!
!! arguments:                                                            !!
!!   sat (in) -- per-atom AO overlap matrix (numint_sat/tomull/tolow)    !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine fspindec(sat)
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
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

!! parallelization: not done. The O(nbasis^2*natoms^2) iat/ibt/mu/nu     !!
!! loops below (a=3/4 and Davidson decompositions) are the dominant     !!
!! cost, same profile as spincorr's in corr.f -- no codebase precedent  !!
!! yet for an array-accumulate REDUCTION at this size, and every active !!
!! test system here is small enough that it isn't currently a           !!
!! bottleneck.                                                          !!

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

      call print_box('EFFECTIVELY UNPAIRED ELECTRONS')
      write(*,'(1x,a7,a12)') 'Atom','u_A'
      write(*,'(2x,a)') repeat('-',18)
      call vprint(ua,nat,maxat,1)
      write(*,'(2x,a)') repeat('-',18)
      write(*,'(2x,a,f10.5)') 'Sum check N_D = ',sum

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

      call print_box('"FUZZY ATOMS" S^2 DECOMPOSITION (a=3/4)')
      call mprint(rindex,natoms,maxat)
      write(*,*)
      write(*,'(2x,a,f10.5)') 'Sum check <S^2> = ',sum

      do i=1,natoms
        do j=1,natoms
          xlsa(i,j)=rindex(i,j)
        enddo
      enddo

!! Davidson-Lowdin-basis twin -- rindex is safely reused here, xlsa      !!
!! already holds the a=3/4 result saved above. P*S^A products first,    !!
!! then the Ps correction (open-shell only).                            !!
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

!! ********************************************************************* !!
!! subroutine: population_density                                        !!
!! purpose: spin- and total-density atomic-population contraction        !!
!! (P*S^A products), writing qat/qsat via COMMON. Split from the print   !!
!! side (population_print_charges/population_print_overlap) so the       !!
!! idoint-guarded print_int call in main.f can sit between them,         !!
!! exactly where it did before this routine existed.                     !!
!! arguments:                                                            !!
!! sat (in) -- per-atom AO overlap matrix (numint_sat/tomull/tolow)      !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine population_density(sat)
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension sat(nbasis,nbasis,natoms)

      call print_box('DOING POPULATION ANALYSIS')
      write(*,'(2x,a)') 'Partial atomic charges'
      write(*,'(2x,a)') 'Atomic spin densities'
      write(*,'(2x,a)') 'Bond orders and Valences'

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

!! ********************************************************************* !!
!! subroutine: print_population_table                                    !!
!! purpose: shared title/column-header/data/sum/fragment-breakdown       !!
!! layout for the three 2-column (apost3d, Mulliken) atomic tables       !!
!! (electron populations, atomic charges, spin populations) -- replaces  !!
!! three near-identical copies of this block, previously inconsistent    !!
!! in small ways (a stray "apost3D" instead of "apost3d" in one of the   !!
!! three, an extra blank line before the fragment breakdown in two of    !!
!! the three but not the first). Header/separator/sum-line widths are    !!
!! sized to match vprint's own fixed data-row format exactly (32 columns !!
!! for the default 6-decimal case, 48 for FULLPRECISION) -- previously a  !!
!! hardcoded 29/35-column layout that didn't cover the actual numbers.    !!
!! arguments:                                                            !!
!! title    (in) -- print_box title, e.g. 'ELECTRON POPULATIONS'         !!
!! fraglabel(in) -- fragment-breakdown label, e.g. 'Electron populations'!!
!! arr      (in) -- the (maxat,2) table to print (qat/dummyvec/qsat)     !!
!! tc1,tc2  (in) -- apost3d/Mulliken column sums, computed by the caller !!
!!                  (the summation itself differs per table)             !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine print_population_table(title,fraglabel,arr,tc1,tc2)
      use input_options_mod, only: idofr
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /printout/iaccur
      character*(*) title,fraglabel
      dimension arr(maxat,2)
      character*80 line

      call print_box(title)
      if(iaccur.eq.0) then
        write(*,'(1x,a7,2a12)') 'Atom','apost3d','Mulliken'
        write(*,'(2x,a)') repeat('-',30)
      else
        write(*,'(1x,a7,2a20)') 'Atom','apost3d','Mulliken'
        write(*,'(2x,a)') repeat('-',46)
      end if
      call vprint(arr,nat,maxat,2)
      if(iaccur.eq.0) then
        write(*,'(2x,a)') repeat('-',30)
        write(*,162) tc1,tc2
      else
        write(*,'(2x,a)') repeat('-',46)
        write(*,172) tc1,tc2
      end if

      if (idofr.eq.1) then
        line ='   FRAGMENT ANALYSIS : '//fraglabel
        call group_by_frag_vec(2,line,arr)
      end if

      return

 162  format(1x,'    Sum',2f12.6)
 172  format(1x,'    Sum',2f20.13)

      end

!! ********************************************************************* !!
!! subroutine: print_valence_table                                       !!
!! purpose: shared title/column-header/data layout for the three         !!
!! single-column atomic valence tables printed by fborder (V_A/VB_A/     !!
!! F_A) -- same consolidation print_population_table already applies to  !!
!! the electron/charge/spin population tables, and for the same reason:  !!
!! the header/separator width now matches vprint's own fixed data-row    !!
!! format (20 columns) instead of a hardcoded value that didn't cover    !!
!! the numbers. Also prints a per-fragment breakdown (DOFRAGS) via        !!
!! group_by_frag_vec, matching the population/bond-order tables --        !!
!! previously the only ones of the six population-analysis tables         !!
!! without one.                                                           !!
!! arguments:                                                            !!
!! title    (in) -- print_box title, e.g. 'TOTAL VALENCES'               !!
!! subtitle (in) -- optional explanatory line under the title, blank     !!
!!                   (' ') if none, e.g. '(SUM OF BOND ORDERS)'          !!
!! fraglabel(in) -- fragment-breakdown label, e.g. 'Total valences'      !!
!! label    (in) -- column header, e.g. 'V_A'                            !!
!! arr      (in) -- the (maxat) vector to print (diag/tindex)            !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine print_valence_table(title,subtitle,fraglabel,label,arr)
      use input_options_mod, only: idofr
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      character*(*) title,subtitle,fraglabel,label
      dimension arr(maxat)
      character*80 line

      call print_box(title)
      if(len_trim(subtitle).gt.0) then
        write(*,'(1x,a)') trim(subtitle)
        write(*,*)
      end if
      write(*,'(1x,a7,a12)') 'Atom',label
      write(*,'(2x,a)') repeat('-',18)
      call vprint(arr,nat,maxat,1)
      write(*,'(2x,a)') repeat('-',18)

      if (idofr.eq.1) then
        line=' FRAGMENT ANALYSIS: '//fraglabel
        call group_by_frag_vec(1,line,arr)
      end if

      return
      end

!! ********************************************************************* !!
!! subroutine: population_print_charges                                  !!
!! purpose: prints electron populations and partial atomic charges       !!
!! (plus per-fragment breakdowns where DOFRAGS is set) from qat, already !!
!! filled by population_density. The ELCOUNT/NCTAIM call that used to    !!
!! sit between these two prints stays in main.f (like idoint/print_int)  !!
!! -- NCTAIM lives in subroutines_mmo.f, which apost3d-eos's link list   !!
!! deliberately excludes, so pop.f itself must not reference it.         !!
!! arguments: none (all via qat/COMMON)                                  !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine population_print_charges()
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension dummyvec(maxat,2)

!! electron populations -- straight sum of qat !!
      tc1=0.d0
      tc2=0.d0
      do i=1,nat
        tc1=tc1+qat(i,1)
        tc2=tc2+qat(i,2)
      enddo
      call print_population_table('ELECTRON POPULATIONS',
     +  'Electron populations',qat,tc1,tc2)

!! partial charges -- nuclear charge minus qat, atom by atom !!
      tc1=0.d0
      tc2=0.d0
      do i=1,nat
        tc1=tc1-qat(i,1)+zn(i)
        tc2=tc2-qat(i,2)+zn(i)
        dummyvec(i,1)=zn(i)-qat(i,1)
        dummyvec(i,2)=zn(i)-qat(i,2)
      enddo
      call print_population_table('TOTAL ATOMIC CHARGES',
     +  'Atomic Charges',dummyvec,tc1,tc2)

      return
      end

!! ********************************************************************* !!
!! subroutine: population_print_overlap                                  !!
!! purpose: prints spin populations (plus per-fragment breakdown) and    !!
!! overlap populations from qsat/op, already filled by                   !!
!! population_density -- op itself is built here via opop/mull_opop.     !!
!! arguments:                                                            !!
!! wp,omp,omp2,rho (in) -- integration grid weights/density, passed      !!
!! straight to opop                                                      !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
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
      character*80 line

!! spin populations -- straight sum of qsat, only for open-shell/       !!
!! correlated-with-DM cases                                              !!
      if(kop.ne.0.OR.(icas.eq.1.and.nalf.ne.nb.and.icorr.ne.0)) then
        tc1=0.d0
        tc2=0.d0
        do i=1,nat
          tc1=tc1+qsat(i,1)
          tc2=tc2+qsat(i,2)
        enddo
        call print_population_table('SPIN POPULATIONS',
     +    'Spin Populations',qsat,tc1,tc2)
      end if

!! overlap populations -- op itself comes from opop (real-space) or     !!
!! mull_opop (Mulliken); iopop=0 just takes the diagonal from qat        !!
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
        call print_box('APOST3D OVERLAP POPULATION MATRIX')
        call mprint2(op,nat,maxat)
        if (idofr.eq.1) then
          line ='   FRAGMENT ANALYSIS : Overlap Populations'
          call group_by_frag_mat(0,line,op)
        end if
      end if

      return
      end

!! ********************************************************************* !!
!! subroutine: bond_order_analysis                                       !!
!! purpose: thin wrapper around fborder -- kept separate from            !!
!! population_density/population_print_charges/population_print_overlap  !!
!! on purpose, so a future second bond-order scheme gets its own         !!
!! equally-thin subroutine instead of being wedged into population math. !!
!! arguments:                                                            !!
!! sat (in) -- per-atom AO overlap matrix, passed straight to fborder    !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine bond_order_analysis(sat)
      use basis_set, only: natoms,nbasis
      use input_options_mod, only: idofr
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      dimension sat(nbasis,nbasis,natoms)
      character*80 line

      call fborder(sat)
      if (idofr.eq.1) then
        line ='   FRAGMENT ANALYSIS : Fuzzy Bond Order'
        call group_by_frag_mat(1,line ,bo)
      end if

      return
      end

!! ********************************************************************* !!
!! subroutine: pca_analysis                                              !!
!! purpose: "fuzzy atoms" PCA -- diagonalizes a covariance-like matrix   !!
!! built from the overlap population (op) and delocalization index (di)  !!
!! matrices, printing eigenvectors/eigenvalues and their projection      !!
!! against qat. Needs di already populated (fborder, called earlier via  !!
!! bond_order_analysis).                                                 !!
!! FIXED 2026-08-15 (M. Gimferrer): scr, the eigenvector-output argument !!
!! to diagonalize, was allocated as scr(nat) -- a vector -- but          !!
!! diagonalize needs a full (M,M) matrix there, same as every other call !!
!! site in this codebase. Wrote past scr's allocation whenever nat>1     !!
!! (heap overflow, never caught -- PCA has zero test coverage). Now      !!
!! allocated scr(nat,nat). Confirmed output-preserving: diagonalize's    !!
!! write order fills scr's first nat elements (column 1) before any      !!
!! out-of-bounds write happens, and every print below only ever read     !!
!! scr(1:nat) -- i.e. column 1 -- so the values printed are identical    !!
!! before and after, only the undefined-behavior overflow is gone.       !!
!! STILL OPEN, not touched -- needs M. Gimferrer + P. Salvador to        !!
!! confirm intent first: the two prints right after diagonalize look     !!
!! swapped relative to its actual contract (pca/A0 comes back with       !!
!! eigenvalues on its diagonal, scr/X holds the eigenvectors) -- but the !!
!! code labels pca "PCA EIGENVECTORS" and prints scr as if it held       !!
!! eigenvalues.                                                          !!
!! arguments: none (all via op/di/qat/COMMON)                            !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
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
      call mprint(pca,nat,nat)
      write(*,*)

      call diagonalize(nat,nat,pca,scr,0)

      call print_box('"FUZZY ATOMS" PCA EIGENVECTORS')
      call mprint(pca,nat,nat)
      write(*,*)
      write(*,'(8f10.4)') (scr(i,1),i=1,nat)
      write(*,*)

      do i=1,nat
        xx=ZERO
        do k=1,nat
          xx=xx+pca(k,i)*qat(k,1)
        end do
        write(*,'(2x,a,i3,a,2f14.6)') 'PC: ',i,' sum: ',xx,xx*scr(i,1)
      end do

      deallocate(pca,scr)

      return
      end

