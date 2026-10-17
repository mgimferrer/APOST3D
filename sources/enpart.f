!! ********************************************************** !!
!! IQA/ENPART ENERGY PARTITIONING SUBROUTINES                 !!
!! one-electron part:                                         !!
!!   numint_one      -- RHF/RKS                               !!
!!   numint_one_uhf  -- UHF/UKS                               !!
!! two-electron part (both call calc_coul/multipolar_mo       !!
!! internally for the Coulomb/multipolar-expansion terms):    !!
!!   numint_two      -- RHF/RKS                               !!
!!   numint_two_uhf  -- UHF/UKS                               !!
!!   calc_coul       -- Coulomb energy term                   !!
!! multipolar-expansion xc estimate for skipped atom pairs:   !!
!!   multipolar      -- stored pair functions (CASSCF/CISD)   !!
!!   multipolar_mo   -- MO products, formed on the fly        !!
!!   mp_init/mp_moments/mp_pairs/mp_sum -- shared pieces      !!
!! polarizability analysis (# POLAR keyword, standalone --    !!
!! not part of the one/two-electron machinery above):         !!
!!   polar                                                    !!
!! ********************************************************** !!

!! ***** !!

!! *********************************************************************** !!
!! subroutine: numint_one                                                  !!
!! purpose: one-electron IQA/ENPART energy partition -- electron-nuclear   !!
!!   attraction, kinetic energy, nuclear repulsion, and (if present) the   !!
!!   external-field dipole terms, each decomposed into atomic/diatomic     !!
!!   contributions. RHF/RKS only -- see numint_one_uhf for open-shell.     !!
!! arguments:                                                              !!
!!   ndim   (in)  -- number of basis functions (leading dim of chp)        !!
!!   itotps (in)  -- total number of grid points                           !!
!!   wp     (in)  -- integration weight of each grid point                 !!
!!   rho    (in)  -- electron density at each grid point                   !!
!!   omp    (in)  -- becke/tfvc weight of each grid point for its own atom !!
!!   omp2   (in)  -- becke/tfvc (or hirshfeld) weight of each point for    !!
!!                   every atom                                            !!
!!   pcoord (in)  -- xyz coordinates of each grid point                    !!
!!   chp    (in)  -- basis-function values at each grid point              !!
!!   eto    (out) -- total one-electron atomic/diatomic energy matrix      !!
!! author: PSalse, MGimf, MMO.                                             !!
!! *********************************************************************** !!
      subroutine numint_one(ndim,itotps,wp,rho,omp,omp2,pcoord,chp,eto)

      use ao_matrices
      use integration_grid
      use output_mod, only: print_value

      implicit real*8(a-h,o-z)

      include 'parameter.h'

      common /nat/nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(200)
      common /actual/iact,jat,icenter
      common /energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common /energ0/ekin0,eelnuc0,evee0,etot0
      common /exchg/exch(maxat,maxat),xmix
      common /efield/field(4),edipole

      integer, intent(in) :: ndim,itotps

      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3)
      dimension omp(itotps),omp2(itotps,nat),rho(itotps)
      dimension eto(maxat,maxat)
      dimension epa(maxat,maxat),ekin(maxat,maxat),dip(maxat,3),enuc(maxat,maxat)
      dimension eecp(maxat),ect(maxat,maxat)
      dimension edip(maxat,maxat)
      dimension tkin(maxat),fdip(maxat,2)
      character*40 cfld
      logical lref
      common /printout/iaccur

      allocatable :: chp2(:,:),scr(:),chpaux(:,:)

      idofr   = iopt(40)
      ifield  = iopt(47)
      ipoints = itotps
      iatps   = nang*nrad
      xxdip   = ZERO
      ndec=7
      if(iaccur.ne.0) ndec=13

      ALLOCATE(chp2(itotps,nocc),scr(itotps))

      call ao_to_mo_grid(itotps,igr,nocc,c,chp,chp2)

      epa=ZERO
      ekin=ZERO
      enuc=ZERO
      eto=ZERO
      ect=ZERO
      edip=ZERO

!! pseudopotential atomic energies (if present), added to the E-N terms !!
      call mga_misc(iecp,eecp)
      if(iecp.eq.1) then
        do ii=1,nat
          epa(ii,ii)=eecp(ii)
        end do
      end if

!! field(:)/ifield were already read by main.f into the /efield/ common  !!
!! and iopt(47) -- no need to re-scan the .inp file here.                !!
      if(ifield.eq.1) then
!! parallel over icenter: each iteration writes only its own dip(icenter,:)/ !!
!! ect(icenter,icenter); x/y/z/zztop/xtot are per-iteration scratch, not     !!
!! running totals, so no reduction is needed. Untested: no test uses an      !!
!! external electric field.                                                  !!
!$OMP PARALLEL DO PRIVATE(icenter,x,y,z,zztop,ifut,distx,disty,distz,xtot)
        do icenter=1,nat
          x=ZERO
          y=ZERO
          z=ZERO
          zztop=ZERO
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            distx=pcoord(ifut,1)
            disty=pcoord(ifut,2)
            distz=pcoord(ifut,3)
            x=x+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*distx
            y=y+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*disty
            z=z+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*distz
            zztop=zztop+wp(ifut)*omp2(ifut,icenter)*rho(ifut)
          end do
          dip(icenter,1)=-x+zztop*coord(1,icenter)
          dip(icenter,2)=-y+zztop*coord(2,icenter)
          dip(icenter,3)=-z+zztop*coord(3,icenter)
          xtot=(coord(1,icenter)*field(2)+coord(2,icenter)*field(3)+coord(3,icenter)*field(4))
          ect(icenter,icenter)=-zztop*xtot
        end do
!$OMP END PARALLEL DO

!! electric field x,y,z components are on field(2-4) !!
        xtot=ZERO
        do i=1,nat
          edip(i,i)=+dip(i,1)*field(2)+dip(i,2)*field(3)+dip(i,3)*field(4)
          xtot=xtot+ect(i,i)+edip(i,i)
        end do
        xxdip=xtot
      end if

!! electron-nuclei attraction integrals; parallel over the (icenter,jcenter) !!
!! atom-pair space (COLLAPSE(2), both bounds constant) -- each pair writes  !!
!! only its own epa(icenter,jcenter), xtot is a genuine running total.      !!
      xtot=ZERO
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,jcenter,x,zz,dx0,dy0,dz0,ifut,distx,disty,distz,rr) REDUCTION(+:xtot)
      do icenter=1,nat
        do jcenter=1,nat
          x=ZERO
          zz=zn(jcenter)
          dx0=coord(1,jcenter)
          dy0=coord(2,jcenter)
          dz0=coord(3,jcenter)
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            distx=pcoord(ifut,1)-dx0
            disty=pcoord(ifut,2)-dy0
            distz=pcoord(ifut,3)-dz0
            rr=dsqrt(distx*distx+disty*disty+distz*distz)
            if(rr.gt.1.0d-10) then
              x=x+zz*wp(ifut)*omp2(ifut,icenter)*rho(ifut)/rr
            end if
          end do
          xtot=xtot+x
          epa(icenter,jcenter)=epa(icenter,jcenter)-x
        end do
      end do
!$OMP END PARALLEL DO

      xtot=ZERO
      do i=1,nat
        xtot=xtot+epa(i,i)
        do j=1,i-1
          epa(i,j)=(epa(i,j)+epa(j,i))
          epa(j,i)=epa(i,j)
          xtot=xtot+epa(i,j)
        end do
      end do
      eelnuc=xtot

      if(iecp.eq.1) then
        call enpart_matrix('ELECTRON-NUCLEAR ATTRACTION','Includes the'//
     +    ' pseudopotential (ECP) energy of each atom on the diagonal.',
     +    ' ',epa,'Electron-nuclear attraction',eelnuc)
      else
        call enpart_matrix('ELECTRON-NUCLEAR ATTRACTION',' ',' ',epa,
     +    'Electron-nuclear attraction',eelnuc)
      end if

!! xxdip is added here (not just eelnuc) since it's excluded from the    !!
!! atomic epa contributions above -- keeps this check apples-to-apples.  !!
      lref=eelnuc0.ne.ZERO
      errk=ZERO
      if(lref) then
        errk=(eelnuc+xxdip-eelnuc0)*tokcal
        call print_value('Integration error',errk,2,'kcal/mol')
      else
        eelnuc0=eelnuc
      end if
      call enpart_frag('ELECTRON-NUCLEAR ATTRACTION',epa)
      call json_enpart_term('electron_nuclear','Electron-nuclear '//
     +  'attraction',epa,eelnuc,lref,errk)

!! second-derivative AO values, needed for the kinetic-energy density below !!
      ALLOCATE(chpaux(itotps,igr))
      call dpoints(chpaux,pcoord)

!! kinetic-energy density at each point: transform the 2nd-derivative AOs !!
!! to MOs and dot with the already-transformed chp2 -- density in scr,   !!
!! off by the usual factor of two. parallel over k: each k only reads    !!
!! its own chpaux(k,:)/chp2(k,:) and writes only scr(k).                 !!
!$OMP PARALLEL DO PRIVATE(k,x,j,xx,i)
      do k=1,ipoints
        x=ZERO
        do j=1,nocc
          xx=ZERO
          do i=1,igr
            xx=xx+c(i,j)*chpaux(k,i)
          end do
          x=x+chp2(k,j)*xx
        end do
        scr(k)=x
      end do
!$OMP END PARALLEL DO

!! contract the kinetic-energy density onto atoms; parallel over icenter, !!
!! each iteration writes only its own ekin(icenter,icenter), xtot is a    !!
!! genuine running total.                                                 !!
      xtot=ZERO
!$OMP PARALLEL DO PRIVATE(icenter,x,ifut) REDUCTION(+:xtot)
      do icenter=1,nat
        x=ZERO
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          x=x+wp(ifut)*scr(ifut)*omp2(ifut,icenter)
        end do
        ekin(icenter,icenter)=ekin(icenter,icenter)+x
        xtot=xtot+x
      end do
!$OMP END PARALLEL DO
      ekinen=xtot

      do i=1,nat
        tkin(i)=ekin(i,i)
      end do
      call enpart_atom_table('KINETIC ENERGY','One-center terms only:'//
     +  ' the kinetic energy has no pair terms.',' ',1,tkin,'Kinetic',
     +  ' ',.false.)
      call print_value('Kinetic energy',ekinen,ndec,'au')
      lref=ekin0.ne.ZERO
      errk=ZERO
      if(lref) then
        errk=(ekinen-ekin0)*tokcal
        call print_value('Integration error',errk,2,'kcal/mol')
      else
        ekin0=ekinen
      end if
      if(idofr.eq.1) call group_by_frag_vec(1,'FRAGMENT KINETIC ENERGY',
     +  tkin,'Kinetic',' ')
!! the results file keeps the kinetic energy as a (diagonal) matrix, !!
!! like the other terms                                              !!
      call json_enpart_term('kinetic','Kinetic energy',ekin,ekinen,lref,
     +  errk)
      DEALLOCATE(chpaux)

      erep=ZERO
      do i=1,nat
        eto(i,i)=ekin(i,i)+epa(i,i)
        do j=i+1,nat
          dist=dsqrt((coord(1,i)-coord(1,j))**TWO+(coord(2,i)-coord(2,j))**TWO+(coord(3,i)-coord(3,j))**TWO)
          eto(i,j)=ekin(i,j)+epa(i,j)+zn(i)*zn(j)/dist
          eto(j,i)=eto(i,j)

!! enuc: nuclear-repulsion pair matrix, kept only for the printout below. !!
          enuc(i,j)=zn(i)*zn(j)/dist
          enuc(j,i)=enuc(i,j)
          erep=erep+zn(i)*zn(j)/dist
        end do
      end do

      call enpart_matrix('NUCLEAR-NUCLEAR REPULSION',' ',' ',enuc,
     +  'Nuclear-nuclear repulsion',erep)
      call enpart_frag('NUCLEAR-NUCLEAR REPULSION',enuc)
      call json_enpart_term('nuclear_repulsion','Nuclear-nuclear '//
     +  'repulsion',enuc,erep,.false.,0.0d0)

      if(ifield.eq.1) then
        xxx=ZERO
        xtot=ZERO
        do i=1,nat
          xxx=zn(i)*(coord(1,i)*field(2)+coord(2,i)*field(3)+coord(3,i)*field(4))
          ect(i,i)=ect(i,i)+xxx
          xtot=xtot+xxx
        end do
        edipole=xtot+xxdip
        call print_box('EXTERNAL ELECTRIC FIELD')
        write(cfld,'(3f10.6)') field(2),field(3),field(4)
        call print_value('Field x, y, z',trim(adjustl(cfld)),'au')
        call print_value('Electronic dipole term',xxdip,ndec,'au')
        call print_value('Nuclear dipole term',xtot,ndec,'au')
        call print_value('Field interaction energy',edipole,ndec,'au')
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Not part of the IQA terms: added to their'//
     +    ' total in the energy check'
        write(*,'(2x,a)') 'of IQA ENERGY TERMS.'
        do i=1,nat
          fdip(i,1)=edip(i,i)
          fdip(i,2)=ect(i,i)
        end do
        call enpart_atom_table('FIELD INTERACTION PER ATOM','Intrinsic'//
     +    ' dipole term (origin independent) and charge-transfer',
     +    'term (origin dependent) of each atom.',2,fdip,'Dipole',
     +    'Charge tr.',.true.)
        if(idofr.eq.1) call group_by_frag_vec(2,
     +    'FRAGMENT FIELD INTERACTION',fdip,'Dipole','Charge tr.')
      end if

      etot=ekinen+eelnuc+erep
      etot0=ekin0+eelnuc0+erep

      DEALLOCATE(chp2,scr)

      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: enpart_matrix
!! Prints one ENPART atom x atom matrix under its own title: optional
!! caption lines, the matrix and its total (diagonal and one triangle).
!! arguments:
!!   title     (in)  section title
!!   cap1,cap2 (in)  caption lines, ' ' for none
!!   A         (in)  (maxat,maxat) matrix
!!   label     (in)  label of the total line (the unit, au, is added)
!!   xtot      (in)  the total
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine enpart_matrix(title,cap1,cap2,A,label,xtot)
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /printout/iaccur
      character*(*) title,cap1,cap2,label
      dimension A(maxat,maxat)

      call print_box(title)
      if(len_trim(cap1).gt.0) then
        write(*,'(2x,a)') trim(cap1)
        if(len_trim(cap2).gt.0) write(*,'(2x,a)') trim(cap2)
        write(*,'(a)') ''
      end if
      call MPRINT2(A,nat,maxat)
      ndec=7
      if(iaccur.ne.0) ndec=13
      call print_value(label,xtot,ndec,'au')
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: enpart_frag
!! With DOFRAGS, prints the fragment matrix of an ENPART matrix, titled
!! "FRAGMENT " + the title of the atom matrix.
!! arguments:
!!   title (in)  title of the atom matrix
!!   A     (in)  (maxat,maxat) atom matrix
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine enpart_frag(title,A)
      use input_options_mod, only: idofr
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      character*(*) title
      dimension A(maxat,maxat)

      if(idofr.eq.1) call group_by_frag_mat(1,'FRAGMENT '//trim(title),A)
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: enpart_atom_table
!! Prints a per-atom ENPART table (one or two columns) under its own
!! title, with optional caption lines and Total row.
!! arguments:
!!   title     (in)  section title
!!   cap1,cap2 (in)  caption lines, ' ' for none
!!   ncol      (in)  number of columns (1 or 2)
!!   arr       (in)  (maxat,ncol) values
!!   h1,h2     (in)  column headings (h2 unused when ncol=1)
!!   ltot      (in)  logical, print a Total row (column sums)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine enpart_atom_table(title,cap1,cap2,ncol,arr,h1,h2,ltot)
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /printout/iaccur
      character*(*) title,cap1,cap2,h1,h2
      logical ltot
      dimension arr(maxat,ncol),tot(2)
      character*30 fmthdr,fmttot

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
      if(len_trim(cap1).gt.0) then
        write(*,'(2x,a)') trim(cap1)
        if(len_trim(cap2).gt.0) write(*,'(2x,a)') trim(cap2)
        write(*,'(a)') ''
      end if
      if(ncol.eq.2) then
        write(*,fmthdr) 'Atom   ',trim(h1),trim(h2)
      else
        write(*,fmthdr) 'Atom   ',trim(h1)
      end if
      write(*,'(2x,a)') repeat('-',7+iw*ncol)
      call vprint(arr,nat,maxat,ncol)
      write(*,'(2x,a)') repeat('-',7+iw*ncol)
      if(ltot) then
        do j=1,ncol
          tot(j)=0.0d0
          do i=1,nat
            tot(j)=tot(j)+arr(i,j)
          end do
        end do
        write(*,fmttot) 'Total  ',(tot(j),j=1,ncol)
      end if
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: enpart_ee_head
!! Opens the TWO-ELECTRON INTEGRATION ERROR section: the electron-electron
!! energy on the first grid, its error against the reference energy of
!! the .fchk and whether the zero-error correction follows, or a note
!! when the .fchk has no reference energy.
!! arguments:
!!   evee     (in)  electron-electron energy (au)
!!   twoelerr (in)  its integration error (kcal/mol)
!!   iref     (in)  1 if the .fchk has the reference energy
!!   tol      (in)  tolerance of the error (kcal/mol)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine enpart_ee_head(evee,twoelerr,iref,tol)
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      common /printout/iaccur

      ndec=7
      if(iaccur.ne.0) ndec=13
      call print_box('TWO-ELECTRON INTEGRATION ERROR')
      call print_value('Electron-electron energy',evee,ndec,'au')
      if(iref.eq.1) then
        call print_value('Integration error',twoelerr,2,'kcal/mol')
        call print_value('Tolerance',tol,2,'kcal/mol')
        if(abs(twoelerr).gt.tol) then
          call print_value('Zero-error correction',
     +      'applied to the one-center terms')
        else
          call print_value('Zero-error correction',
     +      'not needed (error within the tolerance)')
        end if
      else
        write(*,'(a)') ''
        write(*,'(2x,a)') 'No reference electron-electron energy in the'//
     +    ' .fchk file: the error'
        write(*,'(2x,a)') 'is neither estimated nor corrected. To enable'//
     +    ' it, append the'
        write(*,'(2x,a)') 'reference energies to the .fchk'//
     +    ' (utils/get_energy_g16 for Gaussian 16;'
        write(*,'(2x,a)') 'utils/apost3d.py writes them for pySCF).'
      end if
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: enpart_ee_tail
!! Closes the TWO-ELECTRON INTEGRATION ERROR section after the zero-error
!! correction: error on the rotated grid, damping parameter, corrected
!! electron-electron energy and its error.
!! arguments:
!!   err1    (in)  error on the first grid (kcal/mol)
!!   err2    (in)  error on the rotated grid (kcal/mol)
!!   damp    (in)  damping parameter of the interpolation
!!   evee    (in)  corrected electron-electron energy (au)
!!   errc    (in)  its error (kcal/mol)
!!   lphf    (in)  logical, correlated wavefunction (the corrected terms
!!                 are its exchange-correlation, not an exchange)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine enpart_ee_tail(err1,err2,damp,evee,errc,lphf)
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      common /printout/iaccur
      logical lphf

      ndec=7
      if(iaccur.ne.0) ndec=13
      call print_value('Error on the rotated grid',err2,2,'kcal/mol')
      if(err1*err2.gt.0.0d0) write(*,'(2x,a)')
     +  'Same-sign error on both grids: extrapolating, not interpolating'
!! 7 decimals also with FULLPRECISION: a ratio of small differences, !!
!! its further digits depend on the summation order                  !!
      call print_value('Damping parameter',damp,7)
      call print_value('Corrected electron-electron energy',evee,ndec,
     +  'au')
      call print_value('Corrected integration error',errc,2,'kcal/mol')
      write(*,'(a)') ''
      if(lphf) then
        write(*,'(2x,a)') 'The one-center Coulomb and exchange-correlation'//
     +    ' terms of the matrices'
        write(*,'(2x,a)') 'below are the corrected ones; the pair terms'//
     +    ' do not change.'
      else
        write(*,'(2x,a)') 'The one-center Coulomb and exchange terms of the'//
     +    ' matrices below'
        write(*,'(2x,a)') 'are the corrected ones; the pair terms do not'//
     +    ' change.'
      end if
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: enpart_two_mats
!! Prints the two-electron matrices of HF and KS-DFT: Coulomb, and the
!! exchange (HF) or, for a hybrid functional, the exact-exchange share
!! (scaled by the exact-exchange fraction, as it enters the energy) and
!! the complete exchange-correlation; with VERBOSE 1 also the unscaled
!! HF-type exchange. A pure functional has only the Coulomb matrix here
!! (its exchange-correlation is the DFT part, printed before).
!! arguments:
!!   ihf         (in)  1 for Hartree-Fock
!!   xmix        (in)  exact-exchange fraction (0 for a pure functional)
!!   coul,xcoul  (in)  Coulomb matrix and total
!!   exhf,xexhf  (in)  unscaled HF-type exchange matrix and total
!!   exc,xexc    (in)  complete exchange-correlation matrix and total
!!                     (hybrid functionals)
!!   cstage      (in)  appended to the titles (' ' for the final values)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine enpart_two_mats(ihf,xmix,coul,xcoul,exhf,xexhf,exc,
     +  xexc,cstage)
      use input_options_mod, only: iverbose
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension coul(maxat,maxat),exhf(maxat,maxat),exc(maxat,maxat)
      allocatable :: sc(:,:)
      character*(*) cstage
      character*80 c1,c2,ct
      logical lfin

!! the results file gets the final matrices only !!
      lfin=len_trim(cstage).eq.0

      ct='COULOMB ENERGY'//cstage
      call enpart_matrix(ct,' ',' ',coul,'Coulomb energy',xcoul)
      call enpart_frag(ct,coul)
      if(lfin) call json_enpart_term('coulomb','Coulomb energy',coul,
     +  xcoul,.false.,0.0d0)

      if(ihf.eq.1) then
        ct='EXCHANGE ENERGY (HARTREE-FOCK)'//cstage
        call enpart_matrix(ct,' ',' ',exhf,'Exchange energy',xexhf)
        call enpart_frag(ct,exhf)
        if(lfin) call json_enpart_term('xc','Exchange energy (HF)',exhf,
     +    xexhf,.false.,0.0d0)
      else if(xmix.gt.0.0d0) then
        if(iverbose.eq.1) then
          ct='HF-TYPE EXCHANGE, UNSCALED'//cstage
          write(c1,'(a,f5.3,a)') 'The full HF-type exchange: only ',xmix,
     +      ' x this matrix enters'
          c2='the energy (see EXACT-EXCHANGE SHARE).'
          call enpart_matrix(ct,c1,c2,exhf,'HF-type exchange, unscaled',
     +      xexhf)
          call enpart_frag(ct,exhf)
        end if
        ALLOCATE(sc(maxat,maxat))
        sc=0.0d0
        xsc=0.0d0
        do i=1,nat
          do j=1,nat
            sc(i,j)=xmix*exhf(i,j)
          end do
          do j=i,nat
            xsc=xsc+sc(i,j)
          end do
        end do
        ct='EXACT-EXCHANGE SHARE'//cstage
        write(c1,'(f5.3,a)') xmix,' x the HF-type exchange: the exact-'//
     +    'exchange part of the'
        c2='functional, already scaled, as it enters the energy.'
        call enpart_matrix(ct,c1,c2,sc,'Exact-exchange share',xsc)
        call enpart_frag(ct,sc)
        if(lfin) call json_enpart_term('exact_exchange_share',
     +    'Exact-exchange share',sc,xsc,.false.,0.0d0)
        DEALLOCATE(sc)
        ct='EXCHANGE-CORRELATION ENERGY'//cstage
        c1='DFT part (DFT EXCHANGE-CORRELATION ENERGY) + exact-exchange'//
     +    ' share.'
        call enpart_matrix(ct,c1,' ',exc,'Exchange-correlation energy',
     +    xexc)
        call enpart_frag(ct,exc)
        if(lfin) call json_enpart_term('xc','Exchange-correlation '//
     +    'energy',exc,xexc,.false.,0.0d0)
      else
!! pure functional: its exchange-correlation is the DFT part !!
        if(lfin) call json_enpart_term('xc','Exchange-correlation '//
     +    'energy',exc,xexc,.false.,0.0d0)
      end if
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: json_enpart_term
!! Writes one ENPART term to the results file (inside the enpart block):
!! an object with the atom matrix and, with DOFRAGS, the fragment matrix
!! (self terms on the diagonal, whole pair terms off it, totals over the
!! diagonal and one triangle), and the integration error of the term
!! when the output prints one.
!! arguments:
!!   key   (in)  key of the term ('total' and 'classical_interaction'
!!               get their own nouns for the diagonal and off-diagonal)
!!   title (in)  title of the term for the index
!!   A     (in)  (maxat,maxat) atom matrix
!!   xtot  (in)  its total
!!   lerr  (in)  logical, write errk
!!   errk  (in)  integration error of the term (kcal/mol)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_enpart_term(key,title,A,xtot,lerr,errk)
      use input_options_mod, only: idofr
      use json_mod
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /printout/iaccur
      character*(*) key,title
      dimension A(maxat,maxat)
      logical lerr
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort
      character*30 cdiag,coff

      call atomic_definition_names(ckey,cname,cshort,cspace)
      ndec=6
      if(iaccur.ne.0) ndec=13
!! nouns for legends: the IQA matrix has self-energies and interaction !!
!! energies; the classical interactions have no self terms              !!
      cdiag='self term'
      coff='pair term'
      if(key.eq.'total') then
        cdiag='self-energy'
        coff='interaction energy'
      else if(key.eq.'classical_interaction') then
        cdiag='none'
        coff='classical interaction energy'
      end if
      call json_begin_object(key)
      call json_atom_matrix('atoms',A,.true.,trim(cdiag),trim(coff),1,
     +  xtot)
      call json_index('atoms','atom_matrix',trim(title)//' ('//
     +  trim(cshort)//')','au',ndec,' ',.true.,' ')
      if(idofr.eq.1) then
        call json_frag_matrix('fragments',1,A,.true.,trim(cdiag),
     +    trim(coff))
        call json_index('fragments','fragment_matrix','Fragment '//
     +    trim(title)//' ('//trim(cshort)//')','au',ndec,' ',.true.,' ')
      end if
      if(lerr) then
        call json_value('integration_error_kcal',errk)
        call json_index('integration_error_kcal','scalar',trim(title)//
     +    ', integration error','kcal/mol',2,' ',.true.,' ')
      end if
      call json_end()
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: json_enpart_ee
!! Writes the two_electron_integration object of the enpart block: the
!! electron-electron energy on the first grid, its error, the tolerance,
!! whether the zero-error correction was applied and its results.
!! arguments:
!!   evee1 (in)  electron-electron energy, first grid (au)
!!   err1  (in)  its integration error (kcal/mol)
!!   iref  (in)  1 if the .fchk has the reference energy
!!   tol   (in)  tolerance (kcal/mol)
!!   izerr (in)  1 if the zero-error correction was applied
!!   err2  (in)  error on the rotated grid (kcal/mol), izerr=1
!!   damp  (in)  damping parameter, izerr=1
!!   evee  (in)  corrected electron-electron energy (au), izerr=1
!!   errc  (in)  its error (kcal/mol), izerr=1
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_enpart_ee(evee1,err1,iref,tol,izerr,err2,damp,
     +  evee,errc)
      use json_mod
      implicit real*8(a-h,o-z)
      common /printout/iaccur

      ndec=7
      if(iaccur.ne.0) ndec=13
      call json_begin_object('two_electron_integration')
      call json_value('energy',evee1)
      call json_index('energy','scalar','Electron-electron energy, '//
     +  'first grid','au',ndec,' ',.false.,' ')
      if(iref.eq.1) then
        call json_value('error_kcal',err1)
        call json_value('tolerance_kcal',tol)
        if(izerr.eq.1) then
          call json_value('zero_error_correction','applied')
        else
          call json_value('zero_error_correction','not needed')
        end if
      else
        call json_null('error_kcal')
        call json_null('tolerance_kcal')
        call json_value('zero_error_correction','no reference energy')
      end if
      call json_index('error_kcal','scalar','Two-electron integration '//
     +  'error, first grid','kcal/mol',2,' ',.true.,' ')
      if(izerr.eq.1) then
        call json_value('rotated_grid_error_kcal',err2)
        call json_value('damping',damp)
        call json_value('corrected_energy',evee)
        call json_value('corrected_error_kcal',errc)
        call json_index('corrected_error_kcal','scalar','Two-electron '//
     +    'integration error, corrected','kcal/mol',2,' ',.true.,' ')
      else
        call json_null('rotated_grid_error_kcal')
        call json_null('damping')
        call json_null('corrected_energy')
        call json_null('corrected_error_kcal')
      end if
      call json_end()
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: json_enpart_pairs
!! Writes the atom-pair threshold of the enpart block and the atom pairs
!! below it, whose exchange(-correlation) term comes from the multipolar
!! expansion.
!! arguments: none (THREBOD from input_options_mod, bond orders /ovpop/)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_enpart_pairs()
      use input_options_mod, only: ithrebod
      use json_mod
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      dimension ipair(2)

      if(ithrebod.lt.1) then
        threbod=0.0d0
      else
        threbod=real(ithrebod)/10000.0d0
      end if
      call json_value('pair_threshold',threbod)
      call json_begin_array('multipolar_pairs')
      do i=1,nat
        do j=i+1,nat
          if(bo(i,j).lt.threbod) then
            ipair(1)=i
            ipair(2)=j
            call json_vector(' ',ipair,2)
          end if
        end do
      end do
      call json_end()
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: enpart_iqa
!! Prints the final IQA matrix (self-energies on the diagonal, interaction
!! energies off it), its total against the energy of the .fchk, and the
!! summary built from the final matrices: additive atomic energies and
!! classical interaction energies.
!! arguments:
!!   eto    (in)  (maxat,maxat) final IQA matrix
!!   etot   (in)  its total
!!   escf   (in)  energy of the .fchk
!!   ifield (in)  1 with an external field
!!   edip   (in)  field interaction energy (ifield=1)
!!   xcfin  (in)  (maxat,maxat) final exchange-correlation matrix
!!   ljson  (in)  logical, write the results-file entries (the caller has
!!                the enpart block open)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine enpart_iqa(eto,etot,escf,ifield,edip,xcfin,ljson)
      use output_mod, only: print_value
      use input_options_mod, only: idofr
      use json_mod
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,
     +  jfrlist(maxat)
      common /printout/iaccur
      dimension eto(maxat,maxat),xcfin(maxat,maxat)
      dimension arr(maxat,2)
      allocatable :: vcl(:,:),bm(:,:)
      dimension bf(maxat,2),bv(maxat,1),x(2),xs(1),fself(maxat)
      logical ljson
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort
      character*80 ct

      ndec=7
      if(iaccur.ne.0) ndec=13
      call atomic_definition_names(ckey,cname,cshort,cspace)
      ct='IQA ENERGY TERMS ('//trim(cshort)//')'
      call enpart_matrix(ct,'Self-energies on the diagonal, interaction'//
     +  ' energies off it.',' ',eto,'Total energy',etot)
      call print_value('Energy in the .fchk file',escf,ndec,'au')
      if(ifield.eq.1) then
        call print_value('Field interaction energy',edip,ndec,'au')
        err=etot-escf+edip
      else
        err=etot-escf
      end if
      call print_value('Integration error',err,ndec,'au')
      call print_value('Integration error',err*tokcal,2,'kcal/mol')
      call enpart_frag(ct,eto)
      if(ljson) then
        call json_enpart_term('total','IQA energy terms',eto,etot,
     +    .false.,0.0d0)
        call json_value('total_energy',etot)
        call json_index('total_energy','scalar','Total energy (IQA)',
     +    'au',ndec,' ',.true.,' ')
        call json_value('reference_energy',escf)
!! a field run: integration_error = total + field_energy - reference !!
        if(ifield.eq.1) then
          call json_value('field_energy',edip)
        else
          call json_null('field_energy')
        end if
        call json_value('integration_error',err)
        call json_index('integration_error','scalar','Integration '//
     +    'error of the total energy','au',ndec,' ',.true.,' ')
        call json_value('integration_error_kcal',err*tokcal)
      end if

!! additive atomic energies: self-energy + half of each interaction !!
      do i=1,nat
        arr(i,1)=eto(i,i)
        xx=0.0d0
        do j=1,nat
          if(j.ne.i) xx=xx+eto(i,j)
        end do
        arr(i,2)=eto(i,i)+xx/2.0d0
      end do
      call enpart_atom_table('ATOMIC ENERGIES','Additive energy = self-'//
     +  'energy + half of each interaction energy of the atom;',
     +  'the additive energies add up to the total energy.',2,arr,
     +  'Self-energy','Additive',.true.)
      if(ifield.eq.1) write(*,'(2x,a)')
     +  'The field interaction energy is not included.'
      ndv=6
      if(iaccur.ne.0) ndv=13
      if(ljson) then
        call json_begin_object('atomic_energies')
        call json_begin_object('atoms')
        call json_vector('self_energy',arr(1:nat,1),nat)
        call json_index('self_energy','atom_vector','Atomic self-'//
     +    'energies','au',ndv,' ',.true.,' ')
        call json_vector('additive_energy',arr(1:nat,2),nat)
        call json_index('additive_energy','atom_vector','Additive '//
     +    'atomic energies','au',ndv,' ',.true.,' ')
        call json_end()
      end if
      if(idofr.eq.1) then
        ALLOCATE(bm(maxat,maxat))
        call frag_sum_mat(1,eto,bm,xx)
        call frag_sum_vec(1,arr(1,2),bv,xs)
        x(1)=0.0d0
        x(2)=0.0d0
        do k=1,icufr
          bf(k,1)=bm(k,k)
          bf(k,2)=bv(k,1)
          x(1)=x(1)+bf(k,1)
          x(2)=x(2)+bf(k,2)
        end do
        call print_frag_table(2,'FRAGMENT ATOMIC ENERGIES',bf,x,
     +    'Self-energy','Additive')
        do k=1,icufr
          fself(k)=bf(k,1)
        end do
        if(ljson) then
          call json_begin_object('fragments')
          call json_vector('self_energy',fself,icufr)
          call json_index('self_energy','fragment_vector','Fragment '//
     +      'self-energies','au',ndv,' ',.true.,' ')
          call json_frag_vector('additive_energy',arr(1,2))
          call json_index('additive_energy','fragment_vector','Additive'//
     +      ' fragment energies','au',ndv,' ',.true.,' ')
          call json_end()
        end if
        DEALLOCATE(bm)
      end if

!! classical interaction energies: interaction minus its xc part !!
      ALLOCATE(vcl(maxat,maxat))
      vcl=0.0d0
      xx=0.0d0
      do i=1,nat
        vcl(i,i)=0.0d0
        do j=1,nat
          if(j.ne.i) vcl(i,j)=eto(i,j)-xcfin(i,j)
          if(j.gt.i) xx=xx+vcl(i,j)
        end do
      end do
      ct='CLASSICAL INTERACTION ENERGIES'
      call enpart_matrix(ct,'Electron-nuclear + nuclear-nuclear +'//
     +  ' Coulomb part of each interaction;','the rest of it is the'//
     +  ' exchange-correlation part. No self terms.',vcl,
     +  'Classical interaction energy',xx)
      if(ljson) call json_end()
      call enpart_frag(ct,vcl)
      if(ljson) call json_enpart_term('classical_interaction',
     +  'Classical interaction energies',vcl,xx,.false.,0.0d0)
      DEALLOCATE(vcl)
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: numint_two                                                  !!
!! purpose: two-electron IQA/ENPART energy partition, RHF/RKS -- Coulomb   !!
!!   (via calc_coul) and, if requested, HF-type exchange (same-center      !!
!!   and atom-pair kernels, skipping pairs below the THREBOD bond-order    !!
!!   threshold via a multipolar approximation instead), each decomposed    !!
!!   into atomic/diatomic contributions. A second pass on a phase-         !!
!!   rotated grid estimates and corrects the numerical-integration error.  !!
!!   (the "zero-error" interpolation) when it exceeds tolerance. See       !!
!!   numint_two_uhf for open-shell.                                        !!
!! arguments:                                                              !!
!!   ndim   (in)  -- number of basis functions (leading dim of chp)        !!
!!   itotps (in)  -- total number of grid points                           !!
!!   wp     (in)  -- integration weight of each grid point                 !!
!!   omp    (in)  -- becke/tfvc weight of each grid point for its own atom !!
!!   omp2   (in)  -- becke/tfvc (or hirshfeld) weight of each point for    !!
!!                   every atom                                            !!
!!   pcoord (in)  -- xyz coordinates of each grid point                    !!
!!   chp    (in)  -- basis-function values at each grid point              !!
!!   rho    (in)  -- electron density at each grid point                   !!
!!   eto    (inout) -- total energy matrix, accumulated on top of the      !!
!!                   one-electron part already in it                       !!
!! author: PSalse, MGimf; loop order for OpenMP: R. Oswald.               !!
!! *********************************************************************** !!
      subroutine numint_two(ndim,itotps,wp,omp,omp2,pcoord,chp,rho,eto)

      use ao_matrices
      use integration_grid
      use OMP_LIB
      use timing_mod
      use output_mod, only: print_value
      use input_options_mod, only: iverbose

      implicit real*8(a-h,o-z)


      include 'parameter.h'
      parameter (maxpass=5)

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /iops/iopt(200)
      common /actual/iact,jat,icenter
      common /energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common /energ0/ekin0,eelnuc0,evee0,etot0
      common /exchg/exch(maxat,maxat),xmix
      common /twoel/twoeltoler
      common /efield/field(4),edipole
!! rotated-grid angles/thresholds, read from the # GRID input section. !!
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22

      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3)
      dimension omp(itotps),omp2(itotps,nat),rho(itotps)
      dimension eto(maxat,maxat),Excmp(maxat,maxat)

      dimension coul0(maxat,2*maxpass)
      allocatable :: coulb(:,:),exhfb(:,:),exchb(:,:)
      dimension exch_hf(maxat,maxat)
      dimension coul(maxat,maxat)
      
      character*100 threadenv

      allocatable :: chp2(:,:)
      allocatable :: chppha(:,:), wppha(:),pcoordpha(:,:),omppha(:)
      allocatable :: chp2pha(:,:),rhopha(:),omp2pha(:,:),ibaspointpha(:)

!! parallelization-related arrays below. !!
      allocatable :: chp2s(:,:),chp2phas(:,:) !! transposed layout (rows/columns swapped) !!
      allocatable :: ijpaircount(:,:),istart(:),iend(:) !! included atom pairs, thread-slicing bookkeeping !!
      allocatable :: exch_hfk(:,:)
      allocatable :: exch_hfij(:,:)
      allocatable :: f3k(:)

      ihf         = Iopt(18) 
      ithrebod    = Iopt(44)
      ifield      = Iopt(47)
      ipolar      = Iopt(46)
      ianalytical = Iopt(91)

      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if
      idoex=0
      if(xmix.gt.ZERO) idoex=1
      iatps=nang*nrad
!! zero-error strategy only possible with a reference electron-electron !!
!! energy in the .fchk (evee0, see mga_misc); twoelerr stays defined.    !!
      iref=0
      if(evee0.ne.ZERO) iref=1
      twoelerr=ZERO

      npass=2
!! HF: the atom pairs integrated and those from the multipolar expansion !!
!! (for KS-DFT the DFT part has already printed them)                   !!
      if(ihf.eq.1.and.ianalytical.eq.0) then
        npairs=0
        nskip=0
        do i=1,nat
          do j=i+1,nat
            if(bo(i,j).ge.threbod) then
              npairs=npairs+1
            else
              nskip=nskip+1
            end if
          end do
        end do
        call print_value('Pairs integrated (exact exchange)',npairs)
        call print_value('Pairs from the multipolar expansion',nskip)
      end if
      do i=1,nat
        do j=1,nat
          coul(i,j)=ZERO
          exch_hf(i,j)=ZERO
        end do
        do j=1,2*npass
          coul0(i,j)=ZERO
        end do
      end do

      ALLOCATE(chp2(itotps,nocc))
      ALLOCATE(chp2s(nocc,itotps))

      call ao_to_mo_grid_t(itotps,igr,nocc,c,chp,chp2,chp2s)

      if(ianalytical.eq.0) then !! skip this whole rotated-grid pass for the analytical two-electron path !!

!! rotated grid, angle controlled by the # GRID section -- see the diatXC !!
        phb=phb12
        pha=ZERO
        ALLOCATE(wppha(itotps),omppha(itotps),omp2pha(itotps,nat))
        ALLOCATE(chppha(itotps,ndim),pcoordpha(itotps,3),ibaspointpha(itotps))
        ALLOCATE(chp2pha(itotps,nocc),rhopha(itotps))
        ALLOCATE(chp2phas(nocc,itotps))

        call prenumint(ndim,itotps,nat,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)

        call ao_to_mo_grid_t(itotps,igr,nocc,c,chppha,chp2pha,chp2phas)
      end if !! non-analytical skip ends here !!

!! Coulomb energy term. !!
      call cpu_time(xtime)
      call get_wall_time(wxtime)
      if(ianalytical.eq.0) then
        call calc_coul(itotps,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,rho,rhopha,coul)
      else
!! chp2b is not declared here: under implicit typing it is a scalar, not !!
!! the rank-2 array calc_twoel_analytical expects (gfortran warns of the  !!
!! rank mismatch). Only reached with ANALYTIC, which no test covers.      !!
        call calc_twoel_analytical(itotps,wp,omp2,pcoord,rho,chp2,chp2b,coul,exch_hf)
      end if
      call cpu_time(xtime2)
      call get_wall_time(wxtime2)
      write(*,'(a)') ''
      call print_timer('ENPART Coulomb',xtime2-xtime,wxtime2-wxtime)
      xtime=xtime2
      wxtime=wxtime2

      if(ianalytical.eq.0) then !! skip HF-type exchange for the analytical path !!
        if(idoex.eq.1) then

!! multipolar approximation, used below for atom pairs skipped by THREBOD !!
!! (closed shell: weight -2 for phi_i^2, -4 for phi_i*phi_j, i<j)        !!
          call multipolar_mo(nocc,itotps,wp,omp2,pcoord,chp2,-TWO,-FOUR,Excmp)

!! loops ordered for the thread splitting below (R. Oswald) !!

!! same-center terms first; determine the core count for thread splitting. !!
          call getenv('OMP_NUM_THREADS',threadenv)
          if(trim(threadenv)=='') then
            ithreadenv=ZERO
          else
            read(unit=threadenv,FMT='(I4)') ithreadenv
          end if
!! a parallel region only to query the thread and processor counts;  !!
!! every thread writes the same values to iprocs/ithreads            !!
        !$OMP PARALLEL
          iprocs=OMP_GET_MAX_THREADS()
          ithreads=INT(OMP_GET_NUM_PROCS())
        !$OMP END PARALLEL
          if(ithreadenv.ne.ZERO) ithreads=ithreadenv

!! per-thread grid-point slicing bookkeeping. !!
          ALLOCATE(f3k(ithreads))
          ALLOCATE(exch_hfk(nat,ithreads))
          ALLOCATE(istart(ithreads),iend(ithreads))
          itilerest=mod(iatps,ithreads)
          ispace=(iatps-itilerest)/ithreads

!! loop kept in this manual-chunking shape rather than a plain OMP loop --  !!
!! avoids an old cap that limited parallelization to 8 cores.               !!
          exch_hfk=ZERO
          do icenter=1,nat
            ioffset=(icenter-1)*iatps
            istart=0
            iend=0
            do ik=1,(ithreads-1)
              istart(ik)=ioffset+((ik-1)*ispace)+1
              iend(ik)=ioffset+(ik*ispace)
            end do
            istart(ithreads)=ioffset+((ithreads-1)*ispace)+1
            iend(ithreads)=icenter*iatps

!! both calls below are required for the thread count to actually take effect. !!
            call omp_set_dynamic(.false.)
            call omp_set_num_threads(ithreads)
!! parallel over the thread slices ik: f3loc is a PRIVATE scalar that     !!
!! absorbs the i/j accumulation; exch_hfk(icenter,ik) is SHARED but        !!
!! index-partitioned by ik and written once per ifut, which avoids the     !!
!! false sharing of writing adjacent ik slots in every inner iteration.    !!
!$OMP PARALLEL DO PRIVATE(ifut,jfut,x0,dx0,dy0,dz0,x1,dx1,dy1,dz1,dist,i,j,f2,f3loc)
            do ik=1,ithreads
              do ifut=istart(ik),iend(ik)
                x0=wp(ifut)*omp2(ifut,icenter)
                dx0=pcoord(ifut,1)
                dy0=pcoord(ifut,2)
                dz0=pcoord(ifut,3)
                f3loc=ZERO
                do jfut=iatps*(icenter-1)+1,iatps*icenter
                  x1=wppha(jfut)*omp2pha(jfut,icenter)
                  dx1=pcoordpha(jfut,1)
                  dy1=pcoordpha(jfut,2)
                  dz1=pcoordpha(jfut,3)
                  dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                  if(dist.gt.1.0d-12) then
                    do i=1,nocc-1
                      f2=x0*chp2s(i,ifut)*chp2s(i,ifut)
                      f3loc=f3loc+f2*x1*chp2phas(i,jfut)*chp2phas(i,jfut)/dist
                      do j=i+1,nocc
                        f2=TWO*x0*chp2s(i,ifut)*chp2s(j,ifut)
                        f3loc=f3loc+f2*x1*chp2phas(i,jfut)*chp2phas(j,jfut)/dist
                      end do
                    end do
                    f2=x0*chp2s(nocc,ifut)*chp2s(nocc,ifut)
                    f3loc=f3loc+f2*x1*chp2phas(nocc,jfut)*chp2phas(nocc,jfut)/dist
                  end if
                end do
                exch_hfk(icenter,ik)=exch_hfk(icenter,ik)+f3loc
              end do
            end do
!$OMP END PARALLEL DO
          end do

!! reduce the per-thread accumulator into the same-center exch_hf term. !!
          do icenter=1,nat
            do isum=1,ithreads
              exch_hf(icenter,icenter)=exch_hf(icenter,icenter)-exch_hfk(icenter,isum)
            end do
          end do

!! atom-pair terms, skipping pairs below the THREBOD bond-order threshold. !!
          iterms=0
          ipaircounter=0
          ALLOCATE(ijpaircount(nat*nat,2))
          do icenter=1,nat
            do jcenter=icenter+1,nat
              bx0=bo(icenter,jcenter)
              if(bx0.ge.threbod) then
                ipaircounter=ipaircounter+1
                ijpaircount(ipaircounter,1)=icenter
                ijpaircount(ipaircounter,2)=jcenter
              else
                iterms=iterms+1
              end if
            end do
          end do

!! per-thread grid-point slicing bookkeeping for the atom-pair loop. !!
          ALLOCATE(exch_hfij(ipaircounter,ithreads))
          itilerest=mod(iatps,ithreads)
          ispace=(iatps-itilerest)/ithreads
          exch_hfij=ZERO
          do numpairnat=1,ipaircounter
            icenter=ijpaircount(numpairnat,1)
            jcenter=ijpaircount(numpairnat,2)
            ioffset=(icenter-1)*iatps
            istart=0
            iend=0
            do ik=1,(ithreads-1)
              istart(ik)=ioffset+((ik-1)*ispace)+1
              iend(ik)=ioffset+(ik*ispace)
            end do
            istart(ithreads)=ioffset+((ithreads-1)*ispace)+1
            iend(ithreads)=icenter*iatps
!! parallel over the thread slices ik, as the same-center block above: !!
!! f3loc is PRIVATE, exch_hfij(numpairnat,ik) written once per ifut.     !!
!$OMP PARALLEL DO PRIVATE(ifut,jfut,x0,dx0,dy0,dz0,x1,dx1,dy1,dz1,dist,
!$OMP&  chp2snifut,chp2snjfut,chp2sijfut,chp2siifut,i,j,f2,f3loc)
            do ik=1,ithreads
              do ifut=istart(ik),iend(ik)
                x0=wp(ifut)*omp2(ifut,icenter)
                dx0=pcoord(ifut,1)
                dy0=pcoord(ifut,2)
                dz0=pcoord(ifut,3)
                chp2snifut=chp2s(nocc,ifut)
                f3loc=ZERO
                do jfut=iatps*(jcenter-1)+1,iatps*jcenter
                  x1=wp(jfut)*omp2(jfut,jcenter)
                  dx1=pcoord(jfut,1)
                  dy1=pcoord(jfut,2)
                  dz1=pcoord(jfut,3)
                  dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                  if(dist.gt.1.0d-12) then
                    chp2snjfut=chp2s(nocc,jfut)
                    do i=1,nocc-1
                      chp2sijfut=chp2s(i,jfut)
                      chp2siifut=chp2s(i,ifut)
                      do j=i+1,nocc
                        f2=TWO*x0*chp2siifut*chp2s(j,ifut)
                        f3loc=f3loc+f2*x1*chp2sijfut*chp2s(j,jfut)/dist
                      end do
                      f2=x0*chp2siifut*chp2siifut
                      f3loc=f3loc+f2*x1*chp2sijfut*chp2sijfut/dist
                    end do
                    f2=x0*chp2snifut*chp2snifut
                    f3loc=f3loc+f2*x1*chp2snjfut*chp2snjfut/dist
                  end if
                end do
                exch_hfij(numpairnat,ik)=exch_hfij(numpairnat,ik)+f3loc
              end do
            end do
!$OMP END PARALLEL DO
          end do

!! reduce the per-thread accumulator into the atom-pair exch_hf term. !!
          do numpairnat=1,ipaircounter
            icenter=ijpaircount(numpairnat,1)
            jcenter=ijpaircount(numpairnat,2)
            do isum=1,ithreads
              exch_hf(icenter,jcenter)=exch_hf(icenter,jcenter)-exch_hfij(numpairnat,isum)
            end do
          end do

!! account for A-B/B-A symmetry (factor of 2) and for atom pairs whose !!
!! term came from the multipolar approximation instead of integration. !!
          exchen_hf=ZERO
          do i=1,nat
            exchen_hf=exchen_hf+exch_hf(i,i)
            do j=i+1,nat
              if(bo(i,j).lt.threbod) then 
                exch_hf(i,j)=Excmp(i,j)
              else 
                exch_hf(i,j)=TWO*(exch_hf(i,j))
              end if 
              exch_hf(j,i)=exch_hf(i,j)
              exchen_hf=exchen_hf+exch_hf(i,j)
            end do
          end do

          call cpu_time(xtime2)
          call get_wall_time(wxtime2)
          call print_timer('ENPART exact exchange',xtime2-xtime,
     +      wxtime2-wxtime)
          xtime=xtime2
          wxtime=wxtime2
        end if
      end if !! non-analytical skip ends here !!

!! electron-electron energy on the first grid and its integration error; !!
!! for a hybrid functional the exact-exchange share joins the DFT part   !!
      if(ihf.eq.1) then
        evee=coulen+exchen_hf
      else
        exchen=ZERO
        do i=1,nat
          exch(i,i)=exch(i,i)+xmix*exch_hf(i,i)
          exchen=exchen+exch(i,i)
          do j=i+1,nat
            exch(i,j)=exch(i,j)+xmix*exch_hf(i,j)
            exch(j,i)=exch(i,j)
            exchen=exchen+exch(i,j)
          end do
        end do
        evee=coulen+exchen
      end if
      if(evee0.ne.ZERO) then
        twoelerr=(evee-evee0)*tokcal
      else
        evee0=evee
      end if

      if(ianalytical.eq.1) twoeltoler=232000
      call enpart_ee_head(evee,twoelerr,iref,twoeltoler)
!! first-grid values, kept for the results file !!
      eveefg=evee
      errfg=twoelerr
      err2=ZERO
      phabest=ZERO

!! the zero-error correction changes only the one-center terms; the      !!
!! first-grid matrices are kept for VERBOSE 1                            !!
      izerr=0
      if(iref.eq.1.and.abs(twoelerr).gt.twoeltoler) izerr=1
      if(izerr.eq.1.and.iverbose.eq.1) then
        ALLOCATE(coulb(maxat,maxat),exhfb(maxat,maxat))
        ALLOCATE(exchb(maxat,maxat))
        do i=1,nat
          do j=1,nat
            coulb(i,j)=coul(i,j)
            exhfb(i,j)=exch_hf(i,j)
            exchb(i,j)=exch(i,j)
          end do
        end do
        xcoulb=coulen
        xexhfb=exchen_hf
        xexchb=exchen
      end if

      if(izerr.eq.1) then

!! rotated grid, angle controlled by the # GRID section -- see the diatXC !!
!! paper for optimized values: phb=0.162d0, later 0.182d0, for 40/146.   !!
        phb=phb22
        pha=ZERO
        nat0=nat
        call prenumint(ndim,itotps,nat0,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)

        call ao_to_mo_grid_t(itotps,igr,nocc,c,chppha,chp2pha,chp2phas)

!! same-center-only recompute of the Coulomb part on the rotated grid, !!
!! for the zero-error interpolation below. parallel over icenter: each !!
!! iteration accumulates into its own coul0(icenter,3), independent of !!
!! every other icenter.                                                !!
!$OMP PARALLEL DO PRIVATE(icenter,ifut,x0,dx0,dy0,dz0,f3,jfut,x1,dx1,dy1,dz1,dist)
        do icenter=1,nat
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            x0=wp(ifut)*omp2(ifut,icenter)
            dx0=pcoord(ifut,1)
            dy0=pcoord(ifut,2)
            dz0=pcoord(ifut,3)

            f3=ZERO
            do jfut=iatps*(icenter-1)+1,iatps*icenter
              x1=wppha(jfut)*omp2pha(jfut,icenter)
              dx1=pcoordpha(jfut,1)
              dy1=pcoordpha(jfut,2)
              dz1=pcoordpha(jfut,3)
              dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
              if(dist.gt.1.0d-12) f3=f3+rho(ifut)*rhopha(jfut)*x1*x0/dist
            end do
            coul0(icenter,3)=coul0(icenter,3)+f3/TWO
          end do
        end do
!$OMP END PARALLEL DO

        if(idoex.eq.1) then

!! same-center recompute, same manual-chunking shape as above. !!
          exch_hfk=ZERO
          do icenter=1,nat
            ioffset=(icenter-1)*iatps
            istart=0
            iend=0
            do ik=1,(ithreads-1)
              istart(ik)=ioffset+((ik-1)*ispace)+1
              iend(ik)=ioffset+(ik*ispace)
            end do
            istart(ithreads)=ioffset+((ithreads-1)*ispace)+1
            iend(ithreads)=icenter*iatps

!! both calls below are required for the thread count to actually take effect. !!
            call omp_set_dynamic(.false.)
            call omp_set_num_threads(ithreads)
!! parallel over the thread slices ik, as the same-center block above: !!
!! f3loc is PRIVATE, exch_hfk(icenter,ik) written once per ifut.         !!
!$OMP PARALLEL DO PRIVATE(ifut,jfut,x0,dx0,dy0,dz0,x1,dx1,dy1,dz1,dist,i,j,f2,f3loc)
            do ik=1,ithreads
              do ifut=istart(ik),iend(ik)
                x0=wp(ifut)*omp2(ifut,icenter)
                dx0=pcoord(ifut,1)
                dy0=pcoord(ifut,2)
                dz0=pcoord(ifut,3)
                f3loc=ZERO
                do jfut=iatps*(icenter-1)+1,iatps*icenter
                  x1=wppha(jfut)*omp2pha(jfut,icenter)
                  dx1=pcoordpha(jfut,1)
                  dy1=pcoordpha(jfut,2)
                  dz1=pcoordpha(jfut,3)
                  dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                if(dist.gt.1.0d-8) then !! 1e-8 here, 1e-12 in the other loops !!
                    do i=1,nocc-1
                      f2=x0*chp2s(i,ifut)*chp2s(i,ifut)
                      f3loc=f3loc+f2*x1*chp2phas(i,jfut)*chp2phas(i,jfut)/dist
                      do j=i+1,nocc
                        f2=TWO*x0*chp2s(i,ifut)*chp2s(j,ifut)
                        f3loc=f3loc+f2*x1*chp2phas(i,jfut)*chp2phas(j,jfut)/dist
                      end do
                    end do
                    f2=x0*chp2s(nocc,ifut)*chp2s(nocc,ifut)
                    f3loc=f3loc+f2*x1*chp2phas(nocc,jfut)*chp2phas(nocc,jfut)/dist
                  end if
                end do
                exch_hfk(icenter,ik)=exch_hfk(icenter,ik)+f3loc
              end do
            end do
!$OMP END PARALLEL DO
          end do

!! reduce the per-thread accumulator into the same-center coul0 term. !!
          do icenter=1,nat
            do isum=1,ithreads
              coul0(icenter,4)=coul0(icenter,4)-exch_hfk(icenter,isum)
            end do
          end do
        end if

        do i=1,nat
          coul0(i,1)=coul(i,i)
          if(idoex.eq.1) coul0(i,2)=exch_hf(i,i)
        end do
        deltaee=ZERO
        do i=1,nat
          deltaee=deltaee+coul0(i,1)-coul0(i,3)
          if(idoex.eq.1) deltaee=deltaee+xmix*(coul0(i,2)-coul0(i,4))
        end do
        deltaee=deltaee*tokcal
        phabest=ONE-(twoelerr/deltaee)
        err1=twoelerr
        err2=twoelerr-deltaee

!! interpolate energies, replacing the old Coulomb/exchange terms !!
        do i=1,nat
          coul(i,i)=coul0(i,1)*phabest+(ONE-phabest)*coul0(i,3)
          if(idoex.eq.1) exch_hf(i,i)=coul0(i,2)*phabest+(ONE-phabest)*coul0(i,4)
        end do

        coulen=ZERO
        do i=1,nat
          do j=i,nat
            coulen=coulen+coul(i,j)
          end do
        end do
        if(idoex.eq.1) then
          exchen_hf=ZERO
          do i=1,nat
            do j=i,nat
              exchen_hf=exchen_hf+exch_hf(i,j)
            end do
          end do
          if(ihf.ne.1) then
            exchen=ZERO
            do i=1,nat
              exch(i,i)=exch(i,i)+xmix*(exch_hf(i,i)-coul0(i,2))
              exchen=exchen+exch(i,i)
              do j=i+1,nat
                exchen=exchen+exch(i,j)
              end do
            end do
          end if
        else
          exchen=ZERO
          do i=1,nat
            exchen=exchen+exch(i,i)
            do j=i+1,nat
              exchen=exchen+exch(i,j)
            end do
          end do
        end if

        if(ihf.eq.0) then
          evee=coulen+exchen
        else
          evee=coulen+exchen_hf
        end if
        twoelerr=(evee-evee0)*tokcal
        call enpart_ee_tail(err1,err2,phabest,evee,twoelerr,.false.)
        write(*,'(a)') ''
        call cpu_time(xtime2)
        call get_wall_time(wxtime2)
        call print_timer('ENPART zero-error correction',xtime2-xtime,
     +    wxtime2-wxtime)
        xtime=xtime2
        wxtime=wxtime2
      end if

      call json_enpart_ee(eveefg,errfg,iref,twoeltoler,izerr,err2,
     +  phabest,evee,twoelerr)

!! final two-electron matrices; with VERBOSE 1 first the first-grid ones !!
!! when the zero-error correction was applied                            !!
      if(izerr.eq.1.and.iverbose.eq.1) call enpart_two_mats(ihf,xmix,
     +  coulb,xcoulb,exhfb,xexhfb,exchb,xexchb,
     +  ', BEFORE THE ZERO-ERROR CORRECTION')
      if(allocated(coulb)) DEALLOCATE(coulb,exhfb,exchb)
      call enpart_two_mats(ihf,xmix,coul,coulen,exch_hf,exchen_hf,exch,
     +  exchen,' ')

!! IQA matrix: one-electron terms (in eto) + Coulomb + exchange(-corr.) !!
      if(ihf.eq.1) then
        xtot=ZERO
        do i=1,nat
          eto(i,i)=eto(i,i)+coul(i,i)+exch_hf(i,i)
          xtot=xtot+eto(i,i)
          do j=i+1,nat
            eto(i,j)=eto(i,j)+coul(i,j)+exch_hf(i,j)
            eto(j,i)=eto(i,j)
            xtot=xtot+eto(i,j)
          end do
        end do
        etot=xtot
        call enpart_iqa(eto,etot,escf,ifield,edipole,exch_hf,.true.)
      else
        xtot=ZERO
        do i=1,nat
          eto(i,i)=eto(i,i)+coul(i,i)+xmix*exch_hf(i,i)
          xtot=xtot+eto(i,i)
          do j=i+1,nat
            eto(i,j)=eto(i,j)+coul(i,j)+xmix*exch_hf(i,j)
            eto(j,i)=eto(i,j)
            xtot=xtot+eto(i,j)
          end do
        end do
        etot=xtot
        call enpart_iqa(eto,etot,escf,ifield,edipole,exch,.true.)
      end if

      if(ianalytical.eq.0) then
        DEALLOCATE(chp2,chp2pha,rhopha,chppha,wppha,pcoordpha,omppha,omp2pha)
        DEALLOCATE(chp2s,chp2phas)
        if(idoex.eq.1) then
          DEALLOCATE(f3k,exch_hfk,istart,iend,ijpaircount,exch_hfij)
        end if
      else
        DEALLOCATE(chp2)
      end if

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: numint_one_uhf                                            !!
!! purpose: one-electron IQA/ENPART energy partition, UHF/UKS twin of    !!
!!   numint_one -- electron-nuclear attraction, kinetic energy, nuclear  !!
!!   repulsion and (if present) the external-field dipole terms, each    !!
!!   decomposed into atomic/diatomic contributions.                      !!
!! arguments: same as numint_one.                                        !!
!! author: PSalse, MGimf.                                                !!
!! ********************************************************************* !!
      subroutine numint_one_uhf(ndim,itotps,wp,rho,omp,omp2,pcoord,chp,eto)

      use ao_matrices
      use integration_grid
      use output_mod, only: print_value

      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(200)
      common/actual/iact,jat,icenter
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/energ0/ekin0,eelnuc0,evee0,etot0
      common/exchg/exch(maxat,maxat),xmix
      common/efield/field(4),edipole
      dimension eto(maxat,maxat)
      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3)
      dimension omp(itotps),omp2(itotps,nat),rho(itotps)

      dimension Epa(maxat,maxat),ekin(maxat,maxat),dip(maxat,3),enuc(maxat,maxat)
      dimension eecp(maxat),ect(maxat,maxat)
      dimension edip(maxat,maxat)
      dimension tkin(maxat),fdip(maxat,2)
      character*40 cfld
      logical lref
      common /printout/iaccur
      allocatable :: chp2(:,:),scr(:),chpaux(:,:), chp3(:,:)

      idofr   =  Iopt(40)
      ifield  =  Iopt(47)
      ipoints =  itotps
      iatps   =  nang*nrad
      xxdip   =  ZERO
      ndec=7
      if(iaccur.ne.0) ndec=13

      ALLOCATE(chp2(itotps,nalf),chp3(itotps,nb))

      call ao_to_mo_grid(itotps,igr,nalf,c,chp,chp2)
      call ao_to_mo_grid(itotps,igr,nb,cb,chp,chp3)

      epa=ZERO
      ekin=ZERO
      eto=ZERO
      enuc=ZERO
      ect=ZERO
      edip=ZERO

!! pseudopotential atomic energies (if present), added to the E-N terms !!
      call mga_misc(iecp,eecp)
      if(iecp.eq.1) then
        do i=1,nat
          epa(i,i)=eecp(i)
        end do
      end if

!! field(:)/ifield were already read by main.f into the /efield/ common  !!
!! and iopt(47) -- no need to re-scan the .inp file here.                !!
      if(ifield.eq.1) then
!! parallel over icenter: each iteration writes only its own dip(icenter,:)/ !!
!! ect(icenter,icenter); x/y/z/zztop/xtot are per-iteration scratch, not     !!
!! running totals, so no reduction is needed. Untested: no test uses an      !!
!! external electric field.                                                  !!
!$OMP PARALLEL DO PRIVATE(icenter,x,y,z,zztop,ifut,distx,disty,distz,xtot)
        do icenter=1,nat
          x=ZERO
          y=ZERO
          z=ZERO
          zztop=ZERO
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            distx=pcoord(ifut,1)
            disty=pcoord(ifut,2)
            distz=pcoord(ifut,3)
            x=x+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*distx
            y=y+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*disty
            z=z+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*distz
            zztop=zztop+wp(ifut)*omp2(ifut,icenter)*rho(ifut)
          end do
          dip(icenter,1)=-x+zztop*coord(1,icenter)
          dip(icenter,2)=-y+zztop*coord(2,icenter)
          dip(icenter,3)=-z+zztop*coord(3,icenter)
          xtot=(coord(1,icenter)*field(2)+coord(2,icenter)*field(3)+coord(3,icenter)*field(4))
          ect(icenter,icenter)=-zztop*xtot
        end do
!$OMP END PARALLEL DO

!! electric field x,y,z components are on field(2-4) !!
        xtot=ZERO
        do i=1,nat
          edip(i,i)=+dip(i,1)*field(2)+dip(i,2)*field(3)+dip(i,3)*field(4)
          xtot=xtot+ect(i,i)+edip(i,i)
        end do
        xxdip=xtot
      end if

!! electron-nuclei attraction integrals; parallel over the (icenter,jcenter) !!
!! atom-pair space (COLLAPSE(2), both bounds constant) -- each pair writes  !!
!! only its own epa(icenter,jcenter), xtot is a genuine running total.      !!
      xtot=ZERO
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,jcenter,x,zz,dx0,dy0,dz0,ifut,distx,disty,distz,rr) REDUCTION(+:xtot)
      do icenter=1,nat
        do jcenter=1,nat
          x=ZERO
          zz=zn(jcenter)
          dx0=coord(1,jcenter)
          dy0=coord(2,jcenter)
          dz0=coord(3,jcenter)
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            distx=pcoord(ifut,1)-dx0
            disty=pcoord(ifut,2)-dy0
            distz=pcoord(ifut,3)-dz0
            rr=dsqrt(distx*distx+disty*disty+distz*distz)
            if(rr.gt.1.0d-12) then
              x=x+zz*wp(ifut)*omp2(ifut,icenter)*rho(ifut)/rr
            end if
          end do
          xtot=xtot+x
          epa(icenter,jcenter)=epa(icenter,jcenter)-x
        end do
      end do
!$OMP END PARALLEL DO
      xtot=ZERO
      do i=1,nat
        xtot=xtot+epa(i,i)
        do j=1,i-1
          epa(i,j)=(epa(i,j)+epa(j,i))
          epa(j,i)=epa(i,j)
          xtot=xtot+epa(i,j)
        end do
      end do
      eelnuc=xtot

      if(iecp.eq.1) then
        call enpart_matrix('ELECTRON-NUCLEAR ATTRACTION','Includes the'//
     +    ' pseudopotential (ECP) energy of each atom on the diagonal.',
     +    ' ',epa,'Electron-nuclear attraction',eelnuc)
      else
        call enpart_matrix('ELECTRON-NUCLEAR ATTRACTION',' ',' ',epa,
     +    'Electron-nuclear attraction',eelnuc)
      end if

!! xxdip is added here (not just eelnuc) since it's excluded from the    !!
!! atomic epa contributions above -- keeps this check apples-to-apples.  !!
      lref=eelnuc0.ne.ZERO
      errk=ZERO
      if(lref) then
        errk=(eelnuc+xxdip-eelnuc0)*tokcal
        call print_value('Integration error',errk,2,'kcal/mol')
      else
        eelnuc0=eelnuc
      end if
      call enpart_frag('ELECTRON-NUCLEAR ATTRACTION',epa)
      call json_enpart_term('electron_nuclear','Electron-nuclear '//
     +  'attraction',epa,eelnuc,lref,errk)

!! second-derivative AO values, needed for the kinetic-energy density below !!
      ALLOCATE(chpaux(itotps,igr))
      call dpoints(chpaux,pcoord)

!! kinetic-energy density at each point: transform the 2nd-derivative AOs !!
!! to MOs (alpha+beta) and dot with the already-transformed chp2/chp3 --  !!
!! density in scr, /TWO like the restricted case's implicit factor of two !!
!! (kept unmodified from the original UHF formula). parallel over k: each !!
!! k only reads its own chpaux(k,:)/chp2(k,:)/chp3(k,:) and writes scr(k). !!
      ALLOCATE(scr(itotps))
!$OMP PARALLEL DO PRIVATE(k,x,j,xx,xxb,i)
      do k=1,ipoints
        x=ZERO
        do j=1,nalf
          xx=ZERO
          xxb=ZERO
          do i=1,igr
            xx=xx+c(i,j)*chpaux(k,i)
            if(j.le.nb) xxb=xxb+cb(i,j)*chpaux(k,i)
          end do
          x=x+chp2(k,j)*xx
          if(j.le.nb) x=x+chp3(k,j)*xxb
        end do
        scr(k)=x/TWO
      end do
!$OMP END PARALLEL DO

!! contract the kinetic-energy density onto atoms; parallel over icenter, !!
!! each iteration writes only its own ekin(icenter,icenter), xtot is a    !!
!! genuine running total.                                                 !!
      xtot=ZERO
!$OMP PARALLEL DO PRIVATE(icenter,x,ifut) REDUCTION(+:xtot)
      do icenter=1,nat
        x=ZERO
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          x=x+wp(ifut)*scr(ifut)*omp2(ifut,icenter)
        end do
        ekin(icenter,icenter)=ekin(icenter,icenter)+x
        xtot=xtot+x
      end do
!$OMP END PARALLEL DO
      ekinen=xtot

      do i=1,nat
        tkin(i)=ekin(i,i)
      end do
      call enpart_atom_table('KINETIC ENERGY','One-center terms only:'//
     +  ' the kinetic energy has no pair terms.',' ',1,tkin,'Kinetic',
     +  ' ',.false.)
      call print_value('Kinetic energy',ekinen,ndec,'au')
      lref=ekin0.ne.ZERO
      errk=ZERO
      if(lref) then
        errk=(ekinen-ekin0)*tokcal
        call print_value('Integration error',errk,2,'kcal/mol')
      else
        ekin0=ekinen
      end if
      if(idofr.eq.1) call group_by_frag_vec(1,'FRAGMENT KINETIC ENERGY',
     +  tkin,'Kinetic',' ')
!! the results file keeps the kinetic energy as a (diagonal) matrix, !!
!! like the other terms                                              !!
      call json_enpart_term('kinetic','Kinetic energy',ekin,ekinen,lref,
     +  errk)

      erep=ZERO
      do i=1,nat
        eto(i,i)=ekin(i,i)+epa(i,i)
        do j=i+1,nat
          dist=dsqrt((coord(1,i)-coord(1,j))**2 +(coord(2,i)-coord(2,j))**2+(coord(3,i)-coord(3,j))**2)
          eto(i,j)=ekin(i,j)+epa(i,j)+zn(i)*zn(j)/dist
          eto(j,i)=eto(i,j)

!! enuc: nuclear-repulsion pair matrix, kept only for the printout below. !!
          enuc(i,j)=zn(i)*zn(j)/dist
          enuc(j,i)=enuc(i,j)
          erep=erep+zn(i)*zn(j)/dist
        end do
      end do

      call enpart_matrix('NUCLEAR-NUCLEAR REPULSION',' ',' ',enuc,
     +  'Nuclear-nuclear repulsion',erep)
      call enpart_frag('NUCLEAR-NUCLEAR REPULSION',enuc)
      call json_enpart_term('nuclear_repulsion','Nuclear-nuclear '//
     +  'repulsion',enuc,erep,.false.,0.0d0)

      if(ifield.eq.1) then
        xxx=ZERO
        xtot=ZERO
        do i=1,nat
          xxx=zn(i)*(coord(1,i)*field(2)+coord(2,i)*field(3)+coord(3,i)*field(4))
          ect(i,i)=ect(i,i)+xxx
          xtot=xtot+xxx
        end do
        edipole=xtot+xxdip
        call print_box('EXTERNAL ELECTRIC FIELD')
        write(cfld,'(3f10.6)') field(2),field(3),field(4)
        call print_value('Field x, y, z',trim(adjustl(cfld)),'au')
        call print_value('Electronic dipole term',xxdip,ndec,'au')
        call print_value('Nuclear dipole term',xtot,ndec,'au')
        call print_value('Field interaction energy',edipole,ndec,'au')
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Not part of the IQA terms: added to their'//
     +    ' total in the energy check'
        write(*,'(2x,a)') 'of IQA ENERGY TERMS.'
        do i=1,nat
          fdip(i,1)=edip(i,i)
          fdip(i,2)=ect(i,i)
        end do
        call enpart_atom_table('FIELD INTERACTION PER ATOM','Intrinsic'//
     +    ' dipole term (origin independent) and charge-transfer',
     +    'term (origin dependent) of each atom.',2,fdip,'Dipole',
     +    'Charge tr.',.true.)
        if(idofr.eq.1) call group_by_frag_vec(2,
     +    'FRAGMENT FIELD INTERACTION',fdip,'Dipole','Charge tr.')
      end if

      etot=ekinen+eelnuc+erep
      etot0=ekin0+eelnuc0+erep

      DEALLOCATE(chp2,scr,chpaux,chp3)
      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: numint_two_uhf                                              !!
!! purpose: two-electron IQA/ENPART energy partition, UHF/UKS twin of      !!
!!   numint_two -- Coulomb (via calc_coul) and, if requested, HF-type      !!
!!   exchange (same-center and atom-pair kernels, alpha+beta contributions !!
!!   summed via the j.le.nb guard), each decomposed into atomic/diatomic   !!
!!   contributions. Same zero-error rotated-grid interpolation strategy as !!
!!   numint_two. See numint_two for restricted/closed-shell.               !!
!! arguments: same as numint_two.                                          !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine numint_two_uhf(ndim,itotps,wp,omp,omp2,pcoord,chp,rho,eto)

      use ao_matrices
      use integration_grid
      use OMP_LIB
      use timing_mod
      use output_mod, only: print_value
      use input_options_mod, only: iverbose

      implicit real*8(a-h,o-z)


      include 'parameter.h'
      parameter(maxpass=5)

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/twoel/twoeltoler
      common /iops/iopt(200)
      common/actual/iact,jat,icenter
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/energ0/ekin0,eelnuc0,evee0,etot0
      common/exchg/exch(maxat,maxat),xmix
      common/efield/field(4),edipole
!! rotated-grid angles/thresholds, read from the # GRID input section. !!
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22

      dimension eto(maxat,maxat),Excmp(maxat,maxat)
      dimension wp(itotps),chp(itotps,ndim),pcoord(itotps,3)
      dimension omp(itotps),omp2(itotps,nat),rho(itotps)

      dimension coul0(maxat,2*maxpass)
      allocatable :: coulb(:,:),exhfb(:,:),exchb(:,:)
      dimension coul(maxat,maxat),exch_hf(maxat,maxat)

      character*100 threadenv

      allocatable :: chppha(:,:), wppha(:),pcoordpha(:,:),omppha(:)
      allocatable :: chp2pha(:,:),omp2pha(:,:),ibaspointpha(:),rhopha(:)
      allocatable :: chp2(:,:), chp2b(:,:), chp2phab(:,:)
      allocatable :: Excmpb(:,:)

!! parallelization-related arrays below. !!
      allocatable :: chp2s(:,:),chp2phas(:,:),chp2bs(:,:),chp2phabs(:,:) !! transposed layout (rows/columns swapped) !!
      allocatable :: ijpaircount(:,:),istart(:),iend(:) !! included atom pairs, thread-slicing bookkeeping !!
      allocatable :: exch_hfk(:,:)
      allocatable :: exch_hfij(:,:)
      allocatable :: f3k(:)

      ihf         = Iopt(18)
      ithrebod    = Iopt(44)
      ipolar      = Iopt(46)
      ifield      = Iopt(47)
      ianalytical = Iopt(91)

      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if
      idoex=0
      if(xmix.gt.ZERO) idoex=1
      iatps=nang*nrad
!! zero-error strategy only possible with a reference electron-electron !!
!! energy in the .fchk (evee0, see mga_misc); twoelerr stays defined.    !!
      iref=0
      if(evee0.ne.ZERO) iref=1
      twoelerr=ZERO

      npass=2
!! HF: the atom pairs integrated and those from the multipolar expansion !!
!! (for KS-DFT the DFT part has already printed them)                   !!
      if(ihf.eq.1.and.ianalytical.eq.0) then
        npairs=0
        nskip=0
        do i=1,nat
          do j=i+1,nat
            if(bo(i,j).ge.threbod) then
              npairs=npairs+1
            else
              nskip=nskip+1
            end if
          end do
        end do
        call print_value('Pairs integrated (exact exchange)',npairs)
        call print_value('Pairs from the multipolar expansion',nskip)
      end if
      do i=1,nat
        do j=1,nat
          coul(i,j)=ZERO
          exch_hf(i,j)=ZERO
        end do
        do j=1,2*npass
          coul0(i,j)=ZERO
        end do
      end do

      ALLOCATE(chp2(itotps,nalf),chp2b(itotps,nb))
      ALLOCATE(chp2s(nalf,itotps),chp2bs(nb,itotps))

      call ao_to_mo_grid_t(itotps,igr,nalf,c,chp,chp2,chp2s)
      call ao_to_mo_grid_t(itotps,igr,nb,cb,chp,chp2b,chp2bs)

      if(ianalytical.eq.0) then !! not for the analytical path !!

!! rotated grid, angle controlled by the # GRID section -- see the diatXC !!
!! paper for optimized values: phb=0.162d0, later 0.182d0, for 40/146.   !!
        phb=phb12
        pha=ZERO
        ALLOCATE(wppha(itotps),omppha(itotps),omp2pha(itotps,nat))
        ALLOCATE(chppha(itotps,ndim),pcoordpha(itotps,3),ibaspointpha(itotps))
        ALLOCATE(chp2pha(itotps,nalf),rhopha(itotps),chp2phab(itotps,nb))
        ALLOCATE(chp2phas(nalf,itotps),chp2phabs(nb,itotps))

        call prenumint(ndim,itotps,nat,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)

        call ao_to_mo_grid_t(itotps,igr,nalf,c,chppha,chp2pha,chp2phas)
        call ao_to_mo_grid_t(itotps,igr,nb,cb,chppha,chp2phab,chp2phabs)

!! rhopha is recomputed from the MO amplitudes (rather than kept as       !!
!! prenumint's own density output) to stay numerically consistent with   !!
!! chp2pha/chp2phab, which the exchange kernel below reads directly.     !!
!! parallel over grid points: each k only reads its own chp2pha(k,:)/    !!
!! chp2phab(k,:) and writes only its own rhopha(k).                      !!
!$OMP PARALLEL DO PRIVATE(k,j,xtot,xtotb)
        do k=1,itotps
          xtot=ZERO
          xtotb=ZERO
          do j=1,nalf
            xtot=xtot+chp2pha(k,j)*chp2pha(k,j)
          end do
          do j=1,nb
            xtotb=xtotb+chp2phab(k,j)*chp2phab(k,j)
          end do
          rhopha(k)=xtot+xtotb
        end do
!$OMP END PARALLEL DO
      end if !! end of the non-analytical part !!

!! Coulomb energy term. !!
      call cpu_time(xtime)
      call get_wall_time(wxtime)
      if(ianalytical.eq.0) then
        call calc_coul(itotps,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,rho,rhopha,coul)
      else
        CALL calc_twoel_analytical(itotps,wp,omp2,pcoord,rho,chp2,chp2b,coul,exch_hf)
      end if
      call cpu_time(xtime2)
      call get_wall_time(wxtime2)
      write(*,'(a)') ''
      call print_timer('ENPART Coulomb',xtime2-xtime,wxtime2-wxtime)
      xtime=xtime2
      wxtime=wxtime2

!! HF-type exchange, computed only if requested (idoex). !!
      if(ianalytical.eq.0) then !! not for the analytical path !!
        if(idoex.eq.1) then

!! multipolar approximation, used below for atom pairs skipped by THREBOD, !!
!! computed for alpha and beta separately then summed. !!
!! (per spin: weight -1 for phi_i^2, -2 for phi_i*phi_j, i<j)           !!
          ALLOCATE(Excmpb(maxat,maxat))
          call multipolar_mo(nalf,itotps,wp,omp2,pcoord,chp2,-ONE,-TWO,Excmp)
          call multipolar_mo(nb,itotps,wp,omp2,pcoord,chp2b,-ONE,-TWO,Excmpb)
          do ii=1,nat
            do jj=ii+1,nat
              Excmp(ii,jj)=Excmp(ii,jj)+Excmpb(ii,jj)
              Excmp(jj,ii)=Excmp(ii,jj)
            end do 
          end do 
          DEALLOCATE(Excmpb)

!! loops reordered for parallelization purposes, mimicking numint_two's !!
!! strategy (see there for the full explanation). !!

!! same-center terms first; determine the core count for thread splitting. !!
          call getenv('OMP_NUM_THREADS',threadenv)
          if(trim(threadenv)=='') then
            ithreadenv=ZERO
          else
            read(unit=threadenv,FMT='(I4)') ithreadenv
          end if
!! a parallel region only to query the thread and processor counts;  !!
!! every thread writes the same values to iprocs/ithreads            !!
        !$OMP PARALLEL
          iprocs=OMP_GET_MAX_THREADS()
          ithreads=INT(OMP_GET_NUM_PROCS())
        !$OMP END PARALLEL
          if(ithreadenv.ne.ZERO) ithreads=ithreadenv

!! per-thread grid-point slicing bookkeeping. !!
          ALLOCATE(f3k(ithreads))
          ALLOCATE(exch_hfk(nat,ithreads))
          ALLOCATE(istart(ithreads),iend(ithreads))
          itilerest=mod(iatps,ithreads)
          ispace=(iatps-itilerest)/ithreads

!! loop kept in this manual-chunking shape rather than a plain OMP loop --  !!
!! avoids an old cap that limited parallelization to 8 cores.               !!
          exch_hfk=ZERO
          do icenter=1,nat
            ioffset=(icenter-1)*iatps
            istart=0
            iend=0
            do ik=1,(ithreads-1)
              istart(ik)=ioffset+((ik-1)*ispace)+1
              iend(ik)=ioffset+(ik*ispace)
            end do
            istart(ithreads)=ioffset+((ithreads-1)*ispace)+1
            iend(ithreads)=icenter*iatps

!! both calls below are required for the thread count to actually take effect. !!
            call omp_set_dynamic(.false.)
            call omp_set_num_threads(ithreads)
!! same real-OMP-plus-false-sharing fix as the RHF twin (numint_two)      !!
!! above -- f3loc is PRIVATE, exch_hfk written once per ifut, not per i/j.!!
!$OMP PARALLEL DO PRIVATE(ifut,jfut,x0,dx0,dy0,dz0,x1,dx1,dy1,dz1,dist,
!$OMP&  i,j,f2,f2b,f3loc)
            do ik=1,ithreads
              do ifut=istart(ik),iend(ik)
                x0=wp(ifut)*omp2(ifut,icenter)
                dx0=pcoord(ifut,1)
                dy0=pcoord(ifut,2)
                dz0=pcoord(ifut,3)
                f3loc=ZERO
                do jfut=iatps*(icenter-1)+1,iatps*icenter
                  x1=wppha(jfut)*omp2pha(jfut,icenter)
                  dx1=pcoordpha(jfut,1)
                  dy1=pcoordpha(jfut,2)
                  dz1=pcoordpha(jfut,3)
                  dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                  if(dist.gt.1.0d-12) then
                    do i=1,nalf-1
                      f2=x0*chp2s(i,ifut)*chp2s(i,ifut)
                      f3loc=f3loc+f2*x1*chp2phas(i,jfut)*chp2phas(i,jfut)/dist
                      if(i.le.nb) then
                        f2b=x0*chp2bs(i,ifut)*chp2bs(i,ifut)
                        f3loc=f3loc+f2b*x1*chp2phabs(i,jfut)*chp2phabs(i,jfut)/dist
                      end if
                      do j=i+1,nalf
                        f2=TWO*x0*chp2s(i,ifut)*chp2s(j,ifut)
                        f3loc=f3loc+f2*x1*chp2phas(i,jfut)*chp2phas(j,jfut)/dist
                        if(j.le.nb) then
                          f2b=TWO*x0*chp2bs(i,ifut)*chp2bs(j,ifut)
                          f3loc=f3loc+f2b*x1*chp2phabs(i,jfut)*chp2phabs(j,jfut)/dist
                        end if
                      end do
                    end do
                    f2=x0*chp2s(nalf,ifut)*chp2s(nalf,ifut)
                    f3loc=f3loc+f2*x1*chp2phas(nalf,jfut)*chp2phas(nalf,jfut)/dist

!! nalf=nb's own highest-index term needs adding separately (loop above !!
!! only pairs i with i+1..nalf, so the beta-side match never triggers). !!
                    if(nalf.eq.nb) then
                      f2b=x0*chp2bs(nb,ifut)*chp2bs(nb,ifut)
                      f3loc=f3loc+f2b*x1*chp2phabs(nb,jfut)*chp2phabs(nb,jfut)/dist
                    end if
                  end if
                end do
                exch_hfk(icenter,ik)=exch_hfk(icenter,ik)+f3loc
              end do
            end do
!$OMP END PARALLEL DO
          end do

!! reduce the per-thread accumulator into the same-center exch_hf term. !!
          do icenter=1,nat
            do isum=1,ithreads
              exch_hf(icenter,icenter)=exch_hf(icenter,icenter)-exch_hfk(icenter,isum)
            end do
          end do

!! atom-pair terms, skipping pairs below the THREBOD bond-order threshold. !!
          iterms=0
          ipaircounter=0
          ALLOCATE(ijpaircount(nat*nat,2))
          do icenter=1,nat
            do jcenter=icenter+1,nat
              bx0=bo(icenter,jcenter)
              if(bx0.ge.threbod) then
                ipaircounter=ipaircounter+1
                ijpaircount(ipaircounter,1)=icenter
                ijpaircount(ipaircounter,2)=jcenter
              else
                iterms=iterms+1
              end if
            end do
          end do

!! per-thread grid-point slicing bookkeeping for the atom-pair loop. !!
          ALLOCATE(exch_hfij(ipaircounter,ithreads))
          itilerest=mod(iatps,ithreads)
          ispace=(iatps-itilerest)/ithreads
          exch_hfij=ZERO
          do numpairnat=1,ipaircounter
            icenter=ijpaircount(numpairnat,1)
            jcenter=ijpaircount(numpairnat,2)
            ioffset=(icenter-1)*iatps
            istart=0
            iend=0
            do ik=1,(ithreads-1)
              istart(ik)=ioffset+((ik-1)*ispace)+1
              iend(ik)=ioffset+(ik*ispace)
            end do
            istart(ithreads)=ioffset+((ithreads-1)*ispace)+1
            iend(ithreads)=icenter*iatps
!! parallel over the thread slices ik, as the same-center block above: !!
!! f3loc is PRIVATE, exch_hfij(numpairnat,ik) written once per ifut.     !!
!$OMP PARALLEL DO PRIVATE(ifut,jfut,x0,dx0,dy0,dz0,x1,dx1,dy1,dz1,dist,
!$OMP&  i,j,f2,f2b,f3loc)
            do ik=1,ithreads
              do ifut=istart(ik),iend(ik)
                x0=wp(ifut)*omp2(ifut,icenter)
                dx0=pcoord(ifut,1)
                dy0=pcoord(ifut,2)
                dz0=pcoord(ifut,3)
                f3loc=ZERO
                do jfut=iatps*(jcenter-1)+1,iatps*jcenter
                  x1=wp(jfut)*omp2(jfut,jcenter)
                  dx1=pcoord(jfut,1)
                  dy1=pcoord(jfut,2)
                  dz1=pcoord(jfut,3)
                  dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                  if(dist.gt.1.0d-12) then
                    do i=1,nalf-1
                      f2=x0*chp2s(i,ifut)*chp2s(i,ifut)
                      f3loc=f3loc+f2*x1*chp2s(i,jfut)*chp2s(i,jfut)/dist
                      if(i.le.nb) then
                        f2b=x0*chp2bs(i,ifut)*chp2bs(i,ifut)
                        f3loc=f3loc+f2b*x1*chp2bs(i,jfut)*chp2bs(i,jfut)/dist
                      end if
                      do j=i+1,nalf
                        f2=TWO*x0*chp2s(i,ifut)*chp2s(j,ifut)
                        f3loc=f3loc+f2*x1*chp2s(i,jfut)*chp2s(j,jfut)/dist
                        if(j.le.nb) then
                          f2b=TWO*x0*chp2bs(i,ifut)*chp2bs(j,ifut)
                          f3loc=f3loc+f2b*x1*chp2bs(i,jfut)*chp2bs(j,jfut)/dist
                        end if
                      end do
                    end do
                    f2=x0*chp2s(nalf,ifut)*chp2s(nalf,ifut)
                    f3loc=f3loc+f2*x1*chp2s(nalf,jfut)*chp2s(nalf,jfut)/dist

!! same nalf=nb edge case as the same-center block above. !!
                    if(nalf.eq.nb) then
                      f2b=x0*chp2bs(nb,ifut)*chp2bs(nb,ifut)
                      f3loc=f3loc+f2b*x1*chp2bs(nb,jfut)*chp2bs(nb,jfut)/dist
                    end if
                  end if
                end do
                exch_hfij(numpairnat,ik)=exch_hfij(numpairnat,ik)+f3loc
              end do
            end do
!$OMP END PARALLEL DO
          end do

!! reduce the per-thread accumulator into the atom-pair exch_hf term. !!
          do numpairnat=1,ipaircounter
            icenter=ijpaircount(numpairnat,1)
            jcenter=ijpaircount(numpairnat,2)
            do isum=1,ithreads
              exch_hf(icenter,jcenter)=exch_hf(icenter,jcenter)-exch_hfij(numpairnat,isum)
            end do
          end do   

!! account for A-B/B-A symmetry (factor of 2) and for atom pairs whose !!
!! term came from the multipolar approximation instead of integration. !!
          exchen_hf=ZERO
          do i=1,nat
            exch_hf(i,i)=exch_hf(i,i)/TWO
            exchen_hf=exchen_hf+exch_hf(i,i)
            do j=i+1,nat
              if(bo(i,j).lt.threbod) exch_hf(i,j)=Excmp(i,j)
              exch_hf(j,i)=exch_hf(i,j)
              exchen_hf=exchen_hf+exch_hf(i,j)
            end do
          end do
          call cpu_time(xtime2)
          call get_wall_time(wxtime2)
          call print_timer('ENPART exact exchange',xtime2-xtime,
     +      wxtime2-wxtime)
          xtime=xtime2
          wxtime=wxtime2
        end if
      end if !! non-analytical skip ends here !!

!! electron-electron energy on the first grid and its integration error; !!
!! for a hybrid functional the exact-exchange share joins the DFT part   !!
      if(ihf.eq.1) then
        evee=coulen+exchen_hf
      else
        exchen=ZERO
        do i=1,nat
          exch(i,i)=exch(i,i)+xmix*exch_hf(i,i)
          exchen=exchen+exch(i,i)
          do j=i+1,nat
            exch(i,j)=exch(i,j)+xmix*exch_hf(i,j)
            exch(j,i)=exch(i,j)
            exchen=exchen+exch(i,j)
          end do
        end do
        evee=coulen+exchen
      end if
      if(evee0.ne.ZERO) then
        twoelerr=(evee-evee0)*tokcal
      else
        evee0=evee
      end if

      if(ianalytical.eq.1) twoeltoler=232000
      call enpart_ee_head(evee,twoelerr,iref,twoeltoler)
!! first-grid values, kept for the results file !!
      eveefg=evee
      errfg=twoelerr
      err2=ZERO
      phabest=ZERO

!! the zero-error correction changes only the one-center terms; the      !!
!! first-grid matrices are kept for VERBOSE 1                            !!
      izerr=0
      if(iref.eq.1.and.abs(twoelerr).gt.twoeltoler) izerr=1
      if(izerr.eq.1.and.iverbose.eq.1) then
        ALLOCATE(coulb(maxat,maxat),exhfb(maxat,maxat))
        ALLOCATE(exchb(maxat,maxat))
        do i=1,nat
          do j=1,nat
            coulb(i,j)=coul(i,j)
            exhfb(i,j)=exch_hf(i,j)
            exchb(i,j)=exch(i,j)
          end do
        end do
        xcoulb=coulen
        xexhfb=exchen_hf
        xexchb=exchen
      end if

      if(izerr.eq.1) then

!! rotated grid, angle controlled by the # GRID section -- see the diatXC !!
!! paper for optimized values: phb=0.162d0, later 0.182d0, for 40/146.   !!
        phb=phb22
        pha=ZERO
        nat0=nat
        call prenumint(ndim,itotps,nat0,wppha,omppha,omp2pha,chppha,rhopha,pcoordpha,ibaspointpha,0)

        call ao_to_mo_grid_t(itotps,igr,nalf,c,chppha,chp2pha,chp2phas)
        call ao_to_mo_grid_t(itotps,igr,nb,cb,chppha,chp2phab,chp2phabs)

!! rhopha recomputed from the MO amplitudes, same reason as the earlier !!
!! rotated-grid pass above. !!
!$OMP PARALLEL DO PRIVATE(k,j,xtot,xtotb)
        do k=1,itotps
          xtot=ZERO
          xtotb=ZERO
          do j=1,nalf
            xtot=xtot+chp2phas(j,k)*chp2phas(j,k)
          end do
          do j=1,nb
            xtotb=xtotb+chp2phabs(j,k)*chp2phabs(j,k)
          end do
          rhopha(k)=xtot+xtotb
        end do
!$OMP END PARALLEL DO

!! same-center-only recompute of the Coulomb part on the rotated grid, !!
!! for the zero-error interpolation below. parallel over icenter: each !!
!! iteration accumulates into its own coul0(icenter,3), independent of !!
!! every other icenter.                                                !!
!$OMP PARALLEL DO PRIVATE(icenter,ifut,x0,dx0,dy0,dz0,f3,jfut,x1,dx1,dy1,dz1,dist)
        do icenter=1,nat
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            x0=wp(ifut)*omp2(ifut,icenter)
            dx0=pcoord(ifut,1)
            dy0=pcoord(ifut,2)
            dz0=pcoord(ifut,3)
            f3=ZERO
            do jfut=iatps*(icenter-1)+1,iatps*icenter
              x1=wppha(jfut)*omp2pha(jfut,icenter)
              dx1=pcoordpha(jfut,1)
              dy1=pcoordpha(jfut,2)
              dz1=pcoordpha(jfut,3)
              dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
              if(dist.gt.1.0d-12) f3=f3+(rho(ifut)*rhopha(jfut))*x1*x0/dist
            end do
            coul0(icenter,3)=coul0(icenter,3)+f3/TWO
          end do
        end do
!$OMP END PARALLEL DO

!! HF-type exchange, computed only if requested (idoex). !!
        if(idoex.eq.1) then

!! same-center recompute, same manual-chunking shape as above. !!
          exch_hfk=ZERO
          do icenter=1,nat
            ioffset=(icenter-1)*iatps
            istart=0
            iend=0
            do ik=1,(ithreads-1)
              istart(ik)=ioffset+((ik-1)*ispace)+1
              iend(ik)=ioffset+(ik*ispace)
            end do
            istart(ithreads)=ioffset+((ithreads-1)*ispace)+1
            iend(ithreads)=icenter*iatps
  
!! both calls below are required for the thread count to actually take effect. !!
            call omp_set_dynamic(.false.)
            call omp_set_num_threads(ithreads)
!! parallel over the thread slices ik, as the blocks above: f3loc is   !!
!! PRIVATE, exch_hfk(icenter,ik) written once per ifut.                 !!
!$OMP PARALLEL DO PRIVATE(ifut,jfut,x0,dx0,dy0,dz0,x1,dx1,dy1,dz1,dist,
!$OMP&  i,j,f2,f2b,f3loc)
            do ik=1,ithreads
              do ifut=istart(ik),iend(ik)
                x0=wp(ifut)*omp2(ifut,icenter)
                dx0=pcoord(ifut,1)
                dy0=pcoord(ifut,2)
                dz0=pcoord(ifut,3)
                f3loc=ZERO
                do jfut=iatps*(icenter-1)+1,iatps*icenter
                  x1=wppha(jfut)*omp2pha(jfut,icenter)
                  dx1=pcoordpha(jfut,1)
                  dy1=pcoordpha(jfut,2)
                  dz1=pcoordpha(jfut,3)
                  dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
                  if(dist.gt.1.0d-8) then !! 1e-8 here, 1e-12 in the other loops !!
                    do i=1,nalf-1
                      f2=x0*chp2s(i,ifut)*chp2s(i,ifut)
                      f3loc=f3loc+f2*x1*chp2phas(i,jfut)*chp2phas(i,jfut)/dist
                      if(i.le.nb) then
                        f2b=x0*chp2bs(i,ifut)*chp2bs(i,ifut)
                        f3loc=f3loc+f2b*x1*chp2phabs(i,jfut)*chp2phabs(i,jfut)/dist
                      end if
                      do j=i+1,nalf
                        f2=TWO*x0*chp2s(i,ifut)*chp2s(j,ifut)
                        f3loc=f3loc+f2*x1*chp2phas(i,jfut)*chp2phas(j,jfut)/dist
                        if(j.le.nb) then
                          f2b=TWO*x0*chp2bs(i,ifut)*chp2bs(j,ifut)
                          f3loc=f3loc+f2b*x1*chp2phabs(i,jfut)*chp2phabs(j,jfut)/dist
                        end if
                      end do
                    end do
                    f2=x0*chp2s(nalf,ifut)*chp2s(nalf,ifut)
                    f3loc=f3loc+f2*x1*chp2phas(nalf,jfut)*chp2phas(nalf,jfut)/dist

!! same nalf=nb edge case as the same-center block above. !!
                    if(nalf.eq.nb) then
                      f2b=x0*chp2bs(nb,ifut)*chp2bs(nb,ifut)
                      f3loc=f3loc+f2b*x1*chp2phabs(nb,jfut)*chp2phabs(nb,jfut)/dist
                    end if
                  end if
                end do
                exch_hfk(icenter,ik)=exch_hfk(icenter,ik)+f3loc
              end do
            end do
!$OMP END PARALLEL DO
          end do

  !! reduce the per-thread accumulator into the same-center coul0 term. !!
          do icenter=1,nat
            do isum=1,ithreads
              coul0(icenter,4)=coul0(icenter,4)-exch_hfk(icenter,isum)/TWO
            end do
          end do
        end if
  
!! compare same-center coul0/exch_hf against the rotated-grid recompute !!
!! to get the damping parameter for the interpolation below. !!
        do i=1,nat
          coul0(i,1)=coul(i,i)
          if(idoex.eq.1) coul0(i,2)=exch_hf(i,i)
        end do
        deltaee=ZERO
        do i=1,nat
          deltaee=deltaee+coul0(i,1)-coul0(i,3)
          if(idoex.eq.1)deltaee=deltaee+xmix*(coul0(i,2)-coul0(i,4))
        end do
        deltaee=deltaee*tokcal
        phabest=ONE-(twoelerr/deltaee)
        err1=twoelerr
        err2=twoelerr-deltaee

!! interpolate energies, replacing the old Coulomb/exchange terms !!
        do i=1,nat
          coul(i,i)=coul0(i,1)*phabest+(ONE-phabest)*coul0(i,3)
          if(idoex.eq.1) exch_hf(i,i)=coul0(i,2)*phabest+(ONE-phabest)*coul0(i,4)
        end do

        coulen=ZERO
        do i=1,nat
          do j=i,nat
            coulen=coulen+coul(i,j)
          end do
        end do
        if(idoex.eq.1) then
          exchen_hf=ZERO
          do i=1,nat
            do j=i,nat
              exchen_hf=exchen_hf+exch_hf(i,j)
            end do
          end do
          if(ihf.ne.1) then
            exchen=ZERO
            do i=1,nat
              exch(i,i)=exch(i,i)+xmix*(exch_hf(i,i)-coul0(i,2))
              exchen=exchen+exch(i,i)
              do j=i+1,nat
                exch(j,i)=exch(i,j)
                exchen=exchen+exch(i,j)
              end do
            end do
          end if
        else
          exchen=ZERO
          do i=1,nat
            exchen=exchen+exch(i,i)
            do j=i+1,nat
              exchen=exchen+exch(i,j)
            end do
          end do
        end if

        if(ihf.eq.0) then
          evee=coulen+exchen
        else
          evee=coulen+exchen_hf
        end if
        twoelerr=(evee-evee0)*tokcal
        call enpart_ee_tail(err1,err2,phabest,evee,twoelerr,.false.)
        write(*,'(a)') ''
        call cpu_time(xtime2)
        call get_wall_time(wxtime2)
        call print_timer('ENPART zero-error correction',xtime2-xtime,
     +    wxtime2-wxtime)
        xtime=xtime2
        wxtime=wxtime2
      end if

      call json_enpart_ee(eveefg,errfg,iref,twoeltoler,izerr,err2,
     +  phabest,evee,twoelerr)

!! final two-electron matrices; with VERBOSE 1 first the first-grid ones !!
!! when the zero-error correction was applied                            !!
      if(izerr.eq.1.and.iverbose.eq.1) call enpart_two_mats(ihf,xmix,
     +  coulb,xcoulb,exhfb,xexhfb,exchb,xexchb,
     +  ', BEFORE THE ZERO-ERROR CORRECTION')
      if(allocated(coulb)) DEALLOCATE(coulb,exhfb,exchb)
      call enpart_two_mats(ihf,xmix,coul,coulen,exch_hf,exchen_hf,exch,
     +  exchen,' ')

!! IQA matrix: one-electron terms (in eto) + Coulomb + exchange(-corr.) !!
      if(ihf.eq.1) then
        xtot=ZERO
        do i=1,nat
          eto(i,i)=eto(i,i)+coul(i,i)+exch_hf(i,i)
          xtot=xtot+eto(i,i)
          do j=i+1,nat
            eto(i,j)=eto(i,j)+coul(i,j)+exch_hf(i,j)
            eto(j,i)=eto(i,j)
            xtot=xtot+eto(i,j)
          end do
        end do
        etot=xtot
        call enpart_iqa(eto,etot,escf,ifield,edipole,exch_hf,.true.)
      else
        xtot=ZERO
        do i=1,nat
          eto(i,i)=eto(i,i)+coul(i,i)+xmix*exch_hf(i,i)
          xtot=xtot+eto(i,i)
          do j=i+1,nat
            eto(i,j)=eto(i,j)+coul(i,j)+xmix*exch_hf(i,j)
            eto(j,i)=eto(i,j)
            xtot=xtot+eto(i,j)
          end do
        end do
        etot=xtot
        call enpart_iqa(eto,etot,escf,ifield,edipole,exch,.true.)
      end if

      if(ianalytical.eq.0) then
        DEALLOCATE(chp2,chp2b,chp2pha,rhopha,chp2phab)
        DEALLOCATE(wppha,omppha,omp2pha,chppha,pcoordpha)
        DEALLOCATE(chp2s,chp2bs,chp2phas,chp2phabs)
        if(idoex.eq.1) then
          DEALLOCATE(f3k,exch_hfk,istart,iend,ijpaircount,exch_hfij)
        end if
      else
        DEALLOCATE(chp2,chp2b)
      end if
      end

!! ***** !!

      subroutine polar(itotps,nat0,wp,omp,omp2,pcoord,rho)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in) :: itotps,nat0
      include 'parameter.h'
      parameter (maxbonds=300)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /atomrad/atr(maxat),dist(maxat,maxat)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(200)
      common/actual/iact,jat,icenter
      common/efield/field(4),edipole
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension wp(itotps),pcoord(itotps,3),rho(itotps)
      dimension omp(itotps),omp2(itotps,nat0)
      dimension dipint(maxat,3)
      dimension eint(maxat),dip(maxat,3),dipct(3),diptot(3)         
c MM
      dimension itracker(maxbonds),jtracker(maxbonds)
      real*8 MAT
      dimension MAT(maxbonds,maxbonds+1),ichi(maxat,maxat)
      dimension csolv(maxbonds),qqnorm(maxat),ct(maxat,3)
      dimension scr(3),iringatoms(20)
      integer unique_qij,B(maxat,maxat)
      logical ilog
c

      idofr   =   iopt(40)
      ifield  =   iopt(47)

      ipoints=itotps
      iatps=nang*nrad

      do i=1,3
        dipct(i)=0.0d0
      end do


       write(*,*)
       print *,' --------------------------------------'
       print *,'       DIPOLE MOMENT DECOMPOSITION'
       print *,' --------------------------------------'
       write(*,*)
       xtot=0.0d0
       do icenter=1,nat
        x=0.0d0
        y=0.0d0
        z=0.0d0
        zztop=0.0d0
        do ifut=iatps*(icenter-1)+1,iatps*icenter
         distx=pcoord(ifut,1)
         disty=pcoord(ifut,2)
         distz=pcoord(ifut,3)
         x=x+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*distx
         y=y+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*disty
         z=z+wp(ifut)*omp2(ifut,icenter)*rho(ifut)*distz
         zztop=zztop+wp(ifut)*omp2(ifut,icenter)*rho(ifut)
        end do
        dip(icenter,1)=-x+zn(icenter)*coord(1,icenter)    
        dip(icenter,2)=-y+zn(icenter)*coord(2,icenter)    
        dip(icenter,3)=-z+zn(icenter)*coord(3,icenter)    
        dipint(icenter,1)=-x+zztop*coord(1,icenter)    
        dipint(icenter,2)=-y+zztop*coord(2,icenter)    
        dipint(icenter,3)=-z+zztop*coord(3,icenter)    
        do j=1,3
         dipct(j)=dipct(j)+(zn(icenter)-zztop)*coord(j,icenter)    
        end do
       enddo
       write(*,*)
       write(*,*) 'Origin-dependent atomic dipole moments'
       call mprint_nlop(dip,nat,maxat)
       write(*,*)
       write(*,*) 'Total dipole moment: '
       do j=1,3  
        diptot(j)=0.0d0
        do i=1,nat
         diptot(j)=diptot(j)+dip(i,j)
        end do
       end do
       write(*,*) '               X                   Y                      Z' 
       write(*,'(3f24.15)') (diptot(i),i=1,3)
       write(*,*)
C
       write(*,*)
       write(*,*) 'Intrinsic atomic dipole moments'
       call mprint_nlop(dipint,nat,maxat)
       write(*,*)
       write(*,*) 'Total intrinsic dipole moment: '
       do j=1,3  
        diptot(j)=0.0d0
        do i=1,nat
         diptot(j)=diptot(j)+dipint(i,j)
        end do
       end do
       write(*,*) '               X                   Y                      Z' 
       write(*,'(3f24.15)') (diptot(i),i=1,3)
       write(*,*)
C
C Generalized Keith algorithm to decompose the charge-transfer term
C Code by MMontilla

C Recalculating geometrical connectivities 
C in case TFVC is not used

       do iatom=1,nat-1
         do jatom=iatom+1,nat
           ichi(iatom,jatom)=1
           rijx2=(coord(1,jatom)+coord(1,iatom))/2.0d0
           rijy2=(coord(2,jatom)+coord(2,iatom))/2.0d0
           rijz2=(coord(3,jatom)+coord(3,iatom))/2.0d0
           disth=dist(iatom,jatom)/2.0d0
           do ii=1,nat
             if(ii.ne.iatom.and.ii.ne.jatom) then
               xdist=sqrt((rijx2-coord(1,ii))**2+(rijy2-coord(2,ii))**2+(rijz2-coord(3,ii))**2)
               if(xdist.lt.disth) then
                 ichi(iatom,jatom)=0
                 go to 123
               end if
             end if
           end do
123        continue
         end do
       end do

       do i=1,nat
         do j=i,nat
           B(i,j)=ichi(i,j)
           B(j,i)=-B(i,j)
         end do
       end do

C Finding unique q_ij contributions, aka Q(A|B) contributions.
       unique_qij=0
       do i=1,Nat
         do j=i+1,Nat
           unique_qij=unique_qij+B(i,j)
         end do
       end do
       write(*,*) 'Unique Q(A|B) contributions found:',unique_qij

C Initializing coeff. matrix.
       do i=1,unique_qij
         do j=1,unique_qij
           MAT(i,j)=0.0d0
         end do
       end do

C Checking for lack of equations (unique_qij=number of unknowns;nat-1=number of independent equations
C ;excess unknowns require additional equations)
       iexcess=unique_qij-nat+1
       if(iexcess.ne.0) then 
         write(*,*) 'WARNING, CYCLIC CONNECTIVITY.'
         write(*,*) 'Will need ',iexcess,'additional equations.'
         write(*,*)

         open(45,FILE="ichi.inp",status="old",iostat=ii)
         if(ii.ne.0) then
           write(*,*) 'Printing connectivity matrix:'
           do i=1,nat
             write(*,'(40i2)') (b(i,j),j=1,nat)
           end do
           call apost_stop('error in polar, ichi.inp file is needed')
         end if
         read(45,'(a80)') 
         read(45,*) N

C Introducing (-N) modifications on the connectiviy matrix.
         if(N.lt.0) then
           write(*,*) 'Removing ',abs(N),' bonds.'
           do k=1,abs(N)
             read(45,*) i,j
             B(i,j)=0
             B(j,i)=0
             unique_qij=unique_qij-1
             write(*,*) 'Removed bond between atoms',i,j            
           end do
         end if
       end if

C Labeling the unique q_ij, or Q(A|B), contributions.
       l=1
       do i=1,Nat
         do j=i+1,Nat
           if(B(i,j).eq.1.0) then
             itracker(l)=i
             jtracker(l)=j
             write(*,*) 'Tracker pair ',l,':',i,j
             l=l+1
           end if
         end do
       end do
       if(l-1.ne.unique_qij) call apost_stop('problem in polar')
C End of modifications. Final connectiviy matrix.

C Additional equations for rings if needed. The same iexcess check is performed here,
C  because the loop is never needed if iexcess is 0.
       iextra=0
       if(iexcess.ne.0) then 
         read(45,*) N
         if(N.ne.0) then
           write(*,*) 'Reading ',N,'additional ring equations.'
           do ii=1,N
             read(45,*) iringsize
             read(45,*) (iringatoms(k),k=1,iringsize)
             iextra=iextra+1
             do i=1,iringsize-1
               j=i+1
               do m=1,unique_qij
                 if(iringatoms(i).eq.itracker(m).and.iringatoms(j).eq.jtracker(m)) then
                   Mat(iextra,m)=1.0d0
                 else if(iringatoms(j).eq.itracker(m).and.iringatoms(i).eq.jtracker(m)) then
                   Mat(iextra,m)=-1.0d0
                 end if
               end do
             end do
C Particular case (connectivity between last and first atoms in the ring). 
             i=iringsize
             j=1
             do m=1,unique_qij
               if(iringatoms(i).eq.itracker(m).and.iringatoms(j).eq.jtracker(m)) then
                 Mat(iextra,m)=1.0d0
               else if(iringatoms(j).eq.itracker(m).and.iringatoms(i).eq.jtracker(m)) then
                 Mat(iextra,m)=-1.0d0
               end if
             end do
C Condition for cycles. 
             Mat(iextra,unique_qij+1)=0.0d0
           end do
         end if
       end if
C End reading.
121    if(iexcess.ne.0) close(45)

C Building the matrix of coefficients.
       do i=1,Nat
         do j=1,Nat
           if(B(i,j).eq.1.or.B(i,j).eq.-1) then
             do m=1,unique_qij
               if(i.eq.itracker(m).and.j.eq.jtracker(m).or.j.eq.itracker(m).and.i.eq.jtracker(m)) then
                 Mat(i+iextra,m)=B(i,j)
               end if
             end do
           end if
         end do
       end do

       if(iextra+nat-1.ne.unique_qij) call apost_stop('problem building coeff matrix in polar')

c renormalizing atomic populations
        xx=0.0d0
        do i=1,nat
         xx=xx+qat(i,1)
        end do
        do i=1,nat
         qqnorm(i)=zn(i)-qat(i,1)*(nalf+nb)/xx
        end do
        
        do i=1,nat       
          Mat(i+iextra,unique_qij+1)=qqnorm(i)
        end do

        write(*,*) 'Echoing matrix of coefficients.'
        do i=1,Nat
            write(*,'(8f10.6)') (Mat(i,j),j=1,unique_qij+1)
        end do

        call SOLVESYSTEM(unique_qij,maxbonds,MAT,csolv,ilog)
        if(ilog.eqv..false.) call apost_stop('error in polar, CT decomposition')


C Echoing the solutions of the system of equations, for testing purposes.
c        do i=1,unique_qij
c            write(*,*) 'Unkown number',i,'is',csolv(i)
c        end do

C Final calculations.
       do i=1,nat
         do j=1,3
           CT(i,j)=0.0d0
           scr(j)=0.0d0
         end do
       end do

       do i=1,Nat
         do j=1,Nat
           if(B(i,j).NE.0) then
             Rabx=(coord(1,i)-coord(1,j))/2.0d0
             Raby=(coord(2,i)-coord(2,j))/2.0d0
             Rabz=(coord(3,i)-coord(3,j))/2.0d0
             do m=1,unique_qij
               if(i.eq.itracker(m).and.j.eq.jtracker(m).or.j.eq.itracker(m).and.i.eq.jtracker(m)) l=m
             end do
             CT(i,1)=CT(i,1)+B(i,j)*csolv(l)*Rabx
             CT(i,2)=CT(i,2)+B(i,j)*csolv(l)*Raby
             CT(i,3)=CT(i,3)+B(i,j)*csolv(l)*Rabz
           end if
         end do
       end do

C Echoing final Keith results.
       write(*,*)
       write(*,*) 'Charge-transfer atomic dipole moments (renormalized charges):'
       call mprint_nlop(ct,nat,maxat)
       write(*,*)
       write(*,*) 'Total CT dipole moment (renormalized charges):'
        do i=1,Nat
         do j=1,3
          scr(j)=scr(j)+CT(i,j)
         end do
        end do
        write(*,'(3f24.15)') (scr(i),i=1,3)
       write(*,*)
       write(*,*) 'Total CT dipole moment (unscaled charges):'
       write(*,*) '               X                   Y                      Z' 
       write(*,'(3f24.15)') (dipct(i),i=1,3)

C PSS
       write(*,*)
       write(*,*) 'Total dipole moment:'
       write(*,*) '               X                   Y                      Z' 
       write(*,'(3f24.15)') (dipct(i)+diptot(i),i=1,3)
       write(*,*)

c  energetics 
       if(ifield.eq.1) then
       dipen=0.0d0
       do icenter=1,nat
        eint(icenter)=dip(icenter,1)*field(2)+dip(icenter,2)*field(3)+dip(icenter,3)*field(4)  
        dipen=dipen+eint(icenter)
       end do
       cten=cten+dipct(1)*field(2)+dipct(2)*field(3)+dipct(3)*field(4)

       write(*,*) 'Intrinsic atomic dipole energy terms :'
       do i=1,nat
        write(*,'(a6,i4,f22.15)') 'Atom: ',i,eint(i)
       end do
       write(*,*)
       write(*,*) 'Total Intrinsic dipole energy term :',dipen

       write(*,*) 'Total Charge-transfer  energy term :',cten
       edipole=cten+dipen
       write(*,*) 'Total dipole energy term :',edipole

       end if

       end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: multipolar                                                  !!
!! purpose: multipolar-expansion estimate of the exchange(-correlation)    !!
!!   energy of every atom pair, used for pairs skipped by the THREBOD     !!
!!   bond-order threshold. Each pair "density" f_ij of the xc density is  !!
!!   expanded in per-atom charge/dipole/quadrupole moments, and the        !!
!!   pairwise terms are summed up to quadrupole-quadrupole. This version  !!
!!   takes the f_ij stored on the grid (CASSCF/CISD: eigenfunctions of     !!
!!   the xc density); see multipolar_mo for products of MOs.               !!
!! arguments:                                                              !!
!!   norb   (in)  -- number of orbitals, norb*(norb+1)/2 pair functions    !!
!!   itotps (in)  -- total number of grid points                          !!
!!   wp     (in)  -- integration weight of each grid point                 !!
!!   omp2   (in)  -- becke/tfvc (or hirshfeld) weight of each point for    !!
!!                   every atom                                            !!
!!   pcoord (in)  -- xyz coordinates of each grid point                    !!
!!   fij    (in)  -- pair functions on the grid, (itotps,norb2)            !!
!!   xocc   (in)  -- their weights (eigenvalues); only the diagonal is read !!
!!   Excmp  (out) -- multipolar estimate per atom pair (i/=j elements)     !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine multipolar(norb,itotps,wp,omp2,pcoord,fij,xocc,Excmp)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop

      dimension :: Excmp(maxat,maxat),wp(itotps),pcoord(itotps,3)
      dimension :: fij(itotps,norb*(norb+1)/2),omp2(itotps,nat)
      dimension :: xocc(norb*(norb+1)/2,norb*(norb+1)/2)

      allocatable :: dip(:,:),quadp(:,:,:),sij(:),atdist(:,:),exk(:,:,:)
      allocatable :: xones(:)

      ALLOCATE(dip(nat,3),quadp(nat,3,3),sij(nat),atdist(nat,nat))
      ALLOCATE(exk(nat,nat,6),xones(itotps))
      xones=ONE
      call mp_init(atdist,exk,Excmp)

!! fij*1 is exact, so these moments are bit-identical to using fij alone. !!
      irun=0
      do i=1,norb
        do j=i,norb
          irun=irun+1
          call mp_moments(itotps,wp,omp2,pcoord,fij(1,irun),xones,sij,dip,quadp)
          call mp_pairs(xocc(irun,irun),sij,dip,quadp,atdist,exk)
        end do
      end do

      call mp_sum(exk,Excmp)
      DEALLOCATE(dip,quadp,sij,atdist,exk,xones)

      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: multipolar_mo                                               !!
!! purpose: same multipolar estimate as multipolar, for a single           !!
!!   determinant, whose pair functions are the MO products phi_i*phi_j.   !!
!!   Each product is formed inside the moment loop rather than stored, so  !!
!!   memory is itotps*norb instead of itotps*norb*(norb+1)/2.              !!
!! arguments:                                                              !!
!!   norb   (in)  -- number of occupied MOs (nocc, nalf or nb)             !!
!!   itotps (in)  -- total number of grid points                          !!
!!   wp, omp2, pcoord (in) -- as in multipolar                             !!
!!   chpmo  (in)  -- MO values at each grid point, (itotps,norb)           !!
!!   ffd    (in)  -- weight of the i=j products (-2 closed shell, -1 per   !!
!!                   spin)                                                 !!
!!   ffo    (in)  -- weight of the i<j products (-4 closed shell, -2 per   !!
!!                   spin)                                                 !!
!!   Excmp  (out) -- multipolar estimate per atom pair (i/=j elements)     !!
!! author: MGimf                                                           !!
!! *********************************************************************** !!
      subroutine multipolar_mo(norb,itotps,wp,omp2,pcoord,chpmo,ffd,ffo,Excmp)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop

      dimension :: Excmp(maxat,maxat),wp(itotps),pcoord(itotps,3)
      dimension :: chpmo(itotps,norb),omp2(itotps,nat)

      allocatable :: dip(:,:),quadp(:,:,:),sij(:),atdist(:,:),exk(:,:,:)

      ALLOCATE(dip(nat,3),quadp(nat,3,3),sij(nat),atdist(nat,nat))
      ALLOCATE(exk(nat,nat,6))
      call mp_init(atdist,exk,Excmp)

      do i=1,norb
        do j=i,norb
          ffact=ffo
          if(i.eq.j) ffact=ffd
          call mp_moments(itotps,wp,omp2,pcoord,chpmo(1,i),chpmo(1,j),sij,dip,quadp)
          call mp_pairs(ffact,sij,dip,quadp,atdist,exk)
        end do
      end do

      call mp_sum(exk,Excmp)
      DEALLOCATE(dip,quadp,sij,atdist,exk)

      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: mp_init                                                     !!
!! purpose: interatomic distances and zeroed accumulators for              !!
!!   multipolar/multipolar_mo.                                             !!
!! arguments:                                                              !!
!!   atdist (out) -- distance between every pair of atoms                  !!
!!   exk    (out) -- the six multipole-order accumulators, zeroed          !!
!!   Excmp  (out) -- off-diagonal elements zeroed                          !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine mp_init(atdist,exk,Excmp)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      dimension :: atdist(nat,nat),exk(nat,nat,6),Excmp(maxat,maxat)

      do iat=1,nat
        atdist(iat,iat)=ZERO
        do jat=iat+1,nat
          xx=ZERO
          do ii=1,3
            xx=xx+(coord(ii,jat)-coord(ii,iat))**2
          end do
          atdist(iat,jat)=dsqrt(xx)
          atdist(jat,iat)=atdist(iat,jat)
          Excmp(iat,jat)=ZERO
          do k=1,6
            exk(iat,jat,k)=ZERO
          end do
        end do
      end do

      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: mp_moments                                                  !!
!! purpose: per-atom charge, dipole and (traceless) quadrupole moments of  !!
!!   one pair function f = fa*fb, each atom integrated on its own grid     !!
!!   points with its own weight.                                           !!
!! arguments:                                                              !!
!!   itotps, wp, omp2, pcoord (in) -- as in multipolar                     !!
!!   fa, fb (in)  -- the two factors of the pair function on the grid     !!
!!   sij    (out) -- charge of each atom                                   !!
!!   dip    (out) -- dipole of each atom                                   !!
!!   quadp  (out) -- quadrupole of each atom                               !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine mp_moments(itotps,wp,omp2,pcoord,fa,fb,sij,dip,quadp)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      dimension :: wp(itotps),pcoord(itotps,3),omp2(itotps,nat)
      dimension :: fa(itotps),fb(itotps)
      dimension :: sij(nat),dip(nat,3),quadp(nat,3,3)

      iatps = nang*nrad

!! parallel over icenter: each iteration writes only its own dip(icenter,:), !!
!! quadp(icenter,:,:) and sij(icenter), all independent across atoms.    !!
!$OMP PARALLEL DO PRIVATE(icenter,xx,yy,zz,xy,xz,yz,x,y,z,ifut,distx,disty,distz,wccij)
      do icenter=1,nat

        xx=ZERO
        yy=ZERO
        zz=ZERO
        xy=ZERO
        xz=ZERO
        yz=ZERO

        x=ZERO
        y=ZERO
        z=ZERO

        sij(icenter)=ZERO
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          distx=pcoord(ifut,1)-coord(1,icenter)
          disty=pcoord(ifut,2)-coord(2,icenter)
          distz=pcoord(ifut,3)-coord(3,icenter)
          wccij=wp(ifut)*omp2(ifut,icenter)*(fa(ifut)*fb(ifut))
          xx=xx+wccij*distx*distx
          yy=yy+wccij*disty*disty
          zz=zz+wccij*distz*distz
          xy=xy+wccij*distx*disty
          xz=xz+wccij*distx*distz
          yz=yz+wccij*disty*distz
          x=x+wccij*distx
          y=y+wccij*disty
          z=z+wccij*distz
          sij(icenter)=sij(icenter)+wccij
        end do
        dip(icenter,1)=x
        dip(icenter,2)=y
        dip(icenter,3)=z
        quadp(icenter,1,1)=xx-(yy+zz)/TWO
        quadp(icenter,1,2)=(THREE/TWO)*xy
        quadp(icenter,2,1)= quadp(icenter,1,2)
        quadp(icenter,1,3)=(THREE/TWO)*xz
        quadp(icenter,3,1)= quadp(icenter,1,3)
        quadp(icenter,2,2)=yy-(xx+zz)/TWO
        quadp(icenter,2,3)=(THREE/TWO)*yz
        quadp(icenter,3,2)= quadp(icenter,2,3)
        quadp(icenter,3,3)=zz-(yy+xx)/TWO
      end do
!$OMP END PARALLEL DO

      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: mp_pairs                                                    !!
!! purpose: adds one pair function's multipole-multipole terms, weighted  !!
!!   by ffact, to the six accumulators of every atom pair.                 !!
!! arguments:                                                              !!
!!   ffact  (in)    -- weight of this pair function                        !!
!!   sij, dip, quadp (in) -- its per-atom moments (from mp_moments)        !!
!!   atdist (in)    -- interatomic distances (from mp_init)                !!
!!   exk    (inout) -- charge-charge, charge-dipole, dipole-dipole,        !!
!!                     charge-quadrupole, dipole-quadrupole and            !!
!!                     quadrupole-quadrupole accumulators                  !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine mp_pairs(ffact,sij,dip,quadp,atdist,exk)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 muamub,muar,mubr,muaqbr,mubqar
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      dimension :: sij(nat),dip(nat,3),quadp(nat,3,3)
      dimension :: atdist(nat,nat),exk(nat,nat,6)
      dimension :: rvect(3)

!! pairwise electrostatic terms, iat==A, jat==B. parallel over iat: each   !!
!! iteration only writes exk(iat,jat>iat,:), a disjoint row per iat.      !!
!! rvect is PRIVATE (plain fixed-size local) so each thread gets its own  !!
!! scratch copy.                                                          !!
!$OMP PARALLEL DO PRIVATE(iat,jat,xx,ii,jj,kk,rvect,muamub,muar,mubr,
!$OMP&  rqar,rqbr,muaqbr,mubqar,xm,qaqb,rqaqbr)
      do iat=1,nat
        do jat=iat+1,nat
          xx=atdist(iat,jat)
          do ii=1,3
            rvect(ii)=coord(ii,jat)-coord(ii,iat)
          end do

!! charge-charge !!
          exk(iat,jat,1)=exk(iat,jat,1)+ffact*sij(iat)*sij(jat)/xx

!! charge-dipole !!
          muamub=ZERO
          muar=ZERO
          mubr=ZERO
          do ii=1,3
            muamub=muamub+dip(iat,ii)*dip(jat,ii)
            muar=muar+dip(iat,ii)*rvect(ii)
            mubr=mubr+dip(jat,ii)*rvect(ii)
          end do
          exk(iat,jat,2)=exk(iat,jat,2)+ffact*(muar*sij(jat)-mubr*sij(iat))/(xx**THREE)

!! dipole-dipole !!
          exk(iat,jat,3)=exk(iat,jat,3)-ffact*(THREE*muar*mubr/(xx**FIVE)-muamub/(xx**THREE))

!! charge-quadrupole !!
          rqar=ZERO
          rqbr=ZERO
          do ii=1,3
            do jj=1,3
              rqar=rqar+rvect(ii)*quadp(iat,ii,jj)*rvect(jj)
              rqbr=rqbr+rvect(ii)*quadp(jat,ii,jj)*rvect(jj)
            end do
          end do
          exk(iat,jat,4)=exk(iat,jat,4)+ffact*(rqar*sij(jat)+rqbr*sij(iat))/(xx**FIVE)

!! dipole-quadrupole !!
          muaqbr=ZERO
          mubqar=ZERO
          do ii=1,3
            do jj=1,3
              muaqbr=muaqbr+rvect(ii)*quadp(jat,ii,jj)*dip(iat,jj)
              mubqar=mubqar+rvect(ii)*quadp(iat,ii,jj)*dip(jat,jj)
            end do
          end do
          xm=-FIVE*(mubr*rqar-muar*rqbr)+TWO*xx*xx*(mubqar-muaqbr)
          exk(iat,jat,5)=exk(iat,jat,5)+ffact*xm/(xx**7.0d0)

!! quadrupole-quadrupole !!
          qaqb=ZERO
          rqaqbr=ZERO
          do ii=1,3
            do jj=1,3
              qaqb=qaqb+quadp(jat,ii,jj)*quadp(iat,ii,jj)
              do kk=1,3
                rqaqbr=rqaqbr+rvect(ii)*quadp(jat,ii,kk)*quadp(iat,kk,jj)*rvect(jj)
              end do
            end do
          end do
          xm=(35.0d0/THREE)*(rqar*rqbr)/(xx**9.0d0)+(TWO/THREE)*qaqb/(xx**FIVE)-(60.0/9.0d0)*rqaqbr/(xx**7.0d0)
          exk(iat,jat,6)=exk(iat,jat,6)+ffact*xm
        end do
      end do
!$OMP END PARALLEL DO

      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: mp_sum                                                      !!
!! purpose: sums the six multipole orders into the symmetric Excmp.        !!
!! arguments:                                                              !!
!!   exk   (in)  -- the six accumulators (from mp_pairs)                   !!
!!   Excmp (out) -- multipolar estimate per atom pair (i/=j elements)      !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine mp_sum(exk,Excmp)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension :: exk(nat,nat,6),Excmp(maxat,maxat)

      do i=1,nat
        do j=i+1,nat
          Excmp(i,j)=exk(i,j,1)+exk(i,j,2)+exk(i,j,3)+exk(i,j,4)+exk(i,j,5)+exk(i,j,6)
          Excmp(j,i)=Excmp(i,j)
        end do
      end do

      end

!! ***** !!

!! *********************************************************************** !!
!! subroutine: calc_coul                                                   !!
!! purpose: Coulomb (electron-electron) energy term for the two-electron   !!
!!   IQA/ENPART part, shared by numint_two and numint_two_uhf. Same-center !!
!!   and atom-pair kernels, phase-rotated grid throughout (all pairs are   !!
!!   computed here, unlike the exchange kernels -- no THREBOD skip).       !!
!! arguments:                                                              !!
!!   itotps (in)  -- total number of grid points                          !!
!!   wp/wppha (in) -- integration weight of each point, original/rotated   !!
!!   omp2/omp2pha (in) -- becke/tfvc weight of each point for every atom,  !!
!!                   original/rotated                                     !!
!!   pcoord/pcoordpha (in) -- xyz coordinates of each point, original/     !!
!!                   rotated                                               !!
!!   rho/rhopha (in) -- electron density at each point, original/rotated   !!
!!   coul   (out) -- Coulomb energy matrix, atomic/diatomic contributions  !!
!! author: PSalse, MGimf.                                                  !!
!! *********************************************************************** !!
      subroutine calc_coul(itotps,wp,wppha,omp2,omp2pha,pcoord,pcoordpha,rho,rhopha,coul)

      use integration_grid
      use OMP_LIB

      IMPLICIT REAL*8(A-H,O-Z)


      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot

      dimension :: wp(itotps),wppha(itotps),pcoord(itotps,3),pcoordpha(itotps,3)
      dimension :: omp2(itotps,nat),omp2pha(itotps,nat),rho(itotps),rhopha(itotps)
      dimension :: coul(maxat,maxat)

      character*100 threadenv

!! parallelization-related arrays below. !!
      allocatable :: ijpaircount(:,:),istart(:),iend(:) !! included atom pairs, thread-slicing bookkeeping !!
      allocatable :: ecoul_k(:,:)
      allocatable :: ecoul_ij(:,:)
      allocatable :: f3k(:)

      iatps = nrad*nang

!! same manual thread-chunking strategy as numint_two's exchange kernels; !!
!! determine the core count for thread splitting. !!
      call getenv('OMP_NUM_THREADS',threadenv)
      if(trim(threadenv)=='') then
        ithreadenv=ZERO
      else
        read(unit=threadenv,FMT='(I4)') ithreadenv
      end if
!! a parallel region only to query the thread and processor counts;  !!
!! every thread writes the same values to iprocs/ithreads            !!
      !$OMP PARALLEL
      iprocs=OMP_GET_MAX_THREADS()
      ithreads=INT(OMP_GET_NUM_PROCS())
      !$OMP END PARALLEL
      if(ithreadenv.ne.ZERO) ithreads=ithreadenv

!! per-thread grid-point slicing bookkeeping. !!
      ALLOCATE(f3k(ithreads))
      ALLOCATE(ecoul_k(nat,ithreads))
      ALLOCATE(istart(ithreads),iend(ithreads))
      itilerest=mod(iatps,ithreads)
      ispace=(iatps-itilerest)/ithreads

!! loop kept in this manual-chunking shape rather than a plain OMP loop --  !!
!! avoids an old cap that limited parallelization to 8 cores.               !!
      coul=ZERO
      ecoul_k=ZERO
      do icenter=1,nat
        ioffset=(icenter-1)*iatps
        istart=0
        iend=0
        do ik=1,(ithreads-1)
          istart(ik)=ioffset+((ik-1)*ispace)+1
          iend(ik)=ioffset+(ik*ispace)
        end do
        istart(ithreads)=ioffset+((ithreads-1)*ispace)+1
        iend(ithreads)=icenter*iatps

!! both calls below are required for the thread count to actually take effect. !!
        call omp_set_dynamic(.false.)
        call omp_set_num_threads(ithreads)
!! parallel over the thread slices ik, as in numint_two: f3loc is      !!
!! PRIVATE, ecoul_k(icenter,ik) written once per ifut.                  !!
!$OMP PARALLEL DO PRIVATE(ifut,jfut,x0,dx0,dy0,dz0,x1,dx1,dy1,dz1,dist)
        do ik=1,ithreads
          do ifut=istart(ik),iend(ik)
            x0=wp(ifut)*omp2(ifut,icenter)
            dx0=pcoord(ifut,1)
            dy0=pcoord(ifut,2)
            dz0=pcoord(ifut,3)
            f3k(ik)=ZERO
            do jfut=iatps*(icenter-1)+1,iatps*icenter
              x1=wppha(jfut)*omp2pha(jfut,icenter)
              dx1=pcoordpha(jfut,1)
              dy1=pcoordpha(jfut,2)
              dz1=pcoordpha(jfut,3)
              dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
              if(dist.gt.1.0d-12) f3k(ik)=f3k(ik)+rho(ifut)*rhopha(jfut)*x1*x0/dist
            end do
            ecoul_k(icenter,ik)=ecoul_k(icenter,ik)+f3k(ik)
          end do
        end do
!$OMP END PARALLEL DO
      end do

!! reduce the per-thread accumulator into the same-center coul term !!
!! (the factor of two is accounted for below, once coul is complete). !!
      do icenter=1,nat
        do isum=1,ithreads
          coul(icenter,icenter)=coul(icenter,icenter)+ecoul_k(icenter,isum)
        end do
      end do

!! atom-pair terms -- every pair is computed here, no THREBOD skip. !!
      ipaircounter=0
      ALLOCATE(ijpaircount(nat*nat,2))
      do icenter=1,nat
        do jcenter=icenter+1,nat
          ipaircounter=ipaircounter+1
          ijpaircount(ipaircounter,1)=icenter
          ijpaircount(ipaircounter,2)=jcenter
        end do
      end do

!! per-thread grid-point slicing bookkeeping for the atom-pair loop. !!
      ALLOCATE(ecoul_ij(ipaircounter,ithreads))
      itilerest=mod(iatps,ithreads)
      ispace=(iatps-itilerest)/ithreads
      ecoul_ij=ZERO
      do numpairnat=1,ipaircounter
        icenter=ijpaircount(numpairnat,1)
        jcenter=ijpaircount(numpairnat,2)
        ioffset=(icenter-1)*iatps
        istart=0
        iend=0
        do ik=1,(ithreads-1)
          istart(ik)=ioffset+((ik-1)*ispace)+1
          iend(ik)=ioffset+(ik*ispace)
        end do
        istart(ithreads)=ioffset+((ithreads-1)*ispace)+1
        iend(ithreads)=icenter*iatps
!! parallel over the thread slices ik, as the same-center block above. !!
!$OMP PARALLEL DO PRIVATE(ifut,jfut,x0,dx0,dy0,dz0,x1,dx1,dy1,dz1,dist)
        do ik=1,ithreads
          do ifut=istart(ik),iend(ik)
            x0=wp(ifut)*omp2(ifut,icenter)
            dx0=pcoord(ifut,1)
            dy0=pcoord(ifut,2)
            dz0=pcoord(ifut,3)
            f3k(ik)=ZERO
            do jfut=iatps*(jcenter-1)+1,iatps*jcenter
              x1=wp(jfut)*omp2(jfut,jcenter)
              dx1=pcoord(jfut,1)
              dy1=pcoord(jfut,2)
              dz1=pcoord(jfut,3)
              dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
              if(dist.gt.1.0d-12) f3k(ik)=f3k(ik)+rho(ifut)*rho(jfut)*x1*x0/dist
            end do
            ecoul_ij(numpairnat,ik)=ecoul_ij(numpairnat,ik)+f3k(ik)
          end do
        end do
!$OMP END PARALLEL DO
      end do

!! reduce the per-thread accumulator into the atom-pair coul term. !!
      do numpairnat=1,ipaircounter
        icenter=ijpaircount(numpairnat,1)
        jcenter=ijpaircount(numpairnat,2)
        do isum=1,ithreads
          coul(icenter,jcenter)=coul(icenter,jcenter)+ecoul_ij(numpairnat,isum)
        end do
      end do

!! same-center coul(i,i) was accumulated twice over (both (ifut,jfut) !!
!! orderings), halve it before summing into coulen. !!
      coulen=ZERO
      do i=1,nat
        coul(i,i)=coul(i,i)/TWO
        coulen=coulen+coul(i,i)
        do j=i+1,nat
          coul(j,i)=coul(i,j)
          coulen=coulen+coul(i,j)
        end do
      end do

      end

!! ***** !!

