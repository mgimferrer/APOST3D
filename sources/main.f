!! ----------------------------------------------------------------------
!! APOST-3D, version 5: real-space and Hilbert-space tools for wave
!! function analysis. P. Salvador, E. Ramos-Cordoba, M. Montilla,
!! L. Pujal and M. Gimferrer, J. Chem. Phys. 2024, 160, 172502.
!! The analyses and their papers are listed in the banner (kiir,
!! print.f) and on the citations page of the documentation.
!! Contact: pedro.salvador@udg.edu, mgimferrer18@gmail.com
!!
!! Written using parts of the program APOST by I. Mayer and A. Hamza
!! (Budapest, 2000-2003). The angular quadrature is the Lebedev-Laikov
!! code from CCL (V. I. Lebedev and D. N. Laikov, Doklady Mathematics
!! 59, 477, 1999); the exchange-correlation functionals come from libxc
!! (https://libxc.gitlab.io). We are grateful for these routines.
!!
!! The main program, in five steps:
!!   1) the job name, the .inp and .fchk files and the results file
!!   2) the keywords of the .inp (read_input, read_input.f)
!!   3) checks of keyword combinations
!!   4) iopt(200): the flags read by the analysis routines
!!   5) the atomic definition (Hilbert-space matrices or the real-space
!!      grid), then each requested analysis in turn: populations, local
!!      spin, EFOs/EOS/GEOS, LOBA, ENPART, OSLO, X-ray scattering
!! ----------------------------------------------------------------------

      use basis_set
      use ao_matrices
      use integration_grid
      use timing_mod
      use input_options_mod
      use json_mod, only: json_open, json_close, json_end,
     +  json_begin_object, json_begin_array, json_value, json_null
      use output_mod, only: print_value
      use omp_lib, only: omp_get_max_threads
      implicit real*8(a-h,o-z)
      include 'parameter.h'
!! general parameters !!
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
!! orbitals and density matrices; atom and fragment lists !!
      common /atlist/iatlist(maxat),icuat
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
!! populations and EOS !!
      common /loba/ oxi(maxat),errsav(maxat),elec(maxat),effpop(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
!! Enpart !!
      common /exchg/exch(maxat,maxat),xmix
!! TFVC features !!
      common /achi/achi(maxat,maxat),ibcp
      common /erf/aerf,ierf
!! MODGRID (diatXC) features !!
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22
!! NLOP features !!
      common /twoel/twoeltoler
      common /efield/field(4),edipole
!! printing and internal options !!
      common /filename/name0
      common /printout/iaccur
      common /iops/iopt(200)
!! for enpart !!
      dimension eto(maxat,maxat)
      character*60 name,name2,namepat,name3,name0,cgridarg
      character*40 cbuf
      character*20 cakey,caspace
      character*60 caname
      character*12 cashort

      allocatable wp(:),omp(:),omp2(:,:),chp(:,:),pcoord(:,:),rho(:)
      allocatable ibaspoint(:)
      allocatable sat(:,:,:)
      allocatable dm1(:,:),dm2(:,:,:,:)


      call cpu_time(time)
      call get_wall_time(wtime)

!! Processing arguments !!
      CALL GETARG(1,name0)
      if(name0.ne."") then
        j=len(name0)     
        do i=1,j
          if(name0(i:i).eq.' ') then
            l=i-1
            go to 10
          end if
        end do
10     name=name0(1:l)//".fchk"
        name2=name0(1:l)//".scr"
        name3=name0(1:l)//".inp"
        namepat=name0(1:l)
      else
        call apost_stop('The required input filename is missing')
      end if

!! print version info and the thread count of this run !!
      call kiir()
      write(*,'(a)') ''
      call print_value('OpenMP threads',omp_get_max_threads())

!! the grid is set in # GRID; a second argument stops with a message !!
      CALL GETARG(2,cgridarg)
      if(cgridarg.ne."") then
        write(*,'(2x,a)') 'The integration grid is no longer given on '//
     +    'the command line: use RADIAL and ANGULAR in # GRID'
        call apost_stop(' Grid on the command line')
      end if

!! the job name is given without extension; both files must exist (an   !!
!! open without status='old' would create them empty)                   !!
      if(index(name0(1:l),'.inp',back=.true.).eq.l-3.and.l.ge.4.or.
     +   index(name0(1:l),'.fchk',back=.true.).eq.l-4.and.l.ge.5) then
        write(*,'(2x,a,a,a)') 'Give the job name without extension: ',
     +    'apost3d <jobname> reads <jobname>.inp and <jobname>.fchk'
        call apost_stop(' Job name given with its extension')
      end if
      open(16,file=name3,status='old',iostat=ios)
      if(ios.ne.0) then
        write(*,'(2x,a,a,a)') 'Input file ',trim(name3),' not found'
        call apost_stop(' Input file not found')
      end if
      open(15,file=name,status='old',iostat=ios)
      if(ios.ne.0) then
        write(*,'(2x,a,a,a)') 'Wavefunction file ',trim(name),' not found'
        call apost_stop(' Wavefunction file not found')
      end if

!! the results file <job>.apost.json: started here, closed at every end !!
!! of the run (normal or apost_stop)                                     !!
      call json_open(namepat,name3,name)

      call read_input()
      call json_skeleton()

!! kop=1 -> unrestricted calculation; iposthf -> correlated             !!
!! calculation; icorr=0 -> no external dm1/2 provided; icas=1 ->        !!
!! CASSCF calculation; icisd -> CISD calculation; idono=1 -> do         !!
!! natural orbitals, idono=0 -> restricted SD calculation.              !!

!! Dependencies and cross-keyword validation !!
      iposthf=0
      idono=0
      if(icas.eq.1.or.icisd.eq.1) iposthf=1
      iopt(65)=iposthf
!! GEOS and EFFAO-U (ieffao=3, effao3d_u) need occ_no/c_no even for a  !!
!! restricted SD wavefunction -- there they're just the canonical MOs   !!
!! (integer occupied), giving a trivially ~0 unpaired channel, same as  !!
!! a plain restricted-wavefunction EFFAO run -- so force natural        !!
!! orbitals on for both regardless of kop.                             !!
      if(iposthf.eq.1.or.kop.eq.1.or.ieffao.eq.3) idono=1

      if(iposthf.eq.1) then
        if(ispin.eq.1.and.icorr.lt.2) call apost_stop(' Local Spin needs dm1 and dm2 for correlated WFs')
        if(ienpart.eq.1.and.icorr.lt.2) call apost_stop(' Enpart needs dm1 and dm2 for correlated WFs')
!! the external-field terms and POLAR are for HF and KS-DFT only !!
        if(ipolar.eq.1) then
          write(*,'(2x,a)') 'POLAR works with HF and KS-DFT wavefunctions only.'
          call apost_stop(' POLAR is not available for correlated wavefunctions')
        end if
        if(ienpart.eq.1.and.ifield.eq.1) then
          write(*,'(2x,a)') 'The .fchk file has an external electric field: ENPART'//
     +      ' of a correlated'
          write(*,'(2x,a)') 'wavefunction in a field is not available.'
          call apost_stop(' ENPART in a field is not available for correlated wavefunctions')
        end if
!! OSLO's iterative localize/assign/deflate algorithm assumes integer-  !!
!! occupied canonical-type orbitals -- it has no defined behavior for a !!
!! multireference wavefunction's fractionally-occupied natural orbitals !!
!! (CASSCF/CISD/FCI). A bare stop's message only reaches stderr, not    !!
!! the saved .apost output, so write it explicitly first.               !!
        if(ioslo.eq.1) then
          write(*,'(2x,a)')
     +      'OSLO cannot be performed for multireference (CASSCF/CISD/FCI) wavefunctions'
          call apost_stop(' OSLO cannot be performed for multireference wavefunctions')
        end if
      end if
!! the RDMs are sized by the CASSCF/CI data of the .fchk: with a single  !!
!! determinant they would be read as empty                              !!
      if(icorr.ne.0.and.iposthf.eq.0) then
        write(*,'(2x,a)') 'DM 1 / DM 2 read the RDMs of a correlated '//
     +    'wavefunction: the .fchk must be a CASSCF, FCI or CI one'
        call apost_stop(' DM needs a correlated-wavefunction .fchk')
      end if
      if(iqtaim.eq.1) call apost_stop('This version can not do QTAIM')

      if(imulli.gt.1.or.iqtaim.eq.1) iopop=0
      if(ihirsh.ne.0.and.idoat.eq.1) call apost_stop('Cant do HIRSH with DOATOMS')
!! gated on the real restricted-SD condition directly, not on idono --  !!
!! idono can now also be forced on by GEOS (see above) on a genuinely   !!
!! restricted wavefunction, where Local Spin still doesn't make sense   !!
!! (confirmed: produces nonsense u_A/N_D values, not a trivial zero).   !!
      if(ispin.eq.1.and.kop.eq.0.and.iposthf.eq.0)  then
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Note: SPIN skipped: a closed-shell single '//
     +    'determinant has no local spins'
        ispin=0
      end if
!! ENPART needs a real-space atomic definition !!
      if(imulli.ne.0.and.ienpart.ne.0)
     +  call apost_stop('ENPART needs a real-space AIM (e.g. TFVC)')
      if(iqtaim.ne.0.and.ienpart.ne.0) call apost_stop('Cant do ENPART with QTAIM')
      if(ieos.eq.1.and.idoat.eq.1) call apost_stop('Cant do EOS with DOATOMS')
!! GEOS/EFFAO-U (ieffao=3) exist only on the real-space fragment path   !!
!! (effao3d_u), not with Hilbert-space atoms or DOATOMS                  !!
      if(ieffao.eq.3.and.imulli.ne.0)
     +  call apost_stop('GEOS/EFFAO-U need a real-space AIM (e.g. TFVC)')
      if(ieffao.eq.3.and.idoat.eq.1) call apost_stop('Cant do GEOS/EFFAO-U with DOATOMS')
!! weighted Lowdin (LOWDIN-W) has no EFO construction of its own        !!
      if(ieffao.ne.0.and.imulli.eq.5) then
        write(*,'(2x,a)') 'EFFAO, UEFFAO and EOS are not available with '//
     +    'LOWDIN-W'
        call apost_stop(' EFFAO/EOS not available with LOWDIN-W')
      end if
!! OSLO integrates its spread matrix on the # METHOD grid; the atoms that !!
!! assign OSLOs to fragments can still be Hilbert-space, set in # OSLO    !!
      if(ioslo.eq.1.and.imulli.ne.0) then
        write(*,'(2x,a)') 'OSLO needs a real-space AIM in # METHOD '//
     +    '(e.g. TFVC), used for the numerical integration'
        write(*,'(2x,a)') 'MULLIKEN, LOWDIN, LOWDIN-DAVIDSON or '//
     +    'NAO-BASIS go in # OSLO, for the fragment populations'
        call apost_stop(' OSLO needs a real-space AIM (e.g. TFVC)')
      end if
!! a restricted wavefunction has a single OSLO channel !!
      if(ioslo.eq.1.and.kop.eq.0.and.nbranch(2)+nbranch(3).gt.0) then
        write(*,'(2x,a)') 'BRANCH_ALPHA and BRANCH_BETA need an '//
     +    'unrestricted wavefunction; use BRANCH_ITERATION'
        call apost_stop(' BRANCH_ALPHA/BRANCH_BETA need an unrestricted '//
     +    'wavefunction')
      end if

      if(kop.eq.0.and.ispin.eq.1) then
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Note: closed-shell correlated wavefunction: '//
     +    'its local spins come from electron correlation'
      end if

!! iopt(1..200), every slot listed: a commented line is a free slot,  !!
!! ready for a new option                                             !!
      iopt(1) = idoint
!      iopt(2) =
!      iopt(3) =
      iopt(4) = idono
      iopt(5) = imulli
      iopt(6) = ihirsh
      iopt(7) = iallpo
      iopt(8) = iopop
!      iopt(9) =   !! ndens0, set in read_input() !!
!      iopt(10) =

!      iopt(11) =
      iopt(12) = ieffao
      iopt(13) = icube
      iopt(14) = ibcp
      iopt(15) = ispin
      iopt(16) = iqtaim
      iopt(17) = ienpart
      iopt(18) = ihf
!      iopt(19) =
      iopt(20) = iexact

      iopt(21) = ihomo
      iopt(22) = idek
      iopt(23) = iionic
      iopt(24) = ieffthr
      iopt(25) = istiff
      iopt(26) = icorr
!      iopt(27) =
!      iopt(28) =
!      iopt(29) =
!      iopt(30) =

      iopt(31) = inewbec
!      iopt(32) =
!      iopt(33) =
      iopt(34) = ilaplacian
      iopt(35) = istep
      iopt(36) = inna
      iopt(37) = imaxdist
      iopt(38) = iscreening
      iopt(39) = ipath
      iopt(40) = idofr

      iopt(41) = jcubthr
      iopt(42) = kcubthr
      iopt(43) = iorca
      iopt(44) = ithrebod
!      iopt(45) =
      iopt(46) = ipolar
      iopt(47) = ifield
      iopt(48) = iradmat
      iopt(49) = ielcount
!      iopt(50) =

!      iopt(51) =
!      iopt(52) =
!      iopt(53) =
!      iopt(54) =
!      iopt(55) =   !! functional family, set in the ENPART dispatch !!
      iopt(56) = idftdm1
      iopt(57) = id_func_dm1 !! dft_dm1.f's own functional-ID read (ifunc=Iopt(57)) !!
!      iopt(58) =
      iopt(59) = iecorr
      iopt(60) = id_xcfunc

      iopt(61) = id_xfunc
      iopt(62) = id_cfunc
      iopt(63) = inatorb_dm1 !! dft_dm1.f's own NATORB read !!
!      iopt(64) =
!      iopt(65) =   !! iposthf, set above with the keyword checks !!
!      iopt(66) =
!      iopt(67) =
!      iopt(68) =
!      iopt(69) =
!      iopt(70) =

!      iopt(71) =
!      iopt(72) =
!      iopt(73) =
!      iopt(74) =
!      iopt(75) =
!      iopt(76) =
!      iopt(77) =
!      iopt(78) =
      iopt(79) = imokit !! MG: Temporarily !!
!      iopt(80) =

!      iopt(81) =
!      iopt(82) =
!      iopt(83) =
!      iopt(84) =
!      iopt(85) =
!      iopt(86) =
      iopt(87) = ipyscf
      iopt(88) = inopop
!      iopt(89) =
      iopt(90) = ieoscent

      iopt(91) = ianalytical
!      iopt(92) =
!      iopt(93) =
!      iopt(94) =
      iopt(95) = iqchem
      iopt(96) = ifolitol
!      iopt(97) =
      iopt(98) = ioslofchk
      iopt(99) = iloba
!      iopt(100) =

!      iopt(101) =
!      iopt(102) =
!      iopt(103) =
!      iopt(104) =
!      iopt(105) =
!      iopt(106) =
!      iopt(107) =
!      iopt(108) =
!      iopt(109) =
!      iopt(110) =

!      iopt(111) =
!      iopt(112) =
!      iopt(113) =
!      iopt(114) =
!      iopt(115) =
!      iopt(116) =
!      iopt(117) =
!      iopt(118) =
!      iopt(119) =
!      iopt(120) =

!      iopt(121) =
!      iopt(122) =
!      iopt(123) =
!      iopt(124) =
!      iopt(125) =
!      iopt(126) =
!      iopt(127) =
!      iopt(128) =
!      iopt(129) =
!      iopt(130) =

!      iopt(131) =
!      iopt(132) =
!      iopt(133) =
!      iopt(134) =
!      iopt(135) =
!      iopt(136) =
!      iopt(137) =
!      iopt(138) =
!      iopt(139) =
!      iopt(140) =

!      iopt(141) =
!      iopt(142) =
!      iopt(143) =
!      iopt(144) =
!      iopt(145) =
!      iopt(146) =
!      iopt(147) =
!      iopt(148) =
!      iopt(149) =
!      iopt(150) =

!      iopt(151) =
!      iopt(152) =
!      iopt(153) =
!      iopt(154) =
!      iopt(155) =
!      iopt(156) =
!      iopt(157) =
!      iopt(158) =
!      iopt(159) =
!      iopt(160) =

!      iopt(161) =
!      iopt(162) =
!      iopt(163) =
!      iopt(164) =
!      iopt(165) =
!      iopt(166) =
!      iopt(167) =
!      iopt(168) =
!      iopt(169) =
!      iopt(170) =

!      iopt(171) =
!      iopt(172) =
!      iopt(173) =
!      iopt(174) =
!      iopt(175) =
!      iopt(176) =
!      iopt(177) =
!      iopt(178) =
!      iopt(179) =
!      iopt(180) =

!      iopt(181) =
!      iopt(182) =
!      iopt(183) =
!      iopt(184) =
!      iopt(185) =
!      iopt(186) =
!      iopt(187) =
!      iopt(188) =
!      iopt(189) =
!      iopt(190) =

!      iopt(191) =
!      iopt(192) =
!      iopt(193) =
!      iopt(194) =
!      iopt(195) =
!      iopt(196) =
!      iopt(197) =
!      iopt(198) =
!      iopt(199) =
!      iopt(200) =
!! End of iop listing !!

!! input summary: every flag is final here !!
      call print_input_summary()

!! Natural orbitals from P-matrix in .fchk !!
      if (idono.eq.1) call gennatural()

      write(*,'(a)') ''
      call cpu_time(time2)
      call get_wall_time(wtime2)
      call print_timer('Initialization',time2-time,wtime2-wtime)
      time=time2
      wtime=wtime2

!! density at the iatdens atom !!
      if(iatdens.ne.0) call atdens_int(Rmax,iatdens)
!! NOPOPU ends the run here (normal end, exit 0); flush so the stop !!
!! message (stderr) comes after the buffered output                 !!
      if(inopop.eq.1) then
        call json_close('normal','')
        flush(6)
        stop 'Normal termination of APOST3D'
      end if

!! set atomic radii just in case needed !!
      call prepar()

!! Hilbert-space atomic overlap matrices (Mulliken, Lowdin, NAO, ...)   !!
      ndim=igr
      if(imulli.ge.1) then
        allocate (sat(ndim,ndim,nat))
        call print_box('SETTING ATOMIC DEFINITION')
        if(imulli.eq.1) then
          call print_value('Hilbert-space atoms','Mulliken')
          call tomull(sat)
        else if(imulli.eq.2.or.imulli.eq.3) then
          if(imulli.eq.2) call print_value('Hilbert-space atoms',
     +      'Lowdin basis')
          if(imulli.eq.3) call print_value('Hilbert-space atoms',
     +      'Lowdin-Davidson basis')
          call tolow(sat)
        else if(imulli.eq.4) then
          call print_value('Hilbert-space atoms','NAO basis')
          call tonao(sat)
        else if(imulli.eq.5) then
          call print_value('Hilbert-space atoms','weighted Lowdin basis')
          call tolow2(sat)
        end if

      else

!! Prepare for numerical integrations !!
        call build_integration_grid(ienpart, ipolar)

!! kept to rebuild this grid after ENPART/DFT-DM1 replace it !!
        nradp=nrad
        nangp=nang
        rr00p=rr00

        iatps=nang*nrad
        itotps=nat*iatps
        allocate (wp(itotps),omp(itotps),omp2(itotps,nat))
        allocate (chp(itotps,ndim),rho(itotps))
        allocate (pcoord(itotps,3),ibaspoint(itotps))

! init populations to zero
        do i=1,nat
          qat(i,1)=ZERO
        end do

!! prepare for Becke-type atomic weighting !!
!! the stiffness of the cell functions does not apply to Hirshfeld !!
        call print_box('SETTING ATOMIC DEFINITION')
        if(ihirsh.eq.0) then
          if(ierf.ne.1) then
            call print_value('Cell stiffness k',istiff)
          else
            write(cbuf,'(f8.4)') aerf
            call print_value('Adjustable cell profile A',
     +        trim(adjustl(cbuf)))
          end if
        end if
        if(ibcp.eq.1) call khi()
        if(ihirsh.eq.1.or.ihirsh.eq.2.or.ielcount.eq.1) call makeatdens
 
        iiter=1
        call prenumint(ndim,itotps,nat,wp,omp,omp2,chp,rho,pcoord,ibaspoint,iiter)

!! OS from localized MOs (centroids) !!
        if(ieoscent.eq.1) then
          write(*,'(2x,a)') 'Oxidation states from the centroids of '//
     +      'localized orbitals (OS-CENTROID)'
          call eos_centroid(itotps,chp,wp,omp,pcoord)
!! normal end of this path: flush so stdout keeps its order, exit 0 !!
          call json_close('normal','')
          flush(6)
          stop
        end if

!! Integrate atomic overlap by default !!
        allocate (sat(ndim,ndim,nat))
        call numint_sat(ndim,itotps,nat,wp,omp,omp2,chp,ibaspoint,sat)

      end if
!! end of preparation for numerical integrations !!

!! atomic-domain partition time, covers both branches above (Hilbert    !!
!! and real-space).                                                     !!
      write(*,'(a)') ''
      call cpu_time(time2)
      call get_wall_time(wtime2)
      call print_timer('Atomic definition',time2-time,wtime2-wtime)
      time=time2
      wtime=wtime2

!! Correlated-WF input, needs dm1 and/or dm2 files !!
!! cas/cisd specifications !!
      nelec=nalf+nb
      if(icorr.ne.0) then
        call print_box('POST-HARTREE-FOCK CALCULATION')

        if(icisd.eq.1) nspinorb=nbasis*2
!! a CASSCF 1-RDM covers the active space: dm1input reports the counts !!
        if(icas.ne.1) call print_value('Spin-orbitals',nspinorb)

        ALLOCATE(dm1(nspinorb,nspinorb))
        call dm1input(dm1)
        if(icorr.eq.2) then
          norb=nspinorb/2
          ALLOCATE(dm2(norb,norb,norb,norb))
          if(iorca.eq.1.or.ipyscf.eq.1) then
            call dm2input_pyscf(dm1,dm2)
          else
            call dm2input_dmn(dm1,dm2)
          end if
        end if
      end if
!! end of correlated-WF input !!

!! Population analysis !!
      call population_density(sat)

      if(idoint.eq.1)  call print_int(ndim,nat,sat,namepat)

      call population_print_charges()

      if(ielcount.eq.1) then
        iatps=nang*nrad
        CALL NCTAIM(iatps,wp,omp,omp2,nat,sat,igr,ibaspoint,chp)
      end if

      call population_print_overlap(wp,omp,omp2,rho)
      call json_populations()

!! Bond orders and valences !!
      call bond_order_analysis(sat)

      write(*,'(a)') ''
      call cpu_time(time2)
      call get_wall_time(wtime2)
      call print_timer('Population analysis',time2-time,wtime2-wtime)
      time=time2
      wtime=wtime2

!! Local spin decomposition, single-determinant WF !!
      if(ispin.eq.1.and.icas.eq.0.and.icisd.eq.0) then
        call print_box('DOING LOCAL SPIN ANALYSIS')
        call print_value('Wavefunction','single determinant')
        call fspindec(sat)
      end if

!! Local spin and DIs for correlated WFs, needs dm1/dm2 from the DMN   !!
!! code (E. Matito) or from pySCF                                      !!
      if((icas.eq.1.or.icisd.eq.1).and.ispin.eq.1)then
        call print_box('DOING LOCAL SPIN ANALYSIS')
        call print_value('Wavefunction','correlated (1- and 2-RDM)')
        call spincorr(sat,dm1,dm2)
        write(*,'(a)') ''
        call cpu_time(time2)
        call get_wall_time(wtime2)
        call print_timer('Local spin analysis',time2-time,wtime2-wtime)
        time=time2
        wtime=wtime2
      end if

!! Nonlinear optical properties (POLAR) !!
      if(ipolar.ne.0) call polar(itotps,nat,wp,omp,omp2,pcoord,rho)

!! EFFAO part                                                          !!
!! ieffao: 0 nothing, 1 eff-AOs, 2 spin-resolved eff-AOs (a must for   !!
!!   EOS), 3 paired/unpaired eff-AOs                                   !!
!! imulli: 1 Mulliken, 2 Lowdin, 3 Lowdin-Davidson, 4 NAO, 5 weighted   !!
!!   Lowdin (no EFOs, stopped at the input checks)                       !!
      if (ieffao.ne.0) then

!! accumulate eos_analysis time separately (see print_timer calls below). !!
!! not wired up for ieffao.eq.3 (GEOS) -- effao3d_u times as EFFAO only. !!
        xeos_cpu=ZERO
        xeos_wall=ZERO

        if(idofr.eq.0) then
          icufr=nat
          do i=1,icufr
            nfrlist(i)=1
            ifrlist(1,i)=i                
          end do
        end if

!! beta EFOs computed separately (open shell, or a correlated          !!
!! wavefunction with a spin density; the test differs between the      !!
!! Hilbert-space and real-space paths)                                 !!
        idobeta=0
        if(imulli.ge.1) then
          if(kop.ne.0.or.(icas.eq.1.and.nalf.ne.nb)) idobeta=1
        else
          if(kop.ne.0.or.(icas.eq.1.and.icass.ne.0)) idobeta=1
        end if
        if(ieffao.ne.3) call efo_section_header(idobeta,
     +    dble(ieffthr)/1000.0d0)
        if(ieffao.ne.3.and.idoat.eq.0) call efo_json_open()

!! Mulliken !!
        if(imulli.eq.1) then
          if(ieffao.eq.1) then
            call ueffaomull_frag(0)
          else if (ieffao.eq.2) then
            call ueffaomull_frag(1)
            if (ieos.eq.1) then
              call cpu_time(xeos1)
              call get_wall_time(wxeos1)
              call eos_analysis(0,1,xthresh)
              call cpu_time(xeos2)
              call get_wall_time(wxeos2)
              xeos_cpu=xeos_cpu+(xeos2-xeos1)
              xeos_wall=xeos_wall+(wxeos2-wxeos1)
            end if
            if(idobeta.eq.1) call ueffaomull_frag(2)
            if (ieos.eq.1) then
              call cpu_time(xeos1)
              call get_wall_time(wxeos1)
              call eos_analysis(idobeta,2,xthresh)
              call cpu_time(xeos2)
              call get_wall_time(wxeos2)
              xeos_cpu=xeos_cpu+(xeos2-xeos1)
              xeos_wall=xeos_wall+(wxeos2-wxeos1)
            end if
          end if

!! Lowdin !!
        else if(imulli.gt.1) then
          if(ieffao.eq.1) then
            call ueffaolow_frag(0)
          else if (ieffao.eq.2) then
            call ueffaolow_frag(1)
            if(ieos.eq.1) then
              call cpu_time(xeos1)
              call get_wall_time(wxeos1)
              call eos_analysis(0,1,xthresh)
              call cpu_time(xeos2)
              call get_wall_time(wxeos2)
              xeos_cpu=xeos_cpu+(xeos2-xeos1)
              xeos_wall=xeos_wall+(wxeos2-wxeos1)
            end if
            if(idobeta.eq.1) call ueffaolow_frag(2)
            if(ieos.eq.1) then
              call cpu_time(xeos1)
              call get_wall_time(wxeos1)
              call eos_analysis(idobeta,2,xthresh)
              call cpu_time(xeos2)
              call get_wall_time(wxeos2)
              xeos_cpu=xeos_cpu+(xeos2-xeos1)
              xeos_wall=xeos_wall+(wxeos2-wxeos1)
            end if
          end if

!! 3D-space !!
        else

!! EFFAO/UEFFAO for selected atoms only (no EOS) !!
          if(idoat.ne.0) then
            if(ieffao.eq.1) then
              call ueffao3d(itotps,ndim,omp,chp,sat,wp,omp2,p,0)
            else
              call ueffao3d(itotps,ndim,omp,chp,sat,wp,omp2,pa,1)
              call ueffao3d(itotps,ndim,omp,chp,sat,wp,omp2,pb,2)
            end if

          else

!! EFFAO/UEFFAO/EFFAO-U for fragments/all atoms !!
            if(ieffao.eq.1) then
              call ueffao3d_frag(itotps,ndim,omp,chp,sat,wp,omp2,p,0)
            else if(ieffao.eq.2) then
              call ueffao3d_frag(itotps,ndim,omp,chp,sat,wp,omp2,pa,1)
              if(ieos.eq.1) then
                call cpu_time(xeos1)
                call get_wall_time(wxeos1)
                call eos_analysis(idobeta,1,xthresh)
                call cpu_time(xeos2)
                call get_wall_time(wxeos2)
                xeos_cpu=xeos_cpu+(xeos2-xeos1)
                xeos_wall=xeos_wall+(wxeos2-wxeos1)
              end if
              if(idobeta.eq.1) then
                call ueffao3d_frag(itotps,ndim,omp,chp,sat,wp,omp2,pb,2)
              end if
              if(ieos.eq.1) then
                call cpu_time(xeos1)
                call get_wall_time(wxeos1)
                call eos_analysis(idobeta,2,xthresh)
                call cpu_time(xeos2)
                call get_wall_time(wxeos2)
                xeos_cpu=xeos_cpu+(xeos2-xeos1)
                xeos_wall=xeos_wall+(wxeos2-wxeos1)
              end if

!! GEOS part: oxidation-state analysis called from inside the routine !!
            else if(ieffao.eq.3) then
              call effao3d_u(itotps,ndim,omp,chp,sat,wp,omp2,iueos) 
            end if
          end if
        end if
        if(ieffao.ne.3.and.idoat.eq.0) call json_end()

        write(*,'(a)') ''
        call cpu_time(time2)
        call get_wall_time(wtime2)
        call print_timer('EFFAO computation',(time2-time)-xeos_cpu,(wtime2-wtime)-xeos_wall)
        call print_timer('EOS analysis',xeos_cpu,xeos_wall)
        time=time2
        wtime=wtime2
      end if

!! Localized orbital bonding analysis (LOBA) !!
      if(iloba.eq.1) then
        call print_box('DOING LOCALIZED ORBITAL BONDING ANALYSIS (LOBA)')

!! Hilbert-space !!
        if(imulli.gt.0) call apost_stop('LOBA needs a real-space '//
     +    'atomic definition')

!! Real-space !!
        if(imulli.eq.0) call eos_loba(sat)
      end if

!! Energy decomposition (ENPART) !!
      if(ienpart.eq.1) then
        call print_box('ENERGY PARTITIONING (ENPART)')
        if(icas.eq.1) then
          cbuf='CASSCF'
        else if(icisd.eq.1) then
          cbuf='CISD'
        else if(ihf.eq.1) then
          cbuf='Hartree-Fock'
        else
          cbuf='KS-DFT'
        end if
        if(kop.eq.1) then
          cbuf=trim(cbuf)//', unrestricted'
        else
          cbuf=trim(cbuf)//', restricted'
        end if
        call print_value('Wavefunction',trim(cbuf))
        call atomic_definition_names(cakey,caname,cashort,caspace)
        call print_value('Atomic definition',trim(cakey)//' (real space)')
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Each matrix below holds the one-center (self) terms'//
     +    ' on the diagonal and'
        write(*,'(2x,a)') 'the whole A-B pair term off it; its total adds the'//
     +    ' diagonal and one'
        write(*,'(2x,a)') 'triangle. Energies in au.'
        if(iverbose.eq.0) write(*,'(2x,a)')
     +    'VERBOSE 1 in # METHOD also prints the intermediate matrices.'

!! results file: the enpart block (closed after the two-electron    !!
!! part); the external-field terms stay out of it                    !!
        call json_begin_object('enpart')
        if(icas.eq.1) then
          call json_value('wavefunction','CASSCF')
        else if(icisd.eq.1) then
          call json_value('wavefunction','CISD')
        else if(ihf.eq.1) then
          call json_value('wavefunction','Hartree-Fock')
        else
          call json_value('wavefunction','KS-DFT')
        end if
        if(kop.eq.1) then
          call json_value('formalism','unrestricted')
        else
          call json_value('formalism','restricted')
        end if

!! One-electron part: CASSCF and CI WFs !!
        if(iposthf.eq.1) then
          call json_null('functional')
          call json_enpart_pairs()
          call numint_one_rphf(ndim,itotps,wp,omp2,pcoord,chp,rho,eto)
          call cpu_time(time2)
          call get_wall_time(wtime2)
          write(*,'(a)') ''
          call print_timer('ENPART one-electron terms',time2-time,
     +  wtime2-wtime)
          time=time2
          wtime=wtime2

!! One-electron part: DFT and HF WFs !!
        else 

!! Initialize DFT functional for info and initial printing !!
          if(id_xfunc.eq.-1) call json_null('functional')
          if(id_xfunc.ne.-1) then
            call json_begin_object('functional')
            if(len_trim(cfuncname).gt.0) then
              call json_value('keyword',trim(cfuncname))
            else
              call json_value('keyword','LIBRARY')
            end if
            call json_begin_array('components')
            ifuncfirst=1
            itype=0
            jtype=0
            if(id_xcfunc.ne.0) then
              call func_info_print(id_xcfunc,itype,ifuncfirst)
              ifuncfirst=0
            end if
            if(id_cfunc.ne.0) then
              call func_info_print(id_cfunc,itype,ifuncfirst)
              ifuncfirst=0
            end if
            if(id_xfunc.ne.0) then
              call func_info_print(id_xfunc,jtype,ifuncfirst)
              ifuncfirst=0
            end if
            if(itype.ge.jtype) iopt(55) = itype
            if(jtype.gt.itype) iopt(55) = jtype
            call json_end()
            call json_value('exact_exchange_fraction',xmix)
            call json_end()
          end if
          call json_enpart_pairs()

!! Restricted case !!
          if(kop.ne.1) then

!! One-electron terms !!
            call numint_one(ndim,itotps,wp,rho,omp,omp2,pcoord,chp,eto)
            call cpu_time(time2)
            call get_wall_time(wtime2)
            write(*,'(a)') ''
            call print_timer('ENPART one-electron terms',time2-time,
     +  wtime2-wtime)
            time=time2
            wtime=wtime2

!! DFT XC term !!
            if(id_xfunc.ne.-1) then
              if(ianalytical.eq.0) then
                call numint_dft(ndim,itotps,wp,rho,omp,omp2,chp,eto,pcoord,sat)
              else
                call numint_dft_analytical(ndim,itotps,wp,omp,omp2,chp,eto,pcoord,sat)
              end if
              call cpu_time(time2)
              call get_wall_time(wtime2)
              write(*,'(a)') ''
              call print_timer('ENPART DFT exchange-correlation',
     +  time2-time,wtime2-wtime)
              time=time2
              wtime=wtime2
            end if

!! Unrestricted case !!
          else

!! One-electron terms !!
            call numint_one_uhf(ndim,itotps,wp,rho,omp,omp2,pcoord,chp,eto)
            call cpu_time(time2)
            call get_wall_time(wtime2)
            write(*,'(a)') ''
            call print_timer('ENPART one-electron terms',time2-time,
     +  wtime2-wtime)
            time=time2
            wtime=wtime2

!! DFT XC term !!
            if(id_xfunc.ne.-1) then
              call numint_dft_uks(ndim,itotps,wp,omp,omp2,chp,pcoord,eto)
              call cpu_time(time2)
              call get_wall_time(wtime2)
              write(*,'(a)') ''
              call print_timer('ENPART DFT exchange-correlation',
     +  time2-time,wtime2-wtime)
              time=time2
              wtime=wtime2
            end if
          end if

!! End of restricted/unrestricted dispatch !!
        end if
        DEALLOCATE(wp,omp,omp2,chp,pcoord,ibaspoint,rho)

!! Two-electron part !!

!! Two-electron integration defaults !!
        call print_box('TWO-ELECTRON GRID')

!! RADIAL_2E/ANGULAR_2E of # GRID (modgrid common, read_grid); default !!
!! 150/590; the rotation angles follow the angular grid                  !!
        nrad=nrad22
        nang=nang22
        rr00=rr0022

!! ANALYTIC: the 70/434 two-electron grid, unless # GRID sets one     !!
!! (igrid2e=1)                                                        !!
        if(ianalytical.eq.1.and.igrid2e.eq.0) then
          write(*,'(2x,a)') 'ANALYTIC: two-electron grid raised to 70 x 434'//
     +      ' (the two-electron part is a one-electron integral)'
          nrad=70
          nang=434
        else if(ianalytical.eq.1) then
          write(*,'(2x,a)') 'ANALYTIC: two-electron grid of # GRID kept'
        end if

!! Printing info !!
        call print_value('Radial points',nrad)
        call print_value('Angular points',nang)
        call print_value('Grid points',nrad*nang*nat)

!! Generating grid for two-electron numerical integrations !!
        iatps=nang*nrad
        itotps=nrad*nang*nat
        call quad(Nrad,Nang) 
        ALLOCATE(wp(itotps),omp(itotps),omp2(itotps,nat),rho(itotps))
        ALLOCATE(pcoord(itotps,3),ibaspoint(itotps),chp(itotps,ndim))
        call prenumint(ndim,itotps,nat,wp,omp,omp2,chp,rho,pcoord,ibaspoint,0)

!! CASSCF and CI WFs !!
        if(iposthf.eq.1) then
          call numint_two_rphf(ndim,itotps,wp,omp2,pcoord,chp,rho,eto,dm1,dm2)
          call cpu_time(time2)
          call get_wall_time(wtime2)
          write(*,'(a)') ''
          call print_timer('ENPART two-electron terms',time2-time,
     +  wtime2-wtime)
          time=time2
          wtime=wtime2
!! closes the enpart block of the results file !!
          call json_end()

!! DFT and HF WFs !!
        else

!! Unrestricted case !!
          if(kop.eq.1) then 
            call numint_two_uhf(ndim,itotps,wp,omp,omp2,pcoord,chp,rho,eto)

!! Restricted case !!
          else
            call numint_two(ndim,itotps,wp,omp,omp2,pcoord,chp,rho,eto)
          end if
          write(*,'(a)') ''
          call cpu_time(time2)
          call get_wall_time(wtime2)
          call print_timer('ENPART two-electron terms',time2-time,
     +  wtime2-wtime)
          time=time2
          wtime=wtime2
!! closes the enpart block of the results file !!
          call json_end()

!! End of WF-type dispatch !!
        end if 

!! Deallocating any existent grid !!
        DEALLOCATE(wp,omp,omp2,chp,pcoord,ibaspoint,rho)

      end if
!! end of energy decomposition !!

!! DFT-DM1 approximate one-particle RDM1.                                !!
!! Runs standalone, own dedicated grid reusing the two-electron-type      !!
!! RADIAL_2E/ANGULAR_2E of # GRID (read_input.f's read_grid); dft_dm1     !!
!! builds its own rotated second grid internally.                        !!
      if(idftdm1.eq.1) then
        call print_box('DOING DFT-DM1 APPROXIMATE ONE-PARTICLE RDM1')
        ndim=igr
        nrad=nrad22
        nang=nang22
        rr00=rr0022
        iatps=nang*nrad
        itotps=nrad*nang*nat
        pha=ZERO
        phb=ZERO
        call quad(nrad,nang)

!! the primary grid (TFVC/Mulliken/etc.) may still be allocated here --  !!
!! ENPART's own two teardown points only fire when ienpart=1, and        !!
!! DFT-DM1 needs its own grid at nrad22/nang22 regardless                !!
        if(allocated(wp))        DEALLOCATE(wp)
        if(allocated(omp))       DEALLOCATE(omp)
        if(allocated(omp2))      DEALLOCATE(omp2)
        if(allocated(pcoord))    DEALLOCATE(pcoord)
        if(allocated(ibaspoint)) DEALLOCATE(ibaspoint)
        if(allocated(chp))       DEALLOCATE(chp)
        if(allocated(rho))       DEALLOCATE(rho)

        ALLOCATE(wp(itotps),omp(itotps),omp2(itotps,nat),rho(itotps))
        ALLOCATE(pcoord(itotps,3),ibaspoint(itotps),chp(itotps,igr))
        call prenumint(ndim,itotps,nat,wp,omp,omp2,chp,rho,pcoord,ibaspoint,0)
        DEALLOCATE(ibaspoint,omp,rho)
        call dft_dm1(itotps,wp,omp2,pcoord,chp)
        DEALLOCATE(wp,omp2,chp,pcoord)
      end if

!! ENPART and DFT-DM1 replace the primary grid with their own and free  !!
!! it: rebuild it, with its original settings, for OSLO and SCATT-FACT  !!
      if((ioslo.eq.1.or.iscattfact.eq.1).and.imulli.eq.0.and.
     +   .not.allocated(wp)) then
        nrad=nradp
        nang=nangp
        rr00=rr00p
        pha=ZERO
        phb=ZERO
        call quad(nrad,nang)
        ndim=igr
        iatps=nang*nrad
        itotps=nat*iatps
        if(allocated(omp))       DEALLOCATE(omp)
        if(allocated(omp2))      DEALLOCATE(omp2)
        if(allocated(chp))       DEALLOCATE(chp)
        if(allocated(rho))       DEALLOCATE(rho)
        if(allocated(pcoord))    DEALLOCATE(pcoord)
        if(allocated(ibaspoint)) DEALLOCATE(ibaspoint)
        ALLOCATE(wp(itotps),omp(itotps),omp2(itotps,nat))
        ALLOCATE(chp(itotps,ndim),rho(itotps))
        ALLOCATE(pcoord(itotps,3),ibaspoint(itotps))
        call prenumint(ndim,itotps,nat,wp,omp,omp2,chp,rho,pcoord,
     +    ibaspoint,1)
      end if

!! OSLO -- variants of the procedure can be found in the dev version !!
      if(ioslo.eq.1) then
        call print_box('OXIDATION STATES FROM LOCALIZED ORBITALS (OSLO)')
        call cpu_time(time)
        call get_wall_time(wtime)

!! Computing sat for Hilbert-space cases !!
        if(ilow2.ne.0) then
          iopt(5)=ilow2 !! MG: trick, not elegant but works !!
          if(ilow2.eq.1) call tomull(sat)
          if(ilow2.eq.2.or.ilow2.eq.3) call tolow(sat)
          if(ilow2.eq.6) call tonao(sat)
        end if

!! General, independently of the AIM scheme !!
        if(kop.eq.0) then
          call rwf_iterative_oslo(sat,itotps,wp,omp2,chp,pcoord)
        else
          call uwf_iterative_oslo(sat,itotps,wp,omp2,chp,pcoord)
        end if

        write(*,'(a)') ''
        call cpu_time(time2)
        call get_wall_time(wtime2)
        call print_timer('OSLO analysis',time2-time,wtime2-wtime)
      end if

!! X-ray scattering factors !!
      if(iscattfact.eq.1) then
        call print_box('EVALUATING X-RAY SCATTERING FACTORS')
        call scattering_factors(itotps,wp,rho,omp2,pcoord)
      end if

!! End printing... !!
      call json_close('normal','')
      write(*,'(a)') ''
      write(*,'(2x,a)') '...Normal Termination of APOST-3D...'

      end

