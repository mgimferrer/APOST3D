c-----------------------------------------------------------------------------
c                                                                                  
c                        Program APOST-3D, Version 5
c                                 22-09-2026
c                       --------------------------------                           
c                                                                                  
c        Real-space and Hilbert-space tools for wave function analysis             
c                                                                                  
c        Available atomic definitions:                                             
c        ----------------------------                                              
c                                                                                  
c        Real space:                                                               
c          Becke, J. Chem. Phys. 88 2547 1988                                      
c          Hirshfeld, Theor. Chim. Acta 44  129 1977                               
c          Hirshfeld-Iterative, J Chem Phys 126 144111 2007                        
c          Topological fuzzy Voronoi cells (TFVC), J Chem Phys 139 071103 2013    
c          QTAIM, J. Comput. Chem 30 1082 2009                                     
c                                                                                  
c        Hilbert-space : Mulliken, Lowdin, Davidson-Lowdin                         
c                                                                                  
c                                                                                  
c        Calculating:                                                              
c        ------------                                                              
c                                                                                  
c          A) Atomic and overlap populations, bond orders and valences            
c             I. Mayer and P. Salvador, Chem. Phys. Lett. 383 368-375 2004	      
c                                                                                  
c          B) Hartree-Fock molecular energy decomposition                          
c             P. Salvador, M. Duran, I.Mayer, J. Chem. Phys. 115 1153-1157 2001    
c             P. Salvador and I. Mayer, J. Chem. Phys. 120 5046-5052 2004          
c                                                                                  
c          C) KS-DFT molecular energy decomposition                                
c             P. Salvador, I. Mayer, J. Chem. Phys. 126 234113 2007	              
c             M. Gimferrer, P. Salvador, J. Chem. Phys. 158 234105 2023
c                                                                                  
c          D) Molecular energy decomposition for CAS/DMRG wavefunctions            
c                                                                                  
c          E) Effective atomic orbitals:                                           
c             I. Mayer, J. Phys. Chem. 100 6249 1996                               
c             I. Mayer and P. Salvador, J. Chem. Phys. 130 234106 2009             
c             E. Ramos-Cordoba et al., J. Chem. Phys. 138 214107 2013              
c                                                                                  
c          F) Local spin analysis                                                  
c             E. Ramos-Cordoba et al., J. Chem. Theory Comput. 8 1270-1279 2012   
c             E. Ramos-Cordoba et al., Phys. Chem. Chem. Phys. 14 15291-15298 2012 
c                                                                                  
c          G) Effective Oxidation states analysis                                  
c             E. Ramos-Cordoba et al., J. Chem. Theory Comput. 11 1501-1508 2015   
c
c          H) Oxidation states from localized orbitals
c             M. Gimferrer, G. Comas-Vila, P. Salvador, Molecules 25 234 2020
c             M. Gimferrer et al., Inorg. Chem. 59 15410-15420 2020
c             M. Gimferrer et al., J. Chem. Theor. Comput. 18 309-322 2022
c
c          I) Decomposition of EDA quantities into one- and two-center IQA terms
c             M. Gimferrer et al., J. Chem. Theory Comput. 19 3469-3485 2023
c
c          J) Origin-independent decomposition of static polarizabilities
c             M. Montilla, et al., J. Chem. Theory Comput. 17, 1098-1105 2021
c                                                                                  
c                                                                                  
c        Cite this program as:                                                     
C        ---------------------                                                     
c          P. Salvador, E. Ramos-Cordoba, M. Montilla, L. Pujal and M. Gimferrer 
c          J. Chem. Phys., 2024, 160, 172502 DOI: 10.1063/5.0206187                     
c                                                                                  
c      e-mail: psalse@gmail.com, mgimferrer18@gmail.com
c                                                                                  
c----------------------------------------------------------------------------------
c      The program has been written by using parts of the program APOST by         
c      I. Mayer and A. Hamza, Budapest, 2000-2003.                                 
C                                                                                  
c      The numerical integration utilizes the subroutines for Lebedev              
c      quadrature downloaded from CCL. The appropriate reference is:               
c      V.I. Lebedev, and D.N. Laikov "A quadrature formula for the sphere of the   
c      131st algebraic order of accuracy" Doklady Mathematics, 59 477-481 1999.    
c                                                                                  
C      The program makes use of libxc library when necessary, using the F90        
c      interfaces provided by the authors.                                         
c                     (see http://www.tddft.org/programs/libxc)                     
c                                                                                  
c      We are extremely grateful for the possibility of using these routines!      
c
c      -----------------------------------------------------------------------------c

!! **************************************************************** !!
!! PROGRAM APOST-3D -- MAIN ENTRY POINT / DRIVER                    !!
!! single unnamed main program, no PROGRAM statement, no internal   !!
!! subroutines -- this file is the program itself. Runs             !!
!! sequentially in five phases:                                     !!
!!   1) argument/file processing -- opens .fchk/.inp                !!
!!   2) .inp keyword parsing -- delegated to read_input() (see that !!
!!      file), filling input_options_mod's ~70 flags                !!
!!   3) cross-keyword validation                                    !!
!!   4) iopt(200) -- each local flag copied into its own hardcoded  !!
!!      iopt(N) slot, read back by every analysis routine           !!
!!   5) setup (Hilbert-space sat matrix or real-space integration   !!
!!      grid) then sequential dispatch per requested analysis:      !!
!!      population/local spin/POLAR/PCA/EFFAO-EOS/LOBA/ENPART       !!
!!      (one- then two-electron)/OSLO/X-ray scattering              !!
!! population/bond-order/PCA analysis are delegated to pop.f        !!
!! (population_density/population_print_charges/                    !!
!! population_print_overlap/bond_order_analysis/pca_analysis) --    !!
!! this is the final shape of the main.f/pop.f/read_input.f split.  !!
!! several sections are commented out, kept pending evaluation      !!
!! (DAFH, EDAIQA state-construction dispatch, DFT-DM1/HIRAO         !!
!! approximations, legacy topology/effao calls) -- catalogued, not  !!
!! deleted.                                                         !!
!! **************************************************************** !!

      use basis_set
      use ao_matrices
      use integration_grid
      use timing_mod
      use input_options_mod
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
!! local spin !!
      common /localspin/xlsa(maxat,maxat),ua(maxat)
!! Enpart !!
      common /exchg/exch(maxat,maxat),xmix
!! TFVC features !!
      common /achi/achi(maxat,maxat),ibcp
      common /erf/aerf,ierf
!! MODGRID (diatXC) features !!
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22
      common /modgrid2/thr3
!! EDAIQA features !!
      common /edaiqa/xen,xcoul,xnn
      common /edaiqa2/i2deda,iipoints,xptxyz(2,3)
!! NLOP features !!
      common /twoel/twoeltoler
      common /efield/field(4),edipole
!! printing and internal options !!
      common /filename/name0
      common /printout/iaccur
      common /iops/iopt(200)
!! for enpart !!
      dimension eto(maxat,maxat)
      character*60 name,name2,namepat,name3,name0
      character*80 line
!! for testing !!
      dimension xhess(3,3)

      allocatable wp(:),omp(:),omp2(:,:),chp(:,:),pcoord(:,:),rho(:)
      allocatable xkdens(:)
      allocatable ibaspoint(:)
      allocatable sss(:,:),sssi(:,:)
      allocatable sat(:,:,:)
      allocatable dm1(:,:),dm2(:,:,:,:)

!! TO CHANGE !!
      allocatable orbpop(:)
      allocatable rho_orb(:,:)
      allocatable rho_at(:,:)
!! !!


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
9999   stop 'The required input filename is missing'
      end if

!! print version info !!
      call kiir()

!! Processing .inp file !!
      open (16,file=name3,err=9999)
      open (15,file=name,err=9999)

      call read_input()

!! kop=1 -> unrestricted calculation; iposthf -> correlated             !!
!! calculation; icorr=0 -> no external dm1/2 provided; icas=1 ->        !!
!! CASSCF calculation; icisd -> CISD calculation; idono=1 -> do         !!
!! natural orbitals, idono=0 -> restricted SD calculation.              !!

!! Dependencies and cross-keyword validation !!
      iposthf=0
      idono=0
      if(icas.eq.1.or.icisd.eq.1) iposthf=1
      iopt(65)=iposthf
!! GEOS (effao3d_u) needs occ_no/c_no even for a restricted SD          !!
!! wavefunction -- there they're just the canonical MOs (integer        !!
!! occupied), giving a trivially ~0 unpaired channel, same as a plain   !!
!! restricted-wavefunction EFFAO run -- so force natural orbitals on    !!
!! for GEOS regardless of kop.                                         !!
      if(iposthf.eq.1.or.kop.eq.1.or.iueos.eq.1) idono=1

      if(iposthf.eq.1) then
        if(ispin.eq.1.and.icorr.lt.2) stop ' Local Spin needs dm1 and dm2 for correlated WFs'
        if(ienpart.eq.1.and.icorr.lt.2) stop ' Enpart needs dm1 and dm2 for correlated WFs'
      end if
      if(iqtaim.eq.1) stop'This version can not do QTAIM'

      if (ipca.eq.1.and.iqtaim.ne.1) iopop=1
      if(imulli.gt.1.or.iqtaim.eq.1) iopop=0
      if(ihirsh.ne.0.and.idoat.eq.1) stop'Cant do HIRSH with DOATOMS'
!! gated on the real restricted-SD condition directly, not on idono --  !!
!! idono can now also be forced on by GEOS (see above) on a genuinely   !!
!! restricted wavefunction, where Local Spin still doesn't make sense   !!
!! (confirmed: produces nonsense u_A/N_D values, not a trivial zero).   !!
      if(ispin.eq.1.and.kop.eq.0.and.iposthf.eq.0)  then
        write(*,*) 'No Local Spin Analysis needed for Restricted SD WFs'
        ispin=0
      end if
      if(imulli.ne.0.and.ienpart.ne.0)  then
        write(*,*) 'Can not do ENPART with Hilbert-space analysis'
        ienpart=0
      end if
      if(iqtaim.ne.0.and.ienpart.ne.0)  then
        write(*,*) 'Can not do ENPART with QTAIM '
        ienpart=0
      end if
      if(ieos.eq.1.and.idoat.eq.1) stop 'Cant do EOS with DOATOMS'

      if(idoint.eq.1) then
        write(*,*) ' Will do atomic overlaps for FCALC'
        if(iwfn.eq.1) then
          write(*,*) ' Will use orbitals from wfn file'
          write(*,*) ' Assuming orbitals are on fort.92'
c possible call system here...
        end if
      end if
      
      if(kop.eq.0.and.ispin.eq.1) then
        write(*,*) 'Warning, closed-shell Local Spin calculation ' 
      end if

!! Start of iopt listing. Full 1-200 listing so unused slots are   !!
!! visible at a glance -- iopt(200). Commented lines are genuinely !!
!!unused slots.                                                    !!
      iopt(1) = idoint
c      iopt(2) =
      iopt(3) = iwfn
      iopt(4) = idono
      iopt(5) = imulli
      iopt(6) = ihirsh
      iopt(7) = iallpo
      iopt(8) = iopop
c      iopt(9) =   !! ndens0, set above in read_input() !!
c      iopt(10) =

c      iopt(11) =
      iopt(12) = ieffao
      iopt(13) = icube
      iopt(14) = ibcp
      iopt(15) = ispin
      iopt(16) = iqtaim
      iopt(17) = ienpart
      iopt(18) = ihf
c      iopt(19) =
      iopt(20) = iexact

      iopt(21) = ihomo
      iopt(22) = idek
      iopt(23) = iionic
      iopt(24) = ieffthr
      iopt(25) = istiff
      iopt(26) = icorr
      iopt(27) = isha
      iopt(28) = ipca
c      iopt(29) =
      iopt(30) = idafh

      iopt(31) = inewbec
c      iopt(32) =
c      iopt(33) =
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
      iopt(45) = ifinegrid
      iopt(46) = ipolar
      iopt(47) = ifield
      iopt(48) = iradmat
      iopt(49) = ielcount !MMO- NCTAIM
c      iopt(50) =

c      iopt(51) =
c      iopt(52) =
c      iopt(53) =
c      iopt(54) =
c      iopt(55) =   !! itype/jtype, set below in ENPART dispatch (dynamic) !!
      iopt(56) = idftdm1
      iopt(57) = id_func_dm1 !! dft_dm1.f's own functional-ID read (ifunc=Iopt(57)) !!
      iopt(58) = itop
      iopt(59) = iecorr
      iopt(60) = id_xcfunc

      iopt(61) = id_xfunc
      iopt(62) = id_cfunc
      iopt(63) = inatorb_dm1 !! dft_dm1.f's own NATORB read !!
c      iopt(64) =
c      iopt(65) =   !! iposthf, set above in DEPENDENCIES & TO DO !!
c      iopt(66) =
c      iopt(67) =
c      iopt(68) =
c      iopt(69) =
c      iopt(70) =

c      iopt(71) =
c      iopt(72) =
c      iopt(73) =
c      iopt(74) =
c      iopt(75) =
c      iopt(76) =
c      iopt(77) =
c      iopt(78) =
      iopt(79) = imokit !! MG: Temporarily !!
c      iopt(80) =

c      iopt(81) =
c      iopt(82) =
c      iopt(83) =
c      iopt(84) =
      iopt(85) = ipairs
      iopt(86) = ietop
      iopt(87) = ipyscf
      iopt(88) = inopop
c      iopt(89) =
      iopt(90) = ieoscent

      iopt(91) = ianalytical
      iopt(92) = iedaiqa
      iopt(93) = iflip
c      iopt(94) =
      iopt(95) = iqchem
      iopt(96) = ifolitol
      iopt(97) = ibranch
      iopt(98) = ioslofchk
      iopt(99) = iloba
c      iopt(100) =

c      iopt(101) =
c      iopt(102) =
c      iopt(103) =
c      iopt(104) =
c      iopt(105) =
c      iopt(106) =
c      iopt(107) =
c      iopt(108) =
c      iopt(109) =
c      iopt(110) =

c      iopt(111) =
c      iopt(112) =
c      iopt(113) =
c      iopt(114) =
c      iopt(115) =
c      iopt(116) =
c      iopt(117) =
c      iopt(118) =
c      iopt(119) =
c      iopt(120) =

c      iopt(121) =
c      iopt(122) =
c      iopt(123) =
c      iopt(124) =
c      iopt(125) =
c      iopt(126) =
c      iopt(127) =
c      iopt(128) =
c      iopt(129) =
c      iopt(130) =

c      iopt(131) =
c      iopt(132) =
c      iopt(133) =
c      iopt(134) =
c      iopt(135) =
c      iopt(136) =
c      iopt(137) =
c      iopt(138) =
c      iopt(139) =
c      iopt(140) =

c      iopt(141) =
c      iopt(142) =
c      iopt(143) =
c      iopt(144) =
c      iopt(145) =
c      iopt(146) =
c      iopt(147) =
c      iopt(148) =
c      iopt(149) =
c      iopt(150) =

c      iopt(151) =
c      iopt(152) =
c      iopt(153) =
c      iopt(154) =
c      iopt(155) =
c      iopt(156) =
c      iopt(157) =
c      iopt(158) =
c      iopt(159) =
c      iopt(160) =

c      iopt(161) =
c      iopt(162) =
c      iopt(163) =
c      iopt(164) =
c      iopt(165) =
c      iopt(166) =
c      iopt(167) =
c      iopt(168) =
c      iopt(169) =
c      iopt(170) =

c      iopt(171) =
c      iopt(172) =
c      iopt(173) =
c      iopt(174) =
c      iopt(175) =
c      iopt(176) =
c      iopt(177) =
c      iopt(178) =
c      iopt(179) =
c      iopt(180) =

c      iopt(181) =
c      iopt(182) =
c      iopt(183) =
c      iopt(184) =
c      iopt(185) =
c      iopt(186) =
c      iopt(187) =
c      iopt(188) =
c      iopt(189) =
c      iopt(190) =

c      iopt(191) =
c      iopt(192) =
c      iopt(193) =
c      iopt(194) =
c      iopt(195) =
c      iopt(196) =
c      iopt(197) =
c      iopt(198) =
c      iopt(199) =
c      iopt(200) =
!! End of iop listing !!

!! digested echo of the .inp file -- every flag above is now final !!
      call print_input_summary()

!! Natural orbitals from P-matrix in .fchk !!
      if (idono.eq.1) call gennatural()

      write(*,*)
      call cpu_time(time2)
      call get_wall_time(wtime2)
      call print_timer('initialization',time2-time,wtime2-wtime)
      time=time2
      wtime=wtime2

!! density at the iatdens atom !!
      if(iatdens.ne.0) call atdens_int(Rmax,iatdens)
      if(inopop.eq.1) stop 'Normal termination of APOST3D'

!! set atomic radii just in case needed !!
      call prepar()

!! Mulliken-type analysis, either int files for FCALC or effao/readint !!
!! files                                                                !!
      ndim=igr
      if(imulli.ge.1) then
        allocate (sat(ndim,ndim,nat))
        if(imulli.eq.1) then
          call tomull(sat)
        else if(imulli.eq.2.or.imulli.eq.3) then
          write(*,*) ' Doing Hilbert-space analysis in Lowdin basis'
          call tolow(sat)
        else if(imulli.eq.4) then
          write(*,*) ' Doing Hilbert-space analysis in NAO basis'
          call tonao(sat)
        else if(imulli.eq.5) then
          write(*,*) ' Doing Hilbert-space analysis in weighted Lowdin ba
     +      sis'
          call tolow2(sat)
        end if

!      else if(iqtaim.eq.2) then
!       ndim=nocc
!       allocate (sat(ndim,ndim,nat))
!       call readintfiles(sat)
       
      else

!! Prepare for numerical integrations !!
        call build_integration_grid(ienpart, ipolar,ifinegrid)

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
        call print_box('SETTING ATOMIC DEFINITION')
        if(ierf.ne.1) then
          write(*,'(2x,a,1x,i0)') 'Using stiffness k:',istiff
        else
          write(*,'(2x,a,1x,f8.4)') 'Using adjustable profile with A =',aerf
        end if
        if(ibcp.eq.1) call khi()
        if(ihirsh.eq.1.or.ihirsh.eq.2.or.ielcount.eq.1) call makeatdens
 
        iiter=1
        call prenumint(ndim,itotps,nat,wp,omp,omp2,chp,rho,pcoord,ibaspoint,iiter)

!! OS from localized MOs (centroids) !!
        if(ieoscent.eq.1) then
          write(*,*) 'Doing OS from centroids of localized orbitals...'
          call eos_centroid(itotps,chp,wp,omp,pcoord)
          stop
        end if

!! Integrate atomic overlap by default !!
        allocate (sat(ndim,ndim,nat))
        call numint_sat(ndim,itotps,nat,wp,omp,omp2,chp,ibaspoint,sat)

      end if
!! end of preparation for numerical integrations !!

!! atomic-domain partition time, covers both branches above (Hilbert    !!
!! and real-space).                                                     !!
      write(*,*)
      call cpu_time(time2)
      call get_wall_time(wtime2)
      call print_timer('atomic definition',time2-time,wtime2-wtime)
      time=time2
      wtime=wtime2

!! Correlated-WF input, needs dm1 and/or dm2 files !!
!! cas/cisd specifications !!
      nelec=nalf+nb
      if(icorr.ne.0) then
        call print_box('POST-HARTREE-FOCK CALCULATION')

        if(icisd.eq.1) nspinorb=nbasis*2

        write(*,*) 'Number of core + active spin-orbitals : ',nspinorb
        write(*,*) 'Number of electrons : ',nelec
        write(*,*) 'Number of basis functions :',nbasis
        write(*,*)

        write(*,*) 'DM1 input starts'
        ALLOCATE(dm1(nspinorb,nspinorb))
        call dm1input(dm1)
        if(icorr.eq.2) then
          write(*,*) 'DM2 input starts'
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

!       write(*,*) 'idoint',idoint
      if(idoint.eq.1)  call print_int(ndim,nat,sat,namepat)

      call population_print_charges()

      if(ielcount.eq.1) then
        iatps=nang*nrad
        CALL NCTAIM(iatps,wp,omp,omp2,nat,sat,igr,ibaspoint,chp)
      end if

      call population_print_overlap(wp,omp,omp2,rho)

!! Bond orders, valences, number of effectively unpaired electrons !!
      call bond_order_analysis(sat)

      write(*,*)
      call cpu_time(time2)
      call get_wall_time(wtime2)
      call print_timer('population analyses',time2-time,wtime2-wtime)
      time=time2
      wtime=wtime2

!! Local spin decomposition, single-determinant WF !!
      if(ispin.eq.1.and.icas.eq.0.and.icisd.eq.0) then
        call print_box('DOING LOCAL SPIN ANALYSIS')
        write(*,'(3x,a)') 'Single-determinant case'
        call fspindec(sat)
        if (idofr.eq.1) then 
          line ='   FRAGMENT ANALYSIS : Local Spin Analysis'
          call group_by_frag_mat(0,line ,xlsa)
          line ='   FRAGMENT ANALYSIS : Num. eff. unpaired elec.'
          call group_by_frag_vec(1,line ,ua)
        end if
      end if

!! Local spin and DIs for correlated WFs, needs dm1/dm2 from the DMN   !!
!! code (E. Matito) or from pySCF                                      !!
      if((icas.eq.1.or.icisd.eq.1).and.ispin.eq.1)then
        call print_box('DOING LOCAL SPIN ANALYSIS')
        write(*,'(3x,a)') 'Localization/delocalization, correlated WF'
        call spincorr(sat,dm1,dm2)
        call cpu_time(time2)
        call get_wall_time(wtime2)
        call print_timer('local spin analysis',time2-time,wtime2-wtime)
        time=time2
        wtime=wtime2
      end if

!! Nonlinear optical properties (POLAR) !!
      if(ipolar.ne.0) call polar(itotps,nat,wp,omp,omp2,pcoord,rho)

!! Entropies and correlation indicators -- unfinished, see numint_sha !!
c      if(isha.ne.0) call numint_sha(ndim,itotps,nat,wp,chp,omp,omp2,ibaspoint)

!! PCA analysis -- needs DI already populated (bond_order_analysis above) !!
      if(ipca.eq.1) call pca_analysis()

!! DAFH part -- needed files produced by external code (R. Ponec), not !!
!! available in this version                                           !!
c          if(idafh.eq.1) then
c           ncactiv=nspinorb
c           if(icorr.ne.0) then
c            call dafh_input(ncactiv,nbasis,nat,sat)
c           else if (kop.eq.1) then
c            call dafh_input_uhf(nbasis,nat,sat)
c           else
c            call dafh_input_rhf(nbasis,nat,sat)
c           end if
c          end if

!! EFFAO part                                                          !!
!! ieffao: 0 nothing, 1 eff-AOs, 2 spin-resolved eff-AOs (a must for   !!
!!   EOS), 3 paired/unpaired eff-AOs                                   !!
!! imulli: 1 Mulliken, 2 Lowdin, 3 Lowdin-Davidson (not implemented),   !!
!!   4 NAO                                                              !!
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
            if(kop.ne.0.or.(icas.eq.1.and.nalf.ne.nb)) then
              call ueffaomull_frag(2)
              idobeta=1
            end if
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
            if(kop.ne.0.or.(icas.eq.1.and.nalf.ne.nb)) then
              call ueffaolow_frag(2)
              idobeta=1
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
!! disabled call sites for effao.f's uefomo (dead, see that file's   !!
!! top-of-file note) and devel.f's mhg/mhg2 (dead, part of that      !!
!! file's own confirmed-dead subroutines, no note there yet).        !!
c             call  uefomo(itotps,ndim,omp,chp,sat,wp,omp2,0)
c             call mhg(itotps,ndim,omp,chp,sat,wp,omp2,pcoord,pa,0)
c             call mhg2(itotps,ndim,omp,chp,sat,wp,omp2,pcoord,p,0)
            else if(ieffao.eq.2) then
              idobeta=0
              write(*,*) '  '
              write(*,*) ' UEFFAO: alpha and beta treated separately'
              write(*,*) '  '
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
              if(kop.ne.0.or.(icas.eq.1.and.icass.ne.0)) then
                idobeta=1
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

        write(*,*)
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
        if(imulli.gt.0) then
          write(*,*) " LOBA NOT IMPLEMENTED FOR HILBERT-SPACE "
          stop
        end if

!! Real-space !!
        if(imulli.eq.0) call eos_loba(sat)
      end if

!! Energy decomposition (ENPART) !!
      if(ienpart.eq.1) then
        call print_box('DOING MOLECULAR ENERGY DECOMPOSITION')

!! One-electron part: CASSCF and CI WFs !!
        if(iposthf.eq.1) then
          call numint_one_rphf(ndim,itotps,wp,omp2,pcoord,chp,rho,eto)
          call cpu_time(time2)
          call get_wall_time(wtime2)
          call print_timer('enpart one-electron',time2-time,wtime2-wtime)
          time=time2
          wtime=wtime2

!! One-electron part: DFT and HF WFs !!
        else 

!! Initialize DFT functional for info and initial printing !!
          if(id_xfunc.ne.-1) then
            ifuncfirst=1
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
            write(*,*) " "
          end if

!! Restricted case !!
          if(kop.ne.1) then 
            ALLOCATE(xkdens(itotps)) ! (TO DO) Rethink how to include it... only used in metaGGA functionals

!! One-electron terms !!
            call numint_one(ndim,itotps,wp,rho,omp,omp2,pcoord,chp,eto)
            call cpu_time(time2)
            call get_wall_time(wtime2)
            call print_timer('enpart one-electron',time2-time,wtime2-wtime)
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
              call print_timer('enpart dft',time2-time,wtime2-wtime)
              time=time2
              wtime=wtime2
            end if
            DEALLOCATE(xkdens)

!! Unrestricted case !!
          else
            ALLOCATE(xkdens(itotps)) ! (TO DO) Rethink how to include it... only used in metaGGA functionals

!! One-electron terms !!
            call numint_one_uhf(ndim,itotps,wp,rho,omp,omp2,pcoord,chp,eto)
            call cpu_time(time2)
            call get_wall_time(wtime2)
            call print_timer('enpart one-electron',time2-time,wtime2-wtime)
            time=time2
            wtime=wtime2

!! DFT XC term !!
            if(id_xfunc.ne.-1) then
              call numint_dft_uks(ndim,itotps,wp,omp,omp2,chp,pcoord,eto)
              call cpu_time(time2)
              call get_wall_time(wtime2)
              call print_timer('enpart dft',time2-time,wtime2-wtime)
              time=time2
              wtime=wtime2
            end if
            DEALLOCATE(xkdens)
          end if 

!! End of restricted/unrestricted dispatch !!
        end if
        DEALLOCATE(wp,omp,omp2,chp,pcoord,ibaspoint,rho)

!! Two-electron part !!

!! Two-electron integration defaults !!
        call print_box('SETTING GRID FOR TWO-ELECTRON NUMERICAL INTEGRATION')

!! Controlled by # GRID option (modgrid common) !!
!! Default grid is now 150/590, can be changed to 40/146 but ensure to also modify pha and phb !!
        nrad=nrad22
        nang=nang22
        rr00=rr0022

!! Analytical case -- only step up to the 70/434 default when the user   !!
!! hasn't configured # GRID themselves (ienpart_gridtwoel.eq.0); this     !!
!! used to unconditionally overwrite nrad/nang, silently discarding an   !!
!! explicit # GRID (and, since the plain default grew to 150/590, this   !!
!! was a downgrade in that case, not the increase the message claims).   !!
        if(ianalytical.eq.1.and.ienpart_gridtwoel.eq.0) then
          write(*,*) " Analytical calculation has been requested: Increasing grid because 2-electron is now 1-electron "
          write(*,*) " "
          nrad=70
          nang=434
        else if(ianalytical.eq.1) then
          write(*,*) " Analytical calculation has been requested: keeping the user-configured # GRID (MOD-GRIDTWOEL) "
          write(*,*) " "
        end if

!! Printing info !!
        write(*,'(2x,a14,x,i4)') "Radial points:",nrad
        write(*,'(2x,a15,x,i4)') "Angular points:",nang
        write(*,'(2x,a12,x,i9)') "Grid points:",nrad*nang*nat

!! Generating grid for two-electron numerical integrations !!
        iatps=nang*nrad
        itotps=nrad*nang*nat
        call quad(Nrad,Nang) 
        ALLOCATE(wp(itotps),omp(itotps),omp2(itotps,nat),rho(itotps))
        ALLOCATE(pcoord(itotps,3),ibaspoint(itotps),chp(itotps,ndim))
        call prenumint(ndim,itotps,nat,wp,omp,omp2,chp,rho,pcoord,ibaspoint,0)

!! CASSCF and CI WFs !!
        if(iposthf.eq.1) then
!        if(itop.eq.1) call top_3d(norb,2,0,iatpairs) !! MG: TOPOLOGY ROUTINES NEEDS A CHECK !!
          call numint_two_rphf(ndim,itotps,wp,omp2,pcoord,chp,rho,eto,dm1,dm2)
          call cpu_time(time2)
          call get_wall_time(wtime2)
          call print_timer('enpart two-electron',time2-time,wtime2-wtime)
          time=time2
          wtime=wtime2

!! DFT and HF WFs !!
        else

!! Unrestricted case !!
          if(kop.eq.1) then 
            call numint_two_uhf(ndim,itotps,wp,omp,omp2,pcoord,chp,rho,eto)

!! Restricted case !!
          else
!          if(itop.eq.1) call top_3d(nocc,1,0,iatpairs) !! MG: TOPOLOGY ROUTINES NEEDS A CHECK !!
            call numint_two(ndim,itotps,wp,omp,omp2,pcoord,chp,rho,eto)
          end if
          write(*,*) " "
          call cpu_time(time2)
          call get_wall_time(wtime2)
          call print_timer('enpart two-electron',time2-time,wtime2-wtime)
          time=time2
          wtime=wtime2

!! End of WF-type dispatch !!
        end if 

!! Deallocating any existent grid !!
        DEALLOCATE(wp,omp,omp2,chp,pcoord,ibaspoint,rho)

      end if
!! end of energy decomposition !!

!! more topology stuff can be found in version 3.1-devel !!

!! EDAIQA -- to do, creating the E(<A^0 B^0>)^AB state !!
!     if(iedaiqa.eq.1) then
!       write(*,*) " "
!       write(*,*) " Entering EDAIQA Section "
!       write(*,*) " "
!       if(kop.eq.0) call rwf_edatoiqafchk()
!       if(kop.eq.1) call uwf_edatoiqafchk()
!       close(51)
!       close(52)
!     end if

!! DFT-DM1 approximate one-particle RDM1 (formerly HIRAO internally).    !!
!! Runs standalone, own dedicated grid reusing the two-electron-type      !!
!! MOD-GRIDTWOEL settings (read_input.f's read_gridtwoel); dft_dm1        !!
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

!! OSLO -- variants of the procedure can be found in the dev version !!
      if(ioslo.eq.1) then
        call print_box('DOING OXIDATION STATES FROM LOCALIZED ORBITALS (OSLO)')
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

        call cpu_time(time2)
        call get_wall_time(wtime2)
        call print_timer('OSLO analysis',time2-time,wtime2-wtime)

        DEALLOCATE(wp,omp,omp2,chp,pcoord)
      end if

!! X-ray scattering factors !!
      if(iscattfact.eq.1) then
        call print_box('EVALUATING X-RAY SCATTERING FACTORS')
        call scattering_factors(itotps,wp,rho,omp2,pcoord)
      end if

!! End printing... !!
      write(*,*)
      write(*,'(2x,a)') '...Normal Termination of APOST-3D...'

      end

