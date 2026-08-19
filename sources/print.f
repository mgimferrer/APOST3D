!! ********************************************************************* !!
!! FILE STATUS (2026-08-19): Subroutine Cleanup Protocol still NOT       !!
!! applied to: cubegen3, rmat, rarr, ival, cubegen3_mhg. All confirmed   !!
!! dead or near-dead -- only reachable from dead or DOATOMS-adjacent     !!
!! callers -- deliberately left untouched pending a consolidation        !!
!! decision, not yet made. print_int_old (zero callers, superseded by    !!
!! the live print_int) was confirmed dead and deleted outright rather    !!
!! than left pending, since nothing referenced it. Every other           !!
!! subroutine in this file is done.                                      !!
!! ********************************************************************* !!

!! ********************************************************************* !!
!! subroutine: print_box                                                 !!
!! purpose: prints text inside a rule auto-sized to fit it.              !!
!! arguments:                                                            !!
!!   text (in) -- title to print, no leading/trailing padding needed     !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      SUBROUTINE print_box(text)
      IMPLICIT NONE
      character(len=*), intent(in) :: text
      integer :: n

      n=len_trim(text)+4
      write(*,'(/,2x,a)')  repeat('-',n)
      write(*,'(2x,2x,a)') trim(text)
      write(*,'(2x,a,/)')  repeat('-',n)

      END SUBROUTINE print_box

!! ***** !!

!! ********************************************************************* !!
!! subroutine: kiir                                                      !!
!! purpose: prints the startup banner -- version/date identity, feature  !!
!! availability and per-method citation directory, acknowledgments.      !!
!! arguments: none                                                       !!
!! author:                                                               !!
!! ********************************************************************* !!
      SUBROUTINE kiir
      IMPLICIT NONE
      write(*,'(a80)') '------------------------------------------------------------------------------'
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '                   ___    ____  ____  ___________   _____ ____                '
      write(*,'(a80)') '                  /   |  / __ \/ __ \/ ___/_  __/  |__  // __ \               '
      write(*,'(a80)') '                 / /| | / /_/ / / / /\__ \ / /_____ /_ </ / / /               '
      write(*,'(a80)') '                / ___ |/ ____/ /_/ /___/ // /_____/__/ / /_/ /                '
      write(*,'(a80)') '               /_/  |_/_/    \____//____//_/     /____/_____/                 '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '                      Version 5   --   22 September 2026                      '
      write(*,'(a80)') '                       --------------------------------                       '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Real-space and Hilbert-space tools for wave function analysis             '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Documentation: https://apost3d.readthedocs.io                             '
      write(*,'(a80)') '    Source code & issue tracker: https://github.com/mgimferrer/APOST3D        '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Cite this program as:                                                     '
      write(*,'(a80)') '    ---------------------                                                     '
      write(*,'(a80)') '      P. Salvador, E. Ramos-Cordoba, M. Montilla, L. Pujal and M. Gimferrer   '
      write(*,'(a80)') '      J. Chem. Phys., 2024, 160, 172502 DOI: 10.1063/5.0206187                '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '  e-mail: psalse@gmail.com, eloy.raco@gmail.com, mgimferrer18@gmail.com       '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Available atomic definitions:                                             '
      write(*,'(a80)') '    ----------------------------                                              '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Real space:                                                               '
      write(*,'(a80)') '      Becke, J. Chem. Phys. 88 2547 1988                                      '
      write(*,'(a80)') '      Hirshfeld, Theor. Chim. Acta 44  129 1977                               '
      write(*,'(a80)') '      Hirshfeld-Iterative, J Chem Phys 126 144111 2007                        '
      write(*,'(a80)') '      Topological fuzzy Voronoi cells (TFVC), J Chem Phys 139 071103 2013     '
      write(*,'(a80)') '      QTAIM, J. Comput. Chem 30 1082 2009                                     '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Hilbert-space : Mulliken, Lowdin, Davidson-Lowdin                         '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Calculating:                                                              '
      write(*,'(a80)') '    ------------                                                              '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      A) Atomic and overlap populations, bond orders and valences             '
      write(*,'(a80)') '         I. Mayer and P. Salvador, Chem. Phys. Lett. 383 368-375 2004         '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      B) Hartree-Fock molecular energy decomposition                          '
      write(*,'(a80)') '         P. Salvador, M. Duran, I.Mayer, J. Chem. Phys. 115 1153-1157 2001    '
      write(*,'(a80)') '         P. Salvador and I. Mayer, J. Chem. Phys. 120 5046-5052 2004          '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      C) KS-DFT molecular energy decomposition                                '
      write(*,'(a80)') '         P. Salvador, I. Mayer, J. Chem. Phys. 126 234113 2007                '
      write(*,'(a80)') '         M. Gimferrer, P. Salvador, J. Chem. Phys. 158 234105 2023            '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      D) Molecular energy decomposition for CAS/DMRG wavefunctions            '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      E) Effective atomic orbitals:                                           '
      write(*,'(a80)') '         I. Mayer, J. Phys. Chem. 100 6249 1996                               '
      write(*,'(a80)') '         I. Mayer and P. Salvador, J. Chem. Phys. 130 234106 2009             '
      write(*,'(a80)') '         E. Ramos-Cordoba et al., J. Chem. Phys. 138 214107 2013              '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      F) Local spin analysis                                                  '
      write(*,'(a80)') '         E. Ramos-Cordoba et al., J. Chem. Theory Comput. 8 1270-1279 2012    '
      write(*,'(a80)') '         E. Ramos-Cordoba et al., Phys. Chem. Chem. Phys. 14 15291-15298 2012 '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      G) Effective Oxidation states analysis                                  '
      write(*,'(a80)') '         E. Ramos-Cordoba et al., J. Chem. Theory Comput. 11 1501-1508 2015   '
      write(*,'(a80)') '         M. Gimferrer and P. Salvador, manuscript in preparation, 2026        '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      H) Oxidation states from localized orbitals                             '
      write(*,'(a80)') '         M. Gimferrer, G. Comas-Vila, P. Salvador, Molecules 25 234 2020      '
      write(*,'(a80)') '         M. Gimferrer et al., Inorg. Chem. 59 15410-15420 2020                '
      write(*,'(a80)') '         M. Gimferrer et al., J. Chem. Theor. Comput. 18 309-322 2022         '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      I) Decomposition of EDA quantities into one- and two-center IQA terms   '
      write(*,'(a80)') '         M. Gimferrer et al., J. Chem. Theory Comput. 19 3469-3485 2023       '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      J) Origin-independent decomposition of static polarizabilities          '
      write(*,'(a80)') '         M. Montilla, et al., J. Chem. Theory Comput. 17, 1098-1105 2021      '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '------------------------------------------------------------------------------'
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '  The program has been written by using parts of the program APOST by         '
      write(*,'(a80)') '  I. Mayer and A. Hamza, Budapest, 2000-2003.                                 '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '  The numerical integration utilizes the subroutines for Lebedev              '
      write(*,'(a80)') '  quadrature downloaded from CCL. The appropriate reference is:               '
      write(*,'(a80)') '  V.I. Lebedev, and D.N. Laikov "A quadrature formula for the sphere of the   '
      write(*,'(a80)') '  131st algebraic order of accuracy" Doklady Mathematics, 59 477-481 1999.    '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '  The program makes use of libxc library when necessary, using the F90        '
      write(*,'(a80)') '  interfaces provided by the authors.                                         '
      write(*,'(a80)') '                 (see https://libxc.gitlab.io/)                               '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '  We are extremely grateful for the possibility of using these routines!      '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '------------------------------------------------------------------------------'
      END SUBROUTINE kiir

!! ***** !!

!! ********************************************************************* !!
!! subroutine: print_input_summary                                       !!
!! purpose: "digested" echo of the .inp file -- prints which keywords    !!
!! were actually active, grouped by .inp section, mirroring              !!
!! tests/keywords.json's registry. Every keyword in the registry is      !!
!! accounted for internally; only active ones print, sections with       !!
!! nothing active are skipped entirely. Called from main.f right after   !!
!! OPTIONS LIST, once every flag below is finalized.                     !!
!! arguments: none (all via input_options_mod + a handful of COMMONs)    !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      SUBROUTINE print_input_summary()
      use input_options_mod
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /achi/achi(maxat,maxat),ibcp
      common /erf/aerf,ierf
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22
      common /atlist/iatlist(maxat),icuat
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,
     +  jfrlist(maxat)
      common /printout/iaccur
      character*200 cbuf
      character*40 cval
      character*40 chdr

      call print_box('INPUT SUMMARY')

!! ----------------------------------------------------------------- !!
!! # METHOD -- one line per active category, "label : value(s)"      !!
!! ----------------------------------------------------------------- !!
      chdr='# METHOD'
      write(*,'(2x,a)') trim(chdr)
      write(*,'(2x,a)') repeat('-',len_trim(chdr))

!! real-space AIM schemes -- several can be active at once, joined   !!
!! into one line; itfvc/ibcp/inewbec aren't mutually exclusive since !!
!! TFVC sets ibcp/inewbec as a side effect (read_input.f), so TFVC   !!
!! is checked first and BECKE-RHO/NEWBEC are only reported when they !!
!! are the reason those flags are set, not a TFVC side effect        !!
      cbuf=' '
      if(itfvc.eq.1) cbuf=trim(cbuf)//' TFVC (Topological Fuzzy Voronoi Cells);'
      if(ibcp.eq.1.and.itfvc.ne.1)
     +  cbuf=trim(cbuf)//' BECKE-RHO;'
      if(inewbec.eq.1.and.itfvc.ne.1)
     +  cbuf=trim(cbuf)//' NEWBEC;'
      if(ihirsh.eq.1) cbuf=trim(cbuf)//' HIRSH (Hirshfeld);'
      if(ihirsh.eq.2) cbuf=trim(cbuf)//' HIRSH-IT (Hirshfeld-Iterative);'
      if(iqtaim.eq.1) cbuf=trim(cbuf)//' QTAIM;'
      if(iqtaim.eq.2) cbuf=trim(cbuf)//' QTAIM (READINT, reusing prior results);'
      cbuf=adjustl(cbuf)
      nlen=len_trim(cbuf)
      if(nlen.gt.0.and.cbuf(nlen:nlen).eq.';') cbuf(nlen:nlen)=' '
      if(len_trim(cbuf).gt.0) write(*,'(2x,a,1x,a)')
     +  'Atomic partitioning (real-space)     :',trim(cbuf)

      if(imulli.ge.1) then
        if(imulli.eq.1) cval='Mulliken'
        if(imulli.eq.2) cval='Lowdin'
        if(imulli.eq.3) cval='Lowdin-Davidson'
        if(imulli.eq.4) cval='NAO-basis'
        if(imulli.eq.5) cval='Lowdin-W'
        write(*,'(2x,a,1x,a)')
     +    'Atomic partitioning (Hilbert-space)  :',trim(cval)
      end if

      if(iopop.eq.1) write(*,'(2x,a,1x,a)')
     +  'Overlap population analysis          :','OPOP'

      if(ieffao.ge.1.and.ieos.ne.1.and.iueos.ne.1) then
        if(ieffao.eq.1) cval='EFFAO'
        if(ieffao.eq.2) cval='UEFFAO'
        if(ieffao.eq.3) cval='EFFAO-U (paired/unpaired only)'
        write(*,'(2x,a,1x,a)')
     +    'Effective atomic orbitals            :',trim(cval)
      end if

      if(ieos.eq.1) write(*,'(2x,a,1x,a)')
     +  'Oxidation states analysis            :',
     +  'EOS (Effective Oxidation States, fragment-based)'
      if(iueos.eq.1) write(*,'(2x,a,1x,a)')
     +  'Oxidation states analysis            :',
     +  'EOS-U (open-shell, unpaired-density-based)'
      if(ieoscent.eq.1) write(*,'(2x,a,1x,a)')
     +  'Oxidation states analysis            :','OS-CENTROID'
      if(ioslo.eq.1) write(*,'(2x,a,1x,a)')
     +  'Oxidation states analysis            :',
     +  'OSLO (oxidation states localized orbitals)'
      if(iloba.eq.1) write(*,'(2x,a,1x,a)')
     +  'Oxidation states analysis            :',
     +  'LOBA (Localized Orbital Bonding Analysis)'

      if(ienpart.eq.1) write(*,'(2x,a,1x,a)')
     +  'Energy decomposition                 :',
     +  'ENPART (one- and two-center IQA terms, see below)'
      if(iedaiqa.eq.1) write(*,'(2x,a,1x,a)')
     +  'Energy decomposition                 :',
     +  'EDAIQA (EDA interaction energies into IQA terms)'

      if(ispin.eq.1) write(*,'(2x,a,1x,a)')
     +  'Local spin analysis                  :','SPIN'
      if(ipolar.eq.1) write(*,'(2x,a,1x,a)')
     +  'Static polarizability                :','POLAR'
      if(itop.eq.1) write(*,'(2x,a,1x,a)')
     +  'Topology analysis                    :',
     +  'TOPOLOGY (see below)'
      if(idafh.eq.1) write(*,'(2x,a,1x,a)')
     +  'Domain-averaged Fermi holes          :',
     +  'DAFH (needs external files, see kiir banner)'
      if(iscattfact.eq.1) write(*,'(2x,a,1x,a)')
     +  'X-ray scattering factors             :','SCATT-FACT'
      if(ilaplacian.eq.1) write(*,'(2x,a,1x,a)')
     +  'Density Laplacian                    :','LAPLACIAN'
      if(isha.eq.1) write(*,'(2x,a,1x,a)')
     +  'Shannon entropy decomposition        :','SHANNON'
      if(ipca.eq.1) write(*,'(2x,a,1x,a)')
     +  'Principal component analysis         :','PCA'
      if(idoint.eq.1) write(*,'(2x,a,1x,a)')
     +  'Integration diagnostics              :','DOINT'
      if(ielcount.eq.1) write(*,'(2x,a,1x,a)')
     +  'Electron counting (NCTAIM)           :','ELCOUNT'

      if(idofr.eq.1) then
        write(cval,'(i0,a)') icufr,' fragments (see below)'
        write(*,'(2x,a,1x,a)')
     +    'Fragment analysis                    :',trim(cval)
      end if
      if(idoat.eq.1) then
        write(cval,'(i0,a)') icuat,' atoms selected'
        write(*,'(2x,a,1x,a)')
     +    'Atom selection                       :',trim(cval)
      end if

      cbuf=' '
      if(iqchem.eq.1) cbuf=trim(cbuf)//' QCHEM;'
      if(imokit.eq.1) cbuf=trim(cbuf)//' MOKIT;'
      if(iwfn.eq.1) cbuf=trim(cbuf)//' WFN;'
      cbuf=adjustl(cbuf)
      nlen=len_trim(cbuf)
      if(nlen.gt.0.and.cbuf(nlen:nlen).eq.';') cbuf(nlen:nlen)=' '
      if(len_trim(cbuf).gt.0) write(*,'(2x,a,1x,a)')
     +  'Wavefunction source                  :',trim(cbuf)

      if(icube.eq.1) write(*,'(2x,a,1x,a)')
     +  'Cube files                           :','CUBE (see below)'
      if(iaccur.eq.1) write(*,'(2x,a,1x,a)')
     +  'Output precision                     :','FULLPRECISION'
      if(inopop.eq.1) write(*,'(2x,a,1x,a)')
     +  'Population output                    :','NOPOPU (suppressed)'
      if(ifinegrid.eq.1) write(*,'(2x,a,1x,a)')
     +  'Fine angular grid                    :',
     +  'FINEGRID (974 points, one-electron integrals)'

!! ----------------------------------------------------------------- !!
!! # ENPART -- only if ENPART itself is active                       !!
!! ----------------------------------------------------------------- !!
      if(ienpart.eq.1) then
        write(*,*)
        chdr='# ENPART'
        write(*,'(2x,a)') trim(chdr)
        write(*,'(2x,a)') repeat('-',len_trim(chdr))

        if(ihf.eq.1) then
          cval='HF'
        else if(icas.eq.1) then
          cval='CASSCF'
        else if(icisd.eq.1) then
          cval='CISD'
        else if(id_xcfunc.eq.402) then
          cval='B3LYP (hybrid)'
        else if(id_xfunc.eq.106.and.id_cfunc.eq.132) then
          cval='BP86 (GGA)'
        else if(id_xfunc.eq.1.and.id_cfunc.eq.0) then
          cval='LDA'
        else if(id_xcfunc.ne.0.or.id_xfunc.ne.0.or.id_cfunc.ne.0) then
          write(cval,'(a,3(1x,i0))') 'custom libxc ids',
     +      id_xcfunc,id_xfunc,id_cfunc
        else
          cval='unspecified'
        end if
        write(*,'(2x,a,1x,a)') 'Functional              :',trim(cval)

        write(cval,'(a,i0)') 'THREBOD = ',ithrebod
        write(*,'(2x,a,1x,a)')
     +    'Atom-pair skip threshold:',trim(cval)

        if(iigrid.eq.1) write(*,'(2x,a,1x,a)')
     +    'Two-electron grid       :',
     +    'MOD-GRIDTWOEL (user-selected integration grid)'

        cbuf=' '
        if(iexact.eq.1) cbuf=trim(cbuf)//' EXACT;'
        if(ihomo.eq.1) cbuf=trim(cbuf)//' HOMO;'
        if(idek.eq.1) cbuf=trim(cbuf)//' DEKIN;'
        if(iionic.eq.1) cbuf=trim(cbuf)//' IONIC;'
        if(ianalytical.eq.1) cbuf=trim(cbuf)//' ANALYTIC;'
        if(iecorr.eq.1) cbuf=trim(cbuf)//' CORRELATION;'
        cbuf=adjustl(cbuf)
        nlen=len_trim(cbuf)
        if(nlen.gt.0.and.cbuf(nlen:nlen).eq.';') cbuf(nlen:nlen)=' '
        if(len_trim(cbuf).gt.0) write(*,'(2x,a,1x,a)')
     +    'Extra options           :',trim(cbuf)

!! MOD-GRIDTWOEL's own # GRID settings -- iigrid is shared with EDAIQA's !!
!! own MOD-GRIDTWOEL check below; if both ENPART and EDAIQA are active   !!
!! in the same run, whichever parses last in read_input.f wins here too  !!
        if(iigrid.eq.1) then
          write(*,*)
          chdr='# GRID'
          write(*,'(2x,a)') trim(chdr)
          write(*,'(2x,a)') repeat('-',len_trim(chdr))
          write(cval,'(i0)') nrad22
          write(*,'(2x,a,1x,a)') 'Radial points              :',
     +      trim(cval)
          write(cval,'(i0)') nang22
          write(*,'(2x,a,1x,a)') 'Angular points             :',
     +      trim(cval)
          write(cval,'(f6.3,a,f6.3)') phb12,' / ',phb22
          write(*,'(2x,a,1x,a)') 'Rotation angles (phb1/phb2):',
     +      trim(cval)
        end if
      end if

!! ----------------------------------------------------------------- !!
!! # EDAIQA -- only if EDAIQA is active                              !!
!! ----------------------------------------------------------------- !!
      if(iedaiqa.eq.1) then
        write(*,*)
        chdr='# EDAIQA'
        write(*,'(2x,a)') trim(chdr)
        write(*,'(2x,a)') repeat('-',len_trim(chdr))
        if(iflip.eq.1) write(*,'(2x,a,1x,a)')
     +    'Flip alpha/beta spins:','FLIPSPIN'
        if(iigrid.eq.1) write(*,'(2x,a,1x,a)')
     +    'Two-electron grid    :',
     +    'MOD-GRIDTWOEL (user-selected integration grid)'
      end if

!! ----------------------------------------------------------------- !!
!! # TOPOLOGY -- only if TOPOLOGY is active                          !!
!! ----------------------------------------------------------------- !!
      if(itop.eq.1) then
        write(*,*)
        chdr='# TOPOLOGY'
        write(*,'(2x,a)') trim(chdr)
        write(*,'(2x,a)') repeat('-',len_trim(chdr))
        if(ietop.eq.1) cval='Exchange'
        if(ietop.eq.2) cval='Correlation'
        if(ietop.eq.3) cval='Exchange-Correlation'
        if(ietop.eq.9) cval='Density'
        write(*,'(2x,a,1x,a)') 'Energy component:',trim(cval)
        if(ipairs.gt.0) then
          write(cval,'(i0,a)') ipairs,' atom pairs'
        else
          cval='entire molecule'
        end if
        write(*,'(2x,a,1x,a)') 'Atom pairs      :',trim(cval)
      end if

!! ----------------------------------------------------------------- !!
!! # QTAIM -- only if QTAIM is active                                !!
!! ----------------------------------------------------------------- !!
      if(iqtaim.ge.1) then
        write(*,*)
        chdr='# QTAIM'
        write(*,'(2x,a)') trim(chdr)
        write(*,'(2x,a)') repeat('-',len_trim(chdr))
        write(cval,'(i0)') istep
        write(*,'(2x,a,1x,a)') 'Convergence step            :',
     +    trim(cval)
        write(cval,'(i0)') inna
        write(*,'(2x,a,1x,a)') 'NNA (non-nuclear attractors):',
     +    trim(cval)
        write(cval,'(i0)') imaxdist
        write(*,'(2x,a,1x,a)') 'Max basin distance          :',
     +    trim(cval)
        write(cval,'(i0)') iscreening
        write(*,'(2x,a,1x,a)') 'Screening                   :',
     +    trim(cval)
        write(cval,'(i0)') ipath
        write(*,'(2x,a,1x,a)') 'Gradient path               :',
     +    trim(cval)
      end if

!! ----------------------------------------------------------------- !!
!! # CUBE -- only if CUBE is active                                  !!
!! ----------------------------------------------------------------- !!
      if(icube.eq.1) then
        write(*,*)
        chdr='# CUBE'
        write(*,'(2x,a)') trim(chdr)
        write(*,'(2x,a)') repeat('-',len_trim(chdr))
        write(cval,'(i0)') jcubthr
        write(*,'(2x,a,1x,a)') 'Max occupation threshold:',trim(cval)
        write(cval,'(i0)') kcubthr
        write(*,'(2x,a,1x,a)') 'Min occupation threshold:',trim(cval)
        write(cval,'(f6.3)') cubespacing
        write(*,'(2x,a,1x,a)') 'Grid spacing (bohr)     :',trim(adjustl(cval))
        write(cval,'(f6.3)') cuberadscale
        write(*,'(2x,a,1x,a)') 'Radius scale            :',trim(adjustl(cval))
      end if

!! ----------------------------------------------------------------- !!
!! # OSLO -- only if OSLO is active                                  !!
!! ----------------------------------------------------------------- !!
      if(ioslo.eq.1) then
        write(*,*)
        chdr='# OSLO'
        write(*,'(2x,a)') trim(chdr)
        write(*,'(2x,a)') repeat('-',len_trim(chdr))
        if(ilow2.eq.1) cval='Mulliken'
        if(ilow2.eq.2) cval='Lowdin'
        if(ilow2.eq.3) cval='Lowdin-Davidson'
        if(ilow2.eq.6) cval='NAO-basis'
        if(ilow2.eq.0) cval='TFVC/real-space (default)'
        write(*,'(2x,a,1x,a)') 'Overlap matrix          :',trim(cval)
        write(cval,'(i0)') ifolitol
        write(*,'(2x,a,1x,a)') 'FOLI tolerance          :',trim(cval)
        write(cval,'(i0)') ibranch
        write(*,'(2x,a,1x,a)') 'Branch iteration        :',trim(cval)
        if(ioslofchk.eq.2) then
          cval='yes'
        else
          cval='no'
        end if
        write(*,'(2x,a,1x,a)') 'Print non-orthogonalized:',trim(cval)
      end if

!! ----------------------------------------------------------------- !!
!! # FRAGMENTS -- only if DOFRAGS is active; replaces the raw        !!
!! 'Fragment: N' + bare atom-index dump read_input.f used to print   !!
!! ----------------------------------------------------------------- !!
      if(idofr.eq.1) then
        write(*,*)
        write(chdr,'(a,i0,a)') '# FRAGMENTS  (',icufr,' fragments)'
        write(*,'(2x,a)') trim(chdr)
        write(*,'(2x,a)') repeat('-',len_trim(chdr))
        do i=1,icufr
          write(*,'(2x,a,i3,a,20i4)') 'Fragment',i,' :',
     +      (ifrlist(k,i),k=1,nfrlist(i))
        end do
      end if

!! ----------------------------------------------------------------- !!
!! # ATOMS -- only if DOATOMS is active                              !!
!! ----------------------------------------------------------------- !!
      if(idoat.eq.1) then
        write(*,*)
        write(chdr,'(a,i0,a)') '# ATOMS  (',icuat,' selected)'
        write(*,'(2x,a)') trim(chdr)
        write(*,'(2x,a)') repeat('-',len_trim(chdr))
        write(*,'(2x,20i4)') (iatlist(i),i=1,icuat)
      end if

!! ----------------------------------------------------------------- !!
!! # DM -- only if correlated-WF density input is active             !!
!! ----------------------------------------------------------------- !!
      if(icorr.ne.0) then
        write(*,*)
        chdr='# DM'
        write(*,'(2x,a)') trim(chdr)
        write(*,'(2x,a)') repeat('-',len_trim(chdr))
        write(cval,'(i0,a)') icorr,'-RDM'
        write(*,'(2x,a,1x,a)') 'DM level:',trim(cval)
        if(ipyscf.eq.1) then
          cval='pySCF'
        else if(iorca.eq.1) then
          cval='ORCA'
        else
          cval='DMRG'
        end if
        write(*,'(2x,a,1x,a)') 'Source  :',trim(cval)
      end if

      END SUBROUTINE print_input_summary

!! ***** !!

!! ********************************************************************* !!
!! subroutine: VPRINT                                                    !!
!! purpose: unbordered print of a per-atom vector (or 2-column matrix,   !!
!!   via jdim) -- atom index, element symbol, then jdim value(s). No     !!
!!   6-column chunking (jdim is always 1 or 2 at call sites).            !!
!! arguments:                                                            !!
!!   H    (in) -- data (ndim,jdim)                                       !!
!!   N    (in) -- number of atoms actually printed                       !!
!!   ndim (in) -- H's declared leading dimension                         !!
!!   jdim (in) -- number of value columns (1 or 2)                       !!
!! author:                                                               !!
!! ********************************************************************* !!
      SUBROUTINE VPRINT(H,N,ndim,jdim)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim,jdim
      common /printout/iaccur
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      DIMENSION H(NDIM,jdim)

      Dimension mend(92)
      data mend/4H  H ,4H He ,4H Li ,4H Be ,4H  B ,4H  C ,4H  N ,4H  O ,
     $ 4H  F ,4H Ne ,4H Na ,4H Mg ,4H Al ,4H Si ,4H  P ,4H  S ,4H Cl ,
     $4H Ar ,4H  K ,4H Ca ,4H Sc ,4H Ti ,4H  V ,4H Cr ,4H Mn ,4H Fe ,
     $4H Co ,4H Ni ,4H Cu ,4H Zn ,4H Ga ,4H Ge ,4H As ,4H Se ,4H Br ,
     $4H Kr ,4H Rb ,4H Sr ,4H  Y ,4H Zr ,4H Nb ,4H Mo ,4H Tc ,4H Ru ,
     $4H Rh ,4H Pd ,4H Ag ,4H Cd ,4H In ,4H Sn ,4H Sb ,4H Te ,4H  I ,
     $4H Xe ,4H Cs ,4H Ba ,4H La ,4H Ce ,4H Pr ,4H Nd ,4H Pm ,4H Sn ,
     $4H Eu ,4H Gd ,4H Tb ,4H Dy ,4H Ho ,4H Er ,4H Tm ,4H Yb ,4H Lu ,
     $4H Hf ,4H Ta ,4H  W ,4H Re ,4H Os ,4H Ir ,4H Pt ,4H Au ,4H Hg ,
     $4H Tl ,4H Pb ,4H Bi ,4H Po ,4H At ,4H Rn ,4H Fr ,4H Ra ,4H Ac ,
     $4H Th ,4H Pa ,4H  U   /

   62 FORMAT(1X,I3,A4,6F12.6)
   63 FORMAT(1X,I3,A4,6F20.13)
      DO 2 I=1,N
        if(iaccur.eq.0) then
          PRINT 62,I,mend(iznuc(i)),(H(I,J),J=1,jdim)
        else
          PRINT 63,I,mend(iznuc(i)),(H(I,J),J=1,jdim)
        end if
   2  CONTINUE
      RETURN
      END

!! ***** !!

!! ********************************************************************* !!
!! subroutine: MPRINT_NLOP                                               !!
!! purpose: unbordered print of a per-atom 3-column (X/Y/Z) matrix, with !!
!!   its own column header line -- used for dipole/charge-transfer       !!
!!   component tables in enpart.f.                                       !!
!! arguments:                                                            !!
!!   H    (in) -- data (ndim,3)                                          !!
!!   N    (in) -- number of atoms actually printed                       !!
!!   ndim (in) -- H's declared leading dimension                         !!
!! author:                                                               !!
!! ********************************************************************* !!
      SUBROUTINE MPRINT_NLOP(H,N,ndim)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim
      common /printout/iaccur
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      DIMENSION H(NDIM,3)

      Dimension mend(92)
      data mend/4H  H ,4H He ,4H Li ,4H Be ,4H  B ,4H  C ,4H  N ,4H  O ,
     $ 4H  F ,4H Ne ,4H Na ,4H Mg ,4H Al ,4H Si ,4H  P ,4H  S ,4H Cl ,
     $4H Ar ,4H  K ,4H Ca ,4H Sc ,4H Ti ,4H  V ,4H Cr ,4H Mn ,4H Fe ,
     $4H Co ,4H Ni ,4H Cu ,4H Zn ,4H Ga ,4H Ge ,4H As ,4H Se ,4H Br ,
     $4H Kr ,4H Rb ,4H Sr ,4H  Y ,4H Zr ,4H Nb ,4H Mo ,4H Tc ,4H Ru ,
     $4H Rh ,4H Pd ,4H Ag ,4H Cd ,4H In ,4H Sn ,4H Sb ,4H Te ,4H  I ,
     $4H Xe ,4H Cs ,4H Ba ,4H La ,4H Ce ,4H Pr ,4H Nd ,4H Pm ,4H Sn ,
     $4H Eu ,4H Gd ,4H Tb ,4H Dy ,4H Ho ,4H Er ,4H Tm ,4H Yb ,4H Lu ,
     $4H Hf ,4H Ta ,4H  W ,4H Re ,4H Os ,4H Ir ,4H Pt ,4H Au ,4H Hg ,
     $4H Tl ,4H Pb ,4H Bi ,4H Po ,4H At ,4H Rn ,4H Fr ,4H Ra ,4H Ac ,
     $4H Th ,4H Pa ,4H  U   /

      NMIN=1
      NNMAX=3
      PRINT *,'_________________X___________________Y__________________Z__________'
      PRINT *,' '
      DO I=1,N
        if(iaccur.eq.0) then
          PRINT 62,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
        else
          PRINT 63,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
        end if
      end do
      RETURN
   62 FORMAT(1X,I3,A4,6F12.6)
   63 FORMAT(1X,I3,A4,3F20.13)
      END

!! ********************************************************************* !!
!! subroutine: MPRINT                                                    !!
!! purpose: unbordered N x N matrix print, chunked 6 columns at a time   !!
!!   (loops back to label 1 for each further chunk of columns until      !!
!!   NMIN exceeds N). No border -- see MPRINT2 for the bordered twin     !!
!!   used throughout ENPART.                                             !!
!! arguments:                                                            !!
!!   H    (in) -- data matrix (ndim,ndim)                                !!
!!   N    (in) -- number of atoms actually printed                       !!
!!   ndim (in) -- H's declared leading dimension                         !!
!! author:                                                               !!
!! ********************************************************************* !!
      SUBROUTINE MPRINT(H,N,ndim)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim
      common /printout/iaccur
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      DIMENSION H(NDIM,NDIM)

      Dimension mend(92)
      data mend/4H  H ,4H He ,4H Li ,4H Be ,4H  B ,4H  C ,4H  N ,4H  O ,
     $ 4H  F ,4H Ne ,4H Na ,4H Mg ,4H Al ,4H Si ,4H  P ,4H  S ,4H Cl ,
     $4H Ar ,4H  K ,4H Ca ,4H Sc ,4H Ti ,4H  V ,4H Cr ,4H Mn ,4H Fe ,
     $4H Co ,4H Ni ,4H Cu ,4H Zn ,4H Ga ,4H Ge ,4H As ,4H Se ,4H Br ,
     $4H Kr ,4H Rb ,4H Sr ,4H  Y ,4H Zr ,4H Nb ,4H Mo ,4H Tc ,4H Ru ,
     $4H Rh ,4H Pd ,4H Ag ,4H Cd ,4H In ,4H Sn ,4H Sb ,4H Te ,4H  I ,
     $4H Xe ,4H Cs ,4H Ba ,4H La ,4H Ce ,4H Pr ,4H Nd ,4H Pm ,4H Sn ,
     $4H Eu ,4H Gd ,4H Tb ,4H Dy ,4H Ho ,4H Er ,4H Tm ,4H Yb ,4H Lu ,
     $4H Hf ,4H Ta ,4H  W ,4H Re ,4H Os ,4H Ir ,4H Pt ,4H Au ,4H Hg ,
     $4H Tl ,4H Pb ,4H Bi ,4H Po ,4H At ,4H Rn ,4H Fr ,4H Ra ,4H Ac ,
     $4H Th ,4H Pa ,4H  U   /

      K=6
      NMIN=1
      NNMAX=MIN0(N,K)
   62 FORMAT(1X,I3,A4,6F12.6)
   63 FORMAT(1X,I3,A4,6F20.13)
   1  if(iaccur.eq.0) then
        PRINT 60, (I, mend(iznuc(i)),I=NMIN,NNMAX)
      else
        PRINT 61, (I, mend(iznuc(i)),I=NMIN,NNMAX)
      end if
   60 FORMAT(10X,6(2X,I3,A4,3X))
   61 FORMAT(10X,6(6X,I3,A4,7X))
      PRINT 64
   64 FORMAT(1X)
      DO 2 I=1,N
        if(iaccur.eq.0) then
          PRINT 62,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
        else
          PRINT 63,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
        end if
   2  CONTINUE
      NMIN=NMIN+6
      K=K+6
      NNMAX=MIN0(N,K)
      IF(NNMAX.GE.NMIN) GOTO 71
      RETURN
   71 PRINT 66
   66 FORMAT(1X///)
      GO TO 1
      END

!! ***** !!

!! ********************************************************************* !!
!! subroutine: MPRINT2                                                   !!
!! purpose: bordered N x N matrix print, chunked 6 columns at a time --  !!
!!   same chunk/loop-back structure as MPRINT, plus an 80-column rule    !!
!!   before/after the header and each row block. The style used         !!
!!   throughout ENPART for atom-pair matrices.                           !!
!! arguments:                                                            !!
!!   H    (in) -- data matrix (ndim,ndim)                                !!
!!   N    (in) -- number of atoms actually printed                       !!
!!   ndim (in) -- H's declared leading dimension                         !!
!! author:                                                               !!
!! ********************************************************************* !!
      SUBROUTINE MPRINT2(H,N,ndim)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim
      character*100 line
      common /printout/iaccur
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      DIMENSION H(NDIM,NDIM)

      Dimension mend(92)
      data mend/4H  H ,4H He ,4H Li ,4H Be ,4H  B ,4H  C ,4H  N ,4H  O ,
     $ 4H  F ,4H Ne ,4H Na ,4H Mg ,4H Al ,4H Si ,4H  P ,4H  S ,4H Cl ,
     $4H Ar ,4H  K ,4H Ca ,4H Sc ,4H Ti ,4H  V ,4H Cr ,4H Mn ,4H Fe ,
     $4H Co ,4H Ni ,4H Cu ,4H Zn ,4H Ga ,4H Ge ,4H As ,4H Se ,4H Br ,
     $4H Kr ,4H Rb ,4H Sr ,4H  Y ,4H Zr ,4H Nb ,4H Mo ,4H Tc ,4H Ru ,
     $4H Rh ,4H Pd ,4H Ag ,4H Cd ,4H In ,4H Sn ,4H Sb ,4H Te ,4H  I ,
     $4H Xe ,4H Cs ,4H Ba ,4H La ,4H Ce ,4H Pr ,4H Nd ,4H Pm ,4H Sn ,
     $4H Eu ,4H Gd ,4H Tb ,4H Dy ,4H Ho ,4H Er ,4H Tm ,4H Yb ,4H Lu ,
     $4H Hf ,4H Ta ,4H  W ,4H Re ,4H Os ,4H Ir ,4H Pt ,4H Au ,4H Hg ,
     $4H Tl ,4H Pb ,4H Bi ,4H Po ,4H At ,4H Rn ,4H Fr ,4H Ra ,4H Ac ,
     $4H Th ,4H Pa ,4H  U   /

      line="--------------------------------------------------------------------------------"

      K=6
      NMIN=1
      NNMAX=MIN0(N,K)
   62 FORMAT(2X,I3,A4,6F12.6)
   63 FORMAT(2X,I3,A4,6F20.13)
   1  if(iaccur.eq.0) then
        PRINT 666,line
        PRINT 60, (I, mend(iznuc(i)),I=NMIN,NNMAX)
      else
        PRINT 666,line
        PRINT 61, (I, mend(iznuc(i)),I=NMIN,NNMAX)
      end if
   60 FORMAT(10X,6(2X,I3,A4,3X))
   61 FORMAT(10X,6(6X,I3,A4,7X))
      PRINT 666,line
      DO 2 I=1,N
        if(iaccur.eq.0) then
          PRINT 62,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
        else
          PRINT 63,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
        end if
   2  CONTINUE
      PRINT 666,line
      NMIN=NMIN+6
      K=K+6
      NNMAX=MIN0(N,K)
      IF(NNMAX.GE.NMIN) GOTO 71
      RETURN
   71 PRINT 66
   66 FORMAT(1X)
      GO TO 1

!! format for the border lines !!
  666 FORMAT(2x,a100)
      END

!! ***** !!

!! ********************************************************************* !!
!! subroutine: mprintnoat                                                !!
!! purpose: bordered, chunked (6 columns per block) numeric matrix       !!
!! print for a fragment x fragment matrix -- same visual style as        !!
!! MPRINT2, minus the atom-symbol column (fragments don't have one).     !!
!! The column header is just the fragment number, no "Frag" label, so a  !!
!! block of numbers copy-pastes cleanly. Used only by group_by_frag_mat  !!
!! arguments:                                                            !!
!! H (in) -- data matrix, M rows x N columns                             !!
!! M,N (in) -- rows/columns actually used                                !!
!! mdim,ndim (in) -- H's declared dimensions                             !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      SUBROUTINE MPRINTNOAT(H,M,N,mdim,ndim)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim,mdim
      DIMENSION H(MDIM,NDIM)
      common /printout/iaccur
      character*100 line
      character*40 hdrfmt

      line="--------------------------------------------------------------------------------"

      K=6
      NMIN=1
      NNMAX=MIN0(N,K)
   62 FORMAT(2X,I3,6X,6F12.6)
   63 FORMAT(2X,I3,6X,6F20.13)
    1 continue
      write(*,666) line
!! repeat count built at runtime from the actual chunk size (NNMAX-NMIN+1,  !!
!! <=6), not a literal '6(...)', to stay correct however many columns the   !!
!! last chunk actually has. Each fragment number is right-justified in a    !!
!! field exactly as wide as its F12.6/F20.13 data column (I12/I20), so it   !!
!! lands flush with that column's last digit instead of drifting left.      !!
      if(iaccur.eq.0) then
        write(hdrfmt,'(a,i0,a)') '(11X,',NNMAX-NMIN+1,'I12)'
      else
        write(hdrfmt,'(a,i0,a)') '(11X,',NNMAX-NMIN+1,'I20)'
      end if
      write(*,hdrfmt) (I,I=NMIN,NNMAX)
      write(*,666) line
      DO 2 I=1,M
      if(iaccur.eq.0) then
      PRINT 62,I,(H(I,J),J=NMIN,NNMAX)
      else
      PRINT 63,I,(H(I,J),J=NMIN,NNMAX)
      end if
   2  CONTINUE
      write(*,666) line
      NMIN=NMIN+6
      K=K+6
      NNMAX=MIN0(N,K)
      IF(NNMAX.GE.NMIN) GOTO 71
      RETURN
   71 write(*,*)
      GO TO 1

!! FORMAT FOR THE BORDER LINE !!
  666 FORMAT(2x,a100)
      END

!! ***** !!

!! ********************************************************************* !!
!! subroutine: group_by_frag_mat                                         !!
!! purpose: sums a per-atom matrix A into a per-fragment matrix B (using !!
!! /frlist/'s atom-to-fragment map) and prints it via mprintnoat,        !!
!! MPRINT2-style (bordered, chunked, per-fragment column header).        !!
!! arguments:                                                            !!
!! ilog (in) -- 0: sum the full matrix, 1: lower-triangular only         !!
!!              (symmetric quantity, e.g. bond order)                    !!
!! line (in) -- title, printed via print_box (leading/trailing blanks    !!
!!              in the caller's string are stripped)                     !!
!! A    (in) -- the (maxat,maxat) per-atom matrix to sum and print       !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine group_by_frag_mat(ilog,line,A)
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /printout/iaccur
      dimension A(maxat,maxat),B(maxat,maxat)
      character*80 line
      integer ilog

!! ilog=0: full matrix, ilog=1: lower triangular (symmetric quantity) !!
      do i=1,nat
       do j=1,nat
        B(i,j)=0.0d0
       end do
      end do

      nnat=nat
      do i=1,nat
       if(ilog.eq.1) nnat=i
       do j=1,nnat
        if(jfrlist(i).ne.0.and.jfrlist(j).ne.0) B(jfrlist(i),jfrlist(j))=B(jfrlist(i),jfrlist(j))+a(i,j)
       end do
      end do

      x=0.0d0
      if(ilog.eq.1) then
       do i=1,icufr
        x=x+b(i,i)
        do j=1,i-1
         B(i,j)=B(i,j)+B(j,i)
         x=x+b(i,j)
         B(j,i)=B(i,j)
        end do
       end do
      else
      do i=1,icufr
       do j=1,icufr
        x=x+b(i,j)
       end do
      end do
      end if

      call print_box(trim(adjustl(line)))
      call mprintnoat(B,icufr,icufr,maxat,maxat)
      if(iaccur.eq.0) then
        write(*,163) x
      else
        write(*,164) x
      end if

      return

  163 format(2x,'   Total:',f12.6)
  164 format(2x,'   Total:',f20.13)

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: group_by_frag_vec                                         !!
!! purpose: sums a per-atom (maxat,ndim) table A into a per-fragment      !!
!! table B (using /frlist/'s atom-to-fragment map) and prints it,        !!
!! bordered to its own actual width -- unlike group_by_frag_mat's        !!
!! fragment x fragment matrix, this table's column count never grows     !!
!! with the number of fragments (always 1 or 2: 3D-space[, Mulliken]),   !!
!! so a fixed-width MPRINT2-style border would just run off past the     !!
!! last number for no reason.                                            !!
!! arguments:                                                            !!
!! ndim (in) -- number of value columns (1 or 2: 3D-space[, Mulliken])   !!
!! line (in) -- title, printed via print_box (leading/trailing blanks    !!
!!              in the caller's string are stripped)                     !!
!! A    (in) -- the (maxat,ndim) per-atom table to sum and print         !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine group_by_frag_vec(ndim,line,A)
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      character*80 line
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /printout/iaccur
      dimension A(maxat,ndim),B(maxat,ndim)
      dimension x(ndim)
      integer colw,totw

      do i=1,nat
       do j=1,ndim
        B(i,j)=0.0d0
       end do
      end do

      do j=1,ndim
       do i=1,nat
        if(jfrlist(i).ne.0) B(jfrlist(i),j)=B(jfrlist(i),j)+a(i,j)
       end do
      end do

      do j=1,ndim
       x(j)=0.0d0
       do i=1,icufr
        x(j)=x(j)+b(i,j)
       end do
      end do

!! 11-column fragment-index margin (2X+I3+6X), matching mprintnoat's own !!
!! -- total width = margin + one value column (12/20-wide) per ndim      !!
      if(iaccur.eq.0) then
        colw=12
      else
        colw=20
      end if
      totw=11+colw*ndim

      call print_box(trim(adjustl(line)))
      if(iaccur.eq.0) then
        if(ndim.eq.2) then
          write(*,'(a11,2a12)') 'Fragment','3D-space','Mulliken'
        else
          write(*,'(a11,a12)') 'Fragment','3D-space'
        end if
      else
        if(ndim.eq.2) then
          write(*,'(a11,2a20)') 'Fragment','3D-space','Mulliken'
        else
          write(*,'(a11,a20)') 'Fragment','3D-space'
        end if
      end if
      write(*,'(2x,a)') repeat('-',totw-2)
      do i=1,icufr
        if(iaccur.eq.0) then
          write(*,'(2x,i3,6x,2f12.6)') i,(B(i,j),j=1,ndim)
        else
          write(*,'(2x,i3,6x,2f20.13)') i,(B(i,j),j=1,ndim)
        end if
      end do
      write(*,'(2x,a)') repeat('-',totw-2)
      if(iaccur.eq.0) then
        write(*,63) (x(j),j=1,ndim)
      else
        write(*,64) (x(j),j=1,ndim)
      end if

      return

   63 format(2x,'   Total:',2f12.6)
   64 format(2x,'   Total:',2f20.13)

      end

!! ********************************************************************* !!
!! subroutine: print_int                                                 !!
!! purpose: writes one AIMPAC/PROAIMV-format ".int" file per atom (plus  !!
!!   a ".files" index) holding that atom's overlap matrix in the MO/NO   !!
!!   basis -- consumed by the external FCALC program to compute atomic   !!
!!   properties from APOST-3D's atomic partition. Opt-in via # METHOD /  !!
!!   DOINT. Also runs a final MO-orthogonality sanity check (sum of all  !!
!!   atomic overlap matrices vs. the identity) once every atom is done.  !!
!! arguments:                                                            !!
!!   nbas,nat0 (in) -- sat's declared dimensions (basis functions, atoms)!!
!!   sat       (in) -- per-atom atomic-orbital overlap matrix            !!
!!   name      (in) -- job name, used as the ".int"/".files" basename    !!
!! author:                                                               !!
!! ********************************************************************* !!
      subroutine print_int(nbas,nat0,sat,name)
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /atlist/iatlist(maxat),icuat
      common/cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /coord/ coord0(3,maxat),zn(maxat),iznuc(maxat)
      dimension sat(nbas,nbas,nat0)
      character name*60
      character*2 mend(92)
      data mend/' H','He','Li','Be',' B',' C',' N',' O',
     $ ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl',
     $ 'Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe',
     $ 'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
     $ 'Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     $ 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I',
     $ 'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sn',
     $ 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
     $ 'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     $ 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     $ 'Th','Pa',' U'  /

      character nameaim*55,charnu*2,charnu2*2,charnu1
      character charnu3*3,ext*4
      character*80 line

      allocatable c3(:,:),c2(:,:),csave(:,:),scr(:,:)

!! iopt flags used by this routine !!
      ihirsh = Iopt(6)
      imulli=Iopt(5)
      icorr=Iopt(26)
      iqtaim =Iopt(16)

      inato=0
      if(icorr.eq.1.or.icas.eq.1.or.icisd.eq.1) inato=1

!! rescan the .fchk (unit 15, already fully read by input()) for the    !!
!! independent-function count -- same rewind+rescan pattern as          !!
!! input2.f's readchar/readint.                                         !!
      rewind(15)
 998  read(15,'(a80)') line
      if(index(line,"Number of independant functions").ne.0) then
        read(line(54:61),'(i8)') ndim
      else if(index(line,"Number of independent functions").ne.0) then
        read(line(54:61),'(i8)') ndim
      else
        go to 998
      end if

      allocate( c3(2*igr,2*igr))
      allocate( c2(igr,igr),csave(igr,igr),scr(igr,igr))
      scr=0.0d0

      do i=1,igr
        do j=1,igr
          scr(i,j)=c(i,j)
        end do
      end do

      if(inato.eq.1) then
        if(icas.eq.1) ndim=norb
        write(*,'(2x,a,1x,i0,1x,a)') 'Using',ndim,'natural orbitals'
        do i=1,igr
          do j=1,ndim
            scr(i,j)=c_no(i,j)
          end do
        end do
      end if

      do i=1,igr
        do j=1,ndim
          csave(i,j)=0.0d0
        end do
      end do

!! building file !!
      naim=88
      naim3=89
      l=len_trim(name)
      if(imulli.eq.1) then
        nameaim=name(1:l)//"mul.files"
        ext="mul_"
      else if(imulli.ge.2) then
        nameaim=name(1:l)//"low.files"
        ext="low_"
      else if(iqtaim.eq.1) then
        nameaim=name(1:l)//"aim.files"
        ext="aim_"
      else if(ihirsh.eq.0) then
        nameaim=name(1:l)//"fuz.files"
        ext="fuz_"
      else if(ihirsh.eq.1) then
        nameaim=name(1:l)//"hir.files"
        ext="hir_"
      else if(ihirsh.eq.2) then
        nameaim=name(1:l)//"ihi.files"
        ext="ihi_"
      end if
      open(file=nameaim,unit=naim3,status="unknown")
      rewind(naim3)
      do jjat=1,icuat
        if(icuat.ne.nat) then
          jat=iatlist(jjat)
        else
          jat=jjat
        end if
        read(mend(iznuc(jat)),'(A2)')charnu
        charnu=adjustl(charnu)
        l1=len_trim(name)
!! atom index formatted 1/2/3 digits wide -- assumes up to 999 atoms !!
        if(jat.lt.10) then
          write(charnu1,'(i1)')jat
          nameaim=name(1:l1)//ext//trim(charnu)//charnu1
        else if(jat.lt.100) then
          write(charnu2,'(i2)')jat
          nameaim=name(1:l1)//ext//trim(charnu)//charnu2
        else
          write(charnu3,'(i3)')jat
          nameaim=name(1:l1)//ext//trim(charnu)//charnu3
        end if
        j=len(nameaim)
        do i=1,j
          if(nameaim(i:i).eq.' ') then
            l=i-1
            go to 10
          end if
        end do
  10   continue
        nameaim=nameaim(1:l)//".int"
        write(naim3,'(a40)') nameaim
        open(file=nameaim,unit=naim,status="unknown")
        rewind(naim)

!! the following .int-file content is read by the external FCALC       !!
!! program -- format is AIMPAC/PROAIMV-compatible, do not restyle. AE   !!
!! and BK are placeholders (see header) kept only for format            !!
!! compatibility, not computed from the real SCF energy.                !!
        BK=1.0d0
        AE=-1.0d0
        write(naim,'(a20)') 'Created by APOST3D '
        if(iopt(5).eq.1) then
          write(naim,*) 'Using Mulliken atomic definition '
        else if(iopt(5).ge.2) then
          write(naim,*) 'Using Lowdin atomic definition '
        else if(iopt(6).eq.1) then
          write(naim,*) 'Using Hirshfeld atomic definition '
        else if(iopt(6).eq.2) then
          write(naim,*) 'Using Hirshfeld-Iterative atomic definition '
        else if(iopt(14).eq.1.and.iopt(31).eq.0) then
          write(naim,*) 'Using Becke-rho atomic definition '
        else if(iopt(14).eq.1.and.iopt(31).eq.1) then
          write(naim,*) 'Using TFVC atomic definition '
        else if(iopt(16).eq.1) then
          write(naim,*) 'Using QTAIM atomic definition '
        else
          write(naim,*) 'Using Becke atomic definition '
        end if
        if(iopt(5).eq.0) write(naim,*) 'Stiffness parameter k = ',iopt(25)
        if(inato.eq.1)  then
          write(naim,*) 'Correlated wave function '
          if(icorr.eq.0) then
            write(naim,*) '*Warning*, no RDM1 provided.'
            write(naim,*) 'Can not do alpha and beta populations separatedly.'
          end if
        else
          write(naim,*) 'Single-determinant wave function '
          if(kop.eq.1) write(naim,*) 'Unrestricted wave function '
        end if

        write(naim,'(a30,F20.11)')' MOLECULAR SCF ENERGY (AU)  = ',AE
        write(naim,*)''
        write(naim,'(A25,A4,i5)')' INTEGRATION IS OVER ATOM',mend(iznuc(jat)),
     1  jat
      !write(naim2) mend(iznuc(jat)),jat
        write(naim,'(a27)') ' RESULTS OF THE INTEGRATION'
        write(naim,'(a17,e21.14,a14,e21.14)')'              N  ',
     1   qat(jjat,1),'    NET CHARGE',-qat(jjat,1)+iznuc(jjat)
        !write(naim2) qat(jjat,1)
        write(naim,*)'             G'
        write(naim,'(a17,e21.14,a17,e21.14)')'              K  ',BK,
     1  '        E(ATOM)  ',AE
        write(naim,'(a17,e21.14)')'              L  ',0.0d0
        write(naim,*)''
        write(naim,'(a35)')'          The Atomic Overlap Matrix'
        write(naim,*)''

        if(inato.eq.1) then
          write(naim,'(a36)')'  Correlated  Wavefunction'
          write(naim,*)''
        else if(kop.eq.1) then
          write(naim,'(a36)')'Unrestricted  Wavefunction'
          write(naim,*)''
        else
          write(naim,'(a36)')'Restricted Closed-Shell Wavefunction'
          write(naim,*)''
        end if

        do i=1,2*igr
          do j=1,2*igr
            c3(i,j)=0.0d0
          end do
        end do

!! transform sat matrix into the AO basis !!
        numorb=nalf
        if(inato.eq.1) numorb=ndim
        do i=1,numorb
          do j=1,igr
            xx=0.0d0
            do k=1,igr
              xx=xx+scr(k,i)*sat(k,j,jjat)
            end do
            c2(i,j)=xx
          end do
        end do
        do i=1,numorb
          do j=1,numorb
            xx=0.0d0
            do k=1,igr
              xx=xx+c2(i,k)*scr(k,j)
            end do
            c3(i,j)=xx
          end do
        end do
!! open-shell !!
        if(kop.eq.1.and.inato.eq.0)  then
          do i=1,nb
            do j=1,igr
              xx=0.0d0
              do k=1,igr
                xx=xx+cb(k,i)*sat(k,j,jjat)
              end do
              c2(i,j)=xx
            end do
          end do
          do i=1,nb
            do j=1,nb
              xx=0.0d0
              do k=1,igr
                xx=xx+c2(i,k)*cb(k,j)
              end do
              c3(i+nalf,j+nalf)=xx
            end do
          end do
        end if

        do i=1,igr
          do j=1,igr
            csave(i,j)=csave(i,j)+c3(i,j)
          end do
        end do

!! unrestricted single-determinant !!
        if(kop.eq.1.and.inato.eq.0) then
          if (imulli.ne.1) then
            do i=1,nalf+nb
              write(naim,*) (c3(i,j),j=1,i)
            end do
          else
            do i=1,nalf+nb
              write(naim,*) (c3(i,j),j=1,nalf+nb)
            end do
          end if
!! recalculate spin populations !!
          qalf=0.0d0
          do i=1,nalf
            qalf=qalf+c3(i,i)
          end do
          qbet=0.0d0
          do i=1,nb
            qbet=qbet+c3(i+nalf,i+nalf)
          end do
          write(naim,*) '  '
          write(naim,'(a41,e21.14)') 'ALPHA ELECTRONS (NA)',qalf
          write(naim,'(a41,e21.14)') 'BETA ELECTRONS (NB)',qbet
!! restricted single-determinant !!
        else if(inato.eq.0) then
          if (imulli.ne.1) then
            do i=1,nocc
              write(naim,*) (c3(i,j),j=1,i)
            end do
          else
            do i=1,nocc
              write(naim,*) (c3(i,j),j=1,nocc)
            end do
          end if
          write(naim,*) '  '
          write(naim,'(a41,e21.14)') 'ALPHA ELECTRONS (NA)',qat(jjat,1)/2.0d0
          write(naim,'(a41,e21.14)') 'BETA ELECTRONS (NB)',qat(jjat,1)/2.0d0
!! correlated wave function !!
        else if(inato.eq.1) then
          if (imulli.ne.1) then
            do i=1,numorb
              write(naim,*) (c3(i,j),j=1,i)
            end do
          else
            do i=1,numorb
              write(naim,*) (c3(i,j),j=1,numorb)
            end do
          end if
          write(naim,*) '  '
          write(naim,'(a41,e21.14)') 'ALPHA ELECTRONS (NA)',qat(jjat,1)/2.0d0
          write(naim,'(a41,e21.14)') 'BETA ELECTRONS (NB)',qat(jjat,1)/2.0d0
        end if


        write(naim,*)' '
        write(naim,*) 'NORMAL TERMINATION OF PROAIMV'

        close(naim)
      end do
      close(naim3)


      if(icuat.eq.nat) then
!! check for orthogonality of the MOs -- sum of all atomic overlap    !!
!! matrices (csave) should equal the identity.                         !!
        xmax=1.0d-2
        xmaxd=1.0d-3
        xmaxns=1.0d-4
        do i=1,igr
          do j=i,igr
            if(i.eq.j) then
              if(abs(csave(i,i))-1.0d0.gt.xmax) then
                xx=abs(csave(i,i))-1.0d0
                imaxd=i
                if(i.le.nalf) then
                  write(*,'(2x,a,f10.6,i3)') 'Dev. from normalization ',xx,imaxd
                end if
              end if
            else
              if(abs(csave(i,j)).gt.xmaxd) then
                xx=abs(csave(i,j))
                imax1=i
                imax2=j
                if(i.le.nalf.and.j.le.nalf) then
                  write(*,'(2x,a,f10.6,2i3)') 'Dev. from orth. in occ. set:',xx,
     1      imax1,imax2
                end if
              end if
              if(j.gt.i) then
                xxx=abs(csave(i,j))-abs(csave(j,i))
                if(xxx.gt.xmaxns) then
                  xx=xxx
                  imaxns1=i
                  imaxns2=j
                  write(*,'(2x,a,f10.6,2i3)') 'Dev. from hermiticity',xx,imaxns1,imaxns2
                end if
              end if
            end if
          end do
        end do
      end if
      deallocate(c3,scr,c2,csave)

      end

CCCCC
C FOR WRITTING CUBE FILES
CCCCC

      subroutine cubegen3(ifrag,icase)
      use ao_matrices
      use integration_grid
      use effao_mod, only: p0,p0net,p0gro,ip0 !! replaces common /effao/ -- see modules.f90 !!
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /iops/iopt(200)
      character*30 name
      common /filename/name

      character*30 name2
      character nameaim*60, charnu*2,atnu*2,charnu1*3
      character*2 mend(92)
      data mend/' H','He','Li','Be',' B',' C',' N',' O',
     $ ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl',
     $ 'Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe',
     $ 'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
     $ 'Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     $ 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I',
     $ 'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sn',
     $ 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
     $ 'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     $ 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     $ 'Th','Pa',' U'  /

      dimension pgrid(maxgrid,3),igrid(3),xgrid(3,3)
      dimension pop(maxat)          
      allocatable c0(:,:),xyz(:,:,:)

      ihirsh=iopt(6)
      imulli= Iopt(5) 
      ibcp=Iopt(14)
      iqtaim = Iopt(16)  
      inewbec = Iopt(31) 
      idofr=Iopt(40)
      jcubthr=iopt(41)
      kcubthr=iopt(42)

      allocate(c0(igr,igr))
      imaxo=ip0(ifrag)
      do i=1,igr
       do j=1,imaxo
        c0(i,j)=p0(i,j)
       end do
      end do

c setting actual effos to print, instead
      if(jcubthr.lt.0) then
       imaxeff=abs(jcubthr)
       imineff=abs(kcubthr)
      else
       imaxeff=0
       xmaxeff=float(jcubthr)*1.0d-3
       if (icase.eq.0) xmaxeff=2.0d0*xmaxeff
1      imaxeff= imaxeff+1
       if(p0net(imaxeff,ifrag).ge.xmaxeff) go to 1
       imineff=imaxo+1
       xmineff=float(kcubthr)*1.0d-3   
       if (icase.eq.0) xmineff=2.0d0*xmineff
2      imineff= imineff - 1
       if(p0net(imineff,ifrag).le.xmineff) go to 2
      end if
c

      if(imaxeff.gt.imineff) then
       write(*,*) 'No eff-AOs in the occupation range' 
       return
      end if
      write(*,22)'Generating cube files for eff-AOs',imaxeff,' to',imineff,' of atom/fragment ',ifrag
22    format (a33,i3,a3,i3,a18,i3)

        if (imulli.eq.1) then
         name2="mulliken"
        else if (imulli.gt.1) then
         name2="lowdin"
        else
        if (ihirsh.eq.0) then
         name2="becke"
         if(ibcp.eq.1) then 
            if(inewbec.eq.0)then
              name2="beckerho"
            else
             name2="tfvc"
            end if
         end if
         if(iqtaim.eq.1) name2="qtaim"
        else if(ihirsh.eq.1) then
         name2="hirsh"
        else if(ihirsh.eq.2) then
         name2="hirsh-it"
         do i=1,nat
          pop(i)=qat(i,1)
         end do
        end if
        end if

C grid points for cube in each dimension
      ngridp=40 

      do i=1,3
       igrid(i)=ngridp
c       if (i.eq.3) igrid(i)=100
c asuming rectangular grid...
       do j=1,3
       xgrid(i,j)=0.0d0                 
       end do
      end do

C Two options: For moleculs/fragments or centered on atoms      
C now active first option as it deals with fragments

      if(1.eq.1) then
      extra=3.2
      volume=1.0d0
c furthest x y z atomic posisiotns of the fragment
      do i=1,3
       xmax=-1.0d8
       xmin=1.0d8
       do jatom=1,nfrlist(ifrag)
        iatom=ifrlist(jatom,ifrag)
        if(coord(i,iatom).lt.xmin) xmin=coord(i,iatom)
        if(coord(i,iatom).gt.xmax) xmax=coord(i,iatom)
       end do
       xmin=xmin-extra
       xmax=xmax+extra
       dist=xmax-xmin
       volume=volume*dist
       xgrid(i,i)=dist/(igrid(i)-1.0d0)
       do j=1,igrid(i)
        pgrid(j,i)=xmin+(j-1)*xgrid(i,i)
       end do
      end do

c centering the grid on the atom
      else
      extra=3.0
      iatom=jjat  
      do i=1,3
       xmin=coord(i,iatom)
       xmax=coord(i,iatom)
       xmin=xmin-extra
       xmax=xmax+extra
       dist=xmax-xmin
       xgrid(i,i)=dist/(igrid(i)-1.0d0)
       do j=1,igrid(i)
        pgrid(j,i)=xmin+(j-1)*xgrid(i,i)
       end do
      end do
      end if

c Now the grid
      allocate ( xyz(igrid(1),igrid(2),igrid(3)))

      do ivec=imaxeff,imineff
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          xabs=pgrid(i,1)
          yabs=pgrid(j,2)
          zabs=pgrid(k,3)

          if(imulli.ne.0) then
           xyz(i,j,k)=orbxyz(c0,ivec,xabs,yabs,zabs)
          else
           ww=0.0d0
           do iatom=1,nfrlist(ifrag)
           jjat=ifrlist(iatom,ifrag)
           if(ihirsh.eq.0.and.iqtaim.eq.0) then
            ww=ww+wat(jjat,xabs,yabs,zabs)
           else if(ihirsh.eq.1) then
            ww=ww+wathirsh(jjat,xabs,yabs,zabs)
           else if(ihirsh.eq.2) then
            ww=ww+wathirsh2(jjat,xabs,yabs,zabs,pop)
           else if(iqtaim.eq.1) then 
c            call cubeqtaim(jjat,xabs,yabs,zabs,ww0)
            ww=ww+ww0
           end if
           end do
           xyz(i,j,k)=orbxyz(c0,ivec,xabs,yabs,zabs)*ww
          end if
         end do
        end do
       end do
c approximate normalization of orbital
       x0=0.0d0
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          x0=x0+xyz(i,j,k)*xyz(i,j,k)
         end do
        end do
       end do
       write(*,'(a25,f7.4)') 'Normalization from cube: ',x0*volume/(ngridp**3.0d0)

c      OUTPUT    

        if(idofr.eq.0) then
        read(mend(iznuc(ifrag)),'(A2)')charnu
         charnu=adjustl(charnu)
        else
         charnu="FR"
        end if   
        l0=len_trim(name)
        l1=len_trim(name2)
        l2=len_trim(charnu)
c assuming up to 99 atoms
        if(ifrag.lt.10) then
           write(atnu,'(i1)')ifrag
        else
           write(atnu,'(i2)')ifrag
        end if
        if(ivec.lt.10) then
         write(charnu1,'(i1)')ivec
        else if (ivec.lt.100) then
         write(charnu1,'(i2)')ivec
        else
         write(charnu1,'(i3)')ivec
        end if
         if(icase.ne.2) then
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)
         else                
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)//"beta"
         end if 
       j=len(nameaim)
       do i=1,j
        if(nameaim(i:i).eq.' ') then
         llen=i-1
         go to 10
        end if
       end do
  10   continue
       nameaim=nameaim(1:llen)//".cube"

       open(44,file=nameaim,status="unknown")
       rewind(44)
       write(44,*)'Cube generated with APOST-3D code '
       if(idofr.eq.0) then
       write(44,41) trim(name),name2,' EFFAO',ivec," for atom",
     + mend(iznuc(ifrag)),"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       else
       write(44,44) trim(name),name2,' EFFAO',ivec," for frag",
     + ifrag,"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       end if
       write(44,42) nat,(pgrid(1,j),j=1,3)
       do i=1,3
        write(44,42) igrid(i),(xgrid(i,j),j=1,3)
       end do
       do i=1,nat
        write(44,43) iznuc(i),zn(i),(coord(j,i),j=1,3)
       end do
       ione=1
       write(44,'(2i5)')ione,ione
       do i=1,igrid(1)
       do j=1,igrid(2)
        write(44,40)(xyz(i,j,k),k=1,igrid(3))
       end do
       end do
       close(44)
      end do
41    format(a8,x,a8,a7,i3,a9,a2,a11,f8.4,a9,f8.4)
44    format(a8,x,a8,a7,i3,a9,i2,a11,f8.4,a9,f8.4)
42    format(i5,3f12.6)
43    format(i5,4f12.6)
40    format(6e13.5)
     
      deallocate(c0,xyz)
      end


CCCCC
C FOR WRITTING FCHK FILES
CCCCC
       subroutine rmat(iunit,key,ival,jval,ndim,rmatrix)
       implicit double precision (a-h,o-z)
       include 'parameter.h'
       dimension rmatrix(ndim,ndim) 
       character(len=*) ::  key
       integer ival,jval
       character(len=43) ::  title
       title=adjustl(key)
       write(iunit,'(A43,A6,I12)') title,"R   N=",ival*jval
       write(iunit,'(5ES16.8)') ((rmatrix(i,j),i=1,jval),j=1,ival)
       end

       subroutine rarr(iunit,key,ival,ndim,rarray)
       implicit double precision (a-h,o-z)
       include 'parameter.h'
       character(len=*) ::  key
       integer ival
       dimension rarray(ndim)
       character(len=43) ::  title
       title=adjustl(key)
       write(iunit,'(A43,A6,I12)') title,"R   N=",ival
       write(iunit,'(5ES16.8)') (rarray(i),i=1,ival)
       end

       subroutine ival(iunit,key,ivalue)
       implicit double precision (a-h,o-z)
       character(len=*) ::  key
       integer ivalue
       character(len=43) ::  title
       title=adjustl(key)
       write(iunit,'(A43,A,I17)') title,"I",ivalue
       end

        SUBROUTINE PRINTMAT(N_orbital,S)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        include 'parameter.h'
        DIMENSION S(N_orbital,N_orbital)

        nblock=N_orbital/5
        if(nblock*5.ne.N_orbital) nblock=nblock+1
        do k=1,nblock
          ii=min0(N_orbital,5*k)
          write(*,'(4x,10(7X,i6,a1))') (j," ",j=5*(k-1)+1,ii)
          do i=5*(k-1)+1,N_orbital
!            ii=5*k
!            if(i.lt.5*k) ii=i
            ii=min0(i,5*k)
            write(*,'(i7,5(x,d13.6))') i,(S(i,j),j=5*(k-1)+1,ii)
          end do
        end do
        END SUBROUTINE

      subroutine cubegen3_mhg(ifrag,icase)
      use ao_matrices
      use integration_grid
      use effao_mod, only: p0,p0net,p0gro,ip0 !! replaces common /effao/ -- see modules.f90 !!
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /iops/iopt(200)
      character*30 name
      common /filename/name

      character*30 name2
      character nameaim*60, charnu*2,atnu*2,charnu1*3
      character*2 mend(92)
      data mend/' H','He','Li','Be',' B',' C',' N',' O',
     $ ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl',
     $ 'Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe',
     $ 'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
     $ 'Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     $ 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I',
     $ 'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sn',
     $ 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
     $ 'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     $ 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     $ 'Th','Pa',' U'  /

      dimension pgrid(maxgrid,3),igrid(3),xgrid(3,3)
      dimension pop(maxat)          
      allocatable c0(:,:),xyz(:,:,:)

      ihirsh=iopt(6)
      imulli= Iopt(5) 
      ibcp=Iopt(14)
      iqtaim = Iopt(16)  
      inewbec = Iopt(31) 
      idofr=Iopt(40)
      jcubthr=iopt(41)
      kcubthr=iopt(42)

      allocate(c0(igr,igr))
      imaxo=ip0(ifrag)
      do i=1,igr
       do j=1,imaxo
        c0(i,j)=p0(i,j)
       end do
      end do

c setting actual effos to print, instead
      if(jcubthr.lt.0) then
       imaxeff=abs(jcubthr)
       imineff=abs(kcubthr)
      else
       imaxeff=0
       xmaxeff=float(jcubthr)*1.0d-3
       if (icase.eq.0) xmaxeff=2.0d0*xmaxeff
1      imaxeff= imaxeff+1
       if(p0net(imaxeff,ifrag).ge.xmaxeff) go to 1
       imineff=imaxo+1
       xmineff=float(kcubthr)*1.0d-3   
       if (icase.eq.0) xmineff=2.0d0*xmineff
2      imineff= imineff - 1
       if(p0net(imineff,ifrag).le.xmineff) go to 2
      end if
c

      if(imaxeff.gt.imineff) then
       write(*,*) 'No eff-AOs in the occupation range' 
       return
      end if
      write(*,22)'Generating cube files for eff-AOs',imaxeff,' to',imineff,' of atom/fragment ',ifrag
22    format (a33,i3,a3,i3,a18,i3)

      name2="mhg"

C grid points for cube in each dimension
      ngridp=40 

      do i=1,3
       igrid(i)=ngridp
c       if (i.eq.3) igrid(i)=100
c asuming rectangular grid...
       do j=1,3
       xgrid(i,j)=0.0d0                 
       end do
      end do

C Two options: For moleculs/fragments or centered on atoms      
C now active first option as it deals with fragments

      if(1.eq.1) then
      extra=5.2
      volume=1.0d0
c furthest x y z atomic posisiotns of the fragment
      do i=1,3
       xmax=-1.0d8
       xmin=1.0d8
       do jatom=1,nfrlist(ifrag)
        iatom=ifrlist(jatom,ifrag)
        if(coord(i,iatom).lt.xmin) xmin=coord(i,iatom)
        if(coord(i,iatom).gt.xmax) xmax=coord(i,iatom)
       end do
       xmin=xmin-extra
       xmax=xmax+extra
       dist=xmax-xmin
       volume=volume*dist
       xgrid(i,i)=dist/(igrid(i)-1.0d0)
       do j=1,igrid(i)
        pgrid(j,i)=xmin+(j-1)*xgrid(i,i)
       end do
      end do

c centering the grid on the atom
      else
      extra=3.0
      iatom=jjat  
      do i=1,3
       xmin=coord(i,iatom)
       xmax=coord(i,iatom)
       xmin=xmin-extra
       xmax=xmax+extra
       dist=xmax-xmin
       xgrid(i,i)=dist/(igrid(i)-1.0d0)
       do j=1,igrid(i)
        pgrid(j,i)=xmin+(j-1)*xgrid(i,i)
       end do
      end do
      end if

c Now the grid
      allocate ( xyz(igrid(1),igrid(2),igrid(3)))

      do ivec=imaxeff,imineff
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          xabs=pgrid(i,1)
          yabs=pgrid(j,2)
          zabs=pgrid(k,3)
           xyz(i,j,k)=orbxyz(c0,ivec,xabs,yabs,zabs)
         end do
        end do
       end do
c approximate normalization of orbital
       x0=0.0d0
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          x0=x0+xyz(i,j,k)*xyz(i,j,k)
         end do
        end do
       end do
       write(*,'(a25,f7.4)') 'Normalization from cube: ',x0*volume/(ngridp**3.0d0)

c      OUTPUT    

        if(idofr.eq.0) then
        read(mend(iznuc(ifrag)),'(A2)')charnu
         charnu=adjustl(charnu)
        else
         charnu="FR"
        end if   
        l0=len_trim(name)
        l1=len_trim(name2)
        l2=len_trim(charnu)
c assuming up to 99 atoms
        if(ifrag.lt.10) then
           write(atnu,'(i1)')ifrag
        else
           write(atnu,'(i2)')ifrag
        end if
        if(ivec.lt.10) then
         write(charnu1,'(i1)')ivec
        else if (ivec.lt.100) then
         write(charnu1,'(i2)')ivec
        else
         write(charnu1,'(i3)')ivec
        end if
         if(icase.ne.2) then
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)
         else                
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)//"beta"
         end if 
       j=len(nameaim)
       do i=1,j
        if(nameaim(i:i).eq.' ') then
         llen=i-1
         go to 10
        end if
       end do
  10   continue
       nameaim=nameaim(1:llen)//".cube"

       open(44,file=nameaim,status="unknown")
       rewind(44)
       write(44,*)'Cube generated with APOST-3D code '
       if(idofr.eq.0) then
       write(44,41) trim(name),name2,' EFFAO',ivec," for atom",
     + mend(iznuc(ifrag)),"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       else
       write(44,44) trim(name),name2,' EFFAO',ivec," for frag",
     + ifrag,"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       end if
       write(44,42) nat,(pgrid(1,j),j=1,3)
       do i=1,3
        write(44,42) igrid(i),(xgrid(i,j),j=1,3)
       end do
       do i=1,nat
        write(44,43) iznuc(i),zn(i),(coord(j,i),j=1,3)
       end do
       ione=1
       write(44,'(2i5)')ione,ione
       do i=1,igrid(1)
       do j=1,igrid(2)
        write(44,40)(xyz(i,j,k),k=1,igrid(3))
       end do
       end do
       close(44)
      end do
41    format(a8,x,a8,a7,i3,a9,a2,a11,f8.4,a9,f8.4)
44    format(a8,x,a8,a7,i3,a9,i2,a11,f8.4,a9,f8.4)
42    format(i5,3f12.6)
43    format(i5,4f12.6)
40    format(6e13.5)
     
      deallocate(c0,xyz)
      end

!! ***** !!

!! ********************************************************************** !!
!! subroutine: cubegen4                                                   !!
!! purpose: writes one Gaussian-style .cube file per requested EFO        !!
!! (imaxeff..imineff, thresholded by MAX_OCC/MIN_OCC) of fragment/atom    !!
!! ifrag -- adaptive grid, fixed point spacing (# CUBE SPACING, bohr)     !!
!! with padding scaled by RADIUS_SCALE times the extremal atom's          !!
!! covalent radius in each direction, so point count (not spacing)        !!
!! grows with fragment size. icase selects RHF/UHF-alpha/UHF-beta/UEOS    !!
!! paired/unpaired naming; imulli selects orbital-value output (Mulliken/ !!
!! Lowdin) vs AIM-weighted density (Becke/TFVC/Hirshfeld/QTAIM).          !!
!! arguments: ifrag (in) -- fragment/atom index, icase (in) -- 0-4,       !!
!! see above                                                              !!
!! author: PSalse, MGimf                                                  !!
!! ********************************************************************** !!
      subroutine cubegen4(ifrag,icase)
      use ao_matrices
      use integration_grid
      use effao_mod, only: p0,p0net,p0gro,ip0 !! replaces common /effao/ -- see modules.f90 !!
      use input_options_mod, only: cubespacing,cuberadscale
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /atomrad/atr(maxat),dist(maxat,maxat)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /iops/iopt(200)
      character*30 name
      common /filename/name

      character*30 name2
      character nameaim*60, charnu*2,atnu*2,charnu1*3
      character*2 mend(92)
      data mend/' H','He','Li','Be',' B',' C',' N',' O',
     $ ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl',
     $ 'Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe',
     $ 'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
     $ 'Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     $ 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I',
     $ 'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sn',
     $ 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
     $ 'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     $ 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     $ 'Th','Pa',' U'  /

      dimension pgrid(maxgrid,3),igrid(3),xgrid(3,3)
      dimension pop(maxat)          
      allocatable c0(:,:),xyz(:,:,:)

      ihirsh=iopt(6)
      imulli= Iopt(5) 
      ibcp=Iopt(14)
      iqtaim = Iopt(16)  
      inewbec = Iopt(31) 
      idofr=Iopt(40)
      jcubthr=iopt(41)
      kcubthr=iopt(42)


      imaxo=ip0(ifrag)
      ALLOCATE(c0(igr,igr))
      do i=1,igr
        do j=1,imaxo
          c0(i,j)=p0(i,j)
        end do
      end do

!! setting actual effos to print, instead !!
      if(jcubthr.lt.0) then
       imaxeff=abs(jcubthr)
       imineff=abs(kcubthr)
      else
       imaxeff=0
       xmaxeff=float(jcubthr)*1.0d-3
       if (icase.eq.0.or.icase.eq.3) xmaxeff=2.0d0*xmaxeff

1      imaxeff= imaxeff+1
       if(p0net(imaxeff,ifrag).ge.xmaxeff) go to 1
       imineff=imaxo+1
       xmineff=float(kcubthr)*1.0d-3   
       if (icase.eq.0.or.icase.eq.3) xmineff=2.0d0*xmineff

2      imineff= imineff - 1
       if(p0net(imineff,ifrag).le.xmineff) go to 2
      end if

      if(imaxeff.gt.imineff) then
       write(*,*) ' No eff-AOs in the occupation range'
       write(*,*) " "
       deallocate(c0)
       return
      end if

      write(*,22)'Generating cube files for eff-AOs',imaxeff,' to',imineff,' of atom/fragment ',ifrag
22    format (a33,i3,a3,i3,a18,i3)

        if (imulli.eq.1) then
         name2="mulliken"
        else if (imulli.gt.1) then
         name2="lowdin"
        else
        if (ihirsh.eq.0) then
         name2="becke"
         if(ibcp.eq.1) then 
            if(inewbec.eq.0)then
              name2="beckerho"
            else
             name2="tfvc"
            end if
         end if
         if(iqtaim.eq.1) name2="qtaim"
        else if(ihirsh.eq.1) then
         name2="hirsh"
        else if(ihirsh.eq.2) then
         name2="hirsh-it"
         do i=1,nat
          pop(i)=qat(i,1)
         end do
        end if
        end if
        if(icase.eq.3) name2=trim(name2)//"_paired"
        if(icase.eq.4) name2=trim(name2)//"_unpaired"

!! assuming rectangular grid !!
       xgrid=0.0d0

!! adaptive-size cube: fixed point spacing (# CUBE SPACING), padding      !!
!! scaled by RADIUS_SCALE times the extremal atom's covalent radius       !!
      xmesh=cubespacing
      rrmax=cuberadscale
      volume=1.0d0
!! furthest x y z atomic positions of the fragment !!
      do i=1,3
       xmax=-1.0d8
       xmin=1.0d8
       do jatom=1,nfrlist(ifrag)
        iatom=ifrlist(jatom,ifrag)
        if(coord(i,iatom).lt.xmin) then
         xmin=coord(i,iatom)
         iiatom=iatom
        end if
        if(coord(i,iatom).gt.xmax) then
         xmax=coord(i,iatom)
         iiiatom=iatom
        end if
       end do
       xmin=xmin-rrmax*atr(iiatom)
       xmax=xmax+rrmax*atr(iiiatom)
       dist0=xmax-xmin
       xgrid(i,i)=xmesh
       igrid(i)=int(dist0/xmesh)+1
       do j=1,igrid(i)
        pgrid(j,i)=xmin+(j-1)*xgrid(i,i)
       end do
       volume=volume*dist0
      end do

!! now the grid !!
      write(*,'(a21,3i4)')'Size of cube (x,y,z):',(igrid(i),i=1,3)
      allocate ( xyz(igrid(1),igrid(2),igrid(3)))

      do ivec=imaxeff,imineff

!! each (i,j,k) grid point is independent, writing only its own xyz        !!
!! slot; orbxyz/wat/wathirsh(2) are pure functions of their arguments      !!
!! plus read-only shared state (c0/coord/COMMON), safe to call in          !!
!! parallel (see wat.f)                                                    !!
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k,xabs,yabs,zabs,ww,iatom,jjat)
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          xabs=pgrid(i,1)
          yabs=pgrid(j,2)
          zabs=pgrid(k,3)

          if(imulli.ne.0) then
           xyz(i,j,k)=orbxyz(c0,ivec,xabs,yabs,zabs)
          else
           ww=0.0d0
           do iatom=1,nfrlist(ifrag)
            jjat=ifrlist(iatom,ifrag)
            if(ihirsh.eq.0.and.iqtaim.eq.0) then
              ww=ww+wat(jjat,xabs,yabs,zabs)
            else if(ihirsh.eq.1) then
              ww=ww+wathirsh(jjat,xabs,yabs,zabs)
            else if(ihirsh.eq.2) then
              ww=ww+wathirsh2(jjat,xabs,yabs,zabs,pop)
            else if(iqtaim.eq.1) then
!! QTAIM cube weighting was never wired up (cubeqtaim doesn't exist       !!
!! codebase-wide) -- ww stays 0 here, so a QTAIM-weighted cube would be   !!
!! all zeros rather than erroring. Unreachable in practice regardless:    !!
!! main.f:218 hard-stops the whole run at startup when iqtaim=1.          !!
            end if
           end do
           xyz(i,j,k)=orbxyz(c0,ivec,xabs,yabs,zabs)*ww
          end if
         end do
        end do
       end do
!$OMP END PARALLEL DO

!! approximate normalization of orbital -- same independence argument,   !!
!! reduction on x0                                                        !!
       x0=ZERO
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k) REDUCTION(+:x0)
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          x0=x0+xyz(i,j,k)*xyz(i,j,k)
         end do
        end do
       end do
!$OMP END PARALLEL DO
       write(*,'(a25,f7.4)') 'Normalization from cube: ',x0*volume/(igrid(1)*igrid(2)*igrid(3))
       write(*,*) " "

!! output !!

        if(idofr.eq.0) then
        read(mend(iznuc(ifrag)),'(A2)')charnu
         charnu=adjustl(charnu)
        else
         charnu="FR"
        end if   
!! assuming up to 99 atoms !!
        if(ifrag.lt.10) then
           write(atnu,'(i1)')ifrag
        else
           write(atnu,'(i2)')ifrag
        end if
        if(ivec.lt.10) then
         write(charnu1,'(i1)')ivec
        else if (ivec.lt.100) then
         write(charnu1,'(i2)')ivec
        else
         write(charnu1,'(i3)')ivec
        end if
         if(icase.ne.2) then
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)
         else                
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)//"beta"
         end if 
       j=len(nameaim)
       do i=1,j
        if(nameaim(i:i).eq.' ') then
         llen=i-1
         go to 10
        end if
       end do
  10   continue
       nameaim=nameaim(1:llen)//".cube"

       open(44,file=nameaim,status="unknown")
       rewind(44)
       write(44,*)'Cube generated with APOST-3D code '
       if(idofr.eq.0) then
       write(44,41) trim(name),name2,' EFFAO',ivec," for atom",
     + mend(iznuc(ifrag)),"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       else
       write(44,44) trim(name),name2,' EFFAO',ivec," for frag",
     + ifrag,"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       end if
       write(44,42) nat,(pgrid(1,j),j=1,3)
       do i=1,3
        write(44,42) igrid(i),(xgrid(i,j),j=1,3)
       end do
       do i=1,nat
        write(44,43) iznuc(i),zn(i),(coord(j,i),j=1,3)
       end do
       ione=1
       write(44,'(2i5)')ione,ione
       do i=1,igrid(1)
       do j=1,igrid(2)
        write(44,40)(xyz(i,j,k),k=1,igrid(3))
       end do
       end do
       close(44)
      end do
41    format(a8,x,a8,a7,i3,a9,a2,a11,f8.4,a9,f8.4)
44    format(a8,x,a8,a7,i3,a9,i2,a11,f8.4,a9,f8.4)
42    format(i5,3f12.6)
43    format(i5,4f12.6)
40    format(6e13.5)
     
      deallocate(c0,xyz)
      end

!! ***** !!