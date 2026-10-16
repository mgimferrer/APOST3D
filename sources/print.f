!! *********************************************************************** !!
!! ALL FORMATTED OUTPUT -- startup/input banners, shared matrix printers,  !!
!! fragment-analysis printers, the AIMPAC/PROAIMV .int-file writer, the    !!
!! cube-file generator                                                     !!
!! Startup/input banners:                                                  !!
!!   kiir                -- startup banner (version, citations, features)  !!
!!   print_input_summary -- "digested" echo of the active .inp keywords    !!
!! Shared building block, used by every banner-style print_box call below: !!
!!   print_box           -- auto-sized rule around a title, no manual      !!
!!                           blank-line management needed by callers       !!
!!   print_subbox        -- same, no leading blank -- for a box that       !!
!!                           immediately follows another print_box/        !!
!!                           print_subbox with nothing printed between     !!
!! Ending a run on an error:                                               !!
!!   apost_stop          -- "STOP <msg>" as the last output line, exit 1   !!
!! Per-atom vector/matrix printers (no border):                            !!
!!   VPRINT               -- per-atom vector or 2-column matrix            !!
!!   MPRINT_NLOP           -- per-atom 3-column (X/Y/Z) matrix             !!
!!   MPRINT                -- N x N matrix, chunked 6 columns at a time    !!
!! Bordered N x N atom-pair matrix printer, the style used throughout      !!
!! ENPART:                                                                 !!
!!   MPRINT2               -- same chunking as MPRINT, plus an 80-column   !!
!!                             rule around the header and each row block   !!
!! Fragment values (a per-atom quantity summed over the atoms of each     !!
!! fragment, /frlist/ map) and their printers:                             !!
!!   frag_sum_mat / frag_sum_vec -- fragment matrix / table, no printing   !!
!!   mprintnoat            -- fragment x fragment matrix, as MPRINT2       !!
!!   group_by_frag_mat     -- section with the fragment matrix and total   !!
!!   group_by_frag_vec     -- section with the fragment table and total    !!
!!   print_frag_table      -- section with a given fragment table          !!
!! AIMPAC/PROAIMV interface (opt-in, # METHOD / DOINT):                    !!
!!   print_int              -- writes one .int file per atom (plus a        !!
!!                             .files index) holding its MO/NO-basis        !!
!!                             overlap matrix, for external programs       !!
!!                             (ESI-3D, FCALC)                             !!
!! Cube-file generation (opt-in, # CUBE):                                  !!
!!   cubegen_new            -- writes one Gaussian-style .cube file per      !!
!!                             requested EFO of a fragment/atom, adaptive   !!
!!                             grid sized from SPACING/RADIUS_SCALE          !!
!! *********************************************************************** !!

!! ***** !!

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
!! subroutine: print_subbox                                              !!
!! purpose: same as print_box, but without the leading blank line -- for !!
!!   a box that immediately follows another print_box call with nothing  !!
!!   printed in between, so the two don't stack into a double gap.       !!
!!   (a separate subroutine rather than an optional argument on          !!
!!   print_box: this codebase has no explicit interfaces, and gfortran   !!
!!   requires one for optional arguments at every one of print_box's     !!
!!   ~30 call sites.)                                                    !!
!! arguments:                                                            !!
!!   text (in) -- title to print, no leading/trailing padding needed     !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      SUBROUTINE print_subbox(text)
      IMPLICIT NONE
      character(len=*), intent(in) :: text
      integer :: n

      n=len_trim(text)+4
      write(*,'(2x,a)')    repeat('-',n)
      write(*,'(2x,2x,a)') trim(text)
      write(*,'(2x,a,/)')  repeat('-',n)

      END SUBROUTINE print_subbox

!! ***** !!

!! ********************************************************************* !!
!! subroutine: apost_stop                                                !!
!! purpose: ends the run on an error, in place of a stop statement:      !!
!!   writes "STOP <msg>" to stdout and flushes it, so the message is the !!
!!   last line of the output even when stdout is redirected (a stop      !!
!!   statement writes to unbuffered stderr, above the buffered output),  !!
!!   then exits with code 1 so job scripts can tell an error from a      !!
!!   normal end. An empty msg writes no STOP line (the caller has        !!
!!   printed the reason).                                                !!
!! arguments:                                                            !!
!!   msg (in) -- the message, as the stop statement had it               !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine apost_stop(msg)
      use json_mod, only: json_close
      character(len=*) msg

      call json_close('stopped',msg)
      if(len_trim(msg).gt.0) write(*,'(a,a)') 'STOP ',msg
      flush(6)
      call exit(1)

      end subroutine apost_stop

!! ***** !!

!! ********************************************************************* !!
!! subroutine: kiir                                                      !!
!! purpose: prints the startup banner -- version/date identity, feature  !!
!! availability and per-method citation directory, acknowledgments.      !!
!! arguments: none                                                       !!
!! author:                                                               !!
!! ********************************************************************* !!
      SUBROUTINE kiir
      use version_mod, only: version_number, version_date
      IMPLICIT NONE
      character(len=80) :: vline

      write(*,'(a)') '  ------------------------------------------------------------------------------'
      write(*,*)
      write(*,'(a)') '                     ___    ____  ____  ___________   _____ ____'
      write(*,'(a)') '                    /   |  / __ \/ __ \/ ___/_  __/  |__  // __ \'
      write(*,'(a)') '                   / /| | / /_/ / / / /\__ \ / /_____ /_ </ / / /'
      write(*,'(a)') '                  / ___ |/ ____/ /_/ /___/ // /_____/__/ / /_/ /'
      write(*,'(a)') '                 /_/  |_/_/    \____//____//_/     /____/_____/'
      write(*,*)
      vline=repeat(' ',24)//'Version '//version_number//'   --   '//
     +  version_date
      write(*,'(a)') trim(vline)
      write(*,'(a)') '                         --------------------------------'
      write(*,*)
      write(*,'(a)') '      Real-space and Hilbert-space tools for wave function analysis'
      write(*,*)
      write(*,'(a)') '      Documentation: https://apost3d.readthedocs.io'
      write(*,'(a)') '      Source code & issue tracker: https://github.com/mgimferrer/APOST3D'
      write(*,*)
      write(*,'(a)') '      Cite this program as:'
      write(*,'(a)') '      ---------------------'
      write(*,'(a)') '        P. Salvador, E. Ramos-Cordoba, M. Montilla, L. Pujal and M. Gimferrer'
      write(*,'(a)') '        J. Chem. Phys., 2024, 160, 172502 DOI: 10.1063/5.0206187'
      write(*,*)
      write(*,'(a)') '      e-mail: pedro.salvador@udg.edu, mgimferrer18@gmail.com'
      write(*,*)
      write(*,'(a)') '      Available atomic definitions:'
      write(*,'(a)') '      ----------------------------'
      write(*,*)
      write(*,'(a)') '      Real space:'
      write(*,'(a)') '        Becke, J. Chem. Phys. 88 2547 1988'
      write(*,'(a)') '        Hirshfeld, Theor. Chim. Acta 44 129 1977'
      write(*,'(a)') '        Hirshfeld-Iterative, J. Chem. Phys. 126 144111 2007'
      write(*,'(a)') '        Topological fuzzy Voronoi cells (TFVC), J. Chem. Phys. 139 071103 2013'
      write(*,*)
      write(*,'(a)') '      Hilbert space: Mulliken, Lowdin, Lowdin-Davidson, natural atomic'
      write(*,'(a)') '        orbitals (NAO)'
      write(*,*)
      write(*,'(a)') '      Analyses and their references:'
      write(*,'(a)') '      ------------------------------'
      write(*,*)
      write(*,'(a)') '        A) Atomic and overlap populations, bond orders and valences'
      write(*,'(a)') '           I. Mayer and P. Salvador, Chem. Phys. Lett. 383 368-375 2004'
      write(*,*)
      write(*,'(a)') '        B) Hartree-Fock molecular energy decomposition'
      write(*,'(a)') '           P. Salvador, M. Duran and I. Mayer, J. Chem. Phys. 115 1153-1157 2001'
      write(*,'(a)') '           P. Salvador and I. Mayer, J. Chem. Phys. 120 5046-5052 2004'
      write(*,*)
      write(*,'(a)') '        C) KS-DFT molecular energy decomposition'
      write(*,'(a)') '           P. Salvador and I. Mayer, J. Chem. Phys. 126 234113 2007'
      write(*,'(a)') '           M. Gimferrer and P. Salvador, J. Chem. Phys. 158 234105 2023'
      write(*,*)
      write(*,'(a)') '        D) Molecular energy decomposition for CASSCF/CISD wavefunctions'
      write(*,*)
      write(*,'(a)') '        E) Effective atomic orbitals'
      write(*,'(a)') '           I. Mayer, J. Phys. Chem. 100 6249 1996'
      write(*,'(a)') '           I. Mayer and P. Salvador, J. Chem. Phys. 130 234106 2009'
      write(*,'(a)') '           E. Ramos-Cordoba et al., J. Chem. Phys. 138 214107 2013'
      write(*,'(a)') '           Review: M. Gimferrer, Theor. Chem. Acc. 145 74 2026'
      write(*,*)
      write(*,'(a)') '        F) Local spin analysis'
      write(*,'(a)') '           E. Ramos-Cordoba et al., J. Chem. Theory Comput. 8 1270-1279 2012'
      write(*,'(a)') '           E. Ramos-Cordoba et al., Phys. Chem. Chem. Phys. 14 15291-15298 2012'
      write(*,*)
      write(*,'(a)') '        G) Effective oxidation states analysis'
      write(*,'(a)') '           E. Ramos-Cordoba et al., J. Chem. Theory Comput. 11 1501-1508 2015'
      write(*,'(a)') '           M. Gimferrer and P. Salvador, manuscript in preparation, 2026'
      write(*,*)
      write(*,'(a)') '        H) Oxidation states from localized orbitals'
      write(*,'(a)') '           M. Gimferrer, G. Comas-Vila and P. Salvador, Molecules 25 234 2020'
      write(*,'(a)') '           M. Gimferrer et al., Inorg. Chem. 59 15410-15420 2020'
      write(*,'(a)') '           M. Gimferrer et al., J. Chem. Theory Comput. 18 309-322 2022'
      write(*,*)
      write(*,'(a)') '        I) Decomposition of EDA quantities into one- and two-center IQA terms'
      write(*,'(a)') '           M. Gimferrer et al., J. Chem. Theory Comput. 19 3469-3485 2023'
      write(*,'(a)') '           (under reimplementation in version 5)'
      write(*,*)
      write(*,'(a)') '        J) Origin-independent decomposition of static polarizabilities'
      write(*,'(a)') '           M. Montilla et al., J. Chem. Theory Comput. 17 1098-1105 2021'
      write(*,*)
      write(*,'(a)') '  ------------------------------------------------------------------------------'
      write(*,*)
      write(*,'(a)') '      The program has been written by using parts of the program APOST by'
      write(*,'(a)') '      I. Mayer and A. Hamza, Budapest, 2000-2003.'
      write(*,*)
      write(*,'(a)') '      The numerical integration utilizes the subroutines for Lebedev'
      write(*,'(a)') '      quadrature downloaded from CCL. The appropriate reference is:'
      write(*,'(a)') '      V.I. Lebedev, and D.N. Laikov "A quadrature formula for the sphere of'
      write(*,'(a)') '      the 131st algebraic order of accuracy" Doklady Mathematics, 59 477-481'
      write(*,'(a)') '      1999.'
      write(*,*)
      write(*,'(a)') '      The program uses the libxc library when necessary, through the'
      write(*,'(a)') '      Fortran interfaces provided by its authors (https://libxc.gitlab.io/).'
      write(*,*)
      write(*,'(a)') '      We are extremely grateful for the possibility of using these routines!'
      write(*,*)
      write(*,'(a)') '  ------------------------------------------------------------------------------'
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
      use output_mod, only: print_value
      use integration_grid, only: nrad_in,nang_in
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
      character*60 cval
      character*20 cnum,cfmt

      call print_box('INPUT SUMMARY')

!! ----------------------------------------------------------------- !!
!! # METHOD -- one line per active category                          !!
!! ----------------------------------------------------------------- !!
      write(*,'(2x,a)') '# METHOD'

!! real-space atomic definitions: TFVC also sets the BECKE-RHO/NEWBEC !!
!! flags (read_input.f), so those two are named only without TFVC     !!
      cbuf=' '
      if(itfvc.eq.1.and.itfvcdef.eq.0) cbuf=trim(cbuf)//
     +  '; TFVC (topological fuzzy Voronoi cells)'
      if(itfvcdef.eq.1) cbuf=trim(cbuf)//
     +  '; TFVC (topological fuzzy Voronoi cells, none given: default)'
      if(ibecke.eq.1) cbuf=trim(cbuf)//'; BECKE (fixed empirical radii)'
      if(ibcp.eq.1.and.itfvc.ne.1) cbuf=trim(cbuf)//
     +  '; BECKE-RHO (radii from the density)'
      if(inewbec.eq.1.and.itfvc.ne.1) cbuf=trim(cbuf)//'; NEWBEC'
      if(ihirsh.eq.1) cbuf=trim(cbuf)//'; HIRSH (Hirshfeld)'
      if(ihirsh.eq.2) cbuf=trim(cbuf)//'; HIRSH-IT (iterative Hirshfeld)'
      if(iqtaim.eq.1) cbuf=trim(cbuf)//'; QTAIM'
      if(len_trim(cbuf).gt.0) call print_value('Atomic definition',
     +  trim(adjustl(cbuf(3:))))
      if(imulli.ge.1) then
        if(imulli.eq.1) cval='MULLIKEN (Hilbert space)'
        if(imulli.eq.2) cval='LOWDIN (Hilbert space)'
        if(imulli.eq.3) cval='LOWDIN-DAVIDSON (Hilbert space)'
        if(imulli.eq.4) cval='NAO-BASIS (Hilbert space)'
        if(imulli.eq.5) cval='LOWDIN-W (Hilbert space)'
        call print_value('Atomic definition',trim(cval))
      end if
      if(ndens0.gt.1) call print_value('Density (DENS)',ndens0)
      if(iopop.eq.1) call print_value('Overlap populations','OPOP')

      if(ieffao.ge.1.and.ieos.ne.1.and.iueos.ne.1) then
        if(ieffao.eq.1) cval='EFFAO'
        if(ieffao.eq.2) cval='UEFFAO'
        if(ieffao.eq.3) cval='EFFAO-U (paired/unpaired densities)'
        call print_value('Effective atomic orbitals',trim(cval))
      end if
      if(ieos.eq.1) call print_value('Oxidation states',
     +  'EOS (effective oxidation states)')
      if(iueos.eq.1) call print_value('Oxidation states',
     +  'GEOS (generalized effective oxidation states)')
      if(ieoscent.eq.1) call print_value('Oxidation states',
     +  'OS-CENTROID')
      if(ioslo.eq.1) call print_value('Oxidation states',
     +  'OSLO (oxidation states from localized orbitals)')
      if(iloba.eq.1) call print_value('Oxidation states',
     +  'LOBA (localized orbital bonding analysis)')
      if(ienpart.eq.1) call print_value('Energy decomposition',
     +  'ENPART (IQA terms, see # ENPART)')
      if(ispin.eq.1) call print_value('Local spin','SPIN')
      if(ipolar.eq.1) call print_value('Static polarizability','POLAR')
      if(iscattfact.eq.1) call print_value('X-ray scattering factors',
     +  'SCATT-FACT')
      if(ilaplacian.eq.1) call print_value('Density Laplacian',
     +  'LAPLACIAN')
      if(ipca.eq.1) call print_value('Principal component analysis',
     +  'PCA')
      if(idoint.eq.1) call print_value('Atomic overlaps (ESI-3D, FCALC)',
     +  'DOINT (one .int file per atom)')
      if(ielcount.eq.1) call print_value('Electron counting (NCTAIM)',
     +  'ELCOUNT')
      if(idofr.eq.1) then
        write(cval,'(i0,a)') icufr,' (see # FRAGMENTS)'
        call print_value('Fragments',trim(cval))
      end if
      if(idoat.eq.1) then
        write(cval,'(i0,a)') icuat,' (see # ATOMS)'
        call print_value('Selected atoms',trim(cval))
      end if
      if(iqchem.eq.1) call print_value('Wavefunction source','QCHEM')
      if(imokit.eq.1) call print_value('Wavefunction source','MOKIT')
      if(icube.eq.1) then
        call print_value('Cube files','CUBE (see # CUBE)')
        if(ieffao.eq.0) then
          write(*,'(2x,a)') 'WARNING: no cube file is written: no'//
     +      ' analysis in # METHOD makes EFOs (EFFAO, UEFFAO,'
          write(*,'(2x,a)') 'EFFAO-U, EOS, GEOS)'
        end if
      end if
      if(iaccur.eq.1) call print_value('Output precision',
     +  'FULLPRECISION (13 decimals)')
      if(iverbose.eq.1) call print_value('Output detail',
     +  'full (VERBOSE 1)')
      if(inopop.eq.1) call print_value('Population output',
     +  'NOPOPU (suppressed)')
      if(nrad_in.gt.0.or.nang_in.gt.0) then
        cval='default'
        if(nrad_in.gt.0) write(cval,'(i0)') nrad_in
        if(nang_in.gt.0) then
          write(cnum,'(i0)') nang_in
        else
          cnum='default'
        end if
        call print_value('One-electron grid (# GRID)',
     +    trim(cval)//' x '//trim(cnum))
      end if

!! ----------------------------------------------------------------- !!
!! # ENPART -- only if ENPART itself is active                       !!
!! ----------------------------------------------------------------- !!
      if(ienpart.eq.1) then
        write(*,*)
        write(*,'(2x,a)') '# ENPART'
        if(ihf.eq.1) then
          cval='HF'
        else if(icas.eq.1) then
          cval='CASSCF'
        else if(icisd.eq.1) then
          cval='CISD'
        else
          cval=' '
          if(id_xcfunc.ne.0) then
            write(cbuf,'(a,i0)') 'libxc ',id_xcfunc
          else
            write(cbuf,'(a,i0,a,i0)') 'libxc exchange ',id_xfunc,
     +        ', correlation ',id_cfunc
            if(id_cfunc.eq.0) write(cbuf,'(a,i0)') 'libxc exchange ',
     +        id_xfunc
          end if
          if(len_trim(cfuncname).gt.0) then
            cval=trim(cfuncname)//' ('//trim(cbuf)//')'
          else
            cval=trim(cbuf)
          end if
        end if
        call print_value('Functional',trim(cval))

!! THREBOD is given in units of 1e-4 of a bond order; below 1, every !!
!! atom pair is integrated                                          !!
        write(cnum,'(i0)') ithrebod
        if(ithrebod.lt.1) then
          cval='0 (THREBOD '//trim(cnum)//': every atom pair)'
        else
          write(cval,'(f6.4)') real(ithrebod)/10000.0d0
          cval=trim(adjustl(cval))//' (THREBOD '//trim(cnum)//')'
        end if
        call print_value('Atom-pair threshold (bond order)',
     +    trim(cval))
        call grid2e_summary()

        cbuf=' '
        if(iexact.eq.1) cbuf=trim(cbuf)//' EXACT'
        if(ihomo.eq.1) cbuf=trim(cbuf)//' HOMO'
        if(idek.eq.1) cbuf=trim(cbuf)//' DEKIN'
        if(iionic.eq.1) cbuf=trim(cbuf)//' IONIC'
        if(ianalytical.eq.1) cbuf=trim(cbuf)//' ANALYTIC'
        if(iecorr.eq.1) cbuf=trim(cbuf)//' CORRELATION'
        if(len_trim(cbuf).gt.0) call print_value('Other options',
     +    trim(adjustl(cbuf)))
      end if

!! ----------------------------------------------------------------- !!
!! # QTAIM -- only if QTAIM is active                                !!
!! ----------------------------------------------------------------- !!
      if(iqtaim.ge.1) then
        write(*,*)
        write(*,'(2x,a)') '# QTAIM'
        call print_value('Convergence step',istep)
        call print_value('Non-nuclear attractors',inna)
        call print_value('Max basin distance',imaxdist)
        call print_value('Screening',iscreening)
        call print_value('Gradient path',ipath)
      end if

!! ----------------------------------------------------------------- !!
!! # CUBE -- only if CUBE is active; thresholds are net occupations   !!
!! times 1000 in the input                                           !!
!! ----------------------------------------------------------------- !!
      if(icube.eq.1) then
        write(*,*)
        write(*,'(2x,a)') '# CUBE'
        write(cnum,'(i0)') jcubthr
        write(cval,'(f6.3)') real(jcubthr)/1000.0d0
        call print_value('Largest net occupation',
     +    trim(adjustl(cval))//' (MAX_OCC '//trim(cnum)//')')
        write(cnum,'(i0)') kcubthr
        write(cval,'(f6.3)') real(kcubthr)/1000.0d0
        call print_value('Smallest net occupation',
     +    trim(adjustl(cval))//' (MIN_OCC '//trim(cnum)//')')
        write(cval,'(f8.3)') cubespacing
        call print_value('Grid spacing',trim(adjustl(cval)),'bohr')
        write(cval,'(f8.3)') cuberadscale
        call print_value('Radius scale',trim(adjustl(cval)))
        if(inegefos.eq.1) then
          write(cnum,'(i0)') inegcubthr
          write(cval,'(f7.3)') -real(inegcubthr)/1000.0d0
          call print_value('Negative-occupation EFO cubes',
     +      'net occupation <= '//trim(adjustl(cval))//' (NEG_EFOS '//
     +      trim(cnum)//', GEOS paired only)')
        end if
      end if

!! ----------------------------------------------------------------- !!
!! # OSLO -- only if OSLO is active; FOLI_TOLERANCE n means 10**(-n)  !!
!! ----------------------------------------------------------------- !!
      if(ioslo.eq.1) then
        write(*,*)
        write(*,'(2x,a)') '# OSLO'
        if(ilow2.eq.1) cval='MULLIKEN (Hilbert space)'
        if(ilow2.eq.2) cval='LOWDIN (Hilbert space)'
        if(ilow2.eq.3) cval='LOWDIN-DAVIDSON (Hilbert space)'
        if(ilow2.eq.6) cval='NAO-BASIS (Hilbert space)'
        if(ilow2.eq.0) cval='real space (as in # METHOD)'
        call print_value('Overlap matrix',trim(cval))
        write(cnum,'(i0)') ifolitol
        write(cfmt,'(a,i0,a,i0,a)') '(f',ifolitol+3,'.',ifolitol,')'
        write(cval,cfmt) 10.0d0**(-ifolitol)
        call print_value('FOLI tolerance',trim(adjustl(cval))//
     +    ' (FOLI_TOLERANCE '//trim(cnum)//')')
        if(nbranch(1)+nbranch(2)+nbranch(3).eq.0)
     +    call print_value('Branching','none')
        if(nbranch(1).gt.0) then
          write(cbuf,'(20(i0,1x))') (ibranchit(k,1),k=1,nbranch(1))
          call print_value('Branching at iteration',trim(cbuf))
        end if
        if(nbranch(2).gt.0) then
          write(cbuf,'(20(i0,1x))') (ibranchit(k,2),k=1,nbranch(2))
          call print_value('Branching at iteration (alpha)',trim(cbuf))
        end if
        if(nbranch(3).gt.0) then
          write(cbuf,'(20(i0,1x))') (ibranchit(k,3),k=1,nbranch(3))
          call print_value('Branching at iteration (beta)',trim(cbuf))
        end if
        if(ioslofchk.eq.2) then
          cval='yes'
        else
          cval='no'
        end if
        call print_value('Non-orthogonalized OSLOs (.fchk)',trim(cval))
      end if

!! ----------------------------------------------------------------- !!
!! # FRAGMENTS -- only if DOFRAGS is active; more than 20 atoms wrap  !!
!! onto aligned continuation lines                                   !!
!! ----------------------------------------------------------------- !!
      if(idofr.eq.1) then
        write(*,*)
        write(cval,'(a,i0,a)') '# FRAGMENTS (',icufr,' fragments)'
        write(*,'(2x,a)') trim(cval)
        do i=1,icufr
          write(*,'(2x,a,i3,a,20i4:/(15x,20i4))') 'Fragment',i,' :',
     +      (ifrlist(k,i),k=1,nfrlist(i))
        end do
      end if

!! ----------------------------------------------------------------- !!
!! # ATOMS -- only if DOATOMS is active                              !!
!! ----------------------------------------------------------------- !!
      if(idoat.eq.1) then
        write(*,*)
        write(cval,'(a,i0,a)') '# ATOMS (',icuat,' selected)'
        write(*,'(2x,a)') trim(cval)
        write(*,'(2x,20i4)') (iatlist(i),i=1,icuat)
      end if

!! ----------------------------------------------------------------- !!
!! # DM -- only if correlated-wavefunction density matrices are read  !!
!! ----------------------------------------------------------------- !!
      if(icorr.ne.0) then
        write(*,*)
        write(*,'(2x,a)') '# DM'
        write(cval,'(i0,a)') icorr,'-RDM'
        call print_value('Density matrices',trim(cval))
        if(ipyscf.eq.1) then
          cval='pySCF'
        else if(iorca.eq.1) then
          cval='ORCA'
        else
          cval='DMn program (binary)'
        end if
        call print_value('Source',trim(cval))
        call print_value('1-RDM file',trim(dmfile1))
        if(icorr.eq.2) call print_value('2-RDM file',trim(dmfile2))
      end if

      END SUBROUTINE print_input_summary

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: atomic_definition_names
!! Names of the atomic definition the run uses (read_input allows only
!! one), from the # METHOD flags. TFVC also sets the BECKE-RHO and NEWBEC
!! flags, so it is tested before them.
!! arguments:
!!   ckey   (out)  input keyword: 'TFVC', 'MULLIKEN', ...
!!   cname  (out)  descriptive name: 'Topological fuzzy Voronoi cells'
!!   cshort (out)  column heading, at most 11 characters: 'TFVC',
!!                 'Lowdin-Dav.', ...
!!   cspace (out)  'real' or 'hilbert'
!! author: MGimf
!! ----------------------------------------------------------------------
      SUBROUTINE atomic_definition_names(ckey,cname,cshort,cspace)
      use input_options_mod, only: imulli,ihirsh,itfvc,inewbec,ibecke
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /achi/achi(maxat,maxat),ibcp
      character*(*) ckey,cname,cshort,cspace

      cspace='real'
      ckey=' '
      cname=' '
      cshort=' '
      if(imulli.ge.1) then
        cspace='hilbert'
        if(imulli.eq.1) then
          ckey='MULLIKEN'
          cname='Mulliken'
          cshort='Mulliken'
        else if(imulli.eq.2) then
          ckey='LOWDIN'
          cname='Lowdin'
          cshort='Lowdin'
        else if(imulli.eq.3) then
          ckey='LOWDIN-DAVIDSON'
          cname='Davidson-Lowdin'
          cshort='Lowdin-Dav.'
        else if(imulli.eq.4) then
          ckey='NAO-BASIS'
          cname='Natural atomic orbitals'
          cshort='NAO'
        else if(imulli.eq.5) then
          ckey='LOWDIN-W'
          cname='Weighted Lowdin'
          cshort='Lowdin-W'
        end if
      else if(ihirsh.eq.1) then
        ckey='HIRSH'
        cname='Hirshfeld'
        cshort='Hirshfeld'
      else if(ihirsh.eq.2) then
        ckey='HIRSH-IT'
        cname='Iterative Hirshfeld'
        cshort='Hirshfeld-I'
      else if(itfvc.eq.1) then
        ckey='TFVC'
        cname='Topological fuzzy Voronoi cells'
        cshort='TFVC'
      else if(ibcp.eq.1) then
        ckey='BECKE-RHO'
        cname='Becke, radii from the density'
        cshort='Becke-rho'
      else if(inewbec.eq.1) then
        ckey='NEWBEC'
        cname='Becke, NEWBEC radii'
        cshort='Becke-new'
      else if(ibecke.eq.1) then
        ckey='BECKE'
        cname='Becke, fixed empirical radii'
        cshort='Becke'
      end if
      END SUBROUTINE atomic_definition_names

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: json_skeleton
!! Writes the input, molecule and fragments blocks of the results file
!! (json_mod) and opens its analyses object, which each analysis then
!! fills. Called from main.f once read_input has set every flag.
!! arguments: none (input_options_mod, COMMONs, the .inp file on unit 16)
!! author: MGimf
!! ----------------------------------------------------------------------
      SUBROUTINE json_skeleton()
      use input_options_mod
      use json_mod
      use elements_mod
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,
     +  jfrlist(maxat)
      character(len=1024) :: line
      character(len=:), allocatable :: text
      character(len=2) :: sym(maxat)
      character(len=20) :: ckey
      character(len=60) :: cname
      character(len=12) :: cshort
      character(len=8) :: cspace
      integer :: ios,i,ifr

!! the input file as written !!
      rewind(16)
      text=''
      do
        read(16,'(a)',iostat=ios) line
        if(ios.ne.0) exit
        text=text//trim(line)//new_line('a')
      end do
      rewind(16)
      call json_begin_object('input')
      call json_value('text',text)

      call atomic_definition_names(ckey,cname,cshort,cspace)
      call json_begin_object('atomic_definition')
      call json_value('keyword',trim(ckey))
      call json_value('name',trim(cname))
      call json_value('short_name',trim(cshort))
      call json_value('space',trim(cspace))
      call json_value('default',itfvcdef.eq.1)
      call json_end()
      call json_end()

!! geometry (bohr) and nuclear charges (reduced by an ECP) !!
      call json_begin_object('molecule')
      call json_value('n_atoms',nat)
      call json_vector('atomic_numbers',iznuc,nat)
      do i=1,nat
        sym(i)=adjustl(element_symbol(iznuc(i)))
      end do
      call json_vector('symbols',sym,nat)
      call json_vector('nuclear_charges',zn,nat)
      call json_begin_array('coordinates_bohr')
      do i=1,nat
        call json_vector('',coord(1:3,i),3)
      end do
      call json_end()
      call json_end()

!! fragments as given in # FRAGMENTS (empty without DOFRAGS) !!
      call json_begin_array('fragments')
      if(idofr.eq.1) then
        do ifr=1,icufr
          call json_begin_object('')
          call json_value('id',ifr)
          call json_vector('atoms',ifrlist(1:nfrlist(ifr),ifr),
     +      nfrlist(ifr))
          call json_end()
        end do
      end if
      call json_end()

      call json_begin_object('analyses')
      END SUBROUTINE json_skeleton

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: grid2e_summary
!! INPUT SUMMARY lines of the two-electron grid (# GRID's RADIAL_2E/
!! ANGULAR_2E or the default) and of its rotation angles, with where
!! they come from.
!! arguments: none (/modgrid/, input_options_mod)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine grid2e_summary()
      use input_options_mod, only: igrid2e,irot2e
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      character*60 cval
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22

      if(igrid2e.eq.1) then
        write(cval,'(i0,a,i0,a)') nrad22,' x ',nang22,' (# GRID)'
      else
        write(cval,'(i0,a,i0,a)') nrad22,' x ',nang22,' (default)'
      end if
      call print_value('Two-electron grid',trim(cval))
      if(irot2e.eq.1) then
        write(cval,'(f5.3,a,f5.3,a)') phb12,' / ',phb22,' (ROTATION_2E)'
      else if(irot2e.eq.2) then
        write(cval,'(f5.3,a,f5.3,a)') phb12,' / ',phb22,
     +    ' (not calibrated for this grid)'
      else
        write(cval,'(f5.3,a,f5.3,a)') phb12,' / ',phb22,' (calibrated)'
      end if
      call print_value('Rotation angles (two-electron grid)',trim(cval))
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: VPRINT
!! Rows of a per-atom table: atom number, element and jdim values, in the
!! columns of the matrix printer (MPRINT2). 6 decimals, 13 with
!! FULLPRECISION. The caller prints the header and the rules.
!! arguments:
!!   H    (in)  (ndim,jdim) values, one row per atom
!!   N    (in)  number of atoms
!!   ndim (in)  leading dimension of H
!!   jdim (in)  number of value columns (at most 6)
!! author:
!! ----------------------------------------------------------------------
      SUBROUTINE VPRINT(H,N,ndim,jdim)
      use elements_mod
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim,jdim
      common /printout/iaccur
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      DIMENSION H(NDIM,jdim)

      do i=1,N
        if(iaccur.eq.0) then
          write(*,'(2x,i3,a4,6f12.6)') i,element_label(iznuc(i)),
     +      (H(i,j),j=1,jdim)
        else
          write(*,'(2x,i3,a4,6f20.13)') i,element_label(iznuc(i)),
     +      (H(i,j),j=1,jdim)
        end if
      end do
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
      use elements_mod
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim
      common /printout/iaccur
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      DIMENSION H(NDIM,3)


      NMIN=1
      NNMAX=3
      PRINT *,'_________________X___________________Y__________________Z__________'
      write(*,*)
      DO I=1,N
        if(iaccur.eq.0) then
          PRINT 62,I,element_label(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
        else
          PRINT 63,I,element_label(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
        end if
      end do
      RETURN
   62 FORMAT(1X,I3,A4,6F12.6)
   63 FORMAT(1X,I3,A4,3F20.13)
      END

!! ----------------------------------------------------------------------
!! subroutine: MPRINT
!! Atom x atom matrix, printed exactly as MPRINT2 (kept so that its
!! callers need no change).
!! arguments:
!!   H    (in)  (ndim,ndim) matrix
!!   N    (in)  number of atoms
!!   ndim (in)  leading dimension of H
!! author:
!! ----------------------------------------------------------------------
      SUBROUTINE MPRINT(H,N,ndim)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, intent(in) :: ndim
      DIMENSION H(NDIM,NDIM)

      call MPRINT2(H,N,ndim)
      END

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: MPRINT2
!! Atom x atom matrix in blocks of 6 columns: a rule, the column header
!! (atom number and element, ending over the last digit of its column),
!! a rule, one row per atom, a rule; rules as wide as the block. 6
!! decimals in 12-character columns, 13 in 20-character columns with
!! FULLPRECISION.
!! arguments:
!!   H    (in)  (ndim,ndim) matrix
!!   N    (in)  number of atoms
!!   ndim (in)  leading dimension of H
!! author:
!! ----------------------------------------------------------------------
      SUBROUTINE MPRINT2(H,N,ndim)
      use elements_mod
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim
      common /printout/iaccur
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      DIMENSION H(NDIM,NDIM)
      character*40 fmthdr,fmtrow
      integer iw,j0,j1,ncol

!! row label "  1  O " (7 characters after the 2-blank margin) !!
      if(iaccur.eq.0) then
        iw=12
        fmtrow='(2x,i3,a4,6f12.6)'
      else
        iw=20
        fmtrow='(2x,i3,a4,6f20.13)'
      end if

      do j0=1,N,6
        j1=min(N,j0+5)
        ncol=j1-j0+1
        if(j0.gt.1) write(*,*)
        write(fmthdr,'(a,i0,a,i0,a)') '(9x,',ncol,'(',iw-6,
     +    'x,i3,1x,a2))'
        write(*,'(2x,a)') repeat('-',7+iw*ncol)
        write(*,fmthdr) (j,element_symbol(iznuc(j)),j=j0,j1)
        write(*,'(2x,a)') repeat('-',7+iw*ncol)
        do i=1,N
          write(*,fmtrow) i,element_label(iznuc(i)),(H(i,j),j=j0,j1)
        end do
        write(*,'(2x,a)') repeat('-',7+iw*ncol)
      end do
      END

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: MPRINTNOAT
!! Fragment x fragment matrix, laid out as MPRINT2 with the fragment
!! number as the only label. Used by group_by_frag_mat.
!! arguments:
!!   H         (in)  (mdim,ndim) matrix
!!   M, N      (in)  rows and columns printed
!!   mdim,ndim (in)  dimensions of H
!! author: MGimf
!! ----------------------------------------------------------------------
      SUBROUTINE MPRINTNOAT(H,M,N,mdim,ndim)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, intent(in) :: ndim,mdim
      DIMENSION H(MDIM,NDIM)
      common /printout/iaccur
      character*40 fmthdr,fmtrow
      integer iw,j0,j1,ncol

!! row label: the fragment number in 3 characters and 6 blanks !!
      if(iaccur.eq.0) then
        iw=12
        fmtrow='(2x,i3,6x,6f12.6)'
      else
        iw=20
        fmtrow='(2x,i3,6x,6f20.13)'
      end if

      do j0=1,N,6
        j1=min(N,j0+5)
        ncol=j1-j0+1
        if(j0.gt.1) write(*,*)
        write(fmthdr,'(a,i0,a,i0,a)') '(11x,',ncol,'i',iw,')'
        write(*,'(2x,a)') repeat('-',9+iw*ncol)
        write(*,fmthdr) (j,j=j0,j1)
        write(*,'(2x,a)') repeat('-',9+iw*ncol)
        do i=1,M
          write(*,fmtrow) i,(H(i,j),j=j0,j1)
        end do
        write(*,'(2x,a)') repeat('-',9+iw*ncol)
      end do
      END

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: frag_sum_mat
!! Fragment x fragment matrix from an atom x atom one, with the atom-to-
!! fragment map of /frlist/ (atoms in no fragment are left out).
!! arguments:
!!   ilog (in)   0: every element is added (B is A summed by blocks);
!!               1: A is symmetric with its total over the diagonal and
!!               one triangle (bond orders, IQA terms), so each atom pair
!!               is counted once and B is symmetric in the same sense
!!   A    (in)   (maxat,maxat) atom x atom matrix
!!   B    (out)  (maxat,maxat) fragment x fragment matrix
!!   x    (out)  total of B: all elements (ilog=0), or the diagonal and
!!               one triangle (ilog=1)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine frag_sum_mat(ilog,A,B,x)
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,
     +  jfrlist(maxat)
      dimension A(maxat,maxat),B(maxat,maxat)
      integer ilog

      do i=1,nat
        do j=1,nat
          B(i,j)=0.0d0
        end do
      end do

!! ilog=1: atom pairs from the lower triangle and the diagonal only; a !!
!! pair whose fragments come in the other order lands above the        !!
!! diagonal of B and is folded back below                              !!
      nnat=nat
      do i=1,nat
        if(ilog.eq.1) nnat=i
        do j=1,nnat
          if(jfrlist(i).ne.0.and.jfrlist(j).ne.0)
     +      B(jfrlist(i),jfrlist(j))=B(jfrlist(i),jfrlist(j))+A(i,j)
        end do
      end do

      x=0.0d0
      if(ilog.eq.1) then
        do i=1,icufr
          x=x+B(i,i)
          do j=1,i-1
            B(i,j)=B(i,j)+B(j,i)
            x=x+B(i,j)
            B(j,i)=B(i,j)
          end do
        end do
      else
        do i=1,icufr
          do j=1,icufr
            x=x+B(i,j)
          end do
        end do
      end if
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: frag_sum_vec
!! Per-fragment table from a per-atom one (sum over the atoms of each
!! fragment, /frlist/ map; atoms in no fragment are left out).
!! arguments:
!!   ndim (in)   number of columns
!!   A    (in)   (maxat,ndim) per-atom table
!!   B    (out)  (maxat,ndim) per-fragment table
!!   x    (out)  (ndim) column totals over the fragments
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine frag_sum_vec(ndim,A,B,x)
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,
     +  jfrlist(maxat)
      dimension A(maxat,ndim),B(maxat,ndim),x(ndim)

      do i=1,nat
        do j=1,ndim
          B(i,j)=0.0d0
        end do
      end do

      do j=1,ndim
        do i=1,nat
          if(jfrlist(i).ne.0) B(jfrlist(i),j)=B(jfrlist(i),j)+A(i,j)
        end do
      end do

      do j=1,ndim
        x(j)=0.0d0
        do i=1,icufr
          x(j)=x(j)+B(i,j)
        end do
      end do
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: group_by_frag_mat
!! Prints the fragment x fragment matrix of an atom x atom quantity
!! (frag_sum_mat) under its own section title, with its total.
!! arguments:
!!   ilog  (in)  0: total over all elements; 1: symmetric quantity,
!!               total over the diagonal and one triangle (frag_sum_mat)
!!   title (in)  section title (leading blanks are ignored)
!!   A     (in)  (maxat,maxat) atom x atom matrix
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine group_by_frag_mat(ilog,title,A)
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,
     +  jfrlist(maxat)
      common /printout/iaccur
      dimension A(maxat,maxat),B(maxat,maxat)
      character*(*) title
      integer ilog

      call frag_sum_mat(ilog,A,B,x)
      call print_box(trim(adjustl(title)))
      call mprintnoat(B,icufr,icufr,maxat,maxat)
      ndec=6
      if(iaccur.ne.0) ndec=13
      if(ilog.eq.1) then
        call print_value('Total (diagonal and one triangle)',x,ndec)
      else
        call print_value('Total (all elements)',x,ndec)
      end if
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: group_by_frag_vec
!! Prints the per-fragment table of a per-atom quantity (frag_sum_vec)
!! under its own section title: one row per fragment and a Total row.
!! arguments:
!!   ndim  (in)  number of value columns (1 or 2)
!!   title (in)  section title (leading blanks are ignored)
!!   A     (in)  (maxat,ndim) per-atom table
!!   h1,h2 (in)  column headings (h2 unused when ndim=1)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine group_by_frag_vec(ndim,title,A,h1,h2)
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      character*(*) title,h1,h2
      dimension A(maxat,ndim),B(maxat,ndim),x(ndim)

      call frag_sum_vec(ndim,A,B,x)
      call print_frag_table(ndim,title,B,x,h1,h2)
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: print_frag_table
!! Prints a per-fragment table under its own section title: one row per
!! fragment and a Total row with the given column totals.
!! arguments:
!!   ndim  (in)  number of value columns (1 or 2)
!!   title (in)  section title (leading blanks are ignored)
!!   B     (in)  (maxat,ndim) per-fragment table
!!   x     (in)  (ndim) column totals
!!   h1,h2 (in)  column headings (h2 unused when ndim=1)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine print_frag_table(ndim,title,B,x,h1,h2)
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,
     +  jfrlist(maxat)
      common /printout/iaccur
      character*(*) title,h1,h2
      dimension B(maxat,ndim),x(ndim)
      integer colw
      character*20 fmtrow,fmttot,fmthdr

!! fragment number in 3 characters and 6 blanks, as in mprintnoat; value !!
!! columns 12 wide (6 decimals), 20 wide (13) with FULLPRECISION         !!
      if(iaccur.eq.0) then
        colw=12
        fmtrow='(2x,i3,6x,2f12.6)'
        fmttot='(2x,a,2f12.6)'
        fmthdr='(2x,a,2a12)'
      else
        colw=20
        fmtrow='(2x,i3,6x,2f20.13)'
        fmttot='(2x,a,2f20.13)'
        fmthdr='(2x,a,2a20)'
      end if

      call print_box(trim(adjustl(title)))
      if(ndim.eq.2) then
        write(*,fmthdr) 'Frag.    ',trim(h1),trim(h2)
      else
        write(*,fmthdr) 'Frag.    ',trim(h1)
      end if
      write(*,'(2x,a)') repeat('-',9+colw*ndim)
      do i=1,icufr
        write(*,fmtrow) i,(B(i,j),j=1,ndim)
      end do
      write(*,'(2x,a)') repeat('-',9+colw*ndim)
      write(*,fmttot) 'Total    ',(x(j),j=1,ndim)
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: json_atom_matrix
!! Writes an atom x atom matrix to the results file as a matrix object:
!! values (full square), symmetric, what the diagonal and off-diagonal
!! elements are, and the total when the output prints one.
!! arguments:
!!   key   (in)  key of the object
!!   A     (in)  (maxat,maxat) matrix
!!   lsym  (in)  logical, the matrix is symmetric
!!   cdiag (in)  meaning of the diagonal elements
!!   coff  (in)  meaning of the off-diagonal elements
!!   itot  (in)  -1: no total; 0: total over all elements; 1: over the
!!               diagonal and one triangle
!!   xtot  (in)  the total (unused when itot=-1)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_atom_matrix(key,A,lsym,cdiag,coff,itot,xtot)
      use json_mod
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      character*(*) key,cdiag,coff
      dimension A(maxat,maxat)
      logical lsym

      call json_begin_object(key)
      call json_matrix('values',A,nat,nat,maxat)
      call json_value('symmetric',lsym)
      call json_value('diagonal',cdiag)
      call json_value('off_diagonal',coff)
      if(itot.ge.0) then
        call json_value('total',xtot)
        if(itot.eq.1) then
          call json_value('total_rule','diagonal and one triangle')
        else
          call json_value('total_rule','all elements')
        end if
      end if
      call json_end()
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: json_frag_matrix
!! Writes the fragment x fragment matrix of an atom x atom quantity
!! (frag_sum_mat) as a matrix object with its total, the values that
!! group_by_frag_mat prints.
!! arguments:
!!   key   (in)  key of the object
!!   ilog  (in)  as in frag_sum_mat
!!   A     (in)  (maxat,maxat) atom x atom matrix
!!   lsym  (in)  logical, the fragment matrix is symmetric
!!   cdiag (in)  meaning of the diagonal elements
!!   coff  (in)  meaning of the off-diagonal elements
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_frag_matrix(key,ilog,A,lsym,cdiag,coff)
      use json_mod
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,
     +  jfrlist(maxat)
      character*(*) key,cdiag,coff
      dimension A(maxat,maxat),B(maxat,maxat)
      logical lsym

      call frag_sum_mat(ilog,A,B,x)
      call json_begin_object(key)
      call json_matrix('values',B,icufr,icufr,maxat)
      call json_value('symmetric',lsym)
      call json_value('diagonal',cdiag)
      call json_value('off_diagonal',coff)
      call json_value('total',x)
      if(ilog.eq.1) then
        call json_value('total_rule','diagonal and one triangle')
      else
        call json_value('total_rule','all elements')
      end if
      call json_end()
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: json_frag_vector
!! Writes the per-fragment values of a per-atom quantity (frag_sum_vec),
!! the values that group_by_frag_vec prints.
!! arguments:
!!   key (in)  key of the list
!!   A   (in)  (maxat) per-atom values
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_frag_vector(key,A)
      use json_mod
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,
     +  jfrlist(maxat)
      character*(*) key
      dimension A(maxat),B(maxat),x(1)

      call frag_sum_vec(1,A,B,x)
      call json_vector(key,B,icufr)
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: print_efo_block
!! Prints the EFOs of one fragment (or atom): the net occupation of the
!! EFOs above the cutoff and of each of them, the same for the gross
!! occupations and the fit of each reprojected EFO when the scheme has
!! them, and what the cutoff leaves out. Shared by EFFAO, EOS and GEOS.
!! No empty line before or after: the caller separates the blocks.
!! arguments:
!!   clabel  (in)  'Fragment' or 'Atom'
!!   ifrag   (in)  its number
!!   n       (in)  number of EFOs above the cutoff
!!   xnet    (in)  their total net occupation
!!   occnet  (in)  (n) net occupation of each EFO
!!   lgross  (in)  logical: print the gross occupations
!!   xgro    (in)  their total gross occupation (lgross only)
!!   occgro  (in)  (n) gross occupation of each EFO (lgross only)
!!   lfit    (in)  logical: print the fit of each EFO
!!   fitpct  (in)  (n) fit of each EFO, in % (lfit only)
!!   xleftn  (in)  net occupation left out by the cutoff
!!   xleftg  (in)  gross occupation left out by the cutoff (lgross only)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine print_efo_block(clabel,ifrag,n,xnet,occnet,lgross,
     +  xgro,occgro,lfit,fitpct,xleftn,xleftg)
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      character*(*) clabel
      dimension occnet(*),occgro(*),fitpct(*)
      logical lgross,lfit

      write(*,'(2x,a,1x,i0)') clabel,ifrag
      call print_value('Net occupation',xnet,5)
      if(n.gt.0) then
        call print_efo_row('Net',occnet,n,4)
      else
        write(*,'(2x,a)') 'No EFO above the cutoff'
      end if
      if(lgross) then
        call print_value('Gross occupation',xgro,5)
        if(n.gt.0) call print_efo_row('Gross',occgro,n,4)
      end if
      if(lfit.and.n.gt.0) call print_efo_row('Fit (%)',fitpct,n,2)
      call print_value('Left out by the cutoff (net)',xleftn,5)
      if(lgross) call print_value('Left out by the cutoff (gross)',
     +  xleftg,5)
      end

!! ---- !!
!! one value per EFO, 8 per line, the label on the first line only     !!
      subroutine print_efo_row(clabel,v,n,ndec)
      implicit real*8(a-h,o-z)
      character*(*) clabel
      dimension v(n)
      character*20 fmt1,fmt2
      character*7 lab

      write(fmt1,'(a,i0,a)') '(2x,a7,8f9.',ndec,')'
      write(fmt2,'(a,i0,a)') '(9x,8f9.',ndec,')'
      lab=clabel
      write(*,fmt1) lab,(v(i),i=1,min(n,8))
      do i0=9,n,8
        write(*,fmt2) (v(i),i=i0,min(n,i0+7))
      end do
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: print_eos_frontier
!! Prints the result of the electron assignment of one spin channel:
!! the electrons of each fragment, the gross occupations of its last
!! occupied and first unoccupied EFOs ('-' when it has none), and the
!! reliability index R. Shared by EOS and GEOS.
!! arguments:
!!   n     (in)  number of fragments
!!   elec  (in)  (n) electrons of each fragment
!!   xlo   (in)  (n) occupation of its last occupied EFO
!!   ilo   (in)  (n) 1 when that EFO exists
!!   xfu   (in)  (n) occupation of its first unoccupied EFO
!!   ifu   (in)  (n) 1 when that EFO exists
!!   rel   (in)  reliability index R, in %
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine print_eos_frontier(n,elec,xlo,ilo,xfu,ifu,rel)
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      dimension elec(n),xlo(n),xfu(n),ilo(n),ifu(n)
      character*11 clo
      character*14 cfu

      write(*,'(2x,a)') 'Frag.  Electrons  Last occ.  First unocc.'
      write(*,'(2x,a)') repeat('-',41)
      do i=1,n
        if(ilo(i).eq.1) then
          write(clo,'(f11.3)') xlo(i)
        else
          clo='          -'
        end if
        if(ifu(i).eq.1) then
          write(cfu,'(f14.3)') xfu(i)
        else
          cfu='             -'
        end if
        write(*,'(2x,i5,f11.2,a11,a14)') i,elec(i),clo,cfu
      end do
      write(*,'(2x,a)') repeat('-',41)
      call print_value('Reliability index R (%)',rel,3)
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: print_oxidation_states
!! Prints the section with the oxidation state of each fragment and
!! their total, and the overall reliability index when the method has
!! one. Shared by EOS, GEOS and OSLO.
!! arguments:
!!   n     (in)  number of fragments
!!   oxi   (in)  (n) oxidation states
!!   xtot  (in)  their total
!!   lrel  (in)  logical: print the overall reliability index
!!   rel   (in)  overall reliability index R, in % (lrel only)
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine print_oxidation_states(n,oxi,xtot,lrel,rel)
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      dimension oxi(n)
      logical lrel

      call print_box('FRAGMENT OXIDATION STATES')
      write(*,'(2x,a)') 'Frag.  Oxidation state'
      write(*,'(2x,a)') repeat('-',22)
      do i=1,n
        write(*,'(2x,i5,f17.2)') i,oxi(i)
      end do
      write(*,'(2x,a)') repeat('-',22)
      write(*,'(2x,a5,f17.2)') 'Total',xtot
      if(lrel) call print_value('Overall reliability index R (%)',rel,
     +  3)
      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: print_int
!! DOINT: writes the atomic overlap matrix of each atom in the basis of
!! the occupied (or natural) orbitals, one AIMPAC-format .int file per
!! atom plus a .files list, for programs such as ESI-3D. With every atom,
!! checks that the matrices add up to the identity.
!! arguments:
!!   nbas,nat0 (in)  declared dimensions of sat (basis functions, atoms)
!!   sat       (in)  (nbas,nbas,nat0) atomic overlap matrices, AO basis
!!   name      (in)  job name, start of the file names
!! author:
!! ----------------------------------------------------------------------
      subroutine print_int(nbas,nat0,sat,name)
      use elements_mod
      use output_mod
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

      character nameaim*256,charnu*2,atnum*12,ext*4,namefiles*256
      character*40 cval
      logical lfound

      allocatable c3(:,:),c2(:,:),csave(:,:),scr(:,:)

      ihirsh = Iopt(6)
      imulli=Iopt(5)
      icorr=Iopt(26)
      iqtaim =Iopt(16)

      inato=0
      if(icorr.eq.1.or.icas.eq.1.or.icisd.eq.1) inato=1

!! independent-function count from the .fchk (unit 15); one per basis  !!
!! function when the file doesn't give it                               !!
      ndim=int_locate(15,"Number of independ",lfound)
      if(.not.lfound) ndim=igr

      allocate( c3(2*igr,2*igr),csave(2*igr,2*igr))
      allocate( c2(igr,igr),scr(igr,igr))
      scr=0.0d0
      csave=0.0d0

      do i=1,igr
        do j=1,igr
          scr(i,j)=c(i,j)
        end do
      end do

      if(inato.eq.1) then
        if(icas.eq.1) ndim=norb
        do i=1,igr
          do j=1,ndim
            scr(i,j)=c_no(i,j)
          end do
        end do
      end if

!! the .files list and its extension by atomic definition !!
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
      namefiles=nameaim
      open(file=nameaim,unit=naim3,status="unknown")
      rewind(naim3)
      do jjat=1,icuat
        if(icuat.ne.nat) then
          jat=iatlist(jjat)
        else
          jat=jjat
        end if
        charnu=element_symbol(iznuc(jat))
        charnu=adjustl(charnu)
        write(atnum,'(i0)')jat
        nameaim=trim(name)//ext//trim(charnu)//trim(atnum)//".int"
        write(naim3,'(a)') trim(nameaim)
        open(file=nameaim,unit=naim,status="unknown")
        rewind(naim)

!! the .int layout is that of AIMPAC (PROAIMV), read by other programs: !!
!! keep it as it is. AE and BK are placeholders, not computed energies   !!
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
        write(naim,'(A25,A4,i5)')' INTEGRATION IS OVER ATOM',element_symbol(iznuc(jat)),
     1  jat
      !write(naim2) element_symbol(iznuc(jat)),jat
        write(naim,'(a27)') ' RESULTS OF THE INTEGRATION'
        write(naim,'(a17,e21.14,a14,e21.14)')'              N  ',
     1   qat(jat,1),'    NET CHARGE',zn(jat)-qat(jat,1)
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

!! S^A in the orbital basis, C^T S^A C !!
        numorb=nalf
        if(inato.eq.1) numorb=ndim
        do i=1,numorb
          do j=1,igr
            xx=0.0d0
            do k=1,igr
              xx=xx+scr(k,i)*sat(k,j,jat)
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
!! beta orbitals of an unrestricted wavefunction, after the alpha ones !!
        if(kop.eq.1.and.inato.eq.0)  then
          do i=1,nb
            do j=1,igr
              xx=0.0d0
              do k=1,igr
                xx=xx+cb(k,i)*sat(k,j,jat)
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

!! orbitals written: occupied alpha then beta (unrestricted), occupied !!
!! (restricted) or natural orbitals (correlated)                         !!
        nwrit=nocc
        if(kop.eq.1.and.inato.eq.0) nwrit=nalf+nb
        if(inato.eq.1) nwrit=numorb
        do i=1,nwrit
          do j=1,nwrit
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
!! alpha and beta populations from the diagonal !!
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
          write(naim,'(a41,e21.14)') 'ALPHA ELECTRONS (NA)',qat(jat,1)/2.0d0
          write(naim,'(a41,e21.14)') 'BETA ELECTRONS (NB)',qat(jat,1)/2.0d0
!! correlated wavefunction !!
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
          write(naim,'(a41,e21.14)') 'ALPHA ELECTRONS (NA)',qat(jat,1)/2.0d0
          write(naim,'(a41,e21.14)') 'BETA ELECTRONS (NB)',qat(jat,1)/2.0d0
        end if


        write(naim,*)' '
        write(naim,*) 'NORMAL TERMINATION OF PROAIMV'

        close(naim)
      end do
      close(naim3)
      write(cval,'(i0,a)') icuat,' (one per atom)'
      call print_value('Atomic overlap files (.int)',trim(cval))
      call print_value('List of the .int files',trim(namefiles))
      if(inato.eq.1) call print_value('Natural orbitals in the .int files',
     +  nwrit)

!! with every atom, the atomic overlap matrices add up to the identity: !!
!! largest deviation on the diagonal, off it, and asymmetry             !!
      if(icuat.eq.nat) then
        devd=0.0d0
        devo=0.0d0
        devs=0.0d0
        do i=1,nwrit
          devd=max(devd,abs(csave(i,i)-1.0d0))
          do j=i+1,nwrit
            devo=max(devo,abs(csave(i,j)),abs(csave(j,i)))
            devs=max(devs,abs(csave(i,j)-csave(j,i)))
          end do
        end do
        call print_value('Sum of the overlaps, largest deviation',
     +    max(devd,devo),6)
        if(devd.gt.1.0d-2.or.devo.gt.1.0d-3.or.devs.gt.1.0d-4) then
          write(*,'(2x,a)') 'WARNING: the atomic overlap matrices do'//
     +      ' not add up to the identity (diagonal, off-diagonal,'
          write(*,'(2x,a,3es10.2,a)') 'asymmetry:',devd,devo,devs,
     +      '); check the integration grid'
        end if
      end if
      deallocate(c3,scr,c2,csave)

      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: cubegen_new
!! CUBE: writes one cube file for each EFO of fragment (or atom) ifrag
!! whose net occupation is within MAX_OCC/MIN_OCC, and prints a table
!! with the norm of each cube. Real-space schemes write the EFO times the
!! fragment's weight, Hilbert-space schemes the EFO itself.
!! arguments:
!!   ifrag (in)  fragment (atom without DOFRAGS)
!!   icase (in)  density: 0 closed shell, 1 alpha, 2 beta, 3 paired,
!!               4 unpaired, 5 negative paired EFOs (all written)
!! author: PSalse, MGimf
!! ----------------------------------------------------------------------
      subroutine cubegen_new(ifrag,icase)
      use elements_mod
      use output_mod
      use ao_matrices
      use integration_grid
      use effao_mod, only: p0,p0net,p0gro,ip0
      use input_options_mod, only: cubespacing,cuberadscale
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /atomrad/atr(maxat),dist(maxat,maxat)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /iops/iopt(200)
      character*60 name
      common /filename/name

      character*30 name2
      character nameaim*256,charnu*2,atnu*12,charnu1*12,cwho*20
      character*40 cval
      character*256, allocatable :: cfile(:)

      dimension pgrid(maxgrid,3),igrid(3),xgrid(3,3)
      allocatable c0(:,:),xyz(:,:,:),xnrm(:)

      ihirsh=iopt(6)
      imulli= Iopt(5)
      ibcp=Iopt(14)
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

!! EFOs to plot: icase=5 plots every staged EFO; otherwise the first EFO !!
!! with net occ. < MAX_OCC/1000 through the last with net occ. >         !!
!! MIN_OCC/1000 (doubled for occupations up to 2), both scans capped to  !!
!! this fragment's own 1..imaxo EFOs                                     !!
      if(icase.eq.5) then
        imaxeff=1
        imineff=imaxo
      else
        xmaxeff=float(jcubthr)*1.0d-3
        if (icase.eq.0.or.icase.eq.3) xmaxeff=2.0d0*xmaxeff
        imaxeff=1
        do while(imaxeff.le.imaxo)
          if(p0net(imaxeff,ifrag).lt.xmaxeff) exit
          imaxeff=imaxeff+1
        end do

        xmineff=float(kcubthr)*1.0d-3
        if (icase.eq.0.or.icase.eq.3) xmineff=2.0d0*xmineff
        imineff=imaxo
        do while(imineff.ge.1)
          if(p0net(imineff,ifrag).gt.xmineff) exit
          imineff=imineff-1
        end do
      end if

      write(atnu,'(i0)')ifrag
      if(idofr.eq.0) then
        cwho='atom '//trim(atnu)//' ('//
     +    trim(element_symbol(iznuc(ifrag)))//')'
      else
        cwho='fragment '//trim(atnu)
      end if
      if(imaxeff.gt.imineff) then
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Cube files, '//trim(cwho)//
     +    ': no EFO in the MAX_OCC/MIN_OCC range'
        deallocate(c0)
        return
      end if

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
        else if(ihirsh.eq.1) then
          name2="hirsh"
        else if(ihirsh.eq.2) then
          name2="hirsh-it"
        end if
      end if
      if(icase.eq.3) name2=trim(name2)//"_paired"
      if(icase.eq.4) name2=trim(name2)//"_unpaired"
      if(icase.eq.5) name2=trim(name2)//"_paired_neg"

      xgrid=0.0d0

!! box: the fragment's atoms plus RADIUS_SCALE times the covalent      !!
!! radius of the outermost atom on each side, points SPACING apart     !!
      xmesh=cubespacing
      rrmax=cuberadscale
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
      end do

      allocate ( xyz(igrid(1),igrid(2),igrid(3)))
      allocate(xnrm(imaxeff:imineff),cfile(imaxeff:imineff))

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
                  if(ihirsh.eq.0) then
                    ww=ww+wat(jjat,xabs,yabs,zabs)
                  else if(ihirsh.eq.1) then
                    ww=ww+wathirsh(jjat,xabs,yabs,zabs)
                  else if(ihirsh.eq.2) then
                    ww=ww+wathirsh2(jjat,xabs,yabs,zabs)
                  end if
                end do
                xyz(i,j,k)=orbxyz(c0,ivec,xabs,yabs,zabs)*ww
              end if
            end do
          end do
        end do
!$OMP END PARALLEL DO

!! norm in the cube, sum of the squared values times the voxel volume !!
!! (same independence argument, reduction on x0)                       !!
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
        xnrm(ivec)=x0*xmesh**3

!! file name and the cube file !!

        if(idofr.eq.0) then
          charnu=element_symbol(iznuc(ifrag))
          charnu=adjustl(charnu)
        else
          charnu="FR"
        end if
        write(charnu1,'(i0)')ivec
        if(icase.ne.2) then
          nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)
        else
          nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)//"beta"
        end if
        nameaim=trim(nameaim)//".cube"
        cfile(ivec)=nameaim

        open(44,file=nameaim,status="unknown")
        rewind(44)
        write(44,'(a)') 'Cube generated by APOST-3D'
        write(44,'(a,1x,a,a,i0,a,a,a,f8.4,a,f8.4)') trim(name),
     +    trim(name2),' EFO ',ivec,' of ',trim(cwho),': gross occ.',
     +    p0gro(ivec,ifrag),', net occ.',p0net(ivec,ifrag)
        write(44,42) nat,(pgrid(1,j),j=1,3)
        do i=1,3
          write(44,42) igrid(i),(xgrid(i,j),j=1,3)
        end do
        do i=1,nat
          write(44,43) iznuc(i),zn(i),(coord(j,i),j=1,3)
        end do
        do i=1,igrid(1)
          do j=1,igrid(2)
            write(44,40)(xyz(i,j,k),k=1,igrid(3))
          end do
        end do
        close(44)
      end do

!! one table per call: occupation, norm captured by the box, file      !!
      write(*,'(a)') ''
      if(icase.eq.5) then
        write(*,'(2x,a)') 'Cube files, '//trim(cwho)//
     +    ' (negative paired EFOs)'
      else
        write(*,'(2x,a)') 'Cube files, '//trim(cwho)
      end if
      write(cval,'(i0,a,i0,a,i0)') igrid(1),' x ',igrid(2),' x ',
     +  igrid(3)
      call print_value('Grid points (x, y, z)',trim(cval))
      write(*,'(a)') ''
      lmax=4
      do ivec=imaxeff,imineff
        lmax=max(lmax,len_trim(cfile(ivec)))
      end do
      lrule=34+lmax
      write(*,'(2x,a)') '  EFO   Net occ.   Norm in cube   File'
      write(*,'(2x,a)') repeat('-',lrule)
      nlow=0
      do ivec=imaxeff,imineff
        write(*,'(2x,i5,f11.4,f15.4,3x,a)') ivec,p0net(ivec,ifrag),
     +    xnrm(ivec),trim(cfile(ivec))
        if(xnrm(ivec).lt.0.9d0) nlow=nlow+1
      end do
      write(*,'(2x,a)') repeat('-',lrule)
!! each EFO has norm 1 over all space, so the cube holds the share     !!
!! given in the table                                                  !!
      if(nlow.gt.0) then
        write(*,'(2x,a,i0,a)') 'WARNING: ',nlow,' cube file(s) with a'//
     +    ' norm below 0.9; in # CUBE, raise RADIUS_SCALE'
        write(*,'(2x,a)') '(box too small) or lower SPACING (grid too'//
     +    ' coarse)'
      end if
      deallocate(xnrm,cfile)
42    format(i5,3f12.6)
43    format(i5,4f12.6)
40    format(6e13.5)
     
      deallocate(c0,xyz)
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: skip_fchk_orbblock                                        !!
!! purpose: skips exactly the original data of an Alpha/Beta Orbital     !!
!!   Energies or MO coefficients block just overwritten by the caller,   !!
!!   landing on whatever line follows -- used by rwf_effao_orbprint/     !!
!!   uwf_effao_orbprint in place of searching for the next section's     !!
!!   header. The .fchk format always writes an R-type block as exactly   !!
!!   ceil(nval/5) lines of 5 values each immediately after its header,   !!
!!   regardless of which program wrote the file, so this is robust even  !!
!!   when the source .fchk omits optional sections a marker search would !!
!!   otherwise hang or run off the end of the file looking for (MOKIT/   !!
!!   pySCF-written files routinely omit "Orthonormal basis", for one).   !!
!! arguments:                                                            !!
!!   iunit (in)    -- unit number of the source .fchk, already           !!
!!                    positioned right after the block's header line     !!
!!   nval  (in)    -- number of values in the block (e.g. indepigr or    !!
!!                    igr*indepigr)                                      !!
!!   line  (inout) -- on return, holds the line immediately following    !!
!!                    the skipped block                                  !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine skip_fchk_orbblock(iunit,nval,line)

      implicit double precision(a-h,o-z)
      character*80 line

      nskip=(nval+4)/5
      do ii=1,nskip
        read(iunit,'(a80)') line
      end do
      read(iunit,'(a80)') line

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: skip_fchk_block                                           !!
!! purpose: skips the data of the R-type .fchk block whose header was    !!
!!   just read into line, taking its size from the header's own "N="     !!
!!   field (independent of the program that wrote the file), and lands  !!
!!   on the line that follows the block.                                 !!
!! arguments:                                                            !!
!!   iunit (in)    -- unit number of the source .fchk, positioned right  !!
!!                    after the header line                              !!
!!   line  (inout) -- in: the block's header line; out: the line right   !!
!!                    after the block                                    !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine skip_fchk_block(iunit,line)

      implicit double precision(a-h,o-z)
      character*80 line

      ipos=index(line,"N=")
      read(line(ipos+2:),*) nval
      call skip_fchk_orbblock(iunit,nval,line)

      end

!! ***** !!

!! ----------------------------------------------------------------------
!! subroutine: write_fchk_title
!! Title line of an orbital .fchk file written by APOST-3D: what the
!! orbitals are and the run id of the results file, so that a reader can
!! tell the file belongs to that run (the source title is dropped).
!! arguments:
!!   iu    (in)  unit of the .fchk file being written
!!   ctype (in)  file name suffix: "-EOS-EFOs", "-GEOS-EFOs", "-OSLOs"
!!               or "-OSLOs-preortho"
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine write_fchk_title(iu,ctype)
      use json_mod, only: json_run_id
      implicit none
      integer, intent(in) :: iu
      character(len=*), intent(in) :: ctype
      character(len=60) :: what

      if(trim(ctype).eq.'-EOS-EFOs') then
        what='APOST-3D EOS EFOs'
      else if(trim(ctype).eq.'-GEOS-EFOs') then
        what='APOST-3D GEOS EFOs'
      else if(trim(ctype).eq.'-OSLOs-preortho') then
        what='APOST-3D OSLOs before orthogonalization'
      else
        what='APOST-3D OSLOs'
      end if
      if(len(json_run_id()).gt.0) then
        write(iu,'(a)') trim(what)//', run '//json_run_id()
      else
        write(iu,'(a)') trim(what)
      end if
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: rwf_effao_orbprint                                        !!
!! purpose: writes a restricted-wavefunction pooled-EFO .fchk -- splices !!
!!   the original .fchk's structure, replacing "Alpha Orbital Energies"  !!
!!   and "Alpha MO coefficients" with pooled EFO coefficients (from      !!
!!   ueos.f's ueos_analysis for GEOS, or effao.f's eos_analysis for      !!
!!   plain EOS), sorted by decreasing gross occupation and used in place !!
!!   of a real orbital energy, copying everything else through           !!
!!   unchanged (including Total SCF Density -- these are visualization   !!
!!   orbitals, not a real wavefunction). Reuses the same splice pattern  !!
!!   as OSLO's oslo_fchk_write (oslo.f), extended to also cover the     !!
!!   Energies block. Called whenever GEOS or EOS runs (DOFRAGS, real-    !!
!!   space path) on a restricted wavefunction -- see uwf_effao_orbprint  !!
!!   for the unrestricted twin.                                          !!
!! arguments:                                                            !!
!!   pcoef (in) -- (igr,igr) pooled+sorted EFO coefficients, zero-padded !!
!!                 beyond the actual count                               !!
!!   ctype (in) -- filename suffix, e.g. "-GEOS-EFOs" or "-EOS-EFOs"     !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine rwf_effao_orbprint(pcoef,ctype)

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
      indepigr = int_locate(15,"Number of independ",ilog)
      norb     = igr*indepigr

!! fake orbital energies: gross occupation of each pooled EFO, already   !!
!! sorted decreasing by the caller, zero beyond the actual count.        !!
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
      call write_fchk_title(69,ctype)
      read(15,'(a80)') line

!! standard/MOKIT .fchk layout: Alpha Orbital Energies precedes Alpha    !!
!! MO coefficients. Q-Chem's is the other way around -- see uwf_effao_   !!
!! orbprint's header for the same iqchem branch used by OSLO's own       !!
!! printers.                                                             !!
      if(iqchem.eq.0) then
        do while(index(line,"Alpha Orbital").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do
        write(69,11) "Alpha Orbital Energies","R","N= ",indepigr
        write(69,13) (energ(ii),ii=1,indepigr)

        call skip_fchk_orbblock(15,indepigr,line)
        write(69,12) "Alpha MO coefficients","R","N= ",norb
        write(69,13) ((pcoef(ii,jj),ii=1,igr),jj=1,indepigr)

        call skip_fchk_orbblock(15,norb,line)
      else
        do while(index(line,"Alpha MO co").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do
        write(69,12) "Alpha MO coefficients","R","N= ",norb
        write(69,13) ((pcoef(ii,jj),ii=1,igr),jj=1,indepigr)

        call skip_fchk_orbblock(15,norb,line)
        write(69,11) "Alpha Orbital Energies","R","N= ",indepigr
        write(69,13) (energ(ii),ii=1,indepigr)

        call skip_fchk_orbblock(15,indepigr,line)
      end if

!! copy everything else through unchanged, Total SCF Density included    !!
      do while(.true.)
        write(69,'(a80)') line
        read(15,'(a80)',end=99) line
      end do
99    continue
      close(69)
      DEALLOCATE(energ)

!! Printing formats !!
11    FORMAT(a22,21x,a1,3x,a3,i11)
12    FORMAT(a21,22x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: uwf_effao_orbprint                                        !!
!! purpose: writes an unrestricted-wavefunction pooled-EFO .fchk --      !!
!!   splices the original .fchk's structure, replacing "Alpha/Beta       !!
!!   Orbital Energies" and "Alpha/Beta MO coefficients" with pooled EFO  !!
!!   coefficients (paired/alpha -> Alpha, unpaired/beta -> Beta channel, !!
!!   from ueos.f's ueos_analysis for GEOS or effao.f's eos_analysis for  !!
!!   plain EOS), sorted by decreasing gross occupation and used in place !!
!!   of a real orbital energy, copying everything else through unchanged !!
!!   (including Total SCF/Spin SCF Density -- these are visualization    !!
!!   orbitals, not a real wavefunction). Reuses the same splice pattern  !!
!!   as OSLO's oslo_fchk_write (oslo.f), extended to also cover the     !!
!!   Energies blocks and to skip past every original block by its known  !!
!!   size (skip_fchk_orbblock) rather than searching for the next        !!
!!   section's header -- robust to MOKIT/pySCF-written .fchk files that  !!
!!   omit optional sections like "Orthonormal basis" a search would      !!
!!   otherwise run off the end of the file looking for. Q-Chem-format    !!
!!   branch is still unverified by an active test (no unrestricted       !!
!!   Q-Chem source in the suite), though structurally symmetric with the !!
!!   standard/MOKIT branch now.                                          !!
!! arguments:                                                            !!
!!   pcoef_a (in) -- (igr,igr) pooled+sorted first-channel (-> Alpha)    !!
!!                   EFO coefficients, zero-padded beyond the actual     !!
!!                   count                                               !!
!!   pcoef_b (in) -- (igr,igr) pooled+sorted second-channel (-> Beta)    !!
!!                   EFO coefficients, zero-padded beyond the actual     !!
!!                   count                                               !!
!!   ctype   (in) -- filename suffix, e.g. "-GEOS-EFOs" or "-EOS-EFOs"   !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine uwf_effao_orbprint(pcoef_a,pcoef_b,ctype)

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
      indepigr = int_locate(15,"Number of independ",ilog)
      norb     = igr*indepigr

!! fake orbital energies: gross occupation of each pooled EFO, already   !!
!! sorted decreasing by the caller, zero beyond the actual count.        !!
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
      call write_fchk_title(69,ctype)
      read(15,'(a80)') line

      if(iqchem.eq.0) then

!! standard/MOKIT layout: both Orbital Energies blocks precede both MO   !!
!! coefficient blocks.                                                   !!
        do while(index(line,"Alpha Orbital").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do
        write(69,11) "Alpha Orbital Energies","R","N= ",indepigr
        write(69,13) (energ_a(ii),ii=1,indepigr)

        call skip_fchk_orbblock(15,indepigr,line)
        write(69,11) "Beta Orbital Energies ","R","N= ",indepigr
        write(69,13) (energ_b(ii),ii=1,indepigr)

        call skip_fchk_orbblock(15,indepigr,line)
        write(69,12) "Alpha MO coefficients","R","N= ",norb
        write(69,13) ((pcoef_a(ii,jj),ii=1,igr),jj=1,indepigr)

        call skip_fchk_orbblock(15,norb,line)
        write(69,12) "Beta MO coefficients ","R","N= ",norb
        write(69,13) ((pcoef_b(ii,jj),ii=1,igr),jj=1,indepigr)

        call skip_fchk_orbblock(15,norb,line)
      else

!! Q-Chem layout: both MO coefficient blocks precede both Orbital        !!
!! Energies blocks -- unverified, see header.                            !!
        do while(index(line,"Alpha MO co").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do
        write(69,12) "Alpha MO coefficients","R","N= ",norb
        write(69,13) ((pcoef_a(ii,jj),ii=1,igr),jj=1,indepigr)

        call skip_fchk_orbblock(15,norb,line)
        write(69,12) "Beta MO coefficients ","R","N= ",norb
        write(69,13) ((pcoef_b(ii,jj),ii=1,igr),jj=1,indepigr)

        call skip_fchk_orbblock(15,norb,line)
        write(69,11) "Alpha Orbital Energies","R","N= ",indepigr
        write(69,13) (energ_a(ii),ii=1,indepigr)

        call skip_fchk_orbblock(15,indepigr,line)
        write(69,11) "Beta Orbital Energies ","R","N= ",indepigr
        write(69,13) (energ_b(ii),ii=1,indepigr)

        call skip_fchk_orbblock(15,indepigr,line)
      end if

!! copy everything else through unchanged, Total/Spin SCF Density        !!
!! included                                                              !!
      do while(.true.)
        write(69,'(a80)') line
        read(15,'(a80)',end=99) line
      end do
99    continue
      close(69)
      DEALLOCATE(energ_a,energ_b)

!! Printing formats !!
11    FORMAT(a22,21x,a1,3x,a3,i11)
12    FORMAT(a21,22x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: rwfu_effao_orbprint                                       !!
!! purpose: writes a pooled-EFO .fchk for a restricted-wavefunction      !!
!!   source (kop=0) whose unpaired/second channel is genuinely populated !!
!!   -- e.g. a restricted-orbital post-HF wavefunction (CASSCF/FCI/CISD) !!
!!   with real fractional natural-orbital occupations, not the trivial   !!
!!   ~0 unpaired channel of an ordinary restricted single-determinant    !!
!!   wavefunction (see rwf_effao_orbprint for that case). The source     !!
!!   .fchk has no Beta MO set to splice into, so this *inserts* a        !!
!!   synthetic Beta Orbital Energies/MO coefficients block (second       !!
!!   channel) immediately next to the real Alpha one (first channel),    !!
!!   rather than searching for and overwriting an existing Beta block    !!
!!   the way uwf_effao_orbprint does for a genuinely unrestricted        !!
!!   source -- nothing needs skipping for the inserted blocks since      !!
!!   nothing existed there before; the original Alpha-block data is      !!
!!   still skipped by its known size (skip_fchk_orbblock), same as       !!
!!   rwf_effao_orbprint. Total SCF Density and everything else is        !!
!!   copied through unchanged, exactly as the other two writers.         !!
!! arguments:                                                            !!
!!   pcoef_a (in) -- (igr,igr) pooled+sorted first-channel (-> Alpha)    !!
!!                   EFO coefficients, zero-padded beyond the actual     !!
!!                   count                                               !!
!!   pcoef_b (in) -- (igr,igr) pooled+sorted second-channel (-> Beta,    !!
!!                   synthetic) EFO coefficients, zero-padded beyond     !!
!!                   the actual count                                    !!
!!   ctype   (in) -- filename suffix, e.g. "-GEOS-EFOs" or "-EOS-EFOs"   !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine rwfu_effao_orbprint(pcoef_a,pcoef_b,ctype)

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
      indepigr = int_locate(15,"Number of independ",ilog)
      norb     = igr*indepigr

!! fake orbital energies: gross occupation of each pooled EFO, already   !!
!! sorted decreasing by the caller, zero beyond the actual count.        !!
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
      call write_fchk_title(69,ctype)
      read(15,'(a80)') line

      if(iqchem.eq.0) then

!! standard/MOKIT layout: Alpha Orbital Energies precedes Alpha MO      !!
!! coefficients in the source -- insert the synthetic Beta block right  !!
!! after each real Alpha one.                                           !!
        do while(index(line,"Alpha Orbital").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do
        write(69,11) "Alpha Orbital Energies","R","N= ",indepigr
        write(69,13) (energ_a(ii),ii=1,indepigr)
        write(69,11) "Beta Orbital Energies ","R","N= ",indepigr
        write(69,13) (energ_b(ii),ii=1,indepigr)

        call skip_fchk_orbblock(15,indepigr,line)
        write(69,12) "Alpha MO coefficients","R","N= ",norb
        write(69,13) ((pcoef_a(ii,jj),ii=1,igr),jj=1,indepigr)
        write(69,12) "Beta MO coefficients ","R","N= ",norb
        write(69,13) ((pcoef_b(ii,jj),ii=1,igr),jj=1,indepigr)

        call skip_fchk_orbblock(15,norb,line)
      else

!! Q-Chem layout: Alpha MO coefficients precedes Alpha Orbital Energies !!
!! in the source -- unverified, no restricted+Q-Chem+post-HF source     !!
!! available to test against, but structurally symmetric with the       !!
!! standard/MOKIT branch above.                                         !!
        do while(index(line,"Alpha MO co").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do
        write(69,12) "Alpha MO coefficients","R","N= ",norb
        write(69,13) ((pcoef_a(ii,jj),ii=1,igr),jj=1,indepigr)
        write(69,12) "Beta MO coefficients ","R","N= ",norb
        write(69,13) ((pcoef_b(ii,jj),ii=1,igr),jj=1,indepigr)

        call skip_fchk_orbblock(15,norb,line)
        write(69,11) "Alpha Orbital Energies","R","N= ",indepigr
        write(69,13) (energ_a(ii),ii=1,indepigr)
        write(69,11) "Beta Orbital Energies ","R","N= ",indepigr
        write(69,13) (energ_b(ii),ii=1,indepigr)

        call skip_fchk_orbblock(15,indepigr,line)
      end if

!! copy everything else through unchanged, Total SCF Density included    !!
      do while(.true.)
        write(69,'(a80)') line
        read(15,'(a80)',end=99) line
      end do
99    continue
      close(69)
      DEALLOCATE(energ_a,energ_b)

!! Printing formats !!
11    FORMAT(a22,21x,a1,3x,a3,i11)
12    FORMAT(a21,22x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))

      end