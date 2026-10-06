!! ********************************************************************* !!
!! subroutine: read_input                                                !!
!! purpose: parses every '# METHOD'/section keyword from the .inp file   !!
!!   (unit 16, opened by the caller) into the flags in input_options_mod,!!
!!   plus a handful of already-COMMON quantities (icas, ibcp, aerf,      !!
!!   nrad22...). Also triggers the .fchk read (call input()) partway     !!
!!   through, once DENS is known -- same ordering as before extraction.  !!
!! arguments: none (unit 16/15 already open, all output via              !!
!!   input_options_mod + existing COMMON blocks)                         !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine read_input()

      use input_options_mod
      use integration_grid, only: check_grid

      implicit real*8(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /atlist/iatlist(maxat),icuat
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /exchg/exch(maxat,maxat),xmix
      common /achi/achi(maxat,maxat),ibcp
      common /erf/aerf,ierf
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22
      common /modgrid2/thr3
      common /dm1opt/densthresh_dm1
      common /edaiqa/xen,xcoul,xnn
      common /edaiqa2/i2deda,iipoints,xptxyz(2,3)
      common /twoel/twoeltoler
      common /efield/field(4),edipole
      common /printout/iaccur
      common /iops/iopt(200)

      dimension navect(maxat),missat(maxat)
      dimension iatpairs(2,maxat)
      character*80 namedm,linedm
      character*80 namefchk1,namefchk2

!! every input needs # METHOD !!
      call locate_block(16,"# METHOD",ii)
      if(ii.eq.0) then
        write(*,'(2x,a)') 'The input has no # METHOD block'
        call apost_stop(' # METHOD block not found')
      end if

!! keywords and blocks renamed in version 5 (no blanks or hyphens      !!
!! inside a name): the old spellings stop instead of being ignored      !!
      call renamed_block("# DFT-DM1","# DFTDM1")
      call renamed_block("# DFT-DM1 FUNCTIONAL","# DFTDM1_FUNCTIONAL")
      call renamed_block("# ATOM PAIRS DEFINITION",
     +  "# ATOM_PAIRS_DEFINITION")
      call renamed_block("# 2D PLOTS","# 2D_PLOTS")
      call renamed_keyword("# METHOD","DFT-DM1","DFTDM1")
      call renamed_keyword("# OSLO","FOLI TOLERANCE","FOLI_TOLERANCE")
      call renamed_keyword("# OSLO","BRANCH ITERATION","BRANCH_ITERATION")
      call renamed_keyword("# OSLO","PRINT NON-ORTHO","PRINT_NONORTHO")
      call renamed_keyword("# EDAIQA","eN pySCF","eN_pySCF")
      call renamed_keyword("# EDAIQA","Coul pySCF","Coul_pySCF")
      call renamed_keyword("# EDAIQA","NN pySCF","NN_pySCF")

!! choose density from fchk file !!
      call readint("# METHOD","DENS",ndens0,1,1)
      iopt(9) = ndens0

!! processing fchk file !!
      call input()

      idoint=0
      iwfn=0

!! look for options !!
      call readchar("# METHOD","WFN",iwfn)
      call readchar("# METHOD","ALLPOINTS",iallpo)
      call readchar("# METHOD","FULLPRECISION",iaccur)

!! atoms in molecules !!
      call readchar("# METHOD","MULLIKEN",imulli)
      call readchar("# METHOD","MULLI",ii)
      if(ii.eq.1) imulli=1
      call readchar("# METHOD","LOWDIN",ilow)
      if(ilow.eq.1) imulli=2
      call readchar("# METHOD","LOWDIN-DAVIDSON",ilow)
      if(ilow.eq.1) imulli=3
      call readchar("# METHOD","NAO-BASIS",ilow)
      if(ilow.eq.1) imulli=4
      call readchar("# METHOD","LOWDIN-W",ilow)
      if(ilow.eq.1) imulli=5
      call readchar("# METHOD","HIRSH",ihirsh)
      call readchar("# METHOD","HIRSH-IT",ihirsh0)
      if(ihirsh0.eq.1) ihirsh=2
      call readchar("# METHOD","BECKE-RHO",ibcp)
      call readchar("# METHOD","NEWBEC",inewbec)
      call readchar("# METHOD","TFVC",itfvc)
      if(itfvc.eq.1) then
        ibcp=1
        inewbec=1
        istiff=4
      end if
      call readint("# METHOD","STIFFNESS",istiff,4,1)
      call readchar("# METHOD","WMATRAD",iradmat)
      call readchar("# METHOD","RMATRAD",iradmat0)
      if(iradmat0.eq.1) iradmat=2

      call readchar("# METHOD","ERF_PROF",ierf)
      call readreal("# METHOD","ERF_PROF",aerf,6.266d0,1)
      if(ierf.eq.1)  istiff=0

!! ERC: QTAIM input module !!
      call readchar("# METHOD","QTAIM",iqtaim)
      call readchar("# METHOD","READINT",ireadint)
      if(ireadint.eq.1) iqtaim=2
!! plain Becke atoms (fixed empirical radii) only on request; "BECKE"   !!
!! is also inside "BECKE-RHO", which selects its own scheme             !!
      call readchar("# METHOD","BECKE",ibecke)
      if(ibcp.eq.1) ibecke=0

!! no atomic definition at all: the TFVC flags, with a warning (itfvcdef) !!
      itfvcdef=0
      if(imulli.eq.0.and.ihirsh.eq.0.and.ibcp.eq.0.and.inewbec.eq.0
     +   .and.itfvc.eq.0.and.iqtaim.eq.0.and.ibecke.eq.0) then
        itfvc=1
        ibcp=1
        inewbec=1
        itfvcdef=1
        write(*,'(2x,a)') 'WARNING: no atomic definition in # METHOD, '//
     +    'TFVC used; add TFVC or the scheme you want'
      end if

      if(iqtaim.eq.1) then
        call readint("# QTAIM","STEP",istep,300,1)
        call readint("# QTAIM","NNA",inna,0,1)
        call readint("# QTAIM","MAXDIST",imaxdist,12,1)
        call readint("# QTAIM","SCREENING",iscreening,1000,1)
        call readint("# QTAIM","PATH",ipath,0,1)
      end if 

!! miscellaneous options !!
      call readchar("# METHOD","OPOP",iopop)
      call readchar("# METHOD","SHANNON",isha)
      call readchar("# METHOD","DOINT",idoint)
      call readchar("# METHOD","PCA",ipca)
      call readchar("# METHOD","LAPLACIAN",ilaplacian)
      call readchar("# METHOD","FINEGRID",ifinegrid)
      call readchar("# METHOD","ELCOUNT",ielcount) !MMO- NCTAIM
      call readint("# METHOD","RHO_CALC_AT",iatdens,0,1)
      call readreal("# METHOD","RHO_CALC_RAD",Rmax,0.0d0,1) ! fixed: was 0 (integer), must be REAL*8
      call readchar("# METHOD","NOPOPU",inopop)

!! eff-AO-s and EOS !!
      call readchar("# METHOD","EFFAO",ieffao)
      call readchar("# METHOD","UEFFAO",idummy)
      if(idummy.eq.1) ieffao=2

!! effAOs, paired and unpaired only !!
      call readchar("# METHOD","EFFAO-U",idummy)
      if(idummy.eq.1) ieffao=3

      call readint("# METHOD","EFF_THRESH",ieffthr,1,1)
      call readchar("# METHOD","CUBE",icube)
      inegefos=0
      inegcubthr=25
      if(icube.eq.1) then
        call locate_block(16,"# CUBE",ii)
        if(ii.eq.0) call apost_stop('Required section # CUBE not found in input file')
        call readint("# CUBE","MAX_OCC",jcubthr,1000,1)
        call readint("# CUBE","MIN_OCC",kcubthr,0,1)
        call readreal("# CUBE","SPACING",cubespacing,0.25d0,1)
        call readreal("# CUBE","RADIUS_SCALE",cuberadscale,2.0d0,1)
        call readchar("# CUBE","NEG_EFOS",inegefos)
        if(inegefos.eq.1)
     +    call readint("# CUBE","NEG_EFOS",inegcubthr,25,1)
!! occupations x1000, never negative -- negative-occupation EFOs have    !!
!! their own NEG_EFOS keyword                                            !!
        if(jcubthr.lt.0.or.kcubthr.lt.0) then
          write(*,'(2x,a)') 'MAX_OCC/MIN_OCC must be >= 0 (occupation x 1000).'
          write(*,'(2x,a)') 'For negative-occupation EFOs use NEG_EFOS instead.'
          call apost_stop('')
        end if
        if(kcubthr.gt.jcubthr) call apost_stop('MIN_OCC cannot be larger than MAX_OCC')
        if(inegcubthr.lt.0) call apost_stop('NEG_EFOS value must be >= 0 (|occupation| x 1000)')
      end if

!! EOS (standard) !!
      call readchar("# METHOD","EOS",ieos)
      if(ieos.eq.1) then 
        iopop=1
        ieffao=2
        call readreal("# METHOD","EOS_THRESH",xthresh,2.5d-3,1)
      end if

!! GEOS: EOS from the paired and unpaired densities !!
      call readchar("# METHOD","GEOS",iueos)
      if(iueos.eq.1) then
        iopop=1
        ieffao=3
!! readchar matches by substring, so the "EOS" scan just above also        !!
!! matched this GEOS line (pre-existing issue, same as EFFAO/EFFAO-U/UEFFAO !!
!! below -- harmless there since ieffao gets overwritten either way, but   !!
!! ieos/iueos are separate flags with their own downstream consumers, so   !!
!! undo that false positive here.                                          !!
        ieos=0
        call readreal("# METHOD","EOS_THRESH",xthresh,2.5d-3,1)
      end if

!! OS from centroids !!
      call readchar("# METHOD","OS-CENTROID",ieoscent)

!! OS from localized orbitals (LOBA) !!
      call readchar("# METHOD","LOBA",iloba)

!! local spin and methods for correlated WFs !!
      call readchar("# METHOD","SPIN",ispin)
      call readint("# METHOD","DM",icorr,0,1)
      if(icorr.eq.2) ispin=1
      call readchar("# METHOD","DAFH",idafh)

!! energy decomposition options !!
      call readchar("# METHOD","ENPART",ienpart )
      if(ienpart.eq.1) then
        call locate_block(16,"# ENPART",ii)
        if(ii.eq.0) then
          write(*,'(2x,a)') 'ENPART needs a # ENPART block'
          call apost_stop(' # ENPART block not found')
        end if
        xmix=ZERO
        id_xcfunc=0
        id_xfunc=0
        id_cfunc=0
        call readchar("# ENPART","LIBRARY",ilib)
        call enpart_functional_keyword(ikwfunc,idkxc,idkx,idkc)
        if(ilib.eq.1) then
          if(ikwfunc.eq.1)
     +      call apost_stop('GIVE EITHER LIBRARY OR A FUNCTIONAL KEYWORD IN # ENPART. REVISE inp')
          call readint("# ENPART","EXC_FUNCTIONAL",id_xcfunc,0,1)
          call readint("# ENPART","EX_FUNCTIONAL",id_xfunc,0,1)
          call readint("# ENPART","EC_FUNCTIONAL",id_cfunc,0,1)
          id_func=id_xfunc+id_cfunc+id_xcfunc
          if (id_func.eq.0) call apost_stop('FUNCTIONAL ID NOT FOUND IN INPUT FILE')
!! xc (enpart_dft.f) and func_info_print's xmix assume one or the other. !!
          if(id_xcfunc.ne.0.and.(id_xfunc.ne.0.or.id_cfunc.ne.0)) then
            write(*,'(2x,a)') 'Give either EXC_FUNCTIONAL alone, or EX_FUNCTIONAL and/or'
            write(*,'(2x,a)') 'EC_FUNCTIONAL -- not both kinds together.'
            call apost_stop('EXC_FUNCTIONAL COMBINED WITH EX/EC_FUNCTIONAL. REVISE inp')
          end if
!! specific keywords for functionals !!
        else
          call readchar("# ENPART","HF",ihf)
          if(ihf.eq.1) then
            id_xfunc=-1
            xmix=1.0d0
            goto 233
          end if
!! predefined functionals (SVWN, BLYP, B3LYP, PBE0, ...), table in      !!
!! enpart_functional_keyword                                            !!
          if(ikwfunc.eq.1) then
            id_xcfunc=idkxc
            id_xfunc=idkx
            id_cfunc=idkc
            go to 233
          end if

!! specific keywords for correlated methods -- use CORRELATION to        !!
!! decompose both X and C, default is decompose XC                       !!
          call readchar("# ENPART","CASSCF",icas)
          call readchar("# ENPART","CISD",icisd)
          call readchar("# ENPART","CORRELATION",iecorr)
          if(icas.eq.0. and.icisd.eq.0) then
            call apost_stop("NO DFT/HF/CASSCF/CISD SELECTED FOR ENPART. REVISE inp")
          end if
233       continue
        end if 

!! extra options !!
        call readint("# ENPART","THREBOD",ithrebod,50,1) !! bond order 0.005; below 1 sets it to zero !!
        call readchar("# ENPART","EXACT",iexact)
        call readchar("# ENPART","HOMO",ihomo)
        call readchar("# ENPART","DEKIN",idek)
        call readchar("# ENPART","IONIC",iionic)
        call readreal("# ENPART","TWOELTOLER",twoeltoler,0.25d0,1)
        call readchar("# ENPART","ANALYTIC",ianalytical)

!! adding grid tuning for two-el integration !!
        call read_gridtwoel("# ENPART",ienpart_gridtwoel)

!! for topology calculation !!
!! MG: needs to be properly checked, done a long time ago                !!
!! MG: an extended version for 2D/3D and more 1D topology exists in       !!
!! apost3.1-devel -- worth checking if merging it in is worthwhile        !!
        itop=0
        ipairs=0
        call readchar("# METHOD","TOPOLOGY",itop)
        if(itop.eq.1) then
          call locate_block(16,"# ATOM_PAIRS_DEFINITION",ii)
          if(ii.eq.0) call apost_stop(" # ATOM_PAIRS_DEFINITION section missing")
          read(16,*) ipairs
          if(ipairs.gt.0) then
            do ii=1,ipairs
              read(16,*) (iatpairs(kk,ii),kk=1,2)
              write(*,*) (iatpairs(kk,ii),kk=1,2)
            end do
          else
            write(*,*) " DOING CUBE OF THE ENTIRE MOLECULAR SYSTEM "
          end if
!! for choosing the energy component to do the topology on !!
          ietop=-1
          call readchar("# TOPOLOGY","EXCHANGE",itop2)
          if(itop2.eq.1) ietop=1
          call readchar("# TOPOLOGY","CORRELATION",itop2)
          if(itop2.eq.1) ietop=2
          call readchar("# TOPOLOGY","EXCHANGE-CORRELATION",itop2)
          if(itop2.eq.1) ietop=3
          call readchar("# TOPOLOGY","DENSITY",itop2)
          if(itop2.eq.1) ietop=9
          if(ietop.eq.-1) call apost_stop(" FUNCTION FOR TOPOLOGY NOT INTRODUCED ")
        end if

!! end of ENPART options !!
      end if

!! DFT-DM1 approximate one-particle RDM1 for UHF/UKS-DFT (formerly       !!
!! referred to internally as HIRAO). Runs standalone (no ENPART          !!
!! required), reusing ENPART's own two-electron grid/defaults if ENPART  !!
!! is also active, else reading its own MOD-GRIDTWOEL under # DFTDM1.    !!
      idftdm1=0
      id_func_dm1=0
      call readchar("# METHOD","DFTDM1",idftdm1)
      if(idftdm1.eq.1) then
!! DFT-DM1 is DFT-functional-only -- dft_dm1.f's RDM1 construction reads !!
!! the local exchange-energy density straight out of libxc, which has no !!
!! meaning for HF (no functional to evaluate); an old ifunc=999 HF        !!
!! placeholder branch existed in dft_dm1.f but was never a real          !!
!! implementation, and was removed 2026-08-30.                           !!
        call readchar("# DFTDM1_FUNCTIONAL","LIBRARY",ilib)
        if(ilib.eq.1) then
          call readint("# DFTDM1_FUNCTIONAL","EX_FUNCTIONAL",id_func_dm1,0,1)
        end if
        if(id_func_dm1.eq.0) call apost_stop("FUNCTIONAL ID NOT FOUND FOR DFT-DM1. REVISE inp")

!! density-threshold pruning for the double loop's O(itotps^2) grid-point !!
!! pairs, checked at the pair's midpoint R (dft_dm1.f's main loop)        !!
        call readreal("# DFTDM1","DENSTHRESH",densthresh_dm1,1.0d-8,1)

!! project the RDM1 onto the AO basis and diagonalize for natural-orbital !!
!! occupations -- opt-in, since this adds an O(igr^2) accumulation on     !!
!! top of every surviving grid-point pair (dft_dm1.f's own separate,      !!
!! unoptimized block, run only if requested)                              !!
        call readchar("# DFTDM1","NATORB",inatorb_dm1)

!! idftdm1grid is a throwaway local -- nothing outside this call needs   !!
!! DFT-DM1's own MOD-GRIDTWOEL flag today, unlike ENPART's (see above).  !!
        if(ienpart.ne.1) call read_gridtwoel("# DFTDM1",idftdm1grid)
      end if

!! EDAIQA options !!
      iedaiqa=0
      iflip=0
      call readchar("# METHOD","EDAIQA",iedaiqa)
      if(iedaiqa.eq.1) then
        ii=0
        call locate_block(16,"# EDAIQA",ii)
        if(ii.eq.1) then
          read(16,'(a80)') namefchk1
          read(16,'(a80)') namefchk2
          open(unit=55,file=namefchk1)
          open(unit=52,file=namefchk2)
          call readchar("# EDAIQA","FLIPSPIN",iflip) !! swaps alpha for beta !!

!! adding pySCF reference values for the electrostatic calculation !!
          call readreal("# EDAIQA","eN_pySCF",xen,0.0d0,1)
          call readreal("# EDAIQA","Coul_pySCF",xcoul,0.0d0,1)
          call readreal("# EDAIQA","NN_pySCF",xnn,0.0d0,1)

!! adding grid tuning for two-el integration -- MG: repetitive with the  !!
!! # ENPART block above, could be consolidated into one                 !!
          call readchar("# EDAIQA","MOD-GRIDTWOEL",iigrid)
          call readint("# GRID","RADIAL",nrad22,40,1)
          call readint("# GRID","ANGULAR",nang22,146,1)
          call readreal("# GRID","rr00",rr0022,0.5d0,1)
          call readreal("# GRID","phb1",phb12,0.162d0,1)
          call readreal("# GRID","phb2",phb22,0.182d0,1)
          call check_grid(nrad22,nang22,"# GRID")

!! options to make 2D plots of electrostatic potentials !!
          i2deda=0
          call locate_block(16,"# 2D_PLOTS",i2deda)
          if(i2deda.eq.1) then
            read(16,*) iipoints
            read(16,*) (xptxyz(1,j),j=1,3)
            read(16,*) (xptxyz(2,j),j=1,3)

!! printing to confirm the values read in -- candidate for removal !!
            write(*,*) " "
            write(*,*) " SOME PRINTING FOR EDAIQA PURPOSES "
            write(*,*) " "
            write(*,*) " Number of points for 2D plots ",iipoints
            write(*,*) " Coordinates of the two ghost atoms "
            write(*,*) (xptxyz(1,j),j=1,3)
            write(*,*) (xptxyz(2,j),j=1,3)
            write(*,*) " "
          end if
        else
          call apost_stop("EDAIQA SECTION MISSING. REVISE inp")
        end if

!! end of EDAIQA !!
      end if

!! nonlinear optical properties (POLAR) !!
      call readchar("# METHOD","POLAR",ipolar )
      if(ipolar.eq.1) then
        iaccur=1
!! MMO: the old $name.scr raw-output-for-post-processing path was       !!
!! removed, no longer used here                                         !!
      end if

      call field_misc(ifield)
      if(ifield.eq.1) then
        write(*,*) 'The system is under a static electric field' 
        write(*,'(3(a4,f8.6))') 'Fx=',field(2), 'Fy=',field(3),'Fz=',field(4) 
      end if

!! do for restricted number of atoms !!
      idoat=0
      call readchar("# METHOD","DOATOMS",idoat)
      if(idoat.eq.1) then
        call locate_block(16,"# ATOMS",ii)
        if(ii.eq.0) call apost_stop('Required section not found in input file')
        read(16,*) icuat
        read(16,*) (iatlist(i),i=1,icuat)
      else
        icuat=nat
        do i=1,icuat
          iatlist(i)=i
        end do
      end if

!! do for fragments !!
      idofr=0
      call readchar("# METHOD","DOFRAGS",idofr)
      if(idofr.eq.1) then
        call locate_block(16,"# FRAGMENTS",ii)
        if(ii.eq.0) call apost_stop('Required section not found in input file')
        read(16,*) icufr
        if(icufr.lt.1.or.icufr.gt.nat) then
          write(*,'(2x,a,i0,a,i0,a)') '# FRAGMENTS: ',icufr,
     +      ' fragments given, the molecule has ',nat,' atoms'
          call apost_stop(' Wrong number of fragments in # FRAGMENTS')
        end if
        do i=1,icufr
          read(16,*) nfrlist(i)
          if(nfrlist(i).eq.-1.and.i.ne.icufr) then
            write(*,'(2x,a,i0,a)') '# FRAGMENTS: -1 (all remaining '//
     +        'atoms) given for fragment ',i,
     +        ', only allowed for the last one'
            call apost_stop(' -1 only allowed for the last fragment in # FRAGMENTS')
          end if
          if(nfrlist(i).eq.-1) then
            do l=1,nat
              navect(l)=0
            end do
            do l=1,(icufr-1)
              do k=1,nfrlist(l)
                navect(ifrlist(k,l))=1
              end do
            end do
            k=0
            do l=1,nat
              if(navect(l).eq.0) then
                k=k+1
                ifrlist(k,icufr)=l
              end if
            end do
            nfrlist(icufr)=k
            if(k.eq.0) then
              write(*,'(2x,a)') '# FRAGMENTS: no atoms left for the '//
     +          'last fragment (-1)'
              call apost_stop(' Empty fragment in # FRAGMENTS')
            end if
          else
            if(nfrlist(i).lt.1.or.nfrlist(i).gt.nat) then
              write(*,'(2x,a,i0,a,i0,a)') '# FRAGMENTS: fragment ',i,
     +          ' has ',nfrlist(i),' atoms'
              call apost_stop(' Wrong number of atoms in # FRAGMENTS')
            end if
            read(16,*) (ifrlist(k,i),k=1,nfrlist(i))
            do k=1,nfrlist(i)
              if(ifrlist(k,i).lt.1.or.ifrlist(k,i).gt.nat) then
                write(*,'(2x,a,i0,a,i0,a,i0,a)') '# FRAGMENTS: atom ',
     +            ifrlist(k,i),' (fragment ',i,') does not exist, '//
     +            'the molecule has ',nat,' atoms'
                call apost_stop(' Atom out of range in # FRAGMENTS')
              end if
            end do
          end if
        end do

!! jfrlist tells which fragment a given atom belongs to (0: none) !!
        do i=1,nat
          jfrlist(i)=0
        end do
        do i=1,icufr
          do k=1,nfrlist(i)
            if(jfrlist(ifrlist(k,i)).eq.i) then
              write(*,'(2x,a,i0,a,i0)') '# FRAGMENTS: atom ',
     +          ifrlist(k,i),' is listed twice in fragment ',i
              call apost_stop(' Atom listed twice in # FRAGMENTS')
            else if(jfrlist(ifrlist(k,i)).ne.0) then
              write(*,'(2x,a,i0,a,i0,a,i0)') '# FRAGMENTS: atom ',
     +          ifrlist(k,i),' is in fragments ',jfrlist(ifrlist(k,i)),
     +          ' and ',i
              call apost_stop(' Atom in two fragments in # FRAGMENTS')
            end if
            jfrlist(ifrlist(k,i))=i
          end do
        end do
!! for compatibility !!
      else
        icufr=nat
        do i=1,icufr
          nfrlist(i)=1
          ifrlist(1,i)=i                
          jfrlist(i)=i                
        end do
      end if

!! extra warnings for EDAIQA !!
      if(iedaiqa.eq.1) then
        if(idofr.eq.0) then
          write(*,*) " FRAGMENT DEFINITION REQUIRED FOR EDAIQA "
          write(*,*) " INFO : FRAGMENT ORDER MUST MATCH THE ISOLATED "
          call apost_stop('')
        end if
        if(idofr.eq.1.and.icufr.ne.2) call apost_stop(" ONLY 2 FRAGMENTS ALLOWED FOR EDAIQA ")
      end if

!! reading DM1 and DM2 !!
      if(icorr.ne.0) then
        call locate_block(16,"# DM",ii)
        if(ii.eq.0) call apost_stop(" # DM section not found in input file ")

!! the file names come first (one line each, read as they are, so a     !!
!! name may hold any character); the format keywords only after them    !!
        read(16,'(a80)') namedm
        dmfile1=adjustl(namedm)
        dmfile2=' '
        if(icorr.eq.2) then
          read(16,'(a80)') namedm
          dmfile2=adjustl(namedm)
        end if
        ipyscf=0
        iorca=0
        idmrg=0
        ii=0
        do while(ii.eq.0)
          read(16,'(a80)',end=520) linedm
          call keyword_value(linedm,"pySCF",ipos)
          if(ipos.gt.0) ipyscf=1
          call keyword_value(linedm,"ORCA",ipos)
          if(ipos.gt.0) iorca=1
          call keyword_value(linedm,"DMRG",ipos)
          if(ipos.gt.0) idmrg=1
          if(index(linedm,"#").ne.0) ii=1
        end do
 520    continue

        if(iorca.eq.1.or.ipyscf.eq.1) then
          open(11,file=dmfile1,status='OLD',iostat=ios)
        else
          open(11,file=dmfile1,FORM='UNFORMATTED',status='OLD',iostat=ios)
        end if
        if(ios.ne.0) then
          write(*,'(2x,a,a,a)') '# DM: file ',trim(dmfile1),' not found'
          call apost_stop(' # DM file not found')
        end if
        if(icorr.eq.2) then
          if(iorca.eq.1.or.ipyscf.eq.1) then
            open(12,file=dmfile2,status='OLD',iostat=ios)
          else
            open(12,file=dmfile2,FORM='UNFORMATTED',status='OLD',
     +        iostat=ios)
          end if
          if(ios.ne.0) then
            write(*,'(2x,a,a,a)') '# DM: file ',trim(dmfile2),' not found'
            call apost_stop(' # DM file not found')
          end if
        end if
      end if

!! OSLO options !!

      ioslo=0
      call readchar("# METHOD","OSLO",ioslo)
      if(ioslo.eq.1) then
        call locate_block(16,"# OSLO",ii)
        if(ii.eq.0) then
          write(*,'(2x,a)') 'OSLO needs a # OSLO block (it may be empty)'
          call apost_stop(' # OSLO block not found')
        end if

!! MG: by default requires the TFVC AIM in # METHOD (numerical           !!
!! integration), but one can ask for OSLOs using Hilbert-space AIMs      !!
!! instead -- the Hilbert-space cases follow                             !!

        ilow2=0
        call readchar("# OSLO","MULLIKEN",ii)
        if(ii.eq.1) ilow2=1
        call readchar("# OSLO","LOWDIN",ii)
        if(ii.eq.1) ilow2=2
        call readchar("# OSLO","LOWDIN-DAVIDSON",ii)
        if(ii.eq.1) ilow2=3
        call readchar("# OSLO","NAO-BASIS",ii)
        if(ii.eq.1) ilow2=6

!! extra options !!

        call readint("# OSLO","FOLI_TOLERANCE",ifolitol,3,1) !! FOLI value tolerance, for selection !!
        call readint("# OSLO","BRANCH_ITERATION",ibranch,0,1) !! iteration to invoke branching at !!
        ioslofchk=1
        call readchar("# OSLO","PRINT_NONORTHO",ii) !! prints non-ortho OSLOs to an extra .fchk file !!
        if(ii.eq.1) ioslofchk=2

      end if

!! to print .fchk files from QCHEM or MOKIT -- MG: their fchk format is  !!
!! different from Gaussian's                                             !!
      iqchem=0
      call readchar("# METHOD","QCHEM",iqchem)
      imokit=0
      call readchar("# METHOD","MOKIT",imokit)

!! end of OSLO options !!

!! fragments: required by the oxidation-state methods (EOS, GEOS, OSLO), !!
!! which together with ENPART also need every atom in a fragment; other  !!
!! analyses may leave atoms out on purpose, so only warn there           !!
      if((ieos.eq.1.or.iueos.eq.1.or.ioslo.eq.1).and.idofr.eq.0) then
        write(*,'(2x,a)') 'EOS, GEOS and OSLO need fragments: add '//
     +    'DOFRAGS to # METHOD and a # FRAGMENTS block'
        write(*,'(2x,a)') '(an atom can be a fragment on its own)'
        call apost_stop(' EOS, GEOS and OSLO need fragments (DOFRAGS)')
      end if
      if(idofr.eq.1) then
        nmiss=0
        do i=1,nat
          if(jfrlist(i).eq.0) then
            nmiss=nmiss+1
            missat(nmiss)=i
          end if
        end do
        if(nmiss.gt.0) then
          if(ieos.eq.1.or.iueos.eq.1.or.ioslo.eq.1.or.ienpart.eq.1) then
            write(*,'(2x,a,*(1x,i0))') '# FRAGMENTS: atoms in no '//
     +        'fragment:',(missat(i),i=1,nmiss)
            write(*,'(2x,a)') 'EOS, GEOS, OSLO and ENPART need every '//
     +        'atom in a fragment (-1 on the last one adds the rest)'
            call apost_stop(' Atoms missing in # FRAGMENTS')
          else
            write(*,'(2x,a,*(1x,i0))') 'WARNING: # FRAGMENTS leaves '//
     +        'out atoms',(missat(i),i=1,nmiss)
          end if
        end if
      end if

!! X-ray scattering factors !!
      call readchar("# METHOD","SCATT-FACT",iscattfact)

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: read_gridtwoel                                            !!
!! purpose: reads the MOD-GRIDTWOEL/# GRID override shared by ENPART's   !!
!!   own two-electron integration grid and DFT-DM1's grid (both default  !!
!!   to 150/590 either way) -- extracted out of read_input()'s ENPART    !!
!!   block so DFT-DM1 can reuse the exact same mechanism/defaults        !!
!!   without requiring ENPART to also be active in the same run.         !!
!! arguments:                                                            !!
!!   section     (in)  -- .inp section to scan MOD-GRIDTWOEL/# GRID      !!
!!     under (e.g. "# ENPART" or "# DFTDM1")                            !!
!!   igridtwoel  (out) -- 1 if MOD-GRIDTWOEL was set for this section,   !!
!!     0 otherwise. Callers that need to know whether the returned       !!
!!     nrad22/nang22 came from the user's own # GRID (vs. the plain      !!
!!     defaults) must use this, not input_options_mod's iigrid -- that   !!
!!     one is EDAIQA's own separate flag, not shared with this call.     !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine read_gridtwoel(section,igridtwoel)
      use integration_grid, only: check_grid
      implicit real*8(a-h,o-z)
      character section*(*)
      character*80 linia
      integer, intent(out) :: igridtwoel
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22
      common /modgrid2/thr3

      call readchar(section,"MOD-GRIDTWOEL",igridtwoel)
      if(igridtwoel.eq.1) then
        call readint("# GRID","RADIAL",nrad22,150,1)
        call readint("# GRID","ANGULAR",nang22,590,1)
        call readreal("# GRID","rr00",rr0022,0.5d0,1)
        call readreal("# GRID","phb1",phb12,0.169d0,1)
        call readreal("# GRID","phb2",phb22,0.170d0,1)
        call readreal("# GRID","THRESH2",thr3,1.0d-12,1)
        call check_grid(nrad22,nang22,"# GRID")

!! defaults, modified for safe integration setup !!
      else
        nrad22=150
        nang22=590
        rr0022=0.5
        phb12=0.169d0
        phb22=0.170d0
        thr3=1.0d-12

!! Warn if a # GRID block exists in the input but is being ignored        !!
!! because MOD-GRIDTWOEL wasn't set. Scan the file directly here rather   !!
!! than via "locate", which always prints a "section not found" message  !!
!! on a miss                                                              !!
        igridpresent=0
        rewind(16)
        iiscan=0
        do while(iiscan.eq.0)
          read(16,"(a80)",end=234) linia
          if(index(linia,"# GRID").ne.0) then
            igridpresent=1
            iiscan=1
          end if
        end do
234     continue
        if(igridpresent.eq.1) then
          write(*,*) " "
          write(*,*) "WARNING: a # GRID block was found in the input, but MOD-GRIDTWOEL"
          write(*,*) "was not set in ",trim(section),". Using the default integration setup instead"
          write(*,*) "Add MOD-GRIDTWOEL to ",trim(section)," to apply your # GRID settings"
          write(*,*) " "
        end if
      end if

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: enpart_functional_keyword                                 !!
!! purpose: looks up a predefined functional keyword in the # ENPART     !!
!!   section and returns its libxc ids. Matches the whole first word of  !!
!!   each line, case-insensitive (readchar's substring match would mix   !!
!!   up SVWN/SVWN5, PBE/PBE0, ...). Each mapping was checked against     !!
!!   Gaussian 16 wavefunctions with ENPART's two-electron error.         !!
!! arguments:                                                            !!
!!   kfound (out) -- 1 if a functional keyword was found, 0 otherwise    !!
!!   idxc, idx, idc (out) -- libxc ids: combined xc, exchange,           !!
!!     correlation (0 where not used)                                    !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine enpart_functional_keyword(kfound,idxc,idx,idc)
      implicit none
      integer, intent(out) :: kfound,idxc,idx,idc
      integer, parameter :: nkw=13
      character*10 kwname(nkw)
      integer kwxc(nkw),kwx(nkw),kwc(nkw)
      character*80 linea
      character*10 tok
      integer ii,k,kk,ic,iend,nfound

!! idx=-2 marks a keyword that is refused with a message (see below).    !!
!! B3P86 is 315, not 403: 403 is ~100 kcal/mol off Gaussian's B3P86.     !!
      data kwname /'SVWN','SVWN5','BLYP','BP86','PBE','PBEPBE','B3LYP',
     +  'B3PW91','B3P86','PBE0','PBE1PBE','BHANDHLYP','LDA'/
      data kwxc /0,0,0,0,0,0,402,401,315,406,406,436,0/
      data kwx  /1,1,106,106,101,101,0,0,0,0,0,0,-2/
      data kwc  /8,7,131,132,130,130,0,0,0,0,0,0,0/

      kfound=0
      idxc=0
      idx=0
      idc=0
      nfound=0
      call locate_block(16,"# ENPART",ii)
      if(ii.eq.0) return
      do
        read(16,'(a80)',end=10) linea
        if(index(linea,"#").ne.0) exit
        linea=adjustl(linea)
        iend=scan(linea,' =')
        if(iend.le.1) cycle
        tok=linea(1:min(iend-1,10))
        do kk=1,len_trim(tok)
          ic=ichar(tok(kk:kk))
          if(ic.ge.ichar('a').and.ic.le.ichar('z')) tok(kk:kk)=char(ic-32)
        end do
        do k=1,nkw
          if(tok.eq.kwname(k)) then
            nfound=nfound+1
            if(kwx(k).eq.-2) then
              write(*,'(2x,a)') 'The LDA keyword is no longer accepted (it meant Slater'
              write(*,'(2x,a)') 'exchange only). Use SVWN or SVWN5, or LIBRARY with'
              write(*,'(2x,a)') 'EX_FUNCTIONAL 1 for exchange only.'
              call apost_stop('LDA KEYWORD REMOVED. REVISE inp')
            end if
            kfound=1
            idxc=kwxc(k)
            idx=kwx(k)
            idc=kwc(k)
          end if
        end do
      end do
10    continue
      if(nfound.gt.1) call apost_stop('MORE THAN ONE FUNCTIONAL KEYWORD IN # ENPART. REVISE inp')

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: renamed_keyword                                           !!
!! purpose: stops the run when an .inp block uses a keyword's old        !!
!!   spelling, naming the new one.                                       !!
!! arguments:                                                            !!
!!   section, old, new (in) -- block header, old and new keyword         !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine renamed_keyword(section,old,new)
      character section*(*), old*(*), new*(*)
      character linea*80
      integer ii,ipos

      call find_block(16,section,ii)
      if(ii.eq.0) return
      ii=0
      do while(ii.eq.0)
        read(16,"(a80)",end=10) linea
        call keyword_value(linea,old,ipos)
        if(ipos.gt.0) then
          write(*,'(2x,a,a,a,a,a,a)') trim(section),': ',old,
     +      ' is now written ',new
          call apost_stop(' Renamed keyword in the input')
        end if
        if(index(linea,"#").ne.0) ii=1
      end do
10    return
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: renamed_block                                             !!
!! purpose: stops the run when the .inp uses a block's old name, naming  !!
!!   the new one.                                                        !!
!! arguments:                                                            !!
!!   old, new (in) -- old and new block header                           !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine renamed_block(old,new)
      character old*(*), new*(*)
      integer ii

      call find_block(16,old,ii)
      if(ii.eq.1) then
        write(*,'(2x,a,a,a)') old,' is now written ',new
        call apost_stop(' Renamed block in the input')
      end if
      end
