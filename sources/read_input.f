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
      common /edaiqa/xen,xcoul,xnn
      common /edaiqa2/i2deda,iipoints,xptxyz(2,3)
      common /twoel/twoeltoler
      common /efield/field(4),edipole
      common /printout/iaccur
      common /iops/iopt(200)

      dimension navect(maxat)
      dimension iatpairs(2,maxat)
      character*80 linia,namedm
      character*80 namefchk1,namefchk2

c choose density from fchk file
      call readint("# METHOD","DENS",ndens0,1,1)
      iopt(9) = ndens0
      
C Processing FChk file
      call input()

      idoint=0
      iwfn=0  

c look for options      
      call readchar("# METHOD","WFN",iwfn)
      call readchar("# METHOD","ALLPOINTS",iallpo)
      call readchar("# METHOD","FULLPRECISION",iaccur)

C Atoms in molecules
      call readchar("# METHOD","MULLI",imulli)
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

c ERC QTAIM input module
      call readchar("# METHOD","QTAIM",iqtaim)
      call readchar("# METHOD","READINT",ireadint)
      if(ireadint.eq.1) iqtaim=2
      if(iqtaim.eq.1) then
        call readint("# QTAIM","STEP",istep,300,1)
        call readint("# QTAIM","NNA",inna,0,1)
        call readint("# QTAIM","MAXDIST",imaxdist,12,1)
        call readint("# QTAIM","SCREENING",iscreening,1000,1)
        call readint("# QTAIM","PATH",ipath,0,1)
      end if 

c Miscellaneous options
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
 
C eff-AO-s and EOS
      call readchar("# METHOD","EFFAO",ieffao)
      call readchar("# METHOD","UEFFAO",idummy)
      if(idummy.eq.1) ieffao=2

!! EFFAOS PAIRED AND UNPAIRED (ONLY) !!
      call readchar("# METHOD","EFFAO-U",idummy)
      if(idummy.eq.1) ieffao=3

      call readint("# METHOD","EFF_THRESH",ieffthr,1,1)
      call readchar("# METHOD","CUBE",icube)
      if(icube.eq.1) then
        call locate(16,"# CUBE",ii)
        if(ii.eq.0) stop'Required section # CUBE not found in input file'
        call readint("# CUBE","MAX_OCC",jcubthr,1000,1)
        call readint("# CUBE","MIN_OCC",kcubthr,0,1)
      end if

!! EOS (STANDARD) !!
      call readchar("# METHOD","EOS",ieos)
      if(ieos.eq.1) then 
        iopop=1
        ieffao=2
        call readreal("# METHOD","EOS_THRESH",xthresh,2.5d-3,1)
      end if

!! EOS FROM THE PAIRED AND UNPAIRED DENSITIES !!
      call readchar("# METHOD","EOS-U",iueos)
      if(iueos.eq.1) then 
        iopop=1
        ieffao=3
        call readreal("# METHOD","EOS_THRESH",xthresh,2.5d-3,1)
      end if

!! OS FROM CENTROIDS !!
      call readchar("# METHOD","OS-CENTROID",ieoscent)

!! OS FROM LOCALIZED ORBITALS (LOBA) !!
      call readchar("# METHOD","LOBA",iloba)

c Local spin and methods for correlated WFs
      call readchar("# METHOD","SPIN",ispin)
      call readint("# METHOD","DM",icorr,0,1)
      if(icorr.eq.2) ispin=1
      call readchar("# METHOD","DAFH",idafh)

c Energy decomposition options   
      call readchar("# METHOD","ENPART",ienpart )
      if(ienpart.eq.1) then
        xmix=ZERO
        id_xcfunc=0
        id_xfunc=0
        id_cfunc=0
        call readchar("# ENPART","LIBRARY ",ilib)
        if(ilib.eq.1) then
          call readint("# ENPART","EXC_FUNCTIONAL",id_xcfunc,0,1)
          call readint("# ENPART","EX_FUNCTIONAL",id_xfunc,0,1)
          call readint("# ENPART","EC_FUNCTIONAL",id_cfunc,0,1)
          id_func=id_xfunc+id_cfunc+id_xcfunc
          if (id_func.eq.0) stop 'FUNCTIONAL ID NOT FOUND IN INPUT FILE'
          go to 233
c specific keywords for functionanls 
        else 
          call readchar("# ENPART","HF ",ihf)
          if(ihf.eq.1) then
            id_xfunc=-1
            xmix=1.0d0
            goto 233
          end if
          call readchar("# ENPART","LDA",ival )
          if(ival.eq.1) then 
            id_xfunc=1
            go to 233
          end if
          call readchar("# ENPART","BP86",ival)
          if(ival.eq.1) then
            id_xfunc=106
            id_cfunc=132
            go to 233
          end if
          call readchar("# ENPART","B3LYP",ival)
          if(ival.eq.1) then
            id_xcfunc=402
            go to 233
          end if

c specific keywords for correlated methods
c use CORRELATION to decompose both X and C. Default is decompose XC.
          call readchar("# ENPART","CASSCF",icas)
          call readchar("# ENPART","CISD",icisd)
          call readchar("# ENPART","CORRELATION",iecorr)
          if(icas.eq.0. and.icisd.eq.0) then
            stop "NO DFT/HF/CASSCF/CISD SELECTED FOR ENPART. REVISE inp"
          end if
233       continue
        end if 

!!  EXTRA OPTIONS !! 
        call readint("# ENPART","THREBOD",ithrebod,100,1) !! SELECTING VALUE LOWER THAN 1 SETS IT TO ZERO !!
        call readchar("# ENPART","EXACT",iexact)
        call readchar("# ENPART","HOMO",ihomo)
        call readchar("# ENPART","DEKIN",idek)
        call readchar("# ENPART","IONIC",iionic)
        call readreal("# ENPART","TWOELTOLER",twoeltoler,0.00d0,1)
        call readchar("# ENPART","ANALYTIC",ianalytical)

!! ADDING GRID TUNNING FOR TWO-EL INTEGRATION !!
        call readchar("# ENPART","MOD-GRIDTWOEL",iigrid)
        if(iigrid.eq.1) then
          call readint("# GRID","RADIAL",nrad22,150,1)
          call readint("# GRID","ANGULAR",nang22,590,1)
          call readreal("# GRID","rr00",rr0022,0.5d0,1)
          call readreal("# GRID","phb1",phb12,0.169d0,1)
          call readreal("# GRID","phb2",phb22,0.170d0,1)
          call readreal("# GRID","THRESH2",thr3,1.0d-12,1)
        
!! DEFAULTS, MODIFIED FOR SAFE INTEGRATION SETUP !! 
        else
          nrad22=150
          nang22=590
          rr0022=0.5
          phb12=0.169d0
          phb22=0.170d0
          thr3=1.0d-12

!! Warn if a # GRID block exists in the input but is being ignored because !!
!! MOD-GRIDTWOEL wasn't set. Scan the file directly here rather than via   !!
!! "locate", which always prints a "section not found" message on a miss   !!
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
234      continue
          if(igridpresent.eq.1) then
            write(*,*) " "
            write(*,*) "WARNING: a # GRID block was found in the input, but MOD-GRIDTWOEL"
            write(*,*) "was not set in # ENPART. Using the default integration setup instead"
            write(*,*) "Add MOD-GRIDTWOEL to # ENPART to apply your # GRID settings"
            write(*,*) " "
          end if
        end if

!! FOR TOPOLOGY CALCULATION !!
!! MG: needs to be properly checked... done long time ago !!
!! MG: extended version for 2d, 3d, and more 1d topology in apost3.1-devel of my user... we should check if worth merging !!
        itop=0
        ipairs=0
        call readchar("# METHOD","TOPOLOGY",itop)
        if(itop.eq.1) then
          call locate(16,"# ATOM PAIRS DEFINITION",ii)
          if(ii.eq.0) stop " ATOM PAIRS DEFINITION SECTION MISSING "
          read(16,*) ipairs
          if(ipairs.gt.0) then
            do ii=1,ipairs
              read(16,*) (iatpairs(kk,ii),kk=1,2)
              write(*,*) (iatpairs(kk,ii),kk=1,2)
            end do
          else
            write(*,*) " DOING CUBE OF THE ENTIRE MOLECULAR SYSTEM "
          end if
!! FOR CHOOSING THE ENERGY COMPONENT TO DO THE TOPOLOGY !!
          ietop=-1
          call readchar("# TOPOLOGY","EXCHANGE",itop2)
          if(itop2.eq.1) ietop=1
          call readchar("# TOPOLOGY","CORRELATION",itop2)
          if(itop2.eq.1) ietop=2
          call readchar("# TOPOLOGY","EXCHANGE-CORRELATION",itop2)
          if(itop2.eq.1) ietop=3
          call readchar("# TOPOLOGY","DENSITY",itop2)
          if(itop2.eq.1) ietop=9
          if(ietop.eq.-1) stop " FUNCTION FOR TOPOLOGY NOT INTRODUCED "
        end if

!! END OF ENPART OPTIONS !!
      end if

!! EDAIQA OPTIONS !!
      iedaiqa=0
      iflip=0
      call readchar("# METHOD","EDAIQA",iedaiqa)
      if(iedaiqa.eq.1) then
        ii=0
        call locate(16,"# EDAIQA",ii)
        if(ii.eq.1) then
          read(16,'(a80)') namefchk1
          read(16,'(a80)') namefchk2
          open(unit=55,file=namefchk1)
          open(unit=52,file=namefchk2)
          call readchar("# EDAIQA","FLIPSPIN",iflip) !! ALPHA FOR BETA !!

!! ADDING pySCF REF VALUES FOR ELSTAT CALCULATION !!
          call readreal("# EDAIQA","eN pySCF",xen,0.0d0,1)
          call readreal("# EDAIQA","Coul pySCF",xcoul,0.0d0,1)
          call readreal("# EDAIQA","NN pySCF",xnn,0.0d0,1)

!! ADDING GRID TUNNING FOR TWO-EL INTEGRATION !!
!! Potser es fa repetitiu amb el de la seccio ENPART. Es podria fer un 2x1!!
          call readchar("# EDAIQA","MOD-GRIDTWOEL",iigrid)
          call readint("# GRID","RADIAL",nrad22,40,1)
          call readint("# GRID","ANGULAR",nang22,146,1)
          call readreal("# GRID","rr00",rr0022,0.5d0,1)
          call readreal("# GRID","phb1",phb12,0.162d0,1)
          call readreal("# GRID","phb2",phb22,0.182d0,1)

!! OPTIONS TO MAKE 2D PLOTS ABOUT ELECTROSTATIC POTENTIALS !!
          i2deda=0
          call locate(16,"# 2D PLOTS",i2deda)
          if(i2deda.eq.1) then
            read(16,*) iipoints
            read(16,*) (xptxyz(1,j),j=1,3)
            read(16,*) (xptxyz(2,j),j=1,3)

!! PRINTING TO ENSURE... CAN BE REMOVED !!
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
          stop "EDAIQA SECTION MISSING. REVISE inp"
        end if

!! END OF EDAIQA !!
      end if

C NLOPs                       
      call readchar("# METHOD","POLAR",ipolar )
      if(ipolar.eq.1) then
        iaccur=1
c using file $name.scr as raw output for post-processing
!MMO- deleting everything scr-related as it's no longer used
      end if

      call field_misc(ifield)
      if(ifield.eq.1) then
        write(*,*) 'The system is under a static electric field' 
        write(*,'(3(a4,f8.6))') 'Fx=',field(2), 'Fy=',field(3),'Fz=',field(4) 
      end if

C Do for restricted number of atoms
      idoat=0
      call readchar("# METHOD","DOATOMS",idoat)
      if(idoat.eq.1) then
        call locate(16,"# ATOMS",ii)
        if(ii.eq.0) stop 'Required section not found in input file'
        read(16,*) icuat
        read(16,*) (iatlist(i),i=1,icuat)
      else
        icuat=nat
        do i=1,icuat
          iatlist(i)=i
        end do
      end if

c Do for fragments
      idofr=0
      call readchar("# METHOD","DOFRAGS",idofr)
      if(idofr.eq.1) then
        call locate(16,"# FRAGMENTS",ii)
        if(ii.eq.0) stop 'Required section not found in input file'
        read(16,*) icufr
        do i=1,icufr
          read(16,*) nfrlist(i)
          if(i.eq.icufr.and.nfrlist(i).eq.-1) then 
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
          else     
            read(16,*) (ifrlist(k,i),k=1,nfrlist(i))
          end if
        end do
c
        ixx=0
        do i=1,icufr
          write(*,*)'Fragment: ',i 
          write(*,'(20i4)') (ifrlist(k,i),k=1,nfrlist(i))
          ixx=ixx+nfrlist(i)
        end do
        if(ixx.ne.nat.and.(ieos.eq.1.or.ienpart.eq.1)) then
          stop 'Missing/Additional atoms in fragment definition'
        end if

c  jfrlist tells which fragment a given atom belongs to
        do i=1,icufr
          do k=1,nfrlist(i)
            jfrlist(ifrlist(k,i))=i
          end do
        end do
        do i=1,nat
          if(jfrlist(i).eq.0.and.(ieos.eq.1.or.ienpart.eq.1)) then
            write(*,*) 'Unassigned atom to fragment:',i
            stop
          end if
        end do
c for compatibility
      else
        icufr=nat
        do i=1,icufr
          nfrlist(i)=1
          ifrlist(1,i)=i                
          jfrlist(i)=i                
        end do
      end if

!! EXTRA WARNINGS FOR EDAIQA !!
      if(iedaiqa.eq.1) then
        if(idofr.eq.0) then
          write(*,*) " FRAGMENT DEFINITION REQUIRED FOR EDAIQA "
          write(*,*) " INFO : FRAGMENT ORDER MUST MATCH THE ISOLATED "
          stop
        end if
        if(idofr.eq.1.and.icufr.ne.2) stop " ONLY 2 FRAGMENTS ALLOWED FOR EDAIQA "
      end if

C READING DM1 and DM2  
      if(icorr.ne.0) then 
        call locate(16,"# DM",ii)
        if(ii.eq.0) stop " # DM section not found in input file "
        call readchar("# DM","pySCF",ipyscf)
        call readchar("# DM","ORCA",iorca)
        call readchar("# DM","DMRG",idmrg)

        call locate(16,"# DM",ii)
        read(16,'(a80)') namedm
        if(iorca.eq.1.or.ipyscf.eq.1) then
          open(11,file=namedm)
        else
          open(11,file=namedm,FORM='UNFORMATTED',status='OLD')
        end if
        if(icorr.eq.2) then
          read(16,'(a80)') namedm
          if(iorca.eq.1.or.ipyscf.eq.1) then
            open(12,file=namedm)
          else
            open(12,file=namedm,FORM='UNFORMATTED',status='OLD')
          end if
        end if
      end if

!! OSLO OPTIONS !!

      ioslo=0
      call readchar("# METHOD","OSLO",ioslo)
      if(ioslo.eq.1) then

!! MG: BY DEFAULT REQUIRED THE TFVC AIM IN # METHOD (NUMERICAL INTEGRATION). BUT ONE CAN ASK OSLOs USING HILBERT AIMS !!
!! HILBERT AIMS CASES !!

        ilow2=0
        call readchar("# OSLO","MULLIKEN",ii)
        if(ii.eq.1) ilow2=1
        call readchar("# OSLO","LOWDIN",ii)
        if(ii.eq.1) ilow2=2
        call readchar("# OSLO","LOWDIN-DAVIDSON",ii)
        if(ii.eq.1) ilow2=3
        call readchar("# OSLO","NAO-BASIS",ii)
        if(ii.eq.1) ilow2=6

!! EXTRA OPTIONS !!
        
        call readint("# OSLO","FOLI TOLERANCE",ifolitol,3,1) !! FOLI VALUE TOLERANCE (FOR SELECTION) !!
        call readint("# OSLO","BRANCH ITERATION",ibranch,0,1) !! VALUE OF THE ITERATION TO INVOKE BRANCHING !!
        ioslofchk=1
        call readchar("# OSLO","PRINT NON-ORTHO",ii) !! FOR PRINTING NON-ORTHO OSLOs IN AN EXTRA .fchk FILE !!
        if(ii.eq.1) ioslofchk=2

      end if

!! TO PRINT .fchk FILES FROM QCHEM or MOKIT (MG: FCHK FORMAT IS DIFFERENT THAN GAUSSIAN) !!
      iqchem=0
      call readchar("# METHOD","QCHEM",iqchem)
      imokit=0
      call readchar("# METHOD","MOKIT",imokit)

!! END OF OSLO OPTIONS !!


!! X-RAY SCATTERING FACTORS !!
      call readchar("# METHOD","SCATT-FACT",iscattfact)

      end
