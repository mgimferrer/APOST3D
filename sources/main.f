c-----------------------------------------------------------------------------
c                                                                                  
c                        Program APOST-3D, Version 4                               
c                                 31-07-2023                                       
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
c          Topological fuzzy Voronoi cells (TFVC), J Chem Phys, 139 071103 2013    
c          QTAIM, J. Comput. Chem 30 1082 2009                                     
c                                                                                  
c        Hilbert-space : Mulliken, Lowdin, Davidson-Lowdin                         
c                                                                                  
c                                                                                  
c        Calculating:                                                              
c        ------------                                                              
c                                                                                  
c          A)  Atomic and overlap populations, bond orders and valences            
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
c             E. Ramos-Cordoba et.al., J. Chem. Phys. 138 214107 2013              
c                                                                                  
c          F) Local spin analysis                                                  
c             E. Ramos-Cordoba et.al., J. Chem. Theor. Comput. 8, 1270-1279 2012   
c             E. Ramos-Cordoba et.al., Phys. Chem. Chem. Phys. 14 15291-15298 2012 
c                                                                                  
c          G) Effective Oxidation states analysis                                  
c             E. Ramos-Cordoba et.al., J. Chem. Theor. Comput. 11 1501-1508 2015   
c
c          H) Oxidation states from localized orbitals
c             M. Gimferrer et al., J. Chem. Theor. Comput. 18 309-322 2022
c
c          I) Decomposition of EDA quantities into one- and two-center IQA terms
c             M. Gimferrer et al., J. Chem. Theory Comput. 19 3469-3485 2023
c                                                                                  
c                                                                                  
c        Cite this program as:                                                     
C        ---------------------                                                     
c                P. Salvador, E. Ramos-Cordoba, M. Gimferrer, M. Montilla          
c                Program APOST-3D, Version 4, Girona, 2020                       
c                                                                                  
c      e-mail: psalse@gmail.com, eloy.raco@gmail.com, mgimferrer18@gmail.com       
C                                                                                  
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
c                     (see http://www.tddft.org/programs/libxc                     
c                                                                                  
c      We are extremely grateful for the possibility of using these routines!      
c      -----------------------------------------------------------------------------c
      use basis_set    
      use ao_matrices
      use integration_grid
      implicit real*8(a-h,o-z)
      include 'parameter.h'
c general parameters
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
c orbitals and density matrices
c atom and fragment lists
      common /atlist/iatlist(maxat),icuat
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
c populations and EOS
      common /loba/ oxi(maxat),errsav(maxat),elec(maxat),effpop(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
c local spin
      common /localspin/xlsa(maxat,maxat),ua(maxat)
c Enpart
      common /exchg/exch(maxat,maxat),xmix
c TFVC features
      common /achi/achi(maxat,maxat),ibcp
      common /erf/aerf,ierf
!! MODGRID (diatXC) features !!
      common /modgrid/nrad22,nang22,rr0022,phb12,phb22
      common /modgrid2/thr3
!! EDAIQA features !!
      common /edaiqa/xen,xcoul,xnn
      common /edaiqa2/i2deda,iipoints,xptxyz(2,3)
c NLOP features
      common /twoel/twoeltoler
      common /efield/field(4),edipole
c printing and internal options
      common /filename/name0
      common /printout/iaccur
      common /iops/iopt(100)
c for enpart
      dimension eto(maxat,maxat)
c for topology
      dimension iatpairs(2,maxat)
c auxiliary arrays
      dimension dummyvec(maxat,2),navect(maxat)
      character*60 name,name2,namepat,name3,name0
      character*80 linia ,line, namedm
!! EDAIQA extra !!
      character*80 namefchk1,namefchk2
c for testing
      dimension xhess(3,3)
C
      allocatable wp(:),omp(:),omp2(:,:),chp(:,:),pcoord(:,:),rho(:)
      allocatable xkdens(:)
      allocatable ibaspoint(:)
      allocatable pca(:,:),scr(:)
      allocatable sss(:,:),sssi(:,:)
c
      allocatable sat(:,:,:)
      allocatable dm1(:,:),dm2(:,:,:,:)
C
!! TO CHANGE !!
      allocatable orbpop(:)
      allocatable rho_orb(:,:)
      allocatable rho_at(:,:)
!! !!


      call cpu_time(time)

CCCCCCCCCCCCCCCCCCCCCC
C PROCESSING ARGUMENTS
CCCCCCCCCCCCCCCCCCCCC

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

C print version info 
      call kiir()

CCCCCCCCCCCCCCCCCCCCC
C PROCESSING INP FILE 
CCCCCCCCCCCCCCCCCCCCC

      open (16,file=name3,err=9999)
      open (15,file=name,err=9999)

c choose density from fchk file
      call readint("# METODE","DENS",ndens0,1,1)
      iopt(9) = ndens0
      
C Processing FChk file
      call input()

      idoint=0
      iwfn=0  

c look for options      
      call readchar("# METODE","WFN",iwfn)
      call readchar("# METODE","ALLPOINTS",iallpo)
      call readchar("# METODE","FULLPRECISION",iaccur)

C Atoms in molecules
      call readchar("# METODE","MULLI",imulli)
      call readchar("# METODE","LOWDIN",ilow)
      if(ilow.eq.1) imulli=2
      call readchar("# METODE","LOWDIN-DAVIDSON",ilow)
      if(ilow.eq.1) imulli=3
      call readchar("# METODE","NAO-BASIS",ilow)
      if(ilow.eq.1) imulli=4
      call readchar("# METODE","LOWDIN-W",ilow)
      if(ilow.eq.1) imulli=5
      call readchar("# METODE","HIRSH",ihirsh)
      call readchar("# METODE","HIRSH-IT",ihirsh0)
      if(ihirsh0.eq.1) ihirsh=2
      call readchar("# METODE","BECKE-RHO",ibcp)
      call readchar("# METODE","NEWBEC",inewbec)
      call readchar("# METODE","TFVC",itfvc)
      if(itfvc.eq.1) then
        ibcp=1
        inewbec=1
        istiff=4
      end if
      call readint("# METODE","STIFFNESS",istiff,4,1)
      call readchar("# METODE","WMATRAD",iradmat)
      call readchar("# METODE","RMATRAD",iradmat0)
      if(iradmat0.eq.1) iradmat=2

      call readchar("# METODE","ERF_PROF",ierf)
      call readreal("# METODE","ERF_PROF",aerf,6.266d0,1)
      if(ierf.eq.1)  istiff=0

c ERC QTAIM input module
      call readchar("# METODE","QTAIM",iqtaim)
      call readchar("# METODE","READINT",ireadint)
      if(ireadint.eq.1) iqtaim=2
      if(iqtaim.eq.1) then
       call readint("# QTAIM","STEP",istep,300,1)
       call readint("# QTAIM","NNA",inna,0,1)
       call readint("# QTAIM","MAXDIST",imaxdist,12,1)
       call readint("# QTAIM","SCREENING",iscreening,1000,1)
       call readint("# QTAIM","PATH",ipath,0,1)
      end if 

c Miscellaneous options
      call readchar("# METODE","OPOP",iopop)
      call readchar("# METODE","SHANNON",isha)
      call readchar("# METODE","DOINT",idoint)
      call readchar("# METODE","PCA",ipca)
      call readchar("# METODE","LAPLACIAN",ilaplacian)
      call readchar("# METODE","FINEGRID",ifinegrid)
      call readchar("# METODE","ELCOUNT",ielcount) !MMO- NCTAIM
      call readint("# METODE","RHO_CALC_AT",iatdens,0,1)
      call readreal("# METODE","RHO_CALC_RAD",Rmax,0,1) ! MG: typo here!
      call readchar("# METODE","NOPOPU",inopop)
 
C eff-AO-s and EOS
      call readchar("# METODE","EFFAO",ieffao)
      call readchar("# METODE","UEFFAO",idummy)
      if(idummy.eq.1) ieffao=2
!! TO CHANGE FOR UEOS
      call readchar("# METODE","EFFAO-U",idummy)
      if(idummy.eq.1) ieffao=3
!!
      call readint("# METODE","EFF_THRESH",ieffthr,1,1)
      call readchar("# METODE","CUBE",icube)
      if(icube.eq.1) then
       call locate(16,"# CUBE",ii)
       if(ii.eq.0) stop'Required section # CUBE not found in input file'
       call readint("# CUBE","MAX_OCC",jcubthr,1000,1)
       call readint("# CUBE","MIN_OCC",kcubthr,0,1)
      end if
      call readchar("# METODE","EOS",ieos)
      if(ieos.eq.1) then 
       iopop=1
       call readreal("# METODE","EOS_THRESH",xthresh,2.5d-3,1)
       ieffao=2
      end if
      call readchar("# METODE","OS-CENTROID",ieoscent)

c Local spin and methods for correlated WFs
      call readchar("# METODE","SPIN",ispin)
      call readint("# METODE","DM",icorr,0,1)
      if(icorr.eq.2) ispin=1
      call readchar("# METODE","DAFH",idafh)

c Energy decomposition options   
      call readchar("# METODE","ENPART",ienpart )
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
        call readint("# GRID","RADIAL",nrad22,150,1)
        call readint("# GRID","ANGULAR",nang22,590,1)
        call readreal("# GRID","rr00",rr0022,0.5d0,1)
        call readreal("# GRID","phb1",phb12,0.169d0,1)
        call readreal("# GRID","phb2",phb22,0.170d0,1)
        call readreal("# GRID","THRESH2",thr3,1.0d-12,1)

!! FOR TOPOLOGY CALCULATION !!
!! MG: needs to be properly checked... done long time ago !!
!! MG: extended version for 2d, 3d, and more 1d topology in apost3.1-devel of my user... we should check if worth merging !!
        itop=0
        ipairs=0
        call readchar("# METODE","TOPOLOGY",itop)
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
      call readchar("# METODE","EDAIQA",iedaiqa)
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
      call readchar("# METODE","POLAR",ipolar )
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
      call readchar("# METODE","DOATOMS",idoat)
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
      call readchar("# METODE","DOFRAGS",idofr)
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
      call readchar("# METODE","OSLO",ioslo)
      if(ioslo.eq.1) then

!! MG: BY DEFAULT REQUIRED THE TFVC AIM IN # METODE (NUMERICAL INTEGRATION). BUT ONE CAN ASK OSLOs USING HILBERT AIMS !!
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

!! TO PRINT .fchk FILES FROM QCHEM (MG: FCHK FORMAT IS DIFFERENT THAN GAUSSIAN) !!

      iqchem=0
      call readchar("# METODE","QCHEM",iqchem)

!! END OF OSLO OPTIONS !!


CCCCCCCCCCCCCCCCCCCCCCCCC
C END PROCESSING INP FILE 
CCCCCCCCCCCCCCCCCCCCCCCCC


C kop=1 -> unrestricted calculation
C iposthf -> Correlated calculation
C icorr=0 -> no external dm1/2 provided
C icas=1 -> CASSCF calculation
C icisd -> CISD calculation
c idono=1 -> Do natural orbitals 
c idono=0 -> Restricted SD calculation 

CCCCCCCCCCCCCCCCCCCCCCC
C DEPENDENCIES & TO DO
CCCCCCCCCCCCCCCCCCCCCCC
      iposthf=0
      idono=0
      if(icas.eq.1.or.icisd.eq.1) iposthf=1
      iopt(65)=iposthf
      if(iposthf.eq.1.or.kop.eq.1) idono=1

      if(iposthf.eq.1) then
       if(ispin.eq.1.and.icorr.lt.2) stop ' Local Spin needs dm1 and dm2 for correlated WFs'
       if(ienpart.eq.1.and.icorr.lt.2) stop ' Enpart needs dm1 and dm2 for correlated WFs'
      end if
      if(iqtaim.eq.1) stop'This version can not do QTAIM'

      if (ipca.eq.1.and.iqtaim.ne.1) iopop=1
      if(imulli.gt.1.or.iqtaim.eq.1) iopop=0
      if(ihirsh.ne.0.and.idoatoms.eq.1) stop'Cant do HIRSH with DOATOMS'
      if(ispin.eq.1.and.idono.eq.0)  then
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

CCCCCCCCCCCCC
c OPTIONS LIST
CCCCCCCCCCCCC
      iopt(1) = idoint
c      iopt(2) = 
      iopt(3) = iwfn
      iopt(4) = idono
      iopt(5) = imulli
      iopt(6) = ihirsh
      iopt(7) = iallpo
      iopt(8) = iopop
c      iopt(9) = ndens0
c      iopt(10) = 
c     iopt(11) = 
      iopt(12) = ieffao
      iopt(13) = icube 
      iopt(14) = ibcp  
      iopt(15) = ispin 
      iopt(16) = iqtaim 
      iopt(17) = ienpart
      iopt(18) = ihf
!! !!
      iopt(20) = iexact
      iopt(21) = ihomo
      iopt(22) = idek
      iopt(23) = iionic
      iopt(24) = ieffthr
      iopt(25) = istiff 
      iopt(26) = icorr  
      iopt(27) = isha   
      iopt(28) = ipca 
c      iopt(29) =  ipnof
      iopt(30) = idafh  
      iopt(31) = inewbec
      iopt(34) = ilaplacian
c ERC qtaim 
      iopt(35) = istep 
      iopt(36) = inna
      iopt(37) = imaxdist
      iopt(38) = iscreening
      iopt(39) = ipath
c
      iopt(40) = idofr
      iopt(41) = jcubthr
      iopt(42) = kcubthr
      iopt(43) = iorca   
      iopt(44) = ithrebod
      iopt(45) = ifinegrid
c
      iopt(46) = ipolar
      iopt(47) = ifield
      iopt(48) = iradmat
      iopt(49) = ielcount !MMO- NCTAIM

      iopt(56) = ihirao !! MG: to be done !!
      iopt(58) = itop
      iopt(59) = iecorr
!! LIBXC IDs !!
      iopt(60) = id_xcfunc
      iopt(61) = id_xfunc
      iopt(62) = id_cfunc

      iopt(85) = ipairs
      iopt(86) = ietop
      iopt(87) = ipyscf   
      iopt(88) =  inopop
c      iopt(89) = 
      iopt(90) = ieoscent
      iopt(91) = ianalytical
!! EDAIQA !!
      iopt(92) = iedaiqa
      iopt(93) = iflip
!! IF QCHEM .fchk USED AS INPUT !!
      iopt(95) = iqchem
!! OSLO !!
      iopt(96) = ifolitol
      iopt(97) = ibranch
      iopt(98) = ioslofchk

CCCCCCCCCCCCC
c END OPTIONS LIST
CCCCCCCCCCCCC

c Natural orbitals from P-matrix in FChk
      if (idono.eq.1) call gennatural()

      write(*,*)
      call cpu_time(time2)
      write(*,'(a40,f10.1,a2)')'(Elapsed time :: initialization ',time2-time,'s)' 
      time=time2

C density at the iatdens atom
      if(iatdens.ne.0) call atdens_int(Rmax,iatdens) 
      if(inopop.eq.1) stop 'Normal termination of APOST3D'

c set atomic radii just in case needed
       call prepar()

c Do just Mulliken-type analysis...either int files for fcalc or effao
c or readint files 
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
     +   sis'
        call tolow2(sat)
       end if

!      else if(iqtaim.eq.2) then
!       ndim=nocc
!       allocate (sat(ndim,ndim,nat))
!       call readintfiles(sat)
       
       else
CCCCCCCCCCCCCCCCCCC
C PREPARE FOR NUMERICAL INTEGRATIONS
CCCCCCCCCCCCCCCCCCC

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

! Prepare for Becke
       if(ierf.ne.1) then
         write(*,'(A20,i5)') ' Using stiffness k :',istiff
       else
        write(*,'(a34,f8.4)') 'USING ADJUSTABLE PROFILE WITH A= ',aerf
       end if
       if(ibcp.eq.1) call khi()
       if(ihirsh.eq.1.or.ihirsh.eq.2.or.ielcount.eq.1) call makeatdens
 
       iiter=1
       call prenumint(ndim,itotps,nat,wp,omp,omp2,chp,rho,pcoord,ibaspoint,iiter)

CCCCCCCCCCCCCCC
C DO OS FROM LOCALIZED MOs
CCCCCCCCCCCCCCC

      if(ieoscent.eq.1) then
        write(*,*) 'Doing OS from localized orbitals...'
        call eos_centroid(itotps,chp,wp,omp,pcoord)
        stop
      end if

CCCCCCCCCCCCCCCCC
C INTEGRATE ATOMIC OVERLAP BY DEFAULT
CCCCCCCCCCCCCCCCC
       allocate (sat(ndim,ndim,nat))
       call numint_sat(ndim,itotps,nat,wp,omp,omp2,chp,ibaspoint,sat)

      end if
CCCCCCCCCCCCCCCCCCC
C END OF PREPARATION FOR NUMERICAL INTEGRATIONS
CCCCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCCCC
c CORRELATED WFs INPUT
C  Needed dm1 and/or dm2 files 
CCCCCCCCCCCCCCCCC

c cas cisd specifications
      nelec=nalf+nb
      if(icorr.ne.0) then 
       print *,' -------------------------------'
       print *,'  POST-HARTREE-FOCK CALCULATION '
       print *,' -------------------------------'
       print *,' '

       if(icisd.eq.1) nspinorb=nbasis*2

       print *,'Number of core + active spin-orbitals : ',nspinorb
       print *,'Number of electrons : ',nelec
       print *,'Number of basis functions :',nbasis
       print *,' '

       print *,'DM1 input starts'
       allocate (dm1(nspinorb,nspinorb))
       call dm1input(dm1)
       if (icorr.eq.2) then
        norb=nspinorb/2
        allocate (dm2(norb,norb,norb,norb))
        call dm2input(idmrg,dm1,dm2)
       end if
      end if

CCCCCCCCCCCCCCCCC
c END CORRELATED WFs INPUT 
CCCCCCCCCCCCCCCCC


CCCCCCCCCCCCCCC
C DO POPULATION ANALYSIS
CCCCCCCCCCCCCCC
      print *,' '
      print *,' ---------------------------'
      print *,'  DOING POPULATION ANALYSIS '
      print *,'   Partial atomic charges'
      print *,'   Atomic spin densities '
      print *,'   Bond orders and Valences '
      print *,' ---------------------------'
      print *,' '

c Spin density 
      if(kop.ne.0.or.nalf.ne.nb) then
       tc1=0.d0
       tc2=0.d0
       do kat=1,nat
        x=0.d0
        do mu=1,igr
         do nu=1,igr
          x=x+ps(mu,nu)*sat(mu,nu,kat)
         enddo
        enddo
        qsat(kat,1)=x
        qsat(kat,2)=0.d0
       enddo
       do mu=1,igr
        kat=ihold(mu)
        x=0.d0
        do itau=1,igr
         x=x+ps(mu,itau)*s(itau,mu)
        enddo
        qsat(kat,2)=qsat(kat,2)+x
       enddo
      end if

c Total density 
      tc1=0.d0
      tc2=0.d0
      do kat=1,nat
       x=0.d0
       do mu=1,igr
        do nu=1,igr
         x=x+p(mu,nu)*sat(mu,nu,kat)
        enddo
       enddo
       qat(kat,1)=x
       qat(kat,2)=0.d0
      enddo
      do mu=1,igr
       kat=ihold(mu)
       x=0.d0
       do itau=1,igr
        x=x+p(mu,itau)*s(itau,mu)
       enddo
       qat(kat,2)=qat(kat,2)+x
      enddo


CCCCCCCCCCCCCCC
C Write INT FILES
CCCCCCCCCCCCCCC
!       write(*,*) 'idoint',idoint
      if(idoint.eq.1)  call print_int(ndim,nat,sat,namepat)

CCCCCCCCCCCCCCC
C ELECTRON POPULATIONS
CCCCCCCCCCCCCCC
      tc1=0.d0
      tc2=0.d0
      do i=1,nat
       tc1=tc1+qat(i,1)
       tc2=tc2+qat(i,2)
      enddo
      print *,'  '
      print *,'    ELECTRON POPULATIONS'
      print *,'  '
      print *,'  Atom   apost3d      Mulliken'
      print *,' -----------------------------'
      call vprint(qat,nat,maxat,2)
      print *,' -----------------------------'
      if(iaccur.eq.0) then
       print 162, tc1, tc2
      else
       print 172, tc1, tc2
      end if

      if (idofr.eq.1) then
       line ='   FRAGMENT ANALYSIS : Electron populations'
       call group_by_frag_vec(2,line ,qat)
      end if

CCCCCCCCCCCCCCC
C PARTIAL CHARGES     
CCCCCCCCCCCCCCC
      tc1=0.d0
      tc2=0.d0
      do i=1,nat
       tc1=tc1-qat(i,1)+zn(i)
       tc2=tc2-qat(i,2)+zn(i)
       dummyvec(i,1)=zn(i)-qat(i,1)
       dummyvec(i,2)=zn(i)-qat(i,2)
      enddo
      print *,'  '
      print *,'    TOTAL ATOMIC CHARGES    '
      print *,'  '
      print *,'  Atom   apost3d      Mulliken'
      print *,' -----------------------------'
      call vprint(dummyvec,nat,maxat,2)
      print *,' -----------------------------'
      if(iaccur.eq.0) then
       print 162, tc1, tc2
      else
       print 172, tc1, tc2
      end if
      print *,'  '
      
      if (idofr.eq.1) then
       line ='   FRAGMENT ANALYSIS : Atomic Charges'
       call group_by_frag_vec(2,line ,dummyvec)
      end if

      if(ielcount.eq.1) CALL NCTAIM(iatps,wp,omp,omp2,nat,sat,igr,ibaspoint,chp)

CCCCCCCCCCCCCCC
C SPIN POPULATIONS
CCCCCCCCCCCCCCC
      if(kop.ne.0.OR.(icas.eq.1.and.nalf.ne.nb.and.icorr.ne.0)) then
       tc1=0.d0
       tc2=0.d0
       do i=1,nat
        tc1=tc1+qsat(i,1)
        tc2=tc2+qsat(i,2)
       enddo
       print *,'  '
       print *,'     SPIN POPULATIONS'
       print *,'  '
       print *,'  Atom   apost3D      Mulliken'
       print *,' -----------------------------'
       call vprint(qsat,nat,maxat,2)
       print *,' -----------------------------'
       if(iaccur.eq.0) then
        print 162, tc1, tc2
       else
        print 172, tc1, tc2
       end if
       print *,'  '
      if (idofr.eq.1) then
       line ='   FRAGMENT ANALYSIS : Spin Populations'
       call group_by_frag_vec(2,line ,qsat)
      end if
      end if

CCCCCCCCCCCCCCC
C OVERLAP POPULATIONS
CCCCCCCCCCCCCCC
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
       print *,' '
       print *,'          APOST3D  OVERLAP POPULATION MATRIX'
       print *,' '
       call mprint(op,nat,maxat)
       print *,'  '
      if (idofr.eq.1) then 
       line ='   FRAGMENT ANALYSIS : Overlap Populations'
       call group_by_frag_mat(0,line,op)
      end if
      end if

CCCCCCCCCCCCCCC
c Bond orders, valences, number of eff. unpaired electrons
CCCCCCCCCCCCCCC
       call fborder(sat)
       if (idofr.eq.1) then
        line ='   FRAGMENT ANALYSIS : Fuzzy Bond Order'
        call group_by_frag_mat(1,line ,bo)
       end if
c
      write(*,*)
      call cpu_time(time2)
      write(*,'(a37,f10.1,a2)')'(Elapsed time :: population analyses ',time2-time,'s)' 
      time=time2

CCCCCCCCCCCCCCC
C Local Spin decomp for single-determinant WF
CCCCCCCCCCCCCCC
       if(ispin.eq.1.and.icas.eq.0.and.icisd.eq.0) then
        print *,' '
        print *,' ---------------------------'
        print *,'  DOING LOCAL SPIN ANALYSIS '
        print *,'   Single-determinant case  '
        print *,' ---------------------------'
        print *,' '
        call fspindec(sat)
        if (idofr.eq.1) then 
         line ='   FRAGMENT ANALYSIS : Local Spin Analysis'
         call group_by_frag_mat(0,line ,xlsa)
         line ='   FRAGMENT ANALYSIS : Num. eff. unpaired elec.'
         call group_by_frag_vec(1,line ,ua)
        end if
       end if

CCCCCCCCCCCCCCC
C LOCAL SPIN & DIs for CORRELATED WFs         
C  Needed dm1 and/or dm2 files produced by DMN code (E. Matito)
CCCCCCCCCCCCCCC
      if((icas.eq.1.or.icisd.eq.1).and.ispin.eq.1)then
       print *,' '
       print *,' ---------------------------'
       print *,'  DOING LOCAL SPIN ANALYSIS '
       print *,' Localization/Delocalization' 
       print *,'        Correlated WF       '
       print *,' ---------------------------'
       print *,' '
       call spincorr(sat,dm1,dm2)
       call cpu_time(time2)
       write(*,'(a37,f10.1,a2)')'(Elapsed time :: local spin analysis',time2-time,'s)' 
       time=time2
      end if

CCCCCCCCCCCCCCC
C NONLINEAR OPTICAL PROPERTIES          
CCCCCCCCCCCCCCC

      if(ipolar.ne.0) call polar(itotps,nat,wp,omp,omp2,pcoord,rho)

CCCCCCCCCCCCCCC
C ENTROPIES AND CORRELATION INDICATORS
CCCCCCCCCCCCCCC
c      if(isha.ne.0) call numint_sha(ndim,itotps,nat,wp,chp,omp,omp2,ibaspoint)
CCCCCCCCCCCCCCC
C PCA ANALYSIS...after DIs MAY HAVE BEEN CALCULATED
CCCCCCCCCCCCCCC
      if(ipca.eq.1) then
       allocate (pca(nat,natoms))
       allocate (scr(nat))
       do i=1,nat
        do j=1,nat
         pca(i,j)=op(i,j)-0.50d0*di(i,j)
        end do
       end do
       print *,' '
       print *,'                      "FUZZY ATOMS" COVARIANCE MATRIX'
       print *,' '
       call mprint(pca,nat,natoms)
       print *,'  '
       call diagonalize(nat,natoms,pca,scr,0)
       print *,' '
       print *,'                      "FUZZY ATOMS" PCA EIGENVECTORS '
       print *,' '
       call mprint(pca,nat,natoms)
       print *,' '
       write(*,'(8f10.4)') (scr(i),i=1,nat)
       print *,' '
       do i=1,nat
         xx=ZERO
         do k=1,nat
           xx=xx+pca(k,i)*qat(k,1)
         end do
         write(*,*) 'PC: ',i,' sum: ',xx,xx*scr(i)
       end do

      end if

CCCCCCCCCCCCCCC
C DAFH PART
C  Needed files produced by external code (R. Ponec)
C  Not available in this version
CCCCCCCCCCCCCCC
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

CCCCCCCCCCCCCCC
C EFFAO PART
CCCCCCCCCCCCCCC

c ieffao 0 nothing
c ieffao 1 eff-AOs 
c ieffao 2 spin-resolved eff-AOs (a must for EOS)
c ieffao 3 paired/unpaired eff-AOs

c imulli 1 MULLIKEN                                 
c imulli 2 LOWDIN                                   
c imulli 3 LOWDIN-DAVIDSON (not implemented)                                   
c imulli 3 NAO                                    

       if (ieffao.ne.0) then

        if(idofr.eq.0) then
         icufr=nat
         do i=1,icufr
          nfrlist(i)=1
          ifrlist(1,i)=i                
         end do
        end if

CCCCCCCCCCCCCCC
C Mulliken
CCCCCCCCCCCCCCC
        if(imulli.eq.1) then
          if(ieffao.eq.1) then
            call ueffaomull_frag(0)
          else if (ieffao.eq.2) then
           call ueffaomull_frag(1)
           if (ieos.eq.1) call eos_analysis(0,1,xthresh)

           if(kop.ne.0.or.(icas.eq.1.and.nalf.ne.nb))then
            call ueffaomull_frag(2)
            idobeta=1
           end if
           if (ieos.eq.1) call eos_analysis(idobeta,2,xthresh)
          end if
CCCCCCCCCCCCCCC
C Lowdin      
CCCCCCCCCCCCCCC
        else if(imulli.gt.1) then
          if(ieffao.eq.1) then
           call ueffaolow_frag(0)
          else if (ieffao.eq.2) then
           call ueffaolow_frag(1)
           if(ieos.eq.1) call eos_analysis(0,1,xthresh)

           if(kop.ne.0.or.(icas.eq.1.and.nalf.ne.nb))then
            call ueffaolow_frag(2)
            idobeta=1
           end if
           if (ieos.eq.1) call eos_analysis(idobeta,2,xthresh)
          end if
CCCCCCCCCCCCCCC
C 3D-space    
CCCCCCCCCCCCCCC
        else

C DO EFFAO/UEFFAO for selected atoms only. NO EOS

         if(idoat.ne.0) then
          if(ieffao.eq.1) then
            call ueffao3d(itotps,ndim,omp,chp,sat,wp,omp2,p,0)
          else
            call ueffao3d(itotps,ndim,omp,chp,sat,wp,omp2,pa,1)
            call ueffao3d(itotps,ndim,omp,chp,sat,wp,omp2,pb,2)
          end if

         else

C DO EFFAO/UEFFAO for fragments/all atoms

          if(ieffao.eq.1) then
           call ueffao3d_frag(itotps,ndim,omp,chp,sat,wp,omp2,p,0)
c           call  uefomo(itotps,ndim,omp,chp,sat,wp,omp2,0)
c           call mhg(itotps,ndim,omp,chp,sat,wp,omp2,pcoord,pa,0)
c           call mhg2(itotps,ndim,omp,chp,sat,wp,omp2,pcoord,p,0)
          else if (ieffao.eq.2) then
           idobeta=0
           write(*,*) '  '
           write(*,*) 'UEFFAO: alpha and beta treated separatedly'
           write(*,*) '  '
           call ueffao3d_frag(itotps,ndim,omp,chp,sat,wp,omp2,pa,1)
           if (ieos.eq.1) call eos_analysis(idobeta,1,xthresh)

           if(kop.ne.0.or.(icas.eq.1.and.icass.ne.0))then
            idobeta=1
            call ueffao3d_frag(itotps,ndim,omp,chp,sat,wp,omp2,pb,2)
           end if
           if (ieos.eq.1) call eos_analysis(idobeta,2,xthresh)
          else if (ieffao.eq.3) then
           write(*,*) '  '
           write(*,*) 'EFFAO: from paired and unpaired density '
           write(*,*) '  '
           call effao3d_u(itotps,ndim,omp,chp,sat,wp,omp2)
         end if
        end if
       end if

      write(*,*)
      call cpu_time(time2)
      write(*,'(a37,f10.1,a2)')'(Elapsed time :: effAOs/EOS analysis ',time2-time,'s)' 
      time=time2
      end if

CCCCCCCCCCCCCCC
C DO ENERGY DECOMPOSITION
CCCCCCCCCCCCCCC

      if(ienpart.eq.1) then
      print *,' '
      print *,' --------------------------------------'
      print *,'  DOING MOLECULAR ENERGY DECOMPOSITION '
      print *,' --------------------------------------'
      print *,' '

!! CASSCF AND CISD ENERGY DECOMPOSITION !!

      if(iposthf.eq.1) then
        call numint_one_rphf(ndim,itotps,wp,omp2,pcoord,chp,rho,eto)
        call cpu_time(time2)
        write(*,'(a39,f10.1,a2)')'(Elapsed time :: enpart one-electron ',time2-time,'s)' 
        time=time2

!! DFT AND HF ENERGY DECOMPOSITION !!

      else 

!! INITIALIZE DFT FUNCTIONAL FOR INFO AND INITIAL PRINTING !!
      
      if(id_xfunc.ne.-1) then
        if(id_xcfunc.ne.0) call func_info_print(id_xcfunc,itype)
        if(id_cfunc.ne.0) call func_info_print(id_cfunc,itype)
        if(id_xfunc.ne.0) call func_info_print(id_xfunc,jtype)
        if(itype.ge.jtype) iopt(55) = itype
        if(jtype.gt.itype) iopt(55) = jtype
        write(*,*) "itype,jtype,iop : ",itype,jtype,iopt(55)
      else
        write(*,*) "INFO : HF functional selected"
        write(*,*) "INFO : HF-type exchange coeff : ",xmix
        write(*,*) " "
      end if

!! INTEGRATING ONE-ELECTRON CONTRIBUTIONS !!

      if(kop.ne.1) then 
        ALLOCATE(xkdens(itotps)) ! (TO DO) Rethink how to include it... only used in metaGGA functionals

!! RESTRICTED ONE-ELECTRON PART !!

        call numint_one(ndim,itotps,wp,rho,omp,omp2,pcoord,chp,eto)
        write(*,*) " "
        call cpu_time(time2)
        write(*,'(a37,f10.1,a2)')'(Elapsed time :: enpart one-electron ',time2-time,'s)' 
        time=time2

!! RESTRICTED DFT XC PART !!

        if(id_xfunc.ne.-1) then
          if(ianalytical.eq.0) then
            call numint_dft(ndim,itotps,wp,rho,omp,omp2,chp,eto,pcoord,sat)
          else
            call numint_dft_analytical(ndim,itotps,wp,omp,omp2,chp,eto,pcoord,sat)
          end if
          write(*,*) " "
          call cpu_time(time2)
          write(*,'(a28,f10.1,a2)')'(Elapsed time :: enpart dft ',time2-time,'s)' 
          time=time2
        end if
        DEALLOCATE(xkdens)
      else
        ALLOCATE(xkdens(itotps)) ! (TO DO) Rethink how to include it... only used in metaGGA functionals

!! UNRESTRICTED ONE-ELECTRON PART !! 

        call numint_one_uhf(ndim,itotps,wp,rho,omp,omp2,pcoord,chp,eto)
        write(*,*) " "
        call cpu_time(time2)
        write(*,'(a37,f10.1,a2)')'(Elapsed time :: enpart one-electron ',time2-time,'s)' 
        time=time2

!! UNRESTRICTED DFT XC PART !!

        if(id_xfunc.ne.-1) then
          call numint_dft_uks(ndim,itotps,wp,omp,omp2,chp,pcoord,eto)
          write(*,*) " "
          call cpu_time(time2)
          write(*,'(a37,f10.1,a2)')'(Elapsed time :: enpart one-electron ',time2-time,'s)' 
          time=time2
        end if
        DEALLOCATE(xkdens)
      end if 

!! END IF METHOD TYPE !!

      end if
      DEALLOCATE(wp,omp,omp2,chp,pcoord,ibaspoint,rho)

!! END OF ONE-ELECTRON PART !!
!! INTEGRATING TWO-ELECTRON CONTRIBUTIONS !!
!! TWO-ELECTRON INTEGRATION DEFAULTS !!

      print *,'SETTING GRID FOR TWO-ELECTRON ENERGY INTEGRATION'
!! CONTROLLED BY # GRID OPTION (modgrid common) !!
      nrad=nrad22
      nang=nang22
      rr00=rr0022
!      nrad=40
!      nang=146
!! DEFAULT GRID IS LARGER THAN THIS ONE !!
!! MG: CAN BE CHANGED TO THE 40 146, BUT ENSURE TO CHANGE ALSO DEFAULTS FOR pha AND phb !!
      if(ianalytical.eq.1) then
        write(*,*) 'Analytical calculation has been requested: Increasing grid because 2-electron is now 1-electron.'
        nrad=70
        nang=434
      end if

      print *,'  '
      print *,'Radial points: ', nrad
      print *,'Angular points:',nang    


      iatps=Nang*NRad
      itotps=nrad*nang*nat
      call quad(Nrad,Nang) 
      allocate (wp(itotps),omp(itotps),omp2(itotps,nat),rho(itotps))
      allocate (pcoord(itotps,3),ibaspoint(itotps),chp(itotps,ndim))

!! GENERATING GRID FOR TWO-EL CONTRIBUTIONS !!

      call prenumint(ndim,itotps,nat,wp,omp,omp2,chp,rho,pcoord,ibaspoint,0)

!! POST-HF TWO-ELECTRON PART !!

      if(iposthf.eq.1) then
        if(itop.eq.1) call top_3d(norb,2,0,iatpairs) !! MG: needs a check !!
        call numint_two_rphf(ndim,itotps,wp,omp2,pcoord,chp,rho,eto,dm1,dm2)
        write(*,*) " "
        call cpu_time(time2)
        write(*,'(a37,f10.1,a2)')'(Elapsed time :: enpart two-electron ',time2-time,'s)' 
        time=time2
      else

!! SINGLE-DETERMINANT UNRESTRICTED TWO-ELECTRON PART !!

        if(kop.eq.1) then 
          call numint_two_uhf(ndim,itotps,wp,omp,omp2,pcoord,chp,rho,eto)
        else

!! SINGLE-DETERMINANT RESTRICTED TWO-ELECTRON PART !!

          if(itop.eq.1) call top_3d(nocc,1,0,iatpairs) !! MG: needs a check !!
          call numint_two(ndim,itotps,wp,omp,omp2,pcoord,chp,rho,eto)
        end if 

!! END OF TWO-ELECTRON PART !!

        write(*,*) " "
        call cpu_time(time2)
        write(*,'(a37,f10.1,a2)')'(Elapsed time :: enpart two-electron ',time2-time,'s)' 
        time=time2

!! END IF OF WF-TYPE !!
      end if 

!! DEALLOCATING ANY EXISTENT GRID !!
      DEALLOCATE(wp,omp,omp2,chp,pcoord,ibaspoint,rho)

!! END IF OF IENPART !!
      end if

CCCCCCCCCCCCCCCCCCC
!! END OF ENERGY DECOMPOSITION !!
CCCCCCCCCCCCCCCCCCC

!! MG: MORE TOPOLOGY STUFF CAN BE FOUND IN VERSION 3.1-DEVEL !!

CCCCCCCCCCCC
!! EDAIQA !!
CCCCCCCCCCCC

!! CREATING THE E(< A^0 B^0 >)^AB STATE !!

!     if(iedaiqa.eq.1) then
!       write(*,*) " "
!       write(*,*) " Entering EDAIQA Section "
!       write(*,*) " "
!       if(kop.eq.0) call rwf_edatoiqafchk()
!       if(kop.eq.1) call uwf_edatoiqafchk()
!       close(51)
!       close(52)
!     end if

CCCCCCCCCCCCCCCCCCC
!! END OF EDAIQA !!
CCCCCCCCCCCCCCCCCCC


!! MG: HIRAO PART TO BE DONE... LOCATED IN ANOTHER VERSION !!
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
!! DFT DM1 APPROXIMATIONS !!
CCCCCCCCCCCCCCCCCCCCCCCCCCCC

!! GENERATING FIRST GRID HERE FOR NUMERICAL INTEGRATION !!
!     ndim=igr
!     iatps=nang*nrad
!     itotps=nrad*nang*nat
!     pha=ZERO
!     phb=ZERO
!     call quad(nrad,nang)
!     ALLOCATE(wp(itotps),omp(itotps),omp2(itotps,nat))
!     ALLOCATE(pcoord(itotps,3),ibaspoint(itotps),chp(itotps,igr))
!     call prenumint(ndim,itotps,nat,wp,omp,omp2,chp,pcoord,ibaspoint,0)
!     DEALLOCATE(ibaspoint,omp)
!     if(ihirao.eq.1) call dft_dm1(itotps,wp,omp2,pcoord,chp) 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!! END DFT DM1 APPROXIMATIONS !!
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!! (MG: variants of the oslo procedure can be found in my development version) !!
CCCCCCCCCC
!! OSLO !!
CCCCCCCCCC

      if(ioslo.eq.1) then
        write(*,*) " "
        write(*,*) " ********************* "
        write(*,*) " ENTERING OSLO SECTION "
        write(*,*) " ********************* "
        write(*,*) " "

!! COMPUTING sat FOR HILBERT SPACE CASES !!

        if(ilow2.ne.0) then
          iopt(5)=ilow2 !! MG: trick !!
          if(ilow2.eq.1) call tomull(sat)
          if(ilow2.eq.2.or.ilow2.eq.3) call tolow(sat)
          if(ilow2.eq.6) call tonao(sat)
        end if

!! GENERAL INDEPENDENTLY OF THE AIM !!

        if(kop.eq.0) then 
          call rwf_iterative_oslo(sat,itotps,wp,omp2,chp,pcoord)
        else 
          call uwf_iterative_oslo(sat,itotps,wp,omp2,chp,pcoord)
        end if
        DEALLOCATE(wp,omp,omp2,chp,pcoord)
      end if

CCCCCCCCCCCCCCCCC
!! END OF OSLO !!
CCCCCCCCCCCCCCCCC


      print *,' '
      print *,'...Normal Termination of APOST-3D... '

 162  format(1x,'    Sum  ',2(f10.6,2X))  
 172  format(1x,'    Sum  ',2(f20.13,2X))  
        
      end 

