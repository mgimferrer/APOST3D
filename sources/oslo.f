
!! **************************** !!
!! OSLO CALCULATION SUBROUTINES !!
!! **************************** !!

!! FIRST RESTRICTED SUBROUTINES !!

!! ****** !!

      subroutine rwf_iterative_oslo(sat,itotps,wp,omp2,chp,pcoord)

!! MG: ONLY ITERATIVE OSLO IMPLEMENTED IN THIS VERSION !!
!! MG: NON-ITERATIVE PROCEDURE IMPLEMENTED IN DEVELOPMENT VERSION !!

      use basis_set
      use ao_matrices
      use integration_grid
      implicit double precision(a-h,o-z)
      include 'parameter.h'

      integer,intent(in) :: itotps

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /filename/name0

      character*80 line
      character*60 name0,name1
      character*20 ctype

      dimension sat(igr,igr,nat)
      dimension chp(itotps,igr),pcoord(itotps,3)
      dimension wp(itotps),omp2(itotps,nat)

      allocatable :: Smat(:,:,:)
      allocatable :: S0(:,:),Sm(:,:),Splus(:,:),smh(:,:),eigv(:,:)
      allocatable :: c0(:,:),pp0(:,:)

      allocatable :: SSS(:,:),EEE(:,:)

      allocatable :: cmat(:,:),cfrgoslo(:,:,:)
      allocatable :: orbpop(:),orbpop2(:),frgpop(:,:),frgspr(:,:)
      allocatable :: foslo(:,:),foslo2(:,:)

      allocatable :: infopop(:,:),infooslo(:),poposlo(:)
      allocatable :: fspread(:),scr(:),deloc(:,:),delocoslo(:)

      allocatable :: clindep(:,:)

      allocatable :: pcore(:,:),pnocore(:,:)
      allocatable :: ccore(:,:),ccoreorth(:,:)

      allocatable :: coslo(:,:),cosloorth(:,:)
      allocatable :: poslo(:,:)

      allocatable :: ifrgel(:),iznfrg(:)

  
!! LOADING IOPTs !!
      iqchem   = iopt(95)
      ifolitol = iopt(96)
      ibranch  = iopt(97)
      ifchk    = iopt(98)

!! PARAMETERS REQUIRED !!
      iatps = nrad*nang
      niter = 999 !! MAX NUMBER OF ITERATIONS !!
      folitol = 10.0d0**(-REAL(ifolitol)) !! DEFAULT = 10^-3, CONTROLLED IN .inp !!

!! ALLOCATING GENERAL (REQUIRED) MATRICES !!
      ALLOCATE(SSS(nocc,nocc),EEE(nocc,nocc))
      ALLOCATE(Sm(igr,igr),Splus(igr,igr))
      ALLOCATE(c0(igr,igr),pp0(igr,igr))
      ALLOCATE(ccore(igr,igr),pcore(igr,igr),ccoreorth(igr,igr))
      ALLOCATE(coslo(igr,igr),cosloorth(igr,igr))
      ALLOCATE(clindep(igr,igr))
      ALLOCATE(fspread(nocc),scr(nocc),delocoslo(igr))
      ALLOCATE(infooslo(igr),poposlo(igr))
      ALLOCATE(ifrgel(icufr),iznfrg(icufr))

!! FRAGMENT CHARGE EXTRACTED FROM zn (AVOIDS PROBLEMS WHEN PSEUDOPOTENTIALS ARE USED) !!
      do ifrg=1,icufr
        izn=0
        do icenter=1,nfrlist(ifrg)
          iiat=ifrlist(icenter,ifrg)
          izn=izn+INT(zn(iiat))
        end do
        iznfrg(ifrg)=izn
        ifrgel(ifrg)=0
      end do

!! FRAGMENT LOCALIZATION: IMPORTANT :: R_A VALUE SELECTED IS THE CENTER OF CHARGE BETWEEN FRAGMENT ATOMS !!
      write(*,*) " --------------------------------------- "
      write(*,*) "  CHARGE CENTER (R_F) FOR EACH FRAGMENT  "
      write(*,*) " --------------------------------------- "
      write(*,*) " "
      write(*,*) "  Frag.        Charge center (xyz)     "
      write(*,*) " ------------------------------------- "
      ALLOCATE(Smat(ifrg,igr,igr))
      do ifrg=1,icufr
        xcenter=ZERO
        ycenter=ZERO
        zcenter=ZERO
        xchg=ZERO
        do icenter=1,nfrlist(ifrg)
          iiat=ifrlist(icenter,ifrg)
          xatchg=REAL(iznuc(iiat))
          xcenter=xcenter+xatchg*coord(1,iiat)
          ycenter=ycenter+xatchg*coord(2,iiat)
          zcenter=zcenter+xatchg*coord(3,iiat)
          xchg=xchg+xatchg
        end do
        xcenter=xcenter/xchg
        ycenter=ycenter/xchg
        zcenter=zcenter/xchg
        write(*,110) ifrg,xcenter*angtoau,ycenter*angtoau,zcenter*angtoau

!! COMPUTING S ONLY ONCE !!
        do mu=1,igr
          do nu=1,mu
            xx=ZERO
            do jcenter=1,nat
              do ifut=iatps*(jcenter-1)+1,iatps*jcenter
                distx=pcoord(ifut,1)-xcenter
                disty=pcoord(ifut,2)-ycenter
                distz=pcoord(ifut,3)-zcenter
                xoper=(distx*distx+disty*disty+distz*distz)
                xx=xx+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*omp2(ifut,jcenter)*xoper
              end do
            end do
            Smat(ifrg,mu,nu)=xx
            Smat(ifrg,nu,mu)=xx
          end do
        end do
      end do
      write(*,*) " ------------------------------------- "

!! PREPARING FOR ITERATIVE PROCESS !!
      ALLOCATE(cmat(igr,igr),cfrgoslo(icufr,igr,igr))
      ALLOCATE(frgpop(icufr,nocc))
      ALLOCATE(frgspr(icufr,nocc))
      ALLOCATE(pnocore(igr,igr))

!! ZEROING MATRICES LIKE THIS NOW !!
      coslo=ZERO
      cosloorth=ZERO
      pnocore=pa !! SAVING pa IN pnocore !!

!! INITIAL PRINTING !!
      write(*,*) " "
      write(*,*) " ----------------------------------- "
      write(*,*) "  STARTING ITERATIVE OSLO ALGORITHM  "
      write(*,*) " ----------------------------------- "
      write(*,*) " "
      write(*,'(2x,a50,f10.5)') "Tolerance (in delta-FOLI) used for OSLO selection:",folitol
      write(*,*) " "

      ntotcore=0
      iaddoslo=0
      iaddoslo2=0
      do iiter=1,niter
        write(*,*) " ---------------------- "
        write(*,'(3x,a16,x,i3)') "ITERATION NUMBER",iiter
        write(*,*) " ---------------------- "
        write(*,*) " "

!! ZEROING THE INVOLVED MATRICES !!
        pcore=ZERO
        ccore=ZERO
        ccoreorth=ZERO
        ALLOCATE(deloc(icufr,nocc))
        deloc=ZERO

!! 1) OBTAINING OSLOs FOR ALL FRGS !!
        iaddcore=0
        ALLOCATE(S0(igr,igr))
        do ifrg=1,icufr
          do ii=1,igr
            do jj=1,igr

!! USING pp0 AND S0 TO NOT DESTROY pnocore AND Smat !!
              pp0(ii,jj)=pnocore(ii,jj)
              S0(ii,jj)=Smat(ifrg,ii,jj)
            end do
          end do
          call build_Smp(igr,S0,Sm,Splus,0)
          call to_lowdin_basis(igr,Splus,pp0)
          call diagonalize(igr,igr,pp0,c0,0)
          call to_AO_basis(igr,igr,Sm,c0)
          imaxo=nocc !! MAX NUMBER OF OSLOs IS NOW nocc FOR PRACTICITY !!

!! RECOVERING THE COEFFS, SAVED IN cfrgoslo !!
          do kk=1,imaxo
            do mu=1,igr
              cmat(mu,kk)=c0(mu,kk)*dsqrt(pp0(kk,kk))
              cfrgoslo(ifrg,mu,kk)=cmat(mu,kk)
            end do
          end do

!! COMPUTING PIPEK DELOCALIZATION, REQUIRES FRAGMENT POPULATIONS !!
          ALLOCATE(orbpop(nocc))
          do jfrg=1,icufr
            call rwf_uwf_frg_pop(jfrg,sat,nocc,cmat,orbpop)
            do ii=1,imaxo
              deloc(ifrg,ii)=deloc(ifrg,ii)+orbpop(ii)*orbpop(ii) !! ADDING Q_A**2 INSIDE deloc !!

!! IMPORTANT HERE SAVING ONLY FOR THE OWN FRAGMENT !!
              if(jfrg.eq.ifrg) then
                frgspr(ifrg,ii)=pp0(ii,ii) !! SPREADS (ONLY FOR PRINTING) !!
                frgpop(ifrg,ii)=orbpop(ii) !! FRAGMENT POPULATIONS !!
              end if
            end do
          end do
          DEALLOCATE(orbpop)

!! NOW DOING 1/deloc() !!
          do ii=1,imaxo
            if(deloc(ifrg,ii).gt.1.0d-6) then
              deloc(ifrg,ii)=ONE/deloc(ifrg,ii)
            else
              deloc(ifrg,ii)=100.0d0 !! ABSURT VALUE, AVOIDS PROBLEMS !!
            end if
          end do

!! PRINTING VALUABLE INFORMATION !!
          write(*,*) " -------------------------------------- "
          write(*,'(3x,a32,x,i3)') "ORBITAL INFORMATION FOR FRAGMENT",ifrg
          write(*,*) " -------------------------------------- "
          write(*,*) " "
          write(*,*) "  Orb.    Spread      Frg. Pop.      FOLI   "
          write(*,*) " ------------------------------------------ "
          do ii=1,imaxo
            xx=frgpop(ifrg,ii)
            if(xx.gt.1.0d-5) write(*,111) ii,frgspr(ifrg,ii),xx,dsqrt(deloc(ifrg,ii)/xx)
          end do
          write(*,*) " ------------------------------------------ "
          write(*,*) " "
        end do
        DEALLOCATE(S0)

!! 2) CUTOFF EVALUATION !!
        xcutoff=100.0d0 !! SET HIGH FOR FIRST STEP !!
        iifrg=0
        iiorb=0
        do ifrg=1,icufr
          do ii=1,nocc

!! APPLYING CONDITIONS TO REMOVE ORBITALS !!
            xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
            if(xx.lt.xcutoff) then
              iiorb=ii
              iifrg=ifrg
              xcutoff=xx !! NEW LOWEST FOLI !!
            end if
          end do
        end do

!! NOW FRONTIER (SELECTION BY PACKS USING TOLERANCE) !!
        xfront=100.0d0 !! SET HIGH FOR FIRST STEP !!
        jjfrg=0
        jjorb=0
        do ifrg=1,icufr
          do ii=1,nocc
            xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
            if(xcutoff+folitol.lt.xx.and.xx.lt.xfront) then
              jjorb=ii
              jjfrg=ifrg
              xfront=xx
            end if
          end do
        end do

!! PRINTING !!
        write(*,'(2x,a27,x,i3,f10.5)') "Frg. and Lowest FOLI value:",iifrg,xcutoff
        write(*,'(2x,a56,x,i3,f10.5)') "Frg. and Lowest FOLI value including tolerance (cutoff):",jjfrg,xcutoff
        write(*,*) " "
        write(*,*) " ------------------- "
        write(*,*) "  SELECTED ORBITALS  "
        write(*,*) " ------------------- "
        write(*,*) " "

!! 3) EVALUATING DEGENERACIES !!
!! infopop(i,j): SAVING THE FRAGMENT IN i = 1 AND ORBITAL NUMBER IN i = 2 !!
!! j ALLOCATED AS nocc (FOR PRACTICITY) BUT MAXIMUM WILL BE inewcore !!
        ALLOCATE(infopop(2,nocc))
        inewcore=0

!! BRANCHING (CONTROLLED FROM .inp, DEFAULT = 0) !!
        if(iiter.eq.ibranch) then
          write(*,*) " ********************************************** "
          write(*,*) "  WARNING: BRANCHING INVOKED IN THIS ITERATION  "
          write(*,*) " ********************************************** "
          write(*,*) " "

!! MG: TO DO !!
          write(*,*) " Branching code has to be done "
          stop
!         inewcore=2
!         do iiii=1,inewcore
!           infopop(1,iiii)=1
!           if(iiii.eq.1) infopop(2,iiii)=11
!           if(iiii.eq.2) infopop(2,iiii)=12
!           scr(iiii)=dsqrt(deloc(1,infopop(2,iiii))/frgpop(1,infopop(2,iiii)))
!           write(*,'(a29,i3,a12,i3)') " Orbital Equal/Over Cutoff:",infopop(2,iiii),
!    +      " Fragment : ",infopop(1,iiii)
!         end do
        else
          write(*,*) "  Orb.   Frag.   FOLI  "
          write(*,*) " --------------------- "
          do ifrg=1,icufr
            do iorb=1,nocc
              xx=xcutoff-dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              if(ABS(xx).le.folitol) then

!! APPLYING CONDITIONS TO REMOVE ORBITALS !!
                inewcore=inewcore+1
                infopop(1,inewcore)=ifrg
                infopop(2,inewcore)=iorb
                scr(inewcore)=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
                write(*,112) infopop(2,inewcore),infopop(1,inewcore),scr(inewcore)
              end if
            end do
          end do
          write(*,*) " --------------------- "
          write(*,*) " "
        end if
        write(*,'(2x,a25,x,i3)') "Number of OSLOs selected:",inewcore
        write(*,'(2x,a17,x,f10.5)') "delta-FOLI value:",xfront-xcutoff
        write(*,*) " "

!! SAVING THE FRAG OSLOs (CONSIDERED CORE) IN ccore !!
        do iorb=1,inewcore
          iifrg=infopop(1,iorb)
          iiorb=infopop(2,iorb)
          ifrgel(iifrg)=ifrgel(iifrg)+1 !! ADDING THEM HERE !!
          iaddcore=iaddcore+1
          iaddoslo=iaddoslo+1

!! SAVING IN infooslo THE POPULATION AND FRAGMENT PREORTHOGONALIZATION !!
          poposlo(iaddoslo)=frgpop(iifrg,iiorb)
          infooslo(iaddoslo)=iifrg
          delocoslo(iaddoslo)=scr(iorb)
          fspread(iaddoslo)=frgspr(iifrg,iiorb)
          do mu=1,igr
            ccore(mu,iaddcore)=cfrgoslo(iifrg,mu,iiorb)
            coslo(mu,iaddoslo)=cfrgoslo(iifrg,mu,iiorb)
          end do
        end do
        ntotcore=ntotcore+iaddcore

!! EVALUATING OVERASSIGNMENT !!
        nnelect=0
        do ifrg=1,icufr
          nnelect=nnelect+ifrgel(ifrg)
        end do
        write(*,'(2x,a30,x,i3)') "Electron pairs left to assign:",nocc-nnelect
        write(*,*) " "

        if(nocc-nnelect.lt.0) then
          write(*,*) " *************************************** "
          write(*,'(3x,a34,x,i3)') "WARNING: OVERASSIGNING BY (pairs):",-(nocc-nnelect)
          write(*,*) " *************************************** "
          write(*,*) " "

!! DIRTY TRICK !!
          write(*,*) " Continues by tricking the code (overassigned pairs removed) "
          write(*,*) " Check the final OSs, overassigned electrons have to be afterwards " !! TO DO !!
          write(*,*) " "
          inewcore=inewcore+(nocc-nnelect)
          iaddoslo=iaddoslo+(nocc-nnelect)
          nnelect=nnelect+(nocc-nnelect)
        end if

!! SELECTING THE FIRST OUT FOR EVALUATING LINDEP !!
        clindep=ZERO
        write(*,*) " ------------------------------ "
        write(*,*) "  CHECKING LINEAR DEPENDENCIES  "
        write(*,*) " ------------------------------ "
        write(*,*) " "
        write(*,*) "  Orb.   Frag.   FOLI  "
        write(*,*) " --------------------- "
        iselected=0
        do ifrg=1,icufr
          do iorb=1,nocc
            xx=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))-xfront

!! CRITERIA FOR SELECTION: sqrt(deloc/Q_A) <= frontier + folitol !!
            if(xx.lt.folitol) then
              iselected=iselected+1
              xx2=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              write(*,112) iorb,ifrg,xx2
              do mu=1,igr
                clindep(mu,iselected)=cfrgoslo(ifrg,mu,iorb)
              end do
            end if
          end do
        end do
        write(*,*) " --------------------- "
        write(*,'(2x,a38,x,i3)') "Number of OSLOs for LinDep evaluation:",iselected
        write(*,*) " "

!! LAST ITERATION NO LINDEP EVALUATION !!
        iilindep=0
        if(nocc-nnelect.eq.0.and.inewcore.eq.1) then
          write(*,*) " LinDep not evaluated in last iteration if only 1 orbital is selected "
          write(*,*) " "
          iilindep=1
        end if

!! SOME DEALLOCATES... !!
        DEALLOCATE(deloc)
        DEALLOCATE(infopop)

!! EVALUATE LINDEP !!
        ilindep=0
        if(iilindep.eq.0) then !! IF ONLY 1 THERE IS NOTHING TO EVALUATE !!
          SSS=ZERO
          EEE=ZERO
          do ii=1,iselected
            do jj=1,iselected
              xx=ZERO
              do mu=1,igr
                do nu=1,igr
                  xx=xx+clindep(mu,ii)*clindep(nu,jj)*s(mu,nu)
                end do
              end do
              SSS(ii,jj)=xx
            end do
          end do

!! DIAGONALIZING ALL MATRIX, REST IS ZERO SO NO AFFECTS !!
          call diagonalize(nocc,nocc,SSS,EEE,0)

!! PRINTING SMALLEST EIGENVALUE !!
          xx=10.0d0 !! ABSURT VALUE AGAIN... !!
          do ii=1,iselected
            if(SSS(ii,ii).lt.xx) xx=SSS(ii,ii)
          end do
          if(xx.lt.1.0d-4) ilindep=1 !! THRESHOLD FOR LINIAR DEPENDENCY !!
          write(*,'(2x,a36,x,f10.5)') "Lowest eigenvalue obtained (LinDep):",xx
          write(*,*) " "
        end if
        if(ilindep.eq.1) then
          write(*,*) " *********************************** "
          write(*,*) "  WARNING : LINEAR DEPENDENCY FOUND  "
          write(*,*) " *********************************** "
          write(*,*) " "
          write(*,'(2x,a27,3f10.5)') "FOLI values and delta-FOLI:",xcutoff,xfront,xfront-xcutoff
          write(*,*) " Selecting largest to proceed "
          write(*,*) " RECOMMENDED TO BRANCH (.inp) AND CHECK ALTERNATIVE ASSIGNMENT "
          write(*,*) " "
        end if

!! REMOVING ORBITALS FROM P MATRIX, ONLY IF NOT ALL ARE ASSIGNED !!
!! ORTHOGONALIZING FRAGMENT CORE/SEMICORE ORBITALS !!
        ALLOCATE(S0(inewcore,inewcore),eigv(inewcore,inewcore))
        do ii=1,inewcore
          do jj=1,inewcore
            xx=ZERO
            do mu=1,igr
              do nu=1,igr
                xx=xx+ccore(mu,ii)*ccore(nu,jj)*s(mu,nu)
              end do
            end do
            S0(ii,jj)=xx
          end do
        end do
        
        call diagonalize(inewcore,inewcore,S0,eigv,0)
        
        ALLOCATE(smh(inewcore,inewcore))
        do ii=1,inewcore
          do jj=ii,inewcore
            smh(jj,ii)=ZERO
            do kk=1,inewcore
              if(S0(kk,kk).gt.thresh) then !! thresh = 10^-8 (see parameter.h) !!
                xx=eigv(ii,kk)*eigv(jj,kk)
                ssqrt=dsqrt(S0(kk,kk))
                smh(jj,ii)=smh(jj,ii)+xx/ssqrt
              end if
            end do
            smh(ii,jj)=smh(jj,ii)
          end do
        end do
        DEALLOCATE(eigv,S0)
 
        do ii=1,igr
          do jj=1,inewcore
            xx=ZERO
            do kk=1,inewcore
              xx=xx+smh(jj,kk)*ccore(ii,kk)
            end do
            ccoreorth(ii,jj)=xx
          end do
        end do

!! SAVING ORTHOGONAL ORBITALS HERE !!
        do ii=1,inewcore
          iaddoslo2=iaddoslo2+1
          do mu=1,igr
            cosloorth(mu,iaddoslo2)=ccoreorth(mu,ii)
          end do
        end do
        DEALLOCATE(smh)

!! CONSTRUCTING pcore (AND pnocore BY SUBSTRACTION) !!
        do ii=1,igr
          do jj=1,igr
            xx=ZERO
            do ij=1,inewcore
              xx=xx+ccoreorth(ii,ij)*ccoreorth(jj,ij)
            end do
            pcore(ii,jj)=xx
          end do
        end do
        do ii=1,igr
          do jj=1,igr
            pnocore(ii,jj)=pnocore(ii,jj)-pcore(ii,jj)
          end do
        end do

!! IN CASE OF ALL ASSIGNED !!
        if(nocc-nnelect.eq.0) go to 666

!! END OF ITERATIVE PROCEDURE !!
      end do
666   continue

!! FINAL OS ASSIGNMENT !!
      write(*,*) " --------------------------- "
      write(*,*) "  FRAGMENT OXIDATION STATES  "
      write(*,*) " --------------------------- "
      write(*,*) " "
      write(*,*) "  Frag.  Oxidation State  "
      write(*,*) " ------------------------ "
      do ifrg=1,icufr
        write(*,20) ifrg,REAL(iznfrg(ifrg)-2*ifrgel(ifrg))
      end do 
      write(*,*) " ------------------------ "
      write(*,*) " "

!! PRINTING OF THE .fchk FILES WITH THE OSLOs... TO VISUALIZE !!
!! PREORTHOGONALIZATION OSLOs CAN BE VISUALIZED IF DESIRED (.inp) !!
      ALLOCATE(poslo(igr,igr)) !! REQUIRED poslo FOR .fchk CREATION !!
      if(ifchk.eq.2) then
        do ii=1,igr
          do jj=1,igr
            xx=ZERO
            do ij=1,iaddoslo
              xx=xx+coslo(ii,ij)*coslo(jj,ij)
            end do
            poslo(ii,jj)=TWO*xx
          end do
        end do
        ctype="-OSLOs-preortho"
        call rwf_orbprint(coslo,poslo,ctype)
      end if

!! NOW THE FINAL (ORTHOGONALIZED) ONES !!
      do ii=1,igr
        do jj=1,igr
          xx=ZERO
          do ij=1,iaddoslo
            xx=xx+cosloorth(ii,ij)*cosloorth(jj,ij)
          end do
          poslo(ii,jj)=TWO*xx
        end do
      end do
      ctype="-OSLOs"
      call rwf_orbprint(cosloorth,poslo,ctype)

!! EVALUATING FINAL POPULATIONS TO COMPARE !!
      write(*,*) " ---------------------------------- "
      write(*,*) "  PRINTING FINAL OSLOs INFORMATION  "
      write(*,*) " ---------------------------------- "
      write(*,*) " "
      write(*,*) " ------------------------------------------- "
      write(*,*) "  Summary of the selected OSLOs (pre-ortho)  "
      write(*,*) " ------------------------------------------- "
      write(*,*) " "

!! MADE A BIT TRICKY... SORRY !!
      ALLOCATE(orbpop(nocc),orbpop2(nocc))
      ALLOCATE(foslo(nocc,icufr),foslo2(nocc,icufr))
      foslo=ZERO
      foslo2=ZERO
      do jfrg=1,icufr
        orbpop=ZERO
        orbpop2=ZERO
        call rwf_uwf_frg_pop(jfrg,sat,nocc,coslo,orbpop) !! FOR THE NON-ORTHOGONAL OSLOs (ORIGINAL) !!
        call rwf_uwf_frg_pop(jfrg,sat,nocc,cosloorth,orbpop2) !! FOR THE ORTHOGONALIZED ONES (PRINTING LATER) !!
        do ii=1,nocc
          foslo(ii,jfrg)=orbpop(ii)
          foslo2(ii,jfrg)=orbpop2(ii)
        end do
      end do
      call rwf_uwf_print_OSLO_final(1,nocc,delocoslo,foslo)
      write(*,*) " --------------------------------------- "
      write(*,*) "  Summary of the selected OSLOs (final)  "
      write(*,*) " --------------------------------------- "
      write(*,*) " "
      call rwf_uwf_print_OSLO_final(0,nocc,delocoslo,foslo2) !! FOLI VALUES GIVEN JUST FOR USING SAME ROUTINE !!
      
      DEALLOCATE(orbpop,orbpop2)
      DEALLOCATE(foslo,foslo2)

!! TODO: FINAL DEALLOCATE (IF CONSIDERED RELEVANT) !!

!! PRINTING FORMATS !!
110   FORMAT(3x,i3,3x,3f10.5)
111   FORMAT(3x,i3,2x,f10.5,3x,f10.5,3x,f10.5)
112   FORMAT(3x,i3,4x,i3,x,f9.5)
20    FORMAT(3x,i3,6x,f8.2)

      end

!! ****** !!

      subroutine rwf_uwf_frg_pop(ifrg,sat,norb,corb,frgpop)

!! THIS SUBROUTINE COMPUTES FRAGMENT POPULATION ANALYSIS !!
!! MADE INDEPENDENT OF THE WF-TYPE !!

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      dimension sat(igr,igr,nat)
      dimension corb(igr,igr),frgpop(norb)

      allocatable :: orbpop(:,:)

!! LOADING IOPTs !!
      idofr=iopt(40)

!! ALLOCATING MATRICES !!
      ALLOCATE(orbpop(norb,nat))

!! EVALUATING ORBITAL ATOMIC POPULATIONS !!
      do icenter=1,nat
        do iorb=1,norb
          xx=ZERO
          do jj=1,igr
            do kk=1,igr
              xx=xx+corb(kk,iorb)*sat(kk,jj,icenter)*corb(jj,iorb)
            end do
          end do

!! WITHOUT FACTOR OF 2 FOR BEING GENERAL, TO ADD OUTSIDE IF REQUIRED !!
          orbpop(iorb,icenter)=xx
        end do
      end do

!! GROUPING BY FRAGMENTS (ifrg) !!
      if(idofr.eq.1) then
        do iorb=1,norb
          xx=ZERO
          do icenter=1,nfrlist(ifrg)
            xx=xx+orbpop(iorb,ifrlist(icenter,ifrg))
          end do
          frgpop(iorb)=xx
        end do
      end if

!! DEALLOCATING MATRICES !!
      DEALLOCATE(orbpop)
  
      end 

!! ****** !!

      subroutine rwf_orbprint(cmat,pmat,ctype)

!! THIS SUBROUTINE PRINTS ORBITALS IN .fchk FORMAT !!
!! NOT ELEGANT WAY TO DO IT, BUT WORKS !!

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /filename/name0

      character*80 line
      character*60 name0,name1
      character*20 ctype

      dimension cmat(igr,igr),pmat(igr,igr)

      iqchem   = iopt(95)
      imokit   = iopt(79)
      indepigr = int_locate(15,"Number of independ",ilog)
      norb     = igr*indepigr
      norbt    = igr*(igr+1)/2

!! NAME OF THE .fchk FILE !!
      name1=trim(name0)//trim(ctype)//".fchk"
      open(unit=69,file=name1)
      rewind(69)
      rewind(15)

!! PRINTING UNTIL ALPHA MOs !!
      read(15,'(a80)') line
      do while(index(line,"Alpha MO co").eq.0)
        write(69,'(a80)') line
        read(15,'(a80)') line
      end do

!! PRINTING THE NEW ONES !!
      write(69,11) "Alpha MO coefficients","R","N= ",norb
      write(69,13) ((cmat(ii,jj),ii=1,igr),jj=1,indepigr)

!! NOW LOCATING WHAT IS AFTER IT IN THE ORIGINAL ONE TO CONTINUE !!
      if(iqchem.eq.0.and.imokit.eq.0) then
        do while(index(line,"Orthonormal basis").eq.0)
          read(15,'(a80)') line
        end do
      else if(iqchem.eq.1) then
        do while(index(line,"Alpha Orbital").eq.0)
          read(15,'(a80)') line
        end do
      end if

!! RESTART PRINTING UNTIL NEXT STOP !!
      do while(index(line,"Total SCF Dens").eq.0)
        write(69,'(a80)') line
        read(15,'(a80)') line
      end do

!! PRINTING THE NEW ONE !!
      write(69,12) "Total SCF Density","R","N= ",norbt
      write(69,13) ((pmat(ii,jj),jj=1,ii),ii=1,igr)

!! NOW LOCATING WHAT IS AFTER IT IN THE ORIGINAL ONE TO CONTINUE !!
      if(iqchem.eq.0.and.imokit.eq.0) then
        do while(index(line,"Mulliken Charges").eq.0)
          read(15,'(a80)') line
        end do
      else if(iqchem.eq.1) then
        do while(index(line,"Pure Switching").eq.0)
          read(15,'(a80)') line
        end do
      end if

!! RESTART PRINTING UNTIL THE END !!
      do while(.true.)
        write(69,'(a80)') line
        read(15,'(a80)',end=99) line
      end do
99    continue
      close(69)

!! PRINTING FORMATS !!
11    FORMAT(a21,22x,a1,3x,a3,i11)
12    FORMAT(a17,26x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))

      end

!! ****** !!

      subroutine rwf_uwf_print_OSLO_final(iflag,noslo,foli,frgpop)

!! ROUTINE FOR PRINTING THE FRAGMENT POPULATIONS AND FOLI FOR EACH OSLO!!

      implicit double precision (a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      dimension foli(igr)
      dimension frgpop(noslo,icufr)

!! TRICK OF THE INTEGER ROUNDING FOR NUMBER OF COLUMNS !!
      b=noslo/5
      if(b*5.ne.noslo) b=(noslo/5)+1

!! PRINTING !!
      dd=1
      do k=1,b

!! FOR THE LAST PACK OF COLUMNS !!
        if(k.eq.b) then
          write(*,'(2x,a13,5(i7,3x))') "OSLO Number :",(jj,jj=dd,noslo)
          if(iflag.eq.1) write(*,'(2x,a13,5f10.5)') "FOLI Value  :",(foli(jj),jj=dd,noslo)
          do jfrg=1,icufr
            write(*,'(2x,a9,i3,a1,5f10.5)') "Frg. Pop.",jfrg,":",(frgpop(jj,jfrg),jj=dd,noslo)
          end do
          write(*,*) " "

!! FOR PACKS OF 5 COLUMNS !!
        else
          write(*,'(2x,a13,5(i7,3x))') "OSLO Number :",(jj,jj=dd,dd+4)
          if(iflag.eq.1) write(*,'(2x,a13,5f10.5)') "FOLI Value  :",(foli(jj),jj=dd,dd+4)
          do jfrg=1,icufr
            write(*,'(2x,a9,i3,a1,5f10.5)') "Frg. Pop.",jfrg,":",(frgpop(jj,jfrg),jj=dd,dd+4)
          end do
          write(*,*) " "
          dd=dd+5
        end if
      end do

      end

!! ****** !!

!! NOW UNRESTRICTED SUBROUTINES !!

!! ****** !!
      
      subroutine uwf_iterative_oslo(sat,itotps,wp,omp2,chp,pcoord)

!! MG: ONLY ITERATIVE OSLO IMPLEMENTED IN THIS VERSION !!
!! MG: NON-ITERATIVE PROCEDURE IMPLEMENTED IN DEVELOPMENT VERSION !!

      use basis_set
      use ao_matrices
      use integration_grid
      implicit double precision(a-h,o-z)
      include 'parameter.h'
  
      integer,intent(in) :: itotps
  
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /filename/name0
  
      character*80 line
      character*60 name0,name1
      character*20 ctype
  
      dimension sat(igr,igr,nat)
      dimension chp(itotps,igr),pcoord(itotps,3)
      dimension wp(itotps),omp2(itotps,nat)
  
      allocatable :: Smat(:,:,:)
      allocatable :: S0(:,:),Sm(:,:),Splus(:,:),smh(:,:),eigv(:,:)
      allocatable :: c0(:,:),pp0(:,:)
  
      allocatable :: SSS(:,:),EEE(:,:)
  
      allocatable :: cmat(:,:),cfrgoslo(:,:,:)
      allocatable :: orbpop(:),orbpop2(:),frgpop(:,:),frgspr(:,:)
      allocatable :: foslo(:,:),foslo2(:,:),foli(:)
  
      allocatable :: infopop(:,:),infooslo(:,:),poposlo(:,:)
      allocatable :: fspread(:),fbspread(:),scr(:),scrb(:),deloc(:,:),delocoslo(:,:)
  
      allocatable :: clindep(:,:)
  
      allocatable :: pcore(:,:),pnocore(:,:)
      allocatable :: ccore(:,:),ccoreorth(:,:)
  
      allocatable :: coslo(:,:),cosloorth(:,:)
      allocatable :: cboslo(:,:),cbosloorth(:,:)
      allocatable :: poslo(:,:)
  
      allocatable :: ifrgel(:,:),iznfrg(:)
  
!! LOADING IOPTs !!
      iqchem   = iopt(95)
      ifolitol = iopt(96)
      ibranch  = iopt(97)
      ifchk    = iopt(98)

!! PARAMETERS REQUIRED !!
      iatps = nrad*nang
      niter = 999 !! MAX NUMBER OF ITERATIONS !!
      folitol = 10.0d0**(-REAL(ifolitol)) !! DEFAULT = 10^-3, CONTROLLED IN .inp !!

!! ALLOCATING GENERAL (REQUIRED) MATRICES !!
      ALLOCATE(Sm(igr,igr),Splus(igr,igr))
      ALLOCATE(c0(igr,igr),pp0(igr,igr))
      ALLOCATE(ccore(igr,igr),pcore(igr,igr),ccoreorth(igr,igr))
      ALLOCATE(coslo(igr,igr),cosloorth(igr,igr))
      ALLOCATE(cboslo(igr,igr),cbosloorth(igr,igr),delocoslo(2,igr))
      ALLOCATE(clindep(igr,igr))
      ALLOCATE(fspread(nalf),fbspread(nb),scr(nalf),scrb(nb))
      ALLOCATE(infooslo(2,igr),poposlo(2,igr))
      ALLOCATE(ifrgel(2,icufr),iznfrg(icufr))

!! FRAGMENT CHARGE EXTRACTED FROM zn (AVOIDS PROBLEMS WHEN PSEUDOPOTENTIALS ARE USED) !!
      do ifrg=1,icufr
        izn=0
        do icenter=1,nfrlist(ifrg)
          iiat=ifrlist(icenter,ifrg)
          izn=izn+INT(zn(iiat))
        end do
        iznfrg(ifrg)=izn
        ifrgel(1,ifrg)=0 !! ALPHA IN 1 !!
        ifrgel(2,ifrg)=0 !! BETA IN 2 !!
      end do

!! FRAGMENT LOCALIZATION: IMPORTANT :: R_A VALUE SELECTED IS THE CENTER OF CHARGE BETWEEN FRAGMENT ATOMS !!
      write(*,*) " --------------------------------------- "
      write(*,*) "  CHARGE CENTER (R_F) FOR EACH FRAGMENT  "
      write(*,*) " --------------------------------------- "
      write(*,*) " "
      write(*,*) "  Frag.        Charge center (xyz)     "
      write(*,*) " ------------------------------------- "
      ALLOCATE(Smat(ifrg,igr,igr))
      do ifrg=1,icufr
        xcenter=ZERO
        ycenter=ZERO
        zcenter=ZERO
        xchg=ZERO
        do icenter=1,nfrlist(ifrg)
          iiat=ifrlist(icenter,ifrg)
          xatchg=REAL(iznuc(iiat))
          xcenter=xcenter+xatchg*coord(1,iiat)
          ycenter=ycenter+xatchg*coord(2,iiat)
          zcenter=zcenter+xatchg*coord(3,iiat)
          xchg=xchg+xatchg
        end do
        xcenter=xcenter/xchg
        ycenter=ycenter/xchg
        zcenter=zcenter/xchg
        write(*,110) ifrg,xcenter*angtoau,ycenter*angtoau,zcenter*angtoau

!! COMPUTING S ONLY ONCE !!
        do mu=1,igr
          do nu=1,mu
            xx=ZERO
            do jcenter=1,nat
              do ifut=iatps*(jcenter-1)+1,iatps*jcenter
                distx=pcoord(ifut,1)-xcenter
                disty=pcoord(ifut,2)-ycenter
                distz=pcoord(ifut,3)-zcenter
                xoper=(distx*distx+disty*disty+distz*distz)
                xx=xx+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*omp2(ifut,jcenter)*xoper
              end do
            end do
            Smat(ifrg,mu,nu)=xx
            Smat(ifrg,nu,mu)=xx
          end do
        end do
      end do
      write(*,*) " ------------------------------------- "

!! INITIAL PRINTING !!
      write(*,*) " "
      write(*,*) " ----------------------------------- "
      write(*,*) "  STARTING ITERATIVE OSLO ALGORITHM  "
      write(*,*) " ----------------------------------- "
      write(*,*) " "
      write(*,'(2x,a50,f10.5)') "Tolerance (in delta-FOLI) used for OSLO selection:",folitol
      write(*,*) " "

!! PREPARING FOR ITERATIVE PROCESS !!
      ALLOCATE(cmat(igr,igr),cfrgoslo(icufr,igr,igr))
      ALLOCATE(SSS(nalf,nalf),EEE(nalf,nalf)) !! ALLOCATING HERE FOR ALPHA PART !!
      ALLOCATE(frgpop(icufr,nalf))
      ALLOCATE(frgspr(icufr,nalf))
      ALLOCATE(pnocore(igr,igr))

!! ZEROING MATRICES LIKE THIS NOW !!
      coslo=ZERO
      cosloorth=ZERO
      pnocore=pa !! SAVING pa IN pnocore !!

!! FIRST ALPHA, THEN BETA !!
      write(*,*) " ------------ "
      write(*,*) "  ALPHA PART  "
      write(*,*) " ------------ "
      write(*,*) " "

      ntotcore=0
      iaddoslo=0
      iaddoslo2=0
      do iiter=1,niter
        write(*,*) " ---------------------- "
        write(*,'(3x,a16,x,i3)') "ITERATION NUMBER",iiter
        write(*,*) " ---------------------- "
        write(*,*) " "

!! ZEROING THE INVOLVED MATRICES !!
        pcore=ZERO
        ccore=ZERO
        ccoreorth=ZERO
        ALLOCATE(deloc(icufr,nalf))
        deloc=ZERO

!! 1) OBTAINING OSLOs FOR ALL FRGS !!
        iaddcore=0
        ALLOCATE(S0(igr,igr))
        do ifrg=1,icufr
          do ii=1,igr
            do jj=1,igr

!! USING pp0 AND S0 TO NOT DESTROY pnocore AND Smat !!
              pp0(ii,jj)=pnocore(ii,jj)
              S0(ii,jj)=Smat(ifrg,ii,jj)
            end do
          end do
          call build_Smp(igr,S0,Sm,Splus,0)
          call to_lowdin_basis(igr,Splus,pp0)
          call diagonalize(igr,igr,pp0,c0,0)
          call to_AO_basis(igr,igr,Sm,c0)
          imaxo=nalf !! MAX NUMBER OF OSLOs IS NOW nalf FOR PRACTICITY !!

!! RECOVERING THE COEFFS, SAVED IN cfrgoslo !!
          do kk=1,imaxo
            do mu=1,igr
              cmat(mu,kk)=c0(mu,kk)*dsqrt(pp0(kk,kk))
              cfrgoslo(ifrg,mu,kk)=cmat(mu,kk)
            end do
          end do

!! COMPUTING PIPEK DELOCALIZATION, REQUIRES FRAGMENT POPULATIONS !!
          ALLOCATE(orbpop(nalf))
          do jfrg=1,icufr
            call rwf_uwf_frg_pop(jfrg,sat,nalf,cmat,orbpop)
            do ii=1,imaxo
              deloc(ifrg,ii)=deloc(ifrg,ii)+orbpop(ii)*orbpop(ii) !! ADDING Q_A**2 INSIDE deloc !!

!! IMPORTANT HERE SAVING ONLY FOR THE OWN FRAGMENT !!
              if(jfrg.eq.ifrg) then
                frgspr(ifrg,ii)=pp0(ii,ii) !! SPREADS (ONLY FOR PRINTING) !!
                frgpop(ifrg,ii)=orbpop(ii) !! FRAGMENT POPULATIONS !!
              end if
            end do
          end do
          DEALLOCATE(orbpop)

!! NOW DOING 1/deloc() !!
          do ii=1,imaxo
            if(deloc(ifrg,ii).gt.1.0d-6) then
              deloc(ifrg,ii)=ONE/deloc(ifrg,ii)
            else
              deloc(ifrg,ii)=100.0d0 !! ABSURT VALUE, AVOIDS PROBLEMS !!
            end if
          end do

!! PRINTING VALUABLE INFORMATION !!
          write(*,*) " -------------------------------------- "
          write(*,'(3x,a32,x,i3)') "ORBITAL INFORMATION FOR FRAGMENT",ifrg
          write(*,*) " -------------------------------------- "
          write(*,*) " "
          write(*,*) "  Orb.    Spread      Frg. Pop.      FOLI   "
          write(*,*) " ------------------------------------------ "
          do ii=1,imaxo
            xx=frgpop(ifrg,ii)
            if(xx.gt.1.0d-5) write(*,111) ii,frgspr(ifrg,ii),xx,dsqrt(deloc(ifrg,ii)/xx)
          end do
          write(*,*) " ------------------------------------------ "
          write(*,*) " "
        end do
        DEALLOCATE(S0)

!! 2) CUTOFF EVALUATION !!
        xcutoff=100.0d0 !! SET HIGH FOR FIRST STEP !!
        iifrg=0
        iiorb=0
        do ifrg=1,icufr
          do ii=1,nalf

!! APPLYING CONDITIONS TO REMOVE ORBITALS !!
            xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
            if(xx.lt.xcutoff) then
              iiorb=ii
              iifrg=ifrg
              xcutoff=xx !! NEW LOWEST FOLI !!
            end if
          end do
        end do

!! NOW FRONTIER (SELECTION BY PACKS USING TOLERANCE) !!
        xfront=100.0d0 !! SET HIGH FOR FIRST STEP !!
        jjfrg=0
        jjorb=0
        do ifrg=1,icufr
          do ii=1,nalf
            xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
            if(xcutoff+folitol.lt.xx.and.xx.lt.xfront) then
              jjorb=ii
              jjfrg=ifrg
              xfront=xx
            end if
          end do
        end do

!! PRINTING !!
        write(*,'(2x,a27,x,i3,f10.5)') "Frg. and Lowest FOLI value:",iifrg,xcutoff
        write(*,'(2x,a56,x,i3,f10.5)') "Frg. and Lowest FOLI value including tolerance (cutoff):",jjfrg,xcutoff
        write(*,*) " "
        write(*,*) " ------------------- "
        write(*,*) "  SELECTED ORBITALS  "
        write(*,*) " ------------------- "
        write(*,*) " "

!! 3) EVALUATING DEGENERACIES !!
!! infopop(i,j): SAVING THE FRAGMENT IN i = 1 AND ORBITAL NUMBER IN i = 2 !!
!! j ALLOCATED AS nalf (FOR PRACTICITY) BUT MAXIMUM WILL BE inewcore !!
        ALLOCATE(infopop(2,nalf))
        inewcore=0

!! BRANCHING (CONTROLLED FROM .inp, DEFAULT = 0) !!
        if(iiter.eq.ibranch) then
          write(*,*) " ********************************************** "
          write(*,*) "  WARNING: BRANCHING INVOKED IN THIS ITERATION  "
          write(*,*) " ********************************************** "
          write(*,*) " "

!! MG: TO DO !!
          write(*,*) " Branching code has to be done "
          stop
!         inewcore=2
!         do iiii=1,inewcore
!           infopop(1,iiii)=1
!           if(iiii.eq.1) infopop(2,iiii)=11
!           if(iiii.eq.2) infopop(2,iiii)=12
!           scr(iiii)=dsqrt(deloc(1,infopop(2,iiii))/frgpop(1,infopop(2,iiii)))
!           write(*,'(a29,i3,a12,i3)') " Orbital Equal/Over Cutoff:",infopop(2,iiii),
!    +      " Fragment : ",infopop(1,iiii)
!         end do
        else
          write(*,*) "  Orb.   Frag.   FOLI  "
          write(*,*) " --------------------- "
          do ifrg=1,icufr
            do iorb=1,nalf
              xx=xcutoff-dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              if(ABS(xx).le.folitol) then

!! APPLYING CONDITIONS TO REMOVE ORBITALS !!
                inewcore=inewcore+1
                infopop(1,inewcore)=ifrg
                infopop(2,inewcore)=iorb
                scr(inewcore)=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
                write(*,112) infopop(2,inewcore),infopop(1,inewcore),scr(inewcore)
              end if
            end do
          end do
          write(*,*) " --------------------- "
        end if
        write(*,'(2x,a25,x,i3)') "Number of OSLOs selected:",inewcore
        write(*,'(2x,a17,x,f10.5)') "delta-FOLI value:",xfront-xcutoff
        write(*,*) " "

!! SAVING THE FRAG OSLOs (CONSIDERED CORE) IN ccore !!
        do iorb=1,inewcore
          iifrg=infopop(1,iorb)
          iiorb=infopop(2,iorb)
          ifrgel(1,iifrg)=ifrgel(1,iifrg)+1 !! ADDING THEM HERE !!
          iaddcore=iaddcore+1
          iaddoslo=iaddoslo+1

!! SAVING IN infooslo THE POPULATION AND FRAGMENT PREORTHOGONALIZATION !!
          poposlo(1,iaddoslo)=frgpop(iifrg,iiorb)
          infooslo(1,iaddoslo)=iifrg
          delocoslo(1,iaddoslo)=scr(iorb)
          fspread(iaddoslo)=frgspr(iifrg,iiorb)
          do mu=1,igr
            ccore(mu,iaddcore)=cfrgoslo(iifrg,mu,iiorb)
            coslo(mu,iaddoslo)=cfrgoslo(iifrg,mu,iiorb)
          end do
        end do
        ntotcore=ntotcore+iaddcore

!! EVALUATING OVERASSIGNMENT !!
        nnelect=0
        do ifrg=1,icufr
          nnelect=nnelect+ifrgel(1,ifrg)
        end do
        write(*,'(2x,a31,x,i3)') "Alpha electrons left to assign:",nalf-nnelect
        write(*,*) " "

        if(nalf-nnelect.lt.0) then
          write(*,*) " ******************************** "
          write(*,'(3x,a26,x,i3)') "WARNING: OVERASSIGNING BY:",-(nalf-nnelect)
          write(*,*) " ******************************** "
          write(*,*) " "

!! DIRTY TRICK !!
          write(*,*) " Continues by tricking the code (overassigned electrons removed) "
          write(*,*) " Check the final OSs, overassigned electrons have to be afterwards " !! TO DO !!
          write(*,*) " "
          inewcore=inewcore+(nalf-nnelect)
          iaddoslo=iaddoslo+(nalf-nnelect)
          nnelect=nnelect+(nalf-nnelect)
        end if

!! SELECTING THE FIRST OUT FOR EVALUATING LINDEP !!
        clindep=ZERO
        write(*,*) " ------------------------------ "
        write(*,*) "  CHECKING LINIAR DEPENDENCIES  "
        write(*,*) " ------------------------------ "
        write(*,*) " "
        write(*,*) "  Orb.   Frag.   FOLI  "
        write(*,*) " --------------------- "
        iselected=0
        do ifrg=1,icufr
          do iorb=1,nalf
            xx=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))-xfront

!! CRITERIA FOR SELECTION: sqrt(deloc/Q_A) <= frontier + folitol !!
            if(xx.lt.folitol) then
              iselected=iselected+1
              xx2=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              write(*,112) iorb,ifrg,xx2
              do mu=1,igr
                clindep(mu,iselected)=cfrgoslo(ifrg,mu,iorb)
              end do
            end if
          end do
        end do
        write(*,*) " --------------------- "
        write(*,'(2x,a38,x,i3)') "Number of OSLOs for LinDep evaluation:",iselected
        write(*,*) " "

!! LAST ITERATION NO LINDEP EVALUATION !!
        iilindep=0
        if(nalf-nnelect.eq.0.and.inewcore.eq.1) then
          write(*,*) " LinDep not evaluated in last iteration if only 1 orbital is selected "
          write(*,*) " "
          iilindep=1
        end if

!! SOME DEALLOCATES... !!
        DEALLOCATE(deloc)
        DEALLOCATE(infopop)

!! EVALUATE LINDEP !!
        ilindep=0
        if(iilindep.eq.0) then
          SSS=ZERO
          EEE=ZERO
          do ii=1,iselected
            do jj=1,iselected
              xx=ZERO
              do mu=1,igr
                do nu=1,igr
                  xx=xx+clindep(mu,ii)*clindep(nu,jj)*s(mu,nu)
                end do
              end do
              SSS(ii,jj)=xx
            end do
          end do

!! DIAGONALIZING ALL MATRIX, REST IS ZERO SO NO AFFECTS !!
          call diagonalize(nalf,nalf,SSS,EEE,0)

!! PRINTING SMALLEST EIGENVALUE !!
          xx=10.0d0 !! ABSURT VALUE AGAIN... !!
          do ii=1,iselected
            if(SSS(ii,ii).lt.xx) xx=SSS(ii,ii)
          end do
          if(xx.lt.1.0d-4) ilindep=1 !! THRESHOLD FOR LINIAR DEPENDENCY !!
          write(*,'(2x,a36,x,f10.5)') "Lowest eigenvalue obtained (LinDep):",xx
          write(*,*) " "
        end if
        if(ilindep.eq.1) then
          write(*,*) " *********************************** "
          write(*,*) "  WARNING : LINIAR DEPENDENCY FOUND "
          write(*,*) " *********************************** "
          write(*,*) " "
          write(*,'(2x,a27,3f10.5)') "FOLI values and delta-FOLI:",xcutoff,xfront,xfront-xcutoff
          write(*,*) " Selecting largest to proceed "
          write(*,*) " RECOMMENDED TO BRANCH (.inp) AND CHECK ALTERNATIVE ASSIGNMENT "
          write(*,*) " "
        end if

!! REMOVING ORBITALS FROM P MATRIX, ONLY IF NOT ALL ARE ASSIGNED !!
!! ORTHOGONALIZING FRAGMENT CORE/SEMICORE ORBITALS !!
        ALLOCATE(S0(inewcore,inewcore),eigv(inewcore,inewcore))
        do ii=1,inewcore
          do jj=1,inewcore
            xx=ZERO
            do mu=1,igr
              do nu=1,igr
                xx=xx+ccore(mu,ii)*ccore(nu,jj)*s(mu,nu)
              end do
            end do
            S0(ii,jj)=xx
          end do
        end do
        
        call diagonalize(inewcore,inewcore,S0,eigv,0)
        
        ALLOCATE(smh(inewcore,inewcore))
        do ii=1,inewcore
          do jj=ii,inewcore
            smh(jj,ii)=ZERO
            do kk=1,inewcore
              if(S0(kk,kk).gt.thresh) then
                xx=eigv(ii,kk)*eigv(jj,kk)
                ssqrt=dsqrt(S0(kk,kk))
                smh(jj,ii)=smh(jj,ii)+xx/ssqrt
              end if
            end do
            smh(ii,jj)=smh(jj,ii)
          end do
        end do
        DEALLOCATE(eigv,S0)

        do ii=1,igr
          do jj=1,inewcore
            xx=ZERO
            do kk=1,inewcore
              xx=xx+smh(jj,kk)*ccore(ii,kk)
            end do
            ccoreorth(ii,jj)=xx
          end do
        end do

!! SAVING ORTHOGONAL ORBITALS HERE !!
        do ii=1,inewcore
          iaddoslo2=iaddoslo2+1
          do mu=1,igr
            cosloorth(mu,iaddoslo2)=ccoreorth(mu,ii)
          end do
        end do
        DEALLOCATE(smh)

!! CONSTRUCTING pcore (AND pnocore BY SUBSTRACTION) !!
        do ii=1,igr
          do jj=1,igr
            xx=ZERO
            do ij=1,inewcore
              xx=xx+ccoreorth(ii,ij)*ccoreorth(jj,ij)
            end do
            pcore(ii,jj)=xx
          end do
        end do
        do ii=1,igr
          do jj=1,igr
            pnocore(ii,jj)=pnocore(ii,jj)-pcore(ii,jj)
          end do
        end do

!! IN CASE OF ALL ALPHA ELECTRONS ASSIGNED !!
        if(nalf-nnelect.eq.0) go to 666

!! END OF ALPHA ITERATIVE PROCEDURE !!
      end do
666   continue

!! PRINTING OF THE .fchk FILES WITH THE OSLOs... TO VISUALIZE !!
!! PREORTHOGONALIZATION OSLOs CAN BE VISUALIZED IF DESIRED (.inp) !!
      ALLOCATE(poslo(igr,igr)) !! REQUIRED poslo FOR .fchk CREATION !!
      if(ifchk.eq.2) then
        do ii=1,igr
          do jj=1,igr
            xx=ZERO
            do ij=1,iaddoslo
              xx=xx+coslo(ii,ij)*coslo(jj,ij)
            end do
            poslo(ii,jj)=xx
          end do
        end do
        ctype="-alpha-OSLOs-preortho"
        call uwf_orbprint(0,coslo,poslo,ctype)
      end if

!! NOW THE FINAL (ORTHOGONALIZED) ONES !!
      do ii=1,igr
        do jj=1,igr
          xx=ZERO
          do ij=1,iaddoslo
            xx=xx+cosloorth(ii,ij)*cosloorth(jj,ij)
          end do
          poslo(ii,jj)=xx
        end do
      end do
      ctype="-alpha-OSLOs"
      call uwf_orbprint(0,cosloorth,poslo,ctype)

!! DEALLOCATES FOR STARTING BETA !!
      DEALLOCATE(poslo)
      DEALLOCATE(frgpop,frgspr)
      DEALLOCATE(SSS,EEE)

!! NOW BETA !!
      write(*,*) " ----------- "
      write(*,*) "  BETA PART  "
      write(*,*) " ----------- "
      write(*,*) " "

!! PREPARING SOME MATRICES !!
      cboslo=ZERO
      cbosloorth=ZERO
      pnocore=pb

!! REALLOCATES... !!
      ALLOCATE(frgpop(icufr,nb))
      ALLOCATE(frgspr(icufr,nb))
      ALLOCATE(SSS(nb,nb),EEE(nb,nb))

      ntotcore=0
      iaddoslo=0
      iaddoslo2=0
      do iiter=1,niter
        write(*,*) " ---------------------- "
        write(*,'(3x,a16,x,i3)') "ITERATION NUMBER",iiter
        write(*,*) " ---------------------- "
        write(*,*) " "

!! ZEROING THE INVOLVED MATRICES !!
        pcore=ZERO
        ccore=ZERO
        ccoreorth=ZERO
        ALLOCATE(deloc(icufr,nb))
        deloc=ZERO

!! 1) OBTAINING OSLOs FOR ALL FRGS !!
        iaddcore=0
        ALLOCATE(S0(igr,igr))
        do ifrg=1,icufr
          do ii=1,igr
            do jj=1,igr

!! USING pp0 AND S0 TO NOT DESTROY pnocore AND Smat !!
              pp0(ii,jj)=pnocore(ii,jj)
              S0(ii,jj)=Smat(ifrg,ii,jj)
            end do
          end do
          call build_Smp(igr,S0,Sm,Splus,0)
          call to_lowdin_basis(igr,Splus,pp0)
          call diagonalize(igr,igr,pp0,c0,0)
          call to_AO_basis(igr,igr,Sm,c0)
          imaxo=nb !! MAX NUMBER OF OSLOs IS NOW nb FOR PRACTICITY !!

!! RECOVERING THE COEFFS, SAVED IN cfrgoslo !!
          do kk=1,imaxo
            do mu=1,igr
              cmat(mu,kk)=c0(mu,kk)*sqrt(pp0(kk,kk))
              cfrgoslo(ifrg,mu,kk)=cmat(mu,kk)
            end do
          end do

!! COMPUTING PIPEK DELOCALIZATION, REQUIRES FRAGMENT POPULATIONS !!
          ALLOCATE(orbpop(nb))
          do jfrg=1,icufr
            call rwf_uwf_frg_pop(jfrg,sat,nb,cmat,orbpop)
            do ii=1,imaxo
              deloc(ifrg,ii)=deloc(ifrg,ii)+orbpop(ii)*orbpop(ii) !! ADDING Q_A**2 INSIDE deloc !!
  
!! IMPORTANT HERE SAVING ONLY FOR THE OWN FRAGMENT !!
              if(jfrg.eq.ifrg) then
                frgspr(ifrg,ii)=pp0(ii,ii) !! SPREADS (ONLY FOR PRINTING) !!
                frgpop(ifrg,ii)=orbpop(ii) !! FRAGMENT POPULATIONS !!
              end if
            end do
          end do
          DEALLOCATE(orbpop)

!! NOW DOING 1/deloc() !!
          do ii=1,imaxo
            if(deloc(ifrg,ii).gt.1.0d-6) then
              deloc(ifrg,ii)=ONE/deloc(ifrg,ii)
            else
              deloc(ifrg,ii)=100.0d0 !! ABSURT VALUE, AVOIDS PROBLEMS !!
            end if
          end do

!! PRINTING VALUABLE INFORMATION !!
          write(*,*) " -------------------------------------- "
          write(*,'(3x,a32,x,i3)') "ORBITAL INFORMATION FOR FRAGMENT",ifrg
          write(*,*) " -------------------------------------- "
          write(*,*) " "
          write(*,*) "  Orb.    Spread      Frg. Pop.      FOLI   "
          write(*,*) " ------------------------------------------ "
          do ii=1,imaxo
            xx=frgpop(ifrg,ii)
            if(xx.gt.1.0d-5) write(*,111) ii,frgspr(ifrg,ii),xx,dsqrt(deloc(ifrg,ii)/xx)
          end do
          write(*,*) " ------------------------------------------ "
          write(*,*) " "
        end do
        DEALLOCATE(S0)

!! 2) CUTOFF EVALUATION !!
        xcutoff=100.0d0 !! SET HIGH FOR FIRST STEP !!
        iifrg=0
        iiorb=0
        do ifrg=1,icufr
          do ii=1,nb

!! APPLYING CONDITIONS TO REMOVE ORBITALS !!
            xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
            if(xx.lt.xcutoff) then
              iiorb=ii
              iifrg=ifrg
              xcutoff=xx !! NEW LOWEST FOLI !!
            end if
          end do
        end do

!! NOW FRONTIER (SELECTION BY PACKS USING TOLERANCE) !!
        xfront=100.0d0 !! SET HIGH FOR FIRST STEP !!
        jjfrg=0
        jjorb=0
        do ifrg=1,icufr
          do ii=1,nb
            xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
            if(xcutoff+folitol.lt.xx.and.xx.lt.xfront) then
              jjorb=ii
              jjfrg=ifrg
              xfront=xx
            end if
          end do
        end do

!! PRINTING !!
        write(*,'(2x,a27,x,i3,f10.5)') "Frg. and Lowest FOLI value:",iifrg,xcutoff
        write(*,'(2x,a56,x,i3,f10.5)') "Frg. and Lowest FOLI value including tolerance (cutoff):",jjfrg,xcutoff
        write(*,*) " "
        write(*,*) " ------------------- "
        write(*,*) "  SELECTED ORBITALS  "
        write(*,*) " ------------------- "
        write(*,*) " "

!! 3) EVALUATING DEGENERACIES !!
!! infopop(i,j): SAVING THE FRAGMENT IN i = 1 AND ORBITAL NUMBER IN i = 2 !!
!! j ALLOCATED AS nb (FOR PRACTICITY) BUT MAXIMUM WILL BE inewcore !!
        ALLOCATE(infopop(2,nb))
        inewcore=0

        !! BRANCHING (CONTROLLED FROM .inp, DEFAULT = 0) !!
        if(iiter.eq.ibranch) then
          write(*,*) " ********************************************** "
          write(*,*) "  WARNING: BRANCHING INVOKED IN THIS ITERATION  "
          write(*,*) " ********************************************** "
          write(*,*) " "

!! MG: TO DO !!
          write(*,*) " Branching code has to be done "
          stop
!         inewcore=2
!         do iiii=1,inewcore
!           infopop(1,iiii)=1
!           if(iiii.eq.1) infopop(2,iiii)=11
!           if(iiii.eq.2) infopop(2,iiii)=12
!           scr(iiii)=dsqrt(deloc(1,infopop(2,iiii))/frgpop(1,infopop(2,iiii)))
!           write(*,'(a29,i3,a12,i3)') " Orbital Equal/Over Cutoff:",infopop(2,iiii),
!    +      " Fragment : ",infopop(1,iiii)
!         end do
        else
          write(*,*) "  Orb.   Frag.   FOLI  "
          write(*,*) " --------------------- "
          do ifrg=1,icufr
            do iorb=1,nb
              xx=xcutoff-dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              if(ABS(xx).le.folitol) then

!! APPLYING CONDITIONS TO REMOVE ORBITALS !!

                inewcore=inewcore+1
                infopop(1,inewcore)=ifrg
                infopop(2,inewcore)=iorb
                scrb(inewcore)=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
                write(*,112) infopop(2,inewcore),infopop(1,inewcore),scr(inewcore)
              end if
            end do
          end do
          write(*,*) " --------------------- "
        end if
        write(*,'(2x,a25,x,i3)') "Number of OSLOs selected:",inewcore
        write(*,'(2x,a17,x,f10.5)') "delta-FOLI value:",xfront-xcutoff
        write(*,*) " "

!! SAVING THE FRAG OSLOs (CONSIDERED CORE) IN ccore !!
        do iorb=1,inewcore
          iifrg=infopop(1,iorb)
          iiorb=infopop(2,iorb)
          ifrgel(2,iifrg)=ifrgel(2,iifrg)+1 !! ADDING THEM HERE !!
          iaddcore=iaddcore+1
          iaddoslo=iaddoslo+1

!! SAVING IN infooslo THE POPULATION AND FRAGMENT PREORTHOGONALIZATION !!
          poposlo(2,iaddoslo)=frgpop(iifrg,iiorb)
          infooslo(2,iaddoslo)=iifrg
          delocoslo(2,iaddoslo)=scrb(iorb)
          fbspread(iaddoslo)=frgspr(iifrg,iiorb)
          do mu=1,igr
            ccore(mu,iaddcore)=cfrgoslo(iifrg,mu,iiorb)
            cboslo(mu,iaddoslo)=cfrgoslo(iifrg,mu,iiorb)
          end do
        end do
        ntotcore=ntotcore+iaddcore

!! EVALUATING OVERASSIGNMENT !!
        nnelect=0
        do ifrg=1,icufr
          nnelect=nnelect+ifrgel(2,ifrg)
        end do
        write(*,'(2x,a30,x,i3)') "Beta electrons left to assign:",nb-nnelect
        write(*,*) " "

        if(nb-nnelect.lt.0) then
          write(*,*) " ******************************** "
          write(*,'(3x,a26,x,i3)') "WARNING: OVERASSIGNING BY:",-(nb-nnelect)
          write(*,*) " ******************************** "
          write(*,*) " "

!! DIRTY TRICK !!
          write(*,*) " Continues by tricking the code (overassigned electrons removed) "
          write(*,*) " Check the final OSs, overassigned electrons have to be afterwards " !! TO DO!!
          write(*,*) " "
          inewcore=inewcore+(nb-nnelect)
          iaddoslo=iaddoslo+(nb-nnelect)
          nnelect=nnelect+(nb-nnelect)
        end if

!! SELECTING THE FIRST OUT FOR EVALUATING LINDEP !!
        clindep=ZERO
        write(*,*) " ------------------------------ "
        write(*,*) "  CHECKING LINIAR DEPENDENCIES  "
        write(*,*) " ------------------------------ "
        write(*,*) " "
        write(*,*) "  Orb.   Frag.   FOLI  "
        write(*,*) " --------------------- "
        iselected=0
        do ifrg=1,icufr
          do iorb=1,nb
            xx=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))-xfront

!! CRITERIA FOR SELECTION: sqrt(deloc/Q_A) <= frontier + folitol !!
            if(xx.lt.folitol) then
              iselected=iselected+1
              xx2=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              write(*,112) iorb,ifrg,xx2
              do mu=1,igr
                clindep(mu,iselected)=cfrgoslo(ifrg,mu,iorb)
              end do
            end if
          end do
        end do
        write(*,*) " --------------------- "
        write(*,'(2x,a38,x,i3)') "Number of OSLOs for LinDep evaluation:",iselected
        write(*,*) " "

!! LAST ITERATION NO LINDEP EVALUATION !!
        iilindep=0
        if(nb-nnelect.eq.0.and.inewcore.eq.1) then
          write(*,*) " LinDep not evaluated in last iteration if only 1 orbital is selected "
          write(*,*) " "
          iilindep=1
        end if

!! SOME DEALLOCATES... !!
        DEALLOCATE(deloc)
        DEALLOCATE(infopop)

!! EVALUATE LINDEP !!
        ilindep=0
        if(iilindep.eq.0) then
          SSS=ZERO
          EEE=ZERO
          do ii=1,iselected
            do jj=1,iselected
              xx=ZERO
              do mu=1,igr
                do nu=1,igr
                  xx=xx+clindep(mu,ii)*clindep(nu,jj)*s(mu,nu)
                end do
              end do
              SSS(ii,jj)=xx
            end do
          end do

!! DIAGONALIZING ALL MATRIX, REST IS ZERO SO NO AFFECTS !!
          call diagonalize(nb,nb,SSS,EEE,0)

!! PRINTING SMALLEST EIGENVALUE !!
          xx=10.0d0 !! ABSURT VALUE AGAIN... !!
          do ii=1,iselected
            if(SSS(ii,ii).lt.xx) xx=SSS(ii,ii)
          end do
          if(xx.lt.1.0d-4) ilindep=1 !! THRESHOLD FOR LINIAR DEPENDENCY !!
          write(*,'(2x,a36,x,f10.5)') "Lowest eigenvalue obtained (LinDep):",xx
          write(*,*) " "
        end if
        if(ilindep.eq.1) then
          write(*,*) " *********************************** "
          write(*,*) "  WARNING : LINIAR DEPENDENCY FOUND  "
          write(*,*) " *********************************** "
          write(*,*) " "
          write(*,'(2x,a27,3f10.5)') "FOLI values and delta-FOLI:",xcutoff,xfront,xfront-xcutoff
          write(*,*) " Selecting largest to proceed "
          write(*,*) " RECOMMENDED TO BRANCH (.inp) AND CHECK ALTERNATIVE ASSIGNMENT "
          write(*,*) " "
        end if

!! REMOVING ORBITALS FROM P MATRIX, ONLY IF NOT ALL ARE ASSIGNED !!
!! ORTHOGONALIZING FRAGMENT CORE/SEMICORE ORBITALS !!
        ALLOCATE(S0(inewcore,inewcore),eigv(inewcore,inewcore))
        do ii=1,inewcore
          do jj=1,inewcore
            xx=ZERO
            do mu=1,igr
              do nu=1,igr
                xx=xx+ccore(mu,ii)*ccore(nu,jj)*s(mu,nu)
              end do
            end do
            S0(ii,jj)=xx
          end do
        end do
        
        call diagonalize(inewcore,inewcore,S0,eigv,0)
        
        ALLOCATE(smh(inewcore,inewcore))
        do ii=1,inewcore
          do jj=ii,inewcore
            smh(jj,ii)=ZERO
            do kk=1,inewcore
              if(S0(kk,kk).gt.thresh) then
                xx=eigv(ii,kk)*eigv(jj,kk)
                ssqrt=dsqrt(S0(kk,kk))
                smh(jj,ii)=smh(jj,ii)+xx/ssqrt
              end if
            end do
            smh(ii,jj)=smh(jj,ii)
          end do
        end do
        DEALLOCATE(eigv,S0)

        do ii=1,igr
          do jj=1,inewcore
            xx=ZERO
            do kk=1,inewcore
              xx=xx+smh(jj,kk)*ccore(ii,kk)
            end do
            ccoreorth(ii,jj)=xx
          end do
        end do

!! SAVING ORTHOGONAL ORBITALS HERE !!
        do ii=1,inewcore
          iaddoslo2=iaddoslo2+1
          do mu=1,igr
            cbosloorth(mu,iaddoslo2)=ccoreorth(mu,ii)
          end do
        end do
        DEALLOCATE(smh)

!! CONSTRUCTING pcore (AND pnocore BY SUBSTRACTION) !!
        do ii=1,igr
          do jj=1,igr
            xx=ZERO
            do ij=1,inewcore
              xx=xx+ccoreorth(ii,ij)*ccoreorth(jj,ij)
            end do
            pcore(ii,jj)=xx
          end do
        end do
        do ii=1,igr
          do jj=1,igr
            pnocore(ii,jj)=pnocore(ii,jj)-pcore(ii,jj)
          end do
        end do

!! IN CASE OF ALL ALPHA ELECTRONS ASSIGNED !!
        if(nb-nnelect.eq.0) go to 667

!! END OF BETA ITERATIVE PROCEDURE !!
      end do
667   continue

!! FINAL OS ASSIGNMENT !!
      write(*,*) " --------------------------- "
      write(*,*) "  FRAGMENT OXIDATION STATES  "
      write(*,*) " --------------------------- "
      write(*,*) " "
      write(*,*) "  Frag.  Oxidation State  "
      write(*,*) " ------------------------ "
      do ifrg=1,icufr
        write(*,20) ifrg,REAL(iznfrg(ifrg)-(ifrgel(1,ifrg)+ifrgel(2,ifrg)))
      end do 
      write(*,*) " ------------------------ "
      write(*,*) " "

!! PRINTING OF THE .fchk FILES WITH THE OSLOs... TO VISUALIZE !!
!! PREORTHOGONALIZATION OSLOs CAN BE VISUALIZED IF DESIRED (.inp) !!
      ALLOCATE(poslo(igr,igr))
      if(ifchk.eq.2) then
        do ii=1,igr
          do jj=1,igr
            xx=ZERO
            do ij=1,iaddoslo
              xx=xx+cboslo(ii,ij)*cboslo(jj,ij)
            end do
            poslo(ii,jj)=xx
          end do
        end do
        ctype="-beta-OSLOs-preortho"
        call uwf_orbprint(1,cboslo,poslo,ctype)
      end if

!! NOW THE FINAL (ORTHOGONALIZED) ONES !!
      do ii=1,igr
        do jj=1,igr
          xx=ZERO
          do ij=1,iaddoslo
            xx=xx+cbosloorth(ii,ij)*cbosloorth(jj,ij)
          end do
          poslo(ii,jj)=xx
        end do
      end do
      ctype="-beta-OSLOs"
      call uwf_orbprint(1,cbosloorth,poslo,ctype)
      DEALLOCATE(poslo)

!! EVALUATING FINAL POPULATIONS TO COMPARE !!
      write(*,*) " ---------------------------------- "
      write(*,*) "  PRINTING FINAL OSLOs INFORMATION  "
      write(*,*) " ---------------------------------- "
      write(*,*) " "

!! FIRST ALPHA !!
      write(*,*) " ------------------------------------------------- "
      write(*,*) "  Summary of the selected alpha OSLOs (pre-ortho)  "
      write(*,*) " ------------------------------------------------- "
      write(*,*) " "
      ALLOCATE(orbpop(nalf),orbpop2(nalf))
      ALLOCATE(foslo(nalf,icufr),foslo2(nalf,icufr))
      ALLOCATE(foli(nalf))
      foslo=ZERO
      foslo2=ZERO
      foli=ZERO
      do jfrg=1,icufr
        orbpop=ZERO
        orbpop2=ZERO
        call rwf_uwf_frg_pop(jfrg,sat,nalf,coslo,orbpop) !! FOR THE NON-ORTHOGONAL OSLOs (ORIGINAL) !!
        call rwf_uwf_frg_pop(jfrg,sat,nalf,cosloorth,orbpop2) !! FOR THE ORTHOGONALIZED ONES (PRINTING LATER) !!
        do ii=1,nalf
          foslo(ii,jfrg)=orbpop(ii)
          foslo2(ii,jfrg)=orbpop2(ii)
        end do
      end do

!! DEFINING FOLI MATRIX TO REUSE PRINTING SUBROUTINE !!
      do ii=1,nalf
        foli(ii)=delocoslo(1,ii)
      end do
      call rwf_uwf_print_OSLO_final(1,nalf,foli,foslo)
      write(*,*) " --------------------------------------------- "
      write(*,*) "  Summary of the selected alpha OSLOs (final)  "
      write(*,*) " --------------------------------------------- "
      write(*,*) " "
      call rwf_uwf_print_OSLO_final(0,nalf,foli,foslo2) !! FOLI VALUES GIVEN JUST FOR USING SAME ROUTINE !!
      DEALLOCATE(orbpop,orbpop2)
      DEALLOCATE(foslo,foslo2)
      DEALLOCATE(foli)

!! NOW BETA !!
      write(*,*) " ------------------------------------------------ "
      write(*,*) "  Summary of the selected beta OSLOs (pre-ortho)  "
      write(*,*) " ------------------------------------------------ "
      write(*,*) " "
      ALLOCATE(orbpop(nb),orbpop2(nb))
      ALLOCATE(foslo(nb,icufr),foslo2(nb,icufr))
      ALLOCATE(foli(nb))
      foslo=ZERO
      foslo2=ZERO
      foli=ZERO
      do jfrg=1,icufr
        orbpop=ZERO
        orbpop2=ZERO
        call rwf_uwf_frg_pop(jfrg,sat,nb,cboslo,orbpop) !! FOR THE NON-ORTHOGONAL OSLOs (ORIGINAL) !!
        call rwf_uwf_frg_pop(jfrg,sat,nb,cbosloorth,orbpop2) !! FOR THE ORTHOGONALIZED ONES (PRINTING LATER) !!
        do ii=1,nb
          foslo(ii,jfrg)=orbpop(ii)
          foslo2(ii,jfrg)=orbpop2(ii)
        end do
      end do

!! DEFINING FOLI MATRIX TO REUSE PRINTING SUBROUTINE !!
      do ii=1,nb
        foli(ii)=delocoslo(2,ii)
      end do
      call rwf_uwf_print_OSLO_final(1,nb,foli,foslo)
      write(*,*) " -------------------------------------------- "
      write(*,*) "  Summary of the selected beta OSLOs (final)  "
      write(*,*) " -------------------------------------------- "
      write(*,*) " "
      call rwf_uwf_print_OSLO_final(0,nb,foli,foslo2) !! FOLI VALUES GIVEN JUST FOR USING SAME ROUTINE !!
      DEALLOCATE(orbpop,orbpop2)
      DEALLOCATE(foslo,foslo2)
      DEALLOCATE(foli)

!! TODO: FINAL DEALLOCATE (IF CONSIDERED RELEVANT) !!

!! PRINTING FORMATS !!
110   FORMAT(3x,i3,3x,3f10.5)
111   FORMAT(3x,i3,2x,f10.5,3x,f10.5,3x,f10.5)
112   FORMAT(3x,i3,4x,i3,x,f9.5)
20    FORMAT(3x,i3,6x,f8.2)

      end
  
!! ****** !!

      subroutine uwf_orbprint(iiss,cmat,pmat,ctype)

!! THIS SUBROUTINE PRINTS ORBITALS IN .fchk FORMAT !!
!! NOT ELEGANT WAY TO DO IT, BUT WORKS !!

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /filename/name0

      character*80 line
      character*60 name0,name1
      character*20 ctype

      dimension cmat(igr,igr),pmat(igr,igr)

      iqchem   = iopt(95)
      indepigr = int_locate(15,"Number of independ",ilog)
      norb     = igr*indepigr
      norbt    = igr*(igr+1)/2

!! NAME OF THE .fchk FILE !!
      name1=trim(name0)//trim(ctype)//".fchk"
      open(unit=69,file=name1)
      rewind(69)
      rewind(15)

!! ALPHA CASE !!
      if(iiss.eq.0) then

!! PRINTING UNTIL ALPHA MOs !!
        read(15,'(a80)') line
        do while(index(line,"Alpha MO co").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do

!! PRINTING THE NEW ONES !!
        write(69,11) "Alpha MO coefficients","R","N= ",norb
        write(69,13) ((cmat(ii,jj),ii=1,igr),jj=1,indepigr)

!! NOW LOCATING WHAT IS AFTER IT IN THE ORIGINAL ONE TO CONTINUE !!
        do while(index(line,"Beta MO coef").eq.0)
          read(15,'(a80)') line
        end do

!! BETA CASE !!
      else if(iiss.eq.1) then

!! PRINTING UNTIL ALPHA MOs !!
        read(15,'(a80)') line
        do while(index(line,"Beta MO co").eq.0)
          write(69,'(a80)') line
          read(15,'(a80)') line
        end do

!! PRINTING THE NEW ONES !!
        write(69,11) "Beta MO coefficients ","R","N= ",norb
        write(69,13) ((cmat(ii,jj),ii=1,igr),jj=1,indepigr)

!! NOW LOCATING WHAT IS AFTER IT IN THE ORIGINAL ONE TO CONTINUE !!
        if(iqchem.eq.0) then
          do while(index(line,"Orthonormal basis").eq.0)
            read(15,'(a80)') line
          end do
        else
          do while(index(line,"Alpha Orbital").eq.0)
            read(15,'(a80)') line
          end do
        end if

!! END OF SPIN CASES !! 
      end if

!! RESTART PRINTING UNTIL NEXT STOP !!
!! aquesta part s'hauria de canviar, pero tampoc afecta... !!
!! tambe falta la Spin SCF Dens !!
      do while(index(line,"Total SCF Dens").eq.0)
        write(69,'(a80)') line
        read(15,'(a80)') line
      end do

!! PRINTING THE NEW ONE !!
      write(69,12) "Total SCF Density","R","N= ",norbt
      write(69,13) ((pmat(ii,jj),jj=1,ii),ii=1,igr)

!! NOW LOCATING WHAT IS AFTER IT IN THE ORIGINAL ONE TO CONTINUE !!
      if(iqchem.eq.0) then
        do while(index(line,"Mulliken Charges").eq.0)
          read(15,'(a80)') line
        end do
      else
        do while(index(line,"Pure Switching").eq.0)
          read(15,'(a80)') line
        end do
      end if

!! RESTART PRINTING UNTIL THE END !!
      do while(.true.)
        write(69,'(a80)') line
        read(15,'(a80)',end=99) line
      end do
99    continue
      close(69)

!! PRINTING FORMATS !!
11    FORMAT(a21,22x,a1,3x,a3,i11)
12    FORMAT(a17,26x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))

      end

!! ****** !!

