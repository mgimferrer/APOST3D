
!! **************************** !!
!! OSLO CALCULATION SUBROUTINES !!
!! **************************** !!

      subroutine iterative_oslo_rwf(sat,itotps,wp,omp2,chp,pcoord)

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
      character*5 ccent

      dimension sat(igr,igr,nat)
      dimension chp(itotps,igr),pcoord(itotps,3)
      dimension wp(itotps),omp2(itotps,nat)

      allocatable :: Smat(:,:,:)
      allocatable :: S0(:,:),Sm(:,:),Splus(:,:),smh(:,:),eigv(:,:)
      allocatable :: c0(:,:),pp0(:,:)

      allocatable :: SSS(:,:),EEE(:,:)

      allocatable :: cfrgoslo(:,:,:)
      allocatable :: cmat(:,:),pmat(:,:)
      allocatable :: orbpop(:),orbpop2(:),frgpop(:,:),frgspr(:,:)

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
      folitol = 1.0d0**(-ifolitol) !! DEFAULT = 10^-3, CONTROLLED IN .inp !!

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
      write(*,*) " ------------------------------------- "
      write(*,*) " CHARGE CENTER (R_F) FOR EACH FRAGMENT "
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
        write(*,'(a6,i3,a18,3f10.5)') " Frg: ",ifrg," Chg Cent. (xyz): ",
     + xcenter*angtoau,ycenter*angtoau,zcenter*angtoau

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
      write(*,*) " "

!! PREPARING FOR ITERATIVE PROCESS !!
      ALLOCATE(cmat(igr,igr),pmat(igr,igr))
      ALLOCATE(cfrgoslo(icufr,igr,igr))
      ALLOCATE(frgpop(icufr,nocc))
      ALLOCATE(frgspr(icufr,nocc))
      ALLOCATE(pnocore(igr,igr))

!! to change this part for calling subroutines... !!
      do ii=1,igr
        do jj=1,igr
          pnocore(ii,jj)=pa(ii,jj)
          coslo(ii,jj)=ZERO
          cosloorth(ii,jj)=ZERO
        end do
      end do

!! INITIAL PRINTING !!
      write(*,*) " --------------------------------- "
      write(*,*) " STARTING ITERATIVE OSLO ALGORITHM "
      write(*,*) " --------------------------------- "
      write(*,*) " "
      write(*,*) " Tolerance (in delta-FOLI) used for OSLO selection: ",folitol
      write(*,*) " "

      iaddcore=0 !! CHECK (MG: no recordo perque vaig posar check) !!
      ntotcore=0
      iaddoslo=0
      iaddoslo2=0
      do iiter=1,niter
        write(*,*) " ---------------------- "
        write(*,*) " ITERATION NUMBER ",iiter
        write(*,*) " ---------------------- "
        write(*,*) " "

!! ZEROING THE INVOLVED MATRICES !!
        !! to change also
        do ii=1,igr
          do jj=1,igr
            pcore(ii,jj)=ZERO
            ccore(ii,jj)=ZERO
            ccoreorth(ii,jj)=ZERO
          end do
        end do
        ALLOCATE(deloc(icufr,nocc))
        do ifrg=1,icufr
          do ii=1,nocc
            deloc(ifrg,ii)=ZERO
          end do
        end do

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

!! RECOVERING THE COEFFS, SAVED IN cfrgoslo. cmat AND pmat USED FOR .fchk PRINTING !!
          do kk=1,imaxo
            do mu=1,igr
              cmat(mu,kk)=c0(mu,kk)*dsqrt(pp0(kk,kk))
              cfrgoslo(ifrg,mu,kk)=cmat(mu,kk)
            end do
          end do

!! CONSTRUCTING pmat TO LATER PRINT THE .fchk !!
          do ii=1,igr
            do jj=1,igr
              xx=ZERO
              do ij=1,nocc-ntotcore
                xx=xx+cmat(ii,ij)*cmat(jj,ij)
              end do
              pmat(ii,jj)=TWO*xx
            end do
          end do

!! COMPUTING PIPEK DELOCALIZATION, REQUIRES FRAGMENT POPULATIONS !!
          ALLOCATE(orbpop(nocc))
          do jfrg=1,icufr
            call rwf_frg_pop(jfrg,sat,cmat,orbpop)
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
          write(*,*) " ------------------------------------- "
          write(*,*) " ORBITAL INFORMATION FOR FRAGMENT ",ifrg
          write(*,*) " ------------------------------------- "
          write(*,*) " "
          write(*,*) " Orbital   Spread   Frg. Pop.   FOLI "
          do ii=1,imaxo
            write(*,111) ii,frgspr(ifrg,ii),dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
          end do
          write(*,*) " "
        end do
        DEALLOCATE(S0)
111   FORMAT(x,i3,3x,e6.3,3x,e6.3,3x,e6.3) !! PRINTING FORMAT !!

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
        write(*,*) " Frg. and Lowest FOLI value: ",iifrg,xcutoff
        write(*,*) " Frg. and Lowest FOLI value including tolerance (cutoff): ",jjfrg,xcutoff
        write(*,*) " "
        write(*,*) " ----------------- "
        write(*,*) " SELECTED ORBITALS "
        write(*,*) " ----------------- "
        write(*,*) " "
!! 3) EVALUATING DEGENERACIES !!
!! infopop(i,j): SAVING THE FRAGMENT IN i = 1 AND ORBITAL NUMBER IN i = 2 !!
!! j ALLOCATED AS nocc (FOR PRACTICITY) BUT MAXIMUM WILL BE inewcore !!
        ALLOCATE(infopop(2,nocc))
        inewcore=0

!! BRANCHING (CONTROLLED FROM .inp, DEFAULT = 0) !!
        if(iiter.eq.ibranch) then
          write(*,*) " ******************************************** "
          write(*,*) " WARNING: BRANCHING INVOKED IN THIS ITERATION "
          write(*,*) " ******************************************** "
          write(*,*) " "
!! MG: TO DO !!
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
          do ifrg=1,icufr
            do iorb=1,nocc
              xx=xcutoff-dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              if(ABS(xx).le.folitol) then

!! APPLYING CONDITIONS TO REMOVE ORBITALS !!
                inewcore=inewcore+1
                infopop(1,inewcore)=ifrg
                infopop(2,inewcore)=iorb
                scr(inewcore)=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
                write(*,*) " Orb., Frg., FOLI: ",infopop(2,inewcore),infopop(1,inewcore),scr(inewcore)
              end if
            end do
          end do
        end if
        write(*,*) " Number of OSLOs selected: ",inewcore
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
        write(*,*) " Electron pairs left to assign: ",nocc-nnelect
        write(*,*) " "

        if(nocc-nnelect.lt.0) then
          write(*,*) " *************************************** "
          write(*,*) " WARNING: OVERASSIGNING BY (pairs): ",-(nocc-nnelect)
          write(*,*) " *************************************** "
          write(*,*) " "

!! OS PRINTING OF DESPERATION !!
          write(*,*) " ---------------------- "
          write(*,*) " EOS-like OS ASSIGNMENT "
          write(*,*) " ---------------------- "
          do ifrg=1,icufr
            write(*,'(a12,i3,i3)') " FRAG, OS : ",ifrg,iznfrg(ifrg)-2*ifrgel(ifrg)
          end do
          write(*,*) " "

!! DIRTY TRICK !!
          write(*,*) " CONTINUES BY TRICKING THE CODE (overassigned pairs removed) "
          write(*,*) " "
          inewcore=inewcore+(nocc-nnelect)
          iaddoslo=iaddoslo+(nocc-nnelect)
          nnelect=nnelect+(nocc-nnelect)
        end if

!! SELECTING THE FIRST OUT FOR EVALUATING LINDEP !!
        do mu=1,igr
          do nu=1,igr
            clindep(mu,nu)=ZERO
          end do
        end do
        write(*,*) " ---------------------------- "
        write(*,*) " CHECKING LINIAR DEPENDENCIES "
        write(*,*) " ---------------------------- "
        write(*,*) " "
        iselected=0
        do ifrg=1,icufr
          do iorb=1,nocc
            xx=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))-xfront

!! CRITERIA FOR SELECTION: sqrt(deloc/Q_A) <= frontier + folitol !!
            if(xx.lt.folitol) then
              iselected=iselected+1
              xx2=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              write(*,*) " Orb.,Frag., FOLI: ",iorb,ifrg,xx2
              do mu=1,igr
                clindep(mu,iselected)=cfrgoslo(ifrg,mu,iorb)
              end do
            end if
          end do
        end do
        write(*,*) " Number of OSLOs for LinDep evaluation: ",iselected
        write(*,*) " "

!! SOME DEALLOCATES... !!
        DEALLOCATE(deloc)
        DEALLOCATE(infopop)

!! EVALUATE LINDEP !!
        ilindep=0
        if(iselected.gt.1) then !! IF ONLY 1 THERE IS NOTHING TO EVALUATE !!
          do ii=1,nocc
            do jj=1,nocc
              SSS(ii,jj)=ZERO
              EEE(ii,jj)=ZERO
            end do
          end do
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
          if(xx.lt.1.0d-5) ilindep=1 !! THRESHOLD FOR LINIAR DEPENDENCY !!
          write(*,*) " Lowest eigenvalue obtained (LinDep): ",xx
          write(*,*) " "
        end if
        if(ilindep.eq.1) then
          write(*,*) " ********************************* "
          write(*,*) " WARNING : LINIAR DEPENDENCY FOUND "
          write(*,*) " ********************************* "
          write(*,*) " "
          write(*,*) " FOLI values and delta-FOLI: ",xcutoff,xfront,xcutoff-xfront
          write(*,*) " Selecting largest to proceed "
          write(*,*) " "
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

!! IN CASE OF ALL ASSIGNED !!
        if(nocc-nnelect.eq.0) go to 666

!! END OF ITERATIVE PROCEDURE !!
      end do
666   continue

!! FINAL OS ASSIGNMENT !!
      write(*,*) " ---------------------------- "
      write(*,*) " FINAL EOS-like OS ASSIGNMENT "
      write(*,*) " ---------------------------- "
      write(*,*) " "
      do ifrg=1,icufr
        write(*,'(a12,i3,i3)') " FRAG, OS : ",ifrg,iznfrg(ifrg)-2*ifrgel(ifrg)
      end do
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
        ctype="-non-ortho"
        call rwf_orbprint(0,coslo,poslo,ctype)
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
      ctype="-ortho"
      call rwf_orbprint(0,cosloorth,poslo,ctype)

!! TODO: final printing... !!
!! EVALUATING FINAL POPULATIONS TO COMPARE !!

      write(*,*) " PRINTING FINAL POPULATIONS "
      write(*,*) " "

!! DONE UGLY TO USE THE ROUTINES !!

      ALLOCATE(orbpop(nocc))
      do iorb=1,nocc
        iifrg=infooslo(iorb)
        call rwf_frg_pop2(iifrg,sat,cosloorth,orbpop)
        write(*,*) " Frg., Spread, Pop. (pre), FOLI, Pop. (ortho): ",iifrg,fspread(iorb),poposlo(iorb),delocoslo(iorb),orbpop(iorb)
      end do
      DEALLOCATE(orbpop)
      write(*,*) " "

      write(*,*) " "
      write(*,*) " OSLO populations in all fragments (CHECK)"
      write(*,*) " FOLI values are of the selected one! "
      write(*,*) " "
      ALLOCATE(orbpop(nocc))
      ALLOCATE(orbpop2(nocc))
      do iorb=1,nocc
        write(*,*) " "
        write(*,*) " Orbital",iorb
        write(*,*) " "
        do jfrg=1,icufr
          call rwf_frg_pop2(jfrg,sat,coslo,orbpop)
          call rwf_frg_pop2(jfrg,sat,cosloorth,orbpop2)
          write(*,*) " Frg., Pop. (non-ortho), FOLI, Pop. (ortho): ",jfrg,orbpop(iorb),delocoslo(iorb),orbpop2(iorb)
        end do
        write(*,*) " "
      end do
      DEALLOCATE(orbpop)
      DEALLOCATE(orbpop2)
      write(*,*) " "

!! DEALLOCATING !!
!! TODO !!
!     DEALLOCATE(pcore,pnocore)

      end

!! ****** !!
!! TODO !!
      subroutine iterative_oslo_uwf(sat,itotps,wp,omp2,chp,pcoord)
      use basis_set
      use ao_matrices
      use integration_grid
  
!! fer descripcio !!
  
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
      character*5 ccent
  
      dimension sat(igr,igr,nat) !,ppa(igr,igr),ppb(igr,igr)
      dimension chp(itotps,igr),pcoord(itotps,3)
      dimension wp(itotps),omp2(itotps,nat)
  
      allocatable :: Smat(:,:,:)
      allocatable :: S0(:,:),Sm(:,:),Splus(:,:),smh(:,:),eigv(:,:)
      allocatable :: c0(:,:),pp0(:,:)
  
      allocatable :: SSS(:,:),EEE(:,:)
  
      allocatable :: cfrgoslo(:,:,:)
      allocatable :: cmat(:,:),pmat(:,:)
      allocatable :: orbpop(:),frgpop(:,:),frgspr(:,:)
  
      allocatable :: infopop(:,:),infooslo(:,:),poposlo(:,:)
      allocatable :: fspread(:),fbspread(:),scr(:),scrb(:),deloc(:,:),delocoslo(:,:)
  
      allocatable :: clindep(:,:)
  
      allocatable :: pcore(:,:),pnocore(:,:)
      allocatable :: ccore(:,:),ccoreorth(:,:)
  
      allocatable :: coslo(:,:),cosloorth(:,:)
      allocatable :: cboslo(:,:),cbosloorth(:,:)
      allocatable :: poslo(:,:)
  
      allocatable :: ifrgel(:,:),iznfrg(:)
  
  
      iqchem=iopt(95)
      iatps=nrad*nang
  
      niter=999
      xthresh=1.0d-7 !! DELOCALIZATION THRESHOLD !!
  
!! ALLOCATING THE REQUIRED MATRICES !!
  
      ALLOCATE(Sm(igr,igr),Splus(igr,igr))
      ALLOCATE(c0(igr,igr),pp0(igr,igr))
  
      ALLOCATE(ccore(igr,igr),pcore(igr,igr),ccoreorth(igr,igr))
  
      ALLOCATE(coslo(igr,igr),cosloorth(igr,igr))
      ALLOCATE(cboslo(igr,igr),cbosloorth(igr,igr))
  
      ALLOCATE(clindep(igr,igr))
  
      ALLOCATE(fspread(nalf),fbspread(nb),scr(nalf),scrb(nb),delocoslo(2,igr))
      ALLOCATE(infooslo(2,igr),poposlo(2,igr))
  
      ALLOCATE(ifrgel(2,icufr),iznfrg(icufr))
  
!! FRAGMENT CHARGE FROM zn (AVOIDS PROBLEMS WHEN PSEUDOPOTENTIALS ARE USED) !!
  
      do ifrg=1,icufr
        izn=0
        do icenter=1,nfrlist(ifrg)
          iiat=ifrlist(icenter,ifrg)
          izn=izn+INT(zn(iiat))
        end do
        iznfrg(ifrg)=izn
        ifrgel(1,ifrg)=0
        ifrgel(2,ifrg)=0
      end do
  
!! FRAGMENT LOCALIZATION: IMPORTANT :: R_A VALUE SELECTED IS THE CENTER OF CHARGES BETWEEN FRAGMENT ATOMS !!
!! COMPUTING S ONLY ONCE !!

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
        write(*,*) " "
        write(*,'(a6,i3,a18,3f10.5)') " Frg: ",ifrg," Chg Cent. (xyz): ",xcenter*angtoau,ycenter*angtoau,zcenter*angtoau
        write(*,*) " "

!! MIMICKING PEDROs STRATEGY USING P (FIRST ITERATION WITH CORE) !!

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

!! COMPUTING OSLOs ONCE AND EVALUATING LINDEP FOR THE SET !!

      write(*,*) " "
      write(*,*) " COMPUTING OSLOs AND EVALUATING LINDEP. NO IT. PART "
      write(*,*) " "
      write(*,*) " !!! ALPHA PART !!!"
      write(*,*) " "

      iaddoslo=0
      ntotcore=0

!! ZEROING THE INVOLVED MATRICES !!

      do ii=1,igr
        do jj=1,igr
          pcore(ii,jj)=ZERO
          ccore(ii,jj)=ZERO
          ccoreorth(ii,jj)=ZERO
        end do
      end do
      ALLOCATE(deloc(icufr,nalf))
      do ifrg=1,icufr
        do ii=1,nalf
          deloc(ifrg,ii)=ZERO
        end do
      end do

!! 1) OBTAINING OSLOs FOR ALL FRGS !!

      ALLOCATE(cmat(igr,igr),pmat(igr,igr))
      ALLOCATE(cfrgoslo(icufr,igr,igr))
      ALLOCATE(frgpop(icufr,nalf))
      ALLOCATE(frgspr(icufr,nalf))

      iaddcore=0 !! CHECK
      ALLOCATE(S0(igr,igr))
      do ifrg=1,icufr
        do ii=1,igr
          do jj=1,igr
            pp0(ii,jj)=pa(ii,jj)
            S0(ii,jj)=Smat(ifrg,ii,jj)
          end do
        end do
        call build_Smp(igr,S0,Sm,Splus,0)
        call to_lowdin_basis(igr,Splus,pp0)
        call diagonalize(igr,igr,pp0,c0,0)
        call to_AO_basis(igr,igr,Sm,c0)

!! RECOVERING THE COEFFS, SAVED IN cfrgoslo. cmat AND pmat USED FOR .fchk PRINTING !!

        do kk=1,nalf
          do mu=1,igr
            cmat(mu,kk)=c0(mu,kk)*sqrt(pp0(kk,kk))
            cfrgoslo(ifrg,mu,kk)=cmat(mu,kk)
          end do
        end do
        do ii=1,igr
          do jj=1,igr
            xx=ZERO
            do ij=1,nalf-ntotcore
              xx=xx+cmat(ii,ij)*cmat(jj,ij)
            end do
!           pmat(ii,jj)=xx
          end do
        end do

!! PRINTING SOME INFO !!

        write(*,*) " "
        write(*,'(a10,i3,a15)') " FRAGMENT ",ifrg," ORBITAL INFO "
        write(*,*) " "

!! SPREADS !!

        write(*,'(a10)') " SPREADS: "
        write(*,*) " "
        do ii=1,nalf
          write(*,'(a16,i3,f10.5)') " Orb., Spread : ",ii,pp0(ii,ii)
          frgspr(ifrg,ii)=pp0(ii,ii)
        end do
        write(*,*) " "

!! FRAGMENT POPULATION ANALYSIS + SAVING FOR CUTOFF EVALUATION !!

        write(*,'(a14)') " POPULATIONS: "
        write(*,*) " "
        ALLOCATE(orbpop(nalf))
        call uwf_frg_pop(0,ifrg,sat,cmat,orbpop)
        do ii=1,nalf
          frgpop(ifrg,ii)=orbpop(ii)
        end do
        DEALLOCATE(orbpop)

!! COMPUTING PIPEK DELOCALIZATION !!
!! REQUIRED THE POPULATIONS OF EACH FRAGMENT, DOING IT ONCE AGAIN FOR SIMPLICITY !!

        ALLOCATE(orbpop(nalf))
        do jfrg=1,icufr
          call uwf_frg_pop2(0,jfrg,sat,cmat,orbpop)
          do ii=1,nalf

!! NOW JUST Q_A**2 INSIDE IT !!

            deloc(ifrg,ii)=deloc(ifrg,ii)+orbpop(ii)*orbpop(ii) 
          end do
        end do
        DEALLOCATE(orbpop)

!! DOING 1/deloc() !!
 
        do ii=1,nalf
          if(deloc(ifrg,ii).gt.1.0d-5) then
            deloc(ifrg,ii)=ONE/deloc(ifrg,ii)
          else
            deloc(ifrg,ii)=100.0d0
          end if
        end do

        write(*,*) " OSLO DELOCALIZATION: "
        write(*,*) " "
        do ii=1,nalf
          write(*,*) " Orb., Deloc, Q_A/deloc, sqrt(deloc/Q_A) : ",ii,deloc(ifrg,ii),
     +    frgpop(ifrg,ii)/deloc(ifrg,ii),dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
        end do
        write(*,*) " "
      end do
      DEALLOCATE(S0)
  
!! 2) SELECTING THE FIRST nalf OSLOs WHICH MINIMIZE sqrt(DELOC/Q_A) !!

      xcutoff=100.0d0 !! SET HIGH FOR FIRST STEP !!
      iifrg=0
      iiorb=0
      do ifrg=1,icufr
        do ii=1,nalf
          xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
          if(xx.lt.xcutoff) then
            iiorb=ii
            iifrg=ifrg
            xcutoff=xx
            xxcutoff=frgpop(iifrg,iiorb)/deloc(iifrg,iiorb)
          end if
        end do
      end do
      write(*,*) " Orb. ",iiorb," Frg. ",iifrg," Q_A/Deloc. ",xxcutoff,
     +" sqrt(deloc/Q_A) ",xcutoff
      write(*,*) " "
      write(*,'(a15,f10.5)') " CUTOFF VALUE: ",xcutoff
      write(*,*) " "

!! ADDING ONE BY ONE UNTIL nalf USING THE GO TO STRATEGY !!

123   continue
      xfront=100.0d0 !! SET HIGH FOR FIRST STEP !!
      jjfrg=0
      jjorb=0
      do ifrg=1,icufr
        do ii=1,nalf
          xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
          if(xcutoff+xthresh.lt.xx.and.xx.lt.xfront) then
            jjorb=ii
            jjfrg=ifrg
            xfront=xx
            xxfront=frgpop(jjfrg,jjorb)/deloc(jjfrg,jjorb)
          end if
        end do
      end do
      write(*,*) " 1st out ",jjorb," Frg. ",jjfrg," Q_A/Deloc. ",xxfront,
     +" sqrt(deloc/Q_A) ",xfront
      write(*,*) " "
      write(*,*) " Frontier sqrt(deloc/Q_A) values (in,out) : ",xcutoff,xfront
      write(*,*) " "

!! 3) EVALUATING DEGENERACIES !!
!! infopop(i,j): SAVING THE FRAGMENT IN i = 1 AND ORBITAL NUMBER IN i = 2 !!
!! j ALLOCATED AS nalf BUT MAXIMUM WILL BE inewcore !!

      ALLOCATE(infopop(2,nalf))
      inewcore=0
      do ifrg=1,icufr
        do iorb=1,nalf
          xx=xcutoff-dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
          if(ABS(xx).le.xthresh) then

!! APPLYING CONDITIONS TO REMOVE ORBITALS !!

            inewcore=inewcore+1
            infopop(1,inewcore)=ifrg
            infopop(2,inewcore)=iorb
            scr(inewcore)=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
            write(*,'(a29,i3,a12,i3)') " Orbital Equal/Over Cutoff:",infopop(2,inewcore),
     +      " Fragment : ",infopop(1,inewcore)
          end if
        end do
      end do
      write(*,*) " "
      write(*,'(a39,i3)') " NUMBER OF ORBITALS EQUAL/OVER CUTOFF: ",inewcore
      write(*,*) " "

!! SAVING THE FRAG OSLOs (CONSIDERED CORE) IN ccore !!

      do iorb=1,inewcore
        iifrg=infopop(1,iorb)
        iiorb=infopop(2,iorb)
        ifrgel(1,iifrg)=ifrgel(1,iifrg)+1 !! ADDING THEM HERE !!
        iaddcore=iaddcore+1
        iaddoslo=iaddoslo+1

!! SAVING IN INFOOSLO THE POPULATION AND FRAGMENT PREORTHOGONALIZATION !!

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

!! CHECKING IF WE OVERSELECTED !!

      nnelect=0
      do ifrg=1,icufr
        nnelect=nnelect+ifrgel(1,ifrg)
      end do
      write(*,*) " DUMMY: iaddcore: ",iaddcore
      write(*,*) " Alpha electron left to assign: ",nalf-nnelect
      write(*,*) " "

!! IF WE OVERASSIGN, LINDEP ENSURED !!

      if(nalf-nnelect.lt.0) then
        write(*,*) " !! OVERASSIGNING !! "
        write(*,*) " Overassign by :",-(nalf-nnelect)
        ilindep=1
        DEALLOCATE(infopop) !! posat aqui per si fa falta treureli info !!
        go to 99
      end if

!! ADDING MORE ORBITALS !!

      write(*,*) " "
      write(*,*) " Continue adding orbitals "
      write(*,*) " "
      xcutoff=xfront
      DEALLOCATE(infopop) !! posat aqui per si fa falta treureli info !!
      if(nalf-nnelect.gt.0) then
        go to 123
      else
        write(*,*) " !! ALL ORBITALS SELECTED !! "
        write(*,*) " "
      end if

!! EVALUATING LINIAR DEPENDENCIES !!

      write(*,*) " "
      write(*,*) " LINIAR DEPENDENCY CHECK "
      write(*,*) " "

      ALLOCATE(SSS(nalf,nalf),EEE(nalf,nalf))

      do ii=1,nalf
        do jj=1,nalf
          SSS(ii,jj)=ZERO
          EEE(ii,jj)=ZERO
        end do
      end do
      do ii=1,iaddoslo
        do jj=1,iaddoslo
          xx=ZERO
          do mu=1,igr
            do nu=1,igr
              xx=xx+coslo(mu,ii)*coslo(nu,jj)*s(mu,nu)
            end do
          end do
          SSS(ii,jj)=xx
        end do
      end do

!! DIAGONALIZING ALL MATRIX, REST ARE ZEROS !!

      call diagonalize(nalf,nalf,SSS,EEE,0)

!! ONLY PRINTING THE iaddoslo FIRST, REST FOR SURE ARE ZERO !!

      do ii=1,iaddoslo
        write(*,*) " Orb., Eigenvalue : ",ii,SSS(ii,ii) 
      end do
      write(*,*) " "
      ilindep=0
      do ii=1,iaddoslo
        if(SSS(ii,ii).lt.1.0d-5) then
          write(*,*) " WARNING : LINIAR DEPENDENCY FOUND "
          write(*,*) " sqrt(deloc/Q_A) values (and diff): ",xcutoff,xfront,xcutoff-xfront
          ilindep=1
        end if
      end do
      DEALLOCATE(SSS,EEE)

99    continue

!! NON-ITERATIVE PART WARNINGS !!

      if(ilindep.eq.1) then
        write(*,*) " **************************************** "
        write(*,*) "               ALPHA PART                 "
        write(*,*) " "
        write(*,*) " WARNING: Liniar Dependency Will Be Found "
        write(*,*) " "
        write(*,*) " **************************************** "
      else
        write(*,*) " **************************************** "
        write(*,*) "               ALPHA PART                 "
        write(*,*) " "
        write(*,*) " INFO: NO Liniar Dependency Will Be Found "
        write(*,*) " "
        write(*,*) " **************************************** "
      end if

      DEALLOCATE(deloc)
      DEALLOCATE(frgpop,frgspr)

!! NOW BETA !!

      write(*,*) " "
      write(*,*) " !!! BETA PART !!!"
      write(*,*) " "

      iaddoslo=0
      ntotcore=0

!! ZEROING THE INVOLVED MATRICES !!

      do ii=1,igr
        do jj=1,igr
          pcore(ii,jj)=ZERO
          ccore(ii,jj)=ZERO
          ccoreorth(ii,jj)=ZERO
        end do
      end do
      ALLOCATE(deloc(icufr,nb))
      do ifrg=1,icufr
        do ii=1,nb
          deloc(ifrg,ii)=ZERO
        end do
      end do

!! 1) OBTAINING OSLOs FOR ALL FRGS !!

      ALLOCATE(frgpop(icufr,nb))
      ALLOCATE(frgspr(icufr,nb))

      iaddcore=0 !! CHECK
      ALLOCATE(S0(igr,igr))
      do ifrg=1,icufr
        do ii=1,igr
          do jj=1,igr
            pp0(ii,jj)=pb(ii,jj)
            S0(ii,jj)=Smat(ifrg,ii,jj)
          end do
        end do
        call build_Smp(igr,S0,Sm,Splus,0)
        call to_lowdin_basis(igr,Splus,pp0)
        call diagonalize(igr,igr,pp0,c0,0)
        call to_AO_basis(igr,igr,Sm,c0)

!! RECOVERING THE COEFFS, SAVED IN cfrgoslo. cmat AND pmat USED FOR .fchk PRINTING !!

        do kk=1,nalf
          do mu=1,igr
            cmat(mu,kk)=c0(mu,kk)*sqrt(pp0(kk,kk))
            cfrgoslo(ifrg,mu,kk)=cmat(mu,kk)
          end do
        end do
        do ii=1,igr
          do jj=1,igr
            xx=ZERO
            do ij=1,nb-ntotcore
              xx=xx+cmat(ii,ij)*cmat(jj,ij)
            end do
!           pmat(ii,jj)=xx
          end do
        end do

!! PRINTING SOME INFO !!

        write(*,*) " "
        write(*,'(a10,i3,a15)') " FRAGMENT ",ifrg," ORBITAL INFO "
        write(*,*) " "

!! SPREADS !!

        write(*,'(a10)') " SPREADS: "
        write(*,*) " "
        do ii=1,nb
          write(*,'(a16,i3,f10.5)') " Orb., Spread : ",ii,pp0(ii,ii)
          frgspr(ifrg,ii)=pp0(ii,ii)
        end do
        write(*,*) " "

!! FRAGMENT POPULATION ANALYSIS + SAVING FOR CUTOFF EVALUATION !!

        write(*,'(a14)') " POPULATIONS: "
        write(*,*) " "
        ALLOCATE(orbpop(nb))
        call uwf_frg_pop(1,ifrg,sat,cmat,orbpop)
        do ii=1,nb
          frgpop(ifrg,ii)=orbpop(ii)
        end do
        DEALLOCATE(orbpop)

!! COMPUTING PIPEK DELOCALIZATION !!
!! REQUIRED THE POPULATIONS OF EACH FRAGMENT, DOING IT ONCE AGAIN FOR SIMPLICITY !!

        ALLOCATE(orbpop(nb))
        do jfrg=1,icufr
          call uwf_frg_pop2(1,jfrg,sat,cmat,orbpop)
          do ii=1,nb

!! NOW JUST Q_A**2 INSIDE IT !!

            deloc(ifrg,ii)=deloc(ifrg,ii)+orbpop(ii)*orbpop(ii) 
          end do
        end do
        DEALLOCATE(orbpop)

!! DOING 1/deloc() !!
 
        do ii=1,nb
          if(deloc(ifrg,ii).gt.1.0d-5) then
            deloc(ifrg,ii)=ONE/deloc(ifrg,ii)
          else
            deloc(ifrg,ii)=100.0d0
          end if
        end do

        write(*,*) " OSLO DELOCALIZATION: "
        write(*,*) " "
        do ii=1,nb
          write(*,*) " Orb., Deloc, Q_A/deloc, sqrt(deloc/Q_A) : ",ii,deloc(ifrg,ii),
     +    frgpop(ifrg,ii)/deloc(ifrg,ii),dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
        end do
        write(*,*) " "
      end do
      DEALLOCATE(S0)

!! 2) SELECTING THE FIRST nb OSLOs WHICH MINIMIZE sqrt(DELOC/Q_A) !!

      xcutoff=100.0d0 !! SET HIGH FOR FIRST STEP !!
      iifrg=0
      iiorb=0
      do ifrg=1,icufr
        do ii=1,nb
          xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
          if(xx.lt.xcutoff) then
            iiorb=ii
            iifrg=ifrg
            xcutoff=xx
            xxcutoff=frgpop(iifrg,iiorb)/deloc(iifrg,iiorb)
          end if
        end do
      end do
      write(*,*) " Orb. ",iiorb," Frg. ",iifrg," Q_A/Deloc. ",xxcutoff,
     +" sqrt(deloc/Q_A) ",xcutoff
      write(*,*) " "
      write(*,'(a15,f10.5)') " CUTOFF VALUE: ",xcutoff
      write(*,*) " "

!! ADDING ONE BY ONE UNTIL nb USING THE GO TO STRATEGY !!

122   continue
      xfront=100.0d0 !! SET HIGH FOR FIRST STEP !!
      jjfrg=0
      jjorb=0
      do ifrg=1,icufr
        do ii=1,nb
          xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
          if(xcutoff+xthresh.lt.xx.and.xx.lt.xfront) then
            jjorb=ii
            jjfrg=ifrg
            xfront=xx
            xxfront=frgpop(jjfrg,jjorb)/deloc(jjfrg,jjorb)
          end if
        end do
      end do
      write(*,*) " 1st out ",jjorb," Frg. ",jjfrg," Q_A/Deloc. ",xxfront,
     +" sqrt(deloc/Q_A) ",xfront
      write(*,*) " "
      write(*,*) " Frontier sqrt(deloc/Q_A) values (in,out) : ",xcutoff,xfront
      write(*,*) " "

!! 3) EVALUATING DEGENERACIES !!
!! infopop(i,j): SAVING THE FRAGMENT IN i = 1 AND ORBITAL NUMBER IN i = 2 !!
!! j ALLOCATED AS nb BUT MAXIMUM WILL BE inewcore !!

      ALLOCATE(infopop(2,nb))
      inewcore=0
      do ifrg=1,icufr
        do iorb=1,nb
          xx=xcutoff-dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
          if(ABS(xx).le.xthresh) then

!! APPLYING CONDITIONS TO REMOVE ORBITALS !!

            inewcore=inewcore+1
            infopop(1,inewcore)=ifrg
            infopop(2,inewcore)=iorb
            scrb(inewcore)=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
            write(*,'(a29,i3,a12,i3)') " Orbital Equal/Over Cutoff:",infopop(2,inewcore),
     +      " Fragment : ",infopop(1,inewcore)
          end if
        end do
      end do
      write(*,*) " "
      write(*,'(a39,i3)') " NUMBER OF ORBITALS EQUAL/OVER CUTOFF: ",inewcore
      write(*,*) " "

!! SAVING THE FRAG OSLOs (CONSIDERED CORE) IN ccore !!

      do iorb=1,inewcore
        iifrg=infopop(1,iorb)
        iiorb=infopop(2,iorb)
        ifrgel(2,iifrg)=ifrgel(2,iifrg)+1 !! ADDING THEM HERE !!
        iaddcore=iaddcore+1
        iaddoslo=iaddoslo+1

!! SAVING IN INFOOSLO THE POPULATION AND FRAGMENT PREORTHOGONALIZATION !!

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

!! CHECKING IF WE OVERSELECTED !!

      nnelect=0
      do ifrg=1,icufr
        nnelect=nnelect+ifrgel(2,ifrg)
      end do
      write(*,*) " DUMMY: iaddcore: ",iaddcore
      write(*,*) " Beta electron left to assign: ",nb-nnelect
      write(*,*) " "

!! IF WE OVERASSIGN, LINDEP ENSURED !!

      if(nb-nnelect.lt.0) then
        write(*,*) " !! OVERASSIGNING !! "
        write(*,*) " Overassign by :",-(nb-nnelect)
        ilindep=1
        DEALLOCATE(infopop) !! posat aqui per si fa falta treureli info !!
        go to 98
      end if

!! ADDING MORE ORBITALS !!

      write(*,*) " "
      write(*,*) " Continue adding orbitals "
      write(*,*) " "
      xcutoff=xfront
      DEALLOCATE(infopop) !! posat aqui per si fa falta treureli info !!
      if(nb-nnelect.gt.0) then
        go to 122
      else
        write(*,*) " !! ALL ORBITALS SELECTED !! "
        write(*,*) " "
      end if

!! EVALUATING LINIAR DEPENDENCIES !!

      write(*,*) " "
      write(*,*) " LINIAR DEPENDENCY CHECK "
      write(*,*) " "

      ALLOCATE(SSS(nb,nb),EEE(nb,nb))

      do ii=1,nb
        do jj=1,nb
          SSS(ii,jj)=ZERO
          EEE(ii,jj)=ZERO
        end do
      end do
      do ii=1,iaddoslo
        do jj=1,iaddoslo
          xx=ZERO
          do mu=1,igr
            do nu=1,igr
              xx=xx+coslo(mu,ii)*coslo(nu,jj)*s(mu,nu)
            end do
          end do
          SSS(ii,jj)=xx
        end do
      end do

!! DIAGONALIZING ALL MATRIX, REST ARE ZEROS !!

      call diagonalize(nb,nb,SSS,EEE,0)

!! ONLY PRINTING THE iaddoslo FIRST, REST FOR SURE ARE ZERO !!

      do ii=1,iaddoslo
        write(*,*) " Orb., Eigenvalue : ",ii,SSS(ii,ii) 
      end do
      write(*,*) " "
      ilindep=0
      do ii=1,iaddoslo
        if(SSS(ii,ii).lt.1.0d-5) then
          write(*,*) " WARNING : LINIAR DEPENDENCY FOUND "
          write(*,*) " sqrt(deloc/Q_A) values (and diff): ",xcutoff,xfront,xcutoff-xfront
          ilindep=1
        end if
      end do
      DEALLOCATE(SSS,EEE)

98    continue

!! NON-ITERATIVE PART WARNINGS !!

      if(ilindep.eq.1) then
        write(*,*) " **************************************** "
        write(*,*) "                BETA PART                 "
        write(*,*) " "
        write(*,*) " WARNING: Liniar Dependency Will Be Found "
        write(*,*) " "
        write(*,*) " **************************************** "
      else
        write(*,*) " **************************************** "
        write(*,*) "                BETA PART                 "
        write(*,*) " "
        write(*,*) " INFO: NO Liniar Dependency Will Be Found "
        write(*,*) " "
        write(*,*) " **************************************** "
      end if

      DEALLOCATE(deloc)
      DEALLOCATE(frgpop,frgspr)

!! OS ASSIGNMENT !!

      write(*,*) " "
      write(*,*) " EOS-like OS ASSIGNMENT "
      write(*,*) " "
      do ifrg=1,icufr
        write(*,'(a12,i3,i3)') " FRAG, OS : ",ifrg,iznfrg(ifrg)-(ifrgel(1,ifrg)+ifrgel(2,ifrg))

!! RESETING FOR ITERATIVE PROCESS !!

        ifrgel(1,ifrg)=0
        ifrgel(2,ifrg)=0
      end do
      write(*,*) " "

!! Okay... !!
!! PREPARING FOR ITERATIVE PROCESS !!

      ALLOCATE(frgpop(icufr,nalf))
      ALLOCATE(frgspr(icufr,nalf))
      ALLOCATE(SSS(nalf,nalf),EEE(nalf,nalf))

      ALLOCATE(pnocore(igr,igr))
      do ii=1,igr
        do jj=1,igr
          pnocore(ii,jj)=pa(ii,jj)
          coslo(ii,jj)=ZERO
          cosloorth(ii,jj)=ZERO
        end do
      end do

      write(*,*) " "
      write(*,*) " STARTING ITERATIVE PROCEDURE "
      write(*,*) " "
      write(*,*) " !!! ALPHA PART !!!"
      write(*,*) " "
      write(*,*) " "

!! ITERATIVE PROCESS START !!

      xthresh=1.0d-3 !! SWITCHING TO THE 10^-2 THRESH, MATCHING WITH ABDUL (QCHEM) !!
      ntotcore=0
      iaddoslo=0
      iaddoslo2=0
      do iiter=1,niter
        write(*,*) " "
        write(*,'(a12,i3)') " ITERATION: ",iiter
        write(*,*) " "

!! ZEROING THE INVOLVED MATRICES !!

        do ii=1,igr
          do jj=1,igr
            pcore(ii,jj)=ZERO
            ccore(ii,jj)=ZERO
            ccoreorth(ii,jj)=ZERO
          end do
        end do
        ALLOCATE(deloc(icufr,nalf))
        do ifrg=1,icufr
          do ii=1,nalf
            deloc(ifrg,ii)=ZERO
          end do
        end do

!! 1) OBTAINING OSLOs FOR ALL FRGS !!

        iaddcore=0
        ALLOCATE(S0(igr,igr))
        do ifrg=1,icufr
          do ii=1,igr
            do jj=1,igr
              pp0(ii,jj)=pnocore(ii,jj) !! TO NOT DESTROY pnocore and Smat !!
              S0(ii,jj)=Smat(ifrg,ii,jj)
            end do
          end do
          call build_Smp(igr,S0,Sm,Splus,0)
          call to_lowdin_basis(igr,Splus,pp0)
          call diagonalize(igr,igr,pp0,c0,0)
          call to_AO_basis(igr,igr,Sm,c0)
          imaxo=nalf !! MAX NUMBER OF OSLOs IS NOW nalf BUT SHOULD BE nalf-ntotcore !!

!! RECOVERING THE COEFFS, SAVED IN cfrgoslo. cmat AND pmat USED FOR .fchk PRINTING !!

          do kk=1,imaxo
            do mu=1,igr
              cmat(mu,kk)=c0(mu,kk)*sqrt(pp0(kk,kk))
              cfrgoslo(ifrg,mu,kk)=cmat(mu,kk)
            end do
          end do

!! CONSTRUCTING pmat FOR .fchk CREATING !!

          do ii=1,igr
            do jj=1,igr
              xx=ZERO
              do ij=1,nalf-ntotcore
                xx=xx+cmat(ii,ij)*cmat(jj,ij)
              end do
!             pmat(ii,jj)=xx
            end do
          end do

!! SOME INFO !!

          write(*,*) " "
          write(*,'(a10,i3,a15)') " FRAGMENT ",ifrg," ORBITAL INFO "
          write(*,*) " "

!! SPREADS !!

          write(*,'(a18)') " ORBITAL SPREADS: "
          write(*,*) " "
          do ii=1,imaxo
            write(*,'(a16,i3,f10.5)') " Orb., Spread : ",ii,pp0(ii,ii)
            frgspr(ifrg,ii)=pp0(ii,ii)
          end do
          write(*,*) " "

!! FRAGMENT POPULATION ANALYSIS + SAVING FOR CUTOFF EVALUATION !!

          write(*,'(a22)') " ORBITAL POPULATIONS: "
          write(*,*) " "
          ALLOCATE(orbpop(nalf))
          call uwf_frg_pop(0,ifrg,sat,cmat,orbpop)
          do ii=1,nalf
            frgpop(ifrg,ii)=orbpop(ii)
          end do
          DEALLOCATE(orbpop)

!! COMPUTING PIPEK DELOCALIZATION !!
!! REQUIRED THE POPULATIONS OF EACH FRAGMENT, DOING IT ONCE AGAIN FOR SIMPLICITY !!

          ALLOCATE(orbpop(nalf))
          do jfrg=1,icufr
            call uwf_frg_pop2(0,jfrg,sat,cmat,orbpop)
            do ii=1,imaxo

!! NOW JUST Q_A**2 INSIDE IT !!

              deloc(ifrg,ii)=deloc(ifrg,ii)+orbpop(ii)*orbpop(ii) 
            end do
          end do
          DEALLOCATE(orbpop)

!! DOING 1/deloc() !!
 
          do ii=1,imaxo
            if(deloc(ifrg,ii).gt.1.0d-5) then
              deloc(ifrg,ii)=ONE/deloc(ifrg,ii)
            else
              deloc(ifrg,ii)=100.0d0
            end if
          end do

          write(*,*) " OSLO DELOCALIZATION: "
          write(*,*) " "
          do ii=1,imaxo
            write(*,*) " Orb., Deloc, Q_A/deloc, sqrt(deloc/Q_A) : ",ii,deloc(ifrg,ii),
     +      frgpop(ifrg,ii)/deloc(ifrg,ii),dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
          end do
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
              xcutoff=xx
              xxcutoff=frgpop(iifrg,iiorb)/deloc(iifrg,iiorb)
            end if
          end do
        end do
        write(*,*) " Orb. ",iiorb," Frg. ",iifrg," Q_A/Deloc. ",xxcutoff,
     +  " sqrt(deloc/Q_A) ",xcutoff
        write(*,*) " "
        write(*,'(a15,f10.5)') " CUTOFF VALUE: ",xcutoff
        write(*,*) " "

!! NOW FRONTIER !!

        xfront=100.0d0 !! SET HIGH FOR FIRST STEP !!
        jjfrg=0
        jjorb=0
        do ifrg=1,icufr
          do ii=1,nalf
            xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
            if(xcutoff+xthresh.lt.xx.and.xx.lt.xfront) then
              jjorb=ii
              jjfrg=ifrg
              xfront=xx
              xxfront=frgpop(jjfrg,jjorb)/deloc(jjfrg,jjorb)
            end if
          end do
        end do
        write(*,*) " 1st out ",jjorb," Frg. ",jjfrg," Q_A/Deloc. ",xxfront,
     +  " sqrt(deloc/Q_A) ",xfront
        write(*,*) " "
        write(*,*) " Frontier sqrt(deloc/Q_A) values (in,out) : ",xcutoff,xfront
        write(*,*) " "

!! 3) EVALUATING DEGENERACIES !!
!! infopop(i,j): SAVING THE FRAGMENT IN i = 1 AND ORBITAL NUMBER IN i = 2 !!
!! j ALLOCATED AS nalf BUT MAXIMUM WILL BE inewcore !!

        ALLOCATE(infopop(2,nalf))
        inewcore=0
        do ifrg=1,icufr
          do iorb=1,nalf
            xx=xcutoff-dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
            if(ABS(xx).le.xthresh) then

!! APPLYING CONDITIONS TO REMOVE ORBITALS !!

              inewcore=inewcore+1
              infopop(1,inewcore)=ifrg
              infopop(2,inewcore)=iorb
              scr(inewcore)=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              write(*,'(a29,i3,a12,i3)') " Orbital Equal/Over Cutoff:",infopop(2,inewcore),
     +        " Fragment : ",infopop(1,inewcore)
            end if
          end do
        end do
        write(*,*) " "
        write(*,'(a39,i3)') " NUMBER OF ORBITALS EQUAL/OVER CUTOFF: ",inewcore
        write(*,*) " "

!! SAVING THE FRAG OSLOs (CONSIDERED CORE) IN ccore !!

        do iorb=1,inewcore
          iifrg=infopop(1,iorb)
          iiorb=infopop(2,iorb)
          ifrgel(1,iifrg)=ifrgel(1,iifrg)+1 !! ADDING THEM HERE !!
          iaddcore=iaddcore+1
          iaddoslo=iaddoslo+1

!! SAVING IN INFOOSLO THE POPULATION AND FRAGMENT PREORTHOGONALIZATION !!

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

!! EVALUATING IF WE OVERASSIGN !!

        nnelect=0
        do ifrg=1,icufr
          nnelect=nnelect+ifrgel(1,ifrg)
        end do
        write(*,*) " DUMMY: iaddcore: ",iaddcore
        write(*,*) " Alpha electron left to assign: ",nalf-nnelect
        write(*,*) " "

        if(nalf-nnelect.lt.0) then
          write(*,*) " !! OVERASSIGNING !!"
          write(*,*) " Overassignation by :",-(nalf-nnelect)
!! lets trick the code !!
          write(*,*) "Tricking the code, removing last "
          inewcore=inewcore+(nalf-nnelect)
          iaddoslo=iaddoslo+(nalf-nnelect)
          nnelect=nnelect+(nalf-nnelect)
        end if

!! SELECTING THE FIRST OUT FOR EVALUATING LINDEP !!

        do mu=1,igr
          do nu=1,igr
            clindep(mu,nu)=ZERO
          end do
        end do
        write(*,*) " "
        write(*,*) " SELECTING ORBITALS FOR LINDEP EVALUATION "
        write(*,*) " "
        iselected=0
        do ifrg=1,icufr
          do iorb=1,nalf
            xx=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))-xfront

!! CRITERIA FOR SELECTION: sqrt(deloc/Q_A) <= frontier + xthresh !!

            if(xx.lt.xthresh) then
              iselected=iselected+1
              xx2=frgpop(ifrg,iorb)/deloc(ifrg,iorb)
              xx3=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              write(*,*) " Orbital,Fragment,Q_A/deloc,sqrt(deloc/Q_A): ",iorb,ifrg,xx2,xx3
              do mu=1,igr
                clindep(mu,iselected)=cfrgoslo(ifrg,mu,iorb)
              end do
            end if
          end do
        end do
        write(*,*) " "
        write(*,*) " Selected orbitals for LinDep evaluation: ",iselected
        write(*,*) " "

        DEALLOCATE(deloc)
        DEALLOCATE(infopop)

!! EVALUATE LINDEP !!

        if(iselected.gt.0) then
          write(*,*) " LINIAR DEPENDENCY CHECK "
          write(*,*) " "
          do ii=1,nalf
            do jj=1,nalf
              SSS(ii,jj)=ZERO
              EEE(ii,jj)=ZERO
            end do
          end do
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

!! PRINTING ONLY FIRST iselected !!

          do ii=1,iselected
            write(*,*) " Orb., Eigenvalue : ",ii,SSS(ii,ii) 
          end do
          write(*,*) " "
          ilindep=0
          do ii=1,iselected
            if(SSS(ii,ii).lt.1.0d-5) ilindep=1
          end do
        end if

        if(ilindep.eq.1) then
          write(*,*) " WARNING : LINIAR DEPENDENCY FOUND "
          write(*,*) " sqrt(deloc/Q_A) values (and diff): ",xcutoff,xfront,xcutoff-xfront
          write(*,*) " TO DATE: Selecting largest one and continue "
          write(*,*) " "
        else
          write(*,*) " No lindep found in this step "
          write(*,*) " !! Iterative proceeds Safely !! "
          write(*,*) " "
        end if

!! ONLY IF NOT ALL ARE ASSIGNED !!
!! REMOVING ORBITALS FROM P MATRIX !!

        write(*,*) " "
        write(*,*) " Orbitals to remove in this step: ",inewcore
        write(*,*) " "

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

!! CHECKING BY RECONSTRUCTING PS !!

        write(*,*) " CHECK : COMPUTING trace(PS) "
        write(*,*) " "
        xx=ZERO
        xx2=ZERO
        do ii=1,igr
          do jj=1,igr
            xx=xx+pnocore(ii,jj)*s(jj,ii)
            do iat=1,nat
              xx2=xx2+pnocore(ii,jj)*sat(jj,ii,iat)
            end do
          end do
        end do
        write(*,'(a13,f10.5)') " Trace(PS) : ",xx 
        write(*,'(a24,f10.5)') " Trace(sum_nat(PSat)) : ",xx2 
        write(*,*) " " 

        if(nalf-nnelect.eq.0) then
          write(*,*) " "
          write(*,*) "ALL ALPHA ELECTRONS ASSIGNED"
          write(*,*) " "
          go to 666
        end if

!! END DO iiter !!

      end do
666   continue

      write(*,*) " "
      write(*,*) " (ALPHA) END ITERATIVE PROCEDURE "
      write(*,*) " "

!! PRINTING OSLOs IN fchk !!
!! CONSTRUCTING poslo FOR .fchk CREATING !!

      ALLOCATE(poslo(igr,igr))
      do ii=1,igr
        do jj=1,igr
          xx=ZERO
          do ij=1,iaddoslo
            xx=xx+coslo(ii,ij)*coslo(jj,ij)
          end do
          poslo(ii,jj)=xx
        end do
      end do
      ctype="-alpha-final-nonortho"
      call uwf_orbprint(0,0,coslo,poslo,ctype)

!! PRINTING ORTHO OSLOs IN fchk !!
!! CONSTRUCTING poslo FOR .fchk CREATING !!

      do ii=1,igr
        do jj=1,igr
          xx=ZERO
          do ij=1,iaddoslo
            xx=xx+cosloorth(ii,ij)*cosloorth(jj,ij)
          end do
          poslo(ii,jj)=xx
        end do
      end do
      ctype="-alpha-final-ORTHO"
      call uwf_orbprint(0,0,cosloorth,poslo,ctype)

      DEALLOCATE(poslo)
      DEALLOCATE(frgpop,frgspr)
      DEALLOCATE(SSS,EEE)

!! NOW BETA !!

      ALLOCATE(frgpop(icufr,nb))
      ALLOCATE(frgspr(icufr,nb))
      ALLOCATE(SSS(nb,nb),EEE(nb,nb))

      do ii=1,igr
        do jj=1,igr
          pnocore(ii,jj)=pb(ii,jj)
          cboslo(ii,jj)=ZERO
          cbosloorth(ii,jj)=ZERO
        end do
      end do

      write(*,*) " "
      write(*,*) " STARTING ITERATIVE PROCEDURE "
      write(*,*) " "
      write(*,*) " !!! BETA PART !!!"
      write(*,*) " "
      write(*,*) " "

!! ITERATIVE PROCESS START !!

      ntotcore=0
      iaddoslo=0
      iaddoslo2=0
      do iiter=1,niter
        write(*,*) " "
        write(*,'(a12,i3)') " ITERATION: ",iiter
        write(*,*) " "

!! ZEROING THE INVOLVED MATRICES !!

        do ii=1,igr
          do jj=1,igr
            pcore(ii,jj)=ZERO
            ccore(ii,jj)=ZERO
            ccoreorth(ii,jj)=ZERO
          end do
        end do
        ALLOCATE(deloc(icufr,nb))
        do ifrg=1,icufr
          do ii=1,nb
            deloc(ifrg,ii)=ZERO
          end do
        end do

!! 1) OBTAINING OSLOs FOR ALL FRGS !!

        iaddcore=0
        ALLOCATE(S0(igr,igr))
        do ifrg=1,icufr
          do ii=1,igr
            do jj=1,igr
              pp0(ii,jj)=pnocore(ii,jj) !! TO NOT DESTROY pnocore and Smat !!
              S0(ii,jj)=Smat(ifrg,ii,jj)
            end do
          end do
          call build_Smp(igr,S0,Sm,Splus,0)
          call to_lowdin_basis(igr,Splus,pp0)
          call diagonalize(igr,igr,pp0,c0,0)
          call to_AO_basis(igr,igr,Sm,c0)
          imaxo=nb !! MAX NUMBER OF OSLOs IS NOW nb BUT SHOULD BE nb-ntotcore !!

!! RECOVERING THE COEFFS, SAVED IN cfrgoslo. cmat AND pmat USED FOR .fchk PRINTING !!

          do kk=1,imaxo
            do mu=1,igr
              cmat(mu,kk)=c0(mu,kk)*sqrt(pp0(kk,kk))
              cfrgoslo(ifrg,mu,kk)=cmat(mu,kk)
            end do
          end do

!! CONSTRUCTING pmat FOR .fchk CREATING !!

          do ii=1,igr
            do jj=1,igr
              xx=ZERO
              do ij=1,nb-ntotcore
                xx=xx+cmat(ii,ij)*cmat(jj,ij)
              end do
!             pmat(ii,jj)=xx
            end do
          end do

!! SOME INFO !!

          write(*,*) " "
          write(*,'(a10,i3,a15)') " FRAGMENT ",ifrg," ORBITAL INFO "
          write(*,*) " "

!! SPREADS !!

          write(*,'(a18)') " ORBITAL SPREADS: "
          write(*,*) " "
          do ii=1,imaxo
            write(*,'(a16,i3,f10.5)') " Orb., Spread : ",ii,pp0(ii,ii)
            frgspr(ifrg,ii)=pp0(ii,ii)
          end do
          write(*,*) " "

!! FRAGMENT POPULATION ANALYSIS + SAVING FOR CUTOFF EVALUATION !!

          write(*,'(a22)') " ORBITAL POPULATIONS: "
          write(*,*) " "
          ALLOCATE(orbpop(nb))
          call uwf_frg_pop(1,ifrg,sat,cmat,orbpop)
          do ii=1,nb
            frgpop(ifrg,ii)=orbpop(ii)
          end do
          DEALLOCATE(orbpop)

!! COMPUTING PIPEK DELOCALIZATION !!
!! REQUIRED THE POPULATIONS OF EACH FRAGMENT, DOING IT ONCE AGAIN FOR SIMPLICITY !!

          ALLOCATE(orbpop(nb))
          do jfrg=1,icufr
            call uwf_frg_pop2(1,jfrg,sat,cmat,orbpop)
            do ii=1,imaxo

!! NOW JUST Q_A**2 INSIDE IT !!

              deloc(ifrg,ii)=deloc(ifrg,ii)+orbpop(ii)*orbpop(ii) 
            end do
          end do
          DEALLOCATE(orbpop)

!! DOING 1/deloc() !!
 
          do ii=1,imaxo
            if(deloc(ifrg,ii).gt.1.0d-5) then
              deloc(ifrg,ii)=ONE/deloc(ifrg,ii)
            else
              deloc(ifrg,ii)=100.0d0
            end if
          end do

          write(*,*) " OSLO DELOCALIZATION: "
          write(*,*) " "
          do ii=1,imaxo
            write(*,*) " Orb., Deloc, Q_A/deloc, sqrt(deloc/Q_A) : ",ii,deloc(ifrg,ii),
     +      frgpop(ifrg,ii)/deloc(ifrg,ii),dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
          end do
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
              xcutoff=xx
              xxcutoff=frgpop(iifrg,iiorb)/deloc(iifrg,iiorb)
            end if
          end do
        end do
        write(*,*) " Orb. ",iiorb," Frg. ",iifrg," Q_A/Deloc. ",xxcutoff,
     +  " sqrt(deloc/Q_A) ",xcutoff
        write(*,*) " "
        write(*,'(a15,f10.5)') " CUTOFF VALUE: ",xcutoff
        write(*,*) " "

!! NOW FRONTIER !!

        xfront=100.0d0 !! SET HIGH FOR FIRST STEP !!
        jjfrg=0
        jjorb=0
        do ifrg=1,icufr
          do ii=1,nb
            xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
            if(xcutoff+xthresh.lt.xx.and.xx.lt.xfront) then
              jjorb=ii
              jjfrg=ifrg
              xfront=xx
              xxfront=frgpop(jjfrg,jjorb)/deloc(jjfrg,jjorb)
            end if
          end do
        end do
        write(*,*) " 1st out ",jjorb," Frg. ",jjfrg," Q_A/Deloc. ",xxfront,
     +  " sqrt(deloc/Q_A) ",xfront
        write(*,*) " "
        write(*,*) " Frontier sqrt(deloc/Q_A) values (in,out) : ",xcutoff,xfront
        write(*,*) " "

!! 3) EVALUATING DEGENERACIES !!
!! infopop(i,j): SAVING THE FRAGMENT IN i = 1 AND ORBITAL NUMBER IN i = 2 !!
!! j ALLOCATED AS nb BUT MAXIMUM WILL BE inewcore !!

        ALLOCATE(infopop(2,nb))
        inewcore=0
        do ifrg=1,icufr
          do iorb=1,nb
            xx=xcutoff-dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
            if(ABS(xx).le.xthresh) then

!! APPLYING CONDITIONS TO REMOVE ORBITALS !!

              inewcore=inewcore+1
              infopop(1,inewcore)=ifrg
              infopop(2,inewcore)=iorb
              scrb(inewcore)=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              write(*,'(a29,i3,a12,i3)') " Orbital Equal/Over Cutoff:",infopop(2,inewcore),
     +        " Fragment : ",infopop(1,inewcore)
            end if
          end do
        end do
        write(*,*) " "
        write(*,'(a39,i3)') " NUMBER OF ORBITALS EQUAL/OVER CUTOFF: ",inewcore
        write(*,*) " "

!! SAVING THE FRAG OSLOs (CONSIDERED CORE) IN ccore !!

        do iorb=1,inewcore
          iifrg=infopop(1,iorb)
          iiorb=infopop(2,iorb)
          ifrgel(2,iifrg)=ifrgel(2,iifrg)+1 !! ADDING THEM HERE !!
          iaddcore=iaddcore+1
          iaddoslo=iaddoslo+1

!! SAVING IN INFOOSLO THE POPULATION AND FRAGMENT PREORTHOGONALIZATION !!

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

!! EVALUATING IF WE OVERASSIGN !!

        nnelect=0
        do ifrg=1,icufr
          nnelect=nnelect+ifrgel(2,ifrg)
        end do
        write(*,*) " DUMMY: iaddcore: ",iaddcore
        write(*,*) " Beta electron left to assign: ",nb-nnelect
        write(*,*) " "

        if(nb-nnelect.lt.0) then
          write(*,*) " !! OVERASSIGNING !!"
          write(*,*) " Overassignation by :",-(nb-nnelect)
!! lets trick the code !!
          write(*,*) "Tricking the code, removing last "
          inewcore=inewcore+(nb-nnelect)
          iaddoslo=iaddoslo+(nb-nnelect)
          nnelect=nnelect+(nb-nnelect)

        end if

!! SELECTING THE FIRST OUT FOR EVALUATING LINDEP !!

        do mu=1,igr
          do nu=1,igr
            clindep(mu,nu)=ZERO
          end do
        end do
        write(*,*) " "
        write(*,*) " SELECTING ORBITALS FOR LINDEP EVALUATION "
        write(*,*) " "
        iselected=0
        do ifrg=1,icufr
          do iorb=1,nb
            xx=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))-xfront

!! CRITERIA FOR SELECTION: sqrt(deloc/Q_A) <= frontier + xthresh !!

            if(xx.lt.xthresh) then
              iselected=iselected+1
              xx2=frgpop(ifrg,iorb)/deloc(ifrg,iorb)
              xx3=dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              write(*,*) " Orbital,Fragment,Q_A/deloc,sqrt(deloc/Q_A): ",iorb,ifrg,xx2,xx3
              do mu=1,igr
                clindep(mu,iselected)=cfrgoslo(ifrg,mu,iorb)
              end do
            end if
          end do
        end do
        write(*,*) " "
        write(*,*) " Selected orbitals for LinDep evaluation: ",iselected
        write(*,*) " "

        DEALLOCATE(deloc)
        DEALLOCATE(infopop)

!! EVALUATE LINDEP !!

        if(iselected.gt.0) then
          write(*,*) " LINIAR DEPENDENCY CHECK "
          write(*,*) " "
          do ii=1,nb
            do jj=1,nb
              SSS(ii,jj)=ZERO
              EEE(ii,jj)=ZERO
            end do
          end do
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

!! PRINTING ONLY FIRST iselected !!

          do ii=1,iselected
            write(*,*) " Orb., Eigenvalue : ",ii,SSS(ii,ii) 
          end do
          write(*,*) " "
          ilindep=0
          do ii=1,iselected
            if(SSS(ii,ii).lt.1.0d-5) ilindep=1
          end do
        end if

        if(ilindep.eq.1) then
          write(*,*) " WARNING : LINIAR DEPENDENCY FOUND "
          write(*,*) " sqrt(deloc/Q_A) values (and diff): ",xcutoff,xfront,xcutoff-xfront
          write(*,*) " TO DATE: Selecting largest one and continue "
          write(*,*) " "
        else
          write(*,*) " No lindep found in this step "
          write(*,*) " !! Iterative proceeds Safely !! "
          write(*,*) " "
        end if

!! ONLY IF NOT ALL ARE ASSIGNED !!
!! REMOVING ORBITALS FROM P MATRIX !!

        write(*,*) " "
        write(*,*) " Orbitals to remove in this step: ",inewcore
        write(*,*) " "

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

!! CHECKING BY RECONSTRUCTING PS !!

        write(*,*) " CHECK : COMPUTING trace(PS) "
        write(*,*) " "
        xx=ZERO
        xx2=ZERO
        do ii=1,igr
          do jj=1,igr
            xx=xx+pnocore(ii,jj)*s(jj,ii)
            do iat=1,nat
              xx2=xx2+pnocore(ii,jj)*sat(jj,ii,iat)
            end do
          end do
        end do
        write(*,'(a13,f10.5)') " Trace(PS) : ",xx 
        write(*,'(a24,f10.5)') " Trace(sum_nat(PSat)) : ",xx2 
        write(*,*) " " 

        if(nb-nnelect.eq.0) then
          write(*,*) " "
          write(*,*) "ALL BETA ELECTRONS ASSIGNED"
          write(*,*) " "
          go to 667
        end if

!! END DO iiter !!

      end do
667   continue

      write(*,*) " "
      write(*,*) " (BETA) END ITERATIVE PROCEDURE "
      write(*,*) " "

!! PRINTING OSLOs IN fchk !!
!! CONSTRUCTING poslo FOR .fchk CREATING !!

      ALLOCATE(poslo(igr,igr))
      do ii=1,igr
        do jj=1,igr
          xx=ZERO
          do ij=1,iaddoslo
            xx=xx+cboslo(ii,ij)*cboslo(jj,ij)
          end do
          poslo(ii,jj)=xx
        end do
      end do
      ctype="-beta-final-nonortho"
      call uwf_orbprint(1,0,cboslo,poslo,ctype)

!! PRINTING ORTHO OSLOs IN fchk !!
!! CONSTRUCTING poslo FOR .fchk CREATING !!

      do ii=1,igr
        do jj=1,igr
          xx=ZERO
          do ij=1,iaddoslo
            xx=xx+cbosloorth(ii,ij)*cbosloorth(jj,ij)
          end do
          poslo(ii,jj)=xx
        end do
      end do
      ctype="-beta-final-ORTHO"
      call uwf_orbprint(1,0,cbosloorth,poslo,ctype)
      DEALLOCATE(poslo)

!! OS ASSIGNMENT !!

      write(*,*) " "
      write(*,*) " **************************** "
      write(*,*) " "
      write(*,*) " FINAL EOS-like OS ASSIGNMENT "
      write(*,*) " "
      write(*,*) " **************************** "
      write(*,*) " "
      do ifrg=1,icufr
        write(*,'(a12,i3,i3)') " FRAG, OS : ",ifrg,iznfrg(ifrg)-(ifrgel(1,ifrg)+ifrgel(2,ifrg))
      end do
      write(*,*) " "

!! EVALUATING FINAL POPULATIONS TO COMPARE !!

      write(*,*) " "
      write(*,*) " ******************************** "
      write(*,*) " "
      write(*,*) " PRINTING FINAL ALPHA POPULATIONS "
      write(*,*) " "
      write(*,*) " ******************************** "
      write(*,*) " "
      do iorb=1,nalf
        iifrg=infooslo(1,iorb)
        write(*,*) " Frg., Spread, Pop. (pre), sqrt(deloc/Q_A): ",iifrg,fspread(iorb),poposlo(1,iorb),delocoslo(1,iorb)
      end do
      write(*,*) " "

      write(*,*) " "
      write(*,*) " ******************************* "
      write(*,*) " "
      write(*,*) " PRINTING FINAL BETA POPULATIONS "
      write(*,*) " "
      write(*,*) " ******************************* "
      write(*,*) " "
      do iorb=1,nb
        iifrg=infooslo(2,iorb)
        write(*,*) " Frg., Spread, Pop. (pre), sqrt(deloc/Q_A): ",iifrg,fbspread(iorb),poposlo(2,iorb),delocoslo(2,iorb)
      end do
      write(*,*) " "




!! DEALLOCATING !!

!     DEALLOCATE(pcore,pnocore)

      write(*,*) " "
      write(*,*) " ITERATIVE-NOCORE OSLO calculation COMPLETED "
      write(*,*) " "

      end
  
! ******
  

!! ****** !!

      subroutine rwf_frg_pop(ifrg,sat,corb,frgpop)

!! THIS SUBROUTINE COMPUTES FRAGMENT POPULATION ANALYSIS !!

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      dimension sat(igr,igr,nat)
      dimension corb(igr,igr),frgpop(nocc)

      allocatable :: orbpop(:,:)

!! LOADING IOPTs !!
      idofr=iopt(40)

!! ALLOCATING MATRICES !!
      ALLOCATE(orbpop(nocc,nat))

!! EVALUATING ORBITAL ATOMIC POPULATIONS !!
      do icenter=1,nat
        do iorb=1,nocc
          xx=ZERO
          do jj=1,igr
            do kk=1,igr
              xx=xx+corb(kk,iorb)*sat(kk,jj,icenter)*corb(jj,iorb)
            end do
          end do
          orbpop(iorb,icenter)=xx !! WITHOUT FACTOR OF 2, TO ADD OUTSIDE IF REQUIRED !!
        end do
      end do

!! GROUPING BY FRAGMENTS (ifrg) !!
      if(idofr.eq.1) then
        do iorb=1,nocc
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

      subroutine uwf_frg_pop2(iiss,ifrg,sat,corb,frgpop)

!! THIS SUBROUTINE COMPUTES FRAGMENT POPULATION ANALYSIS !!

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      allocatable orbpop(:,:)

      dimension sat(igr,igr,nat)
      dimension corb(igr,igr),frgpop(nocc)

      idofr=iopt(40)

!! IF ALPHA CASE !!

      if(iiss.eq.0) then

      ALLOCATE(orbpop(nalf,nat))
      do icenter=1,nat
        do iorb=1,nalf
          xx=ZERO
          do jj=1,igr
            do kk=1,igr
              xx=xx+corb(kk,iorb)*sat(kk,jj,icenter)*corb(jj,iorb)
            end do
          end do
          orbpop(iorb,icenter)=xx
        end do
      end do

!! FRAG POPULATIONS !!

      if(idofr.eq.1) then
        do iorb=1,nalf
          xx=ZERO
          do icenter=1,nfrlist(ifrg)
            xx=xx+orbpop(iorb,ifrlist(icenter,ifrg))
          end do
          frgpop(iorb)=xx
        end do
      end if
      DEALLOCATE(orbpop)

!! BETA CASE !!

      else if(iiss.eq.1) then

      ALLOCATE(orbpop(nb,nat))
      do icenter=1,nat
        do iorb=1,nb
          xx=ZERO
          do jj=1,igr
            do kk=1,igr
              xx=xx+corb(kk,iorb)*sat(kk,jj,icenter)*corb(jj,iorb)
            end do
          end do
          orbpop(iorb,icenter)=xx
        end do
      end do

!! FRAG POPULATIONS !!

      if(idofr.eq.1) then
        do iorb=1,nb
          xx=ZERO
          do icenter=1,nfrlist(ifrg)
            xx=xx+orbpop(iorb,ifrlist(icenter,ifrg))
          end do
          frgpop(iorb)=xx
        end do
      end if
      DEALLOCATE(orbpop)

      end if

      end

!! ****** !!

      subroutine uwf_frg_pop(iiss,ifrg,sat,corb,frgpop)

!! THIS SUBROUTINE COMPUTES FRAGMENT POPULATION ANALYSIS !!

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      allocatable orbpop(:,:)

      dimension sat(igr,igr,nat)
      dimension corb(igr,igr),frgpop(nocc)

      idofr=iopt(40)

      write(*,*) " "
      write(*,*) " PERFORMING FRAGMENT POPULATION ANALYSIS "
      write(*,*) " "

!! IF ALPHA CASE !!

      if(iiss.eq.0) then

      ALLOCATE(orbpop(nalf,nat))
      do icenter=1,nat
        do iorb=1,nalf
          xx=ZERO
          do jj=1,igr
            do kk=1,igr
              xx=xx+corb(kk,iorb)*sat(kk,jj,icenter)*corb(jj,iorb)
            end do
          end do
          orbpop(iorb,icenter)=xx
        end do
      end do
      write(*,*) " "

!! FRAG POPULATIONS !!

      if(idofr.eq.1) then
        do iorb=1,nalf
          xx=ZERO
          do icenter=1,nfrlist(ifrg)
            xx=xx+orbpop(iorb,ifrlist(icenter,ifrg))
          end do
          frgpop(iorb)=xx
          write(*,*) " Orbital, Fragment, Population : ",iorb,ifrg,xx
        end do
        write(*,*) " "
      end if
      DEALLOCATE(orbpop)

!! BETA CASE !!

      else if(iiss.eq.1) then

      ALLOCATE(orbpop(nb,nat))
      do icenter=1,nat
        do iorb=1,nb
          xx=ZERO
          do jj=1,igr
            do kk=1,igr
              xx=xx+corb(kk,iorb)*sat(kk,jj,icenter)*corb(jj,iorb)
            end do
          end do
          orbpop(iorb,icenter)=xx
        end do
      end do
      write(*,*) " "

!! FRAG POPULATIONS !!

      if(idofr.eq.1) then
        do iorb=1,nb
          xx=ZERO
          do icenter=1,nfrlist(ifrg)
            xx=xx+orbpop(iorb,ifrlist(icenter,ifrg))
          end do
          frgpop(iorb)=xx
          write(*,*) " Orbital, Fragment, Population : ",iorb,ifrg,xx
        end do
        write(*,*) " "
      end if
      DEALLOCATE(orbpop)

      end if

      end 

!! ****** !!

      subroutine rwf_orbprint(icount,cmat,pmat,ctype)

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
      character*5 ccent

      dimension cmat(igr,igr),pmat(igr,igr)

      iqchem=iopt(95)
      indepigr=int_locate(15,"Number of independ",ilog)
      norb=igr*indepigr
      norbt=igr*(igr+1)/2

!! ASSUMING UP TO 999 ATOMS !!

      if(icount.lt.10) then
        write(ccent,'(i1)') icount
        if(icount.eq.0) ccent="all"
      else if(icount.ge.10.and.icount.lt.100) then
        write(ccent,'(i2)') icount
      else
        write(ccent,'(i3)') icount
      end if

      name1=trim(name0)//trim(ctype)//trim(ccent)//".fchk"
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
      if(iqchem.eq.0) then
        do while(index(line,"Orthonormal basis").eq.0)
          read(15,'(a80)') line
        end do
      else
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

      subroutine uwf_orbprint(iiss,icount,cmat,pmat,ctype)

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
      character*5 ccent

      dimension cmat(igr,igr),pmat(igr,igr)

      iqchem=iopt(95)
      indepigr=int_locate(15,"Number of independ",ilog)
      norb=igr*indepigr
      norbt=igr*(igr+1)/2

!! ASSUMING UP TO 999 ATOMS !!

      if(icount.lt.10) then
        write(ccent,'(i1)') icount
        if(icount.eq.0) ccent="all"
      else if(icount.ge.10.and.icount.lt.100) then
        write(ccent,'(i2)') icount
      else
        write(ccent,'(i3)') icount
      end if

      name1=trim(name0)//trim(ctype)//trim(ccent)//".fchk"
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

