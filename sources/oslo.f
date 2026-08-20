
!! ********************************************************************** !!
!! OSLO calculation subroutines                                           !!
!! shared helpers (spin-agnostic, used by both drivers below):            !!
!!   oslo_build_Smat        -- fragment charge centers + position-        !!
!!                             operator matrices, computed once           !!
!!   oslo_density_from_coeffs -- P = factor*C*C^T from a set of OSLO      !!
!!                             coefficients                               !!
!!   oslo_channel_iterate   -- the iterative localize/assign/deflate      !!
!!                             scheme for one spin channel (RHF's only    !!
!!                             channel, or one of UHF's alpha/beta)       !!
!! restricted driver (prints in this order: .fchk output, per-OSLO       !!
!! summary tables, fragment oxidation states last):                      !!
!!   rwf_iterative_oslo     -- one oslo_channel_iterate call (nocc)       !!
!!   rwf_uwf_orbpop_atom    -- per-atom orbital population, the expensive !!
!!                             part of a fragment-population evaluation   !!
!!                             (shared with the unrestricted driver and   !!
!!                             LOBA too)                                  !!
!!   rwf_uwf_frg_pop        -- cheap per-fragment sum of the above        !!
!!                             (shared with the unrestricted driver too)  !!
!!   rwf_orbprint           -- restricted OSLO .fchk writer               !!
!!   rwf_uwf_print_OSLO_final -- per-OSLO summary table (shared with      !!
!!                             the unrestricted driver too)               !!
!! unrestricted driver (same print order as above, alpha then beta):      !!
!!   uwf_iterative_oslo     -- two oslo_channel_iterate calls (nalf, nb)  !!
!!   uwf_orbprint           -- unrestricted OSLO .fchk writer -- one      !!
!!                             file per stage with proper Alpha/Beta MO   !!
!!                             + Total/Spin SCF Density fields            !!
!! ********************************************************************** !!

!! ****** !!

!! ********************************************************************* !!
!! subroutine: oslo_build_Smat                                           !!
!! purpose: fragment-localization setup, shared by the restricted and    !!
!! unrestricted OSLO drivers -- for each fragment, computes its charge-  !!
!! center R_F (nuclear-charge-weighted average position of its own       !!
!! atoms) and the grid-integrated position-operator matrix               !!
!! sum_r w(r)*chi_mu(r)*chi_nu(r)*omega_A(r)*|r-R_F|^2 used by           !!
!! oslo_channel_iterate's per-fragment diagonalization. Independent of   !!
!! spin/occupation, so computed once and reused for every channel.       !!
!! arguments:                                                             !!
!!   itotps (in)  -- total number of grid points                         !!
!!   wp     (in)  -- integration weight of each grid point                !!
!!   omp2   (in)  -- atomic weight of each grid point, per atom           !!
!!   chp    (in)  -- basis-function values at each grid point             !!
!!   pcoord (in)  -- xyz coordinates of each grid point                   !!
!!   Smat   (out) -- (icufr,igr,igr) fragment position-operator matrices  !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine oslo_build_Smat(itotps,wp,omp2,chp,pcoord,Smat)

      use basis_set
      use ao_matrices
      use integration_grid
      implicit double precision(a-h,o-z)
      include 'parameter.h'

      integer,intent(in) :: itotps

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      dimension chp(itotps,igr),pcoord(itotps,3)
      dimension wp(itotps),omp2(itotps,nat)
      dimension Smat(icufr,igr,igr)

      iatps = nrad*nang

      write(*,*) " --------------------------------------- "
      write(*,*) "  CHARGE CENTER (R_F) FOR EACH FRAGMENT  "
      write(*,*) " --------------------------------------- "
      write(*,*) " "
      write(*,*) "  Frag.        Charge center (xyz)     "
      write(*,*) " ------------------------------------- "
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

!! position-operator matrix for this fragment, computed once. !!
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

110   FORMAT(3x,i3,3x,3f10.5)

      end

!! ****** !!

!! ********************************************************************* !!
!! subroutine: oslo_density_from_coeffs                                  !!
!! purpose: builds a density matrix from a set of OSLO coefficients,     !!
!! P = factor * C * C^T -- factor=2 for a doubly-occupied (restricted)   !!
!! channel, factor=1 for a singly-occupied (unrestricted alpha/beta)     !!
!! channel.                                                               !!
!! arguments:                                                             !!
!!   igr    (in)  -- number of basis functions                           !!
!!   norb   (in)  -- number of columns of cmat to include                !!
!!   cmat   (in)  -- (igr,igr) OSLO coefficients (only the first norb    !!
!!                   columns are read)                                   !!
!!   factor (in)  -- occupation factor (2 restricted, 1 unrestricted)     !!
!!   pmat   (out) -- (igr,igr) resulting density matrix                   !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine oslo_density_from_coeffs(igr,norb,cmat,factor,pmat)
      implicit double precision(a-h,o-z)
      include 'parameter.h'
      dimension cmat(igr,igr),pmat(igr,igr)

      do ii=1,igr
        do jj=1,igr
          xx=ZERO
          do ij=1,norb
            xx=xx+cmat(ii,ij)*cmat(jj,ij)
          end do
          pmat(ii,jj)=factor*xx
        end do
      end do

      end

!! ****** !!

!! ********************************************************************* !!
!! subroutine: oslo_channel_iterate                                      !!
!! purpose: iterative greedy OSLO assignment for a single spin channel   !!
!! -- the closed-shell "only channel", or one of UHF's alpha/beta        !!
!! channels. Each iteration diagonalizes every fragment's position-      !!
!! weighted density against the current remaining density (pnocore),    !!
!! scores every candidate orbital by its Pipek-Mezey-style               !!
!! delocalization index (FOLI), assigns the most-localized one (or a     !!
!! tolerance-based pack) to its fragment, Lowdin-orthogonalizes the      !!
!! newly-assigned orbitals, and deflates pnocore by their contribution   !!
!! before the next iteration. Terminates once nel orbitals are assigned. !!
!! Shared by rwf_iterative_oslo (nel=nocc, one call) and                 !!
!! uwf_iterative_oslo (nel=nalf then nel=nb, one call per channel) --    !!
!! extracted 2026-08-20 from what were three near-identical copies of    !!
!! this logic (RHF, UHF-alpha, UHF-beta); one had drifted and picked up  !!
!! a real print bug (pre-refactor state: git tag                        !!
!! oslo-pre-refactor-2026-08-20).                                        !!
!! arguments:                                                             !!
!!   nel      (in)    -- number of orbitals to assign in this channel     !!
!!                       (nocc for RHF, nalf or nb for UHF)               !!
!!   sat      (in)    -- (igr,igr,nat) per-atom AO overlap, for           !!
!!                       rwf_uwf_frg_pop's fragment population calls      !!
!!   Smat     (in)    -- (icufr,igr,igr) fragment position-operator       !!
!!                       matrices, built once by the caller, shared      !!
!!                       across channels                                 !!
!!   pnocore  (inout) -- (igr,igr) remaining density in this channel;     !!
!!                       caller initializes it (P for RHF, Pa/Pb for      !!
!!                       UHF), deflated in place as orbitals are assigned !!
!!   folitol  (in)    -- FOLI tolerance for orbital selection             !!
!!   ibranch  (in)    -- iteration at which to invoke branching (0 =      !!
!!                       disabled; branching itself is not implemented)  !!
!!   coslo    (out)   -- (igr,igr) pre-orthogonalization OSLO coeffs,     !!
!!                       columns 1..nel populated, rest zero (matches     !!
!!                       the .fchk writers' expected shape)               !!
!!   cosloorth(out)   -- (igr,igr) orthogonalized OSLO coeffs, same       !!
!!                       column convention as coslo                       !!
!!   delocoslo(out)   -- (nel) FOLI value of each OSLO, assignment order  !!
!!   ifrgel   (out)   -- (icufr) number of orbitals assigned per fragment !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine oslo_channel_iterate(nel,sat,Smat,pnocore,folitol,ibranch,
     &  coslo,cosloorth,delocoslo,ifrgel)

      use basis_set
      use ao_matrices
      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      integer,intent(in) :: ibranch
      real*8,intent(in) :: folitol

      dimension sat(igr,igr,nat)
      dimension Smat(icufr,igr,igr),pnocore(igr,igr)
      dimension coslo(igr,igr),cosloorth(igr,igr)
      dimension delocoslo(nel),ifrgel(icufr)

      allocatable :: S0(:,:),Sm(:,:),Splus(:,:),smh(:,:),eigv(:,:)
      allocatable :: c0(:,:),pp0(:,:),cmat(:,:),cfrgoslo(:,:,:)
      allocatable :: orbpop(:),orbpopat(:,:),deloc(:,:),clindep(:,:)
      allocatable :: infopop(:,:),scr(:)
      allocatable :: SSS(:,:),EEE(:,:)
      allocatable :: ccore(:,:),ccoreorth(:,:),pcore(:,:)
      allocatable :: frgpop(:,:),frgspr(:,:)

      niter=999

      ALLOCATE(SSS(nel,nel),EEE(nel,nel))
      ALLOCATE(Sm(igr,igr),Splus(igr,igr))
      ALLOCATE(c0(igr,igr),pp0(igr,igr))
      ALLOCATE(ccore(igr,igr),pcore(igr,igr),ccoreorth(igr,igr))
      ALLOCATE(clindep(igr,igr))
      ALLOCATE(scr(nel))
      ALLOCATE(cmat(igr,igr),cfrgoslo(icufr,igr,igr))
      ALLOCATE(frgpop(icufr,nel),frgspr(icufr,nel))

      ifrgel=0
      coslo=ZERO
      cosloorth=ZERO

      iaddoslo=0
      iaddoslo2=0
      do iiter=1,niter
        write(*,*) " ---------------------- "
        write(*,'(3x,a16,x,i3)') "ITERATION NUMBER",iiter
        write(*,*) " ---------------------- "
        write(*,*) " "

!! Zeroing the involved matrices !!
        pcore=ZERO
        ccore=ZERO
        ccoreorth=ZERO
        ALLOCATE(deloc(icufr,nel))
        deloc=ZERO

!! 1) Obtaining OSLOs for all frgs !!
        iaddcore=0
        ALLOCATE(S0(igr,igr))
        do ifrg=1,icufr
          do ii=1,igr
            do jj=1,igr

!! Using pp0 and S0 to not destroy pnocore and smat !!
              pp0(ii,jj)=pnocore(ii,jj)
              S0(ii,jj)=Smat(ifrg,ii,jj)
            end do
          end do
          call build_Smp(igr,S0,Sm,Splus,0)
          call to_lowdin_basis(igr,Splus,pp0)
          call diagonalize(igr,igr,pp0,c0,0)
          call to_AO_basis(igr,igr,Sm,c0)

!! Recovering the coeffs, saved in cfrgoslo !!
          do kk=1,nel
            do mu=1,igr
              cmat(mu,kk)=c0(mu,kk)*dsqrt(pp0(kk,kk))
              cfrgoslo(ifrg,mu,kk)=cmat(mu,kk)
            end do
          end do

!! Computing pipek delocalization, requires fragment populations. cmat  !!
!! is fixed for this ifrg, so the expensive per-atom population matrix  !!
!! is computed once here, then just summed per fragment below.         !!
          ALLOCATE(orbpop(nel),orbpopat(nel,nat))
          call rwf_uwf_orbpop_atom(sat,nel,cmat,orbpopat)
          do jfrg=1,icufr
            call rwf_uwf_frg_pop(jfrg,nel,orbpopat,orbpop)
            do ii=1,nel
              deloc(ifrg,ii)=deloc(ifrg,ii)+orbpop(ii)*orbpop(ii) !! Adding Q_A**2 inside deloc !!

!! Important here saving only for the own fragment !!
              if(jfrg.eq.ifrg) then
                frgspr(ifrg,ii)=pp0(ii,ii) !! Spreads (only for printing) !!
                frgpop(ifrg,ii)=orbpop(ii) !! Fragment populations !!
              end if
            end do
          end do
          DEALLOCATE(orbpop,orbpopat)

!! Now doing 1/deloc() !!
          do ii=1,nel
            if(deloc(ifrg,ii).gt.1.0d-6) then
              deloc(ifrg,ii)=ONE/deloc(ifrg,ii)
            else
              deloc(ifrg,ii)=100.0d0 !! Absurt value, avoids problems !!
            end if
          end do

!! Printing valuable information !!
          write(*,*) " -------------------------------------- "
          write(*,'(3x,a32,x,i3)') "ORBITAL INFORMATION FOR FRAGMENT",ifrg
          write(*,*) " -------------------------------------- "
          write(*,*) " "
          write(*,*) "  Orb.    Spread      Frg. Pop.      FOLI   "
          write(*,*) " ------------------------------------------ "
          do ii=1,nel
            xx=frgpop(ifrg,ii)
            if(xx.gt.1.0d-5) write(*,111) ii,frgspr(ifrg,ii),xx,dsqrt(deloc(ifrg,ii)/xx)
          end do
          write(*,*) " ------------------------------------------ "
          write(*,*) " "
        end do
        DEALLOCATE(S0)

!! 2) Cutoff evaluation !!
        xcutoff=100.0d0 !! Set high for first step !!
        iifrg=0
        iiorb=0
        do ifrg=1,icufr
          do ii=1,nel

!! Applying conditions to remove orbitals !!
            xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
            if(xx.lt.xcutoff) then
              iiorb=ii
              iifrg=ifrg
              xcutoff=xx !! New lowest FOLI !!
            end if
          end do
        end do

!! Now frontier (selection by packs using tolerance) !!
        xfront=100.0d0 !! Set high for first step !!
        jjfrg=0
        jjorb=0
        do ifrg=1,icufr
          do ii=1,nel
            xx=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
            if(xcutoff+folitol.lt.xx.and.xx.lt.xfront) then
              jjorb=ii
              jjfrg=ifrg
              xfront=xx
            end if
          end do
        end do

!! Printing !!
        write(*,'(2x,a27,x,i3,f10.5)') "Frg. and Lowest FOLI value:",iifrg,xcutoff
        write(*,'(2x,a56,x,i3,f10.5)') "Frg. and Lowest FOLI value including tolerance (cutoff):",jjfrg,xcutoff
        write(*,*) " "
        write(*,*) " ------------------- "
        write(*,*) "  SELECTED ORBITALS  "
        write(*,*) " ------------------- "
        write(*,*) " "

!! 3) Evaluating degeneracies !!
!! Infopop(i,j): saving the fragment in i = 1 and orbital number in i = 2 !!
!! J allocated as nel (for practicity) but maximum will be inewcore !!
        ALLOCATE(infopop(2,nel))
        inewcore=0

!! Branching (controlled from .inp, default = 0) !!
        if(iiter.eq.ibranch) then
          write(*,*) " ********************************************** "
          write(*,*) "  WARNING: BRANCHING INVOKED IN THIS ITERATION  "
          write(*,*) " ********************************************** "
          write(*,*) " "

!! MG: to do !!
          write(*,*) " Branching code has to be done "
          stop
        else
          write(*,*) "  Orb.   Frag.   FOLI  "
          write(*,*) " --------------------- "
          do ifrg=1,icufr
            do iorb=1,nel
              xx=xcutoff-dsqrt(deloc(ifrg,iorb)/frgpop(ifrg,iorb))
              if(ABS(xx).le.folitol) then

!! Applying conditions to remove orbitals !!
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

!! Saving the frag OSLOs (considered core) in ccore !!
        do iorb=1,inewcore
          iifrg=infopop(1,iorb)
          iiorb=infopop(2,iorb)
          ifrgel(iifrg)=ifrgel(iifrg)+1 !! Adding them here !!
          iaddcore=iaddcore+1
          iaddoslo=iaddoslo+1
          delocoslo(iaddoslo)=scr(iorb)
          do mu=1,igr
            ccore(mu,iaddcore)=cfrgoslo(iifrg,mu,iiorb)
            coslo(mu,iaddoslo)=cfrgoslo(iifrg,mu,iiorb)
          end do
        end do

!! Evaluating overassignment !!
        nnelect=0
        do ifrg=1,icufr
          nnelect=nnelect+ifrgel(ifrg)
        end do
        write(*,'(2x,a24,x,i3)') "Orbitals left to assign:",nel-nnelect
        write(*,*) " "

        if(nel-nnelect.lt.0) then
          write(*,*) " *************************************** "
          write(*,'(3x,a34,x,i3)') "WARNING: OVERASSIGNING BY (pairs):",-(nel-nnelect)
          write(*,*) " *************************************** "
          write(*,*) " "

!! Dirty trick !!
          write(*,*) " Continues by tricking the code (overassigned electrons removed) "
          write(*,*) " Check the final OSs, overassigned electrons have to be afterwards " !! To do !!
          write(*,*) " "
          inewcore=inewcore+(nel-nnelect)
          iaddoslo=iaddoslo+(nel-nnelect)
          nnelect=nnelect+(nel-nnelect)
        end if

!! Selecting the first out for evaluating LINDEP !!
        clindep=ZERO
        write(*,*) " ------------------------------ "
        write(*,*) "  CHECKING LINEAR DEPENDENCIES  "
        write(*,*) " ------------------------------ "
        write(*,*) " "
        write(*,*) "  Orb.   Frag.   FOLI  "
        write(*,*) " --------------------- "
        iselected=0
        do ifrg=1,icufr
          do iorb=1,nel
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

!! Last iteration no LINDEP evaluation !!
        iilindep=0
        if(nel-nnelect.eq.0.and.inewcore.eq.1) then
          write(*,*) " LinDep not evaluated in last iteration if only 1 orbital is selected "
          write(*,*) " "
          iilindep=1
        end if

!! Some deallocates... !!
        DEALLOCATE(deloc)
        DEALLOCATE(infopop)

!! Evaluate LINDEP !!
        ilindep=0
        if(iilindep.eq.0) then !! If only 1 there is nothing to evaluate !!
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

!! Diagonalizing all matrix, rest is zero so no affects !!
          call diagonalize(nel,nel,SSS,EEE,0)

!! Printing smallest eigenvalue !!
          xx=10.0d0 !! Absurt value again... !!
          do ii=1,iselected
            if(SSS(ii,ii).lt.xx) xx=SSS(ii,ii)
          end do
          if(xx.lt.1.0d-4) ilindep=1 !! Threshold for liniar dependency !!
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

!! Removing orbitals from P matrix, only if not all are assigned !!
!! Orthogonalizing fragment CORE/SEMICORE orbitals !!
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

!! Saving orthogonal orbitals here !!
        do ii=1,inewcore
          iaddoslo2=iaddoslo2+1
          do mu=1,igr
            cosloorth(mu,iaddoslo2)=ccoreorth(mu,ii)
          end do
        end do
        DEALLOCATE(smh)

!! Constructing pcore (and pnocore by substraction) !!
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

!! In case of all assigned !!
        if(nel-nnelect.eq.0) go to 666

!! End of iterative procedure !!
      end do
666   continue

      DEALLOCATE(SSS,EEE,Sm,Splus,c0,pp0,ccore,pcore,ccoreorth,clindep,scr,cmat,cfrgoslo)
      DEALLOCATE(frgpop,frgspr)

111   FORMAT(3x,i3,2x,f10.5,3x,f10.5,3x,f10.5)
112   FORMAT(3x,i3,4x,i3,x,f9.5)

      end

!! ****** !!

!! ********************************************************************* !!
!! subroutine: rwf_iterative_oslo                                        !!
!! purpose: OSLO oxidation-state assignment, restricted (closed-shell)   !!
!! case -- builds the fragment position-operator matrices                !!
!! (oslo_build_Smat), runs the iterative localize/assign/deflate scheme  !!
!! once for the single (doubly-occupied) channel (oslo_channel_iterate), !!
!! writes the pre-ortho/final OSLO .fchk files, prints the per-OSLO      !!
!! summary tables, then the fragment oxidation states last.              !!
!! MG: only the iterative procedure is implemented in this version; the  !!
!! non-iterative one lives in the development version.                   !!
!! arguments:                                                             !!
!!   sat    (in) -- (igr,igr,nat) per-atom AO overlap                     !!
!!   itotps (in) -- total number of grid points                          !!
!!   wp     (in) -- integration weight of each grid point                 !!
!!   omp2   (in) -- atomic weight of each grid point, per atom            !!
!!   chp    (in) -- basis-function values at each grid point              !!
!!   pcoord (in) -- xyz coordinates of each grid point                    !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine rwf_iterative_oslo(sat,itotps,wp,omp2,chp,pcoord)

      use basis_set
      use ao_matrices
      use integration_grid
      implicit double precision(a-h,o-z)
      include 'parameter.h'

      integer,intent(in) :: itotps

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(200)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      character*20 ctype

      dimension sat(igr,igr,nat)
      dimension chp(itotps,igr),pcoord(itotps,3)
      dimension wp(itotps),omp2(itotps,nat)

      allocatable :: Smat(:,:,:),pnocore(:,:),poslo(:,:)
      allocatable :: coslo(:,:),cosloorth(:,:),delocoslo(:)
      allocatable :: ifrgel(:),iznfrg(:)
      allocatable :: orbpop(:),orbpop2(:),foslo(:,:),foslo2(:,:)
      allocatable :: orbpopat(:,:),orbpopat2(:,:)

!! Loading iopts !!
      ifolitol = iopt(96)
      ibranch  = iopt(97)
      ifchk    = iopt(98)
      folitol = 10.0d0**(-REAL(ifolitol)) !! Default = 10^-3, controlled in .inp !!

!! Fragment charge extracted from zn (avoids problems when pseudopotentials are used) !!
      ALLOCATE(ifrgel(icufr),iznfrg(icufr))
      do ifrg=1,icufr
        izn=0
        do icenter=1,nfrlist(ifrg)
          iiat=ifrlist(icenter,ifrg)
          izn=izn+INT(zn(iiat))
        end do
        iznfrg(ifrg)=izn
      end do

      ALLOCATE(Smat(icufr,igr,igr))
      call oslo_build_Smat(itotps,wp,omp2,chp,pcoord,Smat)

!! Initial printing !!
      write(*,*) " "
      write(*,*) " ----------------------------------- "
      write(*,*) "  STARTING ITERATIVE OSLO ALGORITHM  "
      write(*,*) " ----------------------------------- "
      write(*,*) " "
      write(*,'(2x,a50,f10.5)') "Tolerance (in delta-FOLI) used for OSLO selection:",folitol
      write(*,*) " "

      ALLOCATE(pnocore(igr,igr))
      pnocore=pa !! SAVING pa IN pnocore !!

      ALLOCATE(coslo(igr,igr),cosloorth(igr,igr),delocoslo(nocc))

      call oslo_channel_iterate(nocc,sat,Smat,pnocore,folitol,ibranch,
     &  coslo,cosloorth,delocoslo,ifrgel)

      DEALLOCATE(Smat,pnocore)

!! Printing of the .fchk files with the OSLOs... to visualize !!
!! Preorthogonalization OSLOs can be visualized if desired (.inp) !!
      ALLOCATE(poslo(igr,igr)) !! Required poslo for .fchk creation !!
      if(ifchk.eq.2) then
        call oslo_density_from_coeffs(igr,nocc,coslo,TWO,poslo)
        ctype="-OSLOs-preortho"
        call rwf_orbprint(coslo,poslo,ctype)
      end if

!! Now the final (orthogonalized) ones !!
      call oslo_density_from_coeffs(igr,nocc,cosloorth,TWO,poslo)
      ctype="-OSLOs"
      call rwf_orbprint(cosloorth,poslo,ctype)
      DEALLOCATE(poslo)

!! Evaluating final populations to compare !!
      write(*,*) " ---------------------------------- "
      write(*,*) "  PRINTING FINAL OSLOs INFORMATION  "
      write(*,*) " ---------------------------------- "
      write(*,*) " "
      write(*,*) " ------------------------------------------- "
      write(*,*) "  Summary of the selected OSLOs (pre-ortho)  "
      write(*,*) " ------------------------------------------- "
      write(*,*) " "

!! Made A bit tricky... sorry !!
!! coslo/cosloorth are fixed here, so their expensive per-atom         !!
!! population matrices are each computed once, then just summed per   !!
!! fragment below.                                                     !!
      ALLOCATE(orbpop(nocc),orbpop2(nocc))
      ALLOCATE(orbpopat(nocc,nat),orbpopat2(nocc,nat))
      ALLOCATE(foslo(nocc,icufr),foslo2(nocc,icufr))
      foslo=ZERO
      foslo2=ZERO
      call rwf_uwf_orbpop_atom(sat,nocc,coslo,orbpopat) !! For the non-orthogonal OSLOs (original) !!
      call rwf_uwf_orbpop_atom(sat,nocc,cosloorth,orbpopat2) !! For the orthogonalized ones (printing later) !!
      do jfrg=1,icufr
        orbpop=ZERO
        orbpop2=ZERO
        call rwf_uwf_frg_pop(jfrg,nocc,orbpopat,orbpop)
        call rwf_uwf_frg_pop(jfrg,nocc,orbpopat2,orbpop2)
        do ii=1,nocc
          foslo(ii,jfrg)=orbpop(ii)
          foslo2(ii,jfrg)=orbpop2(ii)
        end do
      end do
      DEALLOCATE(orbpopat,orbpopat2)
      call rwf_uwf_print_OSLO_final(1,nocc,delocoslo,foslo)
      write(*,*) " --------------------------------------- "
      write(*,*) "  Summary of the selected OSLOs (final)  "
      write(*,*) " --------------------------------------- "
      write(*,*) " "
      call rwf_uwf_print_OSLO_final(0,nocc,delocoslo,foslo2) !! FOLI values given just for using the same routine !!

!! Final oxidation-state assignment, printed last -- each assigned      !!
!! orbital is doubly occupied.                                          !!
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

      DEALLOCATE(orbpop,orbpop2)
      DEALLOCATE(foslo,foslo2)
      DEALLOCATE(coslo,cosloorth,delocoslo,ifrgel,iznfrg)

20    FORMAT(3x,i3,6x,f8.2)

      end

!! ****** !!

!! ********************************************************************* !!
!! subroutine: rwf_uwf_orbpop_atom                                       !!
!! purpose: per-atom, per-orbital population (Mulliken-type, via sat)    !!
!! for a fixed set of orbital coefficients -- the expensive O(nat*norb*  !!
!! igr^2) part of the old rwf_uwf_frg_pop, split out so callers that     !!
!! need the same corb scored against every fragment (the usual case --  !!
!! see oslo_channel_iterate's and both drivers' do jfrg=1,icufr loops)   !!
!! compute it once instead of once per fragment.                        !!
!! arguments:                                                             !!
!!   sat    (in)  -- (igr,igr,nat) per-atom AO overlap                    !!
!!   norb   (in)  -- number of orbitals in corb to score                  !!
!!   corb   (in)  -- (igr,igr) candidate orbital coefficients (only the   !!
!!                   first norb columns are read)                        !!
!!   orbpop (out) -- (norb,nat) population of each orbital on each atom   !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine rwf_uwf_orbpop_atom(sat,norb,corb,orbpop)

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop

      dimension sat(igr,igr,nat)
      dimension corb(igr,igr),orbpop(norb,nat)

      do icenter=1,nat
        do iorb=1,norb
          xx=ZERO
          do jj=1,igr
            do kk=1,igr
              xx=xx+corb(kk,iorb)*sat(kk,jj,icenter)*corb(jj,iorb)
            end do
          end do
          orbpop(iorb,icenter)=xx
        end do
      end do

      end

!! ****** !!

!! ********************************************************************* !!
!! subroutine: rwf_uwf_frg_pop                                           !!
!! purpose: sums a precomputed per-atom orbital population (see          !!
!! rwf_uwf_orbpop_atom) onto one fragment -- independent of the          !!
!! wavefunction type, shared by the restricted and unrestricted OSLO     !!
!! drivers (and LOBA), called once per fragment per iteration to score   !!
!! how localized each candidate orbital is. No factor of 2 applied (the  !!
!! caller decides whether an orbital is singly- or doubly-occupied).     !!
!! frgpop is only written if DOFRAGS is set (idofr=1) -- always true in  !!
!! practice since OSLO is meaningless without fragments, but left        !!
!! unwritten (not zeroed) otherwise, matching this codebase's existing   !!
!! convention elsewhere for that flag.                                   !!
!! arguments:                                                             !!
!!   ifrg   (in)  -- fragment to compute the population on                !!
!!   norb   (in)  -- number of orbitals scored                            !!
!!   orbpop (in)  -- (norb,nat) per-atom orbital population, from         !!
!!                   rwf_uwf_orbpop_atom                                 !!
!!   frgpop (out) -- (norb) population of each orbital on fragment ifrg   !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine rwf_uwf_frg_pop(ifrg,norb,orbpop,frgpop)

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      dimension orbpop(norb,nat),frgpop(norb)

      idofr=iopt(40)

      if(idofr.eq.1) then
        do iorb=1,norb
          xx=ZERO
          do icenter=1,nfrlist(ifrg)
            xx=xx+orbpop(iorb,ifrlist(icenter,ifrg))
          end do
          frgpop(iorb)=xx
        end do
      end if

      end

!! ****** !!

!! ********************************************************************* !!
!! subroutine: rwf_orbprint                                              !!
!! purpose: writes the restricted OSLO .fchk -- splices the original     !!
!! wavefunction .fchk's structure, replacing "Alpha MO coefficients"     !!
!! and "Total SCF Density" with the OSLO coefficients/density, copying   !!
!! everything else through unchanged. See uwf_orbprint for the           !!
!! unrestricted (alpha+beta combined) twin.                              !!
!! arguments:                                                             !!
!!   cmat  (in) -- (igr,igr) OSLO coefficients (columns 1..nocc)          !!
!!   pmat  (in) -- (igr,igr) OSLO density (2*C*C^T, doubly-occupied)      !!
!!   ctype (in) -- filename suffix, e.g. "-OSLOs" or "-OSLOs-preortho"    !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine rwf_orbprint(cmat,pmat,ctype)

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
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

!! Name of the .fchk file !!
      name1=trim(name0)//trim(ctype)//".fchk"
      open(unit=69,file=name1)
      rewind(69)
      rewind(15)

!! Printing until alpha MOs !!
      read(15,'(a80)') line
      do while(index(line,"Alpha MO co").eq.0)
        write(69,'(a80)') line
        read(15,'(a80)') line
      end do

!! Printing the new ones !!
      write(69,11) "Alpha MO coefficients","R","N= ",norb
      write(69,13) ((cmat(ii,jj),ii=1,igr),jj=1,indepigr)

!! Now locating what is after it in the original one to continue !!
      if(iqchem.eq.0.and.imokit.eq.0) then
        do while(index(line,"Orthonormal basis").eq.0)
          read(15,'(a80)') line
        end do
      else if(iqchem.eq.1) then
        do while(index(line,"Alpha Orbital").eq.0)
          read(15,'(a80)') line
        end do
      end if

!! Restart printing until next stop !!
      do while(index(line,"Total SCF Dens").eq.0)
        write(69,'(a80)') line
        read(15,'(a80)') line
      end do

!! Printing the new one !!
      write(69,12) "Total SCF Density","R","N= ",norbt
      write(69,13) ((pmat(ii,jj),jj=1,ii),ii=1,igr)

!! Now locating what is after it in the original one to continue !!
      if(iqchem.eq.0.and.imokit.eq.0) then
        do while(index(line,"Mulliken Charges").eq.0)
          read(15,'(a80)') line
        end do
      else if(iqchem.eq.1) then
        do while(index(line,"Pure Switching").eq.0)
          read(15,'(a80)') line
        end do
      end if

!! Restart printing until the end !!
      do while(.true.)
        write(69,'(a80)') line
        read(15,'(a80)',end=99) line
      end do
99    continue
      close(69)

!! Printing formats !!
11    FORMAT(a21,22x,a1,3x,a3,i11)
12    FORMAT(a17,26x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))

      end

!! ****** !!

!! ********************************************************************* !!
!! subroutine: rwf_uwf_print_OSLO_final                                  !!
!! purpose: prints the per-OSLO fragment-population (and, if requested,  !!
!! FOLI) summary table in packs of 5 columns -- shared by the restricted !!
!! and unrestricted OSLO drivers, called once per stage (pre-ortho/      !!
!! final) and, for the unrestricted case, once per spin channel.         !!
!! arguments:                                                             !!
!!   iflag  (in) -- 1 to also print the FOLI Value row, 0 to omit it      !!
!!   noslo  (in) -- number of OSLOs (columns) to print                    !!
!!   foli   (in) -- (igr) FOLI value of each OSLO (only read if iflag=1)  !!
!!   frgpop (in) -- (noslo,icufr) fragment population of each OSLO        !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine rwf_uwf_print_OSLO_final(iflag,noslo,foli,frgpop)

      implicit double precision (a-h,o-z)
      include 'parameter.h'

      integer :: b,dd

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      dimension foli(igr)
      dimension frgpop(noslo,icufr)

!! number of 5-column packs needed, rounded up. !!
      b=noslo/5
      if(b*5.ne.noslo) b=(noslo/5)+1

      dd=1
      do k=1,b

!! last (possibly partial) pack. !!
        if(k.eq.b) then
          write(*,'(2x,a13,5(i7,3x))') "OSLO Number :",(jj,jj=dd,noslo)
          if(iflag.eq.1) write(*,'(2x,a13,5f10.5)') "FOLI Value  :",(foli(jj),jj=dd,noslo)
          do jfrg=1,icufr
            write(*,'(2x,a9,i3,a1,5f10.5)') "Frg. Pop.",jfrg,":",(frgpop(jj,jfrg),jj=dd,noslo)
          end do
          write(*,*) " "

!! full 5-column packs. !!
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

!! ********************************************************************* !!
!! subroutine: uwf_iterative_oslo                                        !!
!! purpose: OSLO oxidation-state assignment, unrestricted (open-shell)   !!
!! case -- builds the fragment position-operator matrices                !!
!! (oslo_build_Smat, shared between spins), runs the iterative           !!
!! localize/assign/deflate scheme once per spin channel                  !!
!! (oslo_channel_iterate, alpha then beta), writes the combined (alpha+  !!
!! beta together) pre-ortho/final OSLO .fchk files, prints the per-OSLO  !!
!! summary tables for each spin, then the fragment oxidation states      !!
!! last.                                                                  !!
!! MG: only the iterative procedure is implemented in this version; the  !!
!! non-iterative one lives in the development version.                   !!
!! arguments:                                                             !!
!!   sat    (in) -- (igr,igr,nat) per-atom AO overlap                     !!
!!   itotps (in) -- total number of grid points                          !!
!!   wp     (in) -- integration weight of each grid point                 !!
!!   omp2   (in) -- atomic weight of each grid point, per atom            !!
!!   chp    (in) -- basis-function values at each grid point              !!
!!   pcoord (in) -- xyz coordinates of each grid point                    !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine uwf_iterative_oslo(sat,itotps,wp,omp2,chp,pcoord)

      use basis_set
      use ao_matrices
      use integration_grid
      implicit double precision(a-h,o-z)
      include 'parameter.h'

      integer,intent(in) :: itotps

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(200)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      character*20 ctype

      dimension sat(igr,igr,nat)
      dimension chp(itotps,igr),pcoord(itotps,3)
      dimension wp(itotps),omp2(itotps,nat)

      allocatable :: Smat(:,:,:),pnocore(:,:)
      allocatable :: coslo_a(:,:),cosloorth_a(:,:),delocoslo_a(:)
      allocatable :: coslo_b(:,:),cosloorth_b(:,:),delocoslo_b(:)
      allocatable :: ifrgel_a(:),ifrgel_b(:),iznfrg(:)
      allocatable :: poslo_a(:,:),poslo_b(:,:)
      allocatable :: orbpop(:),orbpop2(:),foslo(:,:),foslo2(:,:)
      allocatable :: orbpopat(:,:),orbpopat2(:,:)

!! Loading iopts !!
      ifolitol = iopt(96)
      ibranch  = iopt(97)
      ifchk    = iopt(98)
      folitol = 10.0d0**(-REAL(ifolitol)) !! Default = 10^-3, controlled in .inp !!

!! Fragment charge extracted from zn (avoids problems when pseudopotentials are used) !!
      ALLOCATE(ifrgel_a(icufr),ifrgel_b(icufr),iznfrg(icufr))
      do ifrg=1,icufr
        izn=0
        do icenter=1,nfrlist(ifrg)
          iiat=ifrlist(icenter,ifrg)
          izn=izn+INT(zn(iiat))
        end do
        iznfrg(ifrg)=izn
      end do

      ALLOCATE(Smat(icufr,igr,igr))
      call oslo_build_Smat(itotps,wp,omp2,chp,pcoord,Smat)

!! Initial printing !!
      write(*,*) " "
      write(*,*) " ----------------------------------- "
      write(*,*) "  STARTING ITERATIVE OSLO ALGORITHM  "
      write(*,*) " ----------------------------------- "
      write(*,*) " "
      write(*,'(2x,a50,f10.5)') "Tolerance (in delta-FOLI) used for OSLO selection:",folitol
      write(*,*) " "

!! Alpha channel !!
      write(*,*) " ------------ "
      write(*,*) "  ALPHA PART  "
      write(*,*) " ------------ "
      write(*,*) " "

      ALLOCATE(pnocore(igr,igr))
      pnocore=pa
      ALLOCATE(coslo_a(igr,igr),cosloorth_a(igr,igr),delocoslo_a(nalf))
      call oslo_channel_iterate(nalf,sat,Smat,pnocore,folitol,ibranch,
     &  coslo_a,cosloorth_a,delocoslo_a,ifrgel_a)
      DEALLOCATE(pnocore)

!! Beta channel !!
      write(*,*) " ----------- "
      write(*,*) "  BETA PART  "
      write(*,*) " ----------- "
      write(*,*) " "

      ALLOCATE(pnocore(igr,igr))
      pnocore=pb
      ALLOCATE(coslo_b(igr,igr),cosloorth_b(igr,igr),delocoslo_b(nb))
      call oslo_channel_iterate(nb,sat,Smat,pnocore,folitol,ibranch,
     &  coslo_b,cosloorth_b,delocoslo_b,ifrgel_b)
      DEALLOCATE(pnocore,Smat)

!! Printing of the combined .fchk files with the OSLOs (alpha+beta together) !!
!! Preorthogonalization OSLOs can be visualized if desired (.inp) !!
      ALLOCATE(poslo_a(igr,igr),poslo_b(igr,igr))
      if(ifchk.eq.2) then
        call oslo_density_from_coeffs(igr,nalf,coslo_a,ONE,poslo_a)
        call oslo_density_from_coeffs(igr,nb,coslo_b,ONE,poslo_b)
        ctype="-OSLOs-preortho"
        call uwf_orbprint(coslo_a,coslo_b,poslo_a,poslo_b,ctype)
      end if

!! Now the final (orthogonalized) ones !!
      call oslo_density_from_coeffs(igr,nalf,cosloorth_a,ONE,poslo_a)
      call oslo_density_from_coeffs(igr,nb,cosloorth_b,ONE,poslo_b)
      ctype="-OSLOs"
      call uwf_orbprint(cosloorth_a,cosloorth_b,poslo_a,poslo_b,ctype)
      DEALLOCATE(poslo_a,poslo_b)

!! Evaluating final populations to compare !!
      write(*,*) " ---------------------------------- "
      write(*,*) "  PRINTING FINAL OSLOs INFORMATION  "
      write(*,*) " ---------------------------------- "
      write(*,*) " "

!! First alpha !!
      write(*,*) " ------------------------------------------------- "
      write(*,*) "  Summary of the selected alpha OSLOs (pre-ortho)  "
      write(*,*) " ------------------------------------------------- "
      write(*,*) " "
!! coslo_a/cosloorth_a are fixed here, so their expensive per-atom     !!
!! population matrices are each computed once, then just summed per   !!
!! fragment below.                                                     !!
      ALLOCATE(orbpop(nalf),orbpop2(nalf))
      ALLOCATE(orbpopat(nalf,nat),orbpopat2(nalf,nat))
      ALLOCATE(foslo(nalf,icufr),foslo2(nalf,icufr))
      foslo=ZERO
      foslo2=ZERO
      call rwf_uwf_orbpop_atom(sat,nalf,coslo_a,orbpopat) !! For the non-orthogonal OSLOs (original) !!
      call rwf_uwf_orbpop_atom(sat,nalf,cosloorth_a,orbpopat2) !! For the orthogonalized ones (printing later) !!
      do jfrg=1,icufr
        orbpop=ZERO
        orbpop2=ZERO
        call rwf_uwf_frg_pop(jfrg,nalf,orbpopat,orbpop)
        call rwf_uwf_frg_pop(jfrg,nalf,orbpopat2,orbpop2)
        do ii=1,nalf
          foslo(ii,jfrg)=orbpop(ii)
          foslo2(ii,jfrg)=orbpop2(ii)
        end do
      end do
      DEALLOCATE(orbpopat,orbpopat2)
      call rwf_uwf_print_OSLO_final(1,nalf,delocoslo_a,foslo)
      write(*,*) " --------------------------------------------- "
      write(*,*) "  Summary of the selected alpha OSLOs (final)  "
      write(*,*) " --------------------------------------------- "
      write(*,*) " "
      call rwf_uwf_print_OSLO_final(0,nalf,delocoslo_a,foslo2) !! FOLI values given just for using same routine !!
      DEALLOCATE(orbpop,orbpop2)
      DEALLOCATE(foslo,foslo2)

!! Now beta !!
      write(*,*) " ------------------------------------------------ "
      write(*,*) "  Summary of the selected beta OSLOs (pre-ortho)  "
      write(*,*) " ------------------------------------------------ "
      write(*,*) " "
!! coslo_b/cosloorth_b are fixed here, so their expensive per-atom     !!
!! population matrices are each computed once, then just summed per   !!
!! fragment below.                                                     !!
      ALLOCATE(orbpop(nb),orbpop2(nb))
      ALLOCATE(orbpopat(nb,nat),orbpopat2(nb,nat))
      ALLOCATE(foslo(nb,icufr),foslo2(nb,icufr))
      foslo=ZERO
      foslo2=ZERO
      call rwf_uwf_orbpop_atom(sat,nb,coslo_b,orbpopat) !! For the non-orthogonal OSLOs (original) !!
      call rwf_uwf_orbpop_atom(sat,nb,cosloorth_b,orbpopat2) !! For the orthogonalized ones (printing later) !!
      do jfrg=1,icufr
        orbpop=ZERO
        orbpop2=ZERO
        call rwf_uwf_frg_pop(jfrg,nb,orbpopat,orbpop)
        call rwf_uwf_frg_pop(jfrg,nb,orbpopat2,orbpop2)
        do ii=1,nb
          foslo(ii,jfrg)=orbpop(ii)
          foslo2(ii,jfrg)=orbpop2(ii)
        end do
      end do
      DEALLOCATE(orbpopat,orbpopat2)
      call rwf_uwf_print_OSLO_final(1,nb,delocoslo_b,foslo)
      write(*,*) " -------------------------------------------- "
      write(*,*) "  Summary of the selected beta OSLOs (final)  "
      write(*,*) " -------------------------------------------- "
      write(*,*) " "
      call rwf_uwf_print_OSLO_final(0,nb,delocoslo_b,foslo2) !! FOLI values given just for using the same routine !!
      DEALLOCATE(orbpop,orbpop2)
      DEALLOCATE(foslo,foslo2)

!! Final oxidation-state assignment, printed last -- each assigned      !!
!! spin-orbital holds 1 electron.                                       !!
      write(*,*) " --------------------------- "
      write(*,*) "  FRAGMENT OXIDATION STATES  "
      write(*,*) " --------------------------- "
      write(*,*) " "
      write(*,*) "  Frag.  Oxidation State  "
      write(*,*) " ------------------------ "
      do ifrg=1,icufr
        write(*,20) ifrg,REAL(iznfrg(ifrg)-(ifrgel_a(ifrg)+ifrgel_b(ifrg)))
      end do
      write(*,*) " ------------------------ "
      write(*,*) " "

      DEALLOCATE(coslo_a,cosloorth_a,delocoslo_a,ifrgel_a)
      DEALLOCATE(coslo_b,cosloorth_b,delocoslo_b,ifrgel_b)
      DEALLOCATE(iznfrg)

20    FORMAT(3x,i3,6x,f8.2)

      end

!! ****** !!

!! ********************************************************************* !!
!! subroutine: uwf_orbprint                                              !!
!! purpose: writes a single combined .fchk with the unrestricted OSLO    !!
!! orbitals -- splices the original wavefunction .fchk's structure,      !!
!! replacing "Alpha MO coefficients"/"Beta MO coefficients" with the     !!
!! OSLO coefficients for each spin, and "Total SCF Density"/"Spin SCF    !!
!! Density" with Pa_oslo+Pb_oslo / Pa_oslo-Pb_oslo -- everything else    !!
!! copied through unchanged. The Spin SCF Density block is only written  !!
!! if the source .fchk has one to begin with. Replaces the previous two  !!
!! separate per-spin files (each mislabeling its own single-spin density !!
!! as "Total SCF Density", with no Spin SCF Density at all) with one     !!
!! properly-labeled file -- see git tag oslo-pre-refactor-2026-08-20 for !!
!! the old behavior.                                                     !!
!! arguments:                                                             !!
!!   cmat_a (in) -- (igr,igr) alpha OSLO coefficients (columns 1..nalf)   !!
!!   cmat_b (in) -- (igr,igr) beta OSLO coefficients (columns 1..nb)      !!
!!   pmat_a (in) -- (igr,igr) alpha OSLO density (Ca*Ca^T)                !!
!!   pmat_b (in) -- (igr,igr) beta OSLO density (Cb*Cb^T)                 !!
!!   ctype  (in) -- filename suffix, e.g. "-OSLOs" or "-OSLOs-preortho"   !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
      subroutine uwf_orbprint(cmat_a,cmat_b,pmat_a,pmat_b,ctype)

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /filename/name0

      character*80 line
      character*60 name0,name1
      character*20 ctype

      dimension cmat_a(igr,igr),cmat_b(igr,igr)
      dimension pmat_a(igr,igr),pmat_b(igr,igr)

      iqchem   = iopt(95)
      imokit   = iopt(79)
      indepigr = int_locate(15,"Number of independ",ilog)
      norb     = igr*indepigr
      norbt    = igr*(igr+1)/2

!! whether the source .fchk carries a Spin SCF Density block to replace !!
!! -- a plain presence check, doesn't disturb the real splice pass below !!
!! since the file gets rewound before that starts.                       !!
      call locate(15,"Spin SCF Dens",ihasspin)

!! Name of the .fchk file !!
      name1=trim(name0)//trim(ctype)//".fchk"
      open(unit=69,file=name1)
      rewind(69)
      rewind(15)

!! Printing until alpha MOs !!
      read(15,'(a80)') line
      do while(index(line,"Alpha MO co").eq.0)
        write(69,'(a80)') line
        read(15,'(a80)') line
      end do

!! Printing the new alpha block !!
      write(69,11) "Alpha MO coefficients","R","N= ",norb
      write(69,13) ((cmat_a(ii,jj),ii=1,igr),jj=1,indepigr)

!! Skipping the original alpha mo data, up to the beta marker !!
      do while(index(line,"Beta MO coef").eq.0)
        read(15,'(a80)') line
      end do

!! Printing the new beta block !!
      write(69,11) "Beta MO coefficients ","R","N= ",norb
      write(69,13) ((cmat_b(ii,jj),ii=1,igr),jj=1,indepigr)

!! Now locating what is after it in the original one to continue !!
      if(iqchem.eq.0.and.imokit.eq.0) then
        do while(index(line,"Orthonormal basis").eq.0)
          read(15,'(a80)') line
        end do
      else if(iqchem.eq.1) then
        do while(index(line,"Alpha Orbital").eq.0)
          read(15,'(a80)') line
        end do
      end if

!! Restart printing until total scf density !!
      do while(index(line,"Total SCF Dens").eq.0)
        write(69,'(a80)') line
        read(15,'(a80)') line
      end do

!! Printing the new total SCF density (Pa+Pb) !!
      write(69,12) "Total SCF Density","R","N= ",norbt
      write(69,13) ((pmat_a(ii,jj)+pmat_b(ii,jj),jj=1,ii),ii=1,igr)

      if(ihasspin.eq.1) then

!! Skipping the original total density data, up to the spin marker !!
        do while(index(line,"Spin SCF Dens").eq.0)
          read(15,'(a80)') line
        end do

!! Printing the new spin SCF density (Pa-Pb) !!
        write(69,12) "Spin SCF Density ","R","N= ",norbt
        write(69,13) ((pmat_a(ii,jj)-pmat_b(ii,jj),jj=1,ii),ii=1,igr)
      end if

!! Now locating what is after it in the original one to continue !!
      if(iqchem.eq.0.and.imokit.eq.0) then
        do while(index(line,"Mulliken Charges").eq.0)
          read(15,'(a80)') line
        end do
      else if(iqchem.eq.1) then
        do while(index(line,"Pure Switching").eq.0)
          read(15,'(a80)') line
        end do
      end if

!! Restart printing until the end !!
      do while(.true.)
        write(69,'(a80)') line
        read(15,'(a80)',end=99) line
      end do
99    continue
      close(69)

!! Printing formats !!
11    FORMAT(a21,22x,a1,3x,a3,i11)
12    FORMAT(a17,26x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))

      end

!! ****** !!
