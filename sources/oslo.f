
!! ----------------------------------------------------------------------
!! OSLO: oxidation states from localized orbitals.
!!   oslo_build_Smat         -- fragment charge centres and spread matrices
!!   oslo_density_from_coeffs -- P = factor*C*C^T of a set of OSLOs
!!   oslo_channel_iterate    -- the iterative selection for one channel
!!   rwf_iterative_oslo, uwf_iterative_oslo -- restricted / unrestricted
!!                              drivers
!!   oslo_header, oslo_summary, oslo_oxidation_states -- printing shared
!!                              by both drivers
!!   oslo_fchk_write         -- the OSLO .fchk files (FOLI as orbital
!!                              energies)
!!   oslo_json_scheme, oslo_json_file, json_oslo_channel -- results file
!!   rwf_uwf_orbpop_atom, rwf_uwf_frg_pop -- orbital populations per atom
!!                              and per fragment (also used by LOBA)
!! ----------------------------------------------------------------------

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: oslo_build_Smat
!! Charge centre R_F of each fragment (its atoms' positions weighted by
!! atomic number), printed in Angstrom, and its spread matrix, the grid
!! integral of chi_mu chi_nu |r-R_F|^2 weighted by the atomic weights,
!! used by every iteration of every channel.
!! arguments:
!!   itotps (in)   number of grid points
!!   wp     (in)   (itotps) integration weights
!!   omp2   (in)   (itotps,nat) atomic weights of each point
!!   chp    (in)   (itotps,igr) basis functions at the points
!!   pcoord (in)   (itotps,3) coordinates of the points
!!   Smat   (out)  (icufr,igr,igr) spread matrix of each fragment
!! author: MGimf
!! ----------------------------------------------------------------------
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

      write(*,'(a)') ''
      write(*,'(2x,a)') 'Charge centre R_F of each fragment (Angstrom)'
      write(*,'(2x,a)') 'Frag.         x         y         z'
      write(*,'(2x,a)') repeat('-',35)
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
      write(*,'(2x,a)') repeat('-',35)

110   FORMAT(2x,i5,3f10.5)

      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: oslo_density_from_coeffs
!! Density P = factor * C C^T of the first norb columns of a set of
!! OSLOs (factor 2 for a restricted channel, 1 for alpha or beta).
!! arguments:
!!   igr    (in)   number of basis functions
!!   norb   (in)   number of OSLOs
!!   cmat   (in)   (igr,igr) their coefficients
!!   factor (in)   electrons per orbital
!!   pmat   (out)  (igr,igr) density
!! author: MGimf
!! ----------------------------------------------------------------------
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

!! ----------------------------------------------------------------------
!! subroutine: oslo_channel_iterate
!! Iterative OSLO assignment of one channel (the doubly occupied
!! orbitals of a restricted wavefunction, or alpha or beta). Each
!! iteration localizes the remaining density around every fragment,
!! scores the candidates by their FOLI, selects the lowest one with
!! every candidate within the tolerance (or the next group, when
!! branching), orthogonalizes the selection and removes it from the
!! remaining density, until nel orbitals are assigned. A selection that
!! spans fewer orbitals than it holds (the same orbital found from
!! several fragments) gives them to its fragments in shares.
!! arguments:
!!   nel       (in)     number of orbitals to assign
!!   sat       (in)     (igr,igr,nat) atomic overlap matrices
!!   Smat      (in)     (icufr,igr,igr) spread matrices (oslo_build_Smat)
!!   pnocore   (inout)  (igr,igr) remaining density, deflated in place
!!   folitol   (in)     FOLI tolerance of a selection
!!   nbr       (in)     number of branching iterations
!!   ibr       (in)     (nbr) iterations that select the next FOLI group
!!   chkey     (in)     input keyword that branches this channel
!!   coslo     (out)    (igr,igr) OSLOs before orthogonalization
!!   cosloorth (out)    (igr,igr) orthonormal OSLOs (the .fchk orbitals)
!!   delocoslo (out)    (nel) FOLI value of each OSLO
!!   frgel     (out)    (icufr) orbitals assigned to each fragment
!!   iterosl   (out)    (nel) iteration that selected each OSLO
!!   shrosl    (out)    (icufr,nel) share of each fragment in each OSLO
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine oslo_channel_iterate(nel,sat,Smat,pnocore,folitol,nbr,
     &  ibr,chkey,coslo,cosloorth,delocoslo,frgel,iterosl,shrosl)

      use basis_set
      use ao_matrices
      use output_mod, only: print_value
      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      integer,intent(in) :: nbr,ibr(*)
      integer,intent(out) :: iterosl(nel)
      real*8,intent(out) :: shrosl(icufr,nel)
      real*8,intent(in) :: folitol
      character*(*),intent(in) :: chkey

!! a next OSLO whose overlap with the selected ones is above 1-xlindep  !!
!! is nearly one of them; with a delta-FOLI below dfoliwarn it is a      !!
!! close alternative selection                                           !!
      real*8,parameter :: xlindep=0.1d0,dfoliwarn=0.5d0

      character*40 line
      logical lbranch

      dimension sat(igr,igr,nat)
      dimension Smat(icufr,igr,igr),pnocore(igr,igr)
      dimension coslo(igr,igr),cosloorth(igr,igr)
      dimension delocoslo(nel),frgel(icufr)

      allocatable :: S0(:,:),Sm(:,:),Splus(:,:),smh(:,:),eigv(:,:)
      allocatable :: c0(:,:),pp0(:,:),cmat(:,:),cfrgoslo(:,:,:)
      allocatable :: orbpop(:),orbpopat(:,:),deloc(:,:),foli(:,:)
      allocatable :: infopop(:,:),sfoli(:),share(:),st(:)
      allocatable :: ccore(:,:),ccoreorth(:,:),pcore(:,:)
      allocatable :: frgpop(:,:),frgspr(:,:),ibrused(:)

      frgel=ZERO
      coslo=ZERO
      cosloorth=ZERO
      iterosl=0
      shrosl=ZERO

!! a channel with no electrons (e.g. beta of a high-spin H2 triplet) has  !!
!! nothing to localize                                                   !!
      if(nel.eq.0) then
        write(*,'(2x,a)') 'No electrons in this channel: nothing to '//
     +    'localize'
        do k=1,nbr
          write(*,'(2x,a,i0,a)') 'Branching at iteration ',ibr(k),
     +      ': not applied (no electrons in this channel)'
        end do
        return
      end if

      niter=999

      ALLOCATE(Sm(igr,igr),Splus(igr,igr))
      ALLOCATE(c0(igr,igr),pp0(igr,igr))
      ALLOCATE(pcore(igr,igr),st(igr))
      ALLOCATE(cmat(igr,igr),cfrgoslo(icufr,igr,igr))
      ALLOCATE(frgpop(icufr,nel),frgspr(icufr,nel))
      ALLOCATE(ibrused(max(1,nbr)))
      ibrused=0

      nnelect=0
      iadd=0
      do iiter=1,niter
        write(line,'(a,i0)') "ITERATION ",iiter
        call print_box(trim(line))

        pcore=ZERO
        ALLOCATE(deloc(icufr,nel),foli(icufr,nel))
        deloc=ZERO

!! 1) Obtaining OSLOs for all frgs !!
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

!! Now doing 1/deloc(), then the FOLI of each candidate !!
          do ii=1,nel
            if(deloc(ifrg,ii).gt.1.0d-6) then
              deloc(ifrg,ii)=ONE/deloc(ifrg,ii)
            else
              deloc(ifrg,ii)=100.0d0 !! Absurt value, avoids problems !!
            end if
            foli(ifrg,ii)=dsqrt(deloc(ifrg,ii)/frgpop(ifrg,ii))
          end do

!! candidate orbitals of this fragment !!
          if(ifrg.gt.1) write(*,'(a)') ''
          write(*,'(2x,a,i0,a)') 'Fragment ',ifrg,': candidate orbitals'
          write(*,'(2x,a)') 'Orb.     Spread  Frag. pop.      FOLI'
          write(*,'(2x,a)') repeat('-',37)
          do ii=1,nel
            xx=frgpop(ifrg,ii)
            if(xx.gt.1.0d-5) write(*,111) ii,frgspr(ifrg,ii),xx,foli(ifrg,ii)
          end do
          write(*,'(2x,a)') repeat('-',37)
        end do
        DEALLOCATE(S0)

!! 2) Lowest FOLI value, and the next one beyond the tolerance !!
!! 1.0d-10: FOLIs equal up to rounding (symmetry-equivalent fragments) !!
!! keep the first one found, so the choice is the same on any machine  !!
        xcutoff=100.0d0
        iifrg=0
        do ifrg=1,icufr
          do ii=1,nel
            if(foli(ifrg,ii).lt.xcutoff-1.0d-10) then
              iifrg=ifrg
              xcutoff=foli(ifrg,ii)
            end if
          end do
        end do
        if(iifrg.eq.0) then
          write(*,'(2x,a)') 'No candidate orbital with a FOLI value '//
     +      'below 100 in this iteration'
          call apost_stop(' OSLO: no candidate orbital')
        end if
        xfront=100.0d0
        jjfrg=0
        do ifrg=1,icufr
          do ii=1,nel
            xx=foli(ifrg,ii)
            if(xcutoff+folitol.lt.xx.and.xx.lt.xfront-1.0d-10) then
              jjfrg=ifrg
              xfront=xx
            end if
          end do
        end do

        write(*,'(a)') ''
        write(line,'(a,i0,a)') 'Lowest FOLI value (fragment ',iifrg,')'
        call print_value(trim(line),xcutoff,5)
        if(jjfrg.gt.0) then
          write(line,'(a,i0,a)') 'Next FOLI value (fragment ',jjfrg,')'
          call print_value(trim(line),xfront,5)
        else
          call print_value('Next FOLI value','none')
        end if

!! Branching: the group at the next FOLI value is selected instead !!
        lbranch=.false.
        do k=1,nbr
          if(ibr(k).eq.iiter) then
            if(jjfrg.gt.0) then
              lbranch=.true.
              ibrused(k)=1
            else
              ibrused(k)=2
            end if
          end if
        end do
        xsel=xcutoff
        if(lbranch) then
          xsel=xfront
          write(*,'(2x,a,f10.5,a,f10.5)') "Branching: FOLI value",
     +      xfront," selected instead of",xcutoff
        else if(nbr.gt.0) then
          do k=1,nbr
            if(ibr(k).eq.iiter) write(*,'(2x,a)') 'Branching asked '//
     +        'at this iteration, but there is no other candidate: '//
     +        'not applied'
          end do
        end if

!! 3) Selection: the lowest group (every candidate within the tolerance !!
!! of the lowest FOLI value), or with branching the next group (beyond   !!
!! the lowest one, up to the next FOLI value plus the tolerance)         !!
        mm=0
        do ifrg=1,icufr
          do iorb=1,nel
            if(ingroup(foli(ifrg,iorb),lbranch)) mm=mm+1
          end do
        end do
        ALLOCATE(infopop(2,mm),sfoli(mm))
        write(*,'(a)') ''
        write(*,'(2x,a)') 'Selected orbitals'
        write(*,'(2x,a)') 'Orb.  Frag.      FOLI'
        write(*,'(2x,a)') repeat('-',21)
        mm=0
        do ifrg=1,icufr
          do iorb=1,nel
            if(ingroup(foli(ifrg,iorb),lbranch)) then
              mm=mm+1
              infopop(1,mm)=ifrg
              infopop(2,mm)=iorb
              sfoli(mm)=foli(ifrg,iorb)
              write(*,112) iorb,ifrg,sfoli(mm)
            end if
          end do
        end do
        write(*,'(2x,a)') repeat('-',21)
        call print_value('OSLOs selected',mm)
        if(jjfrg.gt.0) then
          call print_value('Delta-FOLI',xfront-xcutoff,5)
        else
          call print_value('Delta-FOLI','none (no other candidate)')
        end if

!! 4) Orthogonalizing the selection: Lowdin, from the eigenvalues of its !!
!! overlap above thresh (10^-8, parameter.h). nr of them = the number   !!
!! of orbitals the selection spans                                      !!
        ALLOCATE(ccore(igr,mm),ccoreorth(igr,mm))
        do iorb=1,mm
          do mu=1,igr
            ccore(mu,iorb)=cfrgoslo(infopop(1,iorb),mu,infopop(2,iorb))
          end do
        end do
        ALLOCATE(S0(mm,mm),eigv(mm,mm))
        do ii=1,mm
          do jj=1,mm
            xx=ZERO
            do mu=1,igr
              do nu=1,igr
                xx=xx+ccore(mu,ii)*ccore(nu,jj)*s(mu,nu)
              end do
            end do
            S0(ii,jj)=xx
          end do
        end do
        call diagonalize(mm,mm,S0,eigv,0)
        nr=0
        do kk=1,mm
          if(S0(kk,kk).gt.thresh) nr=nr+1
        end do
        nleft=nel-nnelect
        if(nr.gt.nleft) then
          write(*,'(2x,a,i0,a,i0,a)') 'The selected OSLOs span ',nr,
     +      ' orbitals, but only ',nleft,' are left to assign'
          call apost_stop(' OSLO: more orbitals selected than left')
        end if

        if(nr.eq.mm) then

!! independent selection: each OSLO goes to its own fragment !!
          ALLOCATE(smh(mm,mm))
          do ii=1,mm
            do jj=ii,mm
              smh(jj,ii)=ZERO
              do kk=1,mm
                if(S0(kk,kk).gt.thresh) then
                  xx=eigv(ii,kk)*eigv(jj,kk)
                  ssqrt=dsqrt(S0(kk,kk))
                  smh(jj,ii)=smh(jj,ii)+xx/ssqrt
                end if
              end do
              smh(ii,jj)=smh(jj,ii)
            end do
          end do
          do ii=1,igr
            do jj=1,mm
              xx=ZERO
              do kk=1,mm
                xx=xx+smh(jj,kk)*ccore(ii,kk)
              end do
              ccoreorth(ii,jj)=xx
            end do
          end do
          DEALLOCATE(smh)
          do iorb=1,mm
            ifrg=infopop(1,iorb)
            frgel(ifrg)=frgel(ifrg)+ONE
            iadd=iadd+1
            shrosl(ifrg,iadd)=ONE
            delocoslo(iadd)=sfoli(iorb)
            do mu=1,igr
              coslo(mu,iadd)=ccore(mu,iorb)
              cosloorth(mu,iadd)=ccoreorth(mu,iorb)
            end do
          end do
        else

!! dependent selection (the same orbital found from several fragments): !!
!! it assigns nr orbitals, one per independent direction. Each selected !!
!! OSLO's share is its squared norm after the Lowdin orthogonalization  !!
!! (the shares add up to nr); the orbitals kept are the canonical ones  !!
!! of the selection, which span the same space                          !!
          ALLOCATE(share(mm))
          do iorb=1,mm
            xx=ZERO
            do kk=1,mm
              if(S0(kk,kk).gt.thresh) xx=xx+eigv(iorb,kk)**2
            end do
            share(iorb)=xx
            ifrg=infopop(1,iorb)
            frgel(ifrg)=frgel(ifrg)+xx
          end do
          ir=0
          do kk=1,mm
            if(S0(kk,kk).gt.thresh) then
              ir=ir+1
              ssqrt=dsqrt(S0(kk,kk))
              do mu=1,igr
                xx=ZERO
                do iorb=1,mm
                  xx=xx+ccore(mu,iorb)*eigv(iorb,kk)
                end do
                ccoreorth(mu,ir)=xx/ssqrt
              end do
              iadd=iadd+1
              delocoslo(iadd)=xsel
              do mu=1,igr
                coslo(mu,iadd)=ccoreorth(mu,ir)
                cosloorth(mu,iadd)=ccoreorth(mu,ir)
              end do
            end if
          end do
          write(*,'(a)') ''
          if(nr.eq.1) then
            write(*,'(2x,a,i0,a)') 'Shared orbital: the ',mm,
     +        ' selected OSLOs span 1 orbital'
          else
            write(*,'(2x,a,i0,a,i0,a)') 'Shared orbitals: the ',mm,
     +        ' selected OSLOs span ',nr,' orbitals'
          end if
          write(*,'(2x,a)') 'Frag.     Share'
          write(*,'(2x,a)') repeat('-',15)
          do ifrg=1,icufr
            xx=ZERO
            ii=0
            do iorb=1,mm
              if(infopop(1,iorb).eq.ifrg) then
                xx=xx+share(iorb)
                ii=1
              end if
            end do
            if(ii.eq.1) write(*,'(2x,i5,f10.5)') ifrg,xx
!! each of the nr orbitals of the selection carries 1/nr of the        !!
!! fragment's share (the results file lists them per orbital)          !!
            if(ii.eq.1) then
              do ir=iadd-nr+1,iadd
                shrosl(ifrg,ir)=xx/DBLE(nr)
              end do
            end if
          end do
          write(*,'(2x,a)') repeat('-',15)
          DEALLOCATE(share)
        end if
        nadd=nr
        nnelect=nnelect+nr
        do ir=iadd-nr+1,iadd
          iterosl(ir)=iiter
        end do
        call print_value('Orbitals left to assign',nel-nnelect)

!! 5) Close alternative: a next OSLO with an overlap of nearly 1 with   !!
!! the selected ones (it is nearly one of them, found from another      !!
!! fragment), with a small delta-FOLI. Not evaluated with one orbital   !!
!! left (every candidate is then that orbital)                          !!
        write(*,'(a)') ''
        if(lbranch) then
          write(*,'(2x,a)') 'Close alternatives: not checked in a '//
     +      'branching iteration'
        else if(nleft.eq.1) then
          write(*,'(2x,a)') 'Close alternatives: not checked, one '//
     +      'orbital left'
        else if(jjfrg.eq.0) then
          write(*,'(2x,a)') 'Close alternatives: no other candidate'
        else
          write(*,'(2x,a)') 'Next OSLOs against the selected ones'
          write(*,'(2x,a)') 'Orb.  Frag.      FOLI    Overlap'
          write(*,'(2x,a)') repeat('-',32)
          xmax=ZERO
          do ifrg=1,icufr
            do iorb=1,nel
              if(ingroup(foli(ifrg,iorb),.true.)) then

!! overlap of this next OSLO with the (orthonormal) selected ones: !!
!! sum over them of |<selected|next>|^2, over <next|next>          !!
                do mu=1,igr
                  xx=ZERO
                  do nu=1,igr
                    xx=xx+s(mu,nu)*cfrgoslo(ifrg,nu,iorb)
                  end do
                  st(mu)=xx
                end do
                xnorm=ZERO
                do mu=1,igr
                  xnorm=xnorm+cfrgoslo(ifrg,mu,iorb)*st(mu)
                end do
                xproj=ZERO
                do jj=1,nadd
                  xx=ZERO
                  do mu=1,igr
                    xx=xx+ccoreorth(mu,jj)*st(mu)
                  end do
                  xproj=xproj+xx*xx
                end do
                xov=min(ONE,xproj/xnorm)
                if(xov.gt.xmax) xmax=xov
                write(*,'(2x,i4,i7,f10.5,f11.5)') iorb,ifrg,
     +            foli(ifrg,iorb),xov
              end if
            end do
          end do
          write(*,'(2x,a)') repeat('-',32)
          call print_value('Largest overlap with the selected OSLOs',
     +      xmax,5)
          if(xmax.gt.ONE-xlindep.and.xfront-xcutoff.lt.dfoliwarn) then
            write(*,'(2x,a)') 'WARNING: close alternative selection. '//
     +        'The next OSLO is nearly one of the'
            write(*,'(2x,a)') 'selected ones, and its FOLI value is '//
     +        'within 0.5 of theirs: an alternative'
            write(*,'(2x,a)') 'assignment is possible. To follow it, '//
     +        'add to # OSLO'
            write(*,'(4x,a,1x,i0)') trim(chkey),iiter
            write(*,'(2x,a)') 'and compare the sums of FOLI values '//
     +        'of both runs.'
          end if
        end if

!! 6) Removing the selected orbitals from the remaining density !!
        do ii=1,igr
          do jj=1,igr
            xx=ZERO
            do ij=1,nadd
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

        DEALLOCATE(deloc,foli,infopop,sfoli,ccore,ccoreorth,S0,eigv)

!! In case of all assigned !!
        if(nel-nnelect.eq.0) go to 666

!! End of iterative procedure !!
      end do
666   continue

      xx=ZERO
      do ii=1,nel
        xx=xx+delocoslo(ii)
      end do
      write(*,'(a)') ''
      call print_value('Sum of the FOLI values',xx,5)
      do k=1,nbr
        if(ibrused(k).eq.0) write(*,'(2x,a,i0,a)')
     +    'Branching at iteration ',ibr(k),': not reached, not applied'
      end do

      DEALLOCATE(Sm,Splus,c0,pp0,pcore,st,cmat,cfrgoslo)
      DEALLOCATE(frgpop,frgspr,ibrused)

111   FORMAT(2x,i4,f11.5,f12.5,f10.5)
112   FORMAT(2x,i4,i7,f10.5)

      contains

!! ---- !!
!! whether a FOLI value belongs to the lowest group (next=.false.) or to !!
!! the next one (next=.true.), for the current xcutoff/xfront           !!
      logical function ingroup(f,next)
      real*8,intent(in) :: f
      logical,intent(in) :: next
      if(next) then
        ingroup=f.gt.xcutoff+folitol.and.f.le.xfront+folitol
      else
        ingroup=ABS(xcutoff-f).le.folitol
      end if
      end function ingroup

      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: rwf_iterative_oslo
!! OSLO oxidation states of a restricted (closed-shell) determinant: one
!! channel of doubly occupied orbitals (oslo_channel_iterate), the OSLO
!! .fchk files, the summary of the selected OSLOs, the oxidation states
!! and the oslo block of the results file.
!! arguments:
!!   sat    (in)  (igr,igr,nat) atomic overlap matrices
!!   itotps (in)  number of grid points
!!   wp     (in)  (itotps) integration weights
!!   omp2   (in)  (itotps,nat) atomic weights of each point
!!   chp    (in)  (itotps,igr) basis functions at the points
!!   pcoord (in)  (itotps,3) coordinates of the points
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine rwf_iterative_oslo(sat,itotps,wp,omp2,chp,pcoord)

      use basis_set
      use ao_matrices
      use integration_grid
      use input_options_mod, only: nbranch,ibranchit
      use json_mod
      implicit double precision(a-h,o-z)
      include 'parameter.h'

      integer,intent(in) :: itotps

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      character*20 ctype

      dimension sat(igr,igr,nat)
      dimension chp(itotps,igr),pcoord(itotps,3)
      dimension wp(itotps),omp2(itotps,nat)

      allocatable :: Smat(:,:,:),pnocore(:,:),poslo(:,:)
      allocatable :: coslo(:,:),cosloorth(:,:),delocoslo(:)
      allocatable :: frgel(:),oxi(:),shrosl(:,:)
      integer, allocatable :: iterosl(:)

      ifolitol = iopt(96)
      ifchk    = iopt(98)
      folitol = 10.0d0**(-REAL(ifolitol))

      call oslo_header(.false.,folitol)
      ALLOCATE(Smat(icufr,igr,igr))
      call oslo_build_Smat(itotps,wp,omp2,chp,pcoord,Smat)

      ALLOCATE(pnocore(igr,igr),frgel(icufr),oxi(icufr))
      pnocore=pa
      ALLOCATE(coslo(igr,igr),cosloorth(igr,igr),delocoslo(nocc))
      ALLOCATE(iterosl(nocc),shrosl(icufr,nocc))
      call oslo_channel_iterate(nocc,sat,Smat,pnocore,folitol,
     &  nbranch(1),ibranchit(1,1),'BRANCH_ITERATION',coslo,cosloorth,
     &  delocoslo,frgel,iterosl,shrosl)
      DEALLOCATE(Smat,pnocore)

!! .fchk files of the OSLOs, FOLI values as orbital energies: before   !!
!! orthogonalization on request, the final ones always                 !!
      ALLOCATE(poslo(igr,igr))
      if(ifchk.eq.2) then
        call oslo_density_from_coeffs(igr,nocc,coslo,TWO,poslo)
        ctype="-OSLOs-preortho"
        call oslo_fchk_write(.false.,coslo,coslo,poslo,poslo,delocoslo,
     +    nocc,delocoslo,0,ctype)
      end if
      call oslo_density_from_coeffs(igr,nocc,cosloorth,TWO,poslo)
      ctype="-OSLOs"
      call oslo_fchk_write(.false.,cosloorth,cosloorth,poslo,poslo,
     +  delocoslo,nocc,delocoslo,0,ctype)
      DEALLOCATE(poslo)

      call print_box('SELECTED OSLOs')
      call oslo_summary(sat,nocc,coslo,cosloorth,delocoslo,iterosl)

!! each assigned orbital holds two electrons !!
      call oslo_oxidation_states(frgel,TWO,frgel,ZERO,oxi)

      call json_begin_object('oslo')
      call oslo_json_scheme()
      call json_oslo_channel('total',nocc,delocoslo,iterosl,shrosl)
      call json_begin_object('fragments')
      call json_vector('oxidation_states',oxi,icufr)
      call json_index('oxidation_states','fragment_vector',
     +  'OSLO oxidation states',' ',2,'total',.true.,' ')
      call json_vector('electrons',TWO*frgel,icufr)
      call json_index('electrons','fragment_vector',
     +  'OSLO electrons',' ',2,'total',.false.,' ')
      call json_end()
      call oslo_json_file()
      call json_end()

      DEALLOCATE(coslo,cosloorth,delocoslo,frgel,oxi,iterosl,shrosl)

      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: uwf_iterative_oslo
!! OSLO oxidation states of an unrestricted determinant: the alpha and
!! beta channels separately (oslo_channel_iterate), with the branching
!! iterations of each (BRANCH_ITERATION plus BRANCH_ALPHA/BRANCH_BETA),
!! one OSLO .fchk file with both, the summaries, the oxidation states and
!! the oslo block of the results file.
!! arguments: as rwf_iterative_oslo
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine uwf_iterative_oslo(sat,itotps,wp,omp2,chp,pcoord)

      use basis_set
      use ao_matrices
      use integration_grid
      use input_options_mod, only: maxbranch,nbranch,ibranchit
      use output_mod, only: print_value
      use json_mod
      implicit double precision(a-h,o-z)
      include 'parameter.h'

      integer,intent(in) :: itotps

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      character*20 ctype

      dimension sat(igr,igr,nat)
      dimension chp(itotps,igr),pcoord(itotps,3)
      dimension wp(itotps),omp2(itotps,nat)

      allocatable :: Smat(:,:,:),pnocore(:,:)
      allocatable :: coslo_a(:,:),cosloorth_a(:,:),delocoslo_a(:)
      allocatable :: coslo_b(:,:),cosloorth_b(:,:),delocoslo_b(:)
      allocatable :: frgel_a(:),frgel_b(:),oxi(:)
      allocatable :: shrosl_a(:,:),shrosl_b(:,:)
      integer, allocatable :: iterosl_a(:),iterosl_b(:)
      allocatable :: poslo_a(:,:),poslo_b(:,:)

!! branching iterations of each channel: BRANCH_ITERATION plus its own !!
      dimension ibra(2*maxbranch),ibrb(2*maxbranch)

      ifolitol = iopt(96)
      ifchk    = iopt(98)
      folitol = 10.0d0**(-REAL(ifolitol))

      nbra=nbranch(1)+nbranch(2)
      nbrb=nbranch(1)+nbranch(3)
      do k=1,nbranch(1)
        ibra(k)=ibranchit(k,1)
        ibrb(k)=ibranchit(k,1)
      end do
      do k=1,nbranch(2)
        ibra(nbranch(1)+k)=ibranchit(k,2)
      end do
      do k=1,nbranch(3)
        ibrb(nbranch(1)+k)=ibranchit(k,3)
      end do

      call oslo_header(.true.,folitol)
      ALLOCATE(Smat(icufr,igr,igr))
      call oslo_build_Smat(itotps,wp,omp2,chp,pcoord,Smat)
      ALLOCATE(frgel_a(icufr),frgel_b(icufr),oxi(icufr))

      call print_box('ALPHA ELECTRONS')
      call print_value('Electrons',nalf)
      ALLOCATE(pnocore(igr,igr))
      pnocore=pa
      ALLOCATE(coslo_a(igr,igr),cosloorth_a(igr,igr),delocoslo_a(nalf))
      ALLOCATE(iterosl_a(nalf),shrosl_a(icufr,nalf))
      call oslo_channel_iterate(nalf,sat,Smat,pnocore,folitol,nbra,
     &  ibra,'BRANCH_ALPHA',coslo_a,cosloorth_a,delocoslo_a,frgel_a,
     &  iterosl_a,shrosl_a)

      call print_box('BETA ELECTRONS')
      call print_value('Electrons',nb)
      pnocore=pb
      ALLOCATE(coslo_b(igr,igr),cosloorth_b(igr,igr),delocoslo_b(nb))
      ALLOCATE(iterosl_b(nb),shrosl_b(icufr,nb))
      call oslo_channel_iterate(nb,sat,Smat,pnocore,folitol,nbrb,
     &  ibrb,'BRANCH_BETA',coslo_b,cosloorth_b,delocoslo_b,frgel_b,
     &  iterosl_b,shrosl_b)
      DEALLOCATE(pnocore,Smat)

!! one .fchk file with both channels, FOLI values as orbital energies !!
      ALLOCATE(poslo_a(igr,igr),poslo_b(igr,igr))
      if(ifchk.eq.2) then
        call oslo_density_from_coeffs(igr,nalf,coslo_a,ONE,poslo_a)
        call oslo_density_from_coeffs(igr,nb,coslo_b,ONE,poslo_b)
        ctype="-OSLOs-preortho"
        call oslo_fchk_write(.true.,coslo_a,coslo_b,poslo_a,poslo_b,
     +    delocoslo_a,nalf,delocoslo_b,nb,ctype)
      end if
      call oslo_density_from_coeffs(igr,nalf,cosloorth_a,ONE,poslo_a)
      call oslo_density_from_coeffs(igr,nb,cosloorth_b,ONE,poslo_b)
      ctype="-OSLOs"
      call oslo_fchk_write(.true.,cosloorth_a,cosloorth_b,poslo_a,
     +  poslo_b,delocoslo_a,nalf,delocoslo_b,nb,ctype)
      DEALLOCATE(poslo_a,poslo_b)

      call print_box('SELECTED ALPHA OSLOs')
      call oslo_summary(sat,nalf,coslo_a,cosloorth_a,delocoslo_a,
     +  iterosl_a)
      call print_box('SELECTED BETA OSLOs')
      call oslo_summary(sat,nb,coslo_b,cosloorth_b,delocoslo_b,
     +  iterosl_b)

!! each assigned spin-orbital holds one electron !!
      call oslo_oxidation_states(frgel_a,ONE,frgel_b,ONE,oxi)

      call json_begin_object('oslo')
      call oslo_json_scheme()
      call json_oslo_channel('alpha',nalf,delocoslo_a,iterosl_a,
     +  shrosl_a)
      call json_oslo_channel('beta',nb,delocoslo_b,iterosl_b,shrosl_b)
      call json_begin_object('fragments')
      call json_vector('oxidation_states',oxi,icufr)
      call json_index('oxidation_states','fragment_vector',
     +  'OSLO oxidation states',' ',2,'total',.true.,' ')
      call json_vector('electrons_alpha',frgel_a,icufr)
      call json_index('electrons_alpha','fragment_vector',
     +  'OSLO electrons, alpha',' ',2,'alpha',.false.,' ')
      call json_vector('electrons_beta',frgel_b,icufr)
      call json_index('electrons_beta','fragment_vector',
     +  'OSLO electrons, beta',' ',2,'beta',.false.,' ')
      call json_end()
      call oslo_json_file()
      call json_end()

      DEALLOCATE(coslo_a,cosloorth_a,delocoslo_a,frgel_a,oxi)
      DEALLOCATE(coslo_b,cosloorth_b,delocoslo_b,frgel_b)
      DEALLOCATE(iterosl_a,iterosl_b,shrosl_a,shrosl_b)

      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: oslo_header
!! Opens the OSLO section: the atomic definition of the spreads, the
!! atoms of the fragment populations, the FOLI tolerance and the
!! densities.
!! arguments:
!!   lunres  (in)  logical: unrestricted (alpha and beta separately)
!!   folitol (in)  FOLI tolerance of a selection
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine oslo_header(lunres,folitol)
      use input_options_mod, only: ilow2
      use output_mod, only: print_value
      implicit double precision(a-h,o-z)
      logical lunres
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort

      call atomic_definition_names(ckey,cname,cshort,cspace)
      call print_value('Atomic definition (spreads)',
     +  trim(ckey)//' (real space)')
      if(ilow2.eq.1) then
        call print_value('Fragment populations','MULLIKEN (Hilbert space)')
      else if(ilow2.eq.2) then
        call print_value('Fragment populations','LOWDIN (Hilbert space)')
      else if(ilow2.eq.3) then
        call print_value('Fragment populations',
     +    'LOWDIN-DAVIDSON (Hilbert space)')
      else if(ilow2.eq.6) then
        call print_value('Fragment populations',
     +    'NAO-BASIS (Hilbert space)')
      else
        call print_value('Fragment populations',trim(ckey))
      end if
      call print_value('FOLI tolerance',folitol,5)
      if(lunres) then
        call print_value('Densities','alpha and beta separately')
      else
        call print_value('Densities','total (closed shell)')
      end if
      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: oslo_summary
!! Prints the selected OSLOs of one channel, five per block: before the
!! orthogonalization (iteration, FOLI value and population of each
!! fragment) and after it (the orbitals of the .fchk file, fragment
!! populations).
!! arguments:
!!   sat       (in)  (igr,igr,nat) atomic overlap matrices
!!   nel       (in)  number of OSLOs
!!   coslo     (in)  (igr,igr) OSLOs before orthogonalization
!!   cosloorth (in)  (igr,igr) final OSLOs
!!   foli      (in)  (nel) FOLI value of each
!!   iterosl   (in)  (nel) iteration that selected each
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine oslo_summary(sat,nel,coslo,cosloorth,foli,iterosl)
      implicit double precision(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      dimension sat(igr,igr,nat),coslo(igr,igr),cosloorth(igr,igr)
      dimension foli(*)
      integer iterosl(*)
      allocatable :: orbpopat(:,:),orbpop(:),foslo(:,:)

      if(nel.eq.0) then
        write(*,'(2x,a)') 'No electrons in this channel'
        return
      end if
      ALLOCATE(orbpopat(nel,nat),orbpop(nel),foslo(nel,icufr))

!! the per-atom populations of a set of orbitals are the expensive part, !!
!! computed once per set and then summed per fragment                    !!
      do iset=1,2
        if(iset.eq.1) then
          call rwf_uwf_orbpop_atom(sat,nel,coslo,orbpopat)
          write(*,'(2x,a)') 'Before orthogonalization'
        else
          call rwf_uwf_orbpop_atom(sat,nel,cosloorth,orbpopat)
          write(*,'(a)') ''
          write(*,'(2x,a)') 'After orthogonalization (the orbitals of '//
     +      'the .fchk file)'
        end if
        do jfrg=1,icufr
          orbpop=ZERO
          call rwf_uwf_frg_pop(jfrg,nel,orbpopat,orbpop)
          do ii=1,nel
            foslo(ii,jfrg)=orbpop(ii)
          end do
        end do
        do k0=1,nel,5
          k1=min(nel,k0+4)
          ncol=k1-k0+1
          if(k0.gt.1) write(*,'(a)') ''
          write(*,'(2x,a9,5i10)') 'OSLO     ',(k,k=k0,k1)
          if(iset.eq.1) then
            write(*,'(2x,a9,5i10)') 'Iteration',(iterosl(k),k=k0,k1)
            write(*,'(2x,a9,5f10.5)') 'FOLI     ',(foli(k),k=k0,k1)
          end if
          write(*,'(2x,a)') repeat('-',9+10*ncol)
          do jfrg=1,icufr
            write(*,'(2x,a6,i3,5f10.5)') 'Frag. ',jfrg,
     +        (foslo(k,jfrg),k=k0,k1)
          end do
          write(*,'(2x,a)') repeat('-',9+10*ncol)
        end do
      end do

      DEALLOCATE(orbpopat,orbpop,foslo)
      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: oslo_oxidation_states
!! Oxidation state of each fragment, nuclear charge (reduced by an ECP)
!! minus the electrons of the OSLOs assigned to it, printed with the
!! shared table; says when shared orbitals make some fractional.
!! arguments:
!!   frgel_a, frgel_b (in)  (icufr) orbitals assigned to each fragment,
!!                          per channel (the same array twice, restricted)
!!   occ_a, occ_b     (in)  electrons per orbital of each (restricted:
!!                          2 and 0)
!!   oxi              (out) (icufr) oxidation states
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine oslo_oxidation_states(frgel_a,occ_a,frgel_b,occ_b,oxi)
      implicit double precision(a-h,o-z)
      include 'parameter.h'
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      dimension frgel_a(icufr),frgel_b(icufr),oxi(icufr)

!! the nuclear charge counts the electrons an ECP replaces as absent !!
      ishared=0
      xtot=ZERO
      do ifrg=1,icufr
        izn=0
        do icenter=1,nfrlist(ifrg)
          izn=izn+INT(zn(ifrlist(icenter,ifrg)))
        end do
        xos=DBLE(izn)-(occ_a*frgel_a(ifrg)+occ_b*frgel_b(ifrg))
        if(ABS(xos).lt.5.0d-7) xos=ZERO
        if(ABS(frgel_a(ifrg)-ANINT(frgel_a(ifrg))).gt.1.0d-6.or.
     +    ABS(frgel_b(ifrg)-ANINT(frgel_b(ifrg))).gt.1.0d-6) ishared=1
        oxi(ifrg)=xos
        xtot=xtot+xos
      end do
      call print_oxidation_states(icufr,oxi,xtot,.false.,ZERO)
      if(ishared.eq.1) write(*,'(2x,a)') 'Fractional values: orbitals '//
     +  'shared between fragments (see the iterations)'
      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: oslo_fchk_write
!! Writes <jobname><ctype>.fchk: a copy of the source .fchk with the
!! OSLOs as orbitals, their FOLI values in place of the orbital energies
!! (0 after the last OSLO) and their density. Blocks are replaced
!! wherever they come (Q-Chem writes the energies after the
!! coefficients) and skipped by their own N= size, so any .fchk layout
!! works (G16, G09, pySCF, Q-Chem).
!! arguments:
!!   lunres (in)  logical: unrestricted (Beta blocks and spin density)
!!   cmat_a, cmat_b (in)  (igr,igr) OSLOs of each channel (restricted:
!!                        cmat_a only)
!!   pmat_a, pmat_b (in)  (igr,igr) their densities (restricted: P in
!!                        pmat_a; unrestricted: total Pa+Pb, spin Pa-Pb)
!!   foli_a, foli_b (in)  FOLI values of the OSLOs of each channel
!!   nfa, nfb       (in)  number of OSLOs in each channel
!!   ctype          (in)  file name suffix: "-OSLOs" or "-OSLOs-preortho"
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine oslo_fchk_write(lunres,cmat_a,cmat_b,pmat_a,pmat_b,
     +  foli_a,nfa,foli_b,nfb,ctype)

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /filename/name0

      logical lunres,lfound
      character*80 line
      character*60 name0,name1
      character*20 ctype
      dimension cmat_a(igr,igr),cmat_b(igr,igr)
      dimension pmat_a(igr,igr),pmat_b(igr,igr)
      dimension foli_a(*),foli_b(*)
      allocatable :: ea(:),eb(:)

      indepigr = int_locate(15,"Number of independ",lfound)
      norb     = igr*indepigr
      norbt    = igr*(igr+1)/2
      ALLOCATE(ea(indepigr),eb(indepigr))
      ea=ZERO
      eb=ZERO
      do ii=1,min(nfa,indepigr)
        ea(ii)=foli_a(ii)
      end do
      do ii=1,min(nfb,indepigr)
        eb(ii)=foli_b(ii)
      end do

      name1=trim(name0)//trim(ctype)//".fchk"
      open(unit=69,file=name1)
      rewind(69)
      rewind(15)

      read(15,'(a80)',end=99) line
      call write_fchk_title(69,ctype)
      read(15,'(a80)',end=99) line
      do while(.true.)
        if(index(line,"Alpha Orbital Energies").eq.1) then
          write(69,10) "Alpha Orbital Energies","R","N= ",indepigr
          write(69,13) (ea(ii),ii=1,indepigr)
        else if(lunres.and.index(line,"Beta Orbital Energies").eq.1)
     +    then
          write(69,10) "Beta Orbital Energies ","R","N= ",indepigr
          write(69,13) (eb(ii),ii=1,indepigr)
        else if(index(line,"Alpha MO coefficients").eq.1) then
          write(69,11) "Alpha MO coefficients","R","N= ",norb
          write(69,13) ((cmat_a(ii,jj),ii=1,igr),jj=1,indepigr)
        else if(lunres.and.index(line,"Beta MO coefficients").eq.1)
     +    then
          write(69,11) "Beta MO coefficients ","R","N= ",norb
          write(69,13) ((cmat_b(ii,jj),ii=1,igr),jj=1,indepigr)
        else if(index(line,"Total SCF Density").eq.1) then
          write(69,12) "Total SCF Density","R","N= ",norbt
          if(lunres) then
            write(69,13) ((pmat_a(ii,jj)+pmat_b(ii,jj),jj=1,ii),ii=1,igr)
          else
            write(69,13) ((pmat_a(ii,jj),jj=1,ii),ii=1,igr)
          end if
        else if(lunres.and.index(line,"Spin SCF Density").eq.1) then
          write(69,12) "Spin SCF Density ","R","N= ",norbt
          write(69,13) ((pmat_a(ii,jj)-pmat_b(ii,jj),jj=1,ii),ii=1,igr)
        else
          write(69,'(a80)') line
          read(15,'(a80)',end=99) line
          cycle
        end if

!! a replaced block: skip its data (N= values, 5 per line) !!
        ipos=index(line,"N=")
        read(line(ipos+2:),*) nval
        do ii=1,(nval+4)/5
          read(15,'(a80)',end=99) line
        end do
        read(15,'(a80)',end=99) line
      end do
99    continue
      close(69)
      DEALLOCATE(ea,eb)

10    FORMAT(a22,21x,a1,3x,a3,i11)
11    FORMAT(a21,22x,a1,3x,a3,i11)
12    FORMAT(a17,26x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))

      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: oslo_json_scheme
!! Results file: the atomic definition of the spreads and the atoms of
!! the fragment populations, in the current (oslo) object.
!! arguments: none
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine oslo_json_scheme()
      use input_options_mod, only: ilow2
      use json_mod
      implicit double precision(a-h,o-z)
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort

      call atomic_definition_names(ckey,cname,cshort,cspace)
      call json_value('scheme',trim(ckey))
      if(ilow2.eq.1) then
        call json_value('fragment_populations','MULLIKEN')
      else if(ilow2.eq.2) then
        call json_value('fragment_populations','LOWDIN')
      else if(ilow2.eq.3) then
        call json_value('fragment_populations','LOWDIN-DAVIDSON')
      else if(ilow2.eq.6) then
        call json_value('fragment_populations','NAO-BASIS')
      else
        call json_value('fragment_populations',trim(ckey))
      end if
      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: oslo_json_file
!! Results file: the name of the OSLO .fchk file (no folder).
!! arguments: none
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine oslo_json_file()
      use json_mod
      implicit double precision(a-h,o-z)
      common /filename/name0
      character*60 name0
      character*80 cfile

      cfile=trim(name0)//'-OSLOs.fchk'
      cfile=cfile(index(cfile,'/',back=.true.)+1:)
      call json_value('orbital_file',trim(cfile))
      call json_value('orbital_file_id',json_run_id())
      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: json_oslo_channel
!! Results file: one channel of the oslo block, the sum of its FOLI
!! values and the "oslos" table, one row per OSLO in the order of the
!! channel's orbitals in the .fchk file: iteration, FOLI value, the
!! fragments it was assigned to and their shares, and its column.
!! arguments:
!!   ckey    (in)  'total' (restricted), 'alpha' or 'beta'
!!   nel     (in)  number of OSLOs
!!   foli    (in)  (nel) FOLI values
!!   iterosl (in)  (nel) iterations
!!   shrosl  (in)  (icufr,nel) share of each fragment in each OSLO
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine json_oslo_channel(ckey,nel,foli,iterosl,shrosl)
      use json_mod
      implicit double precision(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      character*(*) ckey
      dimension foli(*),shrosl(icufr,*)
      integer iterosl(*)
      integer, allocatable :: ifr(:),icol(:)
      real*8, allocatable :: xs(:)
      logical, allocatable :: lcol(:)
      logical lfound

      ncolf=min(igr,int_locate(15,"Number of independ",lfound))
      xsum=ZERO
      do k=1,nel
        xsum=xsum+foli(k)
      end do
      call json_begin_object(ckey)
      call json_value('sum_foli',xsum)
      call json_index('sum_foli','scalar','Sum of the FOLI values, '//
     +  ckey,' ',5,ckey,.false.,' ')
      call json_begin_object('oslos')
      call json_value('rows','orbital')
      call json_value('n_rows',nel)
!! the restricted channel is the alpha block of the .fchk file !!
      if(ckey.eq.'beta') then
        call json_value('fchk_orbitals','beta')
      else
        call json_value('fchk_orbitals','alpha')
      end if
      call json_begin_array('columns')
      call json_column('iteration','Iteration',' ',0,.false.)
      call json_column('foli','FOLI',' ',5,.false.)
      call json_column('fragments','Fragments',' ',0,.false.)
      call json_column('shares','Fragment shares',' ',5,.false.)
      call json_column('fchk_column','Column in the .fchk file',' ',0,
     +  .false.)
      call json_end()
      call json_begin_object('data')
      call json_vector('iteration',iterosl(1:nel),nel)
      call json_vector('foli',foli(1:nel),nel)
      ALLOCATE(ifr(icufr),xs(icufr))
      call json_begin_array('fragments')
      do k=1,nel
        n=0
        do ifrg=1,icufr
          if(shrosl(ifrg,k).gt.1.0d-12) then
            n=n+1
            ifr(n)=ifrg
          end if
        end do
        call json_vector('',ifr(1:n),n)
      end do
      call json_end()
      call json_begin_array('shares')
      do k=1,nel
        n=0
        do ifrg=1,icufr
          if(shrosl(ifrg,k).gt.1.0d-12) then
            n=n+1
            xs(n)=shrosl(ifrg,k)
          end if
        end do
        call json_vector('',xs(1:n),n)
      end do
      call json_end()
      ALLOCATE(icol(nel),lcol(nel))
      do k=1,nel
        icol(k)=k
        lcol(k)=k.le.ncolf
      end do
      call json_vector_or_null('fchk_column',icol,lcol,nel)
      call json_end()
      call json_end()
      call json_index('oslos','table','OSLOs, '//ckey,' ',-1,ckey,
     +  .false.,' ')
      call json_end()
      DEALLOCATE(ifr,xs,icol,lcol)
      end

!! ****** !!

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: rwf_uwf_orbpop_atom
!! Population of each of a set of orbitals on each atom, through the
!! atomic overlap matrices (the expensive part of a fragment population,
!! computed once per set of orbitals).
!! arguments:
!!   sat    (in)   (igr,igr,nat) atomic overlap matrices
!!   norb   (in)   number of orbitals
!!   corb   (in)   (igr,igr) their coefficients
!!   orbpop (out)  (norb,nat) population of each orbital on each atom
!! author: MGimf
!! ----------------------------------------------------------------------
      subroutine rwf_uwf_orbpop_atom(sat,norb,corb,orbpop)

      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop

      dimension sat(igr,igr,nat)
      dimension corb(igr,igr),orbpop(norb,nat)

!! corb/sat read-only; each (icenter,iorb) writes only its own slot. !!
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,iorb,jj,kk,xx)
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
!$OMP END PARALLEL DO

      end

!! ****** !!

!! ----------------------------------------------------------------------
!! subroutine: rwf_uwf_frg_pop
!! Population of each orbital on one fragment, the sum of its atoms'
!! populations (rwf_uwf_orbpop_atom).
!! arguments:
!!   ifrg   (in)   fragment
!!   norb   (in)   number of orbitals
!!   orbpop (in)   (norb,nat) populations on each atom
!!   frgpop (out)  (norb) populations on the fragment
!! author: MGimf
!! ----------------------------------------------------------------------
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
