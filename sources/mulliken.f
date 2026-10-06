!! ********************************************************************* !!
!! subroutine: tomull                                                    !!
!! purpose: builds the per-atom Mulliken "overlap" matrix sat(:,:,i) --  !!
!!   a masked copy of S keeping only the rows belonging to atom i, zero  !!
!!   elsewhere. Feeds the same downstream population/bond-order/EFFAO    !!
!!   machinery as every other AIM scheme's sat.                          !!
!! arguments:                                                            !!
!!   sat (out) -- per-atom AO overlap matrix                             !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine tomull(sat)
      use basis_set
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      dimension sat(nbasis,nbasis,natoms)

      do i=1,natoms
        do nu=1,nbasis
          do mu=1,nbasis
            if(ihold(mu).eq.i) then
              sat(nu,mu,i)=s(nu,mu)
            else
              sat(nu,mu,i)=0.0d0
            endif
          enddo
        enddo
      enddo

      end

!! ********************************************************************* !!
!! subroutine: tolow                                                     !!
!! purpose: builds the per-atom Lowdin (imulli=2) or Davidson-Lowdin     !!
!!   (imulli=3) "overlap" matrix sat(:,:,i) -- symmetric S^-1/2          !!
!!   orthogonalization for conventional Lowdin; Davidson-Lowdin instead  !!
!!   orthogonalizes within each atom's own block first, then does the    !!
!!   symmetric step (davidson_lowdin).                                   !!
!! arguments:                                                            !!
!!   sat (out) -- per-atom AO overlap matrix                             !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine tolow(sat)
      use basis_set
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /iops/ iopt(200)
      dimension sat(nbasis,nbasis,natoms)
      allocatable ss(:,:),s12(:,:),x(:,:)

      allocate(ss(nbasis,nbasis),s12(nbasis,nbasis),x(nbasis,nbasis))

      imulli=iopt(5)
      idav=0
      if(imulli.eq.3) idav=1

!! not parallelized: dominated by diagonalize's LAPACK call; the O(nbasis^3) !!
!! matrix builds around it are a plausible future OMP target (same        !!
!! independent-(i,j)-accumulate-over-k pattern already parallelized       !!
!! elsewhere, e.g. ao_to_mo_grid), not done yet in this pass.             !!
      if(idav.eq.1) then
!! Davidson: ss = S_at^1/2 S'^1/2 (see davidson_lowdin) !!
        call davidson_lowdin(nbasis,ss,x)
      else
!! conventional Lowdin: ss = S^1/2 !!
        do i=1,nbasis
          do j=1,nbasis
            ss(i,j)=s(i,j)
          enddo
        enddo
        call diagonalize(nbasis,nbasis,ss,x,0)
        do i=1,nbasis
          do j=1,nbasis
            s12(i,j)=0.0d0
            do k=1,nbasis
              s12(i,j)=s12(i,j)+x(i,k)*dsqrt(ss(k,k))*x(j,k)
            enddo
          enddo
        enddo
        do i=1,nbasis
          do j=1,nbasis
            ss(i,j)=s12(i,j)
          enddo
        enddo
      endif

      do i=1,natoms
        do nu=1,nbasis
          do nu1=1,nbasis
            sat(nu,nu1,i)=0.0d0
            do mu=1,nbasis
              if(ihold(mu).eq.i) sat(nu,nu1,i)=sat(nu,nu1,i)+ss(nu,mu)*ss(nu1,mu)
            enddo
          enddo
        enddo
      enddo
      deallocate(ss,s12,x)

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: davidson_lowdin                                           !!
!! purpose: Davidson-Lowdin orthogonalization T = S_at^-1/2 S'^-1/2:     !!
!!   each atom's basis functions are orthonormalized first (S_at, the    !!
!!   atom-block-diagonal overlap), then a symmetric Lowdin step on        !!
!!   S' = S_at^-1/2 S S_at^-1/2. Shared by tolow (populations) and        !!
!!   ueffaolow_frag (EFOs).                                              !!
!! arguments:                                                            !!
!!   n  (in)  -- number of basis functions                               !!
!!   tp (out) -- S_at^1/2 S'^1/2 = T^-T: tp^T P tp is P in the            !!
!!               orthogonal basis                                        !!
!!   tm (out) -- T = S_at^-1/2 S'^-1/2: orthogonal-basis to AO            !!
!!               coefficients                                            !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine davidson_lowdin(n,tp,tm)
      use basis_set, only: s,ihold,natoms
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      integer, intent(in) :: n
      dimension tp(n,n),tm(n,n)
      allocatable ss(:,:),x(:,:),sm12(:,:),sp12(:,:),s12(:,:),s12m(:,:)

      allocate(ss(n,n),x(n,n),sm12(n,n),sp12(n,n),s12(n,n),s12m(n,n))

!! S_at, and its S^-1/2 (sm12) and S^1/2 (sp12) !!
      do i=1,n
        do j=1,n
          ss(i,j)=0.0d0
        enddo
      enddo
      do i=1,natoms
        do mu=1,n
          if(ihold(mu).eq.i) then
            do nu=1,n
              if(ihold(nu).eq.i) ss(mu,nu)=s(mu,nu)
            enddo
          endif
        enddo
      enddo
      call diagonalize(n,n,ss,x,0)
      do i=1,n
        do j=1,n
          sm12(i,j)=0.0d0
          sp12(i,j)=0.0d0
          do k=1,n
            sm12(i,j)=sm12(i,j)+x(i,k)*x(j,k)/dsqrt(ss(k,k))
            sp12(i,j)=sp12(i,j)+x(i,k)*x(j,k)*dsqrt(ss(k,k))
          enddo
        enddo
      enddo

!! S' = sm12^T S sm12, and its S'^1/2 (s12) and S'^-1/2 (s12m) !!
      do i=1,n
        do j=1,n
          ss(i,j)=0.0d0
          do k=1,n
            ss(i,j)=ss(i,j)+sm12(k,i)*s(k,j)
          enddo
        enddo
      enddo
      do i=1,n
        do j=1,n
          s12(i,j)=0.0d0
          do k=1,n
            s12(i,j)=s12(i,j)+ss(i,k)*sm12(k,j)
          enddo
        enddo
      enddo
      do i=1,n
        do j=1,n
          ss(i,j)=s12(i,j)
        enddo
      enddo
      call diagonalize(n,n,ss,x,0)
      do i=1,n
        do j=1,n
          s12(i,j)=0.0d0
          s12m(i,j)=0.0d0
          do k=1,n
            s12(i,j)=s12(i,j)+x(i,k)*dsqrt(ss(k,k))*x(j,k)
            s12m(i,j)=s12m(i,j)+x(i,k)*x(j,k)/dsqrt(ss(k,k))
          enddo
        enddo
      enddo

!! tp = sp12 S'^1/2, tm = sm12 S'^-1/2 !!
      do i=1,n
        do j=1,n
          tp(i,j)=0.0d0
          tm(i,j)=0.0d0
          do k=1,n
            tp(i,j)=tp(i,j)+sp12(i,k)*s12(k,j)
            tm(i,j)=tm(i,j)+sm12(i,k)*s12m(k,j)
          enddo
        enddo
      enddo

      deallocate(ss,x,sm12,sp12,s12,s12m)

      end


!! ********************************************************************* !!
!! subroutine: mull_opop                                                 !!
!! purpose: Mulliken atom-pair overlap population matrix (op), built     !!
!!   directly from P*S restricted to each atom's own basis-function      !!
!!   range (llim/iulim) -- the Hilbert-space twin of opop's real-space   !!
!!   numerical-integration scheme.                                       !!
!! arguments: none (P/S via ao_matrices, output via COMMON /ovpop/)      !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine mull_opop()
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /ovpop/ op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq

      do iatom=1,natoms
        do jatom=iatom,natoms
          xx=0.0d0
          do i=llim(iatom),iulim(iatom)
            do j=llim(jatom),iulim(jatom)
              xx=xx+p(i,j)*s(j,i)
            enddo
          enddo
          op(iatom,jatom)=xx
          if(iatom.ne.jatom) op(jatom,iatom)=op(iatom,jatom)
        enddo
      enddo
      return
      end
!
!! ********************************************************************* !!
!! subroutine: tonao                                                     !!
!! purpose: builds the per-atom NAO-basis "overlap" matrix sat(:,:,i)    !!
!!   from an externally-generated NAO-to-AO transform, read from         !!
!!   <jobname>.nao (produced by NBO's `$NBO AONAO=W $END` keyword).      !!
!!   Zero active-test coverage (NAO-BASIS keyword untested).             !!
!! arguments:                                                            !!
!!   sat (out) -- per-atom AO overlap matrix                             !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine tonao(sat)
      use basis_set
      use nao_mod, only: unao,ssnao !! replaces common /nao/ -- see modules.f90 !!
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /filename/ name0
      character*60 name0,name3
      dimension sat(nbasis,nbasis,natoms)

      name3=trim(name0)//'.nao'
      open(33,file=name3,status='old')
      write(*,*) 'Reading NAO to AO matrix'
      read(33,*)
      read(33,*)
      read(33,*)
      do j=1,nbasis
        read(33,*) (unao(i,j),i=1,nbasis)
      enddo
      close(33)
      write(*,*) 'NAO to AO matrix read'
      do nu=1,nbasis
        do mu=1,nbasis
          ssnao(nu,mu)=0.0d0
          do k=1,nbasis
            ssnao(nu,mu)=ssnao(nu,mu)+unao(k,nu)*s(k,mu)
          enddo
        enddo
      enddo

      do i=1,natoms
        do nu=1,nbasis
          do mu=1,nbasis
            sat(nu,mu,i)=0.0d0
            do k=1,nbasis
              if(ihold(k).eq.i) sat(nu,mu,i)=sat(nu,mu,i)+ssnao(k,nu)*ssnao(k,mu)
            enddo
          enddo
        enddo
      enddo

      end

!! ********************************************************************* !!
!! subroutine: tolow2                                                    !!
!! purpose: weighted-Lowdin per-atom "overlap" matrix sat(:,:,i) -- each  !!
!!   basis function is scaled by an occupation-dependent weight (mapped   !!
!!   from its total P*S overlap onto [1,Rmax]) before the usual           !!
!!   symmetric S^-1/2 orthogonalization, then unscaled again. Zero        !!
!!   active-test coverage (LOWDIN-W keyword untested); experimental       !!
!!   ("testing weigthed-Lowdin" in the original author's own comment).    !!
!! arguments:                                                            !!
!!   sat (out) -- per-atom AO overlap matrix                             !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine tolow2(sat)
      use basis_set
      use ao_matrices, only: p
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      dimension sat(nbasis,nbasis,natoms)
      allocatable ss(:,:),s12(:,:),x(:,:),wlow(:)

      allocate(ss(nbasis,nbasis),s12(nbasis,nbasis))
      allocate(x(nbasis,nbasis),wlow(nbasis))

!! per-basis-function weight, mapped from its total P*S overlap onto !!
!! [1,Rmax] !!
      do mu=1,nbasis
        xx=0.0d0
        do nu=1,nbasis
          xx=xx+p(mu,nu)*s(nu,mu)
        enddo
        rmax=20.0d0
        aa=2.0d0/rmax-1.0d0
        if(xx.lt.0.0d0) then
          xx2=1.0d0
        else if(xx.gt.2.0d0) then
          xx2=rmax
        else
          xx2=(2.0d0+xx)/(2.0d0+aa*xx)
        endif
        wlow(mu)=xx2
      enddo

      do i=1,nbasis
        do j=1,nbasis
          ss(i,j)=s(i,j)*wlow(i)*wlow(j)
        enddo
      enddo

      call diagonalize(nbasis,nbasis,ss,x,0)
      do i=1,nbasis
        do j=1,nbasis
          s12(i,j)=0.0d0
          do k=1,nbasis
            s12(i,j)=s12(i,j)+x(i,k)*dsqrt(ss(k,k))*x(j,k)
          enddo
        enddo
      enddo

!! unscale by the inverse weights !!
      do i=1,nbasis
        do j=1,nbasis
          ss(i,j)=s12(i,j)/wlow(j)
        enddo
      enddo

      do i=1,natoms
        do nu=1,nbasis
          do mu=1,nbasis
            sat(nu,mu,i)=0.0d0
            do k=1,nbasis
              if(ihold(k).eq.i) sat(nu,mu,i)=sat(nu,mu,i)+ss(k,nu)*ss(k,mu)
            enddo
          enddo
        enddo
      enddo
      deallocate(ss,s12,x,wlow)

      end
