!! *********************************************************************** !!
!! LINEAR ALGEBRA -- natural orbitals, the shared symmetric eigensolver,    !!
!! and the Lowdin-basis transform chain (build_Smp/to_lowdin_basis/         !!
!! to_AO_basis) used throughout OSLO/DFT-DM1-NATORB/gennatural. Only        !!
!! BLAS/LAPACK usage in the codebase (diagonalize, via dsyevd).             !!
!!   gennatural      -- natural orbitals from the density matrix P, run     !!
!!                       once per job when idono=1                         !!
!!   diagonalize     -- symmetric eigensolver (LAPACK dsyevd), used for     !!
!!                       every diagonalization in the program               !!
!!   build_Smp       -- S^-1/2 (and optionally S^1/2) from the overlap      !!
!!                       matrix, via diagonalize                            !!
!!   to_lowdin_basis -- similarity transform into the Lowdin basis          !!
!!   to_AO_basis     -- back-transform from Lowdin basis to AO basis        !!
!!   build_Sinv      -- plain S^-1 (analytic overlap), via diagonalize      !!
!!   reproject_efos_ao -- least-squares projection of EFO*fragment-weight   !!
!!                       onto the AO basis, for the GEOS/EOS EFO .fchk      !!
!!   move_to         -- matrix copy (B=A), used by build_Smp                !!
!!   solvesystem     -- linear system by Gauss-Jordan elimination (POLAR)   !!
!! *********************************************************************** !!

!! ***** !!

!! ********************************************************************* !!
!! subroutine: gennatural                                                !!
!! purpose: builds natural orbitals from the density matrix P (S^1/2     !!
!! transform, diagonalize, back-transform to AO basis); prints their     !!
!! occupation numbers and a sum-of-overlap sanity check. Called from     !!
!! main.f when idono=1 (post-HF or open-shell/kop=1 runs).               !!
!! arguments: none (P/S via ao_matrices, output via occ_no/c_no)         !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine gennatural()
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension dm1(igr,igr)

      allocatable :: s0(:,:),sm(:,:),splus(:,:)

      nbasis=igr

      allocate(s0(igr,igr),sm(igr,igr),splus(igr,igr))

!! make S^1/2 !!
      s0=s
      call build_Smp(igr,s0,sm,splus,0)

!! transform P with S^1/2, diagonalize, back-transform to AO basis !!
      do i=1,igr
        do j=1,igr
          occ_no(i,j)=p(i,j)
        end do
      end do
      call to_lowdin_basis(igr,splus,occ_no)
      call diagonalize(igr,igr,occ_no,c_no,0)
      call to_AO_basis(igr,igr,sm,c_no)

!! i-1 = number of natural orbitals with occupation above 1e-6 (all of !!
!! them when none is below)                                           !!
      i=igr+1
      do ii=1,igr
        if(abs(occ_no(ii,ii)).le.1.0d-6) then
          i=ii
          exit
        end if
      end do

      call print_box('NATURAL ORBITALS')
      write(*,'(2x,a)') 'Occupation numbers (>1.0e-6):'
      write(*,'(2x,8f10.5)') (occ_no(ii,ii),ii=1,i-1)

!! the natural orbitals must be orthonormal: C^T S C = 1 (reusing s0).   !!
!! The largest deviation is round-off (about 1e-14); a warning is        !!
!! printed only above 1e-5. O(igr^3) but run once per job.               !!
      dmax=0.0d0
      do i=1,igr
        do j=1,igr
          xx=0.0d0
          do k=1,igr
            xx=xx+c_no(k,i)*s(k,j)
          end do
          s0(i,j)=xx
        end do
      end do
      do i=1,igr
        do j=1,igr
          xx=0.0d0
          do k=1,igr
            xx=xx+s0(i,k)*c_no(k,j)
          end do
          if(i.eq.j) xx=xx-1.0d0
          dmax=max(dmax,abs(xx))
        end do
      end do
      if(dmax.gt.1.0d-5) write(*,'(2x,a,es9.2)') 'WARNING: natural '//
     +  'orbitals not orthonormal, largest deviation:',dmax

      deallocate(s0,sm,splus)
      return
      end
       
 
! *****
!! ********************************************************************* !!
!! subroutine: build_Smp                                                 !!
!! purpose: builds S^-1/2 (sm) and, if ip.ne.1, S^1/2 (sp) from the       !!
!!   overlap matrix s0, via diagonalization -- the shared Lowdin-basis    !!
!!   transform used throughout OSLO/DFT-DM1-NATORB/gennatural.           !!
!! arguments:                                                            !!
!!   n  (in)  -- matrix order                                            !!
!!   s0 (in)  -- overlap matrix to decompose                             !!
!!   sm (out) -- S^-1/2                                                  !!
!!   sp (out) -- S^1/2 (only if ip.ne.1)                                 !!
!!   ip (in)  -- 1 skips computing sp                                    !!
!! author:                                                                !!
!! ********************************************************************* !!
        Subroutine build_Smp(n,s0,sm,sp,ip)
        implicit double precision(a-h,o-z)
        include 'parameter.h'
        integer, intent(in) :: n
        dimension s0(n,n),sm(n,n),sp(n,n)
        allocatable :: x(:,:),s(:,:)

        allocate(x(n,n),s(n,n))

!! diagonalize S: x holds eigenvectors, s's diagonal holds eigenvalues   !!
        call move_to(N,s0,s)
        call diagonalize(n,n,s,x,0)

!! parallelization: not done. This O(n^3) triple loop, plus the calling  !!
!! chain's own to_lowdin_basis/to_AO_basis, is called once per fragment  !!
!! per iteration from OSLO/DFT-DM1-NATORB -- fixing it needs a threaded/ !!
!! BLAS eigensolver or fewer diagonalizations, a bigger, separately-     !!
!! scoped task, not a local fix here.                                    !!
        do i=1,n
         do j=i,n
           Sm(j,i)=0.0d0
           if (ip.ne.1) Sp(j,i)=0.0d0
           do k=1,n
            if(s(k,k).gt.1.0d-8) then
             xx=x(i,k)*x(j,k)
             xs=sqrt(s(k,k))
             Sm(j,i)=Sm(j,i)+xx/xs
             if(ip.ne.1) Sp(j,i)=Sp(j,i)+xx*xs
            end if
           end do
           Sm(i,j)=Sm(j,i)
           if (ip.ne.1) Sp(i,j)=Sp(j,i)
         end do
        end do
        deallocate(x,s)
        return
        end

! *****
!! ********************************************************************* !!
!! subroutine: to_lowdin_basis                                           !!
!! purpose: similarity transform f -> x^T f x in place, e.g. projecting  !!
!!   a Fock/density matrix (f) into the Lowdin basis defined by x         !!
!!   (typically S^-1/2 from build_Smp).                                  !!
!! arguments:                                                            !!
!!   n (in)    -- matrix order                                           !!
!!   x (in)    -- basis-change matrix                                    !!
!!   f (inout) -- matrix to transform                                    !!
!! author:                                                                !!
!! ********************************************************************* !!
        Subroutine to_lowdin_basis(N,X,F)
        implicit double precision(a-h,o-z)
        include 'parameter.h'
        integer,  intent(in) :: n
        dimension f(n,n), X(n,n)
        allocatable :: f0(:,:)

        allocate (f0(n,n))

!! parallelization: not done -- same O(n^3)/per-call-chain reasoning as  !!
!! build_Smp above (part of the same bigger, separately-scoped task).    !!
        do j=1,n
         do i=1,n
          xx=0.0d0
          do k=1,n
           xx=xx+f(i,k)*x(k,j)
          end do
           f0(i,j)=xx
         end do
        end do
        do j=1,n
          do i=1,n
           f(i,j)=f0(i,j)
          end do
        end do
        do j=1,n
         do i=1,n
          xx=0.0d0
          do k=1,n
           xx=xx+x(k,i)*f(k,j)
          end do
           f0(i,j)=xx
         end do
        end do
        do j=1,n
          do i=1,n
           f(i,j)=f0(i,j)
          end do
        end do
        deallocate (f0)
        return
        end

! *****

!! ********************************************************************* !!
!! subroutine: build_Sinv                                                !!
!! purpose: inverts the analytic AO overlap (basis_set's s) via          !!
!!   eigendecomposition (S=U*L*U^T, S^-1=U*L^-1*U^T) -- same             !!
!!   diagonalize-based pattern as build_Smp's S^-1/2/S^1/2 above, just    !!
!!   the plain inverse. Feeds reproject_efos_ao's least-squares solve.    !!
!! arguments:                                                             !!
!!   igr  (in)  -- basis size                                            !!
!!   Sinv (out) -- S^-1                                                  !!
!! author: MGimf                                                          !!
!! ********************************************************************* !!
        subroutine build_Sinv(igr,Sinv)
        use basis_set, only: s
        implicit real*8(a-h,o-z)
        include 'parameter.h'
        integer, intent(in) :: igr
        dimension Sinv(igr,igr)
        allocatable :: Seig(:,:),Uvec(:,:)

        ALLOCATE(Seig(igr,igr),Uvec(igr,igr))
        Seig=s
        call diagonalize(igr,igr,Seig,Uvec,0)

        do mu=1,igr
          do nu=1,mu
            xx=ZERO
            do kk=1,igr
              if(Seig(kk,kk).gt.1.0d-8)
     +          xx=xx+Uvec(mu,kk)*Uvec(nu,kk)/Seig(kk,kk)
            end do
            Sinv(mu,nu)=xx
            Sinv(nu,mu)=xx
          end do
        end do

        DEALLOCATE(Seig,Uvec)
        end

! *****

!! ********************************************************************* !!
!! subroutine: reproject_efos_ao                                         !!
!! purpose: least-squares projection of EFO(r)*w_F(r) (EFO times its     !!
!!   fragment's real-space weight) onto the AO basis, S*c'=b, for EFO    !!
!!   columns ilo:ihi of c0 -- the coefficients exported to the GEOS/EOS  !!
!!   EFO .fchk, so a viewer approximates the fragment-weighted cube.     !!
!! arguments:                                                            !!
!!   igr    (in)  -- basis size                                          !!
!!   itotps (in)  -- total grid points                                   !!
!!   chp    (in)  -- basis-function values at each grid point            !!
!!   wp     (in)  -- integration weight of each grid point               !!
!!   omp    (in)  -- own-atom weight of each grid point (multi-center    !!
!!                   quadrature factor)                                  !!
!!   scr    (in)  -- fragment weight w_F at each grid point              !!
!!   Sinv   (in)  -- inverse analytic AO overlap (build_Sinv)            !!
!!   c0     (in)  -- (igr,igr) raw EFO coefficients                      !!
!!   ilo,ihi(in)  -- column range to project                             !!
!!   cproj  (out) -- (igr,igr), columns ilo:ihi get the projected        !!
!!                   coefficients; other columns untouched               !!
!!   xfit   (out) -- ||residual||/||target|| per EFO (0 = exact), xfit(1)!!
!!                   for column ilo ... xfit(ihi-ilo+1) for column ihi   !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
        subroutine reproject_efos_ao(igr,itotps,chp,wp,omp,scr,Sinv,
     +    c0,ilo,ihi,cproj,xfit)
        implicit real*8(a-h,o-z)
        include 'parameter.h'
        integer, intent(in) :: igr,itotps,ilo,ihi
        dimension chp(itotps,igr),wp(itotps),omp(itotps),scr(itotps)
        dimension Sinv(igr,igr),c0(igr,igr),cproj(igr,igr),xfit(igr)
        allocatable :: efoval(:,:),bvec(:,:)

        nn=ihi-ilo+1
        if(nn.le.0) return

        ALLOCATE(efoval(itotps,nn),bvec(igr,nn))

!! EFO(r) at every grid point, for every target column in the block.     !!
!! parallel over ifut: each point independent, writes only its own row;  !!
!! c0/chp shared read-only.                                               !!
!$OMP PARALLEL DO PRIVATE(ifut,jj,iorb,mu,xx)
        do ifut=1,itotps
          do jj=1,nn
            iorb=ilo+jj-1
            xx=ZERO
            do mu=1,igr
              xx=xx+c0(mu,iorb)*chp(ifut,mu)
            end do
            efoval(ifut,jj)=xx
          end do
        end do
!$OMP END PARALLEL DO

!! b_mu = integral of chi_mu*EFO*w_F, same multi-center quadrature       !!
!! (wp*omp) as the caller's S0 -- without omp every atom's grid would    !!
!! count all space once, overcounting ~nat-fold.                         !!
!! parallel over mu: each mu writes only its own row of bvec.             !!
!$OMP PARALLEL DO PRIVATE(mu,jj,ifut,xx)
        do mu=1,igr
          do jj=1,nn
            xx=ZERO
            do ifut=1,itotps
              xx=xx+wp(ifut)*omp(ifut)*chp(ifut,mu)*efoval(ifut,jj)*scr(ifut)
            end do
            bvec(mu,jj)=xx
          end do
        end do
!$OMP END PARALLEL DO

        do jj=1,nn
          iorb=ilo+jj-1
          do mu=1,igr
            xx=ZERO
            do nu=1,igr
              xx=xx+Sinv(mu,nu)*bvec(nu,jj)
            end do
            cproj(mu,iorb)=xx
          end do
        end do

!! fit quality per EFO, ||residual||/||target|| on the same quadrature.  !!
!! Threaded over grid points, not EFOs (often only 1-2 per call):        !!
!! xres/xtarget are genuine reductions; cproj/chp/efoval/scr read-only.  !!
        do jj=1,nn
          iorb=ilo+jj-1
          xres=ZERO
          xtarget=ZERO
!$OMP PARALLEL DO PRIVATE(ifut,mu,xproj,xtv) REDUCTION(+:xres,xtarget)
          do ifut=1,itotps
            xproj=ZERO
            do mu=1,igr
              xproj=xproj+cproj(mu,iorb)*chp(ifut,mu)
            end do
            xtv=efoval(ifut,jj)*scr(ifut)
            xres=xres+wp(ifut)*omp(ifut)*(xproj-xtv)**2
            xtarget=xtarget+wp(ifut)*omp(ifut)*xtv**2
          end do
!$OMP END PARALLEL DO
          xfit(jj)=ZERO
          if(xtarget.gt.1.0d-8) xfit(jj)=dsqrt(xres/xtarget)
        end do

        DEALLOCATE(efoval,bvec)
        end


!! ********************************************************************* !!
!! subroutine: diagonalize                                               !!
!! purpose: symmetric eigensolver used everywhere in the program, with  !!
!!   LAPACK's dsyevd; eigenvalues sorted from largest to smallest.       !!
!! arguments:                                                            !!
!!   M     (in)    -- leading declared dimension of a0/x (LDA)          !!
!!   N     (in)    -- actual matrix order to diagonalize (N<=M)          !!
!!   A0    (inout) -- in: matrix to diagonalize. out: zeroed except the  !!
!!                    diagonal, which holds the N eigenvalues           !!
!!   X     (out)   -- eigenvectors as columns, same order as A0's diag   !!
!!   ival  (in)    -- not used                                          !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      Subroutine diagonalize(M,N,A0,X,ival)
      implicit double precision(a-h,o-z)
      include 'parameter.h'
      integer, intent(in) :: M,N,ival
      dimension x(M,M),a0(M,M)
      double precision, allocatable :: w(:),work(:),acopy(:,:)
      integer :: lwork,liwork,info
      integer, allocatable :: iwork(:)

      allocate(acopy(N,N),w(N))
      acopy(1:N,1:N)=a0(1:N,1:N)

!! workspace query, then the real call !!
      allocate(work(1),iwork(1))
      call dsyevd('V','U',N,acopy,N,w,work,-1,iwork,-1,info)
      lwork=int(work(1))
      liwork=iwork(1)
      deallocate(work,iwork)
      allocate(work(max(1,lwork)),iwork(max(1,liwork)))
      call dsyevd('V','U',N,acopy,N,w,work,lwork,iwork,liwork,info)
      if(info.ne.0) call apost_stop('dsyevd failed to converge in diagonalize')

!! dsyevd returns ascending order; every caller expects descending:    !!
!! reverse both eigenvalues and the matching eigenvector columns while  !!
!! copying out.                                                         !!
      a0(1:N,1:N)=0.0d0
      do i=1,N
        a0(i,i)=w(N-i+1)
        x(1:N,i)=acopy(1:N,N-i+1)
      end do

      deallocate(acopy,w,work,iwork)
      return
      end
! *****
!! ********************************************************************* !!
!! subroutine: move_to                                                   !!
!! purpose: copies matrix A into B (B=A), n x n.                         !!
!! arguments: n (in), A (in), B (out)                                    !!
!! author:                                                                !!
!! ********************************************************************* !!
        Subroutine move_to(n,A,B)
        implicit double precision(a-h,o-z)
        include 'parameter.h'
        integer, intent(in) :: n
        dimension A(n,n),B(n,n)

!! parallelization: not worth it -- O(n^2) memory-bound copy, thread     !!
!! overhead would exceed the benefit at any realistic n.                 !!
        do j=1,N
         do i=1,N
          b(i,j)=a(i,j)
         end do
        end do
        return
        end
! *****
!! ********************************************************************* !!
!! subroutine: to_AO_basis                                               !!
!! purpose: back-transforms the first na columns of C from the Lowdin    !!
!!   basis to the AO basis, C(:,1:na) = X . C(:,1:na) (typically X is    !!
!!   S^-1/2 from build_Smp).                                             !!
!! arguments:                                                            !!
!!   n  (in)    -- matrix order                                          !!
!!   na (in)    -- number of columns (orbitals) to transform             !!
!!   x  (in)    -- basis-change matrix                                   !!
!!   c  (inout) -- orbital-coefficient matrix                            !!
!! author:                                                                !!
!! ********************************************************************* !!
        Subroutine to_AO_basis(n,na,X,C)
        implicit double precision(a-h,o-z)
        include 'parameter.h'
        integer, intent(in) :: n
        dimension x(n,n),c(n,n)
        allocatable :: cx(:,:)

        allocate (cx(n,n))

!! parallelization: not done -- same reasoning as build_Smp/            !!
!! to_lowdin_basis above.                                                !!
        do j=1,na
         do i=1,n
           cx(i,j)=0.0d0
            do k=1,n
             cx(i,j)=cx(i,j)+x(i,k)*c(k,j)
            end do
         end do
        end do
 
        do j=1,na
         do i=1,n
          c(i,j)=cx(i,j)
         end do
        end do
 
        deallocate (cx)
        return
       end

!! ********************************************************************* !!
!! subroutine: solvesystem                                               !!
!! purpose: solves the n linear equations stored as the augmented matrix !!
!!   A (n x n+1) by Gauss-Jordan elimination with row pivoting; used by  !!
!!   POLAR. ilog is false (and the system left unsolved) when a pivot is !!
!!   zero.                                                               !!
!! arguments:                                                            !!
!!   n     (in)    -- number of equations                                !!
!!   ndim  (in)    -- leading dimension of A and C                       !!
!!   A     (inout) -- (ndim,ndim+1) augmented matrix, destroyed          !!
!!   C     (out)   -- (ndim) solution                                    !!
!!   ilog  (out)   -- .true. when the system was solved                  !!
!! author:                                                               !!
!! ********************************************************************* !!
       subroutine solvesystem(n,ndim,A,C,ilog)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       integer, intent(in) :: n,ndim
       DIMENSION A(ndim,ndim+1),C(ndim)
       logical ilog
 
       ilog=.true.

      do i=1,N
c pivoting
       do l=i,N
        do m=l,N
         if (abs(A(m,i)).gt.abs(A(l,i))) then
          do kk=1,N+1
           f=A(m,kk)
           A(m,kk)=A(l,kk)
           A(l,kk)=f
          end do
         end if
        end do
       end do
c  final pivotin
c compatibilidad del sistema
       if (Abs(A(i,i)).lt.1e-15)then
        write(*,*) 'sistema no comp. determinat'
        ilog=.false.
        return
       end if
       
       do k=1,N
        if (k.ne.i)then
         alfa=A(k,i)/A(i,i)
         do j=1,N+1
          A(k,j)=A(k,j)-alfa*A(i,j)
         end do
        end if
       end do
      end do
   
      do i=1,N
       P=A(i,i)
       do j=1,N+1
        A(i,j)=A(i,j)/P
       end do
      end do

      do i=1,n   
       c(i)= A(i,n+1)
      end do
      end
