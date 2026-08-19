!! *********************************************************************** !!
!! NUMERICAL INTEGRATION CORE -- grid construction, AO/density evaluation, !!
!! atomic-orbital overlap integration                                      !!
!! Per-iteration grid/density driver:                                      !!
!!   prenumint       -- builds the grid (pcoord), AO values (chp, via      !!
!!                       fpoints/rpoints), electron density (rho), and     !!
!!                       becke/tfvc/hirshfeld atomic weights (omp/omp2)    !!
!! Grid-point evaluators called from prenumint (and, for dpoints, from     !!
!! enpart.f/enpart_phf.f directly):                                        !!
!!   fpoints         -- basis-function values at every grid point          !!
!!   rpoints         -- integration weight of every grid point             !!
!!   dpoints         -- Laplacian of every basis function at every grid    !!
!!                       point                                             !!
!! AO->MO transform helpers (shared by enpart.f/enpart_dft.f, replacing    !!
!! a hardcoded loop once duplicated across both):                          !!
!!   ao_to_mo_grid   -- transforms AO grid values to MO grid values        !!
!!   ao_to_mo_grid_t -- same, plus returns the MO-major transpose          !!
!! Atomic-orbital overlap integration (runs by default for every           !!
!! calculation that isn't Hilbert-space Mulliken/Lowdin/NAO):              !!
!!   numint_sat      -- per-atom AO overlap matrix (sat), three mutually   !!
!!                       exclusive integration schemes (iallpo/iqtaim)     !!
!! Hirshfeld/Hirshfeld-Iterative promolecular density support:             !!
!!   spline          -- cubic-spline second-derivative setup for a         !!
!!                       tabulated free-atom radial density profile,       !!
!!                       consumed by splint() (wat.f)                      !!
!!   atdens_int      -- integrates the density (and its anisotropy) over   !!
!!                       a sphere around one atom, on its own small grid   !!
!! Note: spline's zero-density branch never sets y2(iat,ich,n), unlike     !!
!! the normal branch, which splint() (wat.f) can read whenever its         !!
!! binary search lands khi=n. For the moment... left as-is.                !!
!! *********************************************************************** !!

!! ***** !!

!! ********************************************************************* !!
!! subroutine: prenumint                                                 !!
!! purpose: builds the atom-centered numerical integration grid, the     !!
!!   value of every basis function at every grid point, the electron     !!
!!   density at every grid point, and the becke/tfvc (or hirshfeld)      !!
!!   atomic weight of every grid point for every atom. called once per   !!
!!   iteration, right after build_integration_grid(); feeds numint_sat,  !!
!!   numint_one and numint_two.                                          !!
!! arguments:                                                            !!
!!   ndim      (in)  -- number of basis functions (leading dim of chp)   !!
!!   itotps    (in)  -- total number of grid points (nat*iatps)          !!
!!   nat0      (in)  -- number of atoms (leading dim of omp2)            !!
!!   wp        (out) -- integration weight of each grid point            !!
!!   omp       (out) -- becke/tfvc weight of each point for its own atom !!
!!   omp2      (out) -- becke/tfvc (or hirshfeld) weight of each point   !!
!!                       for every atom                                  !!
!!   chp       (out) -- value of each basis function at each grid point  !!
!!   rho       (out) -- electron density at each grid point              !!
!!   pcoord    (out) -- xyz coordinates of each grid point               !!
!!   ibaspoint (--)   -- qtaim basin assignment; unused here, the qtaim  !!
!!                       grid-building block below is disabled           !!
!!   iiter     (in)  -- scf iteration counter, only read by the          !!
!!                       hirshfeld-iterative (ihirsh=2) branch           !!
!! author: PSalse, ERaco, MGimf                                          !!
!! ********************************************************************* !!
      subroutine prenumint(ndim,itotps,nat0,wp,omp,omp2,chp,rho,pcoord,ibaspoint,iiter)
      use basis_set, only: coord
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim,itotps,nat0
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      dimension wp(itotps),chp(itotps,ndim),omp(itotps),rho(itotps)
      dimension pcoord(itotps,3),ibaspoint(itotps),omp2(itotps,nat0)

!! iopt flags used by this routine !!
      imulli = Iopt(5)
      ihirsh = Iopt(6)
      iallpo = Iopt(7)
      ieffao = Iopt(12)
      icube  = Iopt(13)
      ibcp   = Iopt(14)
      iqtaim = Iopt(16)

      if(iqtaim.eq.1) iallpo=1
      iatps=Nang*NRad

!! building pcoords !!
      ifut=1
      do icenter=1,nat
        do k=1,nrad 
          do i=1,nang 
            thx=th(i)
            fix=ph(i)
            rr=xr(k)
            xabs0=rr*dsin(thx)*dcos(fix)
            yabs0=rr*dsin(thx)*dsin(fix)
            zabs0=rr*dcos(thx)
            xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
            yabs=-dcos(pha)*dsin(phb)*xabs0+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
            zabs=dsin(pha)*dsin(phb)*xabs0-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0
            xx0=xabs+coord(1,icenter)
            yy0=yabs+coord(2,icenter)
            zz0=zabs+coord(3,icenter)
            pcoord(ifut,1)=xx0
            pcoord(ifut,2)=yy0
            pcoord(ifut,3)=zz0
            ifut=ifut+1
          enddo
        enddo
      enddo
      call rpoints(wp)
      call fpoints(chp,pcoord)

!! building rho: parallel over grid points: each ifut only reads the    !!
!! shared, read-only density matrix p and its own row chp(ifut,:), and  !!
!! writes only its own rho(ifut)                                        !!
!$OMP PARALLEL DO PRIVATE(ifut,mu,nu,x)
      do ifut=1,itotps
        x=0.0d0
        do mu=1,igr
          do nu=mu+1,igr
            x=x+p(mu,nu)*chp(ifut,mu)*chp(ifut,nu)*2.0d0
          end do
          x=x+p(mu,mu)*chp(ifut,mu)*chp(ifut,mu)
        end do
        rho(ifut)=x
      end do
!$OMP END PARALLEL DO

!! building aim weights for all gridpoints !!
!! parallel over (icenter,k) grid-point pairs. ifut is now computed            !!
!! directly from icenter/k instead of carried as a serially-incremented        !!
!! counter, since a plain "ifut=ifut+1" is not safe once the loop is           !!
!! split across threads -- the computed form gives the exact same values       !!
!! (icenter=1,k=1..iatps -> ifut=1..iatps, icenter=2 -> ifut=iatps+1..2*iatps, !!
!! etc.) as the original serial counter, so results are unchanged. each        !!
!! iteration writes only its own omp(ifut)/omp2(ifut,:) and reads shared,      !!
!! read-only data (pcoord); wat()/wathirsh() are themselves thread-safe        !!
!! (no shared mutable state -- see wat.f, chi is now a passed argument).       !!
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,k,ifut,xx0,yy0,zz0,jcenter)
      do icenter=1,nat
        do k=1,iatps
          ifut=(icenter-1)*iatps+k
          xx0=pcoord(ifut,1)
          yy0=pcoord(ifut,2)
          zz0=pcoord(ifut,3)
          do jcenter=1,nat
            if(ihirsh.ne.0) then
              omp2(ifut,jcenter)=wathirsh(jcenter,xx0,yy0,zz0)
            else if (iqtaim.ne.1)then
              omp2(ifut,jcenter)=wat(jcenter,xx0,yy0,zz0)
            end if
          end do
          omp(ifut)=wat(icenter,xx0,yy0,zz0)
        enddo
      enddo
!$OMP END PARALLEL DO
      if(ihirsh.eq.2) then
        call wathirshit3(rho,iatps,wp,omp2,nat0,iiter)
      end if

!! QTAIM basin assignment/Laplacian integration is not available in this !!
!! version (QTAIM is hard-disabled, main.f:218) -- the block that used   !!
!! to build ibaspoint and integrate the Laplacian here was removed       !!
!! To implement it properly, we can get from older versions... a TODO.   !!

      return
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: ao_to_mo_grid                                             !!
!! purpose: transforms basis-function (AO) values on the numerical grid  !!
!!   into molecular-orbital (MO) values, chpmo(k,j)=sum_i coef(i,j)*     !!
!!   chpao(k,i). Replaces the hardcoded "TO MOs" loop duplicated across  !!
!!   enpart.f/enpart_dft.f (one call per spin: pass c/nocc for RHF,      !!
!!   c/nalf then cb/nb for UHF).                                         !!
!! arguments:                                                            !!
!!   itotps (in)  -- number of grid points                               !!
!!   igr    (in)  -- number of basis functions (AOs)                     !!
!!   nmo    (in)  -- number of MOs to transform (nocc, nalf or nb)       !!
!!   coef   (in)  -- AO coefficient matrix (igr,igr) -- c or cb          !!
!!   chpao  (in)  -- AO values at each grid point (itotps,igr)           !!
!!   chpmo  (out) -- MO values at each grid point (itotps,nmo)           !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine ao_to_mo_grid(itotps,igr,nmo,coef,chpao,chpmo)

      implicit real*8(a-h,o-z)

      integer, intent(in) :: itotps,igr,nmo
      dimension coef(igr,igr),chpao(itotps,igr),chpmo(itotps,nmo)

!! parallel over grid points: each k only reads its own chpao(k,:) and   !!
!! the shared, read-only coef, and writes only its own chpmo(k,:).       !!
!$OMP PARALLEL DO PRIVATE(k,j,i,xx)
      do k=1,itotps
        do j=1,nmo
          xx=ZERO
          do i=1,igr
            xx=xx+coef(i,j)*chpao(k,i)
          end do
          chpmo(k,j)=xx
        end do
      end do
!$OMP END PARALLEL DO

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: ao_to_mo_grid_t                                           !!
!! purpose: same as ao_to_mo_grid, but also returns the transpose        !!
!!   (chpmot(j,k)=chpmo(k,j)) in one pass -- some downstream loops       !!
!!   (numint_two's exchange kernels) want the MO-major layout for        !!
!!   cache-friendly access.                                              !!
!! arguments:                                                            !!
!!   itotps (in)  -- number of grid points                               !!
!!   igr    (in)  -- number of basis functions (AOs)                     !!
!!   nmo    (in)  -- number of MOs to transform                          !!
!!   coef   (in)  -- AO coefficient matrix (igr,igr) -- c or cb          !!
!!   chpao  (in)  -- AO values at each grid point (itotps,igr)           !!
!!   chpmo  (out) -- MO values at each grid point (itotps,nmo)           !!
!!   chpmot (out) -- transpose of chpmo (nmo,itotps)                     !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine ao_to_mo_grid_t(itotps,igr,nmo,coef,chpao,chpmo,chpmot)

      implicit real*8(a-h,o-z)

      integer, intent(in) :: itotps,igr,nmo
      dimension coef(igr,igr),chpao(itotps,igr),chpmo(itotps,nmo),chpmot(nmo,itotps)

!! parallel over grid points: each k only reads its own chpao(k,:) and   !!
!! the shared, read-only coef, and writes only its own chpmo(k,:)/       !!
!! chpmot(:,k).                                                          !!
!$OMP PARALLEL DO PRIVATE(k,j,i,xx)
      do k=1,itotps
        do j=1,nmo
          xx=ZERO
          do i=1,igr
            xx=xx+coef(i,j)*chpao(k,i)
          end do
          chpmo(k,j)=xx
          chpmot(j,k)=xx
        end do
      end do
!$OMP END PARALLEL DO

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: fpoints                                                   !!
!! purpose: evaluates every basis function at every numerical-grid       !!
!!   point (chp(irun,iact)), contracting each primitive Gaussian over    !!
!!   its shell. called once per iteration from prenumint, right after    !!
!!   the grid coordinates (pcoord) are built.                            !!
!! arguments:                                                            !!
!!   chp    (out) -- value of each basis function at each grid point     !!
!!   pcoord (in)  -- xyz coordinates of each grid point                  !!
!! author: PSalse                                                        !!
!! ********************************************************************* !!
      subroutine fpoints(chp,pcoord)
      use basis_set
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      dimension chp(nrad*nang*natoms,nbasis),pcoord(nrad*nang*natoms,3)

      iatps=nrad*nang
      itotps=iatps*natoms

!! parallel over grid points: each irun only reads shared, read-only     !!
!! basis-set data and its own pcoord(irun,:), and writes only its own    !!
!! chp(irun,:) -- no dependency between iterations.                      !!
!$OMP PARALLEL DO PRIVATE(irun,iact,iactat,x,y,z,rr,f,k,ipr,nn,ll,mm)
      do irun=1,itotps
         do iact=1,nbasis
          iactat=ihold(iact)
          x=pcoord(irun,1)-coord(1,iactat)
          y=pcoord(irun,2)-coord(2,iactat)
          z=pcoord(irun,3)-coord(3,iactat)
          rr=dsqrt(x**2.0d0+y**2.0d0+z**2.0d0)
          f=0.d0
          k=1
          do while(nprimbas(k,iact).ne.0)
           ipr=nprimbas(k,iact)
           nn=nlm(ipr,1)
           ll=nlm(ipr,2)
           mm=nlm(ipr,3)
           f=f+(x**nn)*(y**ll)*(z**mm)*dexp(-expp(ipr)*(rr**2))*coefpb(ipr,iact)
           k=k+1
          enddo
          chp(irun,iact)=f
         enddo
      end do
!$OMP END PARALLEL DO
      return
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: rpoints                                                   !!
!! purpose: builds the integration weight of every numerical-grid point  !!
!!   (radial quadrature weight times angular weight times 4*pi*r^2,      !!
!!   folded into wr). called once per iteration from prenumint.          !!
!! arguments:                                                            !!
!!   wp (out) -- integration weight of each grid point                   !!
!! author: PSalse                                                        !!
!! ********************************************************************* !!
      subroutine rpoints(wp)
      use basis_set
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      dimension wp(nang*nrad*natoms)

      iatps=nrad*nang

      irun=1
      do icenter=1,natoms
        do kk=1,Nrad
          xxr=wr(kk)*xr(kk)*xr(kk)
          do i=1,Nang
            wp(irun)=w(i)*xxr*4.d0*Pi
            irun=irun+1
          end do
        end do
      end do

      return
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: spline                                                    !!
!! purpose: cubic-spline second-derivative setup (Numerical Recipes      !!
!!   spline) for a free-atom radial density profile y(iat,ich,:) tabu-   !!
!!   lated at x(:) -- consumed later by splint() (wat.f) for Hirshfeld/  !!
!!   Hirshfeld-iterative promolecular density lookups. Skips the fit     !!
!!   entirely (zeroes y2) for the H+ special case, where the tabulated   !!
!!   density is identically zero.                                        !!
!! arguments:                                                            !!
!!   iat,ich (in) -- atom-type/charge-state index into y2/y/x            !!
!!   n       (in) -- number of tabulated radial points                   !!
!!   yp1,ypn (in) -- first-derivative boundary conditions at the first/  !!
!!                    last point (natural spline if > 0.99e30)           !!
!! author:                                                               !!
!! ********************************************************************* !!
      SUBROUTINE spline(iat,ich,n,yp1,ypn)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      include 'parameter.h'
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
      common/ hirsh/y2(50,5,150),y(50,5,150),x(150),ieq(maxat),
     1 nrad0,nat0,pop(maxat)

!! special case of H+ atom !!
      if(y(iat,ich,1).ne.0.0d0) then

      if (yp1.gt..99e30) then
        y2(iat,ich,1)=0.d0
        u(1)=0.d0
      else
        y2(iat,ich,1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(iat,ich,2)-y(iat,ich,1))/
     1  (x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(iat,ich,i-1)+2.
        y2(iat,ich,i)=(sig-1.)/p
        u(i)=(6.*((y(iat,ich,i+1)-y(iat,ich,i))/(x(i+1)-x(i))-
     1  (y(iat,ich,i)-y(iat,ich,i-1))/(x(i)-x(i-1)))/(x(i+1)-
     1  x(i-1))-sig*u(i-1))/p
      end do

      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(iat,ich,n)-y(iat,ich,n-1))/
     1  (x(n)-x(n-1)))
      endif

      y2(iat,ich,n)=(un-qn*u(n-1))/(qn*y2(iat,ich,n-1)+1.)

      do k=n-1,1,-1
        y2(iat,ich,k)=y2(iat,ich,k)*y2(iat,ich,k+1)+u(k)
      END DO

      else
      write(*,'(2x,a,2(1x,i0))') 'Set zeroes in spline for species',iat,ich

      do k=n-1,1,-1
        y2(iat,ich,k)=0.0d0                                
      END DO
      end if

      END

!! ***** !!

!! ********************************************************************* !!
!! subroutine: numint_sat                                                !!
!! purpose: integrates the atomic-orbital overlap matrix per atom (sat)  !!
!!   over the numerical grid built by prenumint. runs unconditionally    !!
!!   for every calculation that does not request mulliken/lowdin/nao     !!
!!   analysis (see main.f, "integrate atomic overlap by default") --     !!
!!   ENPART, EOS, OSLO, bond orders and populations all consume sat.     !!
!!   three mutually-exclusive integration schemes, selected by iallpo/   !!
!!   iqtaim (from iopt): same-center only (iallpo=0), full multi-center  !!
!!   (iallpo=1), and qtaim-basin (iqtaim=1).                             !!
!! arguments:                                                            !!
!!   ndim      (in)  -- number of basis functions (leading dim of sat)   !!
!!   itotps    (in)  -- total number of grid points (nat*iatps)          !!
!!   nat0      (in)  -- number of atoms (leading dim of omp2/sat)        !!
!!   wp        (in)  -- integration weight of each grid point            !!
!!   omp       (in)  -- becke/tfvc weight of each point for its own atom !!
!!   omp2      (in)  -- becke/tfvc (or hirshfeld) weight of each point   !!
!!                       for every atom                                  !!
!!   chp       (in)  -- value of each basis function at each grid point  !!
!!   ibaspoint (in)  -- qtaim basin assignment per grid point (qtaim     !!
!!                       branch only)                                    !!
!!   sat       (out) -- per-atom atomic-orbital overlap matrix           !!
!! author: PSalse                                                        !!
!! ********************************************************************* !!
      subroutine numint_sat(ndim,itotps,nat0,wp,omp,omp2,chp,ibaspoint,sat)
      use basis_set, only :s
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      integer, intent(in):: ndim,itotps,nat0
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      dimension wp(itotps),chp(itotps,ndim),omp(itotps)
      dimension ibaspoint(itotps),sat(ndim,ndim,nat0),omp2(itotps,nat0)

!! iopt flags used by this routine !!
      imulli = iopt(5)
      ihirsh = iopt(6)
      iallpo = iopt(7)
      ieffao = iopt(12)
      icube  = iopt(13)
      iqtaim = iopt(16)

      iatps=nrad*nang

      if(iqtaim.eq.1) iallpo=1

!! computing atomic orbital overlap !!
      do mu=1,ndim
        do nu=1,ndim
          do icenter=1,nat
            sat(nu,mu,icenter)=0.0d0               
          end do
        end do
      end do

      if(iallpo.eq.0) then

!! parallel over (icenter,mu): each (icenter,mu,nu) triple accumulates !!
!! its own reduction into the local x and writes only its own          !!
!! sat(mu,nu,icenter)/sat(nu,mu,icenter)                               !!
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,mu,nu,ifut,x)
        do icenter=1,nat
          do mu=1,ndim
            do nu=1,mu
              x=ZERO
              do ifut=iatps*(icenter-1)+1,iatps*icenter
                x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*omp2(ifut,icenter)
              end do
              sat(mu,nu,icenter)=x
              sat(nu,mu,icenter)=x
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
      end if

      if(iallpo.eq.1.and.iqtaim.eq.0) then

!! same reasoning as above -- each triple writes only its own sat entry !!
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,mu,nu,jcenter,jfut,x)
        do icenter=1,nat
          do mu=1,ndim
            do nu=1, mu
              x=0.d0
              do jcenter=1,nat
                do jfut=iatps*(jcenter-1)+1,iatps*jcenter
                  x=x+wp(jfut)*chp(jfut,mu)*chp(jfut,nu)*omp(jfut)*omp2(jfut,icenter)
                end do
              end do
              sat(mu,nu,icenter)=sat(mu,nu,icenter)+x
            end do
          end do
        end do
!$OMP END PARALLEL DO

      else if (iqtaim.eq.1) then

        do jcenter=1,nat
          do jfut=iatps*(jcenter-1)+1,iatps*jcenter
            icenter=ibaspoint(jfut)
            if(icenter.ne.0)then 
              x3=wp(jfut)*omp(jfut)
              do mu=1,ndim
                do nu=1, mu 
                  sat(mu,nu,icenter)=sat(mu,nu,icenter)+chp(jfut,mu)*chp(jfut,nu)*x3
                end do
              end do
            end if
          end do
        end do
      end if 

      do mu=1,ndim
        do nu=1,mu 
          do icenter=1,nat
            sat(nu,mu,icenter)=sat(mu,nu,icenter)
          end do
        end do
      end do

!! sanity check: sum of atomic overlap matrices vs. total overlap S !!
      write(*,*)
      write(*,'(2x,a)') 'Checking sum of AOs overlap matrices'
      write(*,'(2x,a)') '(deviations may not affect overall accuracy of MOs)'
      do i=1,igr
        do j=i,igr
          x=0.0d0
          do iatom=1,nat
            x=x+sat(i,j,iatom)
          end do
          if(abs(s(i,j)-x).gt.1.0d-1) then
            write(*,'(2x,a,2(1x,i0),2(1x,f14.8))') 'Large deviation for element ',i,j,x,s(i,j)
          end if
        end do
      end do

      return 
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: dpoints                                                   !!
!! purpose: evaluates the Laplacian of every basis function at every     !!
!!   numerical-grid point (chp(irun,iact)), contracting each primitive   !!
!!   Gaussian's second derivative over its shell. Same grid/basis loop   !!
!!   structure as fpoints, used by enpart.f/enpart_phf.f's Laplacian-    !!
!!   dependent terms.                                                    !!
!! arguments:                                                            !!
!!   chp    (out) -- Laplacian of each basis function at each grid point !!
!!   pcoord (in)  -- xyz coordinates of each grid point                  !!
!! author: PSalse                                                        !!
!! ********************************************************************* !!
      subroutine dpoints(chp,pcoord)
      use basis_set
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      dimension chp(nrad*nang*natoms,nbasis),pcoord(nrad*nang*natoms,3)

      iatps=nrad*nang

!! parallel over (icenter,ifut): irun computed directly (same convention !!
!! as prenumint), each iteration writes only its own chp(irun,:).        !!
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(icenter,ifut,iact,irun,iactat,x,y,z,
!$OMP&  rr,f,k,ipr,nn,ll,mm,alpha,alpha2,dx2,dy2,dz2,dd)
      do icenter=1,natoms
        do ifut=1,iatps
          irun=(icenter-1)*iatps+ifut
          do iact=1,nbasis
            iactat=ihold(iact)
            x=pcoord(irun,1)-coord(1,iactat)
            y=pcoord(irun,2)-coord(2,iactat)
            z=pcoord(irun,3)-coord(3,iactat)
            rr=x*x+y*y+z*z
            f=0.d0
            k=1
            do while(nprimbas(k,iact).ne.0)
              ipr=nprimbas(k,iact)
              nn=nlm(ipr,1)
              ll=nlm(ipr,2)
              mm=nlm(ipr,3)
              alpha=expp(ipr)
              alpha2=alpha*alpha
              dx2=4.0d0*alpha2*(x**(nn+2))-(2.0d0*alpha*(x**nn)*(2.0*nn+1.0))
              if(nn.ge.2) dx2=dx2+nn*(nn-1)*x**(nn-2)
              dx2=dx2*(y**ll)*(z**mm)
              dy2=4.0d0*alpha2*(y**(ll+2))-(2.0d0*alpha*(y**ll)*(2.0*ll+1.0))
              if(ll.ge.2) dy2=dy2+ll*(ll-1)*y**(ll-2)
              dy2=dy2*(x**nn)*(z**mm)
              dz2=4.0d0*alpha2*(z**(mm+2))-(2.0d0*alpha*(z**mm)*(2.0*mm+1.0))
              if(mm.ge.2) dz2=dz2+mm*(mm-1)*z**(mm-2)
              dz2=dz2*(x**nn)*(y**ll)
              dd=(dx2+dy2+dz2)*dexp(-expp(ipr)*rr)
              f=f+dd*coefpb(ipr,iact)
              k=k+1
            enddo
            chp(irun,iact)=-f
          enddo
        enddo
      end do
!$OMP END PARALLEL DO
      return
      end

!! ********************************************************************* !!
!! subroutine: atdens_int                                                !!
!! purpose: debug/diagnostic utility -- reports the electron density at  !!
!!   one atom's nucleus, then integrates the density (and its angular    !!
!!   anisotropy) over a sphere of radius Rmax centered on that atom, on  !!
!!   its own small Gauss-Legendre/Lebedev grid (independent of the main  !!
!!   integration_grid module). Opt-in via # METHOD / RHO_CALC_AT /       !!
!!   RHO_CALC_RAD (iatdens=0 by default, never called).                  !!
!! arguments:                                                            !!
!!   Rmax    (in) -- integration sphere radius                           !!
!!   iatdens (in) -- index of the atom to integrate around                !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine atdens_int(Rmax,iatdens)
      use basis_set, only: coord
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      allocatable :: pcoord(:,:),wp(:),rho(:)
!! local x/y/z/w/xr/wr, unrelated to integration_grid's module-level     !!
!! arrays of the same name -- this routine never "use"s that module.     !!
      allocatable:: x(:),y(:),z(:),w(:),xr(:),wr(:)

      xxx=functxyz(coord(1,iatdens),coord(2,iatdens),coord(3,iatdens))
      write(*,'(2x,a,1x,i0,a,e16.8)') 'Electron density on coords of atom ',iatdens,' : ',xxx

      nnrad=30
      nnang=110
      allocate(xr(30),wr(30))
      allocate(x(nnang),y(nnang),z(nnang),w(nnang))
      CALL LEGZO(nnrad,XR,WR)
      do i=1,nnrad
        wR(i)=0.50d0*Rmax*wr(i)
        XR(i)=Rmax*(Xr(i)+1.0)/2.0d0
      end do
      CALL LD0110(X,Y,Z,W,nnang)

      allocate(pcoord(nnrad*nnang,3),wp(nnrad*nnang),rho(nnrad*nnang))

!! building pcoords !!
      ifut=0
      icenter=iatdens
      do k=1,nnrad
        rr=xr(k)
        xxr=wr(k)*xr(k)*xr(k)
        do i=1,nnang
          ifut=ifut+1
          pcoord(ifut,1)=coord(1,icenter)+rr*x(i)
          pcoord(ifut,2)=coord(2,icenter)+rr*y(i)
          pcoord(ifut,3)=coord(3,icenter)+rr*z(i)
          rho(ifut)=functxyz(pcoord(ifut,1),pcoord(ifut,2),pcoord(ifut,3))
          wp(ifut)=w(i)*xxr*4.d0*Pi
        enddo
      enddo

      xx=0.0d0
      xanis=0.0d0
      ifut=0
      do k=1,nnrad
        xaver=0.0d0
        xaver2=0.0d0
        do i=1,nnang
          ifut=ifut+1
          xx=xx+rho(ifut)*wp(ifut)
          xaver=xaver+rho(ifut)
          xaver2=xaver2+rho(ifut)**2.0d0
        end do
        xaver=xaver/nnang
        xaver2=xaver2/nnang
        xanis=xanis+xr(k)*xr(k)*wr(k)*(xaver2-xaver*xaver)
      end do
      write(*,'(2x,a42,f6.3,a3,e22.12)') 'Electron density integrated on sphere of R',Rmax,' :',xx
      write(*,'(2x,a43,f6.3,a3,e22.12)') 'Integrated anisotropy of rho on sphere of R',Rmax,' :',xanis

      end

