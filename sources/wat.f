!! ********************************************************************** !!
!! ATOMIC WEIGHT FUNCTIONS -- Becke, Hirshfeld(-I), TFVC                  !!
!! Becke/TFVC weight kernels (thread-safe, called from OMP grid loops     !!
!! elsewhere):                                                            !!
!!   wat             -- normalized atomic weight at a point               !!
!!   pp              -- unnormalized weight (product over cell functions) !!
!!   sbecke          -- Becke/erf cell function                           !!
!! TFVC atomic-radius setup (run once per job, before the grid is used):  !!
!!   khi             -- picks each atom pair's TFVC radius-dividing point !!
!!   prepar          -- Koga empirical atomic radii + distance matrix     !!
!! Hirshfeld / Hirshfeld-Iterative weight evaluators:                     !!
!!   wathirsh        -- standard (non-iterative) Hirshfeld weight         !!
!!   splint          -- cubic-spline lookup on a tabulated radial density !!
!!   wathirsh2       -- Hirshfeld-Iterative weight at converged pop(:)    !!
!!   wathirshit3     -- LIVE self-consistency driver (converges pop(:),   !!
!!                      then rebuilds weights) -- called from numint.f    !!
!!   makeatdens      -- reads/splines the external densoutput file that   !!
!!                      wathirsh/wathirshit3 interpolate against          !!
!!   wathirshit      -- DEAD, zero call sites, superseded by wathirshit3  !!
!!   wathirshit2     -- DEAD, zero call sites, superseded by wathirshit3  !!
!! single-point density/MO evaluators (used by khi's density scan and by  !!
!! print.f's cube generators):                                            !!
!!   functxyz        -- total electron density at a point                !!
!!   dfunctxyz       -- directional derivative of the density at a point  !!
!!   d2functxyz      -- second-derivative density contraction at a point  !!
!!   orbxyz          -- value of one MO at a point                        !!
!! Subroutine Cleanup Protocol: fully done, all 16 routines, 2026-08-18.  !!
!! Known Issues #25 (wathirsh2 arg-count mismatch at its print.f call     !!
!! sites) and #26 (wathirshit/wathirshit2 dead) are open, not fixed --    !!
!! see CLAUDE.md.                                                         !!
!! ********************************************************************** !!

!! ***** !!

!! ********************************************************************* !!
!! function: wat                                                         !!
!! purpose: normalized Becke/TFVC atomic weight of atom ii at point       !!
!! (x,y,z): pp(ii)/sum_j pp(j). Thread-safe (pure function of its         !!
!! arguments plus read-only geometry COMMON) -- called from OMP           !!
!! integration-grid loops in numint.f.                                    !!
!! arguments: ii (atom index), x,y,z (point coordinates)                 !!
!! author:                                                                !!
!! ********************************************************************* !!
      function wat(ii,x,y,z)
      implicit real*8 (a-h,o-z)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop

      sum=0.d0
      do i=1,nat
        p=pp(i,x,y,z)
        sum=sum+p
        if(i.eq.ii)ww=p
      enddo
      wat=ww/sum
      return
      end

!! ********************************************************************* !!
!! function: pp                                                          !!
!! purpose: unnormalized Becke/TFVC weight of atom ii at point (x,y,z) -- !!
!! product over all other atoms j of the cell function sbecke() evaluated !!
!! at the scaled coordinate along the ii-j axis, using the TFVC radius     !!
!! ratio achi(ii,j) (ibcp=1) or the plain atomic-radius ratio otherwise.   !!
!! arguments: ii (atom index), x,y,z (point coordinates)                 !!
!! author:                                                                !!
!! ********************************************************************* !!
      function pp(ii,x,y,z)
      implicit real*8 (a-h,o-z)
      include 'parameter.h'
      common /atomrad/atr(maxat),dist(maxat,maxat)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /achi/achi(maxat,maxat),ibcp

      dimension r(maxat)

      RR(x,y,z)=dsqrt(x**2+y**2+z**2)

      do i=1,nat
        r(i)=RR(x-coord(1,i),y-coord(2,i),z-coord(3,i))
      enddo
      p=1.d0
      do j=1,nat
        if(j.ne.ii) then
          if (ibcp.eq.1) then
            chi=achi(ii,j)
          else
            chi=atr(ii)/atr(j)
          endif
          amu=(r(ii)-r(j))/dist(ii,j)
!! chi is passed explicitly (not via COMMON) so that pp/sbecke stay        !!
!! safe to call from a parallel (OMP) region                               !!
          p=p*sbecke(amu,chi)
        end if
      end do
      pp=p
      return
      end

!! ********************************************************************* !!
!! function: sbecke                                                      !!
!! purpose: Becke cell/boundary function at scaled inter-atomic            !!
!! coordinate amu with radius-ratio parameter chi: the standard iterated   !!
!! polynomial profile (stiffness k=iopt(25), or the shifted/asymmetric     !!
!! form when inewbec=iopt(31)/=0), or an erf-based profile when ierf=1.    !!
!! arguments: amu (scaled coordinate), chi (radius-ratio parameter)      !!
!! author:                                                                !!
!! ********************************************************************* !!
      function sbecke(amu,chi)
      implicit real*8 (a-h,o-z)
      common /iops/iopt(200)
      common /erf/aerf,ierf

!! Becke's profile !!
      p(x)=x*(1.5d0-.5d0*x**2)
!! for erf profile if inewbec=2 and k=1 !!

      k=iopt(25)
      inewbec=iopt(31)

      if(inewbec.eq.0) then
        a=.25d0*(1.d0-chi*chi)/chi
        if(a.gt.0.5d0)a=0.5d0
        if(a.lt.-.5d0)a=-.5d0
        anu=amu+a*(1.d0-amu*amu)
      else
        anu=(1.0d0+amu-chi*(1.0d0-amu))/(1.0d0+amu+chi*(1.0d0-amu))
      end if

      if(ierf.eq.0) then
        do i=1,k
          anu=p(anu)
        enddo
      else
        xx=sqrt(1.0-anu*anu)
        if(abs(xx).gt.1.0d-8) then
          anu=derf(aerf*anu/xx)
        end if
      end if

      sbecke=0.5d0*(1.d0-anu)
      return
      end


!! ********************************************************************* !!
!! subroutine: khi                                                       !!
!! purpose: sets each atom-pair's TFVC radius ratio achi(i,j) by scanning !!
!! the density along the interatomic axis and locating the extremum      !!
!! (minimum, or a maximum when no minimum is consistent) that divides the !!
!! two atoms; short-circuits to a plain atomic-radius ratio when the pair  !!
!! is far apart or (inewbec=1) the midpoint is occupied by a third atom,  !!
!! and can read/write the whole achi matrix from/to radmat.inp            !!
!! (iradmat=iopt(48)=2/1) instead of recomputing it. Prints one summary   !!
!! line per atom pair.                                                    !!
!! arguments: none (nat/coord/achi/atr/dist via COMMON)                  !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine khi
      use basis_set, only : coord
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /atomrad/atr(maxat),dist(maxat,maxat)
      common /achi/achi(maxat,maxat),ibcp
      dimension rho(82),kmin(10),kmax(10)
      character*200 cnote

      inewbec=iopt(31)
      iradmat=iopt(48)

      if(iradmat.eq.2) then 
        open(unit=45,file='radmat.inp')
        read(45,*) nat0
        if(nat0.ne.nat) stop 'inconsistency reading radmat.inp'
        read(45,*)((achi(i,j),j=i+1,nat),i=1,nat) 
        do i=1,nat
          achi(i,i)=0.0d0
          do j=i+1,nat
            achi(j,i)=1.0d0/achi(i,j)
          end do
        end do
        write(*,*)'Atomic radii matrix read from radmat.inp'
        go to 101
      end if

      write(*,'(2x,a)') 'Using density along atom pairs to set atomic radii'
      write(*,*)
!! O(nat^2) atom pairs, each doing a handful of O(igr^2) functxyz() density !!
!! evaluations -- a one-time setup cost, not currently parallelized; a     !!
!! candidate for OMP if it ever shows up as hot on a large system         !!
      do 100 iatom=1,nat-1
        do 100 jatom=iatom+1,nat

          cnote=' '

          if(inewbec.eq.1) then
            rijx2=(coord(1,jatom)+coord(1,iatom))/2.0d0
            rijy2=(coord(2,jatom)+coord(2,iatom))/2.0d0
            rijz2=(coord(3,jatom)+coord(3,iatom))/2.0d0
            disth=dist(iatom,jatom)/2.0d0
            do ii=1,nat
              if(ii.ne.iatom.and.ii.ne.jatom) then
                xdist=sqrt((rijx2-coord(1,ii))**2+(rijy2-coord(2,ii))**2+(rijz2-coord(3,ii))**2)
                if(xdist.lt.disth) then 
                  achi(iatom,jatom)=1.0d0                
                  achi(jatom,iatom)=1.0d0                 
                  goto 100
                end if
              end if
            end do

          else

            if(dist(iatom,jatom).gt.2.0d0*(atr(iatom)+atr(jatom)))then
              achi(iatom,jatom)=atr(iatom)/atr(jatom)
              achi(jatom,iatom)=atr(jatom)/atr(iatom)
              goto 100
            end if

          end if

          rijx=coord(1,jatom)-coord(1,iatom)
          rijy=coord(2,jatom)-coord(2,iatom)
          rijz=coord(3,jatom)-coord(3,iatom)

          ndiv=50

          if(ndiv.gt.80)stop 160
          step=1.d0/(dfloat(ndiv))
      
          do k=1,ndiv+1
            xabs=coord(1,iatom)+step*dfloat(k-1)*rijx
            yabs=coord(2,iatom)+step*dfloat(k-1)*rijy
            zabs=coord(3,iatom)+step*dfloat(k-1)*rijz
            rho(k)=functxyz(xabs,yabs,zabs)
          enddo
c      print *,'pollas',(rho(k),k=1,ndiv+1)

          ncmin=0
          ncmax=0
          maxflg=0
          do k=2,ndiv
            if(rho(k).lt.rho(k-1).and.rho(k).lt.rho(k+1))then 
              ncmin=ncmin+1
              kmin(ncmin)=k
            endif    
            if(rho(k).gt.rho(k-1).and.rho(k).gt.rho(k+1))then 
              ncmax=ncmax+1
              kmax(ncmax)=k
            endif    
          enddo

!! note if there are maxima or no minima along the interatomic axis !!
          if(ncmax.ne.0.or.ncmin.eq.0)then
            write(cnote,'(i0,a,i0,a)') ncmax,' max / ',ncmin,' min extrema'
          endif

!! simplest case: 0 max and 1 min !!
          if(ncmax.eq.0.and.ncmin.eq.1)then
            kc=kmin(1)
            goto 7
          endif
!! flag for pseudopotential !!
          ipseudo=0
          if(int(zn(iatom)).ne.iznuc(iatom).or.int(zn(jatom)).ne.iznuc(jatom))then
            if(len_trim(cnote).gt.0) then
              cnote=trim(cnote)//'; pseudopotential present'
            else
              cnote='pseudopotential present'
            end if
            ipseudo=1
          end if
!! case: there is one more minima than maxima or vice versa !!
          if(ipseudo.eq.0.and.ncmax.ne.0.and.abs(ncmin-ncmax).eq.1) then
            kmc=mod(ncmax,2)
            if(kmc.ne.0) then
              kc=kmax((ncmax+1)/2)
              maxflg=1
              cnote=trim(cnote)//'; using central maximum as dividing point'
            else
              kc=kmin((ncmin+1)/2)
              maxflg=0
              cnote=trim(cnote)//'; using central minimum as dividing point'
            end if
            go to 7
          end if

          if(len_trim(cnote).gt.0) then
            cnote=trim(cnote)//'; inconsistent extrema'
          else
            cnote='inconsistent extrema'
          end if

          ijx=ndiv
          iix=ndiv/2
          do i=1,ncmin
            if(abs(kmin(i)-ndiv/2).lt.ijx) then
              iix=kmin(i)
              ijx=abs(kmin(i)-ndiv/2)
              maxflg=0
            end if
          end do
          do i=1,ncmax
            if(abs(kmax(i)-ndiv/2).lt.ijx) then
              iix=kmax(i)
              ijx=abs(kmax(i)-ndiv/2)
              maxflg=1
            end if
          end do
          kc=iix
          if(maxflg.eq.1) then
            cnote=trim(cnote)//', selecting maximum closest to midpoint'
          else
            cnote=trim(cnote)//', selecting minimum closest to midpoint'
          end if

  7   continue

          nsubdiv=40
          substep=2.d0/dfloat(ndiv*nsubdiv)

          do k=1,nsubdiv+1
            xabs=coord(1,iatom)+(step*dfloat(kc-2)+substep*dfloat(k-1))*rijx
            yabs=coord(2,iatom)+(step*dfloat(kc-2)+substep*dfloat(k-1))*rijy
            zabs=coord(3,iatom)+(step*dfloat(kc-2)+substep*dfloat(k-1))*rijz
            rho(k)=functxyz(xabs,yabs,zabs)
          enddo

          do k=2,nsubdiv
            kc1=k
            if(maxflg.eq.0.and.rho(k).lt.rho(k-1).and.rho(k).lt.rho(k+1))  goto 17
            if(maxflg.eq.1.and.rho(k).gt.rho(k-1).and.rho(k).gt.rho(k+1))  goto 17
          enddo
          if(dist(iatom,jatom).gt.10.0d0) then
            write(*,'(2x,a,i3,a,i3,a)') 'Atom pair ',iatom,' - ',jatom,
     +        ' : no density extremum found, atoms too far apart, using 1:1 ratio'
            achi(iatom,jatom)=1.0d0
            achi(jatom,iatom)=1.0d0
            go to 100
          else
            stop  'No extremum of density found'
          end if

  17  continue
          ri=(step*dfloat(kc-2)+substep*dfloat(kc1-1))
          rj=1.d0-ri
          achi(iatom,jatom)=ri/rj
          achi(jatom,iatom)=rj/ri

          if(len_trim(cnote).gt.0) write(*,'(2x,a,i3,a,i3,a,a)')
     +      'Atom pair ',iatom,' - ',jatom,' : ',trim(adjustl(cnote))
          write(*,'(2x,a,i3,a,i3,a,i0,a,f8.5)') 'Atom pair ',iatom,' - ',
     +      jatom,' : subdivision ',kc1,', ratio ',achi(iatom,jatom)
          if(ri.ge.0.90d0.or.ri.le.0.10) write(*,'(2x,a)')
     +      'WARNING: large polarization for this atom pair, please check'

  100 continue
      if(iradmat.eq.1) then
        open(unit=45,file='radmat.inp')
        write(45,*) nat
        write(45,*)((achi(i,j),j=i+1,nat),i=1,nat)
      end if
  101 continue
      return
      end


!! ********************************************************************* !!
!! function: wathirsh                                                    !!
!! purpose: standard (non-iterative) Hirshfeld weight of atom icenter at  !!
!! point (xabs,yabs,zabs) -- ratio of that atom's tabulated free-atom     !!
!! density (splint() on the radial profile read from densoutput, see      !!
!! makeatdens) to the promolecular sum over all atoms; the sum is         !!
!! skipped (weight set to 0) when the numerator itself is negligible.     !!
!! arguments: icenter (atom index), xabs,yabs,zabs (point coordinates)   !!
!! author:                                                                !!
!! ********************************************************************* !!
      Function wathirsh(icenter,xabs,yabs,zabs)
      use basis_set, only: coord,natoms
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /hirsh/y2(50,5,150),radial(50,5,150),xr2(150),ieq(maxat),nrad2,nat0,pop(maxat)

      toler=1.0d-12

      istate=1
!! istate=1 selects the neutral-atom radial profile (regular Hirshfeld) !!
      ipos=ieq(icenter)

      dist=dsqrt((xabs-coord(1,icenter))**2+(yabs-coord(2,icenter))**2+(zabs-coord(3,icenter))**2)

      x1=splint(ipos,istate, nrad2 , dist)
      if(x1.lt.0) write(*,*) 'error in splint ',x1,ipos,istate,nrad2

      x0=0.0d0
      if(x1.gt.toler) then
        do i=1,natoms
          ipos=ieq(i)
          dist=dsqrt((xabs-coord(1,i))**2+(yabs-coord(2,i))**2+(zabs-coord(3,i))**2)
          x0=x0+splint(ipos,istate, nrad2 , dist)
        end do
        wathirsh=x1/x0
      else
        wathirsh=0.0d0
      end if

      end

!! ********************************************************************* !!
!! function: splint                                                      !!
!! purpose: cubic spline interpolation (Numerical Recipes splint) of a    !!
!! tabulated free-atom radial density profile ya(ipos,istate,:) at        !!
!! y2a(ipos,istate,:) (second derivatives, from spline() in makeatdens)   !!
!! at distance x; binary search locates the bracketing table interval.    !!
!! Clamps a spurious negative result to 0 (interpolation artifact, warned !!
!! above -1.0e-3).                                                        !!
!! arguments: ipos (atom-type index), istate (charge-state index),       !!
!! n (table size), x (evaluation distance)                               !!
!! author:                                                                !!
!! ********************************************************************* !!
      FUNCTION splint(ipos,istate,n,x)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /hirsh/y2a(50,5,150),ya(50,5,150),xa(150),ieq(maxat),nrad2,nat0,pop(maxat)
      INTEGER n
      DOUBLE PRECISION x,y
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n

1     if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
        khi=k
      else
        klo=k
      endif
      goto 1
      endif

      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(ipos,istate,klo)+b*ya(ipos,istate,khi)+((a**3-a)*
     1  y2a(ipos,istate,klo)+(b**3-b)*y2a(ipos,istate,khi))*(h**2)/6.

      if(y.lt.-1.0d-3) write(*,*) 'Large negative
     1  atomic density',y,xa(khi),xa(klo),khi,klo
      if(y.lt.0.0d0) y=0.0d0
      splint=y

      END


!! ********************************************************************* !!
!! function: functxyz                                                    !!
!! purpose: total electron density at point (xabs,yabs,zabs), built by    !!
!! evaluating every AO there and contracting with the density matrix P    !!
!! (rho = sum_mu,nu P(mu,nu) chi_mu chi_nu). Used by khi() to scan the    !!
!! density along an interatomic axis.                                     !!
!! arguments: xabs,yabs,zabs (point coordinates)                         !!
!! author:                                                                !!
!! ********************************************************************* !!
      function functxyz(xabs,yabs,zabs)
      use basis_set
      use ao_matrices
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      allocatable ch(:)

!! this single-point evaluator (and its siblings dfunctxyz/d2functxyz/    !!
!! orbxyz) is itself O(igr^2) but deliberately not OMP-parallelized       !!
!! internally -- like pp()/sbecke(), it must stay callable from an        !!
!! already-parallel caller (e.g. a per-grid-point OMP loop in print.f's   !!
!! cube generators) without nesting a second OMP region                  !!
      allocate (ch(igr))

      do iact=1,igr
        iactat=ihold(iact)
        f=0.d0
        x=xabs-coord(1,iactat)   
        y=yabs-coord(2,iactat)   
        z=zabs-coord(3,iactat)   
        rr=x*x+y*y+z*z
        k=1
        do while(nprimbas(k,iact).ne.0)
          ipr=nprimbas(k,iact)
          nn=nlm(ipr,1)
          ll=nlm(ipr,2)
          mm=nlm(ipr,3)
          f=f+(x**nn)*(y**ll)*(z**mm)*dexp(-expp(ipr)*rr)*coefpb(ipr,iact)
          k=k+1
        enddo
        ch(iact)=f
      end do
      x=0.0d0
      do mu=1,igr
        do nu=1,igr
          x=x+p(mu,nu)*ch(mu)*ch(nu)
        end do
      end do
      deallocate (ch)
      functxyz=x
      end

!! ********************************************************************* !!
!! function: dfunctxyz                                                   !!
!! purpose: directional derivative (d/dixyz) of the electron density at   !!
!! point (xabs,yabs,zabs) -- same AO-contraction structure as functxyz,   !!
!! but also accumulates each AO's derivative along the requested          !!
!! Cartesian direction and contracts P against chi*dchi (off-diagonal     !!
!! terms doubled, matching the mu<nu / mu=mu split used throughout).      !!
!! arguments: ixyz (1/2/3 = d/dx, d/dy, d/dz), xabs,yabs,zabs (point)    !!
!! author:                                                                !!
!! ********************************************************************* !!
      function dfunctxyz(ixyz,xabs,yabs,zabs)
      use basis_set
      use ao_matrices
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      allocatable ch(:),chd(:)

      allocate (ch(igr),chd(igr))
      
      do iact=1,igr
        iactat=ihold(iact)
        x=xabs-coord(1,iactat)   
        y=yabs-coord(2,iactat)   
        z=zabs-coord(3,iactat)   
        rr=x*x+y*y+z*z
        f=0.d0
        fd=0.d0
        k=1
        do while(nprimbas(k,iact).ne.0)
          ipr=nprimbas(k,iact)
          nn=nlm(ipr,1)
          ll=nlm(ipr,2)
          mm=nlm(ipr,3)
          alpha=expp(ipr)
          if(ixyz.eq.1) then
            dx=-2.0d0*alpha*x**(nn+1)
            if(nn.gt.0) dx=dx+nn*x**(nn-1)
            dx=dx*(y**ll)*(z**mm)
          else if(ixyz.eq.2) then
            dx=-2.0d0*alpha*y**(ll+1)
            if(ll.gt.0) dx=dx+ll*y**(ll-1)
            dx=dx*(x**nn)*(z**mm)
          else
            dx=-2.0d0*alpha*z**(mm+1)
            if(mm.gt.0) dx=dx+mm*z**(mm-1)
            dx=dx*(x**nn)*(y**ll)
          end if
          f=f+(x**nn)*(y**ll)*(z**mm)*dexp(-expp(ipr)*rr)*coefpb(ipr,iact)
          fd=fd+dx*dexp(-expp(ipr)*rr)*coefpb(ipr,iact)
          k=k+1
        enddo
        ch(iact)=f
        chd(iact)=fd
      end do
      xx=0.0d0
      do mu=1,igr
        do nu=mu+1,igr
          xx=xx+2.0d0*p(mu,nu)*(ch(mu)*chd(nu)+ch(nu)*chd(mu))
        end do
        xx=xx+2.0d0*p(mu,mu)*ch(mu)*chd(mu)
      end do
      deallocate (ch,chd)
      dfunctxyz=xx
      end


!! ********************************************************************* !!
!! subroutine: wathirshit                                                !!
!! purpose: Hirshfeld-Iterative population/weight update (superseded by   !!
!! wathirshit3, which reuses the same converged-weight logic via          !!
!! wathirsh2() instead of duplicating it inline). Confirmed dead --       !!
!! zero call sites codebase-wide -- left alone pending a deprecation      !!
!! decision, same as old_diagonalize/gennatural_old.                      !!
!! arguments: rho,iatps,wp,whi,nat                                       !!
!! author:                                                                !!
!! ********************************************************************* !!
      SUBROUTINE wathirshit(rho,iatps,wp,whi,nat)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /hirsh/y2(50,5,150),radial(50,5,150),xr2(150),ieq(maxat),nradat,nat0,pop(maxat)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)

      double precision, dimension(maxat) :: znpop
c      double precision, dimension(maxat) :: pop
      double precision, dimension(nat) :: change
      double precision, dimension(nat) :: charge
      double precision, dimension(nat) :: error
      integer, dimension(nat) :: istate
      dimension ipollas(5)
      double precision, dimension(nat*iatps) :: rho
      double precision, dimension(nat*iatps) :: whi
      double precision, dimension(nat*iatps) :: wp

      ipollas(1)=5
      ipollas(2)=3
      ipollas(3)=1
      ipollas(4)=2
      ipollas(5)=4

c iatps: radiale maal angulaire punten

c      WRITE (*,*) "nat", nat
c      WRITE (*,*) "iatps", iatps
c      WRITE (*,*) "nradmol", nradmol
c      WRITE (*,*) "nradat", nradat
c      WRITE (*,*) "nang", nang


c Initialization control matrix
      DO ii=1, nat, 1
        znpop(ii)=zn(ii)
      END DO


      toler=1.0d-8
      nconv=0

c Calcultating atomic populations

      isumcheck=0
      
      DO WHILE (isumcheck<nat.and.nconv.lt.40)
        nconv=nconv+1
        isumcheck=0
        xmaxerr=0.0d0
        DO ii=1, nat, 1
          pop(ii)=0.0d0
          DO ifutt=1,iatps,1

c          write(*,*) 'wp',wp(ifutt+((ii-1)*iatps))
c          write(*,*) 'ifutt',(ifutt+((ii-1)*iatps))
c          write(*,*) 'whi',whi(ifutt+((ii-1)*iatps)),ii,n

            pop(ii)=pop(ii)+wp(ifutt+((ii-1)*iatps))*
     *        rho(ifutt+((ii-1)*iatps))*whi(ifutt+((ii-1)*iatps))
          END DO

c        WRITE(*,*) "pop ii",ii, pop(ii)

          change(ii)=znpop(ii)-pop(ii)

          charge(ii) = zn(ii)-pop(ii)    
        

          error(ii)=abs(znpop(ii)-pop(ii))
c        WRITE(*,*) "change err",ii, change(ii),nconv
c        write(*,*) '##'
c        WRITE(*,*) "CHARGE ",ii, charge(ii)
c        write(*,*) '##'
          znpop(ii)=pop(ii)

c checken op error to stay ithe loop
          IF (error(ii)<0.0005d0) THEN 
            isumcheck=isumcheck+1
          END IF
          if(error(ii).gt.xmaxerr) xmaxerr=error(ii)
        END DO

        write(*,*) 'Iteration :',nconv, 'Max err :',xmaxerr
c       write(*,*) '  Atom  Population   Error'
c       write(*,*) '----------------------------'
        xxp=0.0d0
        do ii=1,nat
          write(*,'(2x,i3,2x,2f10.5)'), ii, znpop(ii),error(ii)
          xxp=xxp+znpop(ii)
        end do
c       write(*,*) '----------------------------'
        write(*,*) 'Tot. Pop. ',xxp 
       

        DO icenter=1, nat,1

          IF (charge(icenter)>0.0d0) THEN
            istate(icenter)=3
            if(charge(icenter)>1.0d0) then
              istate(icenter)=5
            else if(charge(icenter)>2.0d0) then
              write(*,*) 'Atomic charge exceeds 2.0 for atom ',icenter
            end if
          ELSE IF(charge(icenter)<0.0d0) then
            istate(icenter)=2
            if(charge(icenter)<-1.0d0) then
              istate(icenter)=4
            else if(charge(icenter)<-2.0d0) then
              write(*,*) 'Atomic charge exceeds -2.0 for atom ',icenter
            end if
          ELSE
            istate(icenter)=1
          END IF

c           WRITE(*,*) "icenter,istate",icenter,istate(icenter)

        END DO

        

        DO icenter=1, nat,1

          ifuty=0


          do k =1,nrad
            do nn=1,nang

              ilint=INT(pop(icenter))-int(zn(icenter))+3
              iuint=INT(pop(icenter))+1-int(zn(icenter))+3
c            write(*,*) 'Using ',icenter,ipollas(iuint),ipollas(ilint)
              ifuty=ifuty+1
              thx=th(nn)
              fix=ph(nn)
              rr=xr(k)
              xabs0=rr*dsin(thx)*dcos(fix)
              yabs0=rr*dsin(thx)*dsin(fix)
              zabs0=rr*dcos(thx)
              xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
              yabs=-dcos(pha)*dsin(phb)*xabs0+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
              zabs=dsin(pha)*dsin(phb)*xabs0-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0
              xabs=xabs+coord(1,icenter)
              yabs=yabs+coord(2,icenter)
              zabs=zabs+coord(3,icenter)


c Recalculate the hirshfeld weights

      
              ipos=ieq(icenter)
              dist=dsqrt((xabs-coord(1,icenter))**2+
     *          (yabs-coord(2,icenter))**2
     *          +(zabs-coord(3,icenter))**2)
      
        

              x1=(pop(icenter)-INT(pop(icenter)))*splint(ipos,ipollas(iuint)
     *          , nradat , dist)+
     *          ((INT(pop(icenter))+1)-pop(icenter))*splint(ipos,ipollas(ilint)
     *          , nradat,dist) 
!min een

c        WRITE(*,*) 'x1',pop(icenter)-INT(pop(icenter)),
c     *(INT(pop(icenter))+1)-pop(icenter), iuint, ilint

              x0=0.0d0  

              if(x1.gt.toler) then
                do i=1,nat
                  ipos=ieq(i)
                  dist=dsqrt((xabs-coord(1,i))**2+(yabs-coord(2,i))**2
     *              +(zabs-coord(3,i))**2)
                  ilint=INT(pop(i))-int(zn(i))+3
                  iuint=INT(pop(i))+1-int(zn(i))+3
c screening of the point
                  x0=x0+(pop(i)-INT(pop(i)))*
     *              splint(ipos,ipollas(iuint), nradat , dist)+
     *              ((INT(pop(i))+1)-pop(i))*splint(ipos,ipollas(ilint)
     *              , nradat,dist) 

                end do


c        write(*,*)"icenter...",ifuty
c        write(*,*)"icenter...",icenter, k, nn

                whi(((icenter-1)*iatps)+ifuty)=x1/x0
 
              else
                whi(((icenter-1)*iatps)+ifuty)=0.0d0
c        write(*,*)"icenterelse...",ifuty
              end if          
c        write(*,*)'wathirsh',whi(((icenter-1)*iatps)+ifuty),dist,n
        

            END DO
          END DO
        END DO

      END DO

      write(*,*) 'Populations Converged '
      write(*,*) '  Atom    Population   Error    '
      write(*,*) '---------------------------------'
      do ii=1,nat
        write(*,'(2x,i3,2x,2f10.5)'), ii, znpop(ii),error(ii)
      end do
      write(*,*) '---------------------------------'

      end

!! ********************************************************************* !!
!! subroutine: wathirshit2                                               !!
!! purpose: Hirshfeld-Iterative population/weight update, cross-atom      !!
!! whi(:,nat) variant with charge-state clamping (superseded by           !!
!! wathirshit3). Confirmed dead -- zero call sites codebase-wide -- left  !!
!! alone pending a deprecation decision, same as old_diagonalize/          !!
!! gennatural_old.                                                        !!
!! arguments: rho,iatps,wp,whi,nat                                       !!
!! author:                                                                !!
!! ********************************************************************* !!
      SUBROUTINE wathirshit2(rho,iatps,wp,whi,nat)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /hirsh/y2(50,5,150),radial(50,5,150),xr2(150),ieq(maxat),nradat,nat0,pop(maxat)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)

      double precision, dimension(maxat) :: znpop
c      double precision, dimension(maxat) :: pop
      double precision, dimension(nat) :: change
      double precision, dimension(nat) :: charge
      double precision, dimension(nat) :: error
      integer, dimension(nat) :: istate
      dimension ipollas(5)
      double precision, dimension(nat*iatps) :: rho
      double precision, dimension(nat*iatps,nat) :: whi
      double precision, dimension(nat*iatps) :: wp

      ipollas(1)=5
      ipollas(2)=3
      ipollas(3)=1
      ipollas(4)=2
      ipollas(5)=4
 

c iatps: radiale maal angulaire punten


c      WRITE (*,*) "nat", nat
c      WRITE (*,*) "iatps", iatps
c      WRITE (*,*) "nradmol", nradmol
c      WRITE (*,*) "nradat", nradat
c      WRITE (*,*) "nang", nang


c Initialization control matrix
      DO ii=1, nat, 1
        znpop(ii)=zn(ii)
      END DO


      toler=1.0d-8
      nconv=0

c Calcultating atomic populations

      isumcheck=0
      DO WHILE (isumcheck<nat.and.nconv.lt.40)
        nconv=nconv+1
        isumcheck=0
        xmaxerr=0.0d0
        DO ii=1, nat, 1
          pop(ii)=0.0d0
          DO ifutt=1,iatps,1

c          write(*,*) 'wp',wp(ifutt+((ii-1)*iatps))
c          write(*,*) 'ifutt',(ifutt+((ii-1)*iatps))
c          write(*,*) 'whi',whi(ifutt+((ii-1)*iatps)),ii,n

            pop(ii)=pop(ii)+wp(ifutt+((ii-1)*iatps))*
     *        rho(ifutt+((ii-1)*iatps))*whi(ifutt+((ii-1)*iatps),ii)
          END DO

c        WRITE(*,*) "pop ii",ii, pop(ii)

          change(ii)=znpop(ii)-pop(ii)

          charge(ii) = zn(ii)-pop(ii)    
        

          error(ii)=abs(znpop(ii)-pop(ii))
c        WRITE(*,*) "change err",ii, change(ii),nconv
c        write(*,*) '##'
c        WRITE(*,*) "CHARGE ",ii, charge(ii)
c        write(*,*) '##'
          znpop(ii)=pop(ii)

c checken op error to stay ithe loop
          IF (error(ii)<0.0005d0) THEN 
            isumcheck=isumcheck+1
          END IF
          if(error(ii).gt.xmaxerr) xmaxerr=error(ii)
        END DO

        write(*,*) 'Iteration :',nconv,' Max err: ',xmaxerr
c       write(*,*) '  Atom  Population   Error'
c       write(*,*) '----------------------------'
        xxp=0.0d0
        do ii=1,nat
          write(*,'(2x,i3,2x,2f10.5)'), ii, znpop(ii),error(ii)
          xxp=xxp+znpop(ii)
        end do
c       write(*,*) '----------------------------'
        write(*,*) 'Tot. Pop. ',xxp 
       

        DO icenter=1, nat,1

          IF (charge(icenter)>0.0d0) THEN
            istate(icenter)=3
            if(charge(icenter)>1.0d0) then
              istate(icenter)=5
            else if(charge(icenter)>2.0d0) then
              write(*,*) 'Atomic charge exceeds 2.0 for atom ',icenter
            end if
          ELSE IF(charge(icenter)<0.0d0) then
            istate(icenter)=2
            if(charge(icenter)<-1.0d0) then
              istate(icenter)=4
            else if(charge(icenter)<-2.0d0) then
              write(*,*) 'Atomic charge exceeds -2.0 for atom ',icenter
            end if
          ELSE
            istate(icenter)=1
          END IF

c           WRITE(*,*) "icenter,istate",icenter,istate(icenter)

        END DO

        

        DO icenter=1, nat,1

          ifuty=0


          do k =1,nrad
            do nn=1,nang

              ilint=INT(pop(icenter))-int(zn(icenter))+3
              iuint=INT(pop(icenter))+1-int(zn(icenter))+3
              if(ilint.le.0) then
                ilint=1
                iuint=2
              else if(iuint.ge.6) then
                ilint=4
                iuint=5
              end if
c            write(*,*) 'Using ',icenter,ipollas(iuint),ipollas(ilint)
              ifuty=ifuty+1
              thx=th(nn)
              fix=ph(nn)
              rr=xr(k)
              xabs0=rr*dsin(thx)*dcos(fix)
              yabs0=rr*dsin(thx)*dsin(fix)
              zabs0=rr*dcos(thx)
              xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
              yabs=-dcos(pha)*dsin(phb)*xabs0+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
              zabs=dsin(pha)*dsin(phb)*xabs0-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0
              xabs=xabs+coord(1,icenter)
              yabs=yabs+coord(2,icenter)
              zabs=zabs+coord(3,icenter)


c Recalculate the hirshfeld weights

      
              ipos=ieq(icenter)
              dist=dsqrt((xabs-coord(1,icenter))**2+
     *          (yabs-coord(2,icenter))**2
     *          +(zabs-coord(3,icenter))**2)
      
        

              x1=(pop(icenter)-INT(pop(icenter)))*splint(ipos,ipollas(iuint)
     *          , nradat , dist)+
     *          ((INT(pop(icenter))+1)-pop(icenter))*splint(ipos,ipollas(ilint)
     *          , nradat,dist) 
!min een

c        WRITE(*,*) 'x1',pop(icenter)-INT(pop(icenter)),
c     *(INT(pop(icenter))+1)-pop(icenter), iuint, ilint

              x0=0.0d0  

              if(x1.gt.toler) then
                do i=1,nat
                  ipos=ieq(i)
                  dist=dsqrt((xabs-coord(1,i))**2+(yabs-coord(2,i))**2
     *              +(zabs-coord(3,i))**2)
                  ilint=INT(pop(i))-int(zn(i))+3
                  iuint=INT(pop(i))+1-int(zn(i))+3
                  if(ilint.le.0) then
                    ilint=1
                    iuint=2
                  else if(iuint.ge.6) then
                    ilint=4
                    iuint=5
                  end if
c screening of the point
                  x0=x0+(pop(i)-INT(pop(i)))*
     *              splint(ipos,ipollas(iuint), nradat , dist)+
     *              ((INT(pop(i))+1)-pop(i))*splint(ipos,ipollas(ilint)
     *              , nradat,dist) 

                end do


c        write(*,*)"icenter...",ifuty
c        write(*,*)"icenter...",icenter, k, nn

                whi(((icenter-1)*iatps)+ifuty,icenter)=x1/x0
 
              else
                whi(((icenter-1)*iatps)+ifuty,icenter)=0.0d0
c        write(*,*)"icenterelse...",ifuty
              end if          
c        write(*,*)'wathirsh',whi(((icenter-1)*iatps)+ifuty),dist,n
        

            END DO
          END DO
        END DO

      END DO

      write(*,*) 'Iterative Hirshfeld Populations Converged '
      write(*,*) '  Atom       Charge     '
      write(*,*) '------------------------'
      do ii=1,nat
        write(*,'(2x,i3,2x,f10.5)'), ii, zn(ii)-znpop(ii)    
      end do
      write(*,*) '------------------------'

c calcualtion weights in all atoms
      DO icenter=1, nat,1
        ifuty=0
        do k =1,nrad
          do nn=1,nang

c            write(*,*) 'Using ',icenter,ipollas(iuint),ipollas(ilint)
            ifuty=ifuty+1
            thx=th(nn)
            fix=ph(nn)
            rr=xr(k)
            xabs0=rr*dsin(thx)*dcos(fix)
            yabs0=rr*dsin(thx)*dsin(fix)
            zabs0=rr*dcos(thx)
            xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
            yabs=-dcos(pha)*dsin(phb)*xabs0+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
            zabs=dsin(pha)*dsin(phb)*xabs0-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0
            xabs=xabs+coord(1,icenter)
            yabs=yabs+coord(2,icenter)
            zabs=zabs+coord(3,icenter)

            x0=0.0d0  
            do i=1,nat
              ilint=INT(pop(i))-int(zn(i))+3
              iuint=INT(pop(i))+1-int(zn(i))+3
              ipos=ieq(i)
              dist=dsqrt((xabs-coord(1,i))**2+(yabs-coord(2,i))**2
     *          +(zabs-coord(3,i))**2)
              ilint=INT(pop(i))-int(zn(i))+3
              iuint=INT(pop(i))+1-int(zn(i))+3
              x0=x0+(pop(i)-INT(pop(i)))*
     *          splint(ipos,ipollas(iuint), nradat , dist)+
     *          ((INT(pop(i))+1)-pop(i))*splint(ipos,ipollas(ilint)
     *          , nradat,dist) 
            end do

            if(x0.gt.toler) then

              do jcenter=1,nat
                ilint=INT(pop(jcenter))-int(zn(jcenter))+3
                iuint=INT(pop(jcenter))+1-int(zn(jcenter))+3
                ipos=ieq(jcenter)
                dist=dsqrt((xabs-coord(1,jcenter))**2+
     *            (yabs-coord(2,jcenter))**2
     *            +(zabs-coord(3,jcenter))**2)
        
                x1=(pop(jcenter)-INT(pop(jcenter)))*splint(ipos,ipollas(iuint)
     *            , nradat , dist)+
     *            ((INT(pop(jcenter))+1)-pop(jcenter))*splint(ipos,ipollas(ilint)
     *            , nradat,dist) 

                whi(((icenter-1)*iatps)+ifuty,jcenter)=x1/x0
              end do
 
            else
              do jcenter=1,nat
                whi(((icenter-1)*iatps)+ifuty,jcenter)=0.0d0
              end do
            end if          

          END DO
        END DO

      END DO


      end

!! ********************************************************************* !!
!! function: wathirsh2                                                   !!
!! purpose: Hirshfeld-Iterative weight of atom jcenter at point           !!
!! (xabs,yabs,zabs) using the current converged fractional-charge-state   !!
!! populations pop(:) (COMMON /hirsh/) -- interpolates each atom's        !!
!! radial density between its floor/ceiling integer charge states via     !!
!! splint(), same promolecular-ratio structure as wathirsh but with a      !!
!! population-dependent profile per atom instead of a fixed neutral one.  !!
!! arguments: jcenter (atom index), xabs,yabs,zabs (point coordinates)   !!
!! author:                                                                !!
!! ********************************************************************* !!
      function wathirsh2(jcenter,xabs,yabs,zabs)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /hirsh/y2(50,5,150),radial(50,5,150),xr2(150),ieq(maxat),nradat,nat0,pop(maxat)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)

      dimension ipollas(0:6)

      ipollas(0)=5
      ipollas(1)=5
      ipollas(2)=3
      ipollas(3)=1
      ipollas(4)=2
      ipollas(5)=4
      ipollas(6)=4
 
      toler=1.0d-12

      x0=0.0d0  
      do i=1,nat
        ilint=INT(pop(i))-int(zn(i))+3
        iuint=INT(pop(i))+1-int(zn(i))+3
        ipos=ieq(i)
        dist=dsqrt((xabs-coord(1,i))**2+(yabs-coord(2,i))**2
     *    +(zabs-coord(3,i))**2)
        x0=x0+(pop(i)-INT(pop(i)))*
     *    splint(ipos,ipollas(iuint), nradat , dist)+
     *    ((INT(pop(i))+1)-pop(i))*splint(ipos,ipollas(ilint)
     *    , nradat,dist) 
      end do

      if(x0.gt.toler) then

        ilint=INT(pop(jcenter))-int(zn(jcenter))+3
        iuint=INT(pop(jcenter))+1-int(zn(jcenter))+3
        ipos=ieq(jcenter)
        dist=dsqrt((xabs-coord(1,jcenter))**2+
     *    (yabs-coord(2,jcenter))**2
     *    +(zabs-coord(3,jcenter))**2)
        
        x1=(pop(jcenter)-INT(pop(jcenter)))*splint(ipos,ipollas(iuint)
     *    , nradat , dist)+
     *    ((INT(pop(jcenter))+1)-pop(jcenter))*splint(ipos,ipollas(ilint)
     *    , nradat,dist) 

        wathirsh2=x1/x0
 
      else
        wathirsh2=0.d0 
      end if          


      end

!! ********************************************************************* !!
!! subroutine: wathirshit3                                               !!
!! purpose: Hirshfeld-Iterative population cycle -- when iiter=1 (and not !!
!! reading a cached radmat.inp), iterates atomic populations pop(:) to    !!
!! self-consistency (each pop(i) driven by the current per-atom weights   !!
!! whi via wathirsh2()) until every atom's population change drops below  !!
!! tolerance, then recomputes whi at the converged populations; can       !!
!! persist/restore pop(:) via hirshit.pop (iopt(48)=1/2) to skip           !!
!! reconverging on a repeat run. iiter=0 just (re)builds whi from the      !!
!! already-current pop(:) with no iteration.                              !!
!! arguments: rho,iatps,wp,whi,nat,iiter                                 !!
!! author:                                                                !!
!! ********************************************************************* !!
      SUBROUTINE wathirshit3(rho,iatps,wp,whi,nat,iiter)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /iops/iopt(200)
      common /hirsh/y2(50,5,150),radial(50,5,150),xr2(150),ieq(maxat),nradat,nat0,pop(maxat)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)

      double precision, dimension(maxat) :: znpop
      double precision, dimension(nat) :: change
      double precision, dimension(nat) :: charge
      double precision, dimension(nat) :: error
      dimension ipollas(0:6)
      double precision, dimension(nat*iatps) :: rho
      double precision, dimension(nat*iatps,nat) :: whi
      double precision, dimension(nat*iatps) :: wp

      ipollas(0)=5
      ipollas(1)=5
      ipollas(2)=3
      ipollas(3)=1
      ipollas(4)=2
      ipollas(5)=4
      ipollas(6)=4
      
      iradmat= iopt(48)
      
      if(iiter.eq.0.or.iradmat.eq.2) go to 1

!! initialize the population-tracking array with the neutral-atom charges !!
      DO ii=1, nat
        znpop(ii)=zn(ii)
      END DO
      toler=1.0d-8
      nconv=0

!! self-consistency loop: recompute atomic populations from the current   !!
!! weights, then recompute the weights from the new populations          !!
      isumcheck=0
      DO WHILE (isumcheck<nat.and.nconv.lt.40)
        nconv=nconv+1
        isumcheck=0
        xmaxerr=0.0d0
        do ii=1,nat
          xx=0.0d0
          do ifut=iatps*(ii-1)+1,iatps*ii     
            xx=xx+wp(ifut)*rho(ifut)*whi(ifut,ii)
          END DO
          pop(ii)=xx    

          if(zn(ii)-pop(ii).ge.2.0d0) then
            write(*,*) 'Atomic charge exceeds 2.0 for atom ',ii
            pop(ii)=zn(ii)-2.0d0
          ELSE if(zn(ii)-pop(ii).le.-2.0d0) then
            write(*,*) 'Atomic charge exceeds -2.0 for atom ',ii
            pop(ii)=zn(ii)+2.0d0
          end if
        
        END DO

!! check convergence: an atom counts as converged once its population    !!
!! change from the previous iteration drops below tolerance              !!
        do ii=1,nat
          change(ii)=znpop(ii)-pop(ii)
          charge(ii) = zn(ii)-pop(ii)
          error(ii)=abs(znpop(ii)-pop(ii))
          znpop(ii)=pop(ii)
          IF (error(ii)<0.0005d0) isumcheck=isumcheck+1
          if(error(ii).eq.0.0d0) STOP 'Absolute convergence unlikely.'
          if(error(ii).gt.xmaxerr) xmaxerr=error(ii)
        END DO
        write(*,*) 'Iteration :',nconv,' Max err: ',xmaxerr
        xxp=0.0d0
        do ii=1,nat
          write(*,'(2x,i3,2x,2f10.5)') ii, znpop(ii),error(ii)
          xxp=xxp+znpop(ii)
        end do
        write(*,*) 'Tot. Pop. ',xxp



        ifut=0
        DO icenter=1, nat
          do k =1,nrad
            do nn=1,nang
              ifut=ifut+1
              thx=th(nn)
              fix=ph(nn)
              rr=xr(k)
              xabs0=rr*dsin(thx)*dcos(fix)
              yabs0=rr*dsin(thx)*dsin(fix)
              zabs0=rr*dcos(thx)
              xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
              yabs=-dcos(pha)*dsin(phb)*xabs0+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
              zabs=dsin(pha)*dsin(phb)*xabs0-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0
              xabs=xabs+coord(1,icenter)
              yabs=yabs+coord(2,icenter)
              zabs=zabs+coord(3,icenter)
              whi(ifut,icenter)=wathirsh2(icenter,xabs,yabs,zabs) 
            END DO
          END DO
        END DO

      END DO
!! end of self-consistency loop !!

      write(*,*) 'Iterative Hirshfeld Populations Converged '
      write(*,*) '  Atom       Charge     '
      write(*,*) '------------------------'
      do ii=1,nat
        write(*,'(2x,i3,2x,f10.5)') ii, zn(ii)-znpop(ii)
      end do
      write(*,*) '------------------------'

      if(iradmat.eq.1) then
        open(unit=45,file='hirshit.pop')
        write(45,*) nat
        write(45,*) (pop(i),i=1,nat)
      end if

1      continue
      if(iradmat.eq.2.and.iiter.eq.1) then
        open(unit=45,file='hirshit.pop')
        read(45,*) nat0
        if(nat0.ne.nat) stop 'inconsistency reading hirshit.pop'
        read(45,*) (pop(i),i=1,nat)
        write(*,*)'Hirshfeld-Iterative atomic populations read from hirshit.pop'
        close(45)
      end if


!! recompute the weights for every atom at the converged populations.    !!
!! independent per grid point/atom pair -- a candidate for OMP, but this  !!
!! whole Hirshfeld-Iterative path has zero test coverage (no active test  !!
!! exercises HIRSH-IT), so left serial for now, same reasoning as         !!
!! prenumint's untested Hirshfeld branch (see numint.f)                   !!
      ifut=0
      DO icenter=1, nat
        do k =1,nrad
          do nn=1,nang
            ifut=ifut+1
            thx=th(nn)
            fix=ph(nn)
            rr=xr(k)
            xabs0=rr*dsin(thx)*dcos(fix)
            yabs0=rr*dsin(thx)*dsin(fix)
            zabs0=rr*dcos(thx)
            xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
            yabs=-dcos(pha)*dsin(phb)*xabs0+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
            zabs=dsin(pha)*dsin(phb)*xabs0-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0
            xabs=xabs+coord(1,icenter)
            yabs=yabs+coord(2,icenter)
            zabs=zabs+coord(3,icenter)
            do jcenter=1,nat
              whi(ifut,jcenter)=wathirsh2(jcenter,xabs,yabs,zabs) 
            END DO
          END DO
        END DO
      END DO

      end

!! ********************************************************************* !!
!! subroutine: prepar                                                    !!
!! purpose: builds the per-atom Becke/TFVC radius table atr(:) (Koga      !!
!! empirical radii, overridable per-run by an optional radius.inp file)   !!
!! and the full interatomic distance matrix dist(:,:).                    !!
!! arguments: none (nat/coord/iznuc via COMMON, output via atr/dist)     !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine prepar
      implicit real*8 (a-h,o-z)
      include 'parameter.h'
      dimension atrad(92)
      common /atomrad/atr(maxat),dist(maxat,maxat)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
!! radii from Koga; original value for F is 0.633. rare-gas data from  !!
!! Wikipedia                                                              !!
      data atrad/0.327, 0.320, 1.219, 0.911, 0.793, 0.766, 0.699, 0.658,
     +  0.900, 0.690, 1.545, 1.333, 1.199, 1.123, 1.110, 1.071, 1.039,
     +  0.970, 1.978, 1.745, 1.337, 1.274, 1.236, 1.128, 1.180, 1.091,
     +  1.089, 1.077, 1.146, 1.187, 1.199, 1.179, 1.209, 1.201, 1.201,
     +  1.100, 2.217, 1.928, 1.482, 1.377, 1.353, 1.240, 1.287, 1.212,
     +  1.229, 1.240, 1.362, 1.429, 1.385, 1.380, 1.421, 1.400, 1.397,
     +  1.300, 2.442, 2.149, 1.653, 1.600, 1.600, 1.600, 1.600, 1.600,
     +  1.600, 1.600, 1.500, 1.500, 1.500, 1.500, 1.500, 1.500, 1.500,
     +  1.364, 1.346, 1.256, 1.258, 1.222, 1.227, 1.227, 1.273, 1.465,
     +  1.531, 1.434, 1.496, 1.500, 1.500, 1.450, 1.500, 1.500, 1.500,
     +  1.500, 1.500,1.500/ 
      RR(x,y,z)=dsqrt(x**2+y**2+z**2)

!! optional per-run override: read empirical radii from radius.inp if     !!
!! present (MG: Hungarian original read "fill in atomic radii sensibly") !!
      open(unit=20,file="radius.inp",STATUS="OLD",err=100)
      read(20,*,end=100) nn
      write(*,*)
      write(*,'(2x,a)') 'Empirical radii read from radius.inp'
      do i=1,nn
        read(20,*) inn,xnn
        write(*,'(2x,a,1x,i3,2x,a,1x,f6.3)') 'Atomic number:',inn,'radius:',xnn
        atrad(inn)=xnn
      end do
      write(*,*)
 100  continue

      do inn=1,92
        atrad(inn)=atrad(inn)/0.52917721067d0
      end do

      do i=1,nat
        ia=iznuc(i)
        atr(i)=atrad(ia)
      enddo

!! O(nat^2) but a one-time setup cost -- not worth OMP !!
      do i=1,nat-1
        dist(i,i)=0.d0
        do j=i+1,nat
          dist(i,j)=RR(coord(1,i)-coord(1,j),coord(2,i)-coord(2,j),coord(3,i)-coord(3,j))
          dist(j,i)=dist(i,j)
        enddo
      enddo
      return
      end
      
**********************************************************************

!! ********************************************************************* !!
!! subroutine: makeatdens                                                !!
!! purpose: reads the tabulated free-atom radial density profiles (for   !!
!! every charge state of every atom type present) from an external        !!
!! densoutput file, matches each entry to the atoms present in this run   !!
!! via ieq(:), and splines each profile (spline(), for later splint()      !!
!! evaluation by wathirsh/wathirsh2). Feeds the Hirshfeld-family AIM       !!
!! schemes (HIRSH/HIRSH-IT).                                              !!
!! arguments: none (nat/coord/iznuc via COMMON, output via ieq/radial/y2) !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine makeatdens
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /hirsh/y2(50,5,150),radial(50,5,150),xr2(150),ieq(maxat),nrad0,nat0,pop(maxat)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      character*2 atname
      character*4 atname2
      character*80 line    

      Dimension mend(92)
      data mend/4H  H ,4H He ,4H Li ,4H Be ,4H  B ,4H  C ,4H  N ,4H  O , 
     $  4H  F ,4H Ne ,4H Na ,4H Mg ,4H Al ,4H Si ,4H  P ,4H  S ,4H Cl ,
     $  4H Ar ,4H  K ,4H Ca ,4H Sc ,4H Ti ,4H  V ,4H Cr ,4H Mn ,4H Fe ,
     $  4H Co ,4H Ni ,4H Cu ,4H Zn ,4H Ga ,4H Ge ,4H As ,4H Se ,4H Br ,
     $  4H Kr ,4H Rb ,4H Sr ,4H  Y ,4H Zr ,4H Nb ,4H Mo ,4H Tc ,4H Ru ,
     $  4H Rh ,4H Pd ,4H Ag ,4H Cd ,4H In ,4H Sn ,4H Sb ,4H Te ,4H  I ,
     $  4H Xe ,4H Cs ,4H Ba ,4H La ,4H Ce ,4H Pr ,4H Nd ,4H Pm ,4H Sn ,
     $  4H Eu ,4H Gd ,4H Tb ,4H Dy ,4H Ho ,4H Er ,4H Tm ,4H Yb ,4H Lu , 
     $  4H Hf ,4H Ta ,4H  W ,4H Re ,4H Os ,4H Ir ,4H Pt ,4H Au ,4H Hg ,
     $  4H Tl ,4H Pb ,4H Bi ,4H Po ,4H At ,4H Rn ,4H Fr ,4H Ra ,4H Ac ,
     $  4H Th ,4H Pa ,4H  U   /
    
      write(*,'(2x,a)') 'Reading atomic densities from densoutput'

      OPEN (UNIT=51, FILE="densoutput")
      read(51,'(a80)') line
      read(51,*) nat0,nrad0
      if(real(size(xr2)).lt.nrad0) STOP 'Array xr2 is too small to read
     +  densoutput.'
      READ (51,*) (xr2(ii),ii=1,nrad0)

      DO iat=1, nat0, 1
        read(51,'(a2)')atname
        atname=adjustl(atname)
        read(51,*) nch 
        do j=1,92
          write(atname2,'(a4)') mend(j)
          if(index(atname2,atname).ne.0) then
            write(*,'(2x,a,1x,a2,a)') 'Matched element',atname,' in densoutput'
            ipos=j
            do i=1,nat
              if(iznuc(i).eq.ipos) then
                ieq(i)=iat
              end if
            end do
            go to 10
          end if
        end do
10      continue

        DO ich=1, nch, 1
          READ (51,*)  icharge,imult
          READ (51,*) (radial(iat, ich, io),io=1,nrad0) 
          CALL spline(iat,ich,nrad0,-10.0d0,0.0d0)
        END DO
      END DO
 
      CLOSE(51)

      do i=1,nat
        if(ieq(i).eq.0) then
          write(*,'(2x,a,a4,a)') 'ERROR: atom ', mend(iznuc(i)),
     1      ' is missing in densoutput'
          stop
        endif
      end do

      end
!! ********************************************************************* !!
!! function: orbxyz                                                      !!
!! purpose: value of MO nu (from coefficient matrix A) at point           !!
!! (xabs,yabs,zabs) -- evaluates every AO there and contracts with A's    !!
!! nu-th column. Used by the cube-file generators (print.f) to plot an    !!
!! individual orbital rather than a density.                              !!
!! arguments: A (MO coefficient matrix), nu (orbital index), xabs,yabs,  !!
!! zabs (point coordinates)                                              !!
!! author:                                                                !!
!! ********************************************************************* !!
      function orbxyz(A,nu,xabs,yabs,zabs)
      use basis_set
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      dimension a(igr,igr)
      allocatable ch(:)

      allocate (ch(igr))

      do iact=1,igr
        iactat=ihold(iact)
        f=0.d0
        x=xabs-coord(1,iactat)
        y=yabs-coord(2,iactat)
        z=zabs-coord(3,iactat)
        rr=dsqrt(x**2+y**2+z**2)
        k=1
        do while(nprimbas(k,iact).ne.0)
          ipr=nprimbas(k,iact)
          nn=nlm(ipr,1)
          ll=nlm(ipr,2)
          mm=nlm(ipr,3)
          f=f+(x**nn)*(y**ll)*(z**mm)*dexp(-expp(ipr)*(rr**2))*coefpb(ipr,iact)
          k=k+1
        enddo
        ch(iact)=f
      end do

      x=0.0d0
      do mu=1,igr
        x=x+a(mu,nu)*ch(mu)
      end do
      orbxyz=x
      deallocate (ch)
      end


!! ********************************************************************* !!
!! function: d2functxyz                                                  !!
!! purpose: second-derivative contraction of the density at point         !!
!! (xabs,yabs,zabs) along the (ixyz,jxyz) Cartesian pair (dx2/dy2/dz2/    !!
!! dxdy/dxdz/dydz, icode=ixyz*jxyz selects which) -- same AO-contraction   !!
!! structure as functxyz/dfunctxyz, extended to first and second AO       !!
!! derivatives.                                                           !!
!! arguments: ixyz,jxyz (1/2/3 = x/y/z, selects the derivative pair),    !!
!! xabs,yabs,zabs (point coordinates)                                    !!
!! author:                                                                !!
!! ********************************************************************* !!
      function d2functxyz(ixyz,jxyz,xabs,yabs,zabs)
      use basis_set
      use ao_matrices
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      allocatable ch(:),chd(:,:)

      allocate (ch(igr),chd(igr,3))

      icode=ixyz*jxyz
      
      do iact=1,igr
        iactat=ihold(iact)
        f=0.d0
        fdx=0.d0
        fdy=0.d0
        fdz=0.d0
        fd2=0.d0
        x=xabs-coord(1,iactat)   
        y=yabs-coord(2,iactat)   
        z=zabs-coord(3,iactat)   
        rr=x*x+y*y+z*z
        k=1
        do while(nprimbas(k,iact).ne.0)
          ipr=nprimbas(k,iact)
          nn=nlm(ipr,1)
          ll=nlm(ipr,2)
          mm=nlm(ipr,3)
          f=f+(x**nn)*(y**ll)*(z**mm)*dexp(-expp(ipr)*rr)*coefpb(ipr,iact)
          ffd=dexp(-expp(ipr)*rr)*coefpb(ipr,iact)

          if (icode.le.3) then   
            xxx=-2.0d0*expp(ipr)*x**(nn+1)
            if(nn.ge.1) xxx=xxx+nn*x**(nn-1)
            fdx=fdx+xxx*(y**ll)*(z**mm)*ffd
            if (icode.eq.1) then   ! dx2
              xxx=-2.0d0*expp(ipr)*(nn+1)*x**nn+4.0d0*expp(ipr)**2.0d0*x**(nn+2)
              if(nn.gt.0) xxx=xxx-2.0d0*expp(ipr)*nn*x**nn
              if(nn.gt.1) xxx=xxx+nn*(nn-1)*x**(nn-2)
              fd2=fd2+xxx*(y**ll)*(z**mm)*ffd
            else if(icode.eq.2) then !dx dy
              yyy=-2.0d0*expp(ipr)*y**(ll+1)
              if(ll.ge.1) yyy=yyy+ll*y**(ll-1)
              fd2=fd2+xxx*yyy*(z**mm)*ffd
            else if(icode.eq.3) then  !dx dz
              yyy=-2.0d0*expp(ipr)*z**(mm+1)
              if(mm.ge.1) yyy=yyy+mm*z**(mm-1)
              fd2=fd2+xxx*yyy*(y**ll)*ffd
            end if
          end if

          if(mod(icode,2).eq.0) then
            xxx=-2.0d0*expp(ipr)*y**(ll+1)
            if(ll.ge.1) xxx=xxx+ll*y**(ll-1)
            fdy=fdy+xxx*(x**nn)*(z**mm)*ffd
            if (icode.eq.4) then   !dy2
              xxx=-2.0d0*expp(ipr)*(ll+1)*y**ll+4.0d0*expp(ipr)**2.0d0*y**(ll+2)
              if(ll.gt.0) xxx=xxx-2.0d0*expp(ipr)*ll*y**ll
              if(ll.gt.1) xxx=xxx+ll*(ll-1)*y**(ll-2)
              fd2=fd2+xxx*(x**nn)*(z**mm)*ffd
            else if(icode.eq.6) then !dy dz
              yyy=-2.0d0*expp(ipr)*z**(mm+1)
              if(mm.ge.1) yyy=yyy+mm*z**(mm-1)
              fd2=fd2+xxx*yyy*(x**nn)*ffd
            end if
          end if

          if(mod(icode,3).eq.0) then
            xxx=-2.0d0*expp(ipr)*z**(mm+1)
            if(mm.ge.1) xxx=xxx+mm*z**(mm-1)
            fdz=fdz+xxx*(y**ll)*(x**nn)*ffd
            if (icode.eq.9) then  !dz2 
              xxx=-2.0d0*expp(ipr)*(mm+1)*z**mm+4.0d0*expp(ipr)**2.0d0*z**(mm+2)
              if(mm.gt.0) xxx=xxx-2.0d0*expp(ipr)*mm*z**mm
              if(mm.gt.1) xxx=xxx+mm*(mm-1)*z**(mm-2)
              fd2=fd2+xxx*(x**nn)*(y**ll)*ffd
            end if
          end if
          k=k+1
        enddo

        ch(iact)=f
        chd(iact,1)=fd2
        select case(icode)
        case(1)
          chd(iact,2)=fdx
          chd(iact,3)=fdx
        case(2)
          chd(iact,2)=fdx
          chd(iact,3)=fdy
        case(3)
          chd(iact,2)=fdx
          chd(iact,3)=fdz
        case(4)
          chd(iact,2)=fdy
          chd(iact,3)=fdy
        case(6)
          chd(iact,2)=fdy
          chd(iact,3)=fdz
        case(9)
          chd(iact,2)=fdz
          chd(iact,3)=fdz
        end select
      end do
      x=0.0d0
      do mu=1,igr
        do nu=mu+1,igr
          x=x+2.0d0*p(mu,nu)*(chd(mu,1)*ch(nu)+chd(mu,2)*chd(nu,3)+chd(nu,1)*ch(mu)+chd(nu,2)*chd(mu,3))
        end do
        x=x+2.0d0*p(mu,mu)*(chd(mu,1)*ch(mu)+chd(mu,2)*chd(mu,3))
      end do
      deallocate (ch,chd)
      d2functxyz=x
      end
