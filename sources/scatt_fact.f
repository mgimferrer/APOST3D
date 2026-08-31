!! *********************************************************************** !!
!! X-RAY SCATTERING FACTORS (# METHOD/SCATT-FACT) -- one subroutine.       !!
!!   scattering_factors -- per-atom real-space density, spherically        !!
!!                         averaged, Fourier-transformed into XRSF(s)      !!
!!                         values, written to <name>.xrsf                  !!
!! *********************************************************************** !!

!! ***** !!

!! ********************************************************************* !!
!! subroutine: scattering_factors                                       !!
!! purpose: for each atom, extracts its real-space density (rho_at),     !!
!!   spherically averages it over the angular grid (rho_sph), then       !!
!!   Fourier-transforms it at a fixed table of scattering angles (sval)  !!
!!   into X-ray scattering factors (fxA), written to <name>.xrsf.        !!
!! arguments (all read-only):                                            !!
!!   itotps (in) -- total number of grid points (nat*iatps)              !!
!!   wp     (in) -- integration weight of each grid point                !!
!!   rho    (in) -- electron density at each grid point                  !!
!!   omp2   (in) -- becke/tfvc (or hirshfeld) weight of each point for    !!
!!                  every atom                                           !!
!!   pcoord (in) -- xyz coordinates of each grid point (unused today)     !!
!! author:                                                                !!
!! ********************************************************************* !!
      subroutine scattering_factors(itotps,wp,rho,omp2,pcoord)

      use ao_matrices
      use integration_grid

      implicit double precision(a-h,o-z)

      include 'parameter.h'

      common /nat/nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/coord(3,maxat),zn(maxat),iznuc(maxat)
      common /filename/name0

      character*60 name0
      character*80 nameout

      integer, intent(in) :: itotps
      integer :: isval

      dimension :: wp(itotps),pcoord(itotps,3)
      dimension :: omp2(itotps,nat),rho(itotps)

      allocatable :: rho_at(:),rho_sph(:),fxA(:)

!! fixed table of scattering angles (inverse Angstrom) fxA is evaluated at !!
      double precision :: sval(56)
      data sval /0.00d0, 0.01d0, 0.02d0, 0.03d0, 0.04d0, 0.05d0, 0.06d0, 0.07d0, 0.08d0, 0.09d0,
     +           0.10d0, 0.11d0, 0.12d0, 0.13d0, 0.14d0, 0.15d0, 0.16d0, 0.17d0, 0.18d0, 0.19d0,
     +           0.20d0, 0.22d0, 0.24d0, 0.25d0, 0.26d0, 0.28d0, 0.30d0, 0.32d0, 0.34d0, 0.35d0,
     +           0.36d0, 0.38d0, 0.40d0, 0.42d0, 0.44d0, 0.45d0, 0.46d0, 0.48d0, 0.50d0, 0.55d0,
     +           0.60d0, 0.65d0, 0.70d0, 0.80d0, 0.90d0, 1.00d0, 1.10d0, 1.20d0, 1.30d0, 1.40d0,
     +           1.50d0, 1.60d0, 1.70d0, 1.80d0, 1.90d0, 2.00d0/

      iatps = nrad*nang
      isval = 56
      write(*,'(2x,a,1x,i0,1x,a,1x,i0)') 'Radial points:',nrad,'Angular points:',nang
      write(*,*) " "

      nameout=trim(name0)//".xrsf"
      open(69,file=nameout,status="unknown")
      write(*,'(2x,a,1x,a)') 'Printing XRSF information in file',trim(nameout)
      write(*,*) " "

      ALLOCATE(rho_at(iatps))
      ALLOCATE(rho_sph(nrad))
      ALLOCATE(fxA(isval))

!! parallelization: not done. The per-atom work (O(nrad*nang) spherize + !!
!! O(isval*nrad) transform) is a plausible future OMP target if this     !!
!! ever shows up as hot, but SCATT-FACT has zero test coverage and is    !!
!! deprioritized, so not pursued now.                                    !!
      do icenter=1,nat

        rho_at=ZERO
        rho_sph=ZERO
        fxA=ZERO

!! rho^A in real space (non-spherical), from the fuzzy-atom-weighted rho !!
        xrhoA=ZERO
        iifut=1
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          rho_at(iifut)=omp2(ifut,icenter)*rho(ifut)
          xrhoA=xrhoA+wp(ifut)*rho_at(iifut)
          iifut=iifut+1
        end do

        if(iifut-1.ne.iatps) stop 'scattering_factors: vector dimension mismatch (iifut =/ iatps)'

!! spherize rho^A: average over the angular grid at each radial shell    !!
!! (angular weights sum to 1), then integrate the spherized density as a !!
!! cross-check against xrhoA above.                                      !!
        xrhoA2=ZERO
        iifut=1
        do irad=1,nrad
          xxav=ZERO
          do iang=1,nang
            xxav=xxav+w(iang)*rho_at(iifut)
            iifut=iifut+1
          end do
          rho_sph(irad)=xxav
          xrhoA2=xrhoA2+wr(irad)*xr(irad)*xr(irad)*rho_sph(irad)
        end do
        xrhoA2=xrhoA2*FOUR*pi

        write(*,'(2x,a,1x,i0)') 'For atom:',icenter
        write(*,'(2x,a,1x,f11.6)') 'Integrated rho^(A)     =',xrhoA
        write(*,'(2x,a,1x,f11.6)') 'Integrated rho^(A,sph) =',xrhoA2
        write(*,*) " "

!! Fourier-transform the spherized density at each scattering angle sval !!
!! (converted from inverse Angstrom to inverse Bohr); the s=0 limit is   !!
!! just the (spherized) integrated density, sin(x)/x -> 1.               !!
        do iisval=1,isval
          xfxA=ZERO
          xx=FOUR*pi*sval(iisval)
          xx=xx*angtoau
          do irad=1,nrad
            if(sval(iisval).lt.10d-8) then
              xfxA=xfxA+wr(irad)*xr(irad)*xr(irad)*rho_sph(irad)
            else
              xft=DSIN(xx*xr(irad))/(xx*xr(irad))
              xfxA=xfxA+wr(irad)*xr(irad)*xr(irad)*rho_sph(irad)*xft
            end if
          end do
          fxA(iisval)=xfxA*FOUR*pi
        end do

        write(69,'(a,1x,i0)') 'Scattering angle (s) and associated XRSF values for atom',icenter
        write(69,'(a,1x,f11.6)') 'Atomic (electron) population        =',xrhoA
        write(69,'(a,1x,f11.6)') 'Integrated spherized atomic density =',xrhoA2
        write(69,'(a,1x,i0)') 'Atom number                         =',INT(zn(icenter))
        write(69,'(a,1x,f11.6)') 'Atomic (AIM) charge                 =',zn(icenter)-xrhoA
        write(69,*) " "
        write(69,*) " ------------------------- "
        write(69,*) "       s        fxA(s)     "
        write(69,*) " ------------------------- "
        do iisval=1,isval
          write(69,'(2x,2f11.6)') sval(iisval),fxA(iisval)
        end do
        write(69,*) " ------------------------- "
        write(69,*) " "
      end do

      DEALLOCATE(rho_at,rho_sph,fxA)

      close(69)
      end
