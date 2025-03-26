!! ***** !!
!! SET OF SUBROUTINES FOR THE EVALUATION OF X-RAY SCATTERING FACTORS !!
!! ***** !!

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

!! ENTERING s VALUES AS DATA HERE !!
      double precision :: sval(56)
      data sval /0.00d0, 0.01d0, 0.02d0, 0.03d0, 0.04d0, 0.05d0, 0.06d0, 0.07d0, 0.08d0, 0.09d0,
     +           0.10d0, 0.11d0, 0.12d0, 0.13d0, 0.14d0, 0.15d0, 0.16d0, 0.17d0, 0.18d0, 0.19d0,
     +           0.20d0, 0.22d0, 0.24d0, 0.25d0, 0.26d0, 0.28d0, 0.30d0, 0.32d0, 0.34d0, 0.35d0,
     +           0.36d0, 0.38d0, 0.40d0, 0.42d0, 0.44d0, 0.45d0, 0.46d0, 0.48d0, 0.50d0, 0.55d0,
     +           0.60d0, 0.65d0, 0.70d0, 0.80d0, 0.90d0, 1.00d0, 1.10d0, 1.20d0, 1.30d0, 1.40d0,
     +           1.50d0, 1.60d0, 1.70d0, 1.80d0, 1.90d0, 2.00d0/

      iatps = nrad*nang
      isval = 56
      write(*,*) "INFO: nrad, nang: ",nrad,nang
      write(*,*) " "

!! PREPARING OUTPUT FILE !!
      nameout=trim(name0)//".xrsf"
      open(69,file=nameout,status="unknown")
      write(*,*) " Printing XRSF information in file ",trim(nameout)
      write(*,*) " "

!! ALLOCATING RHO VECTORS !!
      ALLOCATE(rho_at(iatps)) !! DIMENSIONATED TO THE SINGLE ATOM GRID SIZE !!
      ALLOCATE(rho_sph(nrad)) !! DIMENSIONATED TO NUMBER OF RADIAL POINTS !!
      ALLOCATE(fxA(isval))

!! DOING FOR EACH ATOM !!
      do icenter=1,nat

!! ZEROING MATRICES !!
        rho_at=ZERO
        rho_sph=ZERO
        fxA=ZERO

!! OBTAINING RHO^A IN REAL-SPACE (NON-SPHERICAL) FROM RHO !!
        xrhoA=ZERO
        iifut=1
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          rho_at(iifut)=omp2(ifut,icenter)*rho(ifut)
          xrhoA=xrhoA+wp(ifut)*rho_at(iifut)
          iifut=iifut+1
        end do

!! CHECK: iifut = iatps !!
        if(iifut-1.ne.iatps) then
          write(*,*) " Vector dimension problem, iifut =/ iatps "
          stop
        end if

!! SPHERIZING RHO^A !!
        xrhoA2=ZERO
        iifut=1
        do irad=1,nrad
          xxav=ZERO
          do iang=1,nang
            xxav=xxav+w(iang)*rho_at(iifut)
            iifut=iifut+1
          end do

!! SAVING AVERAGE DENSITY VALUE IN rho_sph !!
!! REMIND THAT SUM OF ANGULAR WEIGHTS IS 1 !!
          rho_sph(irad)=xxav

!! CHECK: INTEGRATING NOW THE SPHERICAL AVERAGED RHO !!
          xrhoA2=xrhoA2+wr(irad)*xr(irad)*xr(irad)*rho_sph(irad)
        end do

!! CHECK: FINAL VALUE OF xx2 !!
        xrhoA2=xrhoA2*FOUR*pi

!! CHECK: PRINTING BOTH QUANTITIES INTEGRATED !!
        write(*,'(a11,i3)') " For atom: ",icenter
        write(*,'(a26,f11.6)') " Integrated rho^(A)     = ",xrhoA
        write(*,'(a26,f11.6)') " Integrated rho^(A,sph) = ",xrhoA2
        write(*,*) " "

!! FOURIER TRANSFORMATION TIME: USING rho_sph TO EVALUATE fx_A !!
        do iisval=1,isval
          xfxA=ZERO
          xfxA2=ZERO
          xx=FOUR*pi*sval(iisval)
          xx=xx*angtoau !! TRANSFORMING FROM INVERSE ANGSTROM TO INVERSE BOHR !! 
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

!! PRINTING THE fxA VALUES !!
        write(69,'(a58,i3)') " Scattering angle (s) and associated XRSF values for Atom ",icenter
        write(69,'(a38,f11.6)') " Atomic (electron) population        = ",xrhoA
        write(69,'(a38,f11.6)') " Integrated spherized atomic density = ",xrhoA2
        write(69,'(a38,1x,i3)') " Atom number                         = ",INT(zn(icenter))
        write(69,'(a38,f11.6)') " Atomic (AIM) charge                 = ",zn(icenter)-xrhoA
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

!! DEALLOCATING MATRICES !!
      DEALLOCATE(rho_at,rho_sph,fxA)

      close(69)
      end

!! ***** !!


