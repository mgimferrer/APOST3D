      subroutine xc(npt,scr_a,scr,scr2)
!! scr_a = RHO, scr = SIGMA RHO, scr2 = EXC (MULTIPLICATION BY RHO REMAINING) !!
      use xc_f90_types_m
      use xc_f90_lib_m
      implicit real*8(a-h,o-z)
      TYPE(xc_f90_pointer_t) :: xc_func
      TYPE(xc_f90_pointer_t) :: xc_info
      include 'parameter.h'
      integer ideriv,npt
      real*8 scr_a(npt),scr(npt),scr2(npt)
      common /exchg/exch(maxat,maxat),xmix
      common /iops/iopt(100)
! (TO DO) xkdens = kinetic energy density (MG)
! (TO DO) Laplacian!!!
      allocatable :: scr2c(:)

      ideriv=0
      id_xcfunc = iopt(60)
      id_cfunc  = iopt(62)
      id_xfunc  = iopt(61)

      ALLOCATE(scr2c(npt))
      do ii=1,npt
        scr2c(ii)=ZERO
      end do 

!! EXC CALCULATION FROM LIBXC LIBRARIES !!

      if(id_xcfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_xcfunc,XC_UNPOLARIZED)
        select case (xc_f90_info_family(xc_info))
          case(XC_FAMILY_LDA)
            call xc_f90_lda_exc(xc_func,npt,scr_a(1),scr2(1))
          case(XC_FAMILY_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2(1))
          case(XC_FAMILY_HYB_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2(1))
          case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2(1))
          case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! CORRELATION SEPARATED !!

      if(id_cfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_cfunc,XC_UNPOLARIZED)
        select case (xc_f90_info_family(xc_info))
          case(XC_FAMILY_LDA)
            call xc_f90_lda_exc(xc_func,npt,scr_a(1),scr2c(1))
          case(XC_FAMILY_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2c(1))
          case(XC_FAMILY_HYB_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2c(1))
          case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2c(1))
          case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2c(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! EXCHANGE SEPARATED !!

      if(id_xfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_xfunc,XC_UNPOLARIZED)
        select case (xc_f90_info_family(xc_info))
          case(XC_FAMILY_LDA)
            call xc_f90_lda_exc(xc_func,npt,scr_a(1),scr2(1))
          case(XC_FAMILY_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2(1))
          case(XC_FAMILY_HYB_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_a(1),scr(1),scr2(1))
          case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2(1))
          case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_a(1),scr(1),lapl(1),tau(1),scr2(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! MULTIPLYING BY RHO !! 
!! EXCHANGE AND CORRELATION TOGETHER !!

      do ii=1,npt
        scr2(ii)=(scr2c(ii)+scr2(ii))*scr_a(ii)
      end do 
      DEALLOCATE(scr2c)

      end

! *****

      subroutine numint_dft(ndim,itotps,wp,rho,omp,omp2,chp,eto,pcoord,sat)
      use xc_f90_types_m
      use xc_f90_lib_m
      use ao_matrices
      use integration_grid
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(100)
      common /achi/achi(maxat,maxat),ibcp
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/exchg/exch(maxat,maxat),xmix
      dimension eto(maxat,maxat),Excmp(maxat,maxat)
      dimension wp(itotps),chp(itotps,ndim)
      dimension omp(itotps),omp2(itotps,nat),rho(itotps),pcoord(itotps,3)
! (TO DO) xkdens = kinetic energy density (MG)
! (TO DO) Laplacian!!!
      dimension exch2(maxat,maxat)
      dimension sat(igr,igr,nat)
      character*80 line

      allocatable :: chp2(:,:),scr(:),scrx(:,:)
      allocatable :: scr_a(:),rho_a(:),scr2(:)
      allocatable :: sab(:,:,:),chpd(:,:,:)

      idofr    =  Iopt(40)
      ithrebod =  Iopt(44)
      itype    =  Iopt(55)
      iatps    =  nang*nrad

      ALLOCATE(chp2(itotps,nocc),scr(itotps))
      ALLOCATE(scr_a(itotps),rho_a(itotps),scr2(itotps))
      ALLOCATE(sab(nocc,nocc,nat),chpd(itotps,nocc,3))
      ALLOCATE(scrx(igr,nocc))

!! THRESHOLD FOR BODEN CALCULATION USING THE BOND ORDER BETWEEN A PAIR OF ATOMS !!
      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if

!! EXCHANGE AND CORRELATION FOR DFT !!
!! BOND ORDER DENSITY !!
      write(*,*) " "
      write(*,*) " ------------------------------ "
      write(*,*) "  RKS-DFT ENERGY DECOMPOSITION  "
      write(*,*) " ------------------------------ "
      write(*,*) " "
      write(*,*) " USING BOND ORDER DENSITY APPROACH "
      write(*,'(2x,a36,1x,f10.6)') "Threshold for atom pair calculation:",threbod
      write(*,*) " "

      do iatom=1,nat
       do nu=1,igr
        do ii=1,nocc
          xx=ZERO
          do mu=1,igr
            xx=xx+c(mu,ii)*sat(mu,nu,iatom)
          end do
          scrx(nu,ii)=xx
        end do
       end do

       do ii=1,nocc
        do jj=1,nocc
          xx=ZERO
          do nu=1,igr
            xx=xx+c(nu,jj)*scrx(nu,ii)        
          end do
          sab(ii,jj,iatom)=xx
          if(ii.ne.jj) sab(jj,ii,iatom)=xx
        end do
       end do
      end do

!! CHECKING !!
      do ii=1,nocc                 
        do jj=ii,nocc
          x=ZERO
          if(ii.eq.jj) x=-ONE
          do kk=1,nat
            x=x+sab(ii,jj,kk)
          end do
          if(abs(x).gt.1.0d-3) then
            write(*,*) ii,jj,xx
            stop " PROBLEM WITH MOs OVERLAPS "
          end if
        end do
      end do

!! chp2 : MOs, chpd : MO DERIVATIVES !!
      do k=1,itotps
       do j=1,nocc
         xx=ZERO
         do i=1,igr
           xx=xx+c(i,j)*chp(k,i)
         end do
         chp2(k,j)=xx
       end do
      end do

      if(itype.gt.1) call grdrho(pcoord,chpd)

!! CALCULATING XC MATRIX USING BODEN !!
      write(*,*) " --------------------------------------- "
      write(*,*) "  BOND ORDER DENSITY FOR ALL ATOM PAIRS  "
      write(*,*) " --------------------------------------- "
      write(*,*) " "
      write(*,*) " --------------------------- "
      write(*,*) "  Atom   Atom   BODEN value  "
      write(*,*) " --------------------------- "
      x0=ZERO
      do iatom=1,nat
        do jatom=iatom+1,nat

!! COMPUTING BODEN FOR ATOMIC PAIRS (CONTROLLED BY THREBOD)
          bx0=bo(iatom,jatom)
          if(bx0.ge.threbod) then
            x1=ZERO
            do icenter=1,nat
              do jfut=iatps*(icenter-1)+1,iatps*icenter
                x2=wp(jfut)*omp(jfut)
                wa=omp2(jfut,iatom)
                wb=omp2(jfut,jatom)
                ff2=ZERO
                do i=1,nocc
                  do j=1,nocc 
                    ff=sab(i,j,jatom)*wa+sab(j,i,iatom)*wb
                    ff2=ff2+ff*chp2(jfut,i)*chp2(jfut,j)
                  end do
                end do

!! DEFINING A-B BODEN AS A-B + B-A, HENCE FACTOR OF 2
                scr_a(jfut)=ff2
                if(iatom.ne.jatom) scr_a(jfut)=TWO*ff2
                x1=x1+x2*scr_a(jfut)
              end do
            end do
            write(*,'(4x,i3,4x,i3,4x,f10.7)') iatom,jatom,x1

!! GRADIENT OF BODEN CALCULATED FOR GGA FUNCTIONALS !!
!! scr_a : BODEN FOR A GIVEN ATOM PAIR, scr : SIGMA BODEN scr2 : XC FUNCTIONAL VALUE !!
! Laplacian of the BODEN required for the MGGA functionals. NOT IMPLEMENTED
            if(itype.gt.1) call grdboden(itotps,omp2,chp2,chpd,scr,iatom,jatom,sab)
            call xc(itotps,scr_a,scr,scr2) 

!! ALLPOINTS INTEGRATION !!
             x1=ZERO
            do icenter=1,nat
              do ifut=iatps*(icenter-1)+1,iatps*icenter
                x1=x1+wp(ifut)*scr2(ifut)*omp2(ifut,icenter)*omp(ifut)
              end do
            end do
            exch2(iatom,jatom)=x1
          end if
        end do
      end do
      write(*,*) " ---------------------------- "
      write(*,*) " "

!! PRINTING !!
      write(*,*) " --------------------------------------- "
      write(*,*) "  DIATOMIC PURE KS-DFT XC TERMS (BODEN)  "
      write(*,*) " --------------------------------------- "
      write(*,*) " "
      exchen=ZERO
      do ii=1,nat
       do jj=ii,nat
         if(ii.ne.jj) exch2(jj,ii)=exch2(ii,jj)
         exchen=exchen+exch2(ii,jj)
       end do
      end do
      call MPRINT2(exch2,nat,maxat)
      write(*,*) " "

      if(itype.gt.1) then
       ideriv=0
       call sigma(pcoord,chp2,scr)
      end if

!! CALCULATE "EXACT" (ONE-CENTER) XC ENERGY !! 
!! Laplacian for the MGGA functionals can be calculated using LAPLACIAN subroutine (qtaim.f). NOT IMPLEMENTED
      call xc(itotps,rho,scr,scr2) 
      xtot=ZERO
      do icenter=1,nat
        x=ZERO
!! ALLPOINTS INTEGRATION !!
        do ifut=1,itotps 
          x=x+wp(ifut)*scr2(ifut)*omp(ifut)*omp2(ifut,icenter)
        end do
        exch(icenter,icenter)=x
        xtot=xtot+x
      end do

!! PRINTING THE "PURE" ONCE-CENTER TERMS !!
      write(*,*) " ----------------------------------------- "
      write(*,*) "  PURE KS-DFT XC ONE-CENTER TERMS (EXACT)  "
      write(*,*) " ----------------------------------------- "
      write(*,*) " "
      call MPRINT2(exch,nat,maxat)
      write(*,'(2x,a35,x,f14.7)') "KS-DFT exchange-correlation energy:",xtot
      write(*,*) " "

!! REARRANGING ATOMIC XC COMPONENTS !!
      write(*,*) " REARRANGING ATOMIC COMPONENTS "
      write(*,*) " "
      do ii=1,nat
        x0=ZERO
        do jj=1,nat
          if(ii.ne.jj) then
            exch(ii,jj)=exch2(ii,jj)
            x0=x0+exch(ii,jj)
          end if
        end do
        exch(ii,ii)=exch(ii,ii)-x0/TWO
      end do

!! PRINTING !!
      write(*,*) " ---------------------------------------------------------- "
      write(*,*) "  FINAL PURE KS-DFT EXCHANGE-CORRELATION ENERGY COMPONENTS  "
      write(*,*) " ---------------------------------------------------------- "
      write(*,*) " "
      exchen=ZERO
      do ii=1,nat
       do jj=ii,nat
         exchen=exchen+exch(ii,jj)
       end do
      end do
      call MPRINT2(exch,nat,maxat)
      write(*,'(2x,a47,x,f14.7)') "Sum of pure KS-DFT exchange-correlation energy:",exchen
      write(*,*) " "
      if(xmix.gt.ZERO) then
        write(*,*) " WARNING: HF-exchange part missing "
        write(*,*) " "
      end if
      if (idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Exc Decomposition' 
        call group_by_frag_mat(1,line,exch)
      end if

!! ACCUMULATING INTO THE ETO MATRIX !!
      xtot=ZERO
      do ii=1,nat
        eto(ii,ii)=eto(ii,ii)+exch(ii,ii)
        xtot=xtot+eto(ii,ii)
        do jj=ii+1,nat
          eto(ii,jj)=eto(ii,jj)+exch(ii,jj)
          eto(jj,ii)=eto(ii,jj)
          xtot=xtot+eto(ii,jj)
        end do
      end do
      etot=xtot
      DEALLOCATE(chp2,scr,scr_a,rho_a,scr2,sab,chpd,scrx)
      end

!! ***** !!

      subroutine sigma(pcoord,chp2,scr)
      use basis_set
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      dimension chp2(natoms*nrad*nang,nocc),scr(natoms*nang*nrad)
      dimension pcoord(natoms*nang*nrad,3)
      allocatable:: chpd(:,:),chp(:,:) 

      ipoints=natoms*nrad*nang
      ALLOCATE(chpd(ipoints,igr),chp(ipoints,igr))

!! GENERATING GRID FOR 2nd DERIVATIVE OVER AOs !!
      do kk=1,ipoints    
        scr(kk)=ZERO
      end do

      do ixyz=1,3
         do irun=1,ipoints
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
              alpha=expp(ipr)
              if(ixyz.eq.1) then
                dx=-2.0d0*alpha*(x**(nn+1))
                if(nn.ge.1) dx=dx+nn*x**(nn-1)
                dx=dx*(y**ll)*(z**mm)*dexp(-expp(ipr)*(rr**2))
              else if(ixyz.eq.2) then
                dx=-2.0d0*alpha*(y**(ll+1))
                if(ll.ge.1) dx=dx+ll*y**(ll-1)
                dx=dx*(x**nn)*(z**mm)*dexp(-expp(ipr)*(rr**2))
              else
                dx=-2.0d0*alpha*(z**(mm+1))
                if(mm.ge.1) dx=dx+mm*z**(mm-1)
                dx=dx*(x**nn)*(y**ll)*dexp(-expp(ipr)*(rr**2))
              end if
              f=f+dx*coefpb(ipr,iact)
              k=k+1
             enddo
             chp(irun,iact)=f
           enddo
         enddo
        
!! TO MOs !!
         do k=1,ipoints
           do j=1,nocc
             xx=ZERO
             do i=1,igr
               xx=xx+c(i,j)*chp(k,i)
             end do
             chpd(k,j)=xx
           end do
         end do
         do k=1,ipoints
          do i=1,nocc
            do j=1,nocc
             scr(k)=scr(k)+chp2(k,i)*chp2(k,j)*chpd(k,i)*chpd(k,j)
            end do
          end do
         end do
      end do

      do kk=1,ipoints
       scr(kk)=16.0d0*scr(kk)
      end do
      deallocate(chpd,chp)
      end

!! ***** !!

      subroutine grdrho(pcoord,chpd)
      use basis_set
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chpd(natoms*nrad*nang,nocc,3),pcoord(natoms*nrad*nang,3)
      allocatable:: chp(:,:)
      
      itotps=natoms*nrad*nang
      ALLOCATE(chp(itotps,igr)) 

!! GENERATING GRID FOR 2nd DERIVATIVE OVER AOs !!
      do ixyz=1,3
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
             alpha=expp(ipr)
             if(ixyz.eq.1) then
               dx=-2.0d0*alpha*(x**(nn+1))
               if(nn.ge.1) dx=dx+nn*x**(nn-1)
               dx=dx*(y**ll)*(z**mm)*dexp(-expp(ipr)*(rr**2))
             else if(ixyz.eq.2) then
               dx=-2.0d0*alpha*(y**(ll+1))
               if(ll.ge.1) dx=dx+ll*y**(ll-1)
               dx=dx*(x**nn)*(z**mm)*dexp(-expp(ipr)*(rr**2))
             else
               dx=-2.0d0*alpha*(z**(mm+1))
               if(mm.ge.1) dx=dx+mm*z**(mm-1)
               dx=dx*(x**nn)*(y**ll)*dexp(-expp(ipr)*(rr**2))
             end if
             f=f+dx*coefpb(ipr,iact)
             k=k+1
            enddo
            chp(irun,iact)=f
          enddo
         enddo
         
!! TO MOs !!
         do k=1,itotps
          do j=1,nocc
            xx=ZERO
            do i=1,igr
              xx=xx+c(i,j)*chp(k,i)
            end do
            chpd(k,j,ixyz)=xx
          end do
         end do

      end do
      DEALLOCATE(chp) 
      end

!! ***** !!

      subroutine ugrdrho(pcoord,chpd,chpbd)
      use basis_set
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chpd(natoms*nrad*nang,nalf,3),pcoord(natoms*nrad*nang,3)
      dimension chpbd(natoms*nrad*nang,nb,3)
      allocatable:: chp(:,:)
      
      itotps=natoms*nrad*nang
      ALLOCATE(chp(itotps,igr)) 

!! GENERATING GRID FOR 2nd DERIVATIVE OVER AOs !!
      do ixyz=1,3
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
             alpha=expp(ipr)
             if(ixyz.eq.1) then
               dx=-2.0d0*alpha*(x**(nn+1))
               if(nn.ge.1) dx=dx+nn*x**(nn-1)
               dx=dx*(y**ll)*(z**mm)*dexp(-expp(ipr)*(rr**2))
             else if(ixyz.eq.2) then
               dx=-2.0d0*alpha*(y**(ll+1))
               if(ll.ge.1) dx=dx+ll*y**(ll-1)
               dx=dx*(x**nn)*(z**mm)*dexp(-expp(ipr)*(rr**2))
             else
               dx=-2.0d0*alpha*(z**(mm+1))
               if(mm.ge.1) dx=dx+mm*z**(mm-1)
               dx=dx*(x**nn)*(y**ll)*dexp(-expp(ipr)*(rr**2))
             end if
             f=f+dx*coefpb(ipr,iact)
             k=k+1
            enddo
            chp(irun,iact)=f
          enddo
         enddo
         
!! TO MOs !!
         do k=1,itotps
          do j=1,nalf
            xx=ZERO
            do i=1,igr
              xx=xx+c(i,j)*chp(k,i)
            end do
            chpd(k,j,ixyz)=xx
          end do
          do j=1,nb  
            xx=ZERO
            do i=1,igr
              xx=xx+cb(i,j)*chp(k,i)
            end do
            chpbd(k,j,ixyz)=xx
          end do
         end do

      end do
      DEALLOCATE(chp) 
      end

!! ***** !!

      subroutine grdboden(itotps,omp2,chp2,chpd,scr,iat,kat,sab)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chp2(itotps,nocc),scr(itotps)
      dimension chpd(itotps,nocc,3),omp2(itotps,nat)
      dimension sab(nocc,nocc,nat)

!! GENERATING GRID FOR BODEN DERIVATIVE !!
      do kk=1,itotps    
        scr(kk)=ZERO
      end do
      do ixyz=1,3
         do irun=1,itotps
            w1=omp2(irun,iat)
            w2=omp2(irun,kat)
            xx=ZERO
            do ii=1,nocc
              do jj=1,nocc
                gab=sab(ii,jj,iat)*w2+sab(ii,jj,kat)*w1
                xx=xx+gab*chpd(irun,jj,ixyz)*chp2(irun,ii)
              end do
            end do
            scr(irun)=scr(irun)+xx*xx
         end do
      end do

      xfact=FOUR
      if(iat.ne.kat) xfact=xfact*FOUR
      do kk=1,itotps    
        scr(kk)=xfact*scr(kk)
      end do
      end

!! ***** !!

      subroutine xc_uks(npt,scr_ab,scr,scr2)
      use xc_f90_types_m
      use xc_f90_lib_m
      implicit real*8(a-h,o-z)
      TYPE(xc_f90_pointer_t) :: xc_func
      TYPE(xc_f90_pointer_t) :: xc_info
      include 'parameter.h'
      integer ideriv,npt
      common /exchg/exch(maxat,maxat),xmix
      common /iops/iopt(100)
! (TO DO) xkdens = kinetic energy density (MG)
! (TO DO) Laplacian!!!
      dimension scr_ab(2,npt),scr(3,npt),scr2(npt)
      allocatable :: scr2c(:)

      ideriv=0
      id_xcfunc = iopt(60)
      id_cfunc  = iopt(62)
      id_xfunc  = iopt(61)

      ALLOCATE(scr2c(npt))
      do ii=1,npt
        scr2c(ii)=ZERO
      end do 

!! scr_ab = RHO (ORDER: ALPHA, BETA) !!
!! scr = SIGMA RHO (ORDER: UP-UP, UP-DOWN, DOWN-DOWN) !!
!! scr2 = EXC (MULTIPLICATION BY RHO REMAINING) !!
      if(id_xcfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_xcfunc,XC_POLARIZED)
        select case (xc_f90_info_family(xc_info))
          case(XC_FAMILY_LDA)
            call xc_f90_lda_exc(xc_func,npt,scr_ab(1,1),scr2(1))
          case(XC_FAMILY_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2(1))
          case(XC_FAMILY_HYB_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2(1))
          case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2(1))
          case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! CORRELATION SEPARATED !!
      if(id_cfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_cfunc,XC_POLARIZED)
        select case (xc_f90_info_family(xc_info))
          case(XC_FAMILY_LDA)
            call xc_f90_lda_exc(xc_func,npt,scr_ab(1,1),scr2c(1))
          case(XC_FAMILY_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2c(1))
          case(XC_FAMILY_HYB_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2c(1))
          case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2c(1))
          case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2c(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! EXCHANGE SEPARATED !!
      if(id_xfunc.ne.0) then
        call xc_f90_func_init(xc_func,xc_info,id_xfunc,XC_POLARIZED)
        select case (xc_f90_info_family(xc_info))
          case(XC_FAMILY_LDA)
            call xc_f90_lda_exc(xc_func,npt,scr_ab(1,1),scr2(1))
          case(XC_FAMILY_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2(1))
          case(XC_FAMILY_HYB_GGA)
            call xc_f90_gga_exc(xc_func,npt,scr_ab(1,1),scr(1,1),scr2(1))
          case(XC_FAMILY_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2(1))
          case(XC_FAMILY_HYB_MGGA)
!            call xc_f90_mgga_exc(xc_func,npt,scr_ab(1,1),scr(1),lapl(1),tau(1),scr2(1))
        end select
        call xc_f90_func_end(xc_func)
      end if 

!! MULTIPLYING BY RHO !! 
!! EXCHANGE AND CORRELATION TOGETHER !!
      do ii=1,npt
        scr2(ii)=(scr2c(ii)+scr2(ii))*(scr_ab(1,ii)+scr_ab(2,ii))
      end do
      DEALLOCATE(scr2c)
      end 

!! ***** !!

      subroutine numint_dft_uks(ndim,itotps,wp,omp,omp2,chp,pcoord,eto)
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /iops/iopt(100)
      common/actual/iact,jat,icenter
      common /achi/achi(maxat,maxat),ibcp
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/exchg/exch(maxat,maxat),xmix
      character*80 line
      dimension eto(maxat,maxat)
      dimension wp(itotps),chp(itotps,ndim)
      dimension omp(itotps),omp2(itotps,nat),pcoord(itotps,3)
      dimension exch2(maxat,maxat)
! (TO DO) xkdens = kinetic energy density (MG)
! (TO DO) Laplacian!!!
      allocatable :: chp2(:,:),chp3(:,:),rho(:,:),scrall(:,:),scr2(:)
      allocatable :: chpd(:,:,:),sab(:,:,:),sab2(:,:,:),scr_bod(:,:)
      allocatable :: chpbd(:,:,:)

      idofr    = Iopt(40)
      ithrebod = Iopt(44)
      itype    = Iopt(55)
      iatps    = nang*nrad

      if(ithrebod.lt.1) then
        threbod=ZERO
      else
        threbod=real(ithrebod)/10000.0d0
      end if

      ALLOCATE(chp2(itotps,nalf),rho(2,itotps),scr_bod(2,itotps))
      ALLOCATE(chpd(itotps,nalf,3),chpbd(itotps,nb,3),sab2(nb,nb,nat))
      ALLOCATE(scrall(3,itotps),chp3(itotps,nb),sab(nalf,nalf,nat))
      ALLOCATE(scr2(itotps))

!! BUILDING RHO AND ACCUMULATING IN THE SAME VECTOR. (ORDER: ALPHA, BETA) !!
      do kk=1,itotps
        xx0=ZERO
        xx0b=ZERO
        do jj=1,nalf
          xx=ZERO
          xxb=ZERO
          do ii=1,igr
            xx=xx+c(ii,jj)*chp(kk,ii)
            if(jj.le.nb) xxb=xxb+cb(ii,jj)*chp(kk,ii)
          end do
          chp2(kk,jj)=xx
          if(jj.le.nb) then
            chp3(kk,jj)=xxb
            xx0b=xx0b+xxb*xxb
          end if
          xx0=xx0+xx*xx
        end do
        rho(1,kk)=xx0
        rho(2,kk)=xx0b
      end do

!! EXCHANGE AND CORRELATION FOR DFT !!
      write(*,*) " "
      write(*,*) " ------------------------------ "
      write(*,*) "  UKS-DFT ENERGY DECOMPOSITION  "
      write(*,*) " ------------------------------ "
      write(*,*) " "
      write(*,*) " USING BOND ORDER DENSITY APPROACH "
      write(*,'(2x,a36,1x,f10.6)') "Threshold for atom pair calculation:",threbod
      write(*,*) " "

!! sab FOR ALPHA, sab2 FOR BETA !!
      do iatom=1,nat
        do ii=1,nalf
          do jj=ii,nalf
            xx=ZERO
            xxb=ZERO
            do ifut=1,itotps
              xx=xx+chp2(ifut,ii)*chp2(ifut,jj)*wp(ifut)*omp(ifut)*omp2(ifut,iatom)
              if(jj.le.nb.and.ii.le.nb) xxb=xxb+chp3(ifut,ii)*chp3(ifut,jj)*wp(ifut)*omp(ifut)*omp2(ifut,iatom)
            end do
            sab(ii,jj,iatom)=xx
            if(ii.ne.jj) sab(jj,ii,iatom)=xx
            if(jj.le.nb.and.ii.le.nb) then
              sab2(ii,jj,iatom)=xxb
              if(ii.ne.jj) sab2(jj,ii,iatom)=xxb
            end if 
          end do
        end do
      end do

!! CHECKING ALPHA AND BETA OVERLAPS !! 
      do ii=1,nalf
        do jj=ii,nalf
          xa=ZERO
          if(ii.eq.jj) xa=-ONE
          do kk=1,nat
            xa=xa+sab(ii,jj,kk)
          end do
          if(abs(xa).gt.1.0d-3) then
            write(*,*) ii,jj,xa
            stop " PROBLEM WITH ALPHA MOs OVERLAPS "
          end if
        end do
      end do
      do ii=1,nb
        do jj=ii,nb
          xb=ZERO
          if(ii.eq.jj) xb=-ONE
          do kk=1,nat
            xb=xb+sab2(ii,jj,kk)
          end do
          if(abs(xa).gt.1.0d-3) then
            write(*,*) ii,jj,xb
            stop " PROBLEM WITH ALPHA MOs OVERLAPS "
          end if
        end do
      end do

!! CALCULATING RHO GRADIENTS !!
      if(itype.ne.1) call ugrdrho(pcoord,chpd,chpbd)

!! CALCULATING APPROXIMATED XC MATRIX USING BODEN !!
      write(*,*) " --------------------------------------- "
      write(*,*) "  BOND ORDER DENSITY FOR ALL ATOM PAIRS  "
      write(*,*) " --------------------------------------- "
      write(*,*) " "
      write(*,*) " --------------------------- "
      write(*,*) "  Atom   Atom   BODEN value  "
      write(*,*) " --------------------------- "
      x0=ZERO
      do iatom=1,nat
        do jatom=iatom+1,nat

!! COMPUTING BODEN FOR ATOM PAIRS IF BODEN LARGE ENOUGH (CONTROLLED BY THREBOD)

          bx0=bo(iatom,jatom)
          if(bx0.ge.threbod) then
            xx=ZERO
            xxb=ZERO
            do icenter=1,nat
              do jfut=iatps*(icenter-1)+1,iatps*icenter
                x2=wp(jfut)*omp(jfut)
                wa=omp2(jfut,iatom)
                wb=omp2(jfut,jatom)
                ff2=ZERO
                ff2b=ZERO
                do ii=1,nalf
                  do jj=1,nalf
                    ff=sab(ii,jj,jatom)*wa+sab(jj,ii,iatom)*wb
                    ff2=ff2+ff*chp2(jfut,ii)*chp2(jfut,jj)
                    if(jj.le.nb.and.ii.le.nb) then 
                      ffb=sab2(ii,jj,jatom)*wa+sab2(jj,ii,iatom)*wb
                      ff2b=ff2b+ffb*chp3(jfut,ii)*chp3(jfut,jj)
                    end if
                  end do
                end do

!! DEFINING A-B BODEN AS A-B + B-A, HENCE FACTOR OF 2 !!
!! scr_bod = BODEN MATRIX (ORDER: ALPHA, BETA)

                scr_bod(1,jfut)=ff2/TWO
                scr_bod(2,jfut)=ff2b/TWO
                if(iatom.ne.jatom) then
                  scr_bod(1,jfut)=ff2
                  scr_bod(2,jfut)=ff2b
                end if
                xx=xx+x2*scr_bod(1,jfut)
                xxb=xxb+x2*scr_bod(2,jfut)
              end do
            end do
            write(*,'(4x,i3,4x,i3,4x,f10.7)') iatom,jatom,xx+xxb

!! BODEN GRADIENT FOR UNRESTRICTED GGA FUNCTIONALS !!
            if(itype.ne.1) call grdboden_uks(itotps,chp2,chp3,omp2,chpd,chpbd,scrall,iatom,jatom,sab,sab2)
            call xc_uks(itotps,scr_bod,scrall,scr2)
            x1=ZERO
            do icenter=1,nat
              do ifut=iatps*(icenter-1)+1,iatps*icenter
                x1=x1+wp(ifut)*scr2(ifut)*omp2(ifut,icenter)*omp(ifut)
              end do
            end do
            exch2(iatom,jatom)=x1
          end if
        end do
      end do
      write(*,*) " ---------------------------- "
      write(*,*) " "

!! PRINTING !!
      write(*,*) " --------------------------------------- "
      write(*,*) "  DIATOMIC PURE KS-DFT XC TERMS (BODEN)  "
      write(*,*) " --------------------------------------- "
      write(*,*) " "
      exchen=ZERO
      do ii=1,nat
        exchen=exchen+exch2(ii,ii)
        do jj=ii,nat
          if(ii.ne.jj) then
            exch2(jj,ii)=exch2(ii,jj)
            exchen=exchen+exch2(ii,jj)
          end if 
        end do
      end do
      CALL MPRINT2(exch2,nat,maxat)
      write(*,*) " "

!! CALCULATE "EXACT" UNRESTRICTED XC ENERGY !!
      if(itype.ne.1) then
        ideriv=0
        call sigma_uks(pcoord,chp2,chp3,scrall)
      end if
      call xc_uks(itotps,rho,scrall,scr2)
      xtot=ZERO
      do icenter=1,nat
        x=ZERO
        do ifut=1,itotps
          x=x+wp(ifut)*scr2(ifut)*omp(ifut)*omp2(ifut,icenter)
        end do
        exch(icenter,icenter)=x
        xtot=xtot+x
      end do

!! PRINTING !!
      write(*,*) " ----------------------------------------- "
      write(*,*) "  PURE KS-DFT XC ONE-CENTER TERMS (EXACT)  "
      write(*,*) " ----------------------------------------- "
      write(*,*) " "
      call MPRINT2(exch,nat,maxat)
      write(*,'(2x,a35,x,f14.7)') "KS-DFT exchange-correlation energy:",xtot
      write(*,*) " "

!! RECALCULATING ATOMIC XC COMPONENTS !!
      write(*,*) " REARRANGING ATOMIC COMPONENTS "
      write(*,*) " "
      do ii=1,nat
        x0=ZERO
        do jj=1,nat
          if(ii.ne.jj) then
            exch(ii,jj)=exch2(ii,jj)
            x0=x0+exch(ii,jj)
          end if
        end do
        exch(ii,ii)=exch(ii,ii)-x0/TWO
      end do

!! FINAL PRINTING !!
      write(*,*) " ---------------------------------------------------------- "
      write(*,*) "  FINAL PURE KS-DFT EXCHANGE-CORRELATION ENERGY COMPONENTS  "
      write(*,*) " ---------------------------------------------------------- "
      write(*,*) " "
      exchen=ZERO
      do ii=1,nat
        do jj=ii,nat
          exchen=exchen+exch(ii,jj)
        end do
      end do
      call MPRINT2(exch,nat,maxat)
      write(*,'(2x,a47,x,f14.7)') "Sum of pure KS-DFT exchange-correlation energy:",exchen
      write(*,*) " "
      if(xmix.gt.ZERO) then
        write(*,*) " WARNING: HF-exchange part missing "
        write(*,*) " "
      end if
      if (idofr.eq.1) then
        line='   FRAGMENT ANALYSIS: Exc Decomposition' 
        call group_by_frag_mat(1,line,exch)
      end if

!! ACCUMULATING IN ETO !!
      xtot=ZERO
      do ii=1,nat
        eto(ii,ii)=eto(ii,ii)+exch(ii,ii)
        xtot=xtot+eto(ii,ii)
        do jj=ii+1,nat
          eto(ii,jj)=eto(ii,jj)+exch(ii,jj)
          eto(jj,ii)=eto(ii,jj)
          xtot=xtot+eto(ii,jj)
        end do
      end do
      DEALLOCATE(chp2,rho,scr_bod,chpd,chpbd,sab2,scrall,chp3,sab,scr2)
      end 

!! ***** !!

      subroutine grdboden_uks(itotps,chp2,chp3,omp2,chpd,chpbd,scrall,iat,kat,sab,sab2)
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension chp2(itotps,nalf),chp3(itotps,nb),scrall(3,itotps)
      dimension chpd(itotps,nalf,3),chpbd(itotps,nb,3)
      dimension sab(nalf,nalf,nat),sab2(nb,nb,nat),omp2(itotps,nat)

!! GENERATING SIGMA BOND ORDER DENSITY VECTOR !!
      do kk=1,itotps
        scrall(1,kk)=ZERO
        scrall(2,kk)=ZERO
        scrall(3,kk)=ZERO
      end do

!! CALCULATING AND ACCUMULATING SIGMA BODEN (ORDER: UP-UP, UP-DOWN, DOWN-DOWN) !!
      do ixyz=1,3
        do irun=1,itotps
          w1=omp2(irun,iat)
          w2=omp2(irun,kat)
          xxa=ZERO
          xxb=ZERO
          do ii=1,nalf
            do jj=1,nalf
              gab=sab(ii,jj,iat)*w2+sab(ii,jj,kat)*w1
              xxa=xxa+gab*chpd(irun,jj,ixyz)*chp2(irun,ii)
              if(jj.le.nb.and.ii.le.nb) then
               gab2=sab2(ii,jj,iat)*w2+sab2(ii,jj,kat)*w1
               xxb=xxb+gab2*chpbd(irun,jj,ixyz)*chp3(irun,ii)
              end if
            end do
          end do
          scrall(1,irun)=scrall(1,irun)+xxa*xxa
          scrall(2,irun)=scrall(2,irun)+xxa*xxb
          scrall(3,irun)=scrall(3,irun)+xxb*xxb
        end do
      end do
      if(iat.ne.kat) then
        do kk=1,itotps    
          scrall(1,kk)=FOUR*scrall(1,kk)
          scrall(2,kk)=FOUR*scrall(2,kk)
          scrall(3,kk)=FOUR*scrall(3,kk)
        end do
      end if
      end

!! ***** !!

      subroutine func_info_print(id_func,itype)

      use xc_f90_types_m
      use xc_f90_lib_m

      implicit real*8(a-h,o-z)

      TYPE(xc_f90_pointer_t) :: xc_func
      TYPE(xc_f90_pointer_t) :: xc_info

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common/exchg/exch(maxat,maxat),xmix
      common /iops/iopt(100)

      character*120 name_ref        
      character*80 name_func

      if(kop.ne.1) call xc_f90_func_init(xc_func,xc_info,id_func,XC_UNPOLARIZED)
      if(kop.eq.1) call xc_f90_func_init(xc_func,xc_info,id_func,XC_POLARIZED)
      call xc_f90_hyb_exx_coef(xc_func,xmix)
      call xc_f90_info_name(xc_info,name_func)

      write(*,*) " -------------------------------- "
      write(*,*) "  DENSITY FUNCTIONAL INFORMATION  "
      write(*,*) " -------------------------------- "
      write(*,*) " "
      write(*,*) " Functional name -->",trim(name_func)
      ii=0
      call xc_f90_info_refs(xc_info,ii,name_ref)
      do while(ii.ge.0)
        write(*,'(2x,a15,i1,a1,x,a120)') "Reference --> [",ii,"]",name_ref
        call xc_f90_info_refs(xc_info,ii,name_ref)
      end do
      select case(xc_f90_info_kind(xc_info))
        case(XC_EXCHANGE)
          write(*,*) " Functional type --> exchange"
        case(XC_CORRELATION)
          write(*,*) " Functional type --> correlation"
        case(XC_EXCHANGE_CORRELATION)
          write(*,*) " Functional type --> exchange-correlation"
        case(XC_KINETIC)
          write(*,*) " Functional type --> kinetic energy"
        case default
          write(*,*) " Functional type --> unknown"
      end select

!! itype = 1 LDA, 2 GGA, 3 META-GGA. FOR HYB XMIX > 0. itype SAVED IN IOPT(55) !!

      itype=0
      select case(xc_f90_info_family(xc_info))
        case(XC_FAMILY_UNKNOWN)
          write(*,*) " Family of the functional unknown"
          stop
        case(XC_FAMILY_NONE)
          write(*,*) " Non-Family identified functional"
          stop
        case(XC_FAMILY_LDA)
          write(*,*) " LDA functional selected"
          itype=1
        case(XC_FAMILY_GGA)
          write(*,*) " GGA functional selected"
          itype=2
        case(XC_FAMILY_HYB_GGA)
          write(*,*) " HYBRID-GGA functional selected"
          write(*,'(2x,a26,x,f5.3)') "HF-type exchange coeff -->",xmix
          itype=2
        case(XC_FAMILY_MGGA)
          write(*,*) " META-GGA functional selected"
          itype=3
          write(*,*) " META-GGA still in development!!!" !! MG: TO DO LIST !!
          stop
        case(XC_FAMILY_HYB_MGGA)
          write(*,*) " HYBRID-META-GGA functional selected"
          write(*,'(2x,a26,x,f5.3)') "HF-type exchange coeff -->",xmix
          itype=3
          write(*,*) " META-GGA still in development!!!" !! MG: TO DO LIST !!
          stop
      end select
      call xc_f90_func_end(xc_func)

      end 

!! ***** !!

      subroutine sigma_uks(pcoord,chp2,chp3,scr)
      use basis_set
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      dimension chp2(natoms*nrad*nang,nalf),scr(3,natoms*nang*nrad)
      dimension chp3(natoms*nrad*nang,nb)
      dimension pcoord(natoms*nang*nrad,3)
      allocatable:: chpd(:,:),chp(:,:),chpbd(:,:)


      ipoints=natoms*nrad*nang
      ALLOCATE(chpd(ipoints,igr),chp(ipoints,igr),chpbd(ipoints,igr))

!! GENERATING GRID FOR 2nd DERIVATIVE OVER AOs !!
      do kk=1,ipoints
       do jj=1,3
        scr(jj,kk)=ZERO
       end do
      end do

      do ixyz=1,3
         do irun=1,ipoints
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
             alpha=expp(ipr)
             if(ixyz.eq.1) then
               dx=-2.0d0*alpha*(x**(nn+1))
               if(nn.ge.1) dx=dx+nn*x**(nn-1)
               dx=dx*(y**ll)*(z**mm)*dexp(-expp(ipr)*(rr**2))
             else if(ixyz.eq.2) then
               dx=-2.0d0*alpha*(y**(ll+1))
               if(ll.ge.1) dx=dx+ll*y**(ll-1)
               dx=dx*(x**nn)*(z**mm)*dexp(-expp(ipr)*(rr**2))
             else
               dx=-2.0d0*alpha*(z**(mm+1))
               if(mm.ge.1) dx=dx+mm*z**(mm-1)
               dx=dx*(x**nn)*(y**ll)*dexp(-expp(ipr)*(rr**2))
             end if
             f=f+dx*coefpb(ipr,iact)
             k=k+1
            enddo
            chp(irun,iact)=f
          enddo
         enddo
         
!! TO MOs !!  
         do k=1,ipoints
           do j=1,nalf
             xx=ZERO
             xxb=ZERO
             do i=1,igr
               xx=xx+c(i,j)*chp(k,i)
               if(j.le.nb) xxb=xxb+cb(i,j)*chp(k,i)
             end do
             chpd(k,j)=xx
             if(j.le.nb) chpbd(k,j)=xxb
           end do
         end do

!! CALCULATION OF SIGMA FOR COUPLING WITH LIBXC LIBRARIES (ORDER: UP-UP, UP-DOWN, DOWN-DOWN) !!
         do k=1,ipoints
           do i=1,nalf
             do j=1,nalf
               scr(1,k)=scr(1,k)+chp2(k,i)*chp2(k,j)*chpd(k,i)*chpd(k,j)
               if(j.le.nb) scr(2,k)=scr(2,k)+chp2(k,i)*chp3(k,j)*chpd(k,i)*chpbd(k,j)
               if(i.le.nb.and.j.le.nb) scr(3,k)=scr(3,k)+chp3(k,i)*chp3(k,j)*chpbd(k,i)*chpbd(k,j)
             end do
           end do
         end do
      end do
c
      do k=1,ipoints
        do j=1,3
        scr(j,k)=4.0d0*scr(j,k)
       end do 
      end do 
      DEALLOCATE(chpd,chpbd,chp)
      end 

!! ***** !!
