!! ********************************************************************* !!
!! subroutine: spincorr                                                  !!
!! purpose: local spin / delocalization-index decomposition from the     !!
!!   post-HF 1- and 2-RDMs (dm1/dm2, natural-orbital basis). Prints the  !!
!!   bond order matrix, Mayer's "new" bond order (CPL 554, 83 (2012)),   !!
!!   the LI/DI matrix, effectively unpaired electrons (u_A), and the     !!
!!   S^2 decomposition (a=3/4, Ramos-Cordoba/Salvador local spin) --     !!
!!   plus their DOFRAGS fragment-analysis breakdowns when requested.     !!
!! arguments:                                                            !!
!!   sat (in) -- per-atom AO overlap matrix (igr,igr,nat)                !!
!!   dm1 (in) -- one-electron reduced density matrix, spin-orbital       !!
!!               basis (nspinorb,nspinorb)                               !!
!!   dm2 (in) -- spinless two-electron reduced density matrix, natural-  !!
!!               orbital basis (norb,norb,norb,norb)                    !!
!! parallelization: not done. The O(norb^4*nat^2) i/j/k/l/iat/jat loop    !!
!!   (the dominant cost) accumulates into shared S234/DI/BO(iat,jat) on  !!
!!   every i/j/k/l iteration -- a genuine reduction, not a plain private- !!
!!   write split, and this codebase has no existing precedent for an     !!
!!   OMP array-REDUCTION on fixed-size arrays this large. The cheaper    !!
!!   O(norb^2*nat) satmo-transform loops reuse a single shared scr(:,:)  !!
!!   scratch array across iat iterations (not iat-indexed), so paral-   !!
!!   lelizing over iat would race on scr without first giving each      !!
!!   thread its own copy. Both left serial pending a dedicated pass.     !!
!! author:                                                               !!
!! ********************************************************************* !!
      subroutine spincorr(sat,dm1,dm2)
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      parameter (TOL=1.0d-10)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /localspin/xlsa(maxat,maxat),ua(maxat)
      common /iops/iopt(200)
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      dimension sat(igr,igr,nat)
      dimension dm1(nspinorb,nspinorb)
      dimension dm2(norb,norb,norb,norb)
      character*80 line

      dimension S234(maxat,maxat)
      dimension bonew(maxat,maxat)
      allocatable scr(:,:),satmo(:,:,:)
      allocatable sfdm1(:,:), psdm1(:,:),u_no(:)

      idofr= iopt(40)

      allocate (scr(igr,igr))
      allocate (satmo(igr,igr,nat))
      allocate (sfdm1(norb,norb),psdm1(norb,norb),u_no(norb))

!! effective unpaired electrons in the NO basis, same dimension as MOs !!
      do iat=1,nat
        do i=1,norb
          do nu=1,igr
            x=0.0d0
            do mu=1,igr
              x=x+c_no(mu,i)*sat(mu,nu,iat)
            end do
            scr(i,nu)=x
          end do
        end do
        do i=1,norb
          do j=1,norb
            x=0.0d0
            do nu=1,igr
              x=x+c_no(nu,j)*scr(i,nu)
            end do
            satmo(i,j,iat)=x
          end do
        end do
      end do

      xnd=0.0d0
      do i=1,norb
        u_no(i)=2.0d0*occ_no(i,i)-occ_no(i,i)**2.0d0
        xnd=xnd+u_no(i)
      end do

!! atom condensed !!
      do iat=1,nat
        xx=0.0d0
        do i=1,norb
          xx=xx+u_no(i)*satmo(i,i,iat)
        end do
        ua(iat)=xx
      end do

!! satmo(i,j,jat) is read with its indices exchanged (satmo(j,i,jat))    !!
!! throughout this routine relative to the naive iat/jat-symmetric form  !!
!! -- a deliberate correction, kept consistently below.                 !!
      do iat=1,nat
        do jat=iat,nat
          yy=0.0d0
          do i=1,norb
            do j=1,norb
              yy=yy+occ_no(i,i)*occ_no(j,j)*satmo(i,j,iat)*satmo(j,i,jat)
            end do
          end do
          xx=0.0d0
          do i=1,norb
            if(u_no(i).gt.TOL) then
              uno12=dsqrt(u_no(i))
              do j=1,norb
                if(u_no(j).gt.TOL) xx=xx+uno12*dsqrt(u_no(j))*satmo(i,j,iat)*satmo(j,i,jat)
              end do
            end if
          end do
          if(iat.eq.jat) then
            bonew(iat,jat)= (yy+xx)/2.0d0
          else
            bonew(iat,jat)= yy+xx
            bonew(jat,iat)= bonew(iat,jat)
          end if
        end do
      end do

!! LSA: transforming now to the MO basis !!
      do iat=1,nat
        do i=1,norb
          do nu=1,igr
            x=0.0d0
            do mu=1,igr
              x=x+c(mu,i)*sat(mu,nu,iat)
            end do
            scr(i,nu)=x
          end do
        end do
        do i=1,norb
          do j=1,norb
            x=0.0d0
            do nu=1,igr
              x=x+c(nu,j)*scr(i,nu)
            end do
            satmo(i,j,iat)=x
          end do
        end do

      end do

!! orthogonality check !!
      do i=1,norb
        do j=i,norb
          delta=0.0d0
          if(i.eq.j) delta=1.0d0
          do iat=1,nat
            delta=delta-satmo(i,j,iat)
          end do
          if(delta.gt.1.0d-2) then
            write(*,*) 'Large deviation from orthonormality in MOs:',i,j,delta
          end if
        end do
      end do

!! initializing arrays for local spin, DI and U decomposition; local    !!
!! spin formula (-1+2a)*Gamma_ijij - 0.5*Gamma_ijji, a=0 (Alcoba) vs     !!
!! a=3/4 (Ramos-Cordoba) -- this routine implements a=3/4.               !!
      do i=1,nat
        do j=1,nat
          s234(i,j)=0.0d0
          DI(i,j)=0.0d0
          BO(i,j)=0.0d0
        end do
      end do

!! P and Ps in the MO basis !!
      do i=1,norb
        do k=1,norb
          sfdm1(i,k)=dm1((i-1)*2+1,(k-1)*2+1)+dm1((i-1)*2+2,(k-1)*2+2)
          psdm1(i,k)=dm1((i-1)*2+1,(k-1)*2+1)-dm1((i-1)*2+2,(k-1)*2+2)
        end do
      end do

!! contributions from the spinless cumulant of dm2, 1122 form (spin      !!
!! density contributions included in the cumulant, i.e. not removed).    !!
!! satmo indices are exchanged throughout, same correction as above.     !!
      do i=1,norb
        do j=1,norb
          do k=1,norb
            do l=1,norb
              xx1=dm2(i,j,k,l)-sfdm1(i,j)*sfdm1(k,l)+(sfdm1(l,i)*sfdm1(j,k)+psdm1(l,i)*psdm1(j,k))/2.0d0
              xx0=dm2(i,j,k,l)-sfdm1(i,j)*sfdm1(k,l)+(sfdm1(l,i)*sfdm1(j,k))/2.0d0
              xx2=(sfdm1(l,i)*sfdm1(j,k))/2.0d0+(psdm1(l,i)*psdm1(j,k))/2.0d0
              do iat=1,nat
                do jat=iat,nat
                  S234(iat,jat)=s234(iat,jat)+xx0*0.5d0*(satmo(j,i,iat)*satmo(l,k,jat)-satmo(l,i,iat)*satmo(j,k,jat))
                  DI(iat,jat)=DI(iat,jat)-xx1*2.0d0*(satmo(j,i,iat)*satmo(l,k,jat))
                  BO(iat,jat)=BO(iat,jat)+xx2*2.0d0*(satmo(j,i,iat)*satmo(l,k,jat))
                end do
              end do
            end do
          end do
        end do
      end do

!! contributions from dm1 -- simply the density of eff. unpaired elec., !!
!! already calculated above. !!
      do iat=1,nat
        S234(iat,iat)=s234(iat,iat)+0.75d0*ua(iat)
      end do
      deallocate (satmo,sfdm1,psdm1)

!! symmetrizing and computing LI/DI !!
      do i=1,nat
        do j=i+1,nat
          s234(j,i)=s234(i,j)
          bo(j,i)=bo(i,j)
          DI(i,j)=DI(i,j)+bo(i,j)
          DI(j,i)=DI(i,j)
        end do
        BO(i,i)=BO(i,i)/2.0d0
        DI(i,i)=DI(i,i)/2.0d0+bo(i,i)
      end do
!! saving Ramos-Cordoba decomposition !!
      do i=1,nat
        do j=1,nat
          xlsa(i,j)=s234(i,j)
        end do
      end do

!! sum checks !!
      x0=0.d0
      x1=0.d0
      x2=0.d0
      do iat=1,nat
        x0=x0+di(iat,iat)
        x1=x1+bonew(iat,iat)
        do jat=1,nat
          x2=x2+s234(iat,jat)
          if(iat.ne.jat) x0=x0+di(iat,jat)/2.0d0
          if(iat.ne.jat) x1=x1+bonew(iat,jat)/2.0d0
        enddo
      enddo

!! printing !!
      WRITE(*,8)
 8    FORMAT(1x,/21X,'  APOST3D BOND ORDER MATRIX')
      print *,' '
      call MPRINT(bo,nat,maxat)
      print *,' '

      WRITE(*,5)
 5    FORMAT(1x,/21X,'  APOST3D NEW BOND ORDER MATRIX')
      print *,' Improved definition by I. Mayer on CPL 554, 83 (2012)'
      print *,' '
      call MPRINT(bonew,nat,maxat)
      print *,' '
      write(*,'(a13,f10.5)') ' Sum check = ' ,x1
      print *,' '

      WRITE(*,4)
 4    FORMAT(1x,/21X,'    APOST3D LI/DI MATRIX'//)
      call MPRINT(DI,nat,maxat)
      print *,' '
      write(*,'(a13,f10.5)') ' Sum check = ' ,x0
      if (idofr.eq.1) then
        line ='   FRAGMENT ANALYSIS : Deloc. Index'
        call group_by_frag_mat(1,line ,di)
      end if

      print *,'  '
      print *,' EFFECTIVELY UNPAIRED ELECTRONS'
      print *,'  '
      print *,'    Atom     u_A'
      print *,' -----------------'
      call vprint(ua,nat,maxat,1)
      print *,' ------------------'
      write(*,'(a17,f10.5)') ' Sum check N_D = ' ,xnd
      if (idofr.eq.1) then
        line ='   FRAGMENT ANALYSIS : Num. eff. unpaired elec.'
        call group_by_frag_vec(1,line ,ua)
      end if

      WRITE(*,3)
 3    FORMAT(1x,/21X,'    APOST3D S^2 DECOMPOSITION (a=3/4)'//)
      call MPRINT(s234,nat,maxat)
      print *,' '
      write(*,'(a20,f10.5)') 'Sum check  <S^2> = ' ,x2
      if (idofr.eq.1) then
        line ='   FRAGMENT ANALYSIS : Local Spin Analysis'
        call group_by_frag_mat(0,line ,xlsa)
      end if

      end
