!! ----------------------------------------------------------------------
!! subroutine: spincorr
!! Local spin analysis of a correlated wavefunction from its 1- and
!! 2-RDM: the exchange-correlation bond orders (localization and
!! delocalization indices) from the 2-RDM, the effectively unpaired
!! electrons of each atom and the a=3/4 decomposition of <S^2>
!! (E. Ramos-Cordoba et al., J. Chem. Theory Comput. 2012, 8, 1270),
!! per atom and per fragment. The exchange-only bond orders of the
!! 1-RDM (bo, /ovpop/) come from fborder (pop.f); the 2-RDM cumulant
!! is added to them.
!! arguments:
!!   sat (in)  (igr,igr,nat) atomic overlap matrices S^A
!!   dm1 (in)  (nspinorb,nspinorb) 1-RDM, spin-orbital basis
!!   dm2 (in)  (norb,norb,norb,norb) spinless 2-RDM, natural orbitals
!! author: PSalse, ERaco
!! ----------------------------------------------------------------------
      subroutine spincorr(sat,dm1,dm2)
      use ao_matrices
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /printout/iaccur
      character*20 ckey,cspace
      character*60 cname
      character*12 cshort
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /localspin/xlsa(maxat,maxat),ua(maxat)
      common /iops/iopt(200)
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      dimension sat(igr,igr,nat)
      dimension dm1(nspinorb,nspinorb)
      dimension dm2(norb,norb,norb,norb)

      dimension S234(maxat,maxat)
      allocatable scr(:,:),satmo(:,:,:)
      allocatable sfdm1(:,:),psdm1(:,:),u_no(:)

      idofr = iopt(40)

      allocate(scr(igr,igr))
      allocate(satmo(igr,igr,nat))
      allocate(sfdm1(norb,norb),psdm1(norb,norb),u_no(norb))

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
            write(*,'(2x,a,i0,a,i0,a,f10.6)') 'WARNING: orbitals ',i,
     +        ' and ',j,' deviate from orthonormality by ',delta
          end if
        end do
      end do

!! initializing arrays for local spin, DI and U decomposition; local !!
!! spin formula (-1+2a)*Gamma_ijij - 0.5*Gamma_ijji, a=0 (Alcoba) vs !!
!! a=3/4 (Ramos-Cordoba) -- this routine implements a=3/4.           !!
      do i=1,nat
        do j=1,nat
          s234(i,j)=0.0d0
          DI(i,j)=0.0d0
        end do
      end do

!! P and Ps in the MO basis !!
      do i=1,norb
        do k=1,norb
          sfdm1(i,k)=dm1((i-1)*2+1,(k-1)*2+1)+dm1((i-1)*2+2,(k-1)*2+2)
          psdm1(i,k)=dm1((i-1)*2+1,(k-1)*2+1)-dm1((i-1)*2+2,(k-1)*2+2)
        end do
      end do

!! contributions from the spinless cumulant of dm2, 1122 form (spin   !!
!! density contributions included in the cumulant, i.e. not removed). !!
!! satmo is read with its indices exchanged (satmo(j,i,iat)) relative !!
!! to the naive iat/jat-symmetric form, consistently.                 !!
      do i=1,norb
        do j=1,norb
          do k=1,norb
            do l=1,norb
              xx1=dm2(i,j,k,l)-sfdm1(i,j)*sfdm1(k,l)+(sfdm1(l,i)*sfdm1(j,k)+psdm1(l,i)*psdm1(j,k))/2.0d0
              xx0=dm2(i,j,k,l)-sfdm1(i,j)*sfdm1(k,l)+(sfdm1(l,i)*sfdm1(j,k))/2.0d0
              do iat=1,nat
                do jat=iat,nat
                  S234(iat,jat)=s234(iat,jat)+xx0*0.5d0*(satmo(j,i,iat)*satmo(l,k,jat)-satmo(l,i,iat)*satmo(j,k,jat))
                  DI(iat,jat)=DI(iat,jat)-xx1*2.0d0*(satmo(j,i,iat)*satmo(l,k,jat))
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

!! symmetrizing; LI/DI = exchange-only bond orders (bo, localization !!
!! index halved on the diagonal) + cumulant part                      !!
      do i=1,nat
        do j=i+1,nat
          s234(j,i)=s234(i,j)
          DI(i,j)=DI(i,j)+bo(i,j)
          DI(j,i)=DI(i,j)
        end do
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
      x2=0.d0
      do iat=1,nat
        x0=x0+di(iat,iat)
        do jat=1,nat
          x2=x2+s234(iat,jat)
          if(iat.ne.jat) x0=x0+di(iat,jat)/2.0d0
        enddo
      enddo

!! printing !!
      call atomic_definition_names(ckey,cname,cshort,cspace)
      ndec=6
      if(iaccur.ne.0) ndec=13

      call print_box('EXCHANGE-CORRELATION BOND ORDER MATRIX (DI)')
      write(*,'(2x,a)') 'From the 2-RDM: delocalization indices, '//
     +  'localization indices on the diagonal'
      write(*,'(a)') ''
      call mprint2(di,nat,maxat)
      call print_value('Total (diagonal and one triangle)',x0,ndec)
      if(idofr.eq.1) call group_by_frag_mat(1,
     +  'FRAGMENT EXCHANGE-CORRELATION BOND ORDERS (DI)',di)

      call print_atom_table('EFFECTIVELY UNPAIRED ELECTRONS',' ',' ',ua)
      call print_value('Total unpaired electrons N_D',xnd,ndec)
      if(idofr.eq.1) call group_by_frag_vec(1,
     +  'FRAGMENT UNPAIRED ELECTRONS',ua,cshort,' ')

      call print_box('S^2 DECOMPOSITION (a=3/4)')
      call mprint2(s234,nat,maxat)
      call print_value('Total <S^2>',x2,ndec)
      if(idofr.eq.1) call group_by_frag_mat(0,
     +  'FRAGMENT S^2 DECOMPOSITION (a=3/4)',xlsa)

      call json_local_spin(.true.,xnd,x2,di,x0)

      end

!! ----------------------------------------------------------------------
!! subroutine: bond_orders_no
!! Bond orders of a correlated wavefunction from its natural orbitals,
!! I. Mayer, Chem. Phys. Lett. 2012, 554, 83:
!!   B_AB = sum_ij [n_i n_j + sqrt(u_i u_j)] S^A_ij S^B_ji,  u = 2n - n^2
!! (S^A in the natural-orbital basis, half of it on the diagonal), per
!! atom and per fragment. Needs the natural orbitals of the 1-RDM
!! (dm1input).
!! arguments:
!!   sat    (in)  (igr,igr,nat) atomic overlap matrices S^A
!!   bonew  (out) (maxat,maxat) the bond orders
!!   xbonew (out) their total, diagonal and one triangle
!! author: PSalse
!! ----------------------------------------------------------------------
      subroutine bond_orders_no(sat,bonew,xbonew)
      use ao_matrices
      use output_mod, only: print_value
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /printout/iaccur
      parameter (TOL=1.0d-10)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      dimension sat(igr,igr,nat)
      dimension bonew(maxat,maxat)
      allocatable scr(:,:),satmo(:,:,:),u_no(:)

      idofr = iopt(40)

      allocate(scr(igr,igr))
      allocate(satmo(igr,igr,nat))
      allocate(u_no(norb))

!! S^A in the NO basis !!
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

      do i=1,norb
        u_no(i)=2.0d0*occ_no(i,i)-occ_no(i,i)**2.0d0
      end do

!! satmo(i,j,jat) is read with its indices exchanged (satmo(j,i,jat)) !!
!! relative to the naive iat/jat-symmetric form                       !!
      xbonew=0.0d0
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
          xbonew=xbonew+bonew(iat,jat)
        end do
      end do
      deallocate(scr,satmo,u_no)

      ndec=6
      if(iaccur.ne.0) ndec=13
      call print_box('NEW BOND ORDER MATRIX')
      write(*,'(2x,a)') 'From the natural orbitals, I. Mayer, '//
     +  'Chem. Phys. Lett. 2012, 554, 83'
      write(*,'(a)') ''
      call mprint2(bonew,nat,maxat)
      call print_value('Total (diagonal and one triangle)',xbonew,ndec)
      if(idofr.eq.1) call group_by_frag_mat(1,
     +  'FRAGMENT NEW BOND ORDERS',bonew)

      end
