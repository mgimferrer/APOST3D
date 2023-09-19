
      subroutine spincorr(sat,dm1,dm2)
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      parameter (TOL=1.0d-10)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /localspin/xlsa(maxat,maxat),ua(maxat)
      common /iops/iopt(100)
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      dimension sat(igr,igr,nat)
      dimension dm1(nspinorb,nspinorb)
      dimension dm2(norb,norb,norb,norb)
      character*80 line

      dimension S234(maxat,maxat)
      dimension bonew(maxat,maxat)
      allocatable scr(:,:),scr2(:,:),satmo(:,:,:)
      allocatable sfdm1(:,:), psdm1(:,:),u_no(:)

      idofr= iopt(40)

      allocate (scr(igr,igr),scr2(igr,igr))
      allocate (satmo(igr,igr,nat))
      allocate (sfdm1(norb,norb),psdm1(norb,norb),u_no(norb))

c Compute eff. unpaired electrons in NO basis
c assuming same dimension as MOs
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
c print
c       write (*,*) 'atom1'
c       do i=1,norb
c         write(*,*) (satmo(i,j,1),j=1,norb)
c       end do
c       write (*,*) 'atom2'
c       do i=1,norb
c         write(*,*) (satmo(i,j,2),j=1,norb)
c       end do

      xnd=0.0d0
      do i=1,norb
        u_no(i)=2.0d0*occ_no(i,i)-occ_no(i,i)**2.0d0 
        xnd=xnd+u_no(i)
      end do 

c atom condensed
      do iat=1,nat
       xx=0.0d0
       do i=1,norb
         xx=xx+u_no(i)*satmo(i,i,iat)                   
       end do 
       ua(iat)=xx
      end do 

      do iat=1,nat
       do jat=iat,nat
        yy=0.0d0
        do i=1,norb
         do j=1,norb
c exchanged indices
c          yy=yy+occ_no(i,i)*occ_no(j,j)*satmo(i,j,iat)*satmo(i,j,jat)
          yy=yy+occ_no(i,i)*occ_no(j,j)*satmo(i,j,iat)*satmo(j,i,jat)
         end do 
        end do 
        xx=0.0d0
        do i=1,norb
         if(u_no(i).gt.TOL) then
          uno12=dsqrt(u_no(i))
          do j=1,norb
c exchanged indices
c          if(u_no(j).gt.TOL) xx=xx+uno12*dsqrt(u_no(j))*satmo(i,j,iat)*satmo(i,j,jat)
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
C done with eff. unpaired electrons and new bond order

C LSA
C Transforming now to the MO basis
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

C orthogonality check
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

C Initializing arrays for Local Spin, DI and U decomposition
c Local spin formulae:  (-1+2a)*Gamma_ijij - 0.5*Gamma_ijji
C a=0 --> Alcoba
C a=3/4 --> Ramos
c now recalculating BO
      do i=1,nat
       do j=1,nat
        s234(i,j)=0.0d0
        DI(i,j)=0.0d0
        BO(i,j)=0.0d0
       end do 
      end do

C P AND Ps in MO basis
      do i=1,norb
       do k=1,norb
        sfdm1(i,k)=dm1((i-1)*2+1,(k-1)*2+1)+dm1((i-1)*2+2,(k-1)*2+2)
        psdm1(i,k)=dm1((i-1)*2+1,(k-1)*2+1)-dm1((i-1)*2+2,(k-1)*2+2)
       end do 
      end do 

C CONTRIBUTIONS FROM THE SPINLESS CUMULANT OF THE DM2 in 1122 form
C spin density contributions included  in the cumulant (i.e. not removed)
      do i=1,norb
       do j=1,norb
        do k=1,norb
         do l=1,norb
          xx1=dm2(i,j,k,l)-sfdm1(i,j)*sfdm1(k,l)+(sfdm1(l,i)*sfdm1(j,k)+psdm1(l,i)*psdm1(j,k))/2.0d0
          xx0=dm2(i,j,k,l)-sfdm1(i,j)*sfdm1(k,l)+(sfdm1(l,i)*sfdm1(j,k))/2.0d0
          xx2=(sfdm1(l,i)*sfdm1(j,k))/2.0d0+(psdm1(l,i)*psdm1(j,k))/2.0d0
          do iat=1,nat
           do jat=iat,nat
c exchanged indices
c            S234(iat,jat)=s234(iat,jat)+xx0*0.5d0*(satmo(j,i,iat)*satmo(k,l,jat)-satmo(l,i,iat)*satmo(k,j,jat)) 
            S234(iat,jat)=s234(iat,jat)+xx0*0.5d0*(satmo(j,i,iat)*satmo(l,k,jat)-satmo(l,i,iat)*satmo(j,k,jat)) 
c exchanged indices 
c            DI(iat,jat)=DI(iat,jat)-xx1*2.0d0*(satmo(j,i,iat)*satmo(k,l,jat))
            DI(iat,jat)=DI(iat,jat)-xx1*2.0d0*(satmo(j,i,iat)*satmo(l,k,jat))
c exchanged indices 
c            BO(iat,jat)=BO(iat,jat)+xx2*2.0d0*(satmo(j,i,iat)*satmo(k,l,jat))
            BO(iat,jat)=BO(iat,jat)+xx2*2.0d0*(satmo(j,i,iat)*satmo(l,k,jat))
c            if(iat.eq.1.and.jat.eq.2) write(8,'(4i3,6f10.6)') j,i,l,k,satmo(j,i,iat),satmo(l,k,jat),xx2,xx1,xx2*2.0d0*(satmo(j,i,iat)*satmo(l,k,jat)),-xx1*2.0d0*(satmo(j,i,iat)*satmo(l,k,jat))
           end do
          end do
         end do
        end do
       end do
      end do

C CONTRIBUTIONS FROM THE  DM1
C simply from the density of eff. unpaired elec, already calculated

      do iat=1,nat
        S234(iat,iat)=s234(iat,iat)+0.75d0*ua(iat)
      end do
      deallocate (satmo,sfdm1,psdm1)

C symmetrizing and computing LI/DI
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
c saving Ramos decomposition
      do i=1,nat
       do j=1,nat
        xlsa(i,j)=s234(i,j)
       end do
      end do
C
C Sum checks
C
      x0=0.d0
      x1=0.d0
      x2=0.d0
      x3=0.d0
      do iat=1,nat
       x0=x0+di(iat,iat)
       x1=x1+bonew(iat,iat)
       do jat=1,nat
         x2=x2+s234(iat,jat)
         if(iat.ne.jat) x0=x0+di(iat,jat)/2.0d0
         if(iat.ne.jat) x1=x1+bonew(iat,jat)/2.0d0
       enddo
      enddo

C printing
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
      write(*,'(a16,f10.5)') ' Sum check N_D = ' ,xnd
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
