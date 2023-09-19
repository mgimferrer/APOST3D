      subroutine tomull(sat)
      use basis_set
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      dimension sat(nbasis,nbasis,natoms)

c            cross-check
       do i=1,natoms
        do nu=1,nbasis
          do mu=1,nbasis
           if (ihold(mu).eq.i)then
            sat(nu,mu,i)=s(nu,mu)
           else
            sat(nu,mu,i)=0.0d0
           end if
          end do
         end do
       end do

       end

      subroutine tolow(sat)
      use basis_set
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /iops/iopt(100)
      dimension sat(nbasis,nbasis,natoms)
      allocatable  ss(:,:),s12(:,:),sm12(:,:),sp12(:,:),x(:,:) 

      allocate  (ss(nbasis,nbasis),s12(nbasis,nbasis), sm12(nbasis,nbasis))
      allocate  (sp12(nbasis,nbasis),x(nbasis,nbasis))

      imulli= Iopt(5)
      idav=0
      if(imulli.eq.3) idav=1

c doing davidson
      if(idav.eq.1) then
       do i=1,nbasis
        do j=1,nbasis
         ss(i,j)=0.0d0
        end do
       end do
       do i=1,natoms
        do mu=1,nbasis
         if (ihold(mu).eq.i)then
          do nu=1,nbasis
           if (ihold(nu).eq.i)then
           ss(mu,nu)=s(mu,nu)
           end if
          end do
         end if
        end do
       end do

       call diagonalize(nbasis,nbasis,ss,x,0)

       do i=1,nbasis
        do j=1,nbasis
         sm12(i,j)=0.0d0
         sp12(i,j)=0.0d0
         do k=1,nbasis
          sm12(i,j)=sm12(i,j)+x(i,k)*x(j,k)/dsqrt(ss(k,k))
          sp12(i,j)=sp12(i,j)+x(i,k)*x(j,k)*dsqrt(ss(k,k))
         end do
        end do
       end do
c
       do i=1,nbasis
        do j=1,nbasis
         ss(i,j)=0.0d0
         do k=1,nbasis
          ss(i,j)=ss(i,j)+sm12(k,i)*S(k,j)
         end do
        end do
       end do
       do i=1,nbasis
        do j=1,nbasis
         s12(i,j)=0.0d0
         do k=1,nbasis
          s12(i,j)=s12(i,j)+ss(i,k)*sm12(k,j)
         end do
        end do
       end do

       do i=1,nbasis
        do j=1,nbasis
         ss(i,j)=s12(i,j)
        end do
       end do

      else
! conventional Lowdin

       do i=1,nbasis
        do j=1,nbasis
         ss(i,j)=s(i,j)
        end do
       end do
      end if

       call diagonalize(nbasis,nbasis,ss,x,0)
       do i=1,nbasis
        do j=1,nbasis
         s12(i,j)=0.0d0
         do k=1,nbasis
          s12(i,j)=s12(i,j)+x(i,k)*dsqrt(ss(k,k))*x(j,k)
         end do
        end do
       end do

       if(idav.eq.1) then
        do i=1,nbasis
         do j=1,nbasis
          ss(i,j)=0.0d0
          do k=1,nbasis
           ss(i,j)=ss(i,j)+sp12(i,k)*s12(k,j)
          end do
         end do
        end do
       else
        do i=1,nbasis
         do j=1,nbasis
          ss(i,j)=s12(i,j)
          end do
         end do
        end if
c            cross-check
       do i=1,natoms
        do nu=1,nbasis
         do nu1=1,nbasis
          sat(nu,nu1,i)=0.0d0
          do mu=1,nbasis
           if (ihold(mu).eq.i)then
            sat(nu,nu1,i)=sat(nu,nu1,i)+ss(nu,mu)*ss(nu1,mu)
           end if
          end do
         end do
        end do
       end do
      deallocate  (ss,s12,sm12,sp12,x)

       end


      subroutine mull_opop()
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq

      do iatom=1,natoms
       do jatom=iatom,natoms
        xx=0.0d0
        do i=llim(iatom),iulim(iatom)
         do j=llim(jatom),iulim(jatom)
          xx=xx+p(i,j)*s(j,i)
         end do
        end do
        op(iatom,jatom)=xx
        if(iatom.ne.jatom) op(jatom,iatom)=op(iatom,jatom)
       end do
      end do
      return
      end
!
      subroutine tonao(sat)
      use basis_set
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /iops/iopt(100)
      common /filename/name0
      common/nao/unao(nmax,nmax),ssnao(nmax,nmax)
      character*60 name0,name3
      dimension sat(nbasis,nbasis,natoms)

c reading NAOAO matrix from FILE.33
c generate file with $NBO AONAO=W $END keyword

      name3=trim(name0)//'.nao'
      open(33,file=name3,status='old')
      write(*,*) 'Reading NAO to AO matrix'
      read(33,*)
      read(33,*)
      read(33,*)
      do j=1,nbasis
       read(33,*)(unao(i,j),i=1,nbasis)
      end do
      write(*,*) 'NAO to AO matrix read'
      do nu=1,nbasis
       do mu=1,nbasis
         ssnao(nu,mu)=0.0d0
         do k=1,nbasis
          ssnao(nu,mu)=ssnao(nu,mu)+unao(k,nu)*s(k,mu)
         end do
       end do
      end do
c cosntruct sat
      do i=1,natoms
       do nu=1,nbasis
        do mu=1,nbasis
         sat(nu,mu,i)=0.0d0
         do k=1,nbasis
          if(ihold(k).eq.i) sat(nu,mu,i)=sat(nu,mu,i)+ssnao(k,nu)*ssnao(k,mu)
         end do
       end do
       end do
      end do

      end

!
! testing weigthed-Lowdin
!
      subroutine tolow2(sat)
      use basis_set
      use ao_matrices, only:p 
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /iops/iopt(100)
      dimension sat(nbasis,nbasis,natoms)
      allocatable  ss(:,:),s12(:,:),sm12(:,:),sp12(:,:),x(:,:),wlow(:)

      allocate  (ss(nbasis,nbasis),s12(nbasis,nbasis), sm12(nbasis,nbasis))
      allocate  (sp12(nbasis,nbasis),x(nbasis,nbasis),wlow(nbasis))

c weigthed-l lowdin
      do mu=1,nbasis
        xx=0.0d0
        xx1=0.0d0
        iatom=ihold(mu)
        do nu=1,nbasis
         xx=xx+p(mu,nu)*s(nu,mu)
        end do
        do nu=1,nbasis
         if (ihold(nu).eq.iatom)then
          xx1=xx1+p(mu,nu)*s(nu,mu)
         end if
        end do
c mapping from [0,2] to [1,Rmax]
        Rmax=20.0d0
        aa=2.0d0/Rmax-1.0d0
        if(xx.lt.0.0d0) then
         xx2=1.0d0
        else if(xx.gt.2.0d0) then
         xx2=Rmax
        else
         xx2=(2.0d0+xx)/(2.0d0+aa*xx)
        end if
        wlow(mu)=xx2
      end do 

      do i=1,nbasis
       do j=1,nbasis
        ss(i,j)=s(i,j)*wlow(i)*wlow(j)
       end do
      end do

       call diagonalize(nbasis,nbasis,ss,x,0)
       do i=1,nbasis
        do j=1,nbasis
         s12(i,j)=0.0d0
         do k=1,nbasis
         s12(i,j)=s12(i,j)+x(i,k)*dsqrt(ss(k,k))*x(j,k)
        end do
       end do
       end do

c multiply by inverse of weights
       do i=1,nbasis
        do j=1,nbasis
         ss(i,j)=s12(i,j)/wlow(j)
        end do
       end do

C Build proper Sat
       do i=1,natoms
        do nu=1,nbasis
         do mu=1,nbasis
          sat(nu,mu,i)=0.0d0
          do k=1,nbasis
           if (ihold(k).eq.i)then
            sat(nu,mu,i)=sat(nu,mu,i)+ss(k,nu)*ss(k,mu)
           end if
          end do
         end do
        end do
       end do
      deallocate  (ss,s12,sm12,sp12,x)

      end
