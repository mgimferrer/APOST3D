!! *********************************************************************** !!
!! OS-CENTROID: oxidation states from the centroids of localized          !!
!! orbitals, each orbital assigned to the atom closest to its centroid.   !!
!!   eos_centroid -- the whole analysis (main.f stops right after it)     !!
!! *********************************************************************** !!

      subroutine eos_centroid(itotps,chp,wp,omp,pcoord)
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      parameter (toau=0.52917721067d0)
      integer, intent(in) :: itotps
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(200)
      common /atomrad/atr(maxat),distance(maxat,maxat)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /nna/ xnonradi(maxnna),atsphradi(maxat),nna
      dimension chp(itotps,igr), pcoord(itotps,3),elect(maxat,4)
      dimension omp(itotps), wp(itotps)
      real*8 hess(3,3),vecp(3,3),hessi(3,3),xgrad(3),xstep(3)
      character*80 line
      allocatable chp2(:)

      idofr=iopt(40)
      iatps=nang*nrad
c
c closed-shell only
      allocate(chp2(itotps))
c
      do i=1,icufr
       do k=1,4
       elect(i,k)=0.0d0
       end do
      end do
c print atomic coordinates for reference
      write(*,*) 'Molecular coordinates'
      do iat=1,nat
       write(*,'(i3,3f12.8)') iznuc(iat),(toau*coord(i,iat),i=1,3)
      end do
      write(*,*) 

c getting beta-sphere
      call atsphere()

c correct for pseudopotential 
      do iat=1,nat
       if(iznuc(iat).ne.int(zn(iat))) then
        write(*,*) 'Warning: pseudopotential on atom ',iat,atr(iat)
        write(*,*) 'Setting trust sphere to (a.u.) :',atr(iat)/2.0d0
        atsphradi(iat)=atr(iat)/2.0d0
       end if
      end do
      write(*,*) 'Atomic trust sphere'
      call vprint(atsphradi,nat,maxat,1)
      write(*,*) 

      ndim=nocc
      nnelec=2
      if(kop.eq.1) then
        ndim=nalf
        nnelec=1 
      end if
c
c Loop over localized mo's
c
      do iorb=1,ndim
       do k=1,itotps 
        xx=0.0d0
        do i=1,igr
          xx=xx+c(i,iorb)*chp(k,i)
        end do
        chp2(k)=xx
       end do
       x=0.0d0
       y=0.0d0
       z=0.0d0
       xx=0.0d0
       yy=0.0d0
       zz=0.0d0
       zzz=0.0d0
       do icenter=1,nat
        do ifut=iatps*(icenter-1)+1,iatps*icenter
         x=x+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,1)
         y=y+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,2)
         z=z+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,3)
         xx=xx+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,1)**2
         yy=yy+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,2)**2
         zz=zz+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,3)**2
        end do
       enddo
       zzz=sqrt(xx+yy+zz-x**2-y**2-z**2)
       write(*,*)
       write(*,'(a9,i4,a4,f8.4,a3,3f12.4)') 'Centroid',iorb,' : s',zzz*toau,' X',x*toau,y*toau,z*toau
c
c closest atom to centroid
c
       dist0=1.0d99
       iiat=0
       do iat=1,nat
        dist=(x-coord(1,iat))**2.0d0+(y-coord(2,iat))**2.0d0+(z-coord(3,iat))**2.0d0
        if(dist.lt.dist0) then
         dist0=dist
         iiat=iat
        end if
       end do
       elect(jfrlist(iiat),1)=elect(jfrlist(iiat),1)-nnelec
       write(*,'(a16,i3,a6,i2,a17,f8.4)') 'Closest center: ',iiat,
     + ' Frag.',jfrlist(iiat),' Distance (a.u.):',dsqrt(dist0)*toau
c
c Searching for atomic basin of centroid
C
c using steepest ascent from centroid
       iter=1
       maxiter=200
       do while(dsqrt(dist0).gt.atsphradi(iiat).and.iter.le.maxiter)
c calcualte hessian and gradient
        do ixyz=1,3
         do jxyz=ixyz,3
          hess(ixyz,jxyz)=d2functxyz(ixyz,jxyz,x,y,z)
          if(ixyz.ne.jxyz) hess(jxyz,ixyz)=hess(ixyz,jxyz)
         end do
         xgrad(ixyz)=dfunctxyz(ixyz,x,y,z)
         zz=zz+xgrad(ixyz)*xgrad(ixyz)
        end do
c        write(*,'(a10,3f12.5,a6,f12.5)') 'Gradient:',(xgrad(i),i=1,3),
c     + 'Norm:',dsqrt(zz)
c        write(*,*) 'HESSIAN '
c        do i=1,3
c         write(*,'(3f12.6)') (hess(i,j),j=1,3)
c        end do
        call diagonalize(3,3,hess,vecp,0)
c reverse eigenvalue to do steepest ascent
        do i=1,3
         do j=1,3
           hessi(i,j)=0.0d0
           do k=1,3
            hessi(i,j)=hessi(i,j)+vecp(i,k)*vecp(j,k)/abs(hess(k,k))
           end do
         end do
        end do
c        call invert(3,hess)
        zz=0.0d0
        do i=1,3
          xstep(i)=0.0d0
         do j=1,3
          xstep(i)=xstep(i)+hessi(i,j)*xgrad(j) 
c          xstep(i)=xstep(i)+hess(i,j)*xgrad(j) 
         end do
         zz=zz+ xstep(i)*xstep(i)
        end do
        write(*,'(a10,3f12.5,a5,f12.5)') 'Step    :',(xstep(i),i=1,3),
     +  ' |s|:',dsqrt(zz)
        step=1.0d0
        if(dsqrt(zz).gt.0.3) step=0.30d0/dsqrt(zz)
c steepest ascent
        x1=x+xstep(1)*step
        y1=y+xstep(2)*step
        z1=z+xstep(3)*step
c check closest atom
        dist1=1.0d99
        iiat1=0
        do iat=1,nat
         dist=(x1-coord(1,iat))**2.0d0+(y1-coord(2,iat))**2.0d0+(z1-coord(3,iat))**2.0d0
         if(dist.lt.dist1) then
          dist1=dist
          iiat1=iat
         end if
        end do
        if(iiat1.ne.iiat) then
         iiat=iiat1
         write(*,*)'Moving away from closest atom. Actual closest atom:',iiat 
        end if
        dist0=dist1
        x=x1
        y=y1
        z=z1
        iter=iter+1
        write(*,'(a12,3f12.5,a12,i3,f12.5)') 'new pos. X',x*toau,y*toau,z*toau,' dist to at:',iiat,dsqrt(dist0)*toau
       end do
       if(iter.gt.maxiter) then
        write(*,*) 'steepest ascent did not converge'
        call apost_stop('')
       else
        write(*,'(a16,i3,a6,i2,a17,f8.4,a6,i2)') 'QTAIM center:   ',iiat,
     + ' Frag.',jfrlist(iiat),' Distance (a.u.):',dsqrt(dist0)*toau,
     + ' iter ',iter
        elect(jfrlist(iiat),2)=elect(jfrlist(iiat),2)-nnelec
       end if
c End loop over localized mo's
      end do

      if(kop.eq.1) then
      itot1=0
      itot2=0
      print *,'   OXIDATION STATES FROM    '
      print *,' CENTROIDS OF LOCALIZED MOs '
      write(*,*)
      print *,'        ALPHA PART'
      write(*,*)
      print *,'  Frag  Na(geom)  Na(QTAIM) '
      print *,' ---------------------------'
      do i=1,icufr
       itot1=itot1-elect(i,1)
       itot2=itot2-elect(i,2)
       write(*,'(3x,i3,5x,i4,5x,i4)') i,-int(elect(i,1)),-int(elect(i,2))
      end do
      print *,' ---------------------------'
      write(*,'(a6,2i9)') '  Sum ', itot1,itot2
      write(*,*)

c run over beta localized mo's
      do iorb=1,nb   
       do k=1,itotps 
        xx=0.0d0
        do i=1,igr
          xx=xx+cb(i,iorb)*chp(k,i)
        end do
        chp2(k)=xx
       end do
       x=0.0d0
       y=0.0d0
       z=0.0d0
       xx=0.0d0
       yy=0.0d0
       zz=0.0d0
       zzz=0.0d0
       do icenter=1,nat
        do ifut=iatps*(icenter-1)+1,iatps*icenter
         x=x+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,1)
         y=y+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,2)
         z=z+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,3)
         xx=xx+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,1)**2
         yy=yy+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,2)**2
         zz=zz+wp(ifut)*omp(ifut)*chp2(ifut)**2*pcoord(ifut,3)**2
        end do
       enddo
       zzz=sqrt(xx+yy+zz-x**2-y**2-z**2)
       write(*,*)
       write(*,'(a9,i4,a4,f8.4,a3,3f12.4)') 'Centroid',iorb,' : s',zzz*toau,' X',x*toau,y*toau,z*toau

c closest atom to centroid
       dist0=1.0d99
       iiat=0
       do iat=1,nat
        dist=(x-coord(1,iat))**2.0d0+(y-coord(2,iat))**2.0d0+(z-coord(3,iat))**2.0d0
        if(dist.lt.dist0) then
         dist0=dist
         iiat=iat
        end if
       end do
       elect(jfrlist(iiat),3)=elect(jfrlist(iiat),3)-1
       write(*,'(a16,i3,a6,i2,a17,f8.4)') 'Closest center: ',iiat,
     + ' Frag.',jfrlist(iiat),' Distance (a.u.):',dsqrt(dist0)*toau
c using steepest ascent from centroid
       iter=1
       maxiter=200
       do while(dsqrt(dist0).gt.atsphradi(iiat).and.iter.le.maxiter)
c calcualte hessian and gradient
        do ixyz=1,3
         do jxyz=ixyz,3
          hess(ixyz,jxyz)=d2functxyz(ixyz,jxyz,x,y,z)
          if(ixyz.ne.jxyz) hess(jxyz,ixyz)=hess(ixyz,jxyz)
         end do
         xgrad(ixyz)=dfunctxyz(ixyz,x,y,z)
         zz=zz+xgrad(ixyz)*xgrad(ixyz)
        end do
c        write(*,'(a10,3f12.5,a6,f12.5)') 'Gradient:',(xgrad(i),i=1,3),
c     + 'Norm:',dsqrt(zz)
c        write(*,*) 'HESSIAN '
c        do i=1,3
c         write(*,'(3f12.6)') (hess(i,j),j=1,3)
c        end do
        call diagonalize(3,3,hess,vecp,0)
c reverse eigenvalue to do steepest ascent
        do i=1,3
         do j=1,3
           hessi(i,j)=0.0d0
           do k=1,3
            hessi(i,j)=hessi(i,j)+vecp(i,k)*vecp(j,k)/abs(hess(k,k))
           end do
         end do
        end do
c        call invert(3,hess)
        zz=0.0d0
        do i=1,3
          xstep(i)=0.0d0
         do j=1,3
          xstep(i)=xstep(i)+hessi(i,j)*xgrad(j) 
c          xstep(i)=xstep(i)+hess(i,j)*xgrad(j) 
         end do
         zz=zz+ xstep(i)*xstep(i)
        end do
        write(*,'(a10,3f12.5,a5,f12.5)') 'Step    :',(xstep(i),i=1,3),
     +  ' |s|:',dsqrt(zz)
        step=1.0d0
        if(dsqrt(zz).gt.0.3) step=0.30d0/dsqrt(zz)
c steepest ascent
        x1=x+xstep(1)*step
        y1=y+xstep(2)*step
        z1=z+xstep(3)*step
c check closest atom
        dist1=1.0d99
        iiat1=0
        do iat=1,nat
         dist=(x1-coord(1,iat))**2.0d0+(y1-coord(2,iat))**2.0d0+(z1-coord(3,iat))**2.0d0
         if(dist.lt.dist1) then
          dist1=dist
          iiat1=iat
         end if
        end do
        if(iiat1.ne.iiat) then
         iiat=iiat1
         write(*,*)'Moving away from closest atom. Actual closest atom:',iiat 
        end if
        dist0=dist1
        x=x1
        y=y1
        z=z1
        iter=iter+1
        write(*,'(a12,3f12.5,a12,i3,f12.5)') 'new pos. X',x*toau,y*toau,z*toau,' dist to at:',iiat,dsqrt(dist0)*toau
       end do
       if(iter.gt.maxiter) then
        write(*,*) 'steepest ascent did not converge'
        call apost_stop('')
       else
        write(*,'(a16,i3,a6,i2,a17,f8.4,a6,i2)') 'QTAIM center:   ',iiat,
     + ' Frag.',jfrlist(iiat),' Distance (a.u.):',dsqrt(dist0)*toau,
     + ' iter ',iter
        elect(jfrlist(iiat),4)=elect(jfrlist(iiat),4)-nnelec
       end if
c End loop over localized mo's
      end do

      itot1=0
      itot2=0
      write(*,*)
      print *,'   OXIDATION STATES FROM    '
      print *,' CENTROIDS OF LOCALIZED MOs '
      write(*,*)
      print *,'         BETA PART'
      write(*,*)
      print *,'  Frag  Nb(geom)  Nb(QTAIM) '
      print *,' ---------------------------'
      do i=1,icufr
       itot1=itot1-elect(i,3)
       itot2=itot2-elect(i,4)
       write(*,'(3x,i3,5x,i4,5x,i4)') i,-int(elect(i,3)),-int(elect(i,4))
       elect(i,1)=int(elect(i,1))+int(elect(i,3))
       elect(i,2)=int(elect(i,2))+int(elect(i,4))
      end do
      print *,' ---------------------------'
      write(*,'(a6,2i9)') '  Sum ', itot1,itot2
      write(*,*)

      end if

c OS
      do iat=1,nat 
       elect(jfrlist(iat),1)=elect(jfrlist(iat),1)+zn(iat)
       elect(jfrlist(iat),2)=elect(jfrlist(iat),2)+zn(iat)
      end do

      itot1=0
      itot2=0
      write(*,*)
      print *,'   OXIDATION STATES FROM    '
      print *,' CENTROIDS OF LOCALIZED MOs '
      write(*,*)
      print *,'  Frag  OS(geom)  OS(QTAIM) '
      print *,' ---------------------------'
      do i=1,icufr
       itot1=itot1+elect(i,1)
       itot2=itot2+elect(i,2)
       write(*,'(3x,i3,5x,i4,5x,i4)') i,int(elect(i,1)),int(elect(i,2))
      end do
      print *,' ---------------------------'
      write(*,'(a6,2i9)') '  Sum ', itot1,itot2
      write(*,*)

      end




