C WARNING. still using old common's with basis functions
C needs to be updated

      subroutine gpoints(x,y,z,gxchp,gychp,gzchp,aochp)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common/actual/iact,jat,icenter
      common /nat/ nat,igr,ifg,idum(4)
      dimension gxchp(igr),gychp(igr),gzchp(igr),aochp(igr)

      do i=1,igr
       iact=i
       gxchp(iact)=gxfunct(x,y,z)
       gychp(iact)=gyfunct(x,y,z)
       gzchp(iact)=gzfunct(x,y,z)
       aochp(iact)=aofunct(iact,x,y,z)
      end do
      return
      end

      function gxfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /stv/ sp(nmax,nmax),tt(nmax,nmax)
      common /nat/ nat,igr,ifg,idum(4)
      common /data/expp(maxp),n(maxp),l(maxp),m(maxp),iat(maxp), ifill(nmax),ifiul(nmax)
      common /coeff/coeff(maxp,maxc)
      common /lim/ llim(nmax),iulim(nmax),ihold(nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/actual/iact,jact,icenter
      
      iactat=ihold(iact)
      f=0.d0
      xx=x-coord(1,iactat)   
      yy=y-coord(2,iactat)   
      zz=z-coord(3,iactat)   
      rr=dsqrt(xx**2+yy**2+zz**2)
      ngri=ifiul(iact)-ifill(iact)+1
      do ki=1,ngri
      ipr=ifill(iact)+ki-1
      nn=n(ipr)
      ll=l(ipr)
      mm=m(ipr)
c      if (abs(xx).lt.1.E-10.and.nn.eq.0)then 
      if (nn.eq.0)then 
       xprod=0.0d0
      else
       xprod=nn*xx**(nn-1)
      end if
      f=f+(xprod-2.d0*expp(ipr)*xx**(nn+1))*(yy**ll)*(zz**mm)*
     $dexp(-expp(ipr)*(rr**2))*coeff(iact,ki)
c      print *,iact,ki,ipr,expp(ipr),coeff(iact,ki)
      enddo

      gxfunct=f
      
      return
      end

      function aofunct(iact,x,y,z)
      use basis_set
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,idum(4)
      
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
      aofunct=f
      
      return
      end
      function gyfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /stv/ sp(nmax,nmax),tt(nmax,nmax)
      common /nat/ nat,igr,ifg,idum(4)
      common /data/expp(maxp),n(maxp),l(maxp),m(maxp),iat(maxp), ifill(nmax),ifiul(nmax)
      common /coeff/coeff(maxp,maxc)
      common /lim/ llim(nmax),iulim(nmax),ihold(nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/actual/iact,jact,icenter
      
      iactat=ihold(iact)
      f=0.d0
c      write(*,*)'coordfunctbef',x,y,z
      xx=x-coord(1,iactat)   
      yy=y-coord(2,iactat)   
      zz=z-coord(3,iactat)   
      rr=dsqrt(xx**2+yy**2+zz**2)
      ngri=ifiul(iact)-ifill(iact)+1
      do ki=1,ngri
      ipr=ifill(iact)+ki-1
      nn=n(ipr)
      ll=l(ipr)
      mm=m(ipr)
c      write(*,*)'coordfunct',xx,yy,zz
      if (ll.eq.0)then 
       xprod=0.0d0
      else
       xprod=ll*yy**(ll-1)     
      end if
      f=f+(xprod-2.d0*expp(ipr)*yy**(ll+1))*(xx**nn)*(zz**mm)*
     $dexp(-expp(ipr)*(rr**2))*coeff(iact,ki)
c      print *,iact,ki,ipr,expp(ipr),coeff(iact,ki)
c      print *,'f',iact,ki,f
      enddo

c      print *,'f',iact,f
      gyfunct=f
      
      return
      end

      function gzfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /stv/ sp(nmax,nmax),tt(nmax,nmax)
      common /nat/ nat,igr,ifg,idum(4)
      common /data/expp(maxp),n(maxp),l(maxp),m(maxp),iat(maxp),ifill(nmax),ifiul(nmax)
      common /coeff/coeff(maxp,maxc)
      common /lim/ llim(nmax),iulim(nmax),ihold(nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/actual/iact,jact,icenter
      
      iactat=ihold(iact)
      f=0.d0
      xx=x-coord(1,iactat)   
      yy=y-coord(2,iactat)   
      zz=z-coord(3,iactat)   
      rr=dsqrt(xx**2+yy**2+zz**2)
      ngri=ifiul(iact)-ifill(iact)+1
      do ki=1,ngri
      ipr=ifill(iact)+ki-1
      nn=n(ipr)
      ll=l(ipr)
      mm=m(ipr)
      if (mm.eq.0)then 
       xprod=0.0d0
      else
       xprod=mm*zz**(mm-1)
      end if
      f=f+(xprod-2.d0*expp(ipr)*zz**(mm+1))*(yy**ll)*(xx**nn)*
     $dexp(-expp(ipr)*(rr**2))*coeff(iact,ki)
c      print *,iact,ki,ipr,expp(ipr),coeff(iact,ki)
      enddo

      gzfunct=f
      
      return
      end


        subroutine densvalue (aochp,xdens)
        use ao_matrices
        implicit double precision (a-h,o-z)
        include 'parameter.h'
        common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
        dimension aochp(igr)
        
        xdens=0.0d0
        do i=1,igr
         do j=1,igr
          xdens=xdens+P(i,j)*aochp(i)*aochp(j)
         end do 
        end do 
        
        end 
        
        subroutine gradvalue(aochp,gxchp,gychp,gzchp,gx,gy,gz)
        use ao_matrices
        implicit double precision (a-h,o-z)
        include 'parameter.h'
        common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
        dimension aochp(igr),gxchp(igr),gychp(igr),gzchp(igr)
    
        gx=0.0d0
        gy=0.0d0
        gz=0.0d0
        do i=1,igr
         do j=1,igr
          gx=gx+P(i,j)*2*gxchp(i)*aochp(j)
          gy=gy+P(i,j)*2*gychp(i)*aochp(j)
          gz=gz+P(i,j)*2*gzchp(i)*aochp(j)
         end do 
        end do 
        end
        subroutine qtaimgrid(pcoord,nrad,nang,omp,ibaspoint,xlap)
        implicit double precision (a-h,o-z)
        include 'parameter.h'
        common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
        common /iops/iopt(100)
        common /nna/ xnonradi(maxnna),atsphradi(maxat),nna
        common /coordnon/ coordnon(3,maxnna)
        dimension omp(nat*nrad*nang)
        dimension pcoord(nrad*nang*nat,3)
        dimension ibaspoint(nat*nrad*nang)
        dimension xlap(nat*nrad*nang)

        dimension gxchp(:),gychp(:),gzchp(:),aochp(:)
        dimension hxxchp(:),hyychp(:),hzzchp(:)
        dimension hxychp(:),hxzchp(:),hzychp(:)
        allocatable gxchp,gychp,gzchp,aochp
        allocatable hxychp,hxzchp,hzychp
        allocatable hxxchp,hyychp,hzzchp

        dimension x0(3),x1(3),g0(3),H0(3,3),hh0(3,3)
        dimension datdist(nat)
        dimension datdist2(nat)
        real*8 maxdist


c  allocatables
        allocate(gxchp(igr))
        allocate(gychp(igr))
        allocate(gzchp(igr))
        allocate(aochp(igr))
        allocate(hxxchp(igr))
        allocate(hyychp(igr))
        allocate(hzzchp(igr))
        allocate(hxychp(igr))
        allocate(hxzchp(igr))
        allocate(hzychp(igr))

        istep       = Iopt(35) 
        step0=real(istep)/1000.0d0  
        inna        = Iopt(36) 
        imaxdist    = Iopt(37) 
        maxdist=real(imaxdist)
        iscreening  = Iopt(38) 
        screen = real(iscreening)/1D10
        ipath       = Iopt(39) 

        nna=0
       do ipoint=1,nat*nrad*nang
      stepe=step0
      ibas=0
        x0(1)=pcoord(ipoint,1) 
        x0(2)=pcoord(ipoint,2)
        x0(3)=pcoord(ipoint,3)

      ydist=maxdist
      do i=1,nat
       call dist(x0,i,xdist)
       datdist(i)=xdist
       if(xdist.lt.ydist)then 
        ydist=xdist
       end if 
       if (xdist.le.atsphradi(i))then 
        ibaspoint(ipoint)=i
        ibas=1
       end if
      end do
      
      
      if (ydist.ge.maxdist.and.ibas.eq.0) then
       ibaspoint(ipoint)=0
       ibas=1
      else  

      call gpoints(pcoord(ipoint,1),pcoord(ipoint,2),pcoord(ipoint,3),
     $  gxchp,gychp,gzchp,aochp)
       call densvalue(aochp,xdens0)
       xinipoint=xdens0*omp(ipoint)
       call gradvalue(aochp,gxchp,gychp,gzchp,gx,gy,gz)
       gmod=dsqrt(gx**2.0d0+gy**2.0d0+gz**2.0d0)
c        write(*,*) x0(1),x0(2),x0(3)
c        write(*,*)ipoint, xdens0*omp(ipoint), gmod, ydist
c        if (xdens0*omp(ipoint).lt.screen.and.gmod.lt.screen) then 
        if (xdens0*omp(ipoint).lt.screen) then 
c       if (xdens0*omp(ipoint).lt.1E-7) then 
        ibaspoint(ipoint)=0
        ibas=1  
c        xlap(ipoint)=0.0d0
       else 
      
        g0(1)=gx
        g0(2)=gy
        g0(3)=gz

       call lappoints(x0(1),x0(2),x0(3),hxxchp,hyychp,hzzchp)
             call laphessvalue(aochp,gxchp,gychp,gzchp,hxxchp,hyychp,
     $      hzzchp,h0)
          xlap(ipoint)=h0(1,1)+h0(2,2)+h0(3,3)
      
        
        do while(ibas.eq.0)

c NONNUCLEAR ATRACTORS   
c ERC only if necessary           
           if(inna.eq.1)then 
          if(xdens0.gt.0.001d0.and.gmod.le.1e-8)then  
             call hpoints(x0(1),x0(2),x0(3),hxxchp,hyychp,hzzchp,hxychp,
     $            hxzchp,hzychp)
             call hessvalue(aochp,gxchp,gychp,gzchp,hxxchp,hyychp,
     $       hzzchp,hxychp,hxzchp,hzychp,h0)
             call diagonalize(3,3,h0,hh0,0)
           if(h0(1,1).lt.0.0d0.and.h0(2,2).lt.0.0d0.and.h0(3,3).lt.0.0d0)then
            nna=nna+1
            do ij=1,3
             coordnon(ij,nna)=x0(ij)
            end do 
            write(*,*)'NONNUCLEARATRACTOR',nna,'position: ',x0(1),x0(2),x0(3)
            write(*,*)'hessian eigenvalues',h0(1,1),h0(2,2),h0(3,3)
c            call nonsphere(x0(1),x0(2),x0(3),nrad,nang)
            STOP
           end if 
          end if
           end if 
c
c         write(*,*)'x0', x0(1),x0(2),x0(3)
c         write(*,*)'g0', g0(1),g0(2),g0(3)

           if (ipath.eq.0)then 
          call rk4(x0,x1,stepe,g0,gmod)
           else if (ipath.eq.1)then 
          call rk2(x0,x1,stepe,g0,gmod)
           else if (ipath.eq.2)then 
          call euler(x0,x1,stepe,g0,gmod)
           else if (ipath.eq.3)then 
          call lnse(npart,x0,x1,stepe,g0,gmod)
           end if 

         call gpoints(x1(1),x1(2),x1(3),
     $     gxchp,gychp,gzchp,aochp)
         call densvalue(aochp,xdens)
c         write(*,*)'xdens0', xdens0
c         write(*,*)'xdens1', xdens
         call gradvalue(aochp,gxchp,gychp,gzchp,gx,gy,gz)
         g0(1)=gx
         g0(2)=gy
         g0(3)=gz
         gmod=0.0d0
         do l=1,3
          gmod=gmod+g0(l)**2.0d0
         end do 
         gmod=dsqrt(gmod)
           if(xdens.gt.xdens0)then 
          do i=1,3
           x0(i)=x1(i)
          end do
c           write(*,*)ipoint, 'dens', xdens, xdens0
c           write(*,*)ipoint, 'gmod', gx,gy,gz
          xdens0=xdens
          stepe=step0
         else
           stepe=stepe-stepe*0.25
c           write(*,*)'stepe', stepe
c           if (stepe.lt.0.00001d0) then 
c           ibas=1
c           write(*,*)'WARNING, gridpoint', ipoint, 'deleted', xdens, xdens0
c           ibaspoint(ipoint)=0
c           end if
         end if
c         write(*,*)'stepe',stepe   
c         write(*,*)'x1', x0(1),x0(2),x0(3)
         icont=0
         do i=1,nat
          call dist(x1,i,xdist)
          datdist2(i)=xdist
          if (xdist.le.atsphradi(i))then
           ibaspoint(ipoint)=i
           ibas=1
          end if
          if (datdist2(i).ge.datdist(i)) icont=icont+1
          datdist(i)=datdist2(i)
         end do
         if (icont.eq.nat)then 
          write(*,'(a20,i5,a8,E11.4,a2,E11.4,a2,E11.4)')
     $'WARNING,gridpoint',ipoint,' deleted',xinipoint,'  ',xlap(ipoint),
     $' ',gmod
          ibas=1
          ibaspoint(ipoint)=0
         end if 
        end do
       end if  
      end if 
      end do            
        
      end
      
      subroutine vmprod(n,W,V,wv)
      implicit double precision (a-h,o-z)
      dimension w(n,n),v(n),wv(n)
      
      do i=1,n 
        wv(i)=0.0d0
        do k=1,n
         wv(i)=wv(i)+w(i,k)*v(k)
        end do  
      end do 
      end  
      
      subroutine dist(x,iiat,xdist)
      implicit double precision (a-h,o-z)
        include 'parameter.h'
      dimension x(3)
        common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      
c      xdist=dsqrt((x(1)-coord(1,iiat))**2.0d0+(x(2)-coord(2,iiat))**2.0d0+(x(
c     $   3)-coord(3,iiat))**2.0d0)
        xdist=0.0d0
       do i=1,3
      xdist=xdist+(x(i)-coord(i,iiat))*(x(i)-coord(i,iiat))
        end do
        xdist=sqrt(xdist)
      
      end 

      subroutine distnna(x,iiat,xdist)
      implicit double precision (a-h,o-z)
        include 'parameter.h'
      dimension x(3)
        common /coordnna/ coordnna(3,maxnna)
      
c      xdist=dsqrt((x(1)-coordnna(1,iiat))**2.0d0+(x(2)-coordnna(2,iiat))**
c     $  2.0d0+(x(3)-coordnna(3,iiat))**2.0d0)
        xdist=0.0d0
       do i=1,3
      xdist=xdist+(x(i)-coordnna(i,iiat))*(x(i)-coordnna(i,iiat))
        end do
        xdist=sqrt(xdist)
      
      end 

      function hxxfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /stv/ sp(nmax,nmax),tt(nmax,nmax)
      common /nat/ nat,igr,ifg,idum(4)
      common /data/expp(maxp),n(maxp),l(maxp),m(maxp),iat(maxp), ifill(nmax),ifiul(nmax)
      common /coeff/coeff(maxp,maxc)
      common /lim/ llim(nmax),iulim(nmax),ihold(nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/actual/iact,jact,icenter
      
      iactat=ihold(iact)
      f=0.d0
      xx=x-coord(1,iactat)   
      yy=y-coord(2,iactat)   
      zz=z-coord(3,iactat)   
      rr=dsqrt(xx**2+yy**2+zz**2)
      ngri=ifiul(iact)-ifill(iact)+1
      do ki=1,ngri
      ipr=ifill(iact)+ki-1
      nn=n(ipr)
      ll=l(ipr)
      mm=m(ipr)
c      if (abs(xx).lt.1.E-10.and.nn.eq.0)then 
      if (nn.eq.0.or.nn.eq.1)then 
       xprod=0.0d0
      else
       xprod=nn*(nn-1.0d0)*xx**(nn-2.0d0)
      end if
      f=f+(xprod-(4.0d0*expp(ipr)*nn+2.0d0*expp(ipr))*xx**(nn)+4.0d0*
     $expp(ipr)**2*xx**(nn+2.0d0))*(yy**ll)*(zz**mm)*
     $dexp(-expp(ipr)*(rr**2))*coeff(iact,ki)
c      print *,iact,ki,ipr,expp(ipr),coeff(iact,ki)
      enddo

      hxxfunct=f
      
      return
      end
      
      function hyyfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /stv/ sp(nmax,nmax),tt(nmax,nmax)
      common /nat/ nat,igr,ifg,idum(4)
      common /data/expp(maxp),n(maxp),l(maxp),m(maxp),iat(maxp),ifill(nmax),ifiul(nmax)
      common /coeff/coeff(maxp,maxc)
      common /lim/ llim(nmax),iulim(nmax),ihold(nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/actual/iact,jact,icenter
      
      iactat=ihold(iact)
      f=0.d0
      xx=x-coord(1,iactat)   
      yy=y-coord(2,iactat)   
      zz=z-coord(3,iactat)   
      rr=dsqrt(xx**2+yy**2+zz**2)
      ngri=ifiul(iact)-ifill(iact)+1
      do ki=1,ngri
      ipr=ifill(iact)+ki-1
      nn=n(ipr)
      ll=l(ipr)
      mm=m(ipr)
c      if (abs(xx).lt.1.E-10.and.nn.eq.0)then 
      if (ll.eq.0.or.ll.eq.1)then 
       xprod=0.0d0
      else
       xprod=ll*(ll-1.0d0)*yy**(ll-2.0d0)
      end if
      f=f+(xprod-(4.0d0*expp(ipr)*ll+2.0d0*expp(ipr))*yy**(ll)+4.0d0*
     $expp(ipr)**2.0d0*yy**(ll+2.0d0))*(xx**nn)*(zz**mm)*
     $dexp(-expp(ipr)*(rr**2))*coeff(iact,ki)
c      print *,iact,ki,ipr,expp(ipr),coeff(iact,ki)
      enddo

      hyyfunct=f
      
      return
      end


      function hzzfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /stv/ sp(nmax,nmax),tt(nmax,nmax)
      common /nat/ nat,igr,ifg,idum(4)
      common /data/expp(maxp),n(maxp),l(maxp),m(maxp),iat(maxp),ifill(nmax),ifiul(nmax)
      common /coeff/coeff(maxp,maxc)
      common /lim/ llim(nmax),iulim(nmax),ihold(nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/actual/iact,jact,icenter
      
      iactat=ihold(iact)
      f=0.d0
      xx=x-coord(1,iactat)   
      yy=y-coord(2,iactat)   
      zz=z-coord(3,iactat)   
      rr=dsqrt(xx**2+yy**2+zz**2)
      ngri=ifiul(iact)-ifill(iact)+1
      do ki=1,ngri
      ipr=ifill(iact)+ki-1
      nn=n(ipr)
      ll=l(ipr)
      mm=m(ipr)
c      if (abs(xx).lt.1.E-10.and.nn.eq.0)then 
      if (mm.eq.0.or.mm.eq.1)then 
       xprod=0.0d0
      else
       xprod=mm*(mm-1.0d0)*zz**(mm-2.0d0)
      end if
      f=f+(xprod-(4.0d0*expp(ipr)*mm+2.0d0*expp(ipr))*zz**(mm)+4.0d0*
     $expp(ipr)**2.0d0*zz**(mm+2.0d0))*(yy**ll)*(xx**nn)*
     $dexp(-expp(ipr)*(rr**2))*coeff(iact,ki)
c      print *,iact,ki,ipr,expp(ipr),coeff(iact,ki)
      enddo

      hzzfunct=f
      
      return
      end

      subroutine hpoints(x,y,z,hxxchp,hyychp,hzzchp,
     $hxychp,hxzchp,hzychp)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common/actual/iact,jat,icenter
      common /nat/ nat,igr,ifg,idum(4)
      dimension hxxchp(igr),hyychp(igr),hzzchp(igr)
      dimension hxychp(igr),hxzchp(igr),hzychp(igr)

      do i=1,igr
       iact=i
       hxxchp(iact)=hxxfunct(x,y,z)
       hyychp(iact)=hyyfunct(x,y,z)
       hzzchp(iact)=hzzfunct(x,y,z)
       hxychp(iact)=hxyfunct(x,y,z)
       hxzchp(iact)=hxzfunct(x,y,z)
       hzychp(iact)=hzyfunct(x,y,z)
      end do
      return
      end

      subroutine lappoints(x,y,z,hxxchp,hyychp,hzzchp)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common/actual/iact,jat,icenter
      common /nat/ nat,igr,ifg,idum(4)
      dimension hxxchp(igr),hyychp(igr),hzzchp(igr)

      do i=1,igr
       iact=i
       hxxchp(iact)=hxxfunct(x,y,z)
       hyychp(iact)=hyyfunct(x,y,z)
       hzzchp(iact)=hzzfunct(x,y,z)
      end do
      return
      end


      function hxyfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /stv/ sp(nmax,nmax),tt(nmax,nmax)
      common /nat/ nat,igr,ifg,idum(4)
      common /data/expp(maxp),n(maxp),l(maxp),m(maxp),iat(maxp),ifill(nmax),ifiul(nmax)
      common /coeff/coeff(maxp,maxc)
      common /lim/ llim(nmax),iulim(nmax),ihold(nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/actual/iact,jact,icenter
      
      iactat=ihold(iact)
      f=0.d0
      xx=x-coord(1,iactat)   
      yy=y-coord(2,iactat)   
      zz=z-coord(3,iactat)   
      rr=dsqrt(xx**2+yy**2+zz**2)
      ngri=ifiul(iact)-ifill(iact)+1
      do ki=1,ngri
      ipr=ifill(iact)+ki-1
      nn=n(ipr)
      ll=l(ipr)
      mm=m(ipr)
c      if (abs(xx).lt.1.E-10.and.nn.eq.0)then 
      if (nn.eq.0)then 
       xprod=0.0d0
      else
       xprod=nn*xx**(nn-1.0d0)
      end if
      if (ll.eq.0)then 
       yprod=0.0d0
      else
       yprod=ll*yy**(ll-1.0d0)
      end if
      f=f+(yprod-2*expp(ipr)*yy**(ll+1.d0))*(xprod-2*expp(ipr)*xx
     $**(nn+1))*(zz**mm)*dexp(-expp(ipr)*(rr**2))*coeff(iact,ki)
c      print *,iact,ki,ipr,expp(ipr),coeff(iact,ki)
      enddo

      hxyfunct=f
      
      return
      end

      function hxzfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /stv/ sp(nmax,nmax),tt(nmax,nmax)
      common /nat/ nat,igr,ifg,idum(4)
      common /data/expp(maxp),n(maxp),l(maxp),m(maxp),iat(maxp),ifill(nmax),ifiul(nmax)
      common /coeff/coeff(maxp,maxc)
      common /lim/ llim(nmax),iulim(nmax),ihold(nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/actual/iact,jact,icenter
      
      iactat=ihold(iact)
      f=0.d0
      xx=x-coord(1,iactat)   
      yy=y-coord(2,iactat)   
      zz=z-coord(3,iactat)   
      rr=dsqrt(xx**2+yy**2+zz**2)
      ngri=ifiul(iact)-ifill(iact)+1
      do ki=1,ngri
      ipr=ifill(iact)+ki-1
      nn=n(ipr)
      ll=l(ipr)
      mm=m(ipr)
c      if (abs(xx).lt.1.E-10.and.nn.eq.0)then 
      if (nn.eq.0)then 
       xprod=0.0d0
      else
       xprod=nn*xx**(nn-1.0d0)
      end if
      if (mm.eq.0)then 
       yprod=0.0d0
      else
       yprod=mm*zz**(mm-1.0d0)
      end if
      f=f+(yprod-2*expp(ipr)*zz**(mm+1.d0))*(xprod-2*expp(ipr)*xx
     $**(nn+1))*(yy**ll)*dexp(-expp(ipr)*(rr**2))*coeff(iact,ki)
c      print *,iact,ki,ipr,expp(ipr),coeff(iact,ki)
      enddo

      hxzfunct=f
      
      return
      end

      function hzyfunct(x,y,z)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /stv/ sp(nmax,nmax),tt(nmax,nmax)
      common /nat/ nat,igr,ifg,idum(4)
      common /data/expp(maxp),n(maxp),l(maxp),m(maxp),iat(maxp),ifill(nmax),ifiul(nmax)
      common /coeff/coeff(maxp,maxc)
      common /lim/ llim(nmax),iulim(nmax),ihold(nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/actual/iact,jact,icenter
      
      iactat=ihold(iact)
      f=0.d0
      xx=x-coord(1,iactat)   
      yy=y-coord(2,iactat)   
      zz=z-coord(3,iactat)   
      rr=dsqrt(xx**2+yy**2+zz**2)
      ngri=ifiul(iact)-ifill(iact)+1
      do ki=1,ngri
      ipr=ifill(iact)+ki-1
      nn=n(ipr)
      ll=l(ipr)
      mm=m(ipr)
c      if (abs(xx).lt.1.E-10.and.nn.eq.0)then 
      if (mm.eq.0)then 
       xprod=0.0d0
      else
       xprod=mm*zz**(mm-1.0d0)
      end if
      if (ll.eq.0)then 
       yprod=0.0d0
      else
       yprod=ll*yy**(ll-1.0d0)
      end if
      f=f+(yprod-2*expp(ipr)*yy**(ll+1.d0))*(xprod-2*expp(ipr)*zz
     $**(mm+1))*(xx**nn)*dexp(-expp(ipr)*(rr**2))*coeff(iact,ki)
c      print *,iact,ki,ipr,expp(ipr),coeff(iact,ki)
      enddo

      hzyfunct=f
      
      return
      end

        subroutine hessvalue(aochp,gxchp,gychp,gzchp,hxxchp,hyychp,
     $hzzchp,hxychp,hxzchp,hzychp,hess)
      use ao_matrices
      implicit double precision (a-h,o-z)
        include 'parameter.h'
        common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension aochp(igr),gxchp(igr),gychp(igr),gzchp(igr)
        dimension hxxchp(igr),hyychp(igr),hzzchp(igr)
        dimension hxychp(igr),hxzchp(igr),hzychp(igr)
      dimension hess(3,3)
   
        do i=1,3
        do j=1,3
           hess(i,j)=0.0d0
          end do 
        end do 
      do i=1,igr
       do j=1,igr
        hess(1,1)=hess(1,1)+P(i,j)*(2.d0*hxxchp(i)*aochp(j)+2.d0*
     $gxchp(i)*gxchp(j))
        hess(2,2)=hess(2,2)+P(i,j)*(2.d0*hyychp(i)*aochp(j)+2.d0*
     $gychp(i)*gychp(j))
        hess(3,3)=hess(3,3)+P(i,j)*(2.d0*hzzchp(i)*aochp(j)+2.d0*
     $gzchp(i)*gzchp(j))
        hess(2,1)=hess(2,1)+P(i,j)*(2.d0*hxychp(i)*aochp(j)+2.d0*
     $gxchp(i)*gychp(j))
        hess(3,1)=hess(3,1)+P(i,j)*(2.d0*hxzchp(i)*aochp(j)+2.d0*
     $gxchp(i)*gzchp(j))
        hess(3,2)=hess(3,2)+P(i,j)*(2.d0*hzychp(i)*aochp(j)+2.d0*
     $gzchp(i)*gychp(j))
      hess(1,2)=hess(2,1)
      hess(1,3)=hess(3,1)
      hess(2,3)=hess(3,2)
       end do 
      end do

      end

        subroutine laphessvalue(aochp,gxchp,gychp,gzchp,hxxchp,hyychp,
     $hzzchp,hess)
      use ao_matrices
      implicit double precision (a-h,o-z)
        include 'parameter.h'
        common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension aochp(igr),gxchp(igr),gychp(igr),gzchp(igr)
        dimension hxxchp(igr),hyychp(igr),hzzchp(igr)
      dimension hess(3,3)
   
        do i=1,3
           hess(i,i)=0.0d0
        end do 
      do i=1,igr
       do j=1,igr
        hess(1,1)=hess(1,1)+P(i,j)*(2.d0*hxxchp(i)*aochp(j)+2.d0*
     $gxchp(i)*gxchp(j))
        hess(2,2)=hess(2,2)+P(i,j)*(2.d0*hyychp(i)*aochp(j)+2.d0*
     $gychp(i)*gychp(j))
        hess(3,3)=hess(3,3)+P(i,j)*(2.d0*hzzchp(i)*aochp(j)+2.d0*
     $gzchp(i)*gzchp(j))
       end do 
      end do

      end


      subroutine lnse(n,x0,xx,step,zsum,totzsum)
      implicit double precision(a-h,o-z)
        common /iops/iopt(100)
        common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
        include 'parameter.h'
      dimension gxchp(:),gychp(:),gzchp(:),aochp(:)
      allocatable gxchp,gychp,gzchp,aochp
      dimension x0(3),x1(3),xx(3)
             dimension zsum(3) 
        allocate(gxchp(igr))
        allocate(gychp(igr))
        allocate(gzchp(igr))
        allocate(aochp(igr))

      nf=1
      call gpoints(x0(1),x0(2),x0(3),gxchp,gychp,gzchp,aochp)
      call densvalue(aochp,xdensmax)
      xn=float(n)
      do j=1,n
      xj=float(j)
      do i=1,3
       x1(i)=x0(i)+xj*step*zsum(i)/totzsum/xn
      end do
c      write(*,*)'xj/xn',xj*step/xn
      call gpoints(x1(1),x1(2),x1(3),gxchp,gychp,gzchp,aochp)
      call densvalue(aochp,xdens1)
c      write(*,*)'dens',xdens1
      if(xdens1.gt.xdensmax)then
       xdensmax=xdens1
       nf=j
       do k=1,3
        xx(k)=x1(k)
       end do 
      end if 
      end do 
c       write(*,*)'n',nf
      end 
      
      subroutine euler(x0,xx,step,zsum,totzsum)
      implicit double precision(a-h,o-z)
        common /iops/iopt(100)
        common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
        include 'parameter.h'
      dimension gxchp(:),gychp(:),gzchp(:),aochp(:)
      allocatable gxchp,gychp,gzchp,aochp
      dimension x0(3),xx(3)
             dimension zsum(3) 
        allocate(gxchp(igr))
        allocate(gychp(igr))
        allocate(gzchp(igr))
        allocate(aochp(igr))

      do i=1,3
       xx(i)=x0(i)+step*zsum(i)/totzsum
      end do
      do i=1,3
      write(*,*)'disp',step*zsum(i)/totzsum
      end do 
      end 
      
      subroutine rk2(x0,xx,step,zsum,totzsum)
      implicit double precision(a-h,o-z)
        common /iops/iopt(100)
        common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
        include 'parameter.h'
      dimension gxchp(:),gychp(:),gzchp(:),aochp(:)
      dimension xk1(3),xk2(3)
      allocatable gxchp,gychp,gzchp,aochp
      dimension x0(3),xx(3),g(3)
             dimension zsum(3) 
        allocate(gxchp(igr))
        allocate(gychp(igr))
        allocate(gzchp(igr))
        allocate(aochp(igr))

      do i=1,3
       xk1(i)=step*zsum(i)/totzsum
      end do
      
      call gpoints((x0(1)+xk1(1)/2.d0),(x0(2)+xk1(2)/2.d0),(x0(3)+
     $  xk1(3)/2.d0),gxchp,gychp,gzchp,aochp)
      call gradvalue(aochp,gxchp,gychp,gzchp,gx,gy,gz)
      g(1)=gx
      g(2)=gy
      g(3)=gz
      
      gzsum=0.0d0
      do l=1,3
       gzsum=gzsum+g(l)**2.0d0
      end do 
      gzsum=dsqrt(gzsum)
      
      do i=1,3
       xk2(i)=step*g(i)/gzsum
       xx(i)=x0(i)+xk2(i)
      end do 
      
      
      end 

      
      subroutine rk4(x0,xx,step,zsum,totzsum)
      implicit double precision(a-h,o-z)
        common /iops/iopt(100)
        common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
        include 'parameter.h'
      dimension gxchp(:),gychp(:),gzchp(:),aochp(:)
      dimension xk1(3),xk2(3),xk3(3),xk4(3)
      allocatable gxchp,gychp,gzchp,aochp
      dimension x0(3),xx(3),g(3)
             dimension zsum(3) 
        allocate(gxchp(igr))
        allocate(gychp(igr))
        allocate(gzchp(igr))
        allocate(aochp(igr))

      do i=1,3
       xk1(i)=step*zsum(i)/totzsum
      end do
      
      call gpoints((x0(1)+xk1(1)/2.d0),(x0(2)+xk1(2)/2.d0),(x0(3)+
     $  xk1(3)/2.d0),gxchp,gychp,gzchp,aochp)
      call gradvalue(aochp,gxchp,gychp,gzchp,gx,gy,gz)
      g(1)=gx
      g(2)=gy
      g(3)=gz
      
      gzsum=0.0d0
      do l=1,3
       gzsum=gzsum+g(l)**2.0d0
      end do 
      gzsum=dsqrt(gzsum)
      
      do i=1,3
       xk2(i)=step*g(i)/gzsum
      end do 

      call gpoints((x0(1)+xk2(1)/2.d0),(x0(2)+xk2(2)/2.d0),(x0(3)+
     $  xk2(3)/2.d0),gxchp,gychp,gzchp,aochp)
      call gradvalue(aochp,gxchp,gychp,gzchp,gx,gy,gz)
      g(1)=gx
      g(2)=gy
      g(3)=gz
      
      gzsum=0.0d0
      do l=1,3
       gzsum=gzsum+g(l)**2.0d0
      end do 
      gzsum=dsqrt(gzsum)
      
      do i=1,3
       xk3(i)=step*g(i)/gzsum
      end do 
      
      call gpoints((x0(1)+xk3(1)),(x0(2)+xk3(2)),(x0(3)+
     $  xk3(3)),gxchp,gychp,gzchp,aochp)
      call gradvalue(aochp,gxchp,gychp,gzchp,gx,gy,gz)
      g(1)=gx
      g(2)=gy
      g(3)=gz
      
      gzsum=0.0d0
      do l=1,3
       gzsum=gzsum+g(l)**2.0d0
      end do 
      gzsum=dsqrt(gzsum)
      
      do i=1,3
       xk4(i)=step*g(i)/gzsum
      end do 
      
      do l=1,3
       xx(l)=x0(l)+xk1(l)/6.d0+xk2(l)/3.d0+xk3(l)/3.d0+xk4(l)/6.d0
      end do 
      
      do i=1,3
c      write(*,*)'disp',(xk1(i)/6.d0+xk2(i)/3.d0+xk3(i)/3.d0+xk4(i)/6.d0)
c
      end do 
      
      end 

      subroutine nonsphere(cx,cy,cz,nrad,nang)
      use ao_matrices
      implicit real*8 (a-h,o-z)
        include 'parameter.h'
        common /quadrat/th(1000),ph(1000),w(1000),wr(500),Xr(500)
        common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
        common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
        common /iops/iopt(100)
        common /nna/ xnonradi(maxnna),atsphradi(maxat),nna
        common/pha2/pha,phb,rr00
        dimension gxchp(:),gychp(:),gzchp(:),aochp(:)
        allocatable aochp,gxchp,gychp,gzchp
      
        allocate(gxchp(igr))
        allocate(gychp(igr))
        allocate(gzchp(igr))
        allocate(aochp(igr))
c pi=dacos(-1.0d0)
      xangle=cos(45.0d0*pi/180.0d0)
      tres=0.05d0 
      xnonradi(nna)=0.0d0
       do irad=nrad,1,-1
        if(xr(irad).gt.tres)then
         do iang=1,nang
          thx=th(iang)
          fix=ph(iang)
          rr=xr(irad)
      xabs0=rr*dsin(thx)*dcos(fix)
      yabs0=rr*dsin(thx)*dsin(fix)
      zabs0=rr*dcos(thx)
      xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
      yabs=-dcos(pha)*dsin(phb)*xabs0+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
      zabs=dsin(pha)*dsin(phb)*xabs0-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0
                x=xabs+cx
                y=yabs+cy
                z=zabs+cz
                call gpoints(x,y,z,gxchp,gychp,gzchp,aochp)
          call gradvalue(aochp,gxchp,gychp,gzchp,gx,gy,gz)
          if(angle(cx-x,cy-y,cz-z,gx,gy,gz
     $).lt.xangle) go to 44
c          write(*,*) xr(irad),angle(coord(1,i)-x,coord(2,i)-y,coord(3,i)-z,gx
c     $,gy,gz)
         end do  
        end if
       end do   
44       continue
       xnonradi(nna)=xr(irad+1)
c       write(*,*)'resatom',atsphradi(i)
      
      write(*,*) 'nonnuclear atractor sphere radi : ',xnonradi(nna)
      
      
      deallocate(gxchp)
      deallocate(gychp)
      deallocate(gzchp)
      deallocate(aochp)
      
      end 
      



      
      subroutine laplacian(ipoints,chp,nrad,nang,xlap)
      use ao_matrices
      implicit real*8 (A-H,O-Z)
      include 'parameter.h'
      common/actual/iact,jat,icenter
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/ iopt(100)
      dimension chp(ipoints,igr)
      dimension xlap(nat*nrad*nang)

      dimension chpder(:,:),chp2der(:,:)
      allocatable chpder, chp2der

      ipoints=nrad*nang*nat
      allocate(chpder(ipoints,igr),chp2der(ipoints,igr))

      do i=1,ipoints
       xlap(i)=0.0d0
      end do


      do ixyz=1,3
       do ii=1,igr
        irun=1
        iact=ii
        do jat=1,nat
         icenter=jat
         do k=1,Nrad
          do i=1,Nang
           chpder(irun,iact)=drho(i,k,ixyz)
           if(ixyz.eq.1)chp2der(irun,iact)=-dfunct(i,k)
           irun=irun+1
          enddo
         enddo
        enddo
       end do


       do l=1,ipoints
c       do l=1100,1100
        x=0.0d0
        do i=1,igr
         do j=1,igr
          x=x+2.0d0*P(i,j)*chpder(l,j)*chpder(l,i)
c          x=x+2.0d0*P(i,j)*chpder(l,j)*chp(l,i)
          if(ixyz.eq.1)x=x+2.0d0*P(i,j)*chp2der(l,j)*chp(l,i)
         end do
        end do
        xlap(l)=xlap(l)+x
       end do

      end do

c      deallocate(chp)

c      write(*,*)'LAPLACIAN TEST'
      deallocate(chpder,chp2der)

      end
      function dfunct(i,k)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,idum(4)
      common /data/expp(maxp),n(maxp),l(maxp),m(maxp),iat(maxp),ifill(nmax),ifiul(nmax)
      common /coeff/coeff(maxp,maxc)
      common /lim/ llim(nmax),iulim(nmax),ihold(nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /quadrat/th(1000),ph(1000),w(1000),wr(500),Xr(500)
      common/actual/iact,jact,icenter
      common/pha2/pha,phb,rr00


      iactat=ihold(iact)
      f=0.d0
      thx=th(i)
      fix=ph(i)
      rr=xr(k)
      xabs0=rr*dsin(thx)*dcos(fix)
      yabs0=rr*dsin(thx)*dsin(fix)
      zabs0=rr*dcos(thx)
      xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
      yabs=-dcos(pha)*dsin(phb)*xabs0+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
      zabs=dsin(pha)*dsin(phb)*xabs0-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0
      x=xabs-coord(1,iactat)   +coord(1,icenter)
      y=yabs-coord(2,iactat)   +coord(2,icenter)
      z=zabs-coord(3,iactat)   +coord(3,icenter)
      rr=dsqrt(x**2+y**2+z**2)
      ngri=ifiul(iact)-ifill(iact)+1
      do ki=1,ngri
       ipr=ifill(iact)+ki-1
       nn=n(ipr)
       ll=l(ipr)
       mm=m(ipr)
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

       dd=(dx2+dy2+dz2)*dexp(-expp(ipr)*(rr**2))
       f=f+dd*coeff(iact,ki)
      enddo

      dfunct=-f

      return
      end
      function drho(i,k,ixyz)
c angular point i, radial point k
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /data/expp(maxp),n(maxp),l(maxp),m(maxp),iat(maxp),ifill(nmax),ifiul(nmax)
      common /coeff/coeff(maxp,maxc)
      common /lim/ llim(nmax),iulim(nmax),ihold(nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /quadrat/th(1000),ph(1000),w(1000),wr(500),Xr(500)
      common/actual/iact,jact,icenter
      
      iactat=ihold(iact)
      fx=0.d0
      thx=th(i)
      fix=ph(i)
      rr=xr(k)
      xabs=rr*dsin(thx)*dcos(fix)+coord(1,icenter)
      yabs=rr*dsin(thx)*dsin(fix)+coord(2,icenter)
      zabs=rr*dcos(thx)+coord(3,icenter)
      x=xabs-coord(1,iactat)   
      y=yabs-coord(2,iactat)   
      z=zabs-coord(3,iactat)   
      rr=dsqrt(x**2+y**2+z**2)
      ngri=ifiul(iact)-ifill(iact)+1
      do ki=1,ngri
       ipr=ifill(iact)+ki-1
       nn=n(ipr)
       ll=l(ipr)
       mm=m(ipr)
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
       fx=fx+dx*coeff(iact,ki)
      enddo
       drho=fx
      return
      end

      subroutine cubeqtaim(iatinp,xin,yin,zin,xbas)
      implicit double precision (a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /nna/ xnonradi(maxnna),atsphradi(maxat),nna
      logical icontrol
      dimension x0(3),x1(3),g0(3),H0(3,3),hh0(3,3)

      allocatable gxchp(:),gychp(:),gzchp(:),aochp(:)
      allocatable hxxchp(:),hyychp(:),hzzchp(:)
      allocatable hxychp(:),hxzchp(:),hzychp(:)

      allocate(gxchp(igr))
      allocate(gychp(igr))
      allocate(gzchp(igr))
      allocate(aochp(igr))
      allocate(hxxchp(igr))
      allocate(hyychp(igr))
      allocate(hzzchp(igr))
      allocate(hxychp(igr))
      allocate(hxzchp(igr))
      allocate(hzychp(igr))


      icontrol=.TRUE.
c nonnuclear atractor control
      nna=0
c  step
      stepe0=0.3d0
      stepe=stepe0
c atom control
      xbas=0.0d0
c            xbas=-1000.0d0
c coordenates
      x0(1)=xin 
      x0(2)=yin
      x0(3)=zin
c maxdist
      ymaxdist=28.d0
      do i=1,nat
       call dist(x0,i,xdist)
       if (xdist.le.atsphradi(i))then 
        if(i.eq.iatinp)then 
         xbas=1.0d0
         icontrol=.FALSE.
        else
         xbas=0.0d0
c            xbas=-1000.0d0
         icontrol=.FALSE.
        end if
       end if
      end do 
      
      if (xdist.gt.ymaxdist.and.xbas.eq.0) then
       xbas=0.d0
c            xbas=-1000.0d0
       icontrol=.FALSE.
      else  
       call gpoints(x0(1),x0(2),x0(3),
     $      gxchp,gychp,gzchp,aochp)
       call densvalue(aochp,xdens0)
       call gradvalue(aochp,gxchp,gychp,gzchp,gx,gy,gz)
       gmod=dsqrt(gx**2.0d0+gy**2.0d0+gz**2.0d0)
       if(xdens0.lt.1E-15) then 
        xbas=0.0d0  
c            xbas=-1000.0d0
        icontrol=.FALSE.
       else 
        g0(1)=gx
        g0(2)=gy
        g0(3)=gz
        do while(icontrol)
c NONNUCLEAR ATRACTORS   
         if(xdens0.gt.0.001d0.and.gmod.le.1e-8)then  
          call hpoints(x0(1),x0(2),x0(3),hxxchp,hyychp,hzzchp,hxychp,
     $         hxzchp,hzychp)
          call hessvalue(aochp,gxchp,gychp,gzchp,hxxchp,hyychp,
     $                   hzzchp,hxychp,hxzchp,hzychp,h0)
          call diagonalize(3,3,h0,hh0,0)
          if(h0(1,1).lt.0.0d0.and.h0(2,2).lt.0.0d0.and.h0(3,3)
     $                                          .lt.0.0d0)then
           nna=nna+1
           write(*,*)'NONNUCLEARATRACTOR',nna,' position: ',x0(1),x0(2)
     $                                                           ,x0(3)
           write(*,*)'hessian eigenvalues',h0(1,1),h0(2,2),h0(3,3)
           STOP
          end if 
         end if
c step         
         call rk4(x0,x1,stepe,g0,gmod)
         call gpoints(x1(1),x1(2),x1(3),
     $            gxchp,gychp,gzchp,aochp)
         call densvalue(aochp,xdens)
         call gradvalue(aochp,gxchp,gychp,gzchp,gx,gy,gz)
         g0(1)=gx
         g0(2)=gy
         g0(3)=gz
         gmod=0.0d0
         do l=1,3
          gmod=gmod+g0(l)**2.0d0
         end do 
         gmod=dsqrt(gmod)
         if(xdens.gt.xdens0)then 
          do i=1,3
           x0(i)=x1(i)
          end do
          xdens0=xdens
          stepe=stepe0
         else
          stepe=stepe-stepe*0.25
         end if
         do i=1,nat
          call dist(x1,i,xdist2)
          if(i.eq.iatinp)then
           if (xdist2.le.atsphradi(i))then
            xbas=1.0d0
            icontrol=.FALSE.
           else if(xdist2.gt.1.2d0*xdist)then  
            xbas=0.0d0
c            xbas=-1000.0d0
            icontrol=.FALSE.
           end if 
          else 
           if (xdist2.le.atsphradi(i))then
c            xbas=-1000.0d0
            xbas=0.0d0
            icontrol=.FALSE.
           end if   
          end if
         end do 
        end do 
       end if  
      end if  
      deallocate(gxchp)
      deallocate(gychp)
      deallocate(gzchp)
      deallocate(aochp)
      deallocate(hxxchp)
      deallocate(hyychp)
      deallocate(hzzchp)
      deallocate(hxychp)
      deallocate(hxzchp)
      deallocate(hzychp)
      end

      subroutine atsphere()
      use basis_set
      implicit real*8 (a-h,o-z)
      include 'parameter.h'
      common /quadrat/th(1000),ph(1000),w(1000),wr(500),Xr(500)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /nna/ xnonradi(maxnna),atsphradi(maxat),nna
      common /gridp/Nrad,Nang
      common/pha2/pha,phb,rr00

       xangle=cos(45.0d0*pi/180.0d0)
       tres=0.15d0 
       do i=1,nat
        atsphradi(i)=0.0d0
        do irad=nrad,1,-1
         if(xr(irad).gt.tres)then
          do iang=1,nang
           thx=th(iang)
           fix=ph(iang)
           rr=xr(irad)
           xabs0=rr*dsin(thx)*dcos(fix)
           yabs0=rr*dsin(thx)*dsin(fix)
           zabs0=rr*dcos(thx)
           xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
           yabs=-dcos(pha)*dsin(phb)*xabs0+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
           zabs=dsin(pha)*dsin(phb)*xabs0-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0
           x=xabs+coord(1,i)
           y=yabs+coord(2,i)
           z=zabs+coord(3,i)
           gx=dfunctxyz(1,x,y,z)
           gy=dfunctxyz(2,x,y,z)
           gz=dfunctxyz(3,x,y,z)
           if(angle(-xabs,-yabs,-zabs,gx,gy,gz).lt.xangle) go to 44
          end do  
         end if
        end do   
44     continue
       atsphradi(i)=xr(irad+1)
      end do 
      end 

       function angle(x,y,z,xx,yy,zz)
       implicit none
       double precision x,y,z,xx,yy,zz
       double precision angle

       angle=((x*xx+y*yy+z*zz)/(dsqrt((x)**2.+(y)**2.+(z)**2.)*dsqrt((xx)**2.+(yy)**2.+(zz)**2.)))

       end function angle
