!! MG: Pedro, aquesta subrutina ja la pots borrar !!
!! MG: A oslo.f esta la meva versio (incloent assignacio electrons) !!
      subroutine old_effao3d_u(itotps,ndim,omp,chp,sat,wp,omp2)
      use basis_set
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      integer,intent(in) :: itotps,ndim
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
c      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /localspin/xlsa(maxat,maxat),ua(maxat)
      common /iops/iopt(100)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat)
      dimension wp(itotps)
      dimension sat(ndim,ndim,nat)
      character*(5) key
      character*80 line
c
      dimension effrho(maxat),effrhoacc(maxat)
      allocatable s0(:,:), sm(:,:), c0(:,:), splus(:,:), pp0(:,:)
      allocatable s0all(:),is0all(:)
      allocatable ::  scr(:),u0(:,:)

      ihirsh = Iopt(6) 
      iallpo= iopt(7) 
      ieffao=  Iopt(12) 
      icube  = Iopt(13) 
      ieffthr  = Iopt(24) 
      iatps=nang*nrad

      xmaxocc=real(ieffthr)/1000.0d0 

      print *,' '
      print *,' --------------------------------------'
      print *,'             DOING EFFAO-3D '
      print *,'         UNPAIRED/PAIRED DENSITY '
      print *,' --------------------------------------'
      print *,' '

      jocc=0
      xmaxotot=0.d0

      allocate(scr(iatps*nat))
      allocate (s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      allocate (c0(ndim,ndim), pp0(ndim,ndim), is0all(ndim))
      allocate (u0(ndim,ndim))

c build U matrix; U = 2P -PSP : u0
c substract from P to get paired density:p0
      do i=1,igr
       do j=1,igr
        pp0(i,j)=p(i,j)
        s0(i,j)=s(i,j)
       end do
      end do
      call to_lowdin_basis(ndim,pp0,s0)
      do i=1,igr
       do j=1,igr
         u0(i,j)=2.0d0*p(i,j)-s0(i,j)
         pp0(i,j)=p(i,j)-u0(i,j)
       end do
      end do
c Calculate number of efectively unpaired electrons
      xxx=0.0d0
      do jcenter=1,nat
        xx=0.0d0 
        do i=1,igr                        
         do k=1,igr
          xx=xx+u0(i,k)*sat(k,i,jcenter)
         end do
        end do
        ua(jcenter)=xx
        xxx=xxx+xx
      end do
      print *,'  '
      print *,' EFFECTIVELY UNPAIRED ELECTRONS'
      print *,'  '
      print *,'    Atom     u_A'
      print *,' -----------------'
      call vprint(ua,nat,maxat,1)
      print *,' ------------------'
      write(*,'(a16,f10.5)') ' Sum check N_D = ' ,xxx 
      line ='   FRAGMENT ANALYSIS : Num. eff. unpaired elec.'
      call group_by_frag_vec(1,line ,ua)

C loop over fragments
      do ifrag=1,icufr

         do ifut=1,iatps*nat                         
          scr(ifut)=0.0d0
          do icenter=1,nfrlist(ifrag)
           scr(ifut)=scr(ifut)+ omp2(ifut,ifrlist(icenter,ifrag))
          end do
         end do

c Computing Becke atomic NET orbital overlap
C ALLPOINTS should be default here
c a distance based screening would be interesting for very large systems
        iallpo0=1
        if(iallpo0.eq.1) then
        do mu=1,ndim
         do nu=1,mu
           x=0.d0
           do jcenter=1,nat
            do ifut=iatps*(jcenter-1)+1,iatps*jcenter
             x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*scr(ifut)*scr(ifut)*omp(ifut)
            end do
           end do
           s0(mu,nu)=x
           s0(nu,mu)=x
         enddo
        enddo
        else
        do mu=1,ndim
         do nu=1,mu
           x=0.d0
           do icenter=1,nfrlist(ifrag)
            jcenter=ifrlist(icenter,ifrag)
            do ifut=iatps*(jcenter-1)+1,iatps*jcenter
             x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*scr(ifut)*scr(ifut)*omp(ifut)
            end do
           end do
           s0(mu,nu)=x
           s0(nu,mu)=x
         enddo
        enddo
        end if

c make SAA1/2 on Splus
       call build_Smp(igr,S0,Sm,Splus,0)

c doing unpaired density
c transform u0, store on s0
       do i=1,igr
        do j=1,igr
         s0(i,j)=u0(i,j)
        end do
       end do
       call to_lowdin_basis(igr,Splus,s0)
c diagonalize to get effao
       call diagonalize(igr,igr,s0,C0,0)
       call to_AO_basis(igr,igr,Sm,C0)

C But max number of effaos in one center is nocc
        i=1
        imaxo=0
        do while(s0(i,i).ge.xmaxocc.and.i.le.nocc) 
         imaxo=i
         i=i+1
        end do
        xmaxo=0.0d0
        do i=1,igr  
         xmaxo=xmaxo+s0(i,i)
        end do
c
        xx0=0.0d0
        xx1=0.0d0
        do icenter=1,nfrlist(ifrag)
         xx1=xx0+qat(ifrlist(icenter,ifrag),1)
         do jcenter=1,nfrlist(ifrag)
          xx0=xx0+op(ifrlist(icenter,ifrag),ifrlist(jcenter,ifrag))
         end do
        end do
        write(*,*) '                  '
        write(*,'(a,i4,a)') ' ** FRAGMENT ',ifrag,' ** '
        write(*,*) '                  '
        write(*,*) '         UNPAIRED DENSITY '
        write(*,*) '                  '

        write(*,'(a29,i3,f10.5)') ' Net occupation for fragment ',ifrag, xmaxo
        write(*,'(a23,f6.4)') ' Net occupation using >',xmaxocc
        write(*,60) (s0(mu,mu),mu=1,imaxo)

c ...calculate gross occupations from orbitals!!!
        xx0=0.0d0
        do i=1,imaxo
         xxx=0.0d0
         do icenter=1,nfrlist(ifrag)
          jcenter=ifrlist(icenter,ifrag)
          xx=0.0d0 
          do j=1,igr                        
           do k=1,igr
            xx=xx+c0(k,i)*sat(k,j,jcenter)*c0(j,i)
           end do
          end do
          xxx=xxx+xx
         end do
         xxx=xxx*s0(i,i)
         s0all(i)=xxx
         xx0=xx0+xxx
        end do
        write(*,*) '                  '
        write(*,'(a31,i3,f10.5)') ' Gross occupation for fragment ',ifrag, xx0
        write(*,60) (s0all(mu),mu=1,imaxo)

        do k=1,imaxo           
         do mu=1,igr                       
           p0(mu,k)=c0(mu,k)
         end do
         p0net(k,ifrag)=s0(k,k)
         p0gro(k,ifrag)=s0all(k)
        end do
        ip0(ifrag)=imaxo

c  write cube file
        if(icube.eq.1) call cubegen4(ifrag,4)
c clean
        p0=0.0d0
        p0net=0.0d0
        p0gro=0.0d0
        ip0=0

c doing paired density
c transform pp0, store on s0
       do i=1,igr
        do j=1,igr
         s0(i,j)=pp0(i,j)
        end do
       end do
       call to_lowdin_basis(igr,Splus,s0)
c diagonalize to get effao
       call diagonalize(igr,igr,s0,C0,0)
       call to_AO_basis(igr,igr,Sm,C0)

C But max number of effaos in one center is nocc
        i=1
        imaxo=0
        do while(s0(i,i).ge.xmaxocc.and.i.le.nocc) 
         imaxo=i
         i=i+1
        end do
        xmaxo=0.0d0
        do i=1,igr  
         xmaxo=xmaxo+s0(i,i)
        end do
c
        xx0=0.0d0
        xx1=0.0d0
        do icenter=1,nfrlist(ifrag)
         xx1=xx0+qat(ifrlist(icenter,ifrag),1)
         do jcenter=1,nfrlist(ifrag)
          xx0=xx0+op(ifrlist(icenter,ifrag),ifrlist(jcenter,ifrag))
         end do
        end do
        write(*,*) '                  '
        write(*,'(a,i4,a)') ' ** FRAGMENT ',ifrag,' ** '
        write(*,*) '                  '
        write(*,*) '         PAIRED DENSITY '
        write(*,*) '                  '

        write(*,'(a29,i3,f10.5)') ' Net occupation for fragment ',ifrag, xmaxo
        write(*,'(a23,f6.4)') ' Net occupation using >',xmaxocc
        write(*,60) (s0(mu,mu),mu=1,imaxo)

c ...calculate gross occupations from orbitals!!!
        xx0=0.0d0
        do i=1,imaxo
         xxx=0.0d0
         do icenter=1,nfrlist(ifrag)
          jcenter=ifrlist(icenter,ifrag)
          xx=0.0d0 
          do j=1,igr                        
           do k=1,igr
            xx=xx+c0(k,i)*sat(k,j,jcenter)*c0(j,i)
           end do
          end do
          xxx=xxx+xx
         end do
         xxx=xxx*s0(i,i)
         s0all(i)=xxx
         xx0=xx0+xxx
        end do
         
        write(*,*) '                  '
        write(*,'(a31,i3,f10.5)') ' Gross occupation for fragment ',ifrag, xx0
        write(*,60) (s0all(mu),mu=1,imaxo)


        do k=1,imaxo           
         do mu=1,igr                       
           p0(mu,k)=c0(mu,k)
         end do
         p0net(k,ifrag)=s0(k,k)
         p0gro(k,ifrag)=s0all(k)
        end do
        ip0(ifrag)=imaxo

c  write cube file
        if(icube.eq.1) call cubegen4(ifrag,3)

c end outer loop over fragments
        end do
60      format(7h OCCUP.   ,8F9.4)
       deallocate (scr,s0, sm, c0,splus,pp0, s0all,is0all)

         end

      subroutine effao_minbas()
      use basis_set
      use ao_matrices, only: p,c
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /actual/iact,jat,icenter
      common /iops/iopt(100)
      common/effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      dimension illim_eff(maxat),iulim_eff(maxat)
      character*(5) key

      allocatable s0(:,:), sm(:,:), c0(:,:), splus(:,:), pp0(:,:)
      allocatable s0all(:),is0all(:)
      
c for svd
      allocatable scr(:,:),ZZ(:,:),v(:,:),op_eff(:,:),s_eff(:,:)
      allocatable u(:,:),eig(:),p_eff(:,:),dum(:,:),a_eff(:,:)

      ndim=igr

      allocate (s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      allocate (c0(ndim,ndim), pp0(ndim,ndim), is0all(ndim))

      icube=iopt(13)
      ieffthr  = Iopt(24)

      xmaxocc=real(ieffthr)/1000.0d0 

      print *,' '
      print *,' ------------------------'
      print *,'  DOING  MULLIKEN EFFAO '
      print *,' ------------------------'
      print *,' '

       jocc=0
       xmaxotot=0.d0
       k1=0

c Making block-diagonal S  and P matrix (AO)
         do i=1,igr
          do j=1,igr
           s0(i,j)=0.0d0
           pp0(i,j)=0.0d0
          end do
         end do
         do icenter=1,nat
         do mu=llim(icenter),iulim(icenter)
          do nu=llim(icenter),iulim(icenter)
           S0(mu,nu)=s(mu,nu)
           pp0(mu,nu)=p(mu,nu)
          end do
         end do
         end do
c make S1/2
       call build_Smp(igr,S0,Sm,Splus,0)

c block P and S  matrix
c tranform blovck P0 with Splus
      call to_lowdin_basis(igr,Splus,pp0)
      call diagonalize(igr,igr,pp0,C0,0)

c reorder colums to keep block diagonality
       k=0
       do icenter=1,nat
        i=k+1
        kmax=iulim(icenter)
        do while (k.lt.kmax) 
         xx=0.0d0
         do mu=llim(icenter),iulim(icenter)
          xx=xx+abs(c0(mu,i))
         end do
         if(xx.gt.1.0d-3) then
          k=k+1
          if(k.ne.i) then
           xx=pp0(k,k)
           pp0(k,k)=pp0(i,i)
           pp0(i,i)=xx
           do j=1,igr
            xx=c0(j,k) 
            c0(j,k)=c0(j,i)
            c0(j,i)=xx
           end do
          end if
         end if
         i=i+1
        end do
        
c sort occup again for each atom
        do ii=llim(icenter),iulim(icenter)
         xx=pp0(ii,ii)
         iix=ii
         do jj=ii+1,iulim(icenter)
          if(pp0(jj,jj).gt.xx) then
           xx=pp0(jj,jj)
           iix=jj
          end if
         end do
         if(iix.ne.ii) then
          xx=pp0(iix,iix)
          pp0(iix,iix)=pp0(ii,ii)
          pp0(ii,ii)=xx
          do j=1,igr
           xx=c0(j,iix) 
           c0(j,iix)=c0(j,ii)
           c0(j,ii)=xx
          end do
         end if
        end do
       end do

C backtrasnform
C Sm is blcok-diagonal so block-diagonality is conserved. 
c eff-aOs expanded in the atom s basis set
      call to_AO_basis(igr,igr,Sm,C0)

CCCCCCC
C loop over atoms
CCCCCCC
        kk=0
        do icenter=1,nat

        xmaxo=0.0d0
        imaxo=0
        i=llim(icenter)
        do while(pp0(i,i).ge.xmaxocc.and.i.lt.iulim(icenter)) 
         xmaxo=xmaxo+pp0(i,i)
         imaxo=i
         i=i+1
        end do

        write(*,*) '                  '
        write(*,'(a,i4,a)') ' ** ATOM ',icenter,' ** '
        write(*,*) '                  '

        write(*,'(a25,i3,f10.5)') ' Net occupation for atom ',icenter, xmaxo
        if(icase.eq.0) then
        write(*,'(a30,f8.4)') ' Deviation from net population ', xmaxo -op(icenter,icenter)
        end if
        write(*,'(a23,f6.4)') ' Net occupation using >',xmaxocc
        write(*,60) (pp0(mu,mu),mu=llim(icenter),imaxo)
60      format(1x,7h OCCUP.   ,8F9.4)

        do i=llim(icenter),imaxo
         jocc=jocc+1
         s0all(jocc)=pp0(i,i)
         is0all(jocc)=icenter
        end do

C can not calculate gross populations. eff-aos already truncated
c
        do k=llim(icenter),imaxo           
         kk=kk+1
         do mu=llim(icenter),iulim(icenter)
           p0(mu,kk)=c0(mu,k)
         end do
         p0net(kk,icenter)=pp0(k,k)
         p0gro(kk,icenter)=pp0(k,k)
        end do
        ip0(icenter)=imaxo-llim(icenter)+1

c OUTPUT ORBITALS FOR VISUALIZATION
        if(icube.eq.1) call cubegen3(icenter,icase)

c end loop over atoms
        end do

        if(icase.eq.1) then
        open(68,file='eff-aos.mull')
        call ival(68,"Number of alpha electrons",k1)
        call ival(68,"Number of basis functions",igr)
        call rarr(68,"Alpha MO occupations",k1,igr,s0all)
        call rmat(68,"Alpha MO coefficients",k1,igr,nmax,p0)          
        else
        call ival(68,"Number of beta electrons",k1)
        call rarr(68,"Beta MO occupations",k1,igr,s0all)
        call rmat(68,"Beta MO coefficients",k1,igr,nmax,p0)          
        close(68)
        end if

       deallocate (s0, sm, c0,splus,pp0, s0all,is0all)

c DO  MINIMAL BASIS
C NEW HILBERT-SPACE OVER EFF-AO BASIS 
C only if eff-AOs have been computed for all atoms
       ndim2=nocc

       ieffdim=0
       do i=1,nat
        ieffdim=ieffdim+ip0(i)
        write(*,*) 'Using ',ip0(i),'orbitals on atom ',i 
       end do
       write(*,*) ' '
       write(*,*) 'Using ',ieffdim,' eff orbitals as minimal basis' 
       write(*,*) ' '
       if(ieffdim.ne.kk) stop 'inconsistency'

       allocate(scr(igr,igr),s_eff(ieffdim,ieffdim))
       allocate(sm(ieffdim,ieffdim),dum(ieffdim,ieffdim))
      allocate(p_eff(ieffdim,ieffdim))
c First off, orthogonalized effaos
c C_eff <- C_eff S_eff^(-1/2)
c Overlap matrix of eff-aos 
       do i=1,ieffdim
        do j=1,ieffdim
         s_eff(i,j)=0.d0
        end do
       end do
c scr=SC_eff
       do i=1,ieffdim
        do mu=1,igr
         x=0.d0
         do nu=1,igr 
          x=x+s(mu,nu)*p0(nu,i)
         end do
         scr(mu,i)=x
        end do
       end do
       do i=1,ieffdim
        do k=1,ieffdim
         x=0.d0
         do mu=1,igr
          x=x+p0(mu,k)*scr(mu,i)
         end do
         s_eff(k,i)=x
        end do
       end do
C  Calculate s_eff^(-1/2) matrix
      call build_Smp(ieffdim,s_eff,sm,dum,0)
c orthogonalize 
       do mu=1,igr    
        do j=1,ieffdim
         x=0.0d0
         do k=1,ieffdim
          x=x+p0(mu,k)*sm(k,j)
         end do
         scr(mu,j)=x
        end do
       end do
       do mu=1,igr    
        do j=1,ieffdim
         p0(mu,j)=scr(mu,j)
        end do
       end do
c C^-1 = C^T (C C^T)^-1
       allocate(a_eff(igr,igr))
       allocate(u(igr,igr))
       do i=1,igr      
        do j=1,igr  
         x=0.d0
         do k=1,ieffdim     
          x=x+p0(i,k)*p0(j,k)
         end do
         u(i,j)=x
        end do
       end do
       call invert(igr,u)
       do i=1,effdim   
        do j=1,igr     
         x=0.d0
         do k=1,igr         
          x=x+p0(k,i)*u(k,j)
         end do
         a_eff(i,j)=x
        end do
       end do
c P_eff=U^-1 P U^-T
       do i=1,ieffdim  
        do j=1,ieffdim    
         x=0.d0
         do mu=1,igr         
         do nu=1,igr         
          x=x+a_eff(i,mu)*p(mu,nu)*a_eff(j,nu)    
         end do
         end do
         p_eff(i,j)=x
        end do
       end do
c min bas is orthogonal so
       x=0.0d0
       do i=1,ieffdim  
        x=x+p_eff(i,i)
       end do
       write(*,*) 'pollas trace p_eff',x
c
      illim_eff(1)=1
      iulim_eff(1)=ip0(1)
      do iatom=2,nat
       illim_eff(iatom)=iulim_eff(iatom-1)+1
       iulim_eff(iatom)=iulim_eff(iatom-1)+ip0(iatom)
      end do

      x=0.0d0
      do iatom=1,nat   
       xx=0.0d0
       do i=illim_eff(iatom),iulim_eff(iatom)
        xx=xx+p_eff(i,i)
       end do
       write(*,*) 'Population atom :',iatom,xx
       x=x+xx
      end do
      write(*,*) 'Total population :',x


       return
       end

      subroutine mhg(itotps,ndim,omp,chp,sat,wp,omp2,pcoord,pk,icase)
      use basis_set
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      integer,intent(in) :: itotps,ndim
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /atlist/iatlist(maxat),icuat
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat)
      dimension wp(itotps)
      dimension pcoord(itotps,3)
      dimension sat(ndim,ndim,nat)
      dimension pk(ndim,ndim)
      character*(5) key
c
      dimension effrho(maxat),effrhoacc(maxat)
      allocatable s0(:,:), sm(:,:), c0(:,:), splus(:,:), pp0(:,:)
      allocatable s0all(:),is0all(:)
      allocatable ::  scr(:)

      ihirsh = Iopt(6) 
      iallpo= iopt(7) 
      ieffao=  Iopt(12) 
      icube  = Iopt(13) 
      ieffthr  = Iopt(24) 
      iatps=nang*nrad

      xmaxocc=real(ieffthr)/1000.0d0 
      key=" "
      if (icase.eq.1) then
       key="ALPHA"
      else if(icase.eq.2) then
        key="BETA "
      end if

      print *,' '
      print *,' --------------------------------------'
      print *,'   DOING EFFAO-3D GENERAL FORMULATION'
      print *,' --------------------------------------'
      print *,' '
      print *,'     ',key,'  PART  '

      jocc=0
      xmaxotot=0.d0

      allocate(scr(iatps*nat))
      allocate (s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      allocate (c0(ndim,ndim), pp0(ndim,ndim), is0all(ndim))

      xnn=2.0d0
c external llop over all atoms
      do iicenter=1,nat

        do mu=1,ndim
         do nu=1,mu
           x=0.d0
           do jcenter=1,nat
            do ifut=iatps*(jcenter-1)+1,iatps*jcenter
             distx=pcoord(ifut,1)-coord(1,iicenter)
             disty=pcoord(ifut,2)-coord(2,iicenter)
             distz=pcoord(ifut,3)-coord(3,iicenter)
             distzz=(distx*distx+disty*disty+distz*distz)**(xnn/2.0d0)
             x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*omp(ifut)*distzz
            end do
           end do
           s0(mu,nu)=x
           s0(nu,mu)=x
         enddo
        enddo

c 
        call build_Smp(igr,S0,Sm,Splus,0)

c 
        do i=1,igr
         do j=1,igr
         pp0(i,j)=pk(i,j)
         end do
        end do
c
        call to_lowdin_basis(igr,Splus,pp0)
        call diagonalize(igr,igr,pp0,C0,0)
        call to_AO_basis(igr,igr,Sm,C0)

C But max number of effaos in one center is nocc
       
        xmaxo=0.0d0
        do i=1,nocc 
         xmaxo=xmaxo+pp0(i,i)
        end do
        imaxo=nocc
        write(*,*) '                  '
        write(*,'(a,i4,a)') ' ** FRAGMENT ',iicenter,' ** '
        write(*,*) '                  '

        write(*,'(a29,i3,f10.5)') ' Spreads for fragment ',iicenter, xmaxo
        write(*,60) (pp0(mu,mu),mu=1,imaxo)
60      format(7h OCCUP.   ,6F12.7)

c
        do k=1,imaxo           
         do mu=1,igr                       
           c0(mu,k)=c0(mu,k)*sqrt(pp0(k,k))
         end do
        end do
c ...calculate gross occupations from orbitals!!!

        xx0=0.0d0
        do i=1,imaxo
         xx1=0.0d0
         xx=0.0d0 
         do j=1,igr                        
          do k=1,igr
           xx=xx+c0(k,i)*sat(k,j,iicenter)*c0(j,i)
           xx1=xx1+c0(k,i)*s(k,j)*c0(j,i)
          end do
         end do
         write(*,*) 'normalization',iicenter,i,xx1
         s0all(i)=xx*2.0d0
         xx0=xx0+xx*2.0d0
        end do
        do k=1,imaxo           
         do mu=1,igr                       
           p0(mu,k)=c0(mu,k)
         end do
         p0net(k,iicenter)=pp0(k,k)
         p0gro(k,iicenter)=s0all(k)
        end do
        ip0(iicenter)=imaxo
         
        write(*,*) '                  '
        write(*,'(a31,i3,f10.5)') ' Gross occupation for fragment ',iicenter, xx0
        write(*,60) (s0all(mu),mu=1,imaxo)



c  write cube file
        if(icube.eq.1) call cubegen3_mhg(iicenter,icase)


c end outer loop over fragments
        end do
       deallocate (scr,s0, sm, c0,splus,pp0, s0all,is0all)

         end
      subroutine mhg2(itotps,ndim,omp,chp,sat,wp,omp2,pcoord,pk,icase)
      use basis_set
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      integer,intent(in) :: itotps,ndim
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /atlist/iatlist(maxat),icuat
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat)
      dimension wp(itotps)
      dimension pcoord(itotps,3)
      dimension sat(ndim,ndim,nat)
      dimension pk(ndim,ndim)
      character*(5) key
c
      dimension effrho(maxat),effrhoacc(maxat)
      allocatable s0(:,:), sm(:,:), c0(:,:), splus(:,:), pp0(:,:)
      allocatable s0all(:),is0all(:)
      allocatable ::  scr(:)

      ihirsh = Iopt(6) 
      iallpo= iopt(7) 
      ieffao=  Iopt(12) 
      icube  = Iopt(13) 
      ieffthr  = Iopt(24) 
      iatps=nang*nrad

      xmaxocc=real(ieffthr)/1000.0d0 
      key=" "
      if (icase.eq.1) then
       key="ALPHA"
      else if(icase.eq.2) then
        key="BETA "
      end if

      print *,' '
      print *,' --------------------------------------'
      print *,'   DOING EFFAO-3D GENERAL FORMULATION'
      print *,' --------------------------------------'
      print *,' '
      print *,'     ',key,'  PART  '

      jocc=0
      xmaxotot=0.d0

      allocate(scr(iatps*nat))
      allocate (s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      allocate (c0(ndim,ndim), pp0(ndim,ndim), is0all(ndim))

      xnn=2.0d0
c external llop over all atoms
      do mu=1,ndim
       do nu=1,mu
        x=0.d0
         do jcenter=1,nat
          do ifut=iatps*(jcenter-1)+1,iatps*jcenter
           do iicenter=1,nat
             distx=pcoord(ifut,1)-coord(1,iicenter)
             disty=pcoord(ifut,2)-coord(2,iicenter)
             distz=pcoord(ifut,3)-coord(3,iicenter)
             distzz=(distx*distx+disty*disty+distz*distz)**(xnn/2.0d0)
             x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*omp(ifut)/distzz
           end do
          end do
         end do
         s0(mu,nu)=x
         s0(nu,mu)=x
         enddo
        enddo

c make SAA1/2
        call build_Smp(igr,S0,Sm,Splus,0)

c tranform block S with P0
        do i=1,igr
         do j=1,igr
         pp0(i,j)=pk(i,j)
         end do
        end do
        call to_lowdin_basis(igr,Splus,pp0)
        call diagonalize(igr,igr,pp0,C0,0)
        call to_AO_basis(igr,igr,Sm,C0)

C But max number of effaos in one center is nocc
       
        xmaxo=0.0d0
        do i=1,nocc 
         xmaxo=xmaxo+pp0(i,i)
        end do
        imaxo=nocc
        write(*,*) '                  '
        write(*,'(a)') ' ** MOLECULE  ** '
        write(*,*) '                  '

        write(*,'(a29,f10.5)') ' Spreads  ', xmaxo
        write(*,60) (pp0(mu,mu),mu=1,imaxo)
60      format(7h OCCUP.   ,6F12.7)

        do k=1,imaxo           
         do mu=1,igr                       
           c0(mu,k)=c0(mu,k)*sqrt(pp0(k,k)/2.0d0)
         end do
        end do
c ...calculate gross occupations from orbitals!!!

        xx0=0.0d0
        do i=1,imaxo
         xx1=0.0d0
         xx=0.0d0 
c         do iicenter=1,nat
         do j=1,igr                        
          do k=1,igr
c           xx=xx+c0(k,i)*sat(k,j,iicenter)*c0(j,i)
           xx1=xx1+c0(k,i)*s(k,j)*c0(j,i)
          end do
         end do
c         end do
         write(*,*) 'normalization',iicenter,i,xx1
         s0all(i)=xx1*2.0d0
         xx0=xx0+xx*2.0d0
        end do

        iicenter=1
        do k=1,imaxo           
         do mu=1,igr                       
           p0(mu,k)=c0(mu,k)
         end do
         p0net(k,iicenter)=pp0(k,k)
         p0gro(k,iicenter)=s0all(k)
        end do
        ip0(iicenter)=imaxo
         
c        write(*,*) '                  '
c        write(*,'(a31,i3,f10.5)') ' Gross occupation for fragment ',iicenter, xx0
c        write(*,60) (s0all(mu),mu=1,imaxo)



c  write cube file
        if(icube.eq.1) call cubegen3_mhg(iicenter,icase)


       deallocate (scr,s0, sm, c0,splus,pp0, s0all,is0all)

         end
! *****

      subroutine eos_centroid(itotps,chp,wp,omp,pcoord)
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      parameter (toau=0.52917721067d0)
      integer, intent(in) :: itotps
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
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
        stop
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
      print *,'  '
      print *,'        ALPHA PART'
      print *,'  '
      print *,'  Frag  Na(geom)  Na(QTAIM) '
      print *,' ---------------------------'
      do i=1,icufr
       itot1=itot1-elect(i,1)
       itot2=itot2-elect(i,2)
       write(*,'(3x,i3,5x,i4,5x,i4)') i,-int(elect(i,1)),-int(elect(i,2))
      end do
      print *,' ---------------------------'
      write(*,'(a6,2i9)') '  Sum ', itot1,itot2
      print *,'  ' 

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
        stop
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
      print *,'  '
      print *,'   OXIDATION STATES FROM    '
      print *,' CENTROIDS OF LOCALIZED MOs '
      print *,'  '
      print *,'         BETA PART'
      print *,'  '
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
      print *,'  ' 

      end if

c OS
      do iat=1,nat 
       elect(jfrlist(iat),1)=elect(jfrlist(iat),1)+zn(iat)
       elect(jfrlist(iat),2)=elect(jfrlist(iat),2)+zn(iat)
      end do

      itot1=0
      itot2=0
      print *,'  '
      print *,'   OXIDATION STATES FROM    '
      print *,' CENTROIDS OF LOCALIZED MOs '
      print *,'  '
      print *,'  Frag  OS(geom)  OS(QTAIM) '
      print *,' ---------------------------'
      do i=1,icufr
       itot1=itot1+elect(i,1)
       itot2=itot2+elect(i,2)
       write(*,'(3x,i3,5x,i4,5x,i4)') i,int(elect(i,1)),int(elect(i,2))
      end do
      print *,' ---------------------------'
      write(*,'(a6,2i9)') '  Sum ', itot1,itot2
      print *,'  ' 

      end




