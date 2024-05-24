!! **************************** !!
!! CONVENTIONAL EOS SUBROUTINES !!
!! **************************** !!

!! ***** !!

      subroutine ueffao3d_frag(itotps,ndim,omp,chp,sat,wp,omp2,pk,icase)
      
      use integration_grid
      
      implicit real*8(a-h,o-z)
      
      include 'parameter.h'
      
      integer,intent(in) :: itotps,ndim

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common /atlist/iatlist(maxat),icuat
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/qat(maxat,2),qsat(maxat,2)

      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat)
      dimension wp(itotps)
      dimension sat(ndim,ndim,nat)
      dimension pk(ndim,ndim)

      allocatable :: s0(:,:),sm(:,:),c0(:,:),splus(:,:),pp0(:,:)
      allocatable :: s0all(:),is0all(:)
      allocatable :: scr(:)

      ihirsh  = Iopt(6) 
      iallpo  = iopt(7) 
      ieffao  = Iopt(12) 
      icube   = Iopt(13) 
      ieffthr = Iopt(24) 
      iatps   = nang*nrad

      xminocc=REAL(ieffthr)/1000.0d0 

      write(*,*) " "
      write(*,*) " ------------------------------------ "
      write(*,*) "  DOING EFFAO-3D GENERAL FORMULATION  "
      write(*,*) " ------------------------------------ "
      write(*,*) " "
      if(icase.eq.1) then
        write(*,*) " ------------------------------- "
        write(*,*) "  EFFAOs FROM THE ALPHA DENSITY  "
        write(*,*) " ------------------------------- "
      else if(icase.eq.2) then
        write(*,*) " ------------------------------ "
        write(*,*) "  EFFAOs FROM THE BETA DENSITY  "
        write(*,*) " ------------------------------ "
      end if
      write(*,*) " "
      
      jocc=0
      xmaxotot=ZERO

      ALLOCATE(scr(iatps*nat))
      ALLOCATE(s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      ALLOCATE(c0(ndim,ndim),pp0(ndim,ndim),is0all(ndim))
      do iicenter=1,icufr
        do ifut=1,iatps*nat                         
          scr(ifut)=ZERO
          do icenter=1,nfrlist(iicenter)
            scr(ifut)=scr(ifut)+ omp2(ifut,ifrlist(icenter,iicenter))
          end do
        end do

!! Computing Becke atomic NET orbital overlap !!
!! ALLPOINTS should be default here !!
!! a distance based screening would be interesting for very large systems !!
        iallpo0=1
        if(iallpo0.eq.1) then
          do mu=1,ndim
            do nu=1,mu
              x=ZERO
              do jcenter=1,nat
                do ifut=iatps*(jcenter-1)+1,iatps*jcenter
                  x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*scr(ifut)*scr(ifut)*omp(ifut)
                end do
              end do
              s0(mu,nu)=x
              s0(nu,mu)=x
            end do
          end do
        else
          do mu=1,ndim
            do nu=1,mu
              x=ZERO
              do icenter=1,nfrlist(iicenter)
                jcenter=ifrlist(icenter,iicenter)
                do ifut=iatps*(jcenter-1)+1,iatps*jcenter
                  x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*scr(ifut)*scr(ifut)*omp(ifut)
                end do
              end do
              s0(mu,nu)=x
              s0(nu,mu)=x
            end do
          end do
        end if

!! make SAA1/2 !!
        call build_Smp(igr,S0,Sm,Splus,0)

!! transform block S with P0 !!
        do i=1,igr
          do j=1,igr
            pp0(i,j)=pk(i,j)
          end do
        end do
        call to_lowdin_basis(igr,Splus,pp0)
        call diagonalize(igr,igr,pp0,C0,0)
        call to_AO_basis(igr,igr,Sm,C0)

!! But max number of effaos in one center is nocc !!
        i=1
        do while(pp0(i,i).ge.xminocc.and.i.le.igr) 
          imaxo=i
          i=i+1
        end do
        xmaxo=ZERO
        do i=1,igr  
          xmaxo=xmaxo+pp0(i,i)
        end do

        xx0=ZERO
        xx1=ZERO
        do icenter=1,nfrlist(iicenter)
          xx1=xx0+qat(ifrlist(icenter,iicenter),1)
          do jcenter=1,nfrlist(iicenter)
            xx0=xx0+op(ifrlist(icenter,iicenter),ifrlist(jcenter,iicenter))
          end do
        end do
        write(*,'(2x,a11,x,i3,x,a2)') "** FRAGMENT",iicenter,"**"
        write(*,*) " "
        if(icase.eq.0) write(*,'(2x,a29,x,f8.4)') "Deviation from net population",xmaxo-xx0 !! ITS NEVER ZERO, PEDRO RECHECK main.f AND MODIFY IN GENERAL !!
        write(*,'(2x,a27,x,i3,x,f10.5)') "Net occupation for fragment",iicenter,xmaxo
        write(*,'(2x,a22,x,f10.5)') "Net occupation using >",xminocc
        write(*,60) (pp0(mu,mu),mu=1,imaxo)
        write(*,*) " "

!! ...calculate gross occupations from orbitals !!

        xx0=ZERO
        do i=1,imaxo
          xxx=ZERO
          xxy=ZERO
          do icenter=1,nfrlist(iicenter)
            jcenter=ifrlist(icenter,iicenter)
            xx=ZERO
            do j=1,igr                        
              do k=1,igr
                xx=xx+c0(k,i)*sat(k,j,jcenter)*c0(j,i)
              end do
            end do
            xxx=xxx+xx
          end do
          xxx=xxx*pp0(i,i)
          s0all(i)=xxx
          xx0=xx0+xxx
        end do
        write(*,*) " "
        write(*,'(2x,a29,x,i3,f10.5)') "Gross occupation for fragment",iicenter,xx0
        if(icase.eq.0) write(*,'(2x,a31,x,f8.4)') "Deviation from gross population",xx0-xx1
        write(*,60) (s0all(mu),mu=1,imaxo)
        write(*,*) " "

        do k=1,imaxo           
          do mu=1,igr                       
            p0(mu,k)=c0(mu,k)
          end do
          p0net(k,iicenter)=pp0(k,k)
          p0gro(k,iicenter)=s0all(k)
        end do
        ip0(iicenter)=imaxo

!!  write cube file !!
        if(icube.eq.1) call cubegen4(iicenter,icase)

!! end outer loop over fragments !!
      end do

!! DEALLOCATING !!
      DEALLOCATE(scr,s0,sm,c0,splus,pp0,s0all,is0all)

!! PRINTING FORMATS !!
60    FORMAT(8h  OCCUP. ,8f9.4)

      end

!! ****** !!

      subroutine eos_analysis(idobeta,icase,thres)

      implicit real*8(a-h,o-z)
      
      include 'parameter.h'
      
      common /effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /atlist/iatlist(maxat),icuat
      common /iops/iopt(100)
      common /loba/ oxi(maxat),errsav(maxat),elec(maxat),effpop(maxat)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /loba2/occup(nmax,2),iorbat(nmax,2),lorb(2),confi0

      dimension occup2(igr), occupg(igr), errnet(maxat)
      dimension qdev(maxat)

      ieffthr = Iopt(24) 
      icorr   = Iopt(26) 

      if(idobeta.eq.0.and.icase.eq.2) then
        write(*,*) " "
        write(*,*) " SKIPPING EFFAOs FOR BETA ELECTRONS "
        write(*,*) " "
        lorb(2)=lorb(1)
        do i=1,lorb(2)
          occup(i,2)=occup(i,1) 
          iorbat(i,2)= iorbat(i,1)
        end do  
        do i=1,icufr
          errnet(i)=errsav(i)
          effpop(i)=effpop(i)*2.0d0
        end do
        confi=confi0
        go to 99
      else if(icase.eq.2.and.nb.eq.0) then
        write(*,*) " "
        write(*,*) " CALCULATION HAS NO BETA ELECTRONS "
        write(*,*) " "
        do i=1,icufr
          elec(i)=ZERO
        end do
        go to 99
      end if

      iorb=0
      do i=1,icufr
        elec(i)=ZERO
        do k=1,ip0(i)
          iorb=iorb+1
          occup(iorb,icase)=p0net(k,i)
          iorbat(iorb,icase)=i
        end do 
      end do  
      write(*,*) " "
      write(*,'(2x,a38,x,i4)') "Total number of eff-AO-s for analysis:",iorb

      lorb(icase)=iorb

      do i=1,iorb-1
        do j=i+1,iorb
          if (occup(j,icase).gt.occup(i,icase))then
            xkk=occup(j,icase)
            occup(j,icase)=occup(i,icase)
            occup(i,icase)=xkk
            xkk=occupg(j)
            occupg(j)=occupg(i)
            occupg(i)=xkk
            ikk=iorbat(j,icase)
            iorbat(j,icase)=iorbat(i,icase)
            iorbat(i,icase)=ikk 
          end if
        end do  
      end do 

      k=0
      nnn=nalf
      if(icase.eq.2) nnn=nb   
333   k=k+1
      if(dabs(occup(nnn,icase)-occup(nnn+k,icase)).lt.thres) go to 333
      k=k-1
      if(k.eq.0) then
        write(*,*) " EOS: Unambiguous integer electron assignation "
        do i=1,iorb
          if(i.le.nnn) then
            occup2(i)=ONE
          else
            occup2(i)=ZERO
          end if
        end do 
      else
        write(*,*) ' ******************************************'
        write(*,*) ' EOS: Warning, pseudo-degeneracies detected'
        write(*,*) ' ******************************************'
        kk=0
334     kk=kk+1
        if(dabs(occup(nnn,icase)-occup(nnn-kk,icase)).lt.thres) go to 334
        kk=kk-1
        frac=float(kk+1)/float(kk+k+1)
        write(*,'(2x,a12,x,i4,x,a14,x,i4,x,a30)') "Distributing",kk+1,"electrons over",kk+k+1,
     +  "pseudodegenerate atomic orbitals"
        do i=1,iorb
          if(i.lt.nnn-kk) then
            occup2(i)=ONE
          else if(i.gt.nnn+k) then
            occup2(i)=ZERO
          else
            occup2(i)=frac 
          end if
        end do 
      end if
     
      do i=1,iorb
        elec(iorbat(i,icase))=elec(iorbat(i,icase))+occup2(i)
      end do

!! PRINTING INFO !!
      if(icase.eq.1) then
        write(*,*) " "
        write(*,*) " ---------------------------------- "
        write(*,*) "  EOS ANALYSIS FOR ALPHA ELECTRONS  "
        write(*,*) " ---------------------------------- "
        do i=1,icuat
          errsav(i)=errnet(i)
          effpop(i)=ZERO
        end do
      else if (icase.eq.2) then
        write(*,*) " "
        write(*,*) " --------------------------------- "
        write(*,*) "  EOS ANALYSIS FOR BETA ELECTRONS  "
        write(*,*) " --------------------------------- "
      end if
      write(*,*) " "
      write(*,*) "  Frag.  Elect.  Last occ.  First unocc.  "
      write(*,*) " ---------------------------------------- "

      xlast=ONE
      ilast=0
      do i=1,icufr
        nn=int(elec(i))
        if(elec(i)-nn.gt.thresh) nn=nn+1
        if(nn+1.gt.ip0(i)) then 
          write(*,10) i,elec(i),p0net(nn,i)
        else
          if(nn.ne.0) then
            write(*,15) i,elec(i),p0net(nn,i),p0net(nn+1,i)
          else
            write(*,15) i,elec(i),ZERO,p0net(nn+1,i)
          end if
        end if
        if(nn.ne.0) then
          if(p0net(nn,i).lt.xlast) then
            xlast=p0net(nn,i)
            ilast=i
          end if
        end if
        do j=1,nn
          effpop(i)=effpop(i)+p0net(j,i)
        end do
      end do 

      xfirst=ZERO
      do i=1,icufr
        if(i.ne.ilast) then
          nn=int(elec(i))
          if(elec(i)-nn.gt.thresh) nn=nn+1
          if(nn+1.ne.0.and.p0net(nn+1,i).gt.xfirst) xfirst=p0net(nn+1,i)
        end if
      end do
      write(*,*) " ---------------------------------------- "

      confi=100.0*min(1.0d0,xlast-xfirst+0.5d0)
      write(*,'(3x,a24,x,f7.3)') "RELIABILITY INDEX R(%) =",confi
      if(icase.eq.1) confi0=confi

99    continue

      zztot=ZERO
      do i=1,icufr
        if(icase.eq.1) then
          zzn=ZERO
          do j=1,nfrlist(i)
            zzn=zzn+zn(ifrlist(j,i))
          end do
          oxi(i)=zzn-elec(i)
        else
          oxi(i)=oxi(i)-elec(i)
          zztot=zztot+oxi(i)
        end if
      end do

!! FINAL PRINTING !!
      if(icase.eq.2) then
        write(*,*) " "
        write(*,*) " --------------------------- "
        write(*,*) "  FRAGMENT OXIDATION STATES  "
        write(*,*) " --------------------------- "
        write(*,*) " "
        write(*,*) "  Frag.  Oxidation State  "
        write(*,*) " ------------------------ "
        do ifrg=1,icufr
          write(*,20) ifrg,oxi(ifrg)
        end do 
        write(*,*) " ------------------------ "
        write(*,'(3x,a4,x,f4.1)') "Sum:",zztot 
        write(*,*) " "

!! TO BE CHANGED !!
        confi2=dmin1(confi,confi0)
        write(*,'(2x,a32,x,f7.3)') "OVERALL RELIABILITY INDEX R(%) =",confi2
      end if

!! PRINTING FORMATS !!
10    FORMAT(3x,i3,3x,f6.2,4x,f12.3,'    < thresh',2f12.3) !tocheck!
15    FORMAT(3x,i3,3x,f6.2,4x,f6.3,4x,f6.3)
20    FORMAT(3x,i3,6x,f8.2)

      end 

!! ****** !!

      subroutine ueffaolow_frag(icase)

      use basis_set
      use ao_matrices

      implicit real*8(A-H,O-Z)

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /atlist/iatlist(maxat),icuat
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /nao/unao(nmax,nmax),ssnao(nmax,nmax)

      dimension iao_frag(nmax)

      character*(25) nameout

      allocatable :: s0(:,:), sm(:,:), c0(:,:), splus(:,:), pp0(:,:)
      allocatable :: s0all(:),is0all(:),efo(:,:),efo2(:)

      ndim    = igr
      icube   = Iopt(13)
      ieffthr = Iopt(24)
      imulli  = Iopt(5)

      iefo=0
      xminocc=REAL(ieffthr)/1000.0d0

      ALLOCATE(s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      ALLOCATE(c0(ndim,ndim),pp0(ndim,ndim),is0all(ndim))
      ALLOCATE(efo(ndim,ndim),efo2(ndim))

c ao to frag map
      iao_frag=0
      do ifrag=1,icufr
        do icenter=1,nfrlist(ifrag)
          jcenter=ifrlist(icenter,ifrag)
          do mu=llim(jcenter),iulim(jcenter)
            iao_frag(mu)=ifrag
          end do
        end do
      end do

c initializing
      if(icase.eq.0) then
        s0=p
      else if (icase.eq.1) then
        s0=pa
      else if (icase.eq.2) then
        s0=pb
      end if

      if(imulli.eq.4) then
c use NA0 to AO matrix
c warning, transpose...nonsymmetric
        do i=1,igr
          do j=1,igr
            splus(i,j)=ssnao(j,i)
            sm(i,j)=unao(i,j)
          end do
        end do
      else
        splus=s12p
        sm=s12m
      end if
      if(icase.eq.1) then

!! TO PRINT ONLY ONCE !!
        if(imulli.eq.4) then
          write(*,*) " "
          write(*,*) " ----------------------------- "
          write(*,*) "  DOING EFFAO NAO FORMULATION  "
          write(*,*) " ----------------------------- "
          write(*,*) " "
        else
          write(*,*) " "
          write(*,*) " -------------------------------- "
          write(*,*) "  DOING EFFAO LOWDIN FORMULATION  "
          write(*,*) " -------------------------------- "
        end if
        write(*,*) " "
        write(*,*) " ------------------------------- "
        write(*,*) "  EFFAOs FROM THE ALPHA DENSITY  "
        write(*,*) " ------------------------------- "
      else if(icase.eq.2) then
        write(*,*) " "
        write(*,*) " ------------------------------ "
        write(*,*) "  EFFAOs FROM THE BETA DENSITY  "
        write(*,*) " ------------------------------ "
      end if
      write(*,*) " "
     
c tranform P with Splus
C will back trasnform effaos to ao basis later...then expanded in whole basis functions
      call to_lowdin_basis(igr,splus,s0)

C LOOP OVER FRAGMENTS
c easiest way to avois redordering is to perform nat diagonaliztions
c jocc will control the total number of efos to process
      jocc=0
      do ifrag=1,icufr
c Block diagonal P
c clean auxiliary pp0 mat
        pp0=0.0d0
        do mu=1,igr
          do nu=1,igr
            if(iao_frag(mu).eq.ifrag.and.iao_frag(nu).eq.ifrag) pp0(mu,nu)=s0(mu,nu)
          end do
        end do

        call diagonalize(igr,igr,pp0,C0,0)
c backtrasnform. eff.aos expanded over the full basis set
        call to_AO_basis(igr,igr,sm,C0)

c max number of effaos
        imaxeff=0
        do i=1,igr
          if(iao_frag(i).eq.ifrag) imaxeff=imaxeff+1
        end do

c actual number of effaos
        xmaxo=ZERO
        imaxo=0
        i=1
        do while(pp0(i,i).ge.xminocc.and.i.lt.imaxeff) 
          xmaxo=xmaxo+pp0(i,i)
          imaxo=i
          i=i+1
        end do

!! PRINTING !!
        write(*,'(2x,a11,x,i3,x,a2)') "** FRAGMENT",ifrag,"**"
        write(*,*) " "
        write(*,'(2x,a27,x,i3,x,f10.5)') "Net occupation for fragment",ifrag,xmaxo
        write(*,'(2x,a22,x,f10.5)') "Net occupation using >",xminocc
        write(*,60) (pp0(mu,mu),mu=1,imaxo)
        write(*,*) " "

c saving efo info for fragment
c saving all efos and occupations in efo matrix
        do k=1,imaxo           
          iefo=iefo+1
          do mu=1,igr                       
            p0(mu,k)=c0(mu,k)
            efo(mu,iefo)=c0(mu,k)
          end do
          p0net(k,ifrag)=pp0(k,k)
          p0gro(k,ifrag)=pp0(k,k)
          efo2(iefo)=pp0(k,k)
        end do
        ip0(ifrag)=imaxo

c OUTPUT ORBITALS FOR VISUALIZATION
c  write cube file
        if(icube.eq.1) call cubegen4(ifrag,icase)

c end loop over fragments
      end do

!! MG: THIS DESERVES BETTER (TO DO) !!
c print out efo matrix
      !write(*,*) 'printing ',iefo,' efos'
      ival=iefo*igr
      nameout='efo_occ.dat'
      nameout=adjustl(nameout)
      if(icase.eq.1) then
        open(unit=44,file=nameout) 
        write(44,'(A49,I12)') "Alpha Orbital Energies                     R   N=",iefo
        write(44,'(5ES16.8)') (efo2(i),i=1,iefo)
c
        nameout='efo_coeff.dat'
        nameout=adjustl(nameout)
        open(unit=45,file=nameout) 
        write(45,'(A49,I12)') "Alpha MO coefficients                      R   N=",ival
        write(45,'(5ES16.8)') ((efo(i,j),i=1,igr),j=1,iefo)
      else
        write(44,'(A49,I12)') "Beta Orbital Energies                      R   N=",iefo
        write(44,'(5ES16.8)') (efo2(i),i=1,iefo)
        write(45,'(A49,I12)') "Beta MO coefficients                       R   N=",ival
        write(45,'(5ES16.8)') ((efo(i,j),i=1,igr),j=1,iefo)
      end if
      if(icase.eq.2) then
        close(44)
        close(45)
      end if

!! DEALLOCATING !!
      DEALLOCATE(efo,efo2)
      DEALLOCATE(s0,sm,c0,splus,pp0,s0all,is0all)

!! PRINTING FORMATS !!
60    FORMAT(8h  OCCUP. ,8f9.4)

      end

!! ****** !!

      subroutine ueffaomull_frag(icase)

      use basis_set
      use ao_matrices

      implicit real*8(A-H,O-Z)

      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /atlist/iatlist(maxat),icuat
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common  /filename/name0

      character*(60) name0,name
      character*(80) line  

      dimension iao_frag(nmax)

      allocatable :: s0(:,:),sm(:,:),c0(:,:),splus(:,:),pp0(:,:)
      allocatable :: s0all(:),is0all(:),pk(:,:)

      ndim    = igr
      icube   = Iopt(13)
      ieffthr = Iopt(24)

      xminocc=real(ieffthr)/1000.0d0

      ALLOCATE(s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      ALLOCATE(c0(ndim,ndim),pp0(ndim,ndim),is0all(ndim),pk(ndim,ndim)) 

!! INITIALIZING !!
      if(icase.eq.0) then
        pk=p
      else if(icase.eq.1) then
        write(*,*) " "
        write(*,*) " ---------------------------------- "
        write(*,*) "  DOING EFFAO MULLIKEN FORMULATION  "
        write(*,*) " ---------------------------------- "
        write(*,*) " "
        write(*,*) " ------------------------------- "
        write(*,*) "  EFFAOs FROM THE ALPHA DENSITY  "
        write(*,*) " ------------------------------- "
        pk=pa
      else if(icase.eq.2) then
        write(*,*) " ------------------------------ "
        write(*,*) "  EFFAOs FROM THE BETA DENSITY  "
        write(*,*) " ------------------------------ "
        pk=pb
      end if
      write(*,*) " "

c ao to frag map
      iao_frag=0
      do ifrag=1,icufr
        do icenter=1,nfrlist(ifrag)
          jcenter=ifrlist(icenter,ifrag)
          do mu=llim(jcenter),iulim(jcenter)
            iao_frag(mu)=ifrag
          end do
        end do
      end do

C LOOP OVER FRAGMENTS
c easiest way to avois redordering is to perform nat diagonaliztions
c jocc will control the total number of efos to process
      jocc=0
      do ifrag=1,icufr
c Making block-diagonal S  and P matrix (AO)
        pp0=ZERO
        s0=ZERO
        do mu=1,igr
          do nu=1,igr
            if(iao_frag(mu).eq.ifrag.and.iao_frag(nu).eq.ifrag) then 
              pp0(mu,nu)=pk(mu,nu)
              s0(mu,nu)=s(mu,nu)
            end if
          end do
        end do
c make S0^1/2
        call build_Smp(igr,s0,Sm,Splus,0)
c tranform blovck P0 with Splus
        call to_lowdin_basis(igr,Splus,pp0)
        call diagonalize(igr,igr,pp0,C0,0)
c backtrasnform. eff.aos expanded over the full basis set
        call to_AO_basis(igr,igr,Sm,C0)
c max number of effaos
        imaxeff=0
        do i=1,igr
          if(iao_frag(i).eq.ifrag) imaxeff=imaxeff+1
        end do
c actual number of effaos
        xmaxo=ZERO
        imaxo=0
        i=1
        do while(pp0(i,i).ge.xminocc.and.i.lt.imaxeff) 
          xmaxo=xmaxo+pp0(i,i)
          imaxo=i
          i=i+1
        end do

!! PRINTING !!
        write(*,'(2x,a11,x,i3,x,a2)') "** FRAGMENT",ifrag,"**"
        write(*,*) " "
        write(*,'(2x,a27,x,i3,x,f10.5)') "Net occupation for fragment",ifrag,xmaxo
        write(*,'(2x,a22,x,f10.5)') "Net occupation using >",xminocc
        write(*,60) (pp0(mu,mu),mu=1,imaxo)
        write(*,*) " "

c saving efo info for fragment
        do k=1,imaxo           
          do mu=1,igr                       
            p0(mu,k)=c0(mu,k)
          end do
          p0net(k,ifrag)=pp0(k,k)
          p0gro(k,ifrag)=pp0(k,k)
        end do
        ip0(ifrag)=imaxo

c OUTPUT ORBITALS FOR VISUALIZATION
c  write cube file
        if(icube.eq.1) call cubegen4(ifrag,icase)

c end loop over fragments
      end do

!! DEALLOCATING !!
      DEALLOCATE(s0, sm, c0,splus,pp0, s0all,is0all,pk)

!! PRINTING FORMATS !!
60    FORMAT(8h  OCCUP. ,8f9.4)

      end

!! ****** !!
!! MG: NOT TOUCHED FROM HERE !!
!! ****** !!

      subroutine uefomo(itotps,ndim,omp,chp,sat,wp,omp2,icase)
      use integration_grid
      use ao_matrices, only :c,cb
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      integer,intent(in) :: itotps,ndim
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord0(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common /atlist/iatlist(maxat),icuat
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat)
      dimension wp(itotps)
      dimension sat(ndim,ndim,nat)
      character*(15) key
c
      allocatable cc0(:,:), c0(:,:),  pp0(:,:),s0all(:)
      allocatable scr(:),s0(:,:),caux(:,:)

      ihirsh = Iopt(6) 
      iallpo= iopt(7) 
      ieffao=  Iopt(12) 
      icube  = Iopt(13) 
      ieffthr  = Iopt(24) 
      iatps=nang*nrad

      xmaxocc=real(ieffthr)/1000.0d0 
      jocc=0
      xmaxotot=0.d0
      allocate (scr(itotps))
      allocate (s0(igr,igr),caux(igr,igr))


      key=" "
      if (icase.eq.1) then
       key="ALPHA"
       nbas=nalf
       caux=c
      else if(icase.eq.2) then
       key="BETA "
       nbas=nb
       caux=cb
      else
       key="ALPHA + BETA"
       nbas=nocc
       caux=c
      end if

      allocate(pp0(nbas,nbas),cc0(nbas,nbas),c0(igr,nbas),s0all(nbas))

      print *,' '
      print *,' --------------------------------------'
      print *,'   DOING EFFAO-3D GENERAL FORMULATION'
      print *,' --------------------------------------'
      print *,' '
      print *,'     ',key,'  PART  '


c Main loop over fragments
      do iicenter=1,icufr
c
c fragment weights (squared)
         do ifut=1,itotps                            
          scr(ifut)=0.0d0
          do icenter=1,nfrlist(iicenter)
           scr(ifut)=scr(ifut)+ omp2(ifut,ifrlist(icenter,iicenter))
          end do
          scr(ifut)=scr(ifut)*scr(ifut)
         end do

c Computing Becke atomic NET orbital overlap
C ALLPOINTS should be default here
c a distance based screening would be interesting for very large systems
        do mu=1,ndim
         do nu=1,mu
           x=0.d0
           do jcenter=1,nat
            do ifut=iatps*(jcenter-1)+1,iatps*jcenter
             if(scr(ifut).gt.1.0d-8) then
             x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*scr(ifut)*omp(ifut)
             end if
            end do
           end do
           s0(mu,nu)=x
           s0(nu,mu)=x
         enddo
        enddo
c        do mu=1,ndim
c         do nu=1,mu
c           x=0.d0
c           do icenter=1,nfrlist(iicenter)
c            jcenter=ifrlist(icenter,iicenter)
c            do ifut=iatps*(jcenter-1)+1,iatps*jcenter
c             x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*scr(ifut)*omp(ifut)
c            end do
c           end do
c           s0(mu,nu)=x
c           s0(nu,mu)=x
c         enddo
c        enddo

c now to MO
      do i=1,nbas 
       do j=i,nbas
        xx=0.0d0
        do k=1,igr
         do l=1,igr
          xx=xx+caux(k,i)*s0(k,l)*caux(l,j)
         end do
        end do
        pp0(i,j)=xx
        pp0(j,i)=xx
       end do
      end do

      call diagonalize(nbas,nbas,pp0,cc0,0)
c get efos in AO basis
      do i=1,igr
       do j=1,nbas
        c0(i,j)=0.0d0
        do k=1,nbas
         c0(i,j)=c0(i,j)+caux(i,k)*cc0(k,j)
        end do
       end do
      end do

      imaxo=0
      do i=1,nbas
       if(pp0(i,i).gt.xmaxocc) imaxo=imaxo+1
       if(icase.eq.0) pp0(i,i)=pp0(i,i)*2.0d0
      end do

      xmaxo=0.0d0
      do i=1,imaxo
       xmaxo=xmaxo+pp0(i,i)
      end do
c
      xx0=0.0d0
      xx1=0.0d0
      do icenter=1,nfrlist(iicenter)
       xx1=xx0+qat(ifrlist(icenter,iicenter),1)
       do jcenter=1,nfrlist(iicenter)
        xx0=xx0+op(ifrlist(icenter,iicenter),ifrlist(jcenter,iicenter))
       end do
      end do
      write(*,*) '                  '
      write(*,'(a,i4,a)') ' ** FRAGMENT ',iicenter,' ** '
      write(*,*) '                  '

      write(*,'(a29,i3,f10.5)') ' Net occupation for fragment ',iicenter, xmaxo
      if(icase.eq.0) then
       write(*,'(a30,f8.4)') ' Deviation from net population ', xmaxo-xx0 
      end if
      write(*,'(a23,f6.4)') ' Net occupation using >',xmaxocc
      write(*,60) (pp0(mu,mu),mu=1,imaxo)
60    format(7h OCCUP.   ,8F9.4)

c ...calculate gross occupations from orbitals!!!

        xx0=0.0d0
        do i=1,imaxo
         xxx=0.0d0
         do icenter=1,nfrlist(iicenter)
          jcenter=ifrlist(icenter,iicenter)
          xx=0.0d0 
          do j=1,igr                        
           do k=1,igr
            xx=xx+c0(k,i)*sat(k,j,jcenter)*c0(j,i)
           end do
          end do
          xxx=xxx+xx
         end do
         s0all(i)=xxx
         if(icase.eq.0) s0all(i)=2.0d0*s0all(i)
         xx0=xx0+s0all(i)
        end do
         
        write(*,*) '                  '
        write(*,'(a31,i3,f10.5)') ' Gross occupation for fragment ',iicenter, xx0
        if(icase.eq.0) then
         write(*,'(a32,f8.4)') ' Deviation from gross population ', xx0-xx1 
        end if
        write(*,60) (s0all(mu),mu=1,imaxo)

c saving efo info for fragment
        do k=1,imaxo           
         do mu=1,igr                       
           p0(mu,k)=c0(mu,k)
         end do
         p0net(k,iicenter)=pp0(k,k)
         p0gro(k,iicenter)=s0all(k)
        end do
        ip0(iicenter)=imaxo

c  write cube file
        if(icube.eq.1) call cubegen4(iicenter,icase)


c end outer loop over fragments
        end do
       deallocate (cc0,c0,caux,pp0,s0all,s0,scr)

         end
c                              


!!!
!!! OLD VERSIONS
!!!
      subroutine ueffaomull2(p,icase)
      use basis_set
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common/effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /atlist/iatlist(maxat),icuat
!      dimension p(nmax,nmax)
      character*(5) key

      allocatable s0(:,:), sm(:,:), c0(:,:), splus(:,:), pp0(:,:)
      allocatable s0all(:),is0all(:)
      
      ndim=igr

      allocate (s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      allocate (c0(ndim,ndim), pp0(ndim,ndim), is0all(ndim))

      icube=iopt(13)
      ieffthr  = Iopt(24)

      xmaxocc=real(ieffthr)/1000.0d0 
      key=" "
      if (icase.eq.1) then
       key="ALPHA"
      else if(icase.eq.2) then
        key="BETA "
      end if

      print *,' '
      print *,' ------------------------'
      print *,'  DOING  MULLIKEN EFFAO '
      print *,'   GENERAL FORMULATION   '
      print *,' ------------------------'
      print *,' '
      if(icase.ne.0) print *,'     ',key,'  PART  '

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
        kk=0
        do k=llim(icenter),imaxo           
         kk=kk+1
         do mu=1,igr                       
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

       return
       end

      subroutine ueffaolow2(p,icase)
      use basis_set
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /atlist/iatlist(maxat),icuat
      common/effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
!      dimension p(nmax,nmax)
      character*(5) key

      allocatable s0(:,:), sm(:,:), c0(:,:), splus(:,:), pp0(:,:)
      allocatable s0all(:),is0all(:)
      ndim=igr

      allocate (s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      allocate (c0(ndim,ndim), pp0(ndim,ndim), is0all(ndim))

      icube=Iopt(13)
      ieffthr  = Iopt(24)

      xmaxocc=real(ieffthr)/1000.0d0

      key=" "
      if (icase.eq.1) then
       key="ALPHA"
      else if(icase.eq.2) then
        key="BETA "
      end if

      print *,' '
      print *,' ------------------------'
      print *,'  DOING  LOWDIN    EFFAO '
      print *,'   GENERAL FORMULATION   '
      print *,' ------------------------'
      print *,' '
      if(icase.ne.0) print *,'     ',key,'  PART  '

      jocc=0
      xmaxotot=0.d0
      k1=0

        do i=1,igr
         do j=1,igr
          s0(i,j)=p(i,j)
          pp0(i,j)=s(i,j)
         end do
        end do
c make S1/2
      call build_Smp(igr,pp0,Sm,Splus,0)
c clean auxiliary pp0 mat
       do i=1,igr
        do j=1,igr
         pp0(i,j)=0.0d0
        end do
       end do

c tranform P with Splus
C will back trasnform effaos to ao basis later...then expanded in whole basis functions
      call to_lowdin_basis(igr,Splus,s0)
C 
c Block diagonal P
       do icenter=1,nat
        do mu=llim(icenter),iulim(icenter)
         do nu=llim(icenter),iulim(icenter)
          pp0(mu,nu)=s0(mu,nu)
         end do
        end do
       end do

       call diagonalize(igr,igr,pp0,C0,0)

c        call mprintnoat(pp0,igr,igr,igr,igr,'pp0')
c        call mprintnoat(c0,igr,igr,nmax,nmax,'c0 before')
c        call mprintnoat(c0,igr,igr,nmax,nmax,'c0 after')


c sort occup again for each atom
       k=0
       do icenter=1,nat
        do i=1,igr 
         kcenter=0
         xmax=0.0d0
         do jcenter=1,nat
          xx=0.0d0
          do mu=llim(jcenter),iulim(jcenter)
           xx=xx+abs(c0(mu,i))
          end do
          if(xx.gt.xmax) then
            xmax=xx
            kcenter=jcenter
          end if
         end do
         if(kcenter.eq.icenter) then
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
        end do
        end do
        
        do icenter=1,nat
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

c backtrasnform. eff.aos expanded over the full basis set
        call to_AO_basis(igr,igr,Sm,C0)

CCCCCC
C loop over atoms
CCCCCC
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

        kk=0
        do k=llim(icenter),imaxo           
         kk=kk+1
         do mu=1,igr                       
           p0(mu,kk)=c0(mu,k)
         end do
         p0net(kk,icenter)=pp0(k,k)
         p0gro(kk,icenter)=pp0(k,k)
        end do
        ip0(icenter)=imaxo-llim(icenter)+1

c OUTPUT ORBITALS FOR VISUALIZATION
c  write cube file
        if(icube.eq.1) call cubegen3(icenter,icase)

c end loop over atoms
        end do

        if(icase.eq.1) then
        open(68,file='eff-aos.low')
        call ival(68,"Number of alpha electrons",k1)
        call ival(68,"Number of basis functions",igr)
        call rarr(68,"Alpha MO occupations",k1,ndim,s0all)
        call rmat(68,"Alpha MO coefficients",k1,igr,nmax,p0)
        else
        call ival(68,"Number of beta electrons",k1)
        call rarr(68,"Beta MO occupations",k1,ndim,s0all)
        call rmat(68,"Beta MO coefficients",k1,igr,nmax,p0) 
        close(68)
        end if

       deallocate (s0, sm, c0,splus,pp0, s0all,is0all)

        return
      end

      subroutine ueffao3d(itotps,ndim,omp,chp,sat,wp,omp2,pk,icase)
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      integer,intent(in) :: itotps,ndim
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common /atlist/iatlist(maxat),icuat
      common/effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/qat(maxat,2),qsat(maxat,2)
      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat)
      dimension wp(itotps)
      dimension sat(ndim,ndim,nat)
      dimension pk(ndim,ndim)
      character*(5) key

      allocatable s0(:,:), sm(:,:), c0(:,:), splus(:,:), pp0(:,:)
      allocatable s0all(:),is0all(:)

      ihirsh = Iopt(6) 
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
       nat0=nat
       if(icuat.ne.nat) then
        nat0=icuat
       end if

      allocate (s0(ndim,ndim),s0all(ndim),sm(ndim,ndim),splus(ndim,ndim))
      allocate (c0(ndim,ndim), pp0(ndim,ndim), is0all(ndim))

        do iicenter=1,nat0
         icenter=iatlist(iicenter) 
c Computing Becke atomic NET orbital overlap
         do mu=1,ndim
          do nu=1,mu
            x=0.d0
            do jcenter=1,nat
            do ifut=iatps*(jcenter-1)+1,iatps*jcenter
             x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*
     +       omp2(ifut,icenter)*omp2(ifut,icenter)*omp(ifut)
            end do
            end do
            s0(mu,nu)=x
            s0(nu,mu)=x
          enddo
         enddo

c make SAA1/2
        call build_Smp(igr,S0,Sm,Splus,0)

c Save P matrix, only if not mulliken
        do i=1,igr
         do j=1,igr
         pp0(i,j)=pk(i,j)
         end do
        end do

c tranform block S with P0
       call to_lowdin_basis(igr,Splus,pp0)
       call diagonalize(igr,igr,pp0,C0,0)
       call to_AO_basis(igr,igr,Sm,C0)

C But max number of effaos in one center is igr 
        i=1
        do while(pp0(i,i).ge.xmaxocc.and.i.le.igr) 
         imaxo=i
         i=i+1
        end do
        xmaxo=0.0d0
        do i=1,igr  
         xmaxo=xmaxo+pp0(i,i)
        end do

        write(*,*) '                  '
        write(*,'(a,i4,a)') ' ** ATOM ',icenter,' ** '
        write(*,*) '                  '

        write(*,'(a25,i3,f10.5)') ' Net occupation for atom ',icenter, xmaxo
        if(icase.eq.0) then
        write(*,'(a30,f8.4)') ' Deviation from net population ', xmaxo -op(icenter,icenter)
        end if
        write(*,'(a23,f6.4)') ' Net occupation using >',xmaxocc
        write(*,60) (pp0(mu,mu),mu=1,imaxo)
60      format(7h OCCUP.   ,8F9.4)

c ...calculate gross occupations from orbitals!!!

        xx0=0.0d0
        do i=1,imaxo
         xx=0.0d0
         do j=1,igr                        
         do k=1,igr
          xx=xx+c0(k,i)*sat(k,j,icenter)*c0(j,i)
c          xx=xx+c0(k,i)*sp(k,j)*c0(j,i)
         end do
         end do
         xx=xx*pp0(i,i)
         s0all(i)=xx
         xx0=xx0+xx
        end do
         
        write(*,*) ' '
        write(*,'(a27,i3,f10.5)') ' Gross occupation for atom ',icenter, xx0   
        if(icase.eq.0) then
        write(*,'(a32,f8.4)') ' Deviation from gross population ', xx0 -qat(icenter,1)
        end if
        write(*,60) (s0all(mu),mu=1,imaxo)

C FORMATION OF HYBRIDS...renormalization

c        do k=1,imaxo           
c         do mu=1,igr                       
c           c0(mu,k)=c0(mu,k)/sqrt((s0(k,k)/2.0d0))
c          end do
c         end do

        do k=1,imaxo           
         do mu=1,igr                       
           p0(mu,k)=c0(mu,k)
         end do
         p0net(k,icenter)=pp0(k,k)
         p0gro(k,icenter)=s0all(k)
        end do
        ip0(icenter)=imaxo

c  write cube file
        if(icube.eq.1) call cubegen3(icenter,icase)
        
        end do
        deallocate (s0, sm, c0,splus,pp0, s0all,is0all)

         end

