!! **************************** !!
!! UEOS CALCULATION SUBROUTINES !!
!! **************************** !!

!! ***** !!

      subroutine effao3d_u(itotps,ndim,omp,chp,sat,wp,omp2)

!! THIS SUBROUTINE COMPUTES EFFAOs AND THEIR OCC. VALUES FROM THE PAIRED AND UNPAIRED DENSITIES !!
!! TAKATSUKA'S DEFINITION OF THE UNPAIRED DENSITY USED. HEAD-GORDON's IMPLEMENTATION IN DEVEL VERSION !!

      use basis_set
      use ao_matrices
      use integration_grid

      implicit real*8(a-h,o-z)
      include 'parameter.h'

      integer,intent(in) :: itotps,ndim

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

!! COMMON ADDED FOR EFFAO-U AND EOS-U, FIRST COLUMN FOR PAIRED AND SECOND FOR UNPAIRED !!
      common /ueffao/ up0(2,nmax,nmax),up0net(2,nmax,maxat),up0gro(2,nmax,maxat),iup0(2,maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /iops/iopt(100)

      character*80 line
      
      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat),wp(itotps)
      dimension sat(ndim,ndim,nat)

      dimension effrho(maxat),effrhoacc(maxat)

      allocatable :: Pno(:,:),Uno(:,:),UnoHG(:,:)
      allocatable :: SS0(:,:),S0(:,:),Sm(:,:),Splus(:,:),c0(:,:),pp0(:,:)
      allocatable :: d0(:,:),d1(:,:),scr(:),s0all(:)

      icube = iopt(13) 
      iatps = nang*nrad

!! THRESH TO INCLUDE EFFAO OR NOT FOR TOTAL NET POPULATION CALCULATION !!
      xminocc=1.0d-4

      write(*,*) " ------------------------------ "
      write(*,*) " DOING EFFAO-3D FROM U FUNCTION "
      write(*,*) " ------------------------------ "
      write(*,*) " "

!! EVALUATING U and P (AO BASIS) FROM NOs !!
      nnorb=0
      do ii=1,igr
        if(occ_no(ii,ii).ge.thresh) nnorb=nnorb+1
      end do
      write(*,'(2x,a33,x,i3)') "Number of natural orbitals (NOs):",nnorb
      write(*,*) " "

      ALLOCATE(Pno(igr,igr),Uno(igr,igr))
      do mu=1,igr
        do nu=1,igr
          xx=ZERO
          xx2=ZERO
          do ii=1,nnorb
            xx=xx+occ_no(ii,ii)*c_no(mu,ii)*c_no(nu,ii)
            xx2=xx2+occ_no(ii,ii)*(TWO-occ_no(ii,ii))*c_no(mu,ii)*c_no(nu,ii)
          end do
          Pno(mu,nu)=xx
          Uno(mu,nu)=xx2 !! TAKATSUKA's DEFINITION !!
        end do
      end do

!! EVALUATING EFFAOs !!
      ALLOCATE(scr(itotps))
      ALLOCATE(S0(igr,igr),s0all(igr),Sm(igr,igr),Splus(igr,igr))
      ALLOCATE(SS0(igr,igr),c0(igr,igr),pp0(igr,igr))

!! icase = 1 PAIRED EFOs, icase = 2 UNPAIRED EFOs !!
      do icase=1,2
        if(icase.eq.1) then
          write(*,*) " ------------------------------ "
          write(*,*) " EFFAOs FROM THE PAIRED DENSITY "
          write(*,*) " ------------------------------ "
          write(*,*) " "
        else if(icase.eq.2) then
          write(*,*) " -------------------------------- "
          write(*,*) " EFFAOs FROM THE UNPAIRED DENSITY "
          write(*,*) " -------------------------------- "
          write(*,*) " "
        end if
        do iicenter=1,icufr

!! W_A (A = fragA) EVALUATION !!
          scr=ZERO
          do ifut=1,itotps  
            do icenter=1,nfrlist(iicenter)
              scr(ifut)=scr(ifut)+omp2(ifut,ifrlist(icenter,iicenter))
            end do
          end do

!! SAA1/2 !!
          do mu=1,igr
            do nu=1,mu

!! ALLPOINTS INTEGRATION !!
              xx=ZERO
              do jcenter=1,nat
                do ifut=iatps*(jcenter-1)+1,iatps*jcenter
                  xx=xx+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*scr(ifut)*scr(ifut)*omp(ifut)
                end do
              end do
              S0(mu,nu)=xx
              S0(nu,mu)=xx
            end do
          end do
          call build_Smp(igr,S0,Sm,Splus,0)

!! TAKATSUKA DEFINITION !!
          do ii=1,igr
            do jj=1,igr
              if(icase.eq.1) pp0(ii,jj)=Pno(ii,jj)-Uno(ii,jj)
              if(icase.eq.2) pp0(ii,jj)=Uno(ii,jj)
            end do
          end do

          call to_lowdin_basis(igr,Splus,pp0)
          call diagonalize(igr,igr,pp0,c0,0)
          call to_AO_basis(igr,igr,Sm,c0)

!! TRUNCATING MAX NUMBER OF EFOs WITH xminocc !!
          ii=1
          do while(pp0(ii,ii).ge.xminocc.and.ii.le.igr) 
            imaxo=ii
            ii=ii+1
          end do
          xmaxo=ZERO
          do ii=1,igr  
            xmaxo=xmaxo+pp0(ii,ii)
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
          write(*,'(2x,a27,x,i3,x,f10.5)') "Net occupation for fragment",iicenter,xmaxo
          write(*,'(2x,a22,x,f10.5)') "Net occupation using >",xminocc
          write(*,60) (pp0(mu,mu),mu=1,imaxo)
          write(*,*) " "

!! COMPUTING GROSS OCC. !!
          xx0=ZERO
          do ii=1,imaxo
            xxx=ZERO
            do icenter=1,nfrlist(iicenter)
              jcenter=ifrlist(icenter,iicenter)
              xx=ZERO 
              do jj=1,igr                        
                do kk=1,igr
                  xx=xx+c0(kk,ii)*sat(kk,jj,jcenter)*c0(jj,ii)
                end do
              end do
              xxx=xxx+xx
            end do
            xxx=xxx*pp0(ii,ii)
            s0all(ii)=xxx
            xx0=xx0+xxx
          end do
          write(*,'(a31,i3,f10.5)') " Gross occupation for fragment ",iicenter,xx0
          write(*,60) (s0all(mu),mu=1,imaxo)
          write(*,*) " "

!! SAVING EFOs AND OCCUPANCIES TO MAKE CUBES AND ASSIGN OSs !!
          do kk=1,imaxo
            do mu=1,igr
              up0(icase,mu,kk)=c0(mu,kk)
            end do
            up0net(icase,kk,iicenter)=pp0(kk,kk)
            up0gro(icase,kk,iicenter)=s0all(kk)
          end do
          iup0(icase,iicenter)=imaxo

!! WRITE CUBEFILE !! (TO DO)
          if(icase2.eq.0) iicase=3
          if(icase2.eq.1) iicase=4
          if(icube.eq.1) call cubegen4(iicenter,iicase)
        end do
      end do

!! DEALLOCATING MATRICES !!
      DEALLOCATE(SS0,S0,Splus,Sm)
      DEALLOCATE(scr,c0,pp0,s0all)
      DEALLOCATE(Pno,Uno)

!! PRINTING FORMATS !!
60    FORMAT(7h OCCUP.   ,8F9.4)

      end

!! ****** !!

      subroutine ueos_analysis(xthresh)

!! SUBROUTINE TO ASSIGN PAIRED AND UNPAIRED ELECTRONS TO FRAGMENTS, RESULTING INTO FORMAL OSs !!

      implicit real*8(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /atlist/iatlist(maxat),icuat
      common /loba/ oxi(maxat),errsav(maxat),elec(maxat),effpop(maxat)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /loba2/occup(nmax,2),iorbat(nmax,2),lorb(2),confi0
      
      common/effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)

!! COMMON ADDED FOR EFFAO-U AND EOS-U, FIRST COLUMN FOR PAIRED AND SECOND FOR UNPAIRED !!
      common /ueffao/ up0(2,nmax,nmax),up0net(2,nmax,maxat),up0gro(2,nmax,maxat),iup0(2,maxat)
      common /iops/iopt(100)

      dimension occup2(igr),occupg(igr),errnet(maxat)
      dimension qdev(maxat)


      if(idobeta.eq.0.and.icase.eq.2) then
       write(*,*)'  '
       write(*,*)' SKIPPING EFFAOs FOR BETA ELECTRONS'
       write(*,*)'  '
       lorb(2)=lorb(1)
       do i=1,lorb(2)
         occup(i,2)=occup(i,1) 
         iorbat(i,2)= iorbat(i,1)
       end do  
       do i=1,icufr
        errnet(i)=errsav(i)
        effpop(i)=effpop(i)*2.0d0
       end do
       go to 99
      else if(icase.eq.2.and.nb.eq.0) then
       write(*,*)'  '
       write(*,*)' CALCULATION HAS NO BETA ELECTRONS'
       write(*,*)'  '
       do i=1,icufr
        elec(i)=0.0d0
       end do
       go to 99
      end if

      iorb=0
      do i=1,icufr
       elec(i)=0.0d0
       do k=1,ip0(i)
        iorb=iorb+1
        occup(iorb,icase)=p0net(k,i)
        iorbat(iorb,icase)=i
       end do 
      end do  
      write(*,*)'  '
      write(*,*) ' Total number of eff-AO-s for analysis: ',iorb
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
       write(*,*) ' EOS: Unambiguous integer electron assignation'
       do i=1,iorb
        if(i.le.nnn) then
        occup2(i)=1.0d0
        else
        occup2(i)=0.0d0
        end if
       end do 
      else
       write(*,*) ' ******************************************'
       write(*,*) ' EOS: Warning, pseudo-degeneracies detected'
       write(*,*) ' ******************************************'
       kk=0
334    kk=kk+1
       if(dabs(occup(nnn,icase)-occup(nnn-kk,icase)).lt.thres) go to 334
       kk=kk-1
       frac=float(kk+1)/float(kk+k+1)
       write(*,'(a14,i4,a16,i4,a32)') ' Distributing ',kk+1,' electrons 
     +over ',kk+k+1,' pseudodegenerate atomic orbitals'
       do i=1,iorb
        if(i.lt.nnn-kk) then
         occup2(i)=1.0d0
        else if(i.gt.nnn+k) then
         occup2(i)=0.0d0
        else
         occup2(i)=frac 
        end if
       end do 
      end if
     
      do i=1,iorb
       elec(iorbat(i,icase))=elec(iorbat(i,icase))+occup2(i)
      end do

      if(icase.eq.1) then
       write(*,*)' EOS ANALYSIS  FOR ALPHA ELECTRONS'
       do i=1,icuat
        errsav(i)=errnet(i)
        effpop(i)=0.0d0
       end do
      else if (icase.eq.2) then
       write(*,*)' EOS ANALYSIS  FOR BETA ELECTRONS'
      end if

      print *,'  Fragm   Elect  Last occ. First unocc '
      print *,' -------------------------------------'

       xlast=1.0d0
       ilast=0
       do i=1,icufr
        nn=int(elec(i))
        if(elec(i)-nn.gt.thresh) nn=nn+1
        if(nn+1.gt.ip0(i)) then 
         write(*,10) i,elec(i), p0net(nn,i)
        else
         if(nn.ne.0) then
          write(*,15) i,elec(i), p0net(nn,i),p0net(nn+1,i)
         else
          write(*,15) i,elec(i), 0.0d0 , p0net(nn+1,i)
         end if
        end if
        if(nn.ne.0) then
         if(p0net(nn,i).lt.xlast) then
          xlast= p0net(nn,i)
          ilast=i
         end if
        end if
        do j=1,nn
         effpop(i)= effpop(i) + p0net(j,i)
        end do
       end do 

       xfirst=0.0d0
       do i=1,icufr
        if(i.ne.ilast) then
        nn=int(elec(i))
        if(elec(i)-nn.gt.thresh) nn=nn+1
        if(nn+1.ne.0.and.p0net(nn+1,i).gt.xfirst) xfirst= p0net(nn+1,i)
        end if
       end do 

      print *,' -------------------------------------'
      confi=100.0*min(1.0d0,xlast-xfirst+0.5d0)
      write(*,'(a30,f8.3)') 'RELIABILITY INDEX R(%) =',confi


c PSS quadratic deviation of the fractional and integer occupations as new fragment indicator
      print *
      print *,'(warning in case of fractional occupations)'
      print *,'  Fragm    Quad Dev '
      print *,' -------------------'
      do i=1,icufr
       nn=int(elec(i))
       if(elec(i)-nn.gt.thresh) nn=nn+1
       qdev(i)=0.0d0
       do k=1,ip0(i)
        xx=p0gro(k,i)
        if(k.le.nn) xx=1.0d0-p0gro(k,i)
        qdev(i)=qdev(i)+xx*xx
c        write(*,*) 'pollas',i,ip0(i),p0gro(k,i),xx
       end do
       qdev(i)=sqrt(qdev(i))
       write(*,'(I4,5x,f8.4)') i, qdev(i)
      end do
      print *,' -------------------'
C PSS 
      if(icase.eq.1) confi0=confi

99     continue
       zztot=0.0d0
       do i=1,icufr
        if(icase.eq.1) then
         zzn=0.0d0
         do j=1,nfrlist(i)
          zzn=zzn+zn(ifrlist(j,i))
         end do
         oxi(i)=zzn-elec(i)
        else
         oxi(i)=oxi(i)-elec(i)
         zztot=zztot+oxi(i)
        end if
       end do

       if(icase.eq.2) then
      print *,'  '
      print *,'  '
      print *,'  FRAGMENT  OXIDATION STATES '
      print *,'  '
      print *,'  Frag   Oxidation State  '
      print *,' -------------------------'
       do i=1,icufr
        write(*,20) i,oxi(i)
       end do 
      print *,' -------------------------'
      write(*,'(a7,f5.1)') '   Sum ', zztot 
      print *,'  '
      confi=min(confi,confi0)
      write(*,'(a39,f8.3)') ' OVERALL RELIABILITY INDEX R(%) =',confi

c Ignore this part for now
      if(1.eq.0) then
      print *,'  '
      print *,'  '
      print *,'  FRAGMENT  EFFECTIVE CHARGES AND POPULATIONS '
      print *,'  '
      print *,'  Frag     Charge     Population  '  
      print *,' -------------------------------'
      zztot=0.0d0
      do i=1,icufr
       zzn=0.0d0
       do j=1,nfrlist(i)
        zzn=zzn+zn(ifrlist(j,i))
       end do
       write(*,25) i,zzn-effpop(i),effpop(i)
      end do 
      print *,' -------------------------'
      end if
c

      end if

      do i=1,icufr
       zzn=0.0d0
       do j=1,nfrlist(i)
       end do
      end do 

10    format(i4,2x,f6.2,f12.3,'    < thresh',2f12.3)
15    format(i4,2x,f6.2,2f12.3)
20    format(i4,2x,f10.2)
25    format(i4,2x,f10.4,2x,f10.4)

      end 

!! ****** !!
