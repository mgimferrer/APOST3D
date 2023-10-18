!! **************************** !!
!! UEOS CALCULATION SUBROUTINES !!
!! **************************** !!

!! ***** !!

      subroutine effao3d_u(itotps,ndim,omp,chp,sat,wp,omp2,iueos)

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
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common/effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /iops/iopt(100)

      character*80 line
      
      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat),wp(itotps)
      dimension sat(ndim,ndim,nat)

      dimension effrho(maxat),effrhoacc(maxat)

      allocatable :: Pno(:,:),Uno(:,:),UnoHG(:,:)
      allocatable :: SS0(:,:),S0(:,:),Sm(:,:),Splus(:,:),c0(:,:),pp0(:,:)
      allocatable :: d0(:,:),d1(:,:),scr(:),s0all(:)

      allocatable :: iup0(:,:),up0net(:,:,:),up0gro(:,:,:)

      icube = iopt(13) 
      iatps = nang*nrad

!! THRESH TO INCLUDE EFFAO OR NOT FOR TOTAL NET POPULATION CALCULATION !!
      xminocc=1.0d-4

      write(*,*) " "
      write(*,*) " -------------------------------- "
      write(*,*) "  DOING EFFAO-3D FROM U FUNCTION  "
      write(*,*) " -------------------------------- "
      write(*,*) " "
      write(*,*) " EFFAO-U: paired and unpaired densities treated separately "

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

!! ALLOCATING MATRICES FOR EFFAO-U EVALUATION AND OSs IF REQUESTED !!
      ALLOCATE(scr(itotps))
      ALLOCATE(S0(igr,igr),s0all(igr),Sm(igr,igr),Splus(igr,igr))
      ALLOCATE(SS0(igr,igr),c0(igr,igr),pp0(igr,igr))
      if(iueos.eq.1) ALLOCATE(iup0(2,icufr),up0net(2,igr,icufr), up0gro(2,igr,icufr))

!! EVALUATING EFFAOs: icase = 1 PAIRED EFOs, icase = 2 UNPAIRED EFOs !!
      do icase=1,2
        if(icase.eq.1) then
          write(*,*) " -------------------------------- "
          write(*,*) "  EFFAOs FROM THE PAIRED DENSITY  "
          write(*,*) " -------------------------------- "
          write(*,*) " "
        else if(icase.eq.2) then
          write(*,*) " ---------------------------------- "
          write(*,*) "  EFFAOs FROM THE UNPAIRED DENSITY  "
          write(*,*) " ---------------------------------- "
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
          write(*,'(2x,a29,x,i3,f10.5)') "Gross occupation for fragment",iicenter,xx0
          write(*,60) (s0all(mu),mu=1,imaxo)
          write(*,*) " "

!! FOR CUBE CREATION AND TO ASSIGN OSs IF REQUESTED !!
          do kk=1,imaxo
            do mu=1,igr
              p0(mu,kk)=c0(mu,kk)
            end do
            p0net(kk,iicenter)=pp0(kk,kk)
            p0gro(kk,iicenter)=s0all(kk)
            if(iueos.eq.1) then
              up0net(icase,kk,iicenter)=pp0(kk,kk)
              up0gro(icase,kk,iicenter)=s0all(kk)
            end if
          end do
          ip0(iicenter)=imaxo
          if(iueos.eq.1) iup0(icase,iicenter)=imaxo

!! CUBE PRINTING !!
          if(icase.eq.1) iicase=3
          if(icase.eq.2) iicase=4
          if(icube.eq.1) call cubegen4(iicenter,iicase)
        end do
      end do

!! EXTRACTING OSs IF REQUESTED !!
      if(iueos.eq.1) call ueos_analysis(iup0,up0gro) !! USING GROSS POPULATIONS FOR ASSIGNMENT !!

!! DEALLOCATING MATRICES !!
      DEALLOCATE(SS0,S0,Splus,Sm)
      DEALLOCATE(scr,c0,pp0,s0all)
      DEALLOCATE(Pno,Uno)
      DEALLOCATE(iup0,up0net,up0gro)

!! PRINTING FORMATS !!
60    FORMAT(8h  OCCUP. ,8f9.4)

      end

!! ****** !!

      subroutine ueos_analysis(iup0,up0gro)

!! SUBROUTINE TO ASSIGN PAIRED AND UNPAIRED ELECTRONS TO FRAGMENTS, RESULTING INTO FORMAL OSs !!

      implicit real*8(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /loba/ oxi(maxat),errsav(maxat),elec(maxat),effpop(maxat)
      common /loba2/occup(nmax,2),iorbat(nmax,2),lorb(2),confi0

      dimension iup0(2,icufr),up0gro(2,igr,icufr)

      allocatable :: infoelec(:,:)

      ALLOCATE(infoelec(maxat,2))

!! ZEROING ELECTRON COUNTING MATRICES !!
      elec=ZERO
      infoelec=0 !! STORES NUMBER OF PAIR AND UNPAIR ELECTRONS ASSIGNED TO EACH FRAGMENT !!

!! ALL EFOs IN THE SAME MATRIX !!
      do icase=1,2 !! 1 = PAIRED, 2 = UNPAIRED !!
        iorb=0
        do ii=1,icufr
          do kk=1,iup0(icase,ii)
            iorb=iorb+1
            occup(iorb,icase)=up0gro(icase,kk,ii)
            iorbat(iorb,icase)=ii
          end do 
        end do
        lorb(icase)=iorb
      end do
      write(*,'(2x,a38,x,i4)') "Total number of eff-AO-s for analysis:",lorb(1)+lorb(2)

!! ASSIGNING ELECTRONS: NO FRACTIONARY ASSIGNMENT ALLOWED IN THE CODE (YET), WILL PRINT R(%) = 50 !!
      nnn=nalf+nb
      do while(nnn.gt.0)

!! EVEN NUMBER OF ELECTRONS LEFT TO ASSIGN: BOTH PAIRED AND UNPAIRED COMPETE !!
        if(mod(nnn,2).eq.0) then

!! FIRST PAIRED EFOs !!
          imaxocc=0
          xmaxocc=ZERO
          do ii=1,lorb(1)
            xx=occup(ii,1)

!! SAVING LARGEST OCC. !! 
            if(xx.gt.xmaxocc) then
              xmaxocc=xx
              imaxocc=ii
            end if
          end do

!! NOW UNPAIRED EFOs: SAVING TWO LARGEST OCC. !!
          imaxocc2=0
          xmaxocc2=ZERO
          do ii=1,lorb(2)
            xx=occup(ii,2) 
            if(xx.gt.xmaxocc2) then
              xmaxocc2=xx
              imaxocc2=ii
            end if
          end do

!! TRICK FOR IF UNPAIRED IS "VERY UNPAIRED" (HAS TO PLAY ALONE) WITH EVEN NUMBER OF EL. TO ASSIGN !!
          if(occup(imaxocc2,2).gt.0.75d0) then

!! ASSIGNING A SINGLE ELECTRON !!
            elec(iorbat(imaxocc2,2))=elec(iorbat(imaxocc2,2))+ONE
            infoelec(iorbat(imaxocc2,2),2)=infoelec(iorbat(imaxocc2,2),2)+1
            occup(imaxocc2,2)=ZERO !! ZEROING FOR NOT BEING INVOLVED IN NEXT ITERATION !!
            nnn=nnn-1

!! EVALUATING IF EL. PAIR SPLITTING IS REQUIRED !!
          else

!! NOW SECOND ONE !!
            imaxocc3=0
            xmaxocc3=ZERO
            do ii=1,lorb(2)

!! ONLY IF THEY ARE OF DIFFERENT FRAGMENT !!
              if(iorbat(imaxocc2,2).ne.iorbat(ii,2)) then
                xx=occup(ii,2)
                if(xx.gt.xmaxocc3.and.ii.ne.imaxocc2) then
                  imaxocc3=ii
                  xmaxocc3=xx
                end if
              end if
            end do
            if(xmaxocc.gt.(xmaxocc2+xmaxocc3)) then

!! ASSIGNING A PAIR OF ELECTRONS !!
              elec(iorbat(imaxocc,1))=elec(iorbat(imaxocc,1))+TWO
              infoelec(iorbat(imaxocc,1),1)=infoelec(iorbat(imaxocc,1),1)+2
              occup(imaxocc,1)=ZERO !! ZEROING FOR NOT BEING INVOLVED IN NEXT ITERATION !!
            else

!! ASSIGNING ONE ELECTRON TO EACH FRAGMENT !!
              elec(iorbat(imaxocc2,2))=elec(iorbat(imaxocc2,2))+ONE
              elec(iorbat(imaxocc3,2))=elec(iorbat(imaxocc3,2))+ONE
              infoelec(iorbat(imaxocc2,2),2)=infoelec(iorbat(imaxocc2,2),2)+1
              infoelec(iorbat(imaxocc3,2),2)=infoelec(iorbat(imaxocc3,2),2)+1

!! ZEROING FOR NOT BEING INVOLVED IN NEXT ITERATION !!
              occup(imaxocc2,2)=ZERO
              occup(imaxocc3,2)=ZERO
            end if
            nnn=nnn-2
          end if

!! ODD NUMBER OF ELECTRONS LEFT TO ASSIGN: ONLY UNPAIRED !!
        else
          imaxocc=0
          xmaxocc=ZERO
          do ii=1,lorb(2)
            xx=occup(ii,2)

!! SAVING FOR BOTH THE LARGEST AND THE SECOND LARGEST !! 
            if(xx.gt.xmaxocc) then
              xmaxocc=xx
              imaxocc=ii
            end if
          end do

!! ASSIGNING A SINGLE ELECTRON !!
          elec(iorbat(imaxocc,2))=elec(iorbat(imaxocc,2))+ONE
          infoelec(iorbat(imaxocc,2),2)=infoelec(iorbat(imaxocc,2),2)+1
          occup(imaxocc,2)=ZERO !! ZEROING FOR NOT BEING INVOLVED IN NEXT ITERATION !!
          nnn=nnn-1
        end if
      end do

!! PRINTING INFO !!
      do icase=1,2
        write(*,*) " "
        if(icase.eq.1) then
          write(*,*) " ----------------------------------- "
          write(*,*) "  EOS ANALYSIS FOR PAIRED ELECTRONS  "
          write(*,*) " ----------------------------------- "
        else if(icase.eq.2) then
          write(*,*) " ------------------------------------- "
          write(*,*) "  EOS ANALYSIS FOR UNPAIRED ELECTRONS  "
          write(*,*) " ------------------------------------- "
        end if
        write(*,*) " "
        write(*,*) "  Frag.  Elect.  Last occ.  First unocc.  "
        write(*,*) " ---------------------------------------- "

!! EVALUATING FIRST AND LAST OCC. FROM EACH FRAGMENT !!

        xlast=TWO !! PAIRED GROSS POPULATIONS NOT DIVIDED BY TWO NOW !!
        ilast=0
        do ifrg=1,icufr
          if(icase.eq.1) nn=infoelec(ifrg,icase)/2
          if(icase.eq.2) nn=infoelec(ifrg,icase)
          if(nn+1.gt.iup0(icase,ifrg)) then 
            write(*,10) ifrg,REAL(infoelec(ifrg,icase)),up0gro(icase,nn,ifrg)
          else
            if(nn.ne.0) then
              write(*,15) ifrg,REAL(infoelec(ifrg,icase)),up0gro(icase,nn,ifrg),up0gro(icase,nn+1,ifrg)
            else
              write(*,15) ifrg,REAL(infoelec(ifrg,icase)),ZERO,up0gro(icase,nn+1,ifrg)
            end if
          end if
          if(nn.ne.0) then
            if(up0gro(icase,nn,ifrg).lt.xlast) then
              xlast=up0gro(icase,nn,ifrg)
              ilast=ifrg
            end if
          end if
        end do

        xfirst=ZERO
        do ifrg=1,icufr
          if(ifrg.ne.ilast) then
            if(icase.eq.1) nn=infoelec(ifrg,icase)/2
            if(icase.eq.2) nn=infoelec(ifrg,icase)
            if(nn+1.ne.0.and.up0gro(icase,nn+1,ifrg).gt.xfirst) xfirst=up0gro(icase,nn+1,ifrg)
          end if
        end do 
        write(*,*) " ---------------------------------------- "

!! NOW PAIRED AND UNPAIRED ARE DIFFERENT (FACTOR OF 2 IN POPULATIONS) !!
        if(icase.eq.1) confi=100.0*min(1.0d0,(xlast-xfirst)/TWO+0.5d0)
        if(icase.eq.2) confi=100.0*min(1.0d0,xlast-xfirst+0.5d0)
        write(*,'(3x,a24,x,f7.3)') "RELIABILITY INDEX R(%) =",confi
        if(icase.eq.1) confi0=confi
      end do

!! EXTRACTING OSs !!
      zztot=ZERO
      do ifrg=1,icufr
        zzn=ZERO
        do jfrg=1,nfrlist(ifrg)
          zzn=zzn+zn(ifrlist(jfrg,ifrg))
        end do
        oxi(ifrg)=zzn-elec(ifrg)
        zztot=zztot+oxi(ifrg)
      end do

!! FINAL PRINTING !!
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
      confi2=min(confi,confi0)
      write(*,'(2x,a32,x,f7.3)') "OVERALL RELIABILITY INDEX R(%) =",confi2

!! PRINTING FORMATS !!
10    FORMAT(3x,i3,3x,f6.2,4x,f12.3,'    < thresh',2f12.3) !tocheck!
15    FORMAT(3x,i3,3x,f6.2,4x,f6.3,4x,f6.3)
20    FORMAT(3x,i3,6x,f8.2)

      DEALLOCATE(infoelec)

      end 

!! ****** !!
