
!! ***** !!

      subroutine uf_efo3d(itotps,ndim,omp,chp,sat,wp,omp2,icase,icase2)
      use basis_set
      use ao_matrices
      use integration_grid

!! THIS SUBROUTINE COMPUTES THE EFFAOs FROM THE EFFECTIVE NUMBER !!
!! OF UNPAIRED ELECTRONS !!

      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'

      integer,intent(in) :: itotps,ndim

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop

!     common /localspin/xlsa(maxat,maxat),ua(maxat)

      common /iops/iopt(100)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /ovpop/op(maxat,maxat),bo(maxat,maxat),di(maxat,maxat),totq
      common /qat/qat(maxat,2),qsat(maxat,2)

      character*(5) key
      character*80 line

      dimension chp(itotps,ndim),omp(itotps),omp2(itotps,nat),wp(itotps)
      dimension sat(ndim,ndim,nat)

      dimension effrho(maxat),effrhoacc(maxat)

      allocatable :: Pno(:,:),Uno(:,:),UnoHG(:,:)

      allocatable :: SS0(:,:),S0(:,:),Sm(:,:),Splus(:,:),c0(:,:),pp0(:,:)
      allocatable :: d0(:,:),d1(:,:),scr(:),s0all(:)

!     ieffao  = Iopt(12) 
      icube   = Iopt(13) 
      iatps   = nang*nrad

      xminocc=1.0d-3

      write(*,*) " "
      write(*,*) " -----------------------------------"
      write(*,*) "   DOING EFFAO-3D FROM U FUNCTION   "
      write(*,*) " -----------------------------------"
      write(*,*) " "

!! EVALUATING U_munu and P_munu FROM NOs !!

      nnorb=0
      do ii=1,igr
        if(occ_no(ii,ii).ge.thresh) nnorb=nnorb+1
      end do
      write(*,*) " "
      write(*,*) "Number of natural orbitals over thresh ",nnorb
      write(*,*) " "

      ALLOCATE(Pno(igr,igr),Uno(igr,igr),UnoHG(igr,igr))
      do mu=1,igr
        do nu=1,igr
          xx=ZERO
          xx2=ZERO
          xx3=ZERO
          do ii=1,nnorb
            xx=xx+occ_no(ii,ii)*c_no(mu,ii)*c_no(nu,ii)
            xx2=xx2+occ_no(ii,ii)*(TWO-occ_no(ii,ii))*c_no(mu,ii)*c_no(nu,ii)
            xx3=xx3+(ONE-ABS(ONE-occ_no(ii,ii)))*c_no(mu,ii)*c_no(nu,ii)
          end do
          Pno(mu,nu)=xx
          Uno(mu,nu)=xx2
          UnoHG(mu,nu)=xx3
        end do
      end do

!! EOS-U PROCEDURE STARTS HERE !!

      ALLOCATE(scr(itotps))
      ALLOCATE(S0(igr,igr),s0all(igr),Sm(igr,igr),Splus(igr,igr))
      ALLOCATE(SS0(igr,igr),c0(igr,igr),pp0(igr,igr))
      do iicenter=1,icufr

!! W_A (A = fragA) evaluation !!

        do ifut=1,itotps
          scr(ifut)=ZERO
          do icenter=1,nfrlist(iicenter)
            scr(ifut)=scr(ifut)+omp2(ifut,ifrlist(icenter,iicenter))
          end do
        end do

!! SAA1/2 !!

        do mu=1,igr
          do nu=1,mu

!! ALLPOINTS INTEGRATION !!

            x=ZERO
            do jcenter=1,nat
              do ifut=iatps*(jcenter-1)+1,iatps*jcenter
                x=x+wp(ifut)*chp(ifut,mu)*chp(ifut,nu)*scr(ifut)*scr(ifut)*omp(ifut)
              end do
            end do
            S0(mu,nu)=x
            S0(nu,mu)=x
          end do
        end do
        call build_Smp(igr,S0,Sm,Splus,0)

!! TAKATSUKA DEFINITION !!

        if(icase.eq.0) then
          write(*,*) " "
          write(*,*) " TAKATSUKA DEFINITION "
          write(*,*) " "
          do ii=1,igr
            do jj=1,igr
              if(icase2.eq.0) pp0(ii,jj)=Pno(ii,jj)-Uno(ii,jj)
              if(icase2.eq.1) pp0(ii,jj)=Uno(ii,jj)
            end do
          end do

!! HEAD-GORDON DEFINITION !!

        else if(icase.eq.1) then
          write(*,*) " "
          write(*,*) " HEAD-GORDON DEFINITION "
          write(*,*) " "
          do ii=1,igr
            do jj=1,igr
              if(icase2.eq.0) pp0(ii,jj)=Pno(ii,jj)-UnoHG(ii,jj)
              if(icase2.eq.1) pp0(ii,jj)=UnoHG(ii,jj)
            end do
          end do
        end if

!! EVALUATING EFOs !!

        call to_lowdin_basis(igr,Splus,pp0)
        call diagonalize(igr,igr,pp0,c0,0)
        call to_AO_basis(igr,igr,Sm,c0)

!! MAX NUMBER OF EFOs IS igr, BUT TRUNCATED WITH A 10^-4 THRESH !!

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
        write(*,*) " "
        write(*,'(a,i4,a)') " ** FRAGMENT ",iicenter," ** "
        write(*,*) " "

        write(*,'(a29,i3,f10.5)') " Net occupation for fragment ",iicenter,xmaxo
        write(*,'(a23,f6.4)') " Net occupation using >",xminocc
        write(*,60) (pp0(mu,mu),mu=1,imaxo)

!! ...calculate gross occupations from orbitals !!

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
        
        write(*,*) " "
        write(*,'(a31,i3,f10.5)') " Gross occupation for fragment ",iicenter,xx0
        write(*,60) (s0all(mu),mu=1,imaxo)

        do kk=1,imaxo           
          do mu=1,igr                       
            p0(mu,kk)=c0(mu,kk)
          end do
          p0net(kk,iicenter)=pp0(kk,kk)
          p0gro(kk,iicenter)=s0all(kk)
        end do
        ip0(iicenter)=imaxo

!! WRITE CUBEFILE !!

        if(icase2.eq.0) iicase=3
        if(icase2.eq.1) iicase=4
        if(icube.eq.1) call cubegen4(iicenter,iicase)

!! PROVISIONAL PRINTING !!

        write(*,*) " "
        write(*,*) " EFO POPULATIONS FROM FRAGMENT ",iicenter
        write(*,*) " "
        if(icase2.eq.0) write(*,*) " !! PAIRED (DIVIDED BY TWO) !! "
        if(icase2.eq.1) write(*,*) " !! UNPAIRED !! "
        write(*,*) " "
        do ii=1,imaxo
          if(icase2.eq.0) x=pp0(ii,ii)/TWO
          if(icase2.eq.1) x=pp0(ii,ii)
          if(x.ge.xminocc) write(*,*) "EFO, OCC : ",ii,x
        end do
        write(*,*) " "

      end do

      DEALLOCATE(SS0,S0,Splus,Sm)
      DEALLOCATE(scr,c0,pp0,s0all)
      DEALLOCATE(Pno,Uno,UnoHG)

60    FORMAT(7h OCCUP.   ,8F9.4)

      end


!! ***** !!

