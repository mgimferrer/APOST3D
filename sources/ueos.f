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

      allocatable :: tmp_occup(:),tmp_iorbat(:)
      allocatable :: elec2(:,:),elec_id(:,:)
      allocatable :: tmp_elec_id(:,:)
      allocatable :: elec_frg_count(:,:)


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

      !! SORTING BY DECREASING OCC VALUE !!
      do icase=1,2
        write(*,*) " "
        write(*,*) "icase: ",icase

        !! 1) COPYING IN TEMPORARY VECTORS !!
        ALLOCATE(tmp_occup(lorb(icase)))
        ALLOCATE(tmp_iorbat(lorb(icase)))
        do ii=1,lorb(icase)
          tmp_occup(ii)=occup(ii,icase)
          tmp_iorbat(ii)=iorbat(ii,icase)
        end do

        !! 2) ORDERING PAIRED AND UNPAIRED EFOs BY DECREASING OCCUPATION VALUE !!
        do ii=1,lorb(icase)-1
          do jj=ii+1,lorb(icase)
            if(tmp_occup(ii).lt.tmp_occup(jj)) then

              !! SWAPPING OCC. VALUES !!
              tmp_swap=tmp_occup(ii)
              tmp_occup(ii)=tmp_occup(jj)
              tmp_occup(jj)=tmp_swap
        
              !! ALSO REORDERING FRAG. INFO FOR CONSISTENCY !!
              tmp_swap=tmp_iorbat(ii)
              tmp_iorbat(ii)=tmp_iorbat(jj)
              tmp_iorbat(jj)=tmp_swap
            end if
          end do
        end do
    
        !! UPDATING THE MAIN ARRAYS AFTER SORTING !!
        do ii=1,lorb(icase)
          occup(ii,icase)=tmp_occup(ii)
          iorbat(ii,icase)=tmp_iorbat(ii)
        end do
        DEALLOCATE(tmp_occup,tmp_iorbat)

        !! checking
        write(*,*) "occ,frg"
        do ii=1,lorb(icase)
          write(*,*) occup(ii,icase),iorbat(ii,icase)
        end do

      end do

      !! CREATING MATRIX OF IDEAL OCCUPATIONS !!
      ALLOCATE(elec_id(lorb(1),2)) !! ALLOCATED TO lorb(1) FOR SIMPLICITY !!
      ALLOCATE(elec2(lorb(1),2)) !! SAME HERE !!
      elec_id=ZERO
      elec2=occup

      !! ASSIGNING ELECTRONS !!
      nnn=nalf+nb

      !! UNPAIRED IN CASE OF ALPHA =/ BETA !!
      nunp=nalf-nb
      write(*,*) "Number of for sure unpaired electrons: ",nunp
      if(nunp.ne.0) then
        do ii=1,nunp
          elec_id(ii,2)=ONE
          nnn=nnn-1
        end do
      end if

      npair=nnn/2
      !! check
      write(*,*) "NUMBER OF EL. PAIRS TO ASSIGN: ",npair

      !! PAIRED FOR THE REST !!
      do ii=1,npair
        elec_id(ii,1)=TWO
        nnn=nnn-2
      end do

      !! check
      write(*,*) "END NUMBER OF EL. PAIRS TO ASSIGN: ",nnn/2

      !! EVALUATING RMSD BETWEEN OCC. AND IDEAL ONES !!
      xrmsd=ZERO
      do icase=1,2
        do ii=1,lorb(icase)
          xx1=(elec2(ii,icase)-elec_id(ii,icase))
          xx1=xx1*xx1
          xrmsd=xrmsd+xx1
        end do
      end do
      xrmsd=xrmsd/REAL(nalf+nb)
      xrmsd=dsqrt(xrmsd)

      !! check
      write(*,*) "INITIAL RMSD: ",xrmsd

      write(*,*) " "
      write(*,*) " Number of pairs: ",npair
      write(*,*) " "

      !! NOW GIVING THE LAST ELECTRON PAIR FROM THE PAIRED DENSITY TO THE MOST OCC. UNPAIRED, ITERATIVELY !!
      ALLOCATE(tmp_elec_id(lorb(1),2))
      tmp_elec_id=elec_id
      do while(npair.gt.0)

        !! REMOVING LEAST OCCUPIED PAIRED EFO !!
        write(*,*) " Unoccupying paired EFOs with occ value of: ",elec2(npair,1)
        tmp_elec_id(npair,1)=ZERO

        !! ASSIGNING TO THE TWO MOST OCCUPIED UNPAIRED EFOs (BUT AVOIDING THE ALREADY FILLED) !!
        tmp_elec_id(nunp+1,2)=ONE
        tmp_elec_id(nunp+2,2)=ONE
        write(*,*) " Occupying unpaired EFOs with occ value of: ",elec2(nunp+1,2),elec2(nunp+2,2)

        !! EVALUATING AGAIN THE RMSD !!
        xrmsd2=ZERO
        do icase=1,2
          do ii=1,lorb(icase)
            xx1=(elec2(ii,icase)-tmp_elec_id(ii,icase))
            xx1=xx1*xx1
            xrmsd2=xrmsd2+xx1
          end do
        end do
        xrmsd2=xrmsd2/REAL(nalf+nb)
        xrmsd2=dsqrt(xrmsd2)
        
        !check!
        write(*,*) " New RMSD: ",xrmsd2
        write(*,*) " "


        !! IF NEW RMSD IS LOWER, SAVE ASSIGNMENT !!
        if((xrmsd2-xrmsd).lt.ZERO) then
          elec_id=tmp_elec_id
          xrmsd=xrmsd2

        !! IF NEW RMSD IS LARGER, STOP THE PROCESS !!
        else
          write(*,*) " Lowest RMSD found "
          write(*,*) " "
          go to 69
        end if

        !! FOR SHIFTING PROPERLY IN THE VECTOR !!
        nunp=nunp+2
        npair=npair-1
      end do
69    continue
      DEALLOCATE(tmp_elec_id)

      !! NOW WE HAVE THE ELECTRONIC ASIGNATION THAT MINIMIZES RMSD !!
      !! INTRODUCING IT INTO THE "ORIGINAL" EOS MACHINERY !!

      !! EVALUATING HOW MANY ELECTRONS WERE ASSIGNED TO EACH FRAGMENT !!
      !! ALSO COUNTING HOW MANY PAIRED AND UNPAIRED FOR LATER ON !!
      ALLOCATE(elec_frg_count(icufr,2))
      elec=ZERO
      elec_frg_count=ZERO
      do icase=1,2
        do ii=1,lorb(icase)
          if(elec_id(ii,icase).gt.ZERO) then
            ifrg=iorbat(ii,icase)
            elec(ifrg)=elec(ifrg)+elec_id(ii,icase)
            elec_frg_count(ifrg,icase)=elec_frg_count(ifrg,icase)+elec_id(ii,icase)
            write(*,*) "elec_id,ifrg: ",elec_id(ii,icase),ifrg
          end if
        end do
      end do
      write(*,*) " "
      write(*,*) " Printing electrons assigned separately, for check"
      write(*,*) " "
      do ifrg=1,icufr
        write(*,*) "Fragment, paired e, unpaired e: ",ifrg,elec_frg_count(ifrg,1),elec_frg_count(ifrg,2)
      end do
      write(*,*) " "

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

        if(icase.eq.1) xlast=TWO
        if(icase.eq.2) xlast=ONE
        ilast=0
        do ifrg=1,icufr

          !! CONDITIONS HAVE TO APPLY HERE !!
          if(icase.eq.1) then
            nn=INT(elec_frg_count(ifrg,icase))/2
            if(elec_frg_count(ifrg,icase)/TWO-nn.gt.thresh) nn=nn+1
          end if 
          if(icase.eq.2) then
            nn=INT(elec_frg_count(ifrg,icase))
            if(elec_frg_count(ifrg,icase)-nn.gt.thresh) nn=nn+1
          end if

!          write(*,*) "nn,elec_frg_count,iup0: ",nn,elec_frg_count(ifrg,icase),iup0(icase,ifrg)

          if(nn+1.gt.iup0(icase,ifrg)) then
            write(*,10) ifrg,elec_frg_count(ifrg,icase),up0gro(icase,nn,ifrg)
          else
            if(nn.ne.0) then
              write(*,15) ifrg,elec_frg_count(ifrg,icase),up0gro(icase,nn,ifrg),up0gro(icase,nn+1,ifrg)
            else
              write(*,15) ifrg,elec_frg_count(ifrg,icase),ZERO,up0gro(icase,nn+1,ifrg)
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

            !! CONDITIONS HAVE TO APPLY HERE, AGAIN !!
            if(icase.eq.1) then
              nn=INT(elec_frg_count(ifrg,icase))/2
              if(elec_frg_count(ifrg,icase)/TWO-nn.gt.thresh) nn=nn+1
            end if 
            if(icase.eq.2) then
              nn=INT(elec_frg_count(ifrg,icase))
              if(elec_frg_count(ifrg,icase)-nn.gt.thresh) nn=nn+1
            end if
            if(nn+1.ne.0.and.up0gro(icase,nn+1,ifrg).gt.xfirst) xfirst=up0gro(icase,nn+1,ifrg)
          end if
        end do
        write(*,*) " ---------------------------------------- "

        write(*,*) "xlast, xfirst: ",xlast,xfirst

        !! NOW PAIRED AND UNPAIRED ARE DIFFERENT (FACTOR OF 2 IN POPULATIONS) !!
        if(icase.eq.1) confi=100.0*dmin1(1.0d0,(xlast-xfirst)/TWO+0.5d0)
        if(icase.eq.2) confi=100.0*dmin1(1.0d0,xlast-xfirst+0.5d0)
        write(*,'(3x,a24,x,f7.3)') "RELIABILITY INDEX R(%) =",confi
        if(icase.eq.1) then
          write(*,*) "  INFO: (PAIRED) OCC. VALUES HALVED IN R(%) CALCULATION "
          confi0=confi
        end if
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

      confi2=dmin1(confi,confi0)
      write(*,'(2x,a32,x,f7.3)') "OVERALL RELIABILITY INDEX R(%) =",confi2
      write(*,'(2x,a34,x,f7.3)') "ELECTRONIC ASSIGNMENT RMSD VALUE =",xrmsd

      !! DEALLOCATING MATRICES !!
      DEALLOCATE(elec_id,elec2,elec_frg_count)

      !! PRINTING FORMATS !!
10    FORMAT(3x,i3,3x,f6.2,4x,f12.3,'    < thresh',2f12.3)
15    FORMAT(3x,i3,3x,f6.2,4x,f6.3,4x,f6.3)
20    FORMAT(3x,i3,6x,f8.2)

      end 

!! ****** !!
