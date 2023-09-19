!*****

      subroutine rwf_edatoiqafchk()

      implicit real*8(a-h,o-z)
      include 'parameter.h'

      common /hold/ihold(nmax),illim(nmax)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /c/ c(nmax,nmax),p(nmax,nmax)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /lim/ llim(nmax),iulim(nmax)
      common /stv/ sp(nmax,nmax),tt(nmax,nmax)
      common /filename/name0

      logical ilog

      character*60 name0,name1
      character*80 line

      allocatable :: cfr1(:,:),cfr2(:,:),cmol(:,:)
      allocatable :: s0(:,:),eigv(:,:),smh(:,:)
      allocatable :: cnew(:,:),pnew(:,:)
      allocatable :: snew(:,:)
      allocatable :: qmull(:)

!! READING THE COEFFICIENTS FROM THE FRAGMENT (TWO) .fchk FILES !!

      igr1=int_locate(55,"Number of basis",ilog)
      nbas1=int_locate(55,"Number of independ",ilog)
      if(nbas1.ne.igr1) then
        ibas1=1
        write(*,*) "INFO : Some basis from F1 has been removed "
      end if
      igr2=int_locate(52,"Number of basis",ilog)
      nbas2=int_locate(52,"Number of independ",ilog)
      if(nbas2.ne.igr2) then
        ibas2=1
        write(*,*) "INFO : Some basis from F2 has been removed "
      end if
      nocc1=int_locate(55,"Number of alpha electr",ilog)
      nocc2=int_locate(52,"Number of alpha electr",ilog)
      write(*,*) " Basis functions of Frag 1 : ",igr1," Indep : ",nbas1
      write(*,*) " Number of occupied orbitals of Frag 1 : ",nocc1
      write(*,*) " Basis functions of Frag 2 : ",igr2," Indep : ",nbas2
      write(*,*) " Number of occupied orbitals of Frag 2 : ",nocc2
      write(*,*) " "
      if(igr1+igr2.ne.igr) then
        write(*,*) " BASIS FUNCTIONS INCONSISTENCY "
        stop
      end if
      ALLOCATE(cfr1(igr1,igr1),cfr2(igr2,igr2))
      do ii=1,igr1
        do jj=ii,igr1
          cfr1(ii,jj)=ZERO
          cfr1(jj,ii)=ZERO
        end do
      end do
      do ii=1,igr2
        do jj=ii,igr2
          cfr2(ii,jj)=ZERO
          cfr2(jj,ii)=ZERO
        end do
      end do
      dummy=int_locate(55,"Alpha MO co",ilog)
      read(55,*)((cfr1(ii,jj),ii=1,igr1),jj=1,nbas1)
      if(ibas1.eq.1) then
        do ii=nbas1+1,igr1
          do jj=1,igr1
            cfr1(jj,ii)=ZERO
          end do
        end do
      end if
      dummy=int_locate(52,"Alpha MO co",ilog)
      read(52,*)((cfr2(ii,jj),ii=1,igr2),jj=1,nbas2)
      if(ibas2.eq.1) then
        do ii=nbas2+1,igr2
          do jj=1,igr2
            cfr2(jj,ii)=ZERO
          end do
        end do
      end if

!! RECONSTRUCTING THE BLOCK MATRIX Cmol !!

      ALLOCATE(cmol(igr,igr))
      do ii=1,igr
        do jj=ii,igr
          cmol(ii,jj)=ZERO
          cmol(jj,ii)=ZERO
        end do
      end do

!! PROPERLY ORDERED Cmol !!
!! OCCUPIED ORBITALS FIRST !!
!! FIRST FRAGMENT !!

      jmolbas=0
      do jj=1,nocc1
        ifrbas=0
        imolbas=0
        jmolbas=jmolbas+1
        do icenter=1,nat
          nbasi=iulim(icenter)-llim(icenter)+1
          if(jfrlist(icenter).eq.1) then
            do ibasis=1,nbasi
              imolbas=imolbas+1
              ifrbas=ifrbas+1
              cmol(imolbas,jmolbas)=cfr1(ifrbas,jj)
            end do
          else
            imolbas=imolbas+nbasi
          end if
        end do
      end do

!! SECOND FRAGMENT !!

      do jj=1,nocc2
        ifrbas=0
        imolbas=0
        jmolbas=jmolbas+1
        do icenter=1,nat
          nbasi=iulim(icenter)-llim(icenter)+1
          if(jfrlist(icenter).eq.2) then
            do ibasis=1,nbasi
              imolbas=imolbas+1
              ifrbas=ifrbas+1
              cmol(imolbas,jmolbas)=cfr2(ifrbas,jj)
            end do
          else
            imolbas=imolbas+nbasi
          end if
        end do
      end do

!! NOW THE VIRTUALS !!
!! FIRST FRAGMENT !!

      do jj=nocc1+1,igr1
        ifrbas=0
        imolbas=0
        jmolbas=jmolbas+1
        do icenter=1,nat
          nbasi=iulim(icenter)-llim(icenter)+1
          if(jfrlist(icenter).eq.1) then
            do ibasis=1,nbasi
              imolbas=imolbas+1
              ifrbas=ifrbas+1
              cmol(imolbas,jmolbas)=cfr1(ifrbas,jj)
            end do
          else
            imolbas=imolbas+nbasi
          end if
        end do
      end do

!! SECOND FRAGMENT !!

      do jj=nocc2+1,igr2
        ifrbas=0
        imolbas=0
        jmolbas=jmolbas+1
        do icenter=1,nat
          nbasi=iulim(icenter)-llim(icenter)+1
          if(jfrlist(icenter).eq.2) then
            do ibasis=1,nbasi
              imolbas=imolbas+1
              ifrbas=ifrbas+1
              cmol(imolbas,jmolbas)=cfr2(ifrbas,jj)
            end do
          else
            imolbas=imolbas+nbasi
          end if
        end do
      end do

!-----
!! chivatos !!
!     write(*,*) " "
!     write(*,*) " Printing Cfr1 "
!     write(*,*) " "
!     do jj=1,igr1
!       do ii=1,igr1
!         write(*,*) ii,jj,cfr1(ii,jj)
!       end do
!     end do
!     write(*,*) " "
!     write(*,*) " Printing Cfr2 "
!     write(*,*) " "
!     do jj=1,igr2
!       do ii=1,igr2
!         write(*,*) ii,jj,cfr2(ii,jj)
!       end do
!     end do
!     write(*,*) " "
!     write(*,*) " Printing Cmol "
!     write(*,*) " "
!     do jj=1,igr
!       do ii=1,igr
!         write(*,*) ii,jj,cmol(ii,jj)
!       end do
!     end do
!     write(*,*) " "
!-----
!     stop

      DEALLOCATE(cfr1,cfr2)

!! ADDING THE VIRTUALS INTO cnew, AS THEY WILL BE UNTOUCHED !!

      ALLOCATE(cnew(igr,igr))
      do jj=1,igr
        do ii=1,igr
          cnew(ii,jj)=cmol(ii,jj)
        end do
      end do

!! EVALUATING THE ELECTROSTATIC ENERGY BETWEEN FRAGMENTS BEFORE ORTHOGONALIZATION !!

      ndim=igr
      write(*,*) "Computing Charge penetration BEFORE Lowdin ortho"
      call edaiqa_frag_charge(ndim,nocc1,cmol)
      write(*,*) "Done"

!! DOING MULLIKEN CHARGES !!

      write(*,*) " "
      write(*,*) " MULLIKEN POPULATIONS (PRE LOWDIN) "
      write(*,*) " "
      ALLOCATE(pnew(igr,igr),qmull(nat))
      do iat=1,nat
        qmull(iat)=ZERO
      end do
      do ii=1,igr
        do jj=1,igr
          xx=ZERO
          do ij=1,nocc
            xx=xx+cmol(ii,ij)*cmol(jj,ij)
          end do
          pnew(ii,jj)=TWO*xx
        end do
      end do
      do mu=1,igr
        kat=ihold(mu)
        x=0.0d0
        do itau=1,igr
          x=x+pnew(mu,itau)*sp(itau,mu)
        end do
        qmull(kat)=qmull(kat)+x
      end do
      xxa=ZERO
      xxb=ZERO
      do iat=1,nat
        if(jfrlist(iat).eq.1) xxa=xxa+qmull(iat)
        if(jfrlist(iat).eq.2) xxb=xxb+qmull(iat)
      end do
      write(*,*) " Fragment Charges: ",xxa,xxb
      DEALLOCATE(pnew,qmull)
      write(*,*) " "

!! EVALUATING ELSTAT !!

      write(*,*) "Computing elstat"
      call rwf_classic_elect_ab(ndim,nocc1,cmol)
      write(*,*) "Done"

!! LOWDIN ORTHOGONALIZATION OF THE OCCUPIED MOs !!

      ALLOCATE(s0(nocc,nocc),eigv(nocc,nocc))
      do ii=1,nocc
        do jj=1,nocc
          xx=ZERO
          xx2=ZERO
          do mu=1,igr
            do nu=1,igr
              xx=xx+cmol(mu,ii)*cmol(nu,jj)*sp(mu,nu) 
            end do
          end do
          s0(ii,jj)=xx
        end do
      end do

      call diagonalize(nocc,nocc,s0,eigv,0)

      ALLOCATE(smh(nocc,nocc))
      do ii=1,nocc
        do jj=ii,nocc
          smh(jj,ii)=ZERO
          do kk=1,nocc
            if(s0(kk,kk).gt.thresh) then
              xx=eigv(ii,kk)*eigv(jj,kk)
              ssqrt=dsqrt(s0(kk,kk))
              smh(jj,ii)=smh(jj,ii)+xx/ssqrt
            end if
          end do
          smh(ii,jj)=smh(jj,ii)
        end do
      end do
      DEALLOCATE(eigv,s0)

      do ii=1,igr
        do jj=1,nocc
          xx=ZERO 
          do kk=1,nocc
            xx=xx+smh(jj,kk)*cmol(ii,kk)
          end do
          cnew(ii,jj)=xx
        end do
      end do
      DEALLOCATE(smh)

!! RECONSTRUCTING OVERLAP FOR CHECKING !!

      ALLOCATE(snew(nocc,nocc))
      do ii=1,nocc
        do jj=1,nocc
          xx=ZERO
          do mu=1,igr
            do nu=1,igr
              xx=xx+cnew(mu,ii)*cnew(nu,jj)*sp(mu,nu) 
            end do
          end do
          snew(ii,jj)=xx
        end do
      end do
      write(*,*) " "
      write(*,*) " Printing new MO overlaps (only occupied MOs) "
      write(*,*) " "
      do ii=1,nocc
        do jj=1,nocc
          if(snew(ii,jj).gt.thresh) write(*,*) ii,jj,snew(ii,jj)
        end do
      end do
      DEALLOCATE(snew)
!----
      
!! RECONSTRUCTING P MATRIX !!

      ALLOCATE(pnew(igr,igr))
      do ii=1,igr
        do jj=1,igr
          xx=ZERO
          do ij=1,nocc
            xx=xx+cnew(ii,ij)*cnew(jj,ij)
          end do
          pnew(ii,jj)=TWO*xx
        end do
      end do

!! DOING MULLIKEN CHARGES !!

      write(*,*) " "
      write(*,*) " MULLIKEN POPULATIONS (AFTER LOWDIN) "
      write(*,*) " "
      ALLOCATE(qmull(nat))
      do iat=1,nat
        qmull(iat)=ZERO
      end do
      do mu=1,igr
        kat=ihold(mu)
        x=0.0d0
        do itau=1,igr
          x=x+pnew(mu,itau)*sp(itau,mu)
        end do
        qmull(kat)=qmull(kat)+x
      end do
      xxa=ZERO
      xxb=ZERO
      do iat=1,nat
        if(jfrlist(iat).eq.1) xxa=xxa+qmull(iat)
        if(jfrlist(iat).eq.2) xxb=xxb+qmull(iat)
      end do
      write(*,*) " Fragment Charges: ",xxa,xxb
      DEALLOCATE(qmull)
      write(*,*) " "


!! COMPUTING CHARGE PENETRATION AFTER LOWDIN ORTHO !!

      write(*,*) "Computing Charge penetration AFTER Lowdin ortho"
      call edaiqa_frag_charge(ndim,nocc1,cnew)
      write(*,*) "Done"

!! PRINTING IN A NEW .fchk FILE !!
!! NOT ELEGANT WAY TO DO IT, BUT WORKS !!

      write(*,*) " Starting AuB fchk file construction "
      name1=trim(name0)//"-AuB.fchk"
      open(unit=69,file=name1)
      rewind(69)
      rewind(15)

!! PRINTING UNTIL ALPHA MOs !!

      read(15,'(a80)') line
      do while(index(line,"Alpha MO co").eq.0)
        write(69,'(a80)') line
        read(15,'(a80)') line
      end do

!! PRINTING THE NEW ONES !!

      indepigr=nbas1+nbas2
      norb=igr*indepigr
      norbt=igr*(igr+1)/2
      write(69,11) "Alpha MO coefficients","R","N= ",norb
      write(69,13) ((cnew(ii,jj),ii=1,igr),jj=1,indepigr)

!! NOW LOCATING WHAT IS AFTER IT IN THE ORIGINAL ONE TO CONTINUE !!

      do while(index(line,"Orthonormal basis").eq.0)
        read(15,'(a80)') line
      end do

!! RESTART PRINTING UNTIL NEXT STOP !!

      do while(index(line,"Total SCF Dens").eq.0)
        write(69,'(a80)') line
        read(15,'(a80)') line
      end do

!! PRINTING THE NEW ONE !!

      write(69,12) "Total SCF Density","R","N= ",norbt
      write(69,13) ((pnew(ii,jj),jj=1,ii),ii=1,igr)

!! NOW LOCATING WHAT IS AFTER IT IN THE ORIGINAL ONE TO CONTINUE !!

      do while(index(line,"Mulliken Charges").eq.0)
        read(15,'(a80)') line
      end do

!! RESTART PRINTING UNTIL THE END !!

      do while(.true.)
        write(69,'(a80)') line
        read(15,'(a80)',end=99) line
      end do
99    continue
      write(*,*) " AuB fchk file creation completed "
      close(69)

11    FORMAT(a21,22x,a1,3x,a3,i11)
12    FORMAT(a17,26x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))

      end

!*****

      subroutine uwf_edatoiqafchk()

      implicit real*8(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /c/ c(nmax,nmax),p(nmax,nmax)
      common /opensh/ cb(nmax,nmax)   
      common /ps/ps(nmax,nmax),pa(nmax,nmax),pb(nmax,nmax)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /lim/ llim(nmax),iulim(nmax)
      common /stv/ sp(nmax,nmax),tt(nmax,nmax)
      common /iops/iopt(100)

      logical ilog

      character*80 line

      allocatable :: cfr1(:,:),cfr2(:,:),cmol(:,:)
      allocatable :: cfr1b(:,:),cfr2b(:,:),cmolb(:,:)
      allocatable :: s0(:,:),eigv(:,:),smh(:,:)
      allocatable :: s0b(:,:),eigvb(:,:),smhb(:,:)
      allocatable :: cnew(:,:),cnewb(:,:),pnew(:,:),psnew(:,:)
      allocatable :: snew(:,:),snewb(:,:)

      iflip=iopt(93)

!! READING THE COEFFICIENTS FROM THE FRAGMENT (TWO) .fchk FILES !!

      igr1=int_locate(55,"Number of basis",ilog)
      igr2=int_locate(52,"Number of basis",ilog)
      nalf1=int_locate(55,"Number of alpha electr",ilog)
      nb1=int_locate(55,"Number of beta electr",ilog)
      nalf2=int_locate(52,"Number of alpha electr",ilog)
      nb2=int_locate(52,"Number of beta electr",ilog)
      write(*,*) " Number of total alpha and beta electrons : ",nalf,nb
      write(*,*) " Basis functions of Frag 1 : ",igr1
      write(*,*) " Number of alpha occupied orbitals of Frag 1 : ",nalf1
      write(*,*) " Number of beta occupied orbitals of Frag 1 : ",nb1
      if(igr1+igr2.ne.igr) then
        write(*,*) " BASIS FUNCTIONS INCONSISTENCY "
        stop
      end if

!! ASSUMING ALSO THAT TOTAL NUMBER OF ALPHA CANNOT BE SMALLER THAN BETA (GAUSSIAN-LIKE) !!

      if(nalf.lt.nb) then
        write(*,*) " MORE BETA ELECTRONS THAT ALPHA, KILL "
        stop
      end if

      ALLOCATE(cfr1(igr1,igr1),cfr2(igr2,igr2))
      ALLOCATE(cfr1b(igr1,igr1),cfr2b(igr2,igr2))
      do ii=1,igr1
        do jj=ii,igr1
          cfr1(ii,jj)=ZERO
          cfr1(jj,ii)=ZERO
          cfr1b(ii,jj)=ZERO
          cfr1b(jj,ii)=ZERO
        end do
      end do
      do ii=1,igr2
        do jj=ii,igr2
          cfr2(ii,jj)=ZERO
          cfr2(jj,ii)=ZERO
          cfr2b(ii,jj)=ZERO
          cfr2b(jj,ii)=ZERO
        end do
      end do
      dummy=int_locate(55,"Alpha MO co",ilog)
      read(55,*)((cfr1(ii,jj),ii=1,igr1),jj=1,igr1)
      dummy=int_locate(55,"Beta MO co",ilog)
      read(55,*)((cfr1b(ii,jj),ii=1,igr1),jj=1,igr1)

!! FLIPPING (OR NOT) THE SPINS OF THE SECOND FRAGMENT !!

      if(iflip.eq.0) then
        write(*,*) " Basis functions of Frag 2 : ",igr2
        write(*,*) " Number of alpha occupied orbitals of Frag 2 : ",nalf2
        write(*,*) " Number of beta occupied orbitals of Frag 2 : ",nb2
        write(*,*) " "
        dummy=int_locate(52,"Alpha MO co",ilog)
        read(52,*)((cfr2(ii,jj),ii=1,igr2),jj=1,igr2)
        dummy=int_locate(52,"Beta MO co",ilog)
        read(52,*)((cfr2b(ii,jj),ii=1,igr2),jj=1,igr2)

!! ALPHA AND BETA INTERCHANGED !!

      else
        iii=nb2
        nb2=nalf2
        nalf2=iii
        write(*,*) " Basis functions of Frag 2 : ",igr2
        write(*,*) " Number of alpha occupied orbitals of Frag 2 : ",nalf2
        write(*,*) " Number of beta occupied orbitals of Frag 2 : ",nb2
        write(*,*) " "
        write(*,*) " WARNING : FLIPSPIN FLAG ON "
        write(*,*) " Alpha to Beta & Beta to Alpha spinflips for frag2 "
        write(*,*) " "
        dummy=int_locate(52,"Alpha MO co",ilog)
        read(52,*)((cfr2b(ii,jj),ii=1,igr2),jj=1,igr2)
        dummy=int_locate(52,"Beta MO co",ilog)
        read(52,*)((cfr2(ii,jj),ii=1,igr2),jj=1,igr2)
      end if

!! RECONSTRUCTING THE BLOCK MATRICES Cmol AND Cmolb !!

      ALLOCATE(cmol(igr,igr))
      ALLOCATE(cmolb(igr,igr))
      do ii=1,igr
        do jj=ii,igr
          cmol(ii,jj)=ZERO
          cmol(jj,ii)=ZERO
          cmolb(ii,jj)=ZERO
          cmolb(jj,ii)=ZERO
        end do
      end do

!! PROPERLY ORDERED Cmol !!
!! OCCUPIED ORBITALS FIRST !!
!! FIRST FRAGMENT !!

      jmolbas=0
      jmolbasb=0
      do jj=1,nalf1
        ifrbas=0
        imolbas=0
        jmolbas=jmolbas+1
        if(jj.le.nb1) jmolbasb=jmolbasb+1
        do icenter=1,nat
          nbasi=iulim(icenter)-llim(icenter)+1
          if(jfrlist(icenter).eq.1) then
            do ibasis=1,nbasi
              imolbas=imolbas+1
              ifrbas=ifrbas+1
              cmol(imolbas,jmolbas)=cfr1(ifrbas,jj)
              if(jj.le.nb1) cmolb(imolbas,jmolbasb)=cfr1b(ifrbas,jj)
            end do
          else
            imolbas=imolbas+nbasi
          end if
        end do
      end do

!! SECOND FRAGMENT (CONSIDERING THE USE OF FLIPSPIN OR NOT IN THE SAME PART) !!

      norb2=nalf2
      if(nalf2.lt.nb2) norb2=nb2
      do jj=1,norb2
        ifrbas=0
        imolbas=0
        if(jj.le.nalf2) jmolbas=jmolbas+1
        if(jj.le.nb2) jmolbasb=jmolbasb+1
        do icenter=1,nat
          nbasi=iulim(icenter)-llim(icenter)+1
          if(jfrlist(icenter).eq.2) then
            do ibasis=1,nbasi
              imolbas=imolbas+1
              ifrbas=ifrbas+1
              if(jj.le.nalf2) cmol(imolbas,jmolbas)=cfr2(ifrbas,jj)
              if(jj.le.nb2) cmolb(imolbas,jmolbasb)=cfr2b(ifrbas,jj)
            end do
          else
            imolbas=imolbas+nbasi
          end if
        end do
      end do

!! NOW THE VIRTUALS (LOOP STARTING AT NBETA FOR FRAG1 AS .GT. NALPHA)!!
!! FIRST FRAGMENT !!

      do jj=nb1+1,igr1
        ifrbas=0
        imolbas=0
        if(jj.gt.nalf1) jmolbas=jmolbas+1
        jmolbasb=jmolbasb+1
        do icenter=1,nat
          nbasi=iulim(icenter)-llim(icenter)+1
          if(jfrlist(icenter).eq.1) then
            do ibasis=1,nbasi
              imolbas=imolbas+1
              ifrbas=ifrbas+1
              if(jj.gt.nalf1) cmol(imolbas,jmolbas)=cfr1(ifrbas,jj)
              cmolb(imolbas,jmolbasb)=cfr1b(ifrbas,jj)
            end do
          else
            imolbas=imolbas+nbasi
          end if
        end do
      end do

!! SECOND FRAGMENT (CONSIDERING THE USE OF FLIPSPIN OR NOT IN THE SAME PART) !!
!! ONCE AGAIN STARTING THE LOOP FROM THE SMALLER ONE !!

      norb2=nalf2
      if(nalf2.gt.nb2) norb2=nb2
      do jj=nocc2+1,igr2
        ifrbas=0
        imolbas=0
        if(jj2.gt.nalf2) jmolbas=jmolbas+1
        if(jj2.gt.nb2) jmolbasb=jmolbasb+1
        do icenter=1,nat
          nbasi=iulim(icenter)-llim(icenter)+1
          if(jfrlist(icenter).eq.2) then
            do ibasis=1,nbasi
              imolbas=imolbas+1
              ifrbas=ifrbas+1
              if(jj.gt.nalf2) cmol(imolbas,jmolbas)=cfr2(ifrbas,jj)
              if(jj.gt.nb2) cmolb(imolbas,jmolbasb)=cfr2b(ifrbas,jj)
            end do
          else
            imolbas=imolbas+nbasi
          end if
        end do
      end do

!!!!!

!-----
!! chivatos !!
!     write(*,*) " "
!     write(*,*) " Printing Cfr1 and Cfr1b "
!     write(*,*) " "
!     do jj=1,igr1
!       do ii=1,igr1
!         write(*,*) ii,jj,cfr1(ii,jj),cfr1b(ii,jj)
!       end do
!     end do
!     write(*,*) " "
!     write(*,*) " Printing Cfr2 and Cfr2b "
!     write(*,*) " "
!     do jj=1,igr2
!       do ii=1,igr2
!         write(*,*) ii,jj,cfr2(ii,jj),cfr2b(ii,jj)
!       end do
!     end do
!     write(*,*) " "
!     write(*,*) " Printing Cmol amb Cmolb "
!     write(*,*) " "
!     do jj=1,igr
!       do ii=1,igr
!         write(*,*) ii,jj,cmol(ii,jj),cmolb(ii,jj)
!       end do
!     end do
!     write(*,*) " "
!-----

      DEALLOCATE(cfr1,cfr2)
      DEALLOCATE(cfr1b,cfr2b)

!! ADDING THE VIRTUALS INTO COEFF MATRIX AS THEY WILL REMAIN UNCHANGED !!

      ALLOCATE(cnew(igr,igr))
      ALLOCATE(cnewb(igr,igr))
      do jj=1,igr
        do ii=1,igr
          cnew(ii,jj)=cmol(ii,jj)
          cnewb(ii,jj)=cmolb(ii,jj)
        end do
      end do

!! LOWDIN ORTHOGONALIZATION OF THE OCCUPIED MOs !!

      ALLOCATE(s0(nalf,nalf),eigv(nalf,nalf))
      ALLOCATE(s0b(nb,nb),eigvb(nb,nb))
      do ii=1,nalf
        do jj=1,nalf
          xx=ZERO
          xx2=ZERO
          xxb=ZERO
          xx2b=ZERO
          do mu=1,igr
            do nu=1,igr
              xx=xx+cmol(mu,ii)*cmol(nu,jj)*sp(mu,nu) 
              if(ii.le.nb.and.jj.le.nb) xxb=xxb+cmolb(mu,ii)*cmolb(nu,jj)*sp(mu,nu)
            end do
          end do
          s0(ii,jj)=xx
          if(ii.le.nb.and.jj.le.nb) s0b(ii,jj)=xxb
        end do
      end do

      call diagonalize(nalf,nalf,s0,eigv,0)
      call diagonalize(nb,nb,s0b,eigvb,0)

      ALLOCATE(smh(nalf,nalf))
      ALLOCATE(smhb(nb,nb))
      do ii=1,nalf
        do jj=ii,nalf
          smh(jj,ii)=ZERO
          if(ii.le.nb.and.jj.le.nb) smhb(jj,ii)=ZERO
          do kk=1,nalf
            if(s0(kk,kk).gt.thresh) then
              xx=eigv(ii,kk)*eigv(jj,kk)
              ssqrt=dsqrt(s0(kk,kk))
              smh(jj,ii)=smh(jj,ii)+xx/ssqrt
            end if
            if(ii.le.nb.and.jj.le.nb) then
              if(s0b(kk,kk).gt.thresh) then
                xx=eigvb(ii,kk)*eigvb(jj,kk)
                ssqrt=dsqrt(s0b(kk,kk))
                smhb(jj,ii)=smhb(jj,ii)+xx/ssqrt
              end if
            end if
          end do
          smh(ii,jj)=smh(jj,ii)
          if(ii.le.nb.and.jj.le.nb) smhb(ii,jj)=smhb(jj,ii)
        end do
      end do
      DEALLOCATE(eigv,s0)
      DEALLOCATE(eigvb,s0b)

!! NEW COEFFS FOR THE OCCUPIED SPIN-ORBITALS !!

      do ii=1,igr
        do jj=1,nalf
          xx=ZERO 
          xxb=ZERO 
          do kk=1,nalf
            xx=xx+smh(jj,kk)*cmol(ii,kk)
            if(jj.le.nb.and.kk.le.nb) xxb=xxb+smhb(jj,kk)*cmolb(ii,kk)
          end do
          cnew(ii,jj)=xx
          if(jj.le.nb) cnewb(ii,jj)=xxb
        end do
      end do
      DEALLOCATE(smh)
      DEALLOCATE(smhb)
!     DEALLOCATE(smh,cmol)
!     DEALLOCATE(smhb,cmolb)

!! RECONSTRUCTING OVERLAP FOR CHECKING !!

      ALLOCATE(snew(nalf,nalf))
      ALLOCATE(snewb(nb,nb))
      do ii=1,nalf
        do jj=1,nalf
          xx=ZERO
          xxb=ZERO
          do mu=1,igr
            do nu=1,igr
              xx=xx+cnew(mu,ii)*cnew(nu,jj)*sp(mu,nu) 
              if(ii.le.nb.and.jj.le.nb) xxb=xxb+cnewb(mu,ii)*cnewb(nu,jj)*sp(mu,nu) 
            end do
          end do
          snew(ii,jj)=xx
          if(ii.le.nb.and.jj.le.nb) snewb(ii,jj)=xxb
        end do
      end do
      write(*,*) " "
      write(*,*) " Printing new alpha overlaps "
      write(*,*) " "
      do ii=1,nalf
        do jj=1,nalf
          if(snew(ii,jj).gt.thresh) write(*,*) ii,jj,snew(ii,jj)
        end do
      end do
      write(*,*) " "
      write(*,*) " Printing new beta overlaps "
      write(*,*) " "
      do ii=1,nb
        do jj=1,nb
          if(snewb(ii,jj).gt.thresh) write(*,*) ii,jj,snewb(ii,jj)
        end do
      end do
      DEALLOCATE(snew)
      DEALLOCATE(snewb)
!----
      
!! RECONSTRUCTING P AND Ps MATRIX !!

      ALLOCATE(pnew(igr,igr))
      ALLOCATE(psnew(igr,igr))
      do ii=1,igr
        do jj=1,igr
          xx=ZERO
          xxb=ZERO
          do ij=1,nalf
            xx=xx+cnew(ii,ij)*cnew(jj,ij)
            if(ij.le.nb) xxb=xxb+cnewb(ii,ij)*cnewb(jj,ij)
          end do
          pnew(ii,jj)=xx+xxb
          psnew(ii,jj)=xx-xxb
        end do
      end do

!! PRINTING IN .fchk FORMAT !!

      norb=igr*igr
      norbt=igr*(igr+1)/2
      write(*,*) " "
      write(*,*) " ******************************* "
      write(*,*) " PRINTING ALPHA MOs COEFFICIENTS "
      write(*,*) " ******************************* "
      write(*,*) " "
      write(*,11) "Alpha MO coefficients","R","N= ",norb
      write(*,13) ((cnew(ii,jj),ii=1,igr),jj=1,igr)
      write(*,*) " "
      write(*,*) " ****************************** "
      write(*,*) " PRINTING BETA MOs COEFFICIENTS "
      write(*,*) " ****************************** "
      write(*,*) " "
      write(*,14) "Beta MO coefficients","R","N= ",norb
      write(*,13) ((cnewb(ii,jj),ii=1,igr),jj=1,igr)
      write(*,*) " "
      write(*,*) " ************************* "
      write(*,*) " PRINTING LOWER T P MATRIX "
      write(*,*) " ************************* "
      write(*,*) " "
      write(*,12) "Total SCF Density","R","N= ",norbt
      write(*,13) ((pnew(ii,jj),jj=1,ii),ii=1,igr)
      write(*,*) " "
      write(*,*) " ************************* "
      write(*,*) " PRINTING LOWER T Ps MATRIX "
      write(*,*) " ************************* "
      write(*,*) " "
      write(*,15) "Spin SCF Density","R","N= ",norbt
      write(*,13) ((psnew(ii,jj),jj=1,ii),ii=1,igr)
      write(*,*) " "

11    FORMAT(a21,22x,a1,3x,a3,i11)
12    FORMAT(a17,26x,a1,3x,a3,i11)
13    FORMAT(5(1p,e16.8))
14    FORMAT(a20,23x,a1,3x,a3,i11)
15    FORMAT(a16,27x,a1,3x,a3,i11)
      end

!*****

      subroutine rwf_classic_elect_ab(ndim,nocc1,cmol)

      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /c/ c(nmax,nmax),p(nmax,nmax)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /gridp/Nrad,Nang
      common /pha2/pha,phb,rr00
      common /iops/iopt(100)
      common /twoel/twoeltoler
      common /edaiqa/xen,xcoul,xnn
      common /edaiqa2/i2deda,iipoints,xptxyz(2,3)
      common /modgrid/nrad22,nang22,rr0022

      character*80 line

      integer, intent(in) :: ndim
      dimension cmol(ndim,ndim)

      allocatable :: pcoord(:,:),chp(:,:),omp2(:,:),wp(:),omp(:),ibaspoint(:)
      allocatable :: pcoordpha(:,:),chppha(:,:),omp2pha(:,:),wppha(:)
      allocatable :: rhoa(:),rhob(:)
      allocatable :: rhoapha(:),rhobpha(:),coul0(:,:),coul1(:,:)
      allocatable :: eelnuc(:,:),enn(:,:),ecoul(:,:),eelstat(:,:) 
      allocatable :: xnsymcoul(:,:)

      idofr = Iopt(40)
      nat0  = nat

!! DOING ZERO THE ENERGY MATRICES !!

      ALLOCATE(coul0(nat,1),coul1(maxat,maxat))
      ALLOCATE(eelnuc(maxat,maxat),enn(maxat,maxat),ecoul(maxat,maxat),eelstat(maxat,maxat))
      ALLOCATE(xnsymcoul(maxat,maxat))
      do ii=1,nat
        coul0(ii,1)=ZERO
        do jj=ii,nat
          coul1(ii,jj)=ZERO
          coul1(jj,ii)=ZERO
          eelnuc(ii,jj)=ZERO
          eelnuc(jj,ii)=ZERO
          enn(ii,jj)=ZERO
          enn(jj,ii)=ZERO
          ecoul(ii,jj)=ZERO
          ecoul(jj,ii)=ZERO
          eelstat(ii,jj)=ZERO
          eelstat(jj,ii)=ZERO

          xnsymcoul(ii,jj)=ZERO
          xnsymcoul(jj,ii)=ZERO
        end do
      end do

!! INTEGRATION GRID FOR THE ELECTRON-NUCLEI ATTRACTION !!

      write(*,*) " Grid for interfrag EN energy calculation "
      pha=ZERO
      phb=ZERO
      iatps=nang*nrad
      itotps=nrad*nang*nat
      call quad(Nrad,Nang) 
      ALLOCATE(wp(itotps),omp(itotps),omp2(itotps,nat))
      ALLOCATE(pcoord(itotps,3),ibaspoint(itotps),chp(itotps,ndim))
      call prenumint(ndim,itotps,nat0,wp,omp,omp2,chp,pcoord,ibaspoint,0)
      DEALLOCATE(omp,ibaspoint)
      write(*,*) " Radial points : ",nrad
      write(*,*) " Angular points : ",nang
      write(*,*) " rr00 value : ",rr00
      write(*,*) " "

!! COMPUTING RHO FOR FRAGMENTS A AND B !!

      ALLOCATE(rhoa(itotps),rhob(itotps))
      do kk=1,itotps
        xx1=ZERO
        xx2=ZERO
        do jj=1,nocc
          xx=ZERO
          do ii=1,igr
            xx=xx+cmol(ii,jj)*chp(kk,ii)
          end do

!! REMEMBER : FIRST OCCUPIED ORBITALS ARE FROM FRAGMENT 1, HERE TAGGED AS "A" !!

          if(jj.le.nocc1) xx1=xx1+xx*xx
          if(jj.gt.nocc1) xx2=xx2+xx*xx
        end do
        rhoa(kk)=TWO*xx1
        rhob(kk)=TWO*xx2
      end do
      DEALLOCATE(chp)

!! TESTING RHO !!

      write(*,*) " Integrating one-electron density to check "
      write(*,*) " INFO: Supersystem weights used "
      write(*,*) " "
      xxa=ZERO
      xxb=ZERO
      do icenter=1,nat
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          xxa=xxa+wp(ifut)*omp2(ifut,icenter)*rhoa(ifut)
          xxb=xxb+wp(ifut)*omp2(ifut,icenter)*rhob(ifut)
        end do
      end do
      write(*,*) " Frag 1 : ",xxa," Frag 2 : ",xxb
      write(*,*) " "

!! ELECTRON-NUCLEI ATTRACTION INTEGRALS !!

      write(*,*) " Evaluating EN integrals "

!! icenter USED FOR NUMERICAL INTEGRATION, jcenter GOVERNS IF CONTRIBUTES !!

      do icenter=1,nat
        do jcenter=1,nat
          zz=zn(jcenter)
          dx0=coord(1,jcenter)
          dy0=coord(2,jcenter)
          dz0=coord(3,jcenter)
          xx=ZERO
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            distx=pcoord(ifut,1)-dx0
            disty=pcoord(ifut,2)-dy0
            distz=pcoord(ifut,3)-dz0
            rr=dsqrt(distx*distx+disty*disty+distz*distz)
            if(rr.gt.1.0d-10) then
              if(jfrlist(jcenter).eq.2) xx=xx+zz*wp(ifut)*omp2(ifut,icenter)*rhoa(ifut)/rr
              if(jfrlist(jcenter).eq.1) xx=xx+zz*wp(ifut)*omp2(ifut,icenter)*rhob(ifut)/rr
            end if
          end do
          eelnuc(icenter,jcenter)=eelnuc(icenter,jcenter)-xx
        end do
      end do

!! PRINTING TO SEE THE AB AND BA DIFFERENCES !!

      write(*,*) " "
      write(*,*) " PRINTING ELSTAT BEFORE DOING AB = AB + BA "
      write(*,*) " CP CONSIDERED IN ATOMIC TERMS (IF IM NOT WRONG) "
      write(*,*) " "
      xtot=ZERO
      do ii=1,nat
        do jj=1,nat
          xtot=xtot+eelnuc(ii,jj)
        end do
      end do
      write(*,*) " Electron-Nuclei attraction terms "
      CALL Mprint(eelnuc,nat,maxat)
      write(*,*) " Electron-nuclei energy (au) : ",xtot

!! PRINTING INTEGRATION ERROR (CHECK WITH THE SYM) !!

      if(xen.ne.ZERO) then
        write(*,'(a31,f8.2)') "Integration error (kcal/mol): ",(xtot-xen)*23.06*27.212
        write(*,*) " "
      end if

!! GROUP BY FRAG FOR ASYM MATRICES !!

      if(idofr.eq.1) then
        line=' FRAGMENT ANALYSIS: Non-symmetryc Electron-Nuclei energy '
        call group_by_frag_mat(0,line,eelnuc)
      end if

!! AB = AB + BA AND PRINTING !!

      write(*,*) " "
      write(*,*) " PRINTING ELSTAT AFTER AB = AB + BA "
      write(*,*) " "
      xtot=ZERO
      do ii=1,nat
        do jj=1,ii-1  
          eelnuc(ii,jj)=eelnuc(ii,jj)+eelnuc(jj,ii)
          eelnuc(jj,ii)=eelnuc(ii,jj)
          xtot=xtot+eelnuc(ii,jj)
        end do
        xtot=xtot+eelnuc(ii,ii)
      end do
      write(*,*) " Electron-Nuclei attraction terms "
      CALL Mprint(eelnuc,nat,maxat)
      write(*,*) " Electron-nuclei energy (au) : ",xtot

!! PRINTING INTEGRATION ERROR !!

      if(xen.ne.ZERO) then
        write(*,'(a31,f8.2)') 'Integration error (kcal/mol): ',(xtot-xen)*23.06*27.212
        write(*,*) ' '
      end if

!! GROUPING BY FRAG !!

      if(idofr.eq.1) then
        line=' FRAGMENT ANALYSIS: Electron-Nuclei energy '
        call group_by_frag_mat(1,line,eelnuc)
      end if

!! DOING THE 2D PLOTS PART HERE !!

      if(i2deda.eq.1) call edaiqa_2dplots(ndim,itotps,nocc1,cmol,rhoa,rhob,wp,omp2,pcoord)

!! DEALLOCATING FOR TWO-ELECTRON INTEGRALS GRID !!

      DEALLOCATE(rhoa,rhob)
      DEALLOCATE(pcoord,omp2,wp)

!! INTEGRATION GRID FOR THE COULOMB INTERACTION AB TERMS !!
!! NO ROTATION FOR SECOND SET OF POINTS REQUIRED AS ONLY COMPUTING AB TERMS !!

      write(*,*) " "
      write(*,*) " CREATING NEW GRID FOR COULOMB ENERGY CALCULATION "
      write(*,*) " "
      nrad=nrad22
      nang=nang22
      rr00=rr0022
      iatps=nang*nrad
      itotps=nrad*nang*nat
      call quad(nrad,nang) 
      ALLOCATE(wp(itotps),omp(itotps),omp2(itotps,nat))
      ALLOCATE(pcoord(itotps,3),ibaspoint(itotps),chp(itotps,ndim))
      call prenumint(ndim,itotps,nat0,wp,omp,omp2,chp,pcoord,ibaspoint,0)
      write(*,*) " Radial points : ",nrad
      write(*,*) " Angular points : ",nang
      write(*,*) " rr00 value : ",rr00
      write(*,*) " "

!! COMPUTING RHO FOR FRAGMENTS A AND B !!

      ALLOCATE(rhoa(itotps),rhob(itotps))
      do kk=1,itotps
        xx1=ZERO
        xx2=ZERO
        do jj=1,nocc
          xx=ZERO
          do ii=1,igr
            xx=xx+cmol(ii,jj)*chp(kk,ii)
          end do

!! REMEMBER : FIRST OCCUPIED ORBITALS ARE FROM FRAGMENT 1, HERE TAGGED AS "A" !!

          if(jj.le.nocc1) xx1=xx1+xx*xx
          if(jj.gt.nocc1) xx2=xx2+xx*xx
        end do
        rhoa(kk)=TWO*xx1
        rhob(kk)=TWO*xx2
      end do
      DEALLOCATE(chp)

!! TESTING RHO !!

      write(*,*) " "
      write(*,*) " CHECK: INTEGRATING ONE-ELECTRON DENSITY "
      write(*,*) " INFO: SUPERSYSTEM WEIGHTS USED "
      write(*,*) " "
      xxa=ZERO
      xxb=ZERO
      do icenter=1,nat
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          xxa=xxa+wp(ifut)*omp2(ifut,icenter)*rhoa(ifut)
          xxb=xxb+wp(ifut)*omp2(ifut,icenter)*rhob(ifut)
        end do
      end do
      write(*,*) " Frag 1 : ",xxa," Frag 2 : ",xxb
      write(*,*) " "

!! NOW ROTATED GRID, ANGLE OPTIMIZED FOR BECKE 40 146 !!

      phb=0.162d0
      pha=ZERO
      write(*,*) " Rotating for angles : ",pha,phb
      ALLOCATE(wppha(itotps),omp2pha(itotps,nat))
      ALLOCATE(chppha(itotps,ndim),pcoordpha(itotps,3))
      call prenumint(ndim,itotps,nat0,wppha,omp,omp2pha,chppha,pcoordpha,ibaspoint,0)

!! JUST ROTATING THE DENSITY FROM FRAGMENT B !!

      ALLOCATE(rhoapha(itotps),rhobpha(itotps))
      do kk=1,itotps
        xx1=ZERO
        xx2=ZERO
        do jj=1,nocc
          xx=ZERO
          do ii=1,igr
            xx=xx+cmol(ii,jj)*chppha(kk,ii)
          end do
          if(jj.le.nocc1) xx1=xx1+xx*xx
          if(jj.gt.nocc1) xx2=xx2+xx*xx
        end do
        rhoapha(kk)=TWO*xx1
        rhobpha(kk)=TWO*xx2
      end do
      DEALLOCATE(chppha)

!! TESTING RHOPHA !!

      write(*,*) " "
      write(*,*) " CHECK: INTEGRATING ROTATED ONE-ELECTRON DENSITY "
      write(*,*) " INFO: SUPERSYSTEM WEIGHTS USED "
      write(*,*) " "
      xxa=ZERO
      xxb=ZERO
      do icenter=1,nat
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          xxa=xxa+wppha(ifut)*omp2pha(ifut,icenter)*rhoapha(ifut)
          xxb=xxb+wppha(ifut)*omp2pha(ifut,icenter)*rhobpha(ifut)
        end do
      end do
      write(*,*) " Frag 1 : ",xxa," Frag 2 : ",xxb
      write(*,*) " "

!! COULOMB REPULSION INTEGRALS !!

!     do icenter=1,nat
!       do ifut=iatps*(icenter-1)+1,iatps*icenter
!         x0=wp(ifut)*omp2(ifut,icenter)
!         dx0=pcoord(ifut,1)
!         dy0=pcoord(ifut,2)
!         dz0=pcoord(ifut,3)
!         do jcenter=1,nat

!! SAME CENTER !!

!           if(icenter.eq.jcenter) then
!             f3=ZERO
!             do jfut=iatps*(jcenter-1)+1,iatps*jcenter
!               x1=wppha(jfut)*omp2pha(jfut,jcenter)
!               dx1=pcoordpha(jfut,1)
!               dy1=pcoordpha(jfut,2)
!               dz1=pcoordpha(jfut,3)
!               dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
!               if(dist.gt.1.0d-12) f3=f3+rhoa(ifut)*rhobpha(jfut)*x1*x0/dist
!             end do
!           else

!! PAIR OF CENTERS !!

!             f3=ZERO
!             do jfut=iatps*(jcenter-1)+1,iatps*jcenter
!               x1=wppha(jfut)*omp2pha(jfut,jcenter)
!               dx1=pcoord(jfut,1)
!               dy1=pcoord(jfut,2)
!               dz1=pcoord(jfut,3)
!               dx1=pcoordpha(jfut,1)
!               dy1=pcoordpha(jfut,2)
!               dz1=pcoordpha(jfut,3)
!               dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
!               if(dist.gt.1.0d-12) f3=f3+rhoa(ifut)*rhob(jfut)*x1*x0/dist
!               if(dist.gt.1.0d-12) f3=f3+rhoa(ifut)*rhobpha(jfut)*x1*x0/dist
!             end do
!           end if
!           ecoul(icenter,jcenter)=ecoul(icenter,jcenter)+f3
!         end do
!       end do
!     end do


!! TESTING !!

      write(*,*) " "
      write(*,*) " COMPUTING COULOMB ENERGY TERMS "
      write(*,*) " "
      do icenter=1,nat
        do jcenter=1,nat
          f3=ZERO
          f4=ZERO
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            x0=wp(ifut)*omp2(ifut,icenter)
            dx0=pcoord(ifut,1)
            dy0=pcoord(ifut,2)
            dz0=pcoord(ifut,3)
            do jfut=iatps*(jcenter-1)+1,iatps*jcenter
              x1=wppha(jfut)*omp2pha(jfut,jcenter)
              dx1=pcoordpha(jfut,1)
              dy1=pcoordpha(jfut,2)
              dz1=pcoordpha(jfut,3)
              dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
!             if(dist.gt.1.0d-12) f3=f3+rhoa(ifut)*rhobpha(jfut)*x1*x0/dist
              if(dist.gt.1.0d-12) then
                f3=f3+rhoa(ifut)*rhobpha(jfut)*x1*x0/dist

!! HERE I AM DOING THE BA TERMS... LIKE AB = AB + BA !!

                f3=f3+rhob(ifut)*rhoapha(jfut)*x1*x0/dist

!! TRYING TO CONSTRUCT THE PIECES !!

                f4=f4+rhoa(ifut)*rhobpha(jfut)*x1*x0/dist
!               if(jfrlist(jcenter).eq.jfrlist(icenter)) f4=f4+rhoa(ifut)*rhobpha(jfut)*x1*x0/dist
!               if(jfrlist(jcenter).eq.2.and.jfrlist(icenter).eq.1) f4=f4+rhoa(ifut)*rhobpha(jfut)*x1*x0/dist
!               if(jfrlist(jcenter).eq.1.and.jfrlist(icenter).eq.2) f4=f4+rhob(ifut)*rhoapha(jfut)*x1*x0/dist
!               if(jfrlist(jcenter).ne.jfrlist(icenter)) f4=f4+rhoa(ifut)*rhobpha(jfut)*x1*x0/dist

              end if
            end do
          end do
          ecoul(icenter,jcenter)=ecoul(icenter,jcenter)+f3
          xnsymcoul(icenter,jcenter)=xnsymcoul(icenter,jcenter)+f4
        end do
      end do

!! NOW ROTATING rhoA INSTEAD OF rhoB !!

!     do icenter=1,nat
!       do jcenter=1,nat
!         f3=ZERO
!         f4=ZERO
!         do ifut=iatps*(icenter-1)+1,iatps*icenter
!           x0=wp(ifut)*omp2(ifut,icenter)
!           dx0=pcoord(ifut,1)
!           dy0=pcoord(ifut,2)
!           dz0=pcoord(ifut,3)
!             do jfut=iatps*(jcenter-1)+1,iatps*jcenter
!               x1=wppha(jfut)*omp2pha(jfut,jcenter)
!               dx1=pcoordpha(jfut,1)
!               dy1=pcoordpha(jfut,2)
!               dz1=pcoordpha(jfut,3)
!               dist=dsqrt((dx0-dx1)**TWO+(dy0-dy1)**TWO+(dz0-dz1)**TWO)
!               if(dist.gt.1.0d-12) then

!! HERE I AM DOING THE BA TERMS... LIKE AB = AB + BA !!

!                 f3=f3+rhob(ifut)*rhoapha(jfut)*x1*x0/dist

!! BUT NOT FOR THE NON-SYMM, IF NOT I FUCK IT UP !!
! test !
!                 f4=f4+rhoapha(jfut)*rhob(ifut)*x1*x0/dist
!               end if
!             end do
!         end do
!         ecoul(icenter,jcenter)=ecoul(icenter,jcenter)+f3
!         xnsymcoul(icenter,jcenter)=xnsymcoul(icenter,jcenter)+f4
!       end do
!     end do

!! FACTOR OF TWO HERE !!

      xtot=ZERO
      xtot2=ZERO
      do ii=1,nat
        ecoul(ii,ii)=ecoul(ii,ii)/TWO
        do jj=ii,nat
          ecoul(jj,ii)=ecoul(ii,jj)
          xtot=xtot+ecoul(ii,jj)
        end do

!! AND FOR THE NON-SYMM ALSO !!

        do jj=1,nat
!         xnsymcoul(ii,jj)=xnsymcoul(ii,jj)/TWO
          xtot2=xtot2+xnsymcoul(ii,jj)
        end do
      end do

      write(*,*) " "
      write(*,*) " NON-SYMMETRIC COULOMB ELECTRONIC REPULSION TERMS "
      CALL Mprint(xnsymcoul,nat,maxat)
      write(*,*) " Coulomb energy (au) : ",xtot2
      write(*,*) " "

!! ERROR FOR NON-SYMM, TO CHECK !!

      if(xcoul.ne.ZERO) then
        twoelerr=(xtot2-xcoul)*23.06*27.212
        write(*,'(a31,f8.2)') " Integration error (kcal/mol): ",twoelerr
        write(*,*) " "
      end if

!! GROUPING BY FRAG !!

      if(idofr.eq.1) then
        line=' FRAGMENT ANALYSIS: Non-symmetric Coulomb energy '
        call group_by_frag_mat(0,line,xnsymcoul)
      end if

!! NOW FOR THE SYMM COULOMB REPULSION MATRIX !!

      write(*,*) " "
      write(*,*) " COULOMB ELECTRONIC REPULSION TERMS "
      CALL Mprint(ecoul,nat,maxat)
      write(*,*) " Coulomb energy (au) : ",xtot

!! PRINTING INTEGRATION ERROR !!

      if(xcoul.ne.ZERO) then
        twoelerr=(xtot-xcoul)*23.06*27.212
        write(*,'(a31,f8.2)') " Integration error (kcal/mol): ",twoelerr
        write(*,*) ' '
      end if

!! GROUPING BY FRAG !!

      if(idofr.eq.1) then
        line=' FRAGMENT ANALYSIS: Coulomb energy '
        call group_by_frag_mat(1,line,ecoul)
      end if

!! DEALLOCATING !!

      DEALLOCATE(rhoa,rhoapha,rhob,rhobpha)
      DEALLOCATE(pcoord,omp2,wp)
      DEALLOCATE(pcoordpha,omp2pha,wppha)

!! NUCLEAR REPULSION !!

      xtot=ZERO
      do ii=1,nat
        do jj=ii+1,nat
          rr=dsqrt((coord(1,ii)-coord(1,jj))**TWO+(coord(2,ii)-coord(2,jj))**TWO+(coord(3,ii)-coord(3,jj))**TWO)
          if(jfrlist(ii).ne.jfrlist(jj)) then
            enn(ii,jj)=zn(ii)*zn(jj)/rr
            enn(jj,ii)=enn(ii,jj)
            xtot=xtot+enn(ii,jj)
          end if
        end do
      end do
      write(*,*) " Nuclear repulsion energy terms "
      CALL Mprint(enn,nat,maxat)
      write(*,*) " Nuclear repulsion energy (au) : ",xtot

!! PRINTING INTEGRATION ERROR !!

      if(xnn.ne.ZERO) then
        write(*,'(a31,f8.2)') 'Integration error (kcal/mol): ',(xtot-xnn)*23.06*27.212
        write(*,*) ' '
      end if

!! GROUPING BY FRAG !!

      if(idofr.eq.1) then
        line=' FRAGMENT ANALYSIS: Nuclear repulsion energy '
        call group_by_frag_mat(1,line,enn)
      end if

!! ADDING ALL TO RECREATE ELECTROSTATICS !!

      xtot=ZERO
      do ii=1,nat
        do jj=ii,nat
          eelstat(ii,jj)=eelnuc(ii,jj)+ecoul(ii,jj)+enn(ii,jj)
          eelstat(jj,ii)=eelstat(ii,jj)
          xtot=xtot+eelstat(ii,jj)
        end do
      end do
      write(*,*) " Electrostatic energy (kcal/mol) : ",xtot*23.06*27.212
      if(idofr.eq.1) then
        line=' FRAGMENT ANALYSIS: Electrostatic energy '
        call group_by_frag_mat(1,line,eelstat)
      end if

      DEALLOCATE(eelnuc,enn,ecoul,eelstat)

      end

!*****

      subroutine edaiqa_frag_charge(ndim,nocc1,cmol)

      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /gridp/Nrad,Nang
      common /pha2/pha,phb,rr00
      common /iops/iopt(100)

      character*80 line

      integer, intent(in) :: ndim
      dimension cmol(ndim,ndim)

      allocatable :: pcoord(:,:),chp(:,:),omp2(:,:),wp(:),omp(:),ibaspoint(:)
      allocatable :: rhoa(:),rhob(:)
      allocatable :: chga(:,:),chgb(:,:),chgt(:,:)

      idofr = Iopt(40)
      nat0  = nat

!! DOING ZERO THE ENERGY MATRICES !!

      ALLOCATE(chga(maxat,maxat),chgb(maxat,maxat),chgt(maxat,maxat))
      do ii=1,nat
        do jj=ii,nat
          chga(ii,jj)=ZERO
          chga(jj,ii)=ZERO
          chgb(ii,jj)=ZERO
          chgb(jj,ii)=ZERO
          chgt(ii,jj)=ZERO
          chgt(jj,ii)=ZERO
        end do
      end do

!! INTEGRATION GRID FOR THE ELECTRON-NUCLEI ATTRACTION !!

      write(*,*) " Grid for Charge calculation "
      nrad=150
      nang=974
      rr00=0.50d0
      pha=ZERO
      phb=ZERO
      iatps=nang*nrad
      itotps=nrad*nang*nat
      call quad(Nrad,Nang) 
      ALLOCATE(wp(itotps),omp(itotps),omp2(itotps,nat))
      ALLOCATE(pcoord(itotps,3),ibaspoint(itotps),chp(itotps,ndim))
      call prenumint(ndim,itotps,nat0,wp,omp,omp2,chp,pcoord,ibaspoint,0)
      DEALLOCATE(omp,ibaspoint)
      write(*,*) " Radial points : ",nrad
      write(*,*) " Angular points : ",nang
      write(*,*) " rr00 value : ",rr00
      write(*,*) " "

!! COMPUTING RHO FOR FRAGMENTS A AND B !!
!! APPLICABLE TO PRE AND POST LOWDIN ORTHO !!

      ALLOCATE(rhoa(itotps),rhob(itotps))
      do kk=1,itotps
        xx1=ZERO
        xx2=ZERO
        do jj=1,nocc
          xx=ZERO
          do ii=1,igr
            xx=xx+cmol(ii,jj)*chp(kk,ii)
          end do

!! REMEMBER : FIRST OCCUPIED ORBITALS ARE FROM FRAGMENT 1, HERE TAGGED AS "A" !!

          if(jj.le.nocc1) xx1=xx1+xx*xx
          if(jj.gt.nocc1) xx2=xx2+xx*xx
        end do
        rhoa(kk)=TWO*xx1
        rhob(kk)=TWO*xx2
      end do
      DEALLOCATE(chp)

!! COMPUTING CHARGE PENETRATION !!

      write(*,*) " Integrating one-electron density "
      write(*,*) " INFO: Supersystem weights used "
      write(*,*) " "
      do icenter=1,nat
        xxa=ZERO
        xxb=ZERO
        do ifut=iatps*(icenter-1)+1,iatps*icenter
          xxa=xxa+wp(ifut)*omp2(ifut,icenter)*rhoa(ifut)
          xxb=xxb+wp(ifut)*omp2(ifut,icenter)*rhob(ifut)
        end do
        chga(icenter,icenter)=chga(icenter,icenter)+xxa
        chgb(icenter,icenter)=chgb(icenter,icenter)+xxb
      end do

!! ADDING FOR TOTAL CHARGE !!

      do ii=1,nat
        chgt(ii,ii)=chga(ii,ii)+chgb(ii,ii)
      end do

!! PRINTING !!

      if(idofr.eq.1) then
        line=' FRAGMENT ANALYSIS: Charge from rhoA '
        call group_by_frag_mat(1,line,chga)
        line=' FRAGMENT ANALYSIS: Charge from rhoB '
        call group_by_frag_mat(1,line,chgb)
        line=' FRAGMENT ANALYSIS: Total Charge '
        call group_by_frag_mat(1,line,chgt)
      end if

!! DEALLOCATING !!

      DEALLOCATE(rhoa,rhob)
      DEALLOCATE(chga,chgb,chgt)
      DEALLOCATE(pcoord,omp2,wp)

      end

!*****

      subroutine edaiqa_2dplots(ndim,itotps,nocc1,cmol,rhoa,rhob,wp,omp2,pcoord)

      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      parameter (toau=0.52917721067d0)

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /edaiqa2/i2deda,iipoints,xptxyz(2,3)

      integer, intent(in) :: ndim,itotps,nocc1
      dimension cmol(ndim,ndim),rhoa(itotps),rhob(itotps),wp(itotps),omp2(itotps,nat),pcoord(itotps,3)

      allocatable :: rpt(:,:),rdist(:)
      allocatable :: vatot(:),vaeff(:),vbtot(:),vbeff(:)

      iatps=itotps/nat

!! COMPUTING FIRST THE r POINTS !!

      ALLOCATE(rpt(iipoints,3),rdist(iipoints))

!! WORKING DISTANCES AND ATOMIC POSITIONS IN AU !!

      extra=TWO
      call edaiqa_top_pointpos_1d(extra,xptxyz,rpt,rdist)
      write(*,*) " "
      write(*,*) " PRINTING xyz OF THE r POINTS "
      write(*,*) " "
      write(*,*) " Point, Distance, xyz "
      write(*,*) " "
      do ii=1,iipoints

!! PASSANT A ANGSTROMS !!

        x1=rpt(ii,1)*toau
        y1=rpt(ii,2)*toau
        z1=rpt(ii,3)*toau
        rdist(ii)=rdist(ii)*toau
        write(*,'(i5,4f10.5)') ii,rdist(ii),x1,y1,z1
      end do
      write(*,*) " "

!! EVALUATING TOTAL AND EFFECTIVE POTENTIALS !!

      ALLOCATE(vatot(iipoints),vaeff(iipoints),vbtot(iipoints),vbeff(iipoints))
      do ii=1,iipoints
        vaeff(ii)=ZERO
        vbeff(ii)=ZERO
        vatot(ii)=ZERO
        vbtot(ii)=ZERO

!! FIRST INTEGRATING rho(r2)/|r - r2| !!

        do icenter=1,nat
          dx0=rpt(ii,1)
          dy0=rpt(ii,2)
          dz0=rpt(ii,3)
          xxa=ZERO
          xxb=ZERO
          xxatot=ZERO
          xxbtot=ZERO
          do ifut=iatps*(icenter-1)+1,iatps*icenter
            distx=pcoord(ifut,1)-dx0
            disty=pcoord(ifut,2)-dy0
            distz=pcoord(ifut,3)-dz0
            rr=dsqrt(distx*distx+disty*disty+distz*distz)
            if(rr.gt.1.0d-8) then

!! OWN FRAGMENT (NET/EFFECTIVE POTENTIAL) !!

              if(jfrlist(icenter).eq.1) xxa=xxa-wp(ifut)*omp2(ifut,icenter)*rhoa(ifut)/rr
              if(jfrlist(icenter).eq.2) xxb=xxb-wp(ifut)*omp2(ifut,icenter)*rhob(ifut)/rr

!! ALL SPACE (TOTAL POTENTIAL) !!

              xxatot=xxatot-wp(ifut)*omp2(ifut,icenter)*rhoa(ifut)/rr
              xxbtot=xxbtot-wp(ifut)*omp2(ifut,icenter)*rhob(ifut)/rr
            end if
          end do

!! INCLUDING NOW Zi/|r -Ri| !!

          zz=zn(icenter)
          distx=coord(1,icenter)-dx0
          disty=coord(2,icenter)-dy0
          distz=coord(3,icenter)-dz0
          rr=dsqrt(distx*distx+disty*disty+distz*distz)
          if(rr.gt.1.0d-8) then

!! HERE BOTH ADDS ONLY OWN FRAGMENT !!

            if(jfrlist(icenter).eq.1) then
              xxa=xxa+zz/rr
              xxatot=xxatot+zz/rr
            else if(jfrlist(icenter).eq.2) then
              xxb=xxb+zz/rr
              xxbtot=xxbtot+zz/rr
            end if
          end if
          vaeff(ii)=vaeff(ii)+xxa
          vbeff(ii)=vbeff(ii)+xxb
          vatot(ii)=vatot(ii)+xxatot
          vbtot(ii)=vbtot(ii)+xxbtot
        end do 
      end do

!! PRINTING !!

      write(*,*) " "
      write(*,*) " PRINTING POTENTIALS EVALUATED AT EACH r POINT "
      write(*,*) " "
      write(*,*) " FORMAT: dist, vatot, vaeff, vbtot, vbeff"
      write(*,*) " "
      do ii=1,iipoints
        if(vatot(ii).gt.10000) write(*,*) "vatot larger than 10000" 
        if(vatot(ii).lt.-10000) write(*,*) "vatot lower than -10000" 
        if(vbtot(ii).gt.10000) write(*,*) "vbtot larger than 10000" 
        if(vbtot(ii).lt.-10000) write(*,*) "vbtot lower than -10000" 
        write(*,'(5f14.6)') rdist(ii),vatot(ii),vaeff(ii),vbtot(ii),vbeff(ii) 
      end do
      write(*,*) " "

!! NOW EVALUATING AND PRINTING FRAGMENT DENSITIES !!

      write(*,*) " "
      write(*,*) " START OF 2D TOPOLOGY "
      write(*,*) " "
      call edaiqa_top_dens_1d(ndim,nocc1,cmol,rpt,rdist)
      write(*,*) " "
      write(*,*) " END OF 2D TOPOLOGY "
      write(*,*) " "
      call flush

      DEALLOCATE(rpt,rdist)
      DEALLOCATE(vatot,vaeff,vbtot,vbeff)

      end

!*****

      subroutine edaiqa_top_pointpos_1d(extra,atcoord,pcoord,ddist)
      implicit double precision(a-h,o-z)
      include 'parameter.h'
      common /coord/coord(3,maxat),zn(maxat),iznuc(maxat)
      common /edaiqa2/i2deda,iipoints,xptxyz(2,3)

      parameter (toau=0.52917721067d0)
      dimension :: atcoord(2,3),pcoord(iipoints,3),ddist(iipoints)

!! TRANSFORMING POSITIONS FROM .INP TO ATOMIC UNITS !!

      do iat=1,2
        do jj=1,3
          atcoord(iat,jj)=atcoord(iat,jj)/toau
        end do
      end do

!! COMPUTING UNITARY VECTOR BETWEEN A PAIR OF ATOMS !!

      xvect=atcoord(2,1)-atcoord(1,1)
      yvect=atcoord(2,2)-atcoord(1,2)
      zvect=atcoord(2,3)-atcoord(1,3)
      xmod=dsqrt(xvect*xvect+yvect*yvect+zvect*zvect)
!! check
      write(*,*) "pollas,dist",xmod*toau
      xvect=xvect/xmod
      yvect=yvect/xmod
      zvect=zvect/xmod

!! POSITION OF A MARGIN POINT (FROM FIRST ATOM) !!

      xpt1=atcoord(1,1)-extra*xvect
      ypt1=atcoord(1,2)-extra*yvect
      zpt1=atcoord(1,3)-extra*zvect

!! POSITION OF THE CENTRAL POINT (USED AS x-axis WHILE REPRESENTING) !!

      xcent=(atcoord(2,1)+atcoord(1,1))/TWO
      ycent=(atcoord(2,2)+atcoord(1,2))/TWO
      zcent=(atcoord(2,3)+atcoord(1,3))/TWO

      totdist=xmod+TWO*extra
      xx=totdist/(REAL(iipoints)-ONE)
      do ii=1,iipoints

!! EVALUTING THE POSITION OF THE POINTS ALONG AN INTERATOMIC DISTANCE (pcoord) !!

        ffact=xx*(REAL(ii)-ONE)
        pcoord(ii,1)=xpt1+ffact*xvect
        pcoord(ii,2)=ypt1+ffact*yvect
        pcoord(ii,3)=zpt1+ffact*zvect

!! COMPUTING THE DISTANCE FROM THE CENTRAL POINT (ddist)

        x0=pcoord(ii,1)-xcent
        y0=pcoord(ii,2)-ycent
        z0=pcoord(ii,3)-zcent
        dist=(dsqrt(x0*x0+y0*y0+z0*z0))
        if(x0*xvect.le.ZERO.and.y0*yvect.le.ZERO.and.z0*zvect.le.ZERO) then
          ddist(ii)=-dist
        else
          ddist(ii)=dist
        end if
      end do

!! NOW DISTANCE TO EACH OF THE "TWO" ATOMS !!

      x1=atcoord(1,1)-xcent
      y1=atcoord(1,2)-ycent
      z1=atcoord(1,3)-zcent
      xd1=dsqrt(x1*x1+y1*y1+z1*z1)
      if(x1*xvect.le.ZERO.and.y1*yvect.le.ZERO.and.z1*zvect.le.ZERO) xd1=-xd1 
      x1=atcoord(2,1)-xcent
      y1=atcoord(2,2)-ycent
      z1=atcoord(2,3)-zcent
      xd2=dsqrt(x1*x1+y1*y1+z1*z1)
      if(x1*xvect.le.ZERO.and.y1*yvect.le.ZERO.and.z1*zvect.le.ZERO) xd2=-xd2 

!! TRANSFORMINT TO ANGSTROMS FOR PRINTING !!

      write(*,*) " "
      write(*,*) " CHECK:"
      write(*,*) " "
      write(*,*) " DISTANCE TO ATOM 1 :",xd1*toau
      write(*,*) " DISTANCE TO ATOM 2 :",xd2*toau
      write(*,*) " "

      end

!*****

      subroutine edaiqa_top_dens_1d(ndim,nocc1,cmol,pcoord,rdist)
      implicit double precision(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /edaiqa2/i2deda,iipoints,xptxyz(2,3)
      common /iops/iopt(100)
      common /c/ c(nmax,nmax),p(nmax,nmax)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      integer, intent(in) :: ndim

      dimension :: cmol(ndim,ndim),pcoord(iipoints,3),rdist(iipoints)

      allocatable :: gx_ao(:),gy_ao(:),gz_ao(:),eval_ao(:)
      allocatable :: rhoa(:),rhob(:)
      allocatable :: rhot(:),csave(:,:)
      allocatable :: omp(:),rhoaa(:),rhoab(:),rhoba(:),rhobb(:)

      i5d = Iopt(11)

      if(kop.eq.0) iwf=1
      if(kop.eq.1) iwf=2
      write(*,*) " "
      if(iwf.eq.1) write(*,*) " DENSITY TOPOLOGY FOR RWF "
      if(iwf.eq.2) write(*,*) " DENSITY TOPOLOGY FOR UWF "
      write(*,*) " "
      ALLOCATE(gx_ao(igr),gy_ao(igr),gz_ao(igr),eval_ao(igr))
      ALLOCATE(rhoa(iipoints),rhob(iipoints))
      ALLOCATE(rhot(iipoints),csave(igr,igr))

!! for check total density
      do ii=1,igr
        do jj=1,igr
          csave(ii,jj)=c(ii,jj)
          c(ii,jj)=cmol(ii,jj)
        end do
      end do 
!!

      ALLOCATE(omp(nat))
      ALLOCATE(rhoaa(iipoints),rhoab(iipoints),rhoba(iipoints),rhobb(iipoints))
      do ii=1,iipoints
        x0=pcoord(ii,1)
        y0=pcoord(ii,2)
        z0=pcoord(ii,3)
        call gpoints(x0,y0,z0,gx_ao,gy_ao,gz_ao,eval_ao)
        if(i5d.eq.1) call gordermat(eval_ao,igr)

!! EVALUATING ATOMIC WEIGHTS AT THE POINT !!

        do iat=1,nat
          omp(iat)=wat(iat,x0,y0,z0)
        end do


!! COMPUTING RHO AND STORAGING IT !!

        if(iwf.eq.1) then
          call edaiqa_calc_rhf_dens(ndim,eval_ao,nocc1,cmol,rhoa0,rhob0)
          rhoa(ii)=rhoa0
          rhob(ii)=rhob0

!! EVALUATING THE "EFFECTIVE" RHOs !!

          rhoa00=ZERO
          rhoa01=ZERO
          rhob00=ZERO
          rhob01=ZERO
          do iat=1,nat
            if(jfrlist(iat).eq.1) rhoa00=rhoa00+rhoa0*omp(iat) !! A in A !!
            if(jfrlist(iat).eq.2) rhoa01=rhoa01+rhoa0*omp(iat) !! A in B !!
            if(jfrlist(iat).eq.1) rhob01=rhob01+rhob0*omp(iat) !! B in A !!
            if(jfrlist(iat).eq.2) rhob00=rhob00+rhob0*omp(iat) !! B in B !!
          end do
          rhoaa(ii)=rhoa00
          rhoab(ii)=rhoa01
          rhobb(ii)=rhob00
          rhoba(ii)=rhob01

!! total for check !!
          call calc_rhf_dens(eval_ao,rhotot)
          rhot(ii)=rhotot

        else if(iwf.eq.2) then
!         call calc_uhf_dens(eval_ao,rhoa,rhob)
!         rho(ii)=rhoa+rhob
          write(*,*) " NOT IMPLEMENTED YET "
          stop
        end if
      end do
      DEALLOCATE(gx_ao,gy_ao,gz_ao,eval_ao)
      DEALLOCATE(omp)

!! PRINTING !!

      write(*,*) " "
      write(*,*) " PRINTING DENSITIES OF FRAGMENTS AT EACH r POINT "
      write(*,*) " "
      write(*,*) " FORMAT: dist, rhoa, rhob, rhotot "
      write(*,*) " "
      do ii=1,iipoints
        write(*,'(4f12.4)') rdist(ii),rhoa(ii),rhob(ii),rhot(ii)
      end do
      write(*,*) " "

!! PRINTING EFFECTIVE DENSITIES !!

      write(*,*) " PRINTING EFFECTIVE DENSITIES "
      write(*,*) " "
      write(*,*) " FORMAT: rhoaa, rhoab, rhobb, rhoba "
      write(*,*) " "
      do ii=1,iipoints
        write(*,'(4f12.4)') rhoaa(ii),rhoab(ii),rhobb(ii),rhoba(ii)
      end do
      write(*,*) " "


!! recovering c !!
      do ii=1,igr
      do jj=1,igr
      c(ii,jj)=csave(ii,jj)
      end do
      end do 
!!

      DEALLOCATE(rhoa,rhob,rhot,csave)
      DEALLOCATE(rhoaa,rhoab,rhoba,rhobb)

      end

!*****

      subroutine edaiqa_calc_rhf_dens(ndim,eval_ao,nocc1,cmol,rhoa0,rhob0)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop

      integer, intent(in) :: ndim

      dimension :: eval_ao(igr),cmol(ndim,ndim)

!     write(*,*) " for safety, nocc1: ",nocc1
      xxa0=ZERO
      xxb0=ZERO
      do jj=1,nocc
        xx=ZERO
        do ii=1,igr
          xx=xx+cmol(ii,jj)*eval_ao(ii)
        end do

!! REMEMBER : FIRST OCCUPIED ORBITALS ARE FROM FRAGMENT 1, HERE TAGGED AS "A" !!

        if(jj.le.nocc1) xxa0=xxa0+xx*xx
        if(jj.gt.nocc1) xxb0=xxb0+xx*xx
      end do
      rhoa0=TWO*xxa0
      rhob0=TWO*xxb0

      end

!*****

