      SUBROUTINE VPRINT(H,N,ndim,jdim)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim,jdim
      common /printout/iaccur
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      DIMENSION H(NDIM,jdim)

      Dimension mend(92)
      data mend/4H  H ,4H He ,4H Li ,4H Be ,4H  B ,4H  C ,4H  N ,4H  O , 
     $ 4H  F ,4H Ne ,4H Na ,4H Mg ,4H Al ,4H Si ,4H  P ,4H  S ,4H Cl ,
     $4H Ar ,4H  K ,4H Ca ,4H Sc ,4H Ti ,4H  V ,4H Cr ,4H Mn ,4H Fe ,
     $4H Co ,4H Ni ,4H Cu ,4H Zn ,4H Ga ,4H Ge ,4H As ,4H Se ,4H Br ,
     $4H Kr ,4H Rb ,4H Sr ,4H  Y ,4H Zr ,4H Nb ,4H Mo ,4H Tc ,4H Ru ,
     $4H Rh ,4H Pd ,4H Ag ,4H Cd ,4H In ,4H Sn ,4H Sb ,4H Te ,4H  I ,
     $4H Xe ,4H Cs ,4H Ba ,4H La ,4H Ce ,4H Pr ,4H Nd ,4H Pm ,4H Sn ,
     $4H Eu ,4H Gd ,4H Tb ,4H Dy ,4H Ho ,4H Er ,4H Tm ,4H Yb ,4H Lu , 
     $4H Hf ,4H Ta ,4H  W ,4H Re ,4H Os ,4H Ir ,4H Pt ,4H Au ,4H Hg ,
     $4H Tl ,4H Pb ,4H Bi ,4H Po ,4H At ,4H Rn ,4H Fr ,4H Ra ,4H Ac ,
     $4H Th ,4H Pa ,4H  U   /
    
      K=6
      NMIN=1
      NNMAX=1
   62 FORMAT(1X,I3,A4,6F12.6)
   63 FORMAT(1X,I3,A4,6F20.13)
      DO 2 I=1,N
      if(iaccur.eq.0) then
      PRINT 62,I,mend(iznuc(i)),(H(I,J),J=1,jdim)
      else
      PRINT 63,I,mend(iznuc(i)),(H(I,J),J=1,jdim)
      end if
   2  CONTINUE
      RETURN
      END

      SUBROUTINE MPRINT_NLOP(H,N,ndim)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim
      common /printout/iaccur
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      DIMENSION H(NDIM,3)

      Dimension mend(92)
      data mend/4H  H ,4H He ,4H Li ,4H Be ,4H  B ,4H  C ,4H  N ,4H  O , 
     $ 4H  F ,4H Ne ,4H Na ,4H Mg ,4H Al ,4H Si ,4H  P ,4H  S ,4H Cl ,
     $4H Ar ,4H  K ,4H Ca ,4H Sc ,4H Ti ,4H  V ,4H Cr ,4H Mn ,4H Fe ,
     $4H Co ,4H Ni ,4H Cu ,4H Zn ,4H Ga ,4H Ge ,4H As ,4H Se ,4H Br ,
     $4H Kr ,4H Rb ,4H Sr ,4H  Y ,4H Zr ,4H Nb ,4H Mo ,4H Tc ,4H Ru ,
     $4H Rh ,4H Pd ,4H Ag ,4H Cd ,4H In ,4H Sn ,4H Sb ,4H Te ,4H  I ,
     $4H Xe ,4H Cs ,4H Ba ,4H La ,4H Ce ,4H Pr ,4H Nd ,4H Pm ,4H Sn ,
     $4H Eu ,4H Gd ,4H Tb ,4H Dy ,4H Ho ,4H Er ,4H Tm ,4H Yb ,4H Lu , 
     $4H Hf ,4H Ta ,4H  W ,4H Re ,4H Os ,4H Ir ,4H Pt ,4H Au ,4H Hg ,
     $4H Tl ,4H Pb ,4H Bi ,4H Po ,4H At ,4H Rn ,4H Fr ,4H Ra ,4H Ac ,
     $4H Th ,4H Pa ,4H  U   /
    
      NMIN=1
      NNMAX=3
      PRINT *, ('_________________X___________________Y__________________Z__________')
      PRINT *,' '
      DO I=1,N
       if(iaccur.eq.0) then
        PRINT 62,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
       else
        PRINT 63,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
       end if
      end do
      RETURN
      PRINT 66
   62 FORMAT(1X,I3,A4,6F12.6)
   63 FORMAT(1X,I3,A4,3F20.13)
   66 FORMAT(1X///)
      END

      SUBROUTINE MPRINT(H,N,ndim)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim
      common /printout/iaccur
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      DIMENSION H(NDIM,NDIM)


      Dimension mend(92)
      data mend/4H  H ,4H He ,4H Li ,4H Be ,4H  B ,4H  C ,4H  N ,4H  O , 
     $ 4H  F ,4H Ne ,4H Na ,4H Mg ,4H Al ,4H Si ,4H  P ,4H  S ,4H Cl ,
     $4H Ar ,4H  K ,4H Ca ,4H Sc ,4H Ti ,4H  V ,4H Cr ,4H Mn ,4H Fe ,
     $4H Co ,4H Ni ,4H Cu ,4H Zn ,4H Ga ,4H Ge ,4H As ,4H Se ,4H Br ,
     $4H Kr ,4H Rb ,4H Sr ,4H  Y ,4H Zr ,4H Nb ,4H Mo ,4H Tc ,4H Ru ,
     $4H Rh ,4H Pd ,4H Ag ,4H Cd ,4H In ,4H Sn ,4H Sb ,4H Te ,4H  I ,
     $4H Xe ,4H Cs ,4H Ba ,4H La ,4H Ce ,4H Pr ,4H Nd ,4H Pm ,4H Sn ,
     $4H Eu ,4H Gd ,4H Tb ,4H Dy ,4H Ho ,4H Er ,4H Tm ,4H Yb ,4H Lu , 
     $4H Hf ,4H Ta ,4H  W ,4H Re ,4H Os ,4H Ir ,4H Pt ,4H Au ,4H Hg ,
     $4H Tl ,4H Pb ,4H Bi ,4H Po ,4H At ,4H Rn ,4H Fr ,4H Ra ,4H Ac ,
     $4H Th ,4H Pa ,4H  U   /
    
      K=6
      NMIN=1
      NNMAX=MIN0(N,K)
   62 FORMAT(1X,I3,A4,6F12.6)
   63 FORMAT(1X,I3,A4,6F20.13)
   1  if(iaccur.eq.0) then
      PRINT 60, (I, mend(iznuc(i)),I=NMIN,NNMAX)
      else
      PRINT 61, (I, mend(iznuc(i)),I=NMIN,NNMAX)
      end if
   60 FORMAT(10X,6(2X,I3,A4,3X))
   61 FORMAT(10X,6(6X,I3,A4,7X))
      PRINT 64
   64 FORMAT(1X)
      DO 2 I=1,N
      if(iaccur.eq.0) then
      PRINT 62,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
      else
      PRINT 63,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
      end if
   2  CONTINUE
      NMIN=NMIN+6
      K=K+6
      NNMAX=MIN0(N,K)
      IF(NNMAX.GE.NMIN) GOTO 71
      RETURN
   71 PRINT 66
   66 FORMAT(1X///)
      GO TO 1
      END

!! ***** !!

      SUBROUTINE MPRINT2(H,N,ndim)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim
      character*100 line
      common /printout/iaccur
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      DIMENSION H(NDIM,NDIM)

      Dimension mend(92)
      data mend/4H  H ,4H He ,4H Li ,4H Be ,4H  B ,4H  C ,4H  N ,4H  O , 
     $ 4H  F ,4H Ne ,4H Na ,4H Mg ,4H Al ,4H Si ,4H  P ,4H  S ,4H Cl ,
     $4H Ar ,4H  K ,4H Ca ,4H Sc ,4H Ti ,4H  V ,4H Cr ,4H Mn ,4H Fe ,
     $4H Co ,4H Ni ,4H Cu ,4H Zn ,4H Ga ,4H Ge ,4H As ,4H Se ,4H Br ,
     $4H Kr ,4H Rb ,4H Sr ,4H  Y ,4H Zr ,4H Nb ,4H Mo ,4H Tc ,4H Ru ,
     $4H Rh ,4H Pd ,4H Ag ,4H Cd ,4H In ,4H Sn ,4H Sb ,4H Te ,4H  I ,
     $4H Xe ,4H Cs ,4H Ba ,4H La ,4H Ce ,4H Pr ,4H Nd ,4H Pm ,4H Sn ,
     $4H Eu ,4H Gd ,4H Tb ,4H Dy ,4H Ho ,4H Er ,4H Tm ,4H Yb ,4H Lu , 
     $4H Hf ,4H Ta ,4H  W ,4H Re ,4H Os ,4H Ir ,4H Pt ,4H Au ,4H Hg ,
     $4H Tl ,4H Pb ,4H Bi ,4H Po ,4H At ,4H Rn ,4H Fr ,4H Ra ,4H Ac ,
     $4H Th ,4H Pa ,4H  U   /
    
      line="--------------------------------------------------------------------------------"

      K=6
      NMIN=1
      NNMAX=MIN0(N,K)
   62 FORMAT(2X,I3,A4,6F12.6)
   63 FORMAT(2X,I3,A4,6F20.13)
   1  if(iaccur.eq.0) then
      PRINT 666,line
      PRINT 60, (I, mend(iznuc(i)),I=NMIN,NNMAX)
      else
      PRINT 666,line
      PRINT 61, (I, mend(iznuc(i)),I=NMIN,NNMAX)
      end if
   60 FORMAT(10X,6(2X,I3,A4,3X))
   61 FORMAT(10X,6(6X,I3,A4,7X))
      PRINT 666,line
      DO 2 I=1,N
      if(iaccur.eq.0) then
      PRINT 62,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
      else
      PRINT 63,I,mend(iznuc(i)),(H(I,J),J=NMIN,NNMAX)
      end if
   2  CONTINUE
      PRINT 666,line
      NMIN=NMIN+6
      K=K+6
      NNMAX=MIN0(N,K)
      IF(NNMAX.GE.NMIN) GOTO 71
      RETURN
   71 PRINT 66  
   66 FORMAT(1X)
      GO TO 1

!! FORMAT FOR THE LINES !!
  666 FORMAT(2x,a100)
      END

!! ***** !!


      SUBROUTINE MPRINTNOAT(H,M,N,mdim,ndim)
      IMPLICIT REAL*8 (A-H,O-Z)
      include 'parameter.h'
      integer, intent(in) :: ndim,mdim
      DIMENSION H(MDIM,NDIM)
      common /printout/iaccur

      K=6
      NMIN=1
      NNMAX=MIN0(N,K)
   62 FORMAT(1X,I3,6X,6F12.6)
   63 FORMAT(1X,I3,X,6F20.13)
    1 continue                              
      DO 2 I=1,M
      if(iaccur.eq.0) then
      PRINT 62,I,(H(I,J),J=NMIN,NNMAX)
      else
      PRINT 63,I,(H(I,J),J=NMIN,NNMAX)
      end if
   2  CONTINUE
      NMIN=NMIN+6
      K=K+6
      NNMAX=MIN0(N,K)
      IF(NNMAX.GE.NMIN) GOTO 71
      RETURN
   71 PRINT 66
   66 FORMAT(1X//)
      GO TO 1
      END
    
CCCC
C GROUP PRINTING
CCCC
      subroutine group_by_frag_mat(ilog,line,A)
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /printout/iaccur
      dimension A(maxat,maxat),B(maxat,maxat)
      character*80 line
      integer ilog

c ilog=0 do all matrix, ilog=1 lower triangular

      do i=1,nat
       do j=1,nat
        B(i,j)=0.0d0
       end do
      end do

      nnat=nat
      do i=1,nat
       if(ilog.eq.1) nnat=i
       do j=1,nnat
        if(jfrlist(i).ne.0.and.jfrlist(j).ne.0) B(jfrlist(i),jfrlist(j))=B(jfrlist(i),jfrlist(j))+a(i,j)
       end do
      end do

      x=0.0d0
      if(ilog.eq.1) then
       do i=1,icufr
        x=x+b(i,i)
        do j=1,i-1
         B(i,j)=B(i,j)+B(j,i)
         x=x+b(i,j)
         B(j,i)=B(i,j)
        end do
       end do
      else
      do i=1,icufr
       do j=1,icufr
        x=x+b(i,j)
       end do
      end do
      end if

      print 60
      write(*,'(a80)') line 
      write(*,*) ' '
      call mprintnoat(B,icufr,icufr,maxat,maxat)
      write(*,*) ' '
      if(iaccur.eq.0) then
      PRINT 63,x
      else
      PRINT 64,x
      end if
      write(*,*) ' '

      write(*,*) 'Additional printing :'
      do i=1,icufr
       do j=1,icufr
        if(iaccur.eq.0) then
         write(*,'(i3,i3,1X,F12.6)') i,j,b(i,j)
        else
         write(*,'(i3,i3,1X,F20.13)') i,j,b(i,j)
        end if
       end do
      end do
      write(*,*) ' '

   60 FORMAT(1X//)
   63 FORMAT('  Total:',5X,F20.13)
   64 FORMAT('  Total:',5X,F20.13)
      end

      subroutine group_by_frag_vec(ndim,line,A)
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      character*80 line
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /printout/iaccur
      dimension A(maxat,ndim),B(maxat,ndim)
      dimension x(ndim)

      do i=1,nat
       do j=1,ndim
        B(i,j)=0.0d0
       end do
      end do

      do j=1,ndim
       do i=1,nat
        if(jfrlist(i).ne.0) B(jfrlist(i),j)=B(jfrlist(i),j)+a(i,j)
       end do
      end do

      do j=1,ndim
      x(j)=0.0d0
      do i=1,icufr
       x(j)=x(j)+b(i,j)
      end do
      end do

      print 60
      write(*,'(a80)') line 
      write(*,*) ' '
      if(ndim.eq.2) then
      write(*,*) '             3D-space   Mulliken'
      write(*,*) '           ----------------------'
      else
      write(*,*) '             3D-space '
      write(*,*) '           ------------'
      end if
      call mprintnoat(B,icufr,ndim,maxat,ndim)
      if(ndim.eq.2) then
      write(*,*) '           ----------------------'
      else
      write(*,*) '           ------------'
      end if
      if(iaccur.eq.0) then
      write(*,63)(x(j),j=1,ndim)
      else
      write(*,64)(x(j),j=1,ndim)
      end if
   63 FORMAT('  Total:',1X,2(F12.6,1X))
   64 FORMAT('  Total:',1X,2(F20.13,1X))
   60 FORMAT(1X//)
      end

CCCCC
C FOR WRITTING INT FILES
CCCCC
      subroutine print_int_old(nbas,nat0,sat,name)
      use basis_set
      use ao_matrices, only: cb,c_no,occ_no !excluding c because otherwise there's problems
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /atlist/iatlist(maxat),icuat
      common/cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /coord/ coord0(3,maxat),zn(maxat),iznuc(maxat)
      dimension sat(nbas,nbas,nat0)
      character name*60
      character*2 mend(92)
      data mend/' H','He','Li','Be',' B',' C',' N',' O',
     $ ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl',
     $ 'Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe',
     $ 'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
     $ 'Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     $ 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I',
     $ 'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sn',
     $ 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
     $ 'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     $ 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     $ 'Th','Pa',' U'  /

      character nameaim*55,charnu*2,charnu2*2,charnu1
      character charnu3*3,ext*4
      character*80 line

      allocatable c3(:,:),c(:,:),c2(:,:),csave(:,:),scr(:,:)
      allocatable clin(:)


C IOPS
       inofu = Iopt(2)
       ihirsh = Iopt(6)
       imulli=Iopt(5) 
       icorr=Iopt(26)   
       iqtaim =Iopt(16)

       inato=0
       if(icorr.eq.1.or.icas.eq.1.or.icisd.eq.1) inato=1

        igr0=igr
        rewind(15)
 998    read(15,'(a80)') line
        if(index(line,"Number of independant functions").ne.0) then
         read(line(54:61),'(i8)') igr0
        else if(index(line,"Number of independent functions").ne.0) then
         read(line(54:61),'(i8)') igr0
        else
         go to 998
        end if

        allocate( c3(2*igr,2*igr))
        allocate( c(igr,igr),c2(igr,igr),csave(igr,igr),scr(igr,igr))

        if(inofu.eq.0) then

c getting Alpha MO
      rewind(15)
 999    read(15,'(a80)') line
        if(index(line,"Alpha MO coefficients").ne.0) then
         ntriang=igr*igr0
         allocate(clin(ntriang))
         read(15,*)(clin(i),i=1,ntriang)
      else
       go to 999
      end if
      k=1
      do i=1,igr0
         do j=1,igr
        c(j,i)=clin(k)
        k=k+1
       end do
        end do

C check orthogonality of MO 
c transform to MO 
       write(*,*) 'Checking MO overlap matrix '

        do i=1,igr0
         do j=1,igr
          xx=0.0d0
          do k=1,igr
           xx=xx+C(k,i)*S(k,j)
          end do
          c2(i,j)=xx
         end do
        end do
        do i=1,igr0
         do j=1,igr0
          xx=0.0d0
          do k=1,igr
           xx=xx+C2(i,k)*c(k,j)  
          end do
         if(i.eq.j.and.abs(xx-1.0d0).gt.1.0d-5) write(*,*) i,j,xx          
         if(i.ne.j.and.abs(xx).gt.1.0d-5) write(*,*) i,j,xx      
         end do
        end do

        do i=1,igr
         do j=1,igr
         scr(i,j)=c(i,j)
         end do
        end do

c ok..but now it ll get the orbitals from wfn
        if(inato.eq.1) then
c         rewind(92)
c         read(92) imo,((c(i,k),i=1,imo),k=1,imo)
c        if(imo.ne.igr) then
c        write(*,*) 'inconsistency problem',imo,igr
c        stop
c       end if
c         close(92)
c
C making trnafomation matrix from MO to NO
c transform to MO
        do i=1,igr
         do j=1,igr
          xx=0.0d0
          do k=1,igr
           xx=xx+C(k,i)*S(k,j)
          end do
          c2(i,j)=xx
         end do
        end do
c        do i=1,igr
c         do j=1,igr
c          xx=0.0d0
c          do k=1,igr
c           xx=xx+C2(i,k)*c_no(k,j)
c          end do
c          c3(i,j)=xx
c         end do
c        end do
c writting transformation matrix on disk
c         write(69) igr,((c3(i,k),i=1,igr),k=1,igr)

        do i=1,igr
         do j=1,igr
         scr(i,j)=c_no(i,j)
         end do
        end do

       write(*,*) 'Checking NO overlap matrix '
        do i=1,igr
         do j=1,igr
          xx=0.0d0
          do k=1,igr
           xx=xx+C_no(k,i)*S(k,j)
          end do
          c2(i,j)=xx
         end do
        end do
        do i=1,igr
         do j=1,igr
          xx=0.0d0
          do k=1,igr
           xx=xx+C2(i,k)*c_no(k,j)  
          end do
         if(i.eq.j.and.abs(xx-1.0d0).gt.1.0d-5) write(*,*) i,j,xx          
         if(i.ne.j.and.abs(xx).gt.1.0d-5) write(*,*) i,j,xx      
         end do
        end do
        end if

       end if

      do i=1,igr
       do j=1,igr
        csave(i,j)=0.0d0                
       end do
      end do
 
c building file
      naim=88
      naim3=89
      l=len_trim(name)
c      write(*,*) l,'----'
      if(imulli.eq.1) then
       nameaim=name(1:l)//"mul.files"
       ext="mul_"
      else if(imulli.ge.2) then
       nameaim=name(1:l)//"low.files"
       ext="low_"
      else if(iqtaim.eq.1) then
       nameaim=name(1:l)//"aim.files"
       ext="aim_"
      else if(ihirsh.eq.0) then
       nameaim=name(1:l)//"fuz.files"
       ext="fuz_"
      else if(ihirsh.eq.1) then
       nameaim=name(1:l)//"hir.files"
       ext="hir_"
      else if(ihirsh.eq.2) then
       nameaim=name(1:l)//"ihi.files"
       ext="ihi_"
      end if
      open(file=nameaim,unit=naim3,status="unknown")
      rewind(naim3)
      do jjat=1,icuat
        if(icuat.ne.nat) then
         jat=iatlist(jjat)
        else
         jat=jjat
        end if
       read(mend(iznuc(jat)),'(A2)')charnu  
         charnu=adjustl(charnu)
        l1=len_trim(name)
        l2=len_trim(charnu)
        if(jat.lt.10) then 
         write(charnu1,'(i1)')jat
         nameaim=name(1:l1)//ext//trim(charnu)//charnu1
      else if(jat.lt.100) then
         write(charnu2,'(i2)')jat
         nameaim=name(1:l1)//ext//trim(charnu)//charnu2         
      else
         write(charnu3,'(i3)')jat
         nameaim=name(1:l1)//ext//trim(charnu)//charnu3         
c assuming up to 999 atoms
      end if
       j=len(nameaim)
       do i=1,j
        if(nameaim(i:i).eq.' ') then
         l=i-1
         go to 10
        end if
       end do
  10   continue
       nameaim=nameaim(1:l)//".int"
       write(naim3,'(a40)') nameaim
      open(file=nameaim,unit=naim,status="unknown")
      rewind(naim)

c     writting stuff will be looked for by FCALC
        BK=1.0d0
        AE=-1.0d0
      write(naim,'(a20)') 'Created by APOST3D '
        if(iopt(5).eq.1) then
         write(naim,*) 'Using Mulliken atomic definition '
        else if(iopt(5).ge.2) then
         write(naim,*) 'Using Lowdin atomic definition '
        else if(iopt(6).eq.1) then
         write(naim,*) 'Using Hirshfeld atomic definition '
        else if(iopt(6).eq.2) then
         write(naim,*) 'Using Hirshfeld-Iterative atomic definition '
        else if(iopt(14).eq.1.and.iopt(31).eq.0) then
         write(naim,*) 'Using Becke-rho atomic definition '
        else if(iopt(14).eq.1.and.iopt(31).eq.1) then
         write(naim,*) 'Using TFVC atomic definition '
        else if(iopt(16).eq.1) then
         write(naim,*) 'Using QTAIM atomic definition '
        else
         write(naim,*) 'Using Becke atomic definition '
        end if
         if(iopt(5).eq.0) write(naim,*) 'Stiffness parameter k = ',iopt(25)
        if(inato.eq.1)  then
         write(naim,*) 'Correlated wave function '
         if(icorr.eq.0) then
          write(naim,*) '*Warning*, no RDM1 provided.' 
          write(naim,*) 'Can not do alpha and beta populations separatedly.' 
         end if
        else
         write(naim,*) 'Single-determinant wave function '
         if(kop.eq.1) write(naim,*) 'Unrestricted wave function '
        end if

      write(naim,'(a30,F20.11)')' MOLECULAR SCF ENERGY (AU)  = ',AE
      write(naim,*)''
      write(naim,'(A25,A4,i5)')' INTEGRATION IS OVER ATOM',mend(iznuc(jat)),
     1  jat
      !write(naim2) mend(iznuc(jat)),jat
        write(naim,'(a27)') ' RESULTS OF THE INTEGRATION'
        write(naim,'(a17,e21.14,a14,e21.14)')'              N  ',
     1   qat(jjat,1),'    NET CHARGE',-qat(jjat,1)+iznuc(jjat)
        !write(naim2) qat(jjat,1)
        write(naim,*)'             G'
        write(naim,'(a17,e21.14,a17,e21.14)')'              K  ',BK,
     1  '        E(ATOM)  ',AE
        write(naim,'(a17,e21.14)')'              L  ',0.0d0
      write(naim,*)''
        write(naim,'(a35)')'          The Atomic Overlap Matrix'
        write(naim,*)' '
        if(inato.eq.1) then
        write(naim,*)' '
        else if(kop.eq.1) then
        write(naim,'(a36)')'Unrestricted  Wavefunction'
        write(naim,*)' '
        else
        write(naim,'(a36)')'Restricted Closed-Shell Wavefunction'
        write(naim,*)' '
        end if

        do i=1,2*igr
         do j=1,2*igr
           c3(i,j)=0.0d0
         end do
        end do

        if(inofu.eq.0) then
c transform 
        numorb=nalf
        if(inato.eq.1) numorb=igr
      do i=1,numorb
       do j=1,igr
        xx=0.0d0
        do k=1,igr
         xx=xx+scr(k,i)*Sat(k,j,jjat)
        end do
        c2(i,j)=xx
       end do
        end do
      do i=1,numorb
       do j=1,numorb
        xx=0.0d0
        do k=1,igr
         xx=xx+C2(i,k)*scr(k,j)       
        end do
        c3(i,j)=xx
       end do
        end do
C open-shell
        if(kop.eq.1.and.inato.eq.0)  then
      do i=1,nb    
       do j=1,igr
        xx=0.0d0
        do k=1,igr
         xx=xx+Cb(k,i)*Sat(k,j,jjat)
        end do
        c2(i,j)=xx
       end do
        end do
      do i=1,nb   
       do j=1,nb   
        xx=0.0d0
        do k=1,igr
         xx=xx+C2(i,k)*cb(k,j)       
        end do
        c3(i+nalf,j+nalf)=xx
       end do
        end do
        end if
        
        else

C PSS somethign was wrong here using c3  17/7/2015
c      do i=1,nocc
c       do j=1,nocc
c        c3(i,j)=sat(i,j,jjat)
c       end do
c        end do


        end if

c         do i=1,igr 
c          write(naim,*) (sat(i,j,jat),j=1,i)
c         end do

      do i=1,igr
       do j=1,igr
        csave(i,j)=csave(i,j)+c3(i,j)
       end do
      end do

c only occupied
c         do i=1,igr0
c unrestricted single-determinant
         if(kop.eq.1.and.inato.eq.0) then
          if (imulli.ne.1) then
          do i=1,nalf+nb 
           write(naim,*) (c3(i,j),j=1,i)
          end do
          else
          do i=1,nalf+nb 
           write(naim,*) (c3(i,j),j=1,nalf+nb)
          end do
          end if
c recalculate spin populations
          qalf=0.0d0
          do i=1,nalf 
            qalf=qalf+c3(i,i)
          end do
          qbet=0.0d0
          do i=1,nb   
            qbet=qbet+c3(i+nalf,i+nalf)
          end do
          write(naim,*) '  '               
          write(naim,'(a41,e21.14)') 'ALPHA ELECTRONS (NA)',qalf     
          write(naim,'(a41,e21.14)') 'BETA ELECTRONS (NB)',qbet        
c restricted single-determinant
         else if(inato.eq.0) then
          if (imulli.ne.1) then
          do i=1,nocc
           write(naim,*) (c3(i,j),j=1,i)
          end do
          else
          do i=1,nocc
           write(naim,*) (c3(i,j),j=1,nocc)
          end do
          end if
          write(naim,*) '  '               
          write(naim,'(a41,e21.14)') 'ALPHA ELECTRONS (NA)',qat(jjat,
     +    1)/2.0d0
          write(naim,'(a41,e21.14)') 'BETA ELECTRONS (NB)',qat(jjat,
     +    1)/2.0d0
c correlated wave function        
         else if(inato.eq.1) then
          if (imulli.ne.1) then
          do i=1,numorb
           write(naim,*) (c3(i,j),j=1,i)
          end do
          else
          do i=1,numorb
           write(naim,*) (c3(i,j),j=1,numorb)
          end do
          end if
          write(naim,*) '  '               
          write(naim,'(a41,e21.14)') 'ALPHA ELECTRONS (NA)',qat(jjat,
     +    1)/2.0d0
          write(naim,'(a41,e21.14)') 'BETA ELECTRONS (NB)',qat(jjat,
     +    1)/2.0d0
         end if


        write(naim,*)' '
      write(naim,*) 'NORMAL TERMINATION OF PROAIMV'

      close(naim)
        end do
      close(naim3)


        if(icuat.eq.nat) then
c check for orthogonality of the MOs
      xmax=1.0d-2
      xmaxd=1.0d-3
      xmaxns=1.0d-4
c         do i=1,igr
c          write(*,*) (csave(i,j),j=1,igr)
c         end do
c      xmaxg=1.0d-2
      do i=1,igr
       do j=i,igr
          if(i.eq.j) then
         if(abs(csave(i,i))-1.0d0.gt.xmax) then 
          xx=abs(csave(i,i))-1.0d0
          imaxd=i
            if(i.le.nalf) then
          write(*,*) 'Dev. from normalization ',xx,imaxd
            end if
         end if
        else 
         if(abs(csave(i,j)).gt.xmaxd) then 
          xx=abs(csave(i,j))
          imax1=i
          imax2=j
            if(i.le.nalf.and.j.le.nalf) then
          write(*,'(a27,f10.6,2i3)') 'Dev. from orth. in occ. set:',xx,
     1      imax1,imax2
            end if
         end if
         if(j.gt.i) then
          xxx=abs(csave(i,j))-abs(csave(j,i)) 
          if(xxx.gt.xmaxns) then 
           xx=xxx
           imaxns1=i   
           imaxns2=j   
           write(*,*) 'Dev. from hermiticity',xx,imaxns1,imaxns2
          end if
         end if
        end if
       end do
      end do
        end if
         deallocate(c3,c,scr,c2,csave,clin)

      end

      subroutine print_int(nbas,nat0,sat,name)
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /atlist/iatlist(maxat),icuat
      common/cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      common /coord/ coord0(3,maxat),zn(maxat),iznuc(maxat)
      dimension sat(nbas,nbas,nat0)
      character name*60
      character*2 mend(92)
      data mend/' H','He','Li','Be',' B',' C',' N',' O',
     $ ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl',
     $ 'Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe',
     $ 'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
     $ 'Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     $ 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I',
     $ 'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sn',
     $ 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
     $ 'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     $ 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     $ 'Th','Pa',' U'  /

      character nameaim*55,charnu*2,charnu2*2,charnu1
      character charnu3*3,ext*4
      character*80 line

      allocatable c3(:,:),c2(:,:),csave(:,:),scr(:,:)


C IOPS
       ihirsh = Iopt(6)
       imulli=Iopt(5) 
       icorr=Iopt(26)   
       iqtaim =Iopt(16)

       inato=0
       if(icorr.eq.1.or.icas.eq.1.or.icisd.eq.1) inato=1

        rewind(15)
 998    read(15,'(a80)') line
        if(index(line,"Number of independant functions").ne.0) then
         read(line(54:61),'(i8)') ndim
        else if(index(line,"Number of independent functions").ne.0) then
         read(line(54:61),'(i8)') ndim
        else
         go to 998
        end if

        allocate( c3(2*igr,2*igr))
        allocate( c2(igr,igr),csave(igr,igr),scr(igr,igr))
        scr=0.0d0

        do i=1,igr
         do j=1,igr
         scr(i,j)=c(i,j)
         end do
        end do

        if(inato.eq.1) then
         if(icas.eq.1) ndim=norb
         write(*,*) 'Using',ndim,' natural orbitals'
         do i=1,igr
          do j=1,ndim
          scr(i,j)=c_no(i,j)
          end do
         end do
        end if

       do i=1,igr
        do j=1,ndim
         csave(i,j)=0.0d0                
        end do
       end do
 
c building file
      naim=88
      naim3=89
      l=len_trim(name)
c      write(*,*) l,'----'
      if(imulli.eq.1) then
       nameaim=name(1:l)//"mul.files"
       ext="mul_"
      else if(imulli.ge.2) then
       nameaim=name(1:l)//"low.files"
       ext="low_"
      else if(iqtaim.eq.1) then
       nameaim=name(1:l)//"aim.files"
       ext="aim_"
      else if(ihirsh.eq.0) then
       nameaim=name(1:l)//"fuz.files"
       ext="fuz_"
      else if(ihirsh.eq.1) then
       nameaim=name(1:l)//"hir.files"
       ext="hir_"
      else if(ihirsh.eq.2) then
       nameaim=name(1:l)//"ihi.files"
       ext="ihi_"
      end if
      open(file=nameaim,unit=naim3,status="unknown")
      rewind(naim3)
      do jjat=1,icuat
        if(icuat.ne.nat) then
         jat=iatlist(jjat)
        else
         jat=jjat
        end if
       read(mend(iznuc(jat)),'(A2)')charnu  
         charnu=adjustl(charnu)
        l1=len_trim(name)
        l2=len_trim(charnu)
        if(jat.lt.10) then 
         write(charnu1,'(i1)')jat
         nameaim=name(1:l1)//ext//trim(charnu)//charnu1
      else if(jat.lt.100) then
         write(charnu2,'(i2)')jat
         nameaim=name(1:l1)//ext//trim(charnu)//charnu2         
      else
         write(charnu3,'(i3)')jat
         nameaim=name(1:l1)//ext//trim(charnu)//charnu3         
c assuming up to 999 atoms
      end if
       j=len(nameaim)
       do i=1,j
        if(nameaim(i:i).eq.' ') then
         l=i-1
         go to 10
        end if
       end do
  10   continue
       nameaim=nameaim(1:l)//".int"
       write(naim3,'(a40)') nameaim
      open(file=nameaim,unit=naim,status="unknown")
      rewind(naim)

c     writting stuff will be looked for by FCALC
        BK=1.0d0
        AE=-1.0d0
        write(naim,'(a20)') 'Created by APOST3D '
        if(iopt(5).eq.1) then
         write(naim,*) 'Using Mulliken atomic definition '
        else if(iopt(5).ge.2) then
         write(naim,*) 'Using Lowdin atomic definition '
        else if(iopt(6).eq.1) then
         write(naim,*) 'Using Hirshfeld atomic definition '
        else if(iopt(6).eq.2) then
         write(naim,*) 'Using Hirshfeld-Iterative atomic definition '
        else if(iopt(14).eq.1.and.iopt(31).eq.0) then
         write(naim,*) 'Using Becke-rho atomic definition '
        else if(iopt(14).eq.1.and.iopt(31).eq.1) then
         write(naim,*) 'Using TFVC atomic definition '
        else if(iopt(16).eq.1) then
         write(naim,*) 'Using QTAIM atomic definition '
        else
         write(naim,*) 'Using Becke atomic definition '
        end if
         if(iopt(5).eq.0) write(naim,*) 'Stiffness parameter k = ',iopt(25)
        if(inato.eq.1)  then
         write(naim,*) 'Correlated wave function '
         if(icorr.eq.0) then
          write(naim,*) '*Warning*, no RDM1 provided.' 
          write(naim,*) 'Can not do alpha and beta populations separatedly.' 
         end if
        else
         write(naim,*) 'Single-determinant wave function '
         if(kop.eq.1) write(naim,*) 'Unrestricted wave function '
        end if

      write(naim,'(a30,F20.11)')' MOLECULAR SCF ENERGY (AU)  = ',AE
      write(naim,*)''
      write(naim,'(A25,A4,i5)')' INTEGRATION IS OVER ATOM',mend(iznuc(jat)),
     1  jat
      !write(naim2) mend(iznuc(jat)),jat
        write(naim,'(a27)') ' RESULTS OF THE INTEGRATION'
        write(naim,'(a17,e21.14,a14,e21.14)')'              N  ',
     1   qat(jjat,1),'    NET CHARGE',-qat(jjat,1)+iznuc(jjat)
        !write(naim2) qat(jjat,1)
        write(naim,*)'             G'
        write(naim,'(a17,e21.14,a17,e21.14)')'              K  ',BK,
     1  '        E(ATOM)  ',AE
        write(naim,'(a17,e21.14)')'              L  ',0.0d0
      write(naim,*)''
        write(naim,'(a35)')'          The Atomic Overlap Matrix'
        write(naim,*)''

        if(inato.eq.1) then
         write(naim,'(a36)')'  Correlated  Wavefunction'
         write(naim,*)''
        else if(kop.eq.1) then
         write(naim,'(a36)')'Unrestricted  Wavefunction'
         write(naim,*)''
        else
         write(naim,'(a36)')'Restricted Closed-Shell Wavefunction'
         write(naim,*)''
        end if

        do i=1,2*igr
         do j=1,2*igr
           c3(i,j)=0.0d0
         end do
        end do

c transform  sat matrix in AO basis
      numorb=nalf
      if(inato.eq.1) numorb=ndim
      do i=1,numorb
       do j=1,igr
        xx=0.0d0
        do k=1,igr
         xx=xx+scr(k,i)*sat(k,j,jjat)
        end do
        c2(i,j)=xx
       end do
      end do
      do i=1,numorb
       do j=1,numorb
        xx=0.0d0
        do k=1,igr
         xx=xx+c2(i,k)*scr(k,j)       
        end do
        c3(i,j)=xx
       end do
        end do
C open-shell
        if(kop.eq.1.and.inato.eq.0)  then
      do i=1,nb    
       do j=1,igr
        xx=0.0d0
        do k=1,igr
         xx=xx+cb(k,i)*sat(k,j,jjat)
        end do
        c2(i,j)=xx
       end do
        end do
      do i=1,nb   
       do j=1,nb   
        xx=0.0d0
        do k=1,igr
         xx=xx+c2(i,k)*cb(k,j)       
        end do
        c3(i+nalf,j+nalf)=xx
       end do
        end do
        end if
        
c         do i=1,igr 
c          write(naim,*) (sat(i,j,jat),j=1,i)
c         end do

      do i=1,igr
       do j=1,igr
        csave(i,j)=csave(i,j)+c3(i,j)
       end do
      end do

c only occupied
c         do i=1,igr0
c unrestricted single-determinant
         if(kop.eq.1.and.inato.eq.0) then
          if (imulli.ne.1) then
          do i=1,nalf+nb 
           write(naim,*) (c3(i,j),j=1,i)
          end do
          else
          do i=1,nalf+nb 
           write(naim,*) (c3(i,j),j=1,nalf+nb)
          end do
          end if
c recalculate spin populations
          qalf=0.0d0
          do i=1,nalf 
            qalf=qalf+c3(i,i)
          end do
          qbet=0.0d0
          do i=1,nb   
            qbet=qbet+c3(i+nalf,i+nalf)
          end do
          write(naim,*) '  '               
          write(naim,'(a41,e21.14)') 'ALPHA ELECTRONS (NA)',qalf     
          write(naim,'(a41,e21.14)') 'BETA ELECTRONS (NB)',qbet        
c restricted single-determinant
         else if(inato.eq.0) then
          if (imulli.ne.1) then
          do i=1,nocc
           write(naim,*) (c3(i,j),j=1,i)
          end do
          else
          do i=1,nocc
           write(naim,*) (c3(i,j),j=1,nocc)
          end do
          end if
          write(naim,*) '  '               
          write(naim,'(a41,e21.14)') 'ALPHA ELECTRONS (NA)',qat(jjat,1)/2.0d0
          write(naim,'(a41,e21.14)') 'BETA ELECTRONS (NB)',qat(jjat,1)/2.0d0
c correlated wave function        
         else if(inato.eq.1) then
          if (imulli.ne.1) then
          do i=1,numorb
           write(naim,*) (c3(i,j),j=1,i)
          end do
          else
          do i=1,numorb
           write(naim,*) (c3(i,j),j=1,numorb)
          end do
          end if
          write(naim,*) '  '               
          write(naim,'(a41,e21.14)') 'ALPHA ELECTRONS (NA)',qat(jjat,1)/2.0d0
          write(naim,'(a41,e21.14)') 'BETA ELECTRONS (NB)',qat(jjat,1)/2.0d0
         end if


        write(naim,*)' '
       write(naim,*) 'NORMAL TERMINATION OF PROAIMV'

      close(naim)
        end do
      close(naim3)


        if(icuat.eq.nat) then
c check for orthogonality of the MOs
      xmax=1.0d-2
      xmaxd=1.0d-3
      xmaxns=1.0d-4
c         do i=1,igr
c          write(*,*) (csave(i,j),j=1,igr)
c         end do
c      xmaxg=1.0d-2
      do i=1,igr
       do j=i,igr
          if(i.eq.j) then
         if(abs(csave(i,i))-1.0d0.gt.xmax) then 
          xx=abs(csave(i,i))-1.0d0
          imaxd=i
            if(i.le.nalf) then
          write(*,*) 'Dev. from normalization ',xx,imaxd
            end if
         end if
        else 
         if(abs(csave(i,j)).gt.xmaxd) then 
          xx=abs(csave(i,j))
          imax1=i
          imax2=j
            if(i.le.nalf.and.j.le.nalf) then
          write(*,'(a27,f10.6,2i3)') 'Dev. from orth. in occ. set:',xx,
     1      imax1,imax2
            end if
         end if
         if(j.gt.i) then
          xxx=abs(csave(i,j))-abs(csave(j,i)) 
          if(xxx.gt.xmaxns) then 
           xx=xxx
           imaxns1=i   
           imaxns2=j   
           write(*,*) 'Dev. from hermiticity',xx,imaxns1,imaxns2
          end if
         end if
        end if
       end do
      end do
      end if
      deallocate(c3,scr,c2,csave)

      end

CCCCC
C FOR WRITTING CUBE FILES
CCCCC

      subroutine cubegen3(ifrag,icase)
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /iops/iopt(100)
      character*30 name
      common /filename/name

      character*30 name2
      character nameaim*60, charnu*2,atnu*2,charnu1*3
      character*2 mend(92)
      data mend/' H','He','Li','Be',' B',' C',' N',' O',
     $ ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl',
     $ 'Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe',
     $ 'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
     $ 'Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     $ 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I',
     $ 'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sn',
     $ 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
     $ 'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     $ 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     $ 'Th','Pa',' U'  /

      dimension pgrid(maxgrid,3),igrid(3),xgrid(3,3)
      dimension pop(maxat)          
      allocatable c0(:,:),xyz(:,:,:)

      ihirsh=iopt(6)
      imulli= Iopt(5) 
      ibcp=Iopt(14)
      iqtaim = Iopt(16)  
      inewbec = Iopt(31) 
      idofr=Iopt(40)
      jcubthr=iopt(41)
      kcubthr=iopt(42)

      allocate(c0(igr,igr))
      imaxo=ip0(ifrag)
      do i=1,igr
       do j=1,imaxo
        c0(i,j)=p0(i,j)
       end do
      end do

c setting actual effos to print, instead
      if(jcubthr.lt.0) then
       imaxeff=abs(jcubthr)
       imineff=abs(kcubthr)
      else
       imaxeff=0
       xmaxeff=float(jcubthr)*1.0d-3
       if (icase.eq.0) xmaxeff=2.0d0*xmaxeff
1      imaxeff= imaxeff+1
       if(p0net(imaxeff,ifrag).ge.xmaxeff) go to 1
       imineff=imaxo+1
       xmineff=float(kcubthr)*1.0d-3   
       if (icase.eq.0) xmineff=2.0d0*xmineff
2      imineff= imineff - 1
       if(p0net(imineff,ifrag).le.xmineff) go to 2
      end if
c

      if(imaxeff.gt.imineff) then
       write(*,*) 'No eff-AOs in the occupation range' 
       return
      end if
      write(*,22)'Generating cube files for eff-AOs',imaxeff,' to',imineff,' of atom/fragment ',ifrag
22    format (a33,i3,a3,i3,a18,i3)

        if (imulli.eq.1) then
         name2="mulliken"
        else if (imulli.gt.1) then
         name2="lowdin"
        else
        if (ihirsh.eq.0) then
         name2="becke"
         if(ibcp.eq.1) then 
            if(inewbec.eq.0)then
              name2="beckerho"
            else
             name2="tfvc"
            end if
         end if
         if(iqtaim.eq.1) name2="qtaim"
        else if(ihirsh.eq.1) then
         name2="hirsh"
        else if(ihirsh.eq.2) then
         name2="hirsh-it"
         do i=1,nat
          pop(i)=qat(i,1)
         end do
        end if
        end if

C grid points for cube in each dimension
      ngridp=40 

      do i=1,3
       igrid(i)=ngridp
c       if (i.eq.3) igrid(i)=100
c asuming rectangular grid...
       do j=1,3
       xgrid(i,j)=0.0d0                 
       end do
      end do

C Two options: For moleculs/fragments or centered on atoms      
C now active first option as it deals with fragments

      if(1.eq.1) then
      extra=3.2
      volume=1.0d0
c furthest x y z atomic posisiotns of the fragment
      do i=1,3
       xmax=-1.0d8
       xmin=1.0d8
       do jatom=1,nfrlist(ifrag)
        iatom=ifrlist(jatom,ifrag)
        if(coord(i,iatom).lt.xmin) xmin=coord(i,iatom)
        if(coord(i,iatom).gt.xmax) xmax=coord(i,iatom)
       end do
       xmin=xmin-extra
       xmax=xmax+extra
       dist=xmax-xmin
       volume=volume*dist
       xgrid(i,i)=dist/(igrid(i)-1.0d0)
       do j=1,igrid(i)
        pgrid(j,i)=xmin+(j-1)*xgrid(i,i)
       end do
      end do

c centering the grid on the atom
      else
      extra=3.0
      iatom=jjat  
      do i=1,3
       xmin=coord(i,iatom)
       xmax=coord(i,iatom)
       xmin=xmin-extra
       xmax=xmax+extra
       dist=xmax-xmin
       xgrid(i,i)=dist/(igrid(i)-1.0d0)
       do j=1,igrid(i)
        pgrid(j,i)=xmin+(j-1)*xgrid(i,i)
       end do
      end do
      end if

c Now the grid
      allocate ( xyz(igrid(1),igrid(2),igrid(3)))

      do ivec=imaxeff,imineff
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          xabs=pgrid(i,1)
          yabs=pgrid(j,2)
          zabs=pgrid(k,3)

          if(imulli.ne.0) then
           xyz(i,j,k)=orbxyz(c0,ivec,xabs,yabs,zabs)
          else
           ww=0.0d0
           do iatom=1,nfrlist(ifrag)
           jjat=ifrlist(iatom,ifrag)
           if(ihirsh.eq.0.and.iqtaim.eq.0) then
            ww=ww+wat(jjat,xabs,yabs,zabs)
           else if(ihirsh.eq.1) then
            ww=ww+wathirsh(jjat,xabs,yabs,zabs)
           else if(ihirsh.eq.2) then
            ww=ww+wathirsh2(jjat,xabs,yabs,zabs,pop)
           else if(iqtaim.eq.1) then 
c            call cubeqtaim(jjat,xabs,yabs,zabs,ww0)
            ww=ww+ww0
           end if
           end do
           xyz(i,j,k)=orbxyz(c0,ivec,xabs,yabs,zabs)*ww
          end if
         end do
        end do
       end do
c approximate normalization of orbital
       x0=0.0d0
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          x0=x0+xyz(i,j,k)*xyz(i,j,k)
         end do
        end do
       end do
       write(*,'(a25,f7.4)') 'Normalization from cube: ',x0*volume/(ngridp**3.0d0)

c      OUTPUT    

        if(idofr.eq.0) then
        read(mend(iznuc(ifrag)),'(A2)')charnu
         charnu=adjustl(charnu)
        else
         charnu="FR"
        end if   
        l0=len_trim(name)
        l1=len_trim(name2)
        l2=len_trim(charnu)
c assuming up to 99 atoms
        if(ifrag.lt.10) then
           write(atnu,'(i1)')ifrag
        else
           write(atnu,'(i2)')ifrag
        end if
        if(ivec.lt.10) then
         write(charnu1,'(i1)')ivec
        else if (ivec.lt.100) then
         write(charnu1,'(i2)')ivec
        else
         write(charnu1,'(i3)')ivec
        end if
         if(icase.ne.2) then
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)
         else                
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)//"beta"
         end if 
       j=len(nameaim)
       do i=1,j
        if(nameaim(i:i).eq.' ') then
         llen=i-1
         go to 10
        end if
       end do
  10   continue
       nameaim=nameaim(1:llen)//".cube"

       open(44,file=nameaim,status="unknown")
       rewind(44)
       write(44,*)'Cube generated with APOST-3D code '
       if(idofr.eq.0) then
       write(44,41) trim(name),name2,' EFFAO',ivec," for atom",
     + mend(iznuc(ifrag)),"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       else
       write(44,44) trim(name),name2,' EFFAO',ivec," for frag",
     + ifrag,"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       end if
       write(44,42) nat,(pgrid(1,j),j=1,3)
       do i=1,3
        write(44,42) igrid(i),(xgrid(i,j),j=1,3)
       end do
       do i=1,nat
        write(44,43) iznuc(i),zn(i),(coord(j,i),j=1,3)
       end do
       ione=1
       write(44,'(2i5)')ione,ione
       do i=1,igrid(1)
       do j=1,igrid(2)
        write(44,40)(xyz(i,j,k),k=1,igrid(3))
       end do
       end do
       close(44)
      end do
41    format(a8,x,a8,a7,i3,a9,a2,a11,f8.4,a9,f8.4)
44    format(a8,x,a8,a7,i3,a9,i2,a11,f8.4,a9,f8.4)
42    format(i5,3f12.6)
43    format(i5,4f12.6)
40    format(6e13.5)
     
      deallocate(c0,xyz)
      end


CCCCC
C FOR WRITTING FCHK FILES
CCCCC
       subroutine rmat(iunit,key,ival,jval,ndim,rmatrix)
       implicit double precision (a-h,o-z)
       include 'parameter.h'
       dimension rmatrix(ndim,ndim) 
       character(len=*) ::  key
       integer ival,jval
       character(len=43) ::  title
       title=adjustl(key)
       write(iunit,'(A43,A6,I12)') title,"R   N=",ival*jval
       write(iunit,'(5ES16.8)') ((rmatrix(i,j),i=1,jval),j=1,ival)
       end

       subroutine rarr(iunit,key,ival,ndim,rarray)
       implicit double precision (a-h,o-z)
       include 'parameter.h'
       character(len=*) ::  key
       integer ival
       dimension rarray(ndim)
       character(len=43) ::  title
       title=adjustl(key)
       write(iunit,'(A43,A6,I12)') title,"R   N=",ival
       write(iunit,'(5ES16.8)') (rarray(i),i=1,ival)
       end

       subroutine ival(iunit,key,ivalue)
       implicit double precision (a-h,o-z)
       character(len=*) ::  key
       integer ivalue
       character(len=43) ::  title
       title=adjustl(key)
       write(iunit,'(A43,A,I17)') title,"I",ivalue
       end

        SUBROUTINE PRINTMAT(N_orbital,S)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        include 'parameter.h'
        DIMENSION S(N_orbital,N_orbital)

        nblock=N_orbital/5
        if(nblock*5.ne.N_orbital) nblock=nblock+1
        do k=1,nblock
          ii=min0(N_orbital,5*k)
          write(*,'(4x,10(7X,i6,a1))') (j," ",j=5*(k-1)+1,ii)
          do i=5*(k-1)+1,N_orbital
!            ii=5*k
!            if(i.lt.5*k) ii=i
            ii=min0(i,5*k)
            write(*,'(i7,5(x,d13.6))') i,(S(i,j),j=5*(k-1)+1,ii)
          end do
        end do
        END SUBROUTINE

      subroutine cubegen3_mhg(ifrag,icase)
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /iops/iopt(100)
      character*30 name
      common /filename/name

      character*30 name2
      character nameaim*60, charnu*2,atnu*2,charnu1*3
      character*2 mend(92)
      data mend/' H','He','Li','Be',' B',' C',' N',' O',
     $ ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl',
     $ 'Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe',
     $ 'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
     $ 'Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     $ 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I',
     $ 'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sn',
     $ 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
     $ 'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     $ 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     $ 'Th','Pa',' U'  /

      dimension pgrid(maxgrid,3),igrid(3),xgrid(3,3)
      dimension pop(maxat)          
      allocatable c0(:,:),xyz(:,:,:)

      ihirsh=iopt(6)
      imulli= Iopt(5) 
      ibcp=Iopt(14)
      iqtaim = Iopt(16)  
      inewbec = Iopt(31) 
      idofr=Iopt(40)
      jcubthr=iopt(41)
      kcubthr=iopt(42)

      allocate(c0(igr,igr))
      imaxo=ip0(ifrag)
      do i=1,igr
       do j=1,imaxo
        c0(i,j)=p0(i,j)
       end do
      end do

c setting actual effos to print, instead
      if(jcubthr.lt.0) then
       imaxeff=abs(jcubthr)
       imineff=abs(kcubthr)
      else
       imaxeff=0
       xmaxeff=float(jcubthr)*1.0d-3
       if (icase.eq.0) xmaxeff=2.0d0*xmaxeff
1      imaxeff= imaxeff+1
       if(p0net(imaxeff,ifrag).ge.xmaxeff) go to 1
       imineff=imaxo+1
       xmineff=float(kcubthr)*1.0d-3   
       if (icase.eq.0) xmineff=2.0d0*xmineff
2      imineff= imineff - 1
       if(p0net(imineff,ifrag).le.xmineff) go to 2
      end if
c

      if(imaxeff.gt.imineff) then
       write(*,*) 'No eff-AOs in the occupation range' 
       return
      end if
      write(*,22)'Generating cube files for eff-AOs',imaxeff,' to',imineff,' of atom/fragment ',ifrag
22    format (a33,i3,a3,i3,a18,i3)

      name2="mhg"

C grid points for cube in each dimension
      ngridp=40 

      do i=1,3
       igrid(i)=ngridp
c       if (i.eq.3) igrid(i)=100
c asuming rectangular grid...
       do j=1,3
       xgrid(i,j)=0.0d0                 
       end do
      end do

C Two options: For moleculs/fragments or centered on atoms      
C now active first option as it deals with fragments

      if(1.eq.1) then
      extra=5.2
      volume=1.0d0
c furthest x y z atomic posisiotns of the fragment
      do i=1,3
       xmax=-1.0d8
       xmin=1.0d8
       do jatom=1,nfrlist(ifrag)
        iatom=ifrlist(jatom,ifrag)
        if(coord(i,iatom).lt.xmin) xmin=coord(i,iatom)
        if(coord(i,iatom).gt.xmax) xmax=coord(i,iatom)
       end do
       xmin=xmin-extra
       xmax=xmax+extra
       dist=xmax-xmin
       volume=volume*dist
       xgrid(i,i)=dist/(igrid(i)-1.0d0)
       do j=1,igrid(i)
        pgrid(j,i)=xmin+(j-1)*xgrid(i,i)
       end do
      end do

c centering the grid on the atom
      else
      extra=3.0
      iatom=jjat  
      do i=1,3
       xmin=coord(i,iatom)
       xmax=coord(i,iatom)
       xmin=xmin-extra
       xmax=xmax+extra
       dist=xmax-xmin
       xgrid(i,i)=dist/(igrid(i)-1.0d0)
       do j=1,igrid(i)
        pgrid(j,i)=xmin+(j-1)*xgrid(i,i)
       end do
      end do
      end if

c Now the grid
      allocate ( xyz(igrid(1),igrid(2),igrid(3)))

      do ivec=imaxeff,imineff
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          xabs=pgrid(i,1)
          yabs=pgrid(j,2)
          zabs=pgrid(k,3)
           xyz(i,j,k)=orbxyz(c0,ivec,xabs,yabs,zabs)
         end do
        end do
       end do
c approximate normalization of orbital
       x0=0.0d0
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          x0=x0+xyz(i,j,k)*xyz(i,j,k)
         end do
        end do
       end do
       write(*,'(a25,f7.4)') 'Normalization from cube: ',x0*volume/(ngridp**3.0d0)

c      OUTPUT    

        if(idofr.eq.0) then
        read(mend(iznuc(ifrag)),'(A2)')charnu
         charnu=adjustl(charnu)
        else
         charnu="FR"
        end if   
        l0=len_trim(name)
        l1=len_trim(name2)
        l2=len_trim(charnu)
c assuming up to 99 atoms
        if(ifrag.lt.10) then
           write(atnu,'(i1)')ifrag
        else
           write(atnu,'(i2)')ifrag
        end if
        if(ivec.lt.10) then
         write(charnu1,'(i1)')ivec
        else if (ivec.lt.100) then
         write(charnu1,'(i2)')ivec
        else
         write(charnu1,'(i3)')ivec
        end if
         if(icase.ne.2) then
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)
         else                
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)//"beta"
         end if 
       j=len(nameaim)
       do i=1,j
        if(nameaim(i:i).eq.' ') then
         llen=i-1
         go to 10
        end if
       end do
  10   continue
       nameaim=nameaim(1:llen)//".cube"

       open(44,file=nameaim,status="unknown")
       rewind(44)
       write(44,*)'Cube generated with APOST-3D code '
       if(idofr.eq.0) then
       write(44,41) trim(name),name2,' EFFAO',ivec," for atom",
     + mend(iznuc(ifrag)),"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       else
       write(44,44) trim(name),name2,' EFFAO',ivec," for frag",
     + ifrag,"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       end if
       write(44,42) nat,(pgrid(1,j),j=1,3)
       do i=1,3
        write(44,42) igrid(i),(xgrid(i,j),j=1,3)
       end do
       do i=1,nat
        write(44,43) iznuc(i),zn(i),(coord(j,i),j=1,3)
       end do
       ione=1
       write(44,'(2i5)')ione,ione
       do i=1,igrid(1)
       do j=1,igrid(2)
        write(44,40)(xyz(i,j,k),k=1,igrid(3))
       end do
       end do
       close(44)
      end do
41    format(a8,x,a8,a7,i3,a9,a2,a11,f8.4,a9,f8.4)
44    format(a8,x,a8,a7,i3,a9,i2,a11,f8.4,a9,f8.4)
42    format(i5,3f12.6)
43    format(i5,4f12.6)
40    format(6e13.5)
     
      deallocate(c0,xyz)
      end

      subroutine cubegen4(ifrag,icase)
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /atomrad/atr(maxat),dist(maxat,maxat)
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common/effao/p0(nmax,nmax),p0net(nmax,maxat),p0gro(nmax,maxat),ip0(maxat)
      common /qat/qat(maxat,2),qsat(maxat,2)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)
      common /iops/iopt(100)
      character*30 name
      common /filename/name

      character*30 name2
      character nameaim*60, charnu*2,atnu*2,charnu1*3
      character*2 mend(92)
      data mend/' H','He','Li','Be',' B',' C',' N',' O',
     $ ' F','Ne','Na','Mg','Al','Si',' P',' S','Cl',
     $ 'Ar',' K','Ca','Sc','Ti',' V','Cr','Mn','Fe',
     $ 'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',
     $ 'Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     $ 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I',
     $ 'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sn',
     $ 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',
     $ 'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     $ 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     $ 'Th','Pa',' U'  /

      dimension pgrid(maxgrid,3),igrid(3),xgrid(3,3)
      dimension pop(maxat)          
      allocatable c0(:,:),xyz(:,:,:)

      ihirsh=iopt(6)
      imulli= Iopt(5) 
      ibcp=Iopt(14)
      iqtaim = Iopt(16)  
      inewbec = Iopt(31) 
      idofr=Iopt(40)
      jcubthr=iopt(41)
      kcubthr=iopt(42)


      imaxo=ip0(ifrag)
      ALLOCATE(c0(igr,igr))
      do i=1,igr
        do j=1,imaxo
          c0(i,j)=p0(i,j)
        end do
      end do
c setting actual effos to print, instead
      if(jcubthr.lt.0) then
       imaxeff=abs(jcubthr)
       imineff=abs(kcubthr)
      else
       imaxeff=0
       xmaxeff=float(jcubthr)*1.0d-3
       if (icase.eq.0.or.icase.eq.3) xmaxeff=2.0d0*xmaxeff

1      imaxeff= imaxeff+1
       if(p0net(imaxeff,ifrag).ge.xmaxeff) go to 1
       imineff=imaxo+1
       xmineff=float(kcubthr)*1.0d-3   
       if (icase.eq.0.or.icase.eq.3) xmineff=2.0d0*xmineff

2      imineff= imineff - 1
       if(p0net(imineff,ifrag).le.xmineff) go to 2
      end if
c

      if(imaxeff.gt.imineff) then
       write(*,*) ' No eff-AOs in the occupation range' 
       write(*,*) " "
       deallocate(c0)
       return
      end if

      write(*,22)'Generating cube files for eff-AOs',imaxeff,' to',imineff,' of atom/fragment ',ifrag
22    format (a33,i3,a3,i3,a18,i3)

        if (imulli.eq.1) then
         name2="mulliken"
        else if (imulli.gt.1) then
         name2="lowdin"
        else
        if (ihirsh.eq.0) then
         name2="becke"
         if(ibcp.eq.1) then 
            if(inewbec.eq.0)then
              name2="beckerho"
            else
             name2="tfvc"
            end if
         end if
         if(iqtaim.eq.1) name2="qtaim"
        else if(ihirsh.eq.1) then
         name2="hirsh"
        else if(ihirsh.eq.2) then
         name2="hirsh-it"
         do i=1,nat
          pop(i)=qat(i,1)
         end do
        end if
        end if
        if(icase.eq.3) name2=trim(name2)//"_paired"
        if(icase.eq.4) name2=trim(name2)//"_unpaired"

c asuming rectangular grid...
       xgrid=0.0d0                 

C For moleculs/fragments 
c adaptative size cube
c extra set to 3 times the atomic radii
      xmesh=0.2d0
      rrmax=3.0d0
      volume=1.0d0
c furthest x y z atomic positions of the fragment
      do i=1,3
       xmax=-1.0d8
       xmin=1.0d8
       do jatom=1,nfrlist(ifrag)
        iatom=ifrlist(jatom,ifrag)
        if(coord(i,iatom).lt.xmin) then
         xmin=coord(i,iatom)
         iiatom=iatom
        end if
        if(coord(i,iatom).gt.xmax) then
         xmax=coord(i,iatom)
         iiiatom=iatom
        end if
       end do
       xmin=xmin-rrmax*atr(iiatom)
       xmax=xmax+rrmax*atr(iiiatom)
       dist0=xmax-xmin
       xgrid(i,i)=xmesh
       igrid(i)=int(dist0/xmesh)+1
       do j=1,igrid(i)
        pgrid(j,i)=xmin+(j-1)*xgrid(i,i)
       end do
       volume=volume*dist0
      end do

c Now the grid
      write(*,'(a21,3i4)')'Size of cube (x,y,z):',(igrid(i),i=1,3)
      allocate ( xyz(igrid(1),igrid(2),igrid(3)))

      do ivec=imaxeff,imineff
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          xabs=pgrid(i,1)
          yabs=pgrid(j,2)
          zabs=pgrid(k,3)

          if(imulli.ne.0) then
           xyz(i,j,k)=orbxyz(c0,ivec,xabs,yabs,zabs)
          else
           ww=0.0d0
           do iatom=1,nfrlist(ifrag)
            jjat=ifrlist(iatom,ifrag)
            if(ihirsh.eq.0.and.iqtaim.eq.0) then
              ww=ww+wat(jjat,xabs,yabs,zabs)
            else if(ihirsh.eq.1) then
              ww=ww+wathirsh(jjat,xabs,yabs,zabs)
            else if(ihirsh.eq.2) then
              ww=ww+wathirsh2(jjat,xabs,yabs,zabs,pop)
            else if(iqtaim.eq.1) then 
c             call cubeqtaim(jjat,xabs,yabs,zabs,ww0)
c              ww=ww+ww0
            end if
           end do
           xyz(i,j,k)=orbxyz(c0,ivec,xabs,yabs,zabs)*ww
          end if
         end do
        end do
       end do
c approximate normalization of orbital
       x0=0.0d0
       do i=1,igrid(1)
        do j=1,igrid(2)
         do k=1,igrid(3)
          x0=x0+xyz(i,j,k)*xyz(i,j,k)
         end do
        end do
       end do
       write(*,'(a25,f7.4)') 'Normalization from cube: ',x0*volume/(igrid(1)*igrid(2)*igrid(3))
       write(*,*) " "

c      OUTPUT    

        if(idofr.eq.0) then
        read(mend(iznuc(ifrag)),'(A2)')charnu
         charnu=adjustl(charnu)
        else
         charnu="FR"
        end if   
        l0=len_trim(name)
        l1=len_trim(name2)
        l2=len_trim(charnu)
c assuming up to 99 atoms
        if(ifrag.lt.10) then
           write(atnu,'(i1)')ifrag
        else
           write(atnu,'(i2)')ifrag
        end if
        if(ivec.lt.10) then
         write(charnu1,'(i1)')ivec
        else if (ivec.lt.100) then
         write(charnu1,'(i2)')ivec
        else
         write(charnu1,'(i3)')ivec
        end if
         if(icase.ne.2) then
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)
         else                
         nameaim=trim(name)//"_"//trim(name2)//"_"//trim(charnu)//
     +   trim(atnu)//"_"//trim(charnu1)//"beta"
         end if 
       j=len(nameaim)
       do i=1,j
        if(nameaim(i:i).eq.' ') then
         llen=i-1
         go to 10
        end if
       end do
  10   continue
       nameaim=nameaim(1:llen)//".cube"

       open(44,file=nameaim,status="unknown")
       rewind(44)
       write(44,*)'Cube generated with APOST-3D code '
       if(idofr.eq.0) then
       write(44,41) trim(name),name2,' EFFAO',ivec," for atom",
     + mend(iznuc(ifrag)),"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       else
       write(44,44) trim(name),name2,' EFFAO',ivec," for frag",
     + ifrag,"Gross Occ.",p0gro(ivec,ifrag),"Net Occ.",
     + p0net(ivec,ifrag)
       end if
       write(44,42) nat,(pgrid(1,j),j=1,3)
       do i=1,3
        write(44,42) igrid(i),(xgrid(i,j),j=1,3)
       end do
       do i=1,nat
        write(44,43) iznuc(i),zn(i),(coord(j,i),j=1,3)
       end do
       ione=1
       write(44,'(2i5)')ione,ione
       do i=1,igrid(1)
       do j=1,igrid(2)
        write(44,40)(xyz(i,j,k),k=1,igrid(3))
       end do
       end do
       close(44)
      end do
41    format(a8,x,a8,a7,i3,a9,a2,a11,f8.4,a9,f8.4)
44    format(a8,x,a8,a7,i3,a9,i2,a11,f8.4,a9,f8.4)
42    format(i5,3f12.6)
43    format(i5,4f12.6)
40    format(6e13.5)
     
      deallocate(c0,xyz)
      end

