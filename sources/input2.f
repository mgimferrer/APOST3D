      subroutine input  
      use basis_set
      use ao_matrices
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord0(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common /energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      dimension clin(nmax**2)
      character*80 line
      logical ilog
      real*8, allocatable ::vv(:,:)

      ndens0=iopt(9)
     
      irhf=0
      iuhf=0
      icas=0
      icisd=0
      ikop=0
      irohf=0

ccccccccccccccccccccccccccccc
c start processing fchk file
ccccccccccccccccccccccccccccc
      rewind 15
      read(15,'(a80)')line
      read(15,'(a80)')line
      if(index(line(11:11),"R").ne.0) then
       irhf=1
       write(*,*) 'Restricted MO formalism'
      else if(index(line(11:11),"U").ne.0) then
       write(*,*) 'Unrestricted MO formalism'
       iuhf=1
      end if
      if(index(line(11:25),"CASSCF").ne.0) then
       icas=1
       write(*,*) 'CASSCF type calculation'
      else if(index(line(11:25),"HF").ne.0) then
       write(*,*) 'Hartree-Fock type calculation'
      else if(index(line(11:25),"MP").ne.0) then
       write(*,*) 'Moller-Plesset PT type calculation'
      else if(index(line(11:25),"CI").ne.0) then
       icisd=1
       write(*,*) 'CI type calculation'
      else if(index(line(11:25),"CCSD").ne.0) then
       write(*,*) 'Coupled-cluster type calculation'
       icisd=1
      else 
       write(*,*)'KS-DFT type calculation'
       write(*,*)'(single det. built up of the Kohn-Sham orbitals)'
      end if

      call build_basis()

      igr=nbasis
      nat=natoms

      call build_ao_matrices(igr) !MMO- WIP!

      igr0=int_locate(15,"Number of independ",ilog)
      if(igr.ne.igr0) write(*,*) 'WARNING, some basis functions have bee
     +n removed'
   

      nelectr=int_locate(15,"Number of electr",ilog)

      nalf=int_locate(15,"Number of alpha electr",ilog)

      nb=int_locate(15,"Number of beta electr",ilog)

      dummy=int_locate(15,"Atomic numbers",ilog)
      read(15,*)(iznuc(i),i=1,nat)   

      dummy=int_locate(15,"Nuclear charg",ilog)
      read(15,*)(zn(i),i=1,nat)

      dummy=int_locate(15,"Current cartesian",ilog)
      read(15,*)(coord0(1,i),coord0(2,i),coord0(3,i),i=1,nat)

      escf=real_locate(15,"SCF Energ",ilog)
      ss2=real_locate(15,"S**2    ",ilog)

      dummy=int_locate(15,"Alpha MO co",ilog)
      read(15,*)((c(i,j),i=1,igr),j=1,igr0)
      if(igr0.ne.igr) then
       do i=igr0+1,igr
        do j=1,igr
         c(j,i)=0.0d0
        end do
       end do
      end if

      dummy=int_locate(15,"Beta MO co",ilog)
      if(ilog) then
       read(15,*)((cb(i,j),i=1,igr),j=1,igr0)
       if(igr0.ne.igr) then
        do i=igr0+1,igr
         do j=1,igr
          c(j,i)=0.0d0
         end do
        end do
       end if
c just in case not detected above
       iuhf=1
       kop=1
      end if


C Reading P-matrix from fchk
      ntriang=int_locate(15,"otal SCF D",ilog)
      read(15,*)(clin(i),i=1,ntriang)
      ndens=1     

      if(ndens0.gt.1) then
       ndens=-1
       rewind 15
 200   read(15,'(a80)',end=211)line
       if(index(line,"Total").ne.0) ndens=ndens+1
       goto 200
 211   continue                 

       if(ndens.ne.1) then
        print *, ndens,' densities found in the fchk file'   
       else
        print *,' There is only a P-matrix the fchk file - it is used.' 
       end if
       if(ndens0.le.ndens)then
        ndens=ndens0
        print *,'Using the ',ndens,'-d/th density in the fchk file'
       else 
        print *,ndens0,'-d/th density not found in the fchk file'
        stop
       endif

c Getting the n-th Density matrix in fchk instead
       ncou=-1
       rewind 15
       do while (ncou.lt.ndens)
        read(15,'(a80)')line
        if(index(line,"Total").ne.0) ncou=ncou+1   
       end do
       read(15,*)(clin(i),i=1,ntriang)
      end if
c 
      ncou=1
      do i=1,igr
       do j=1,i
        p(i,j)=clin(ncou)
        p(j,i)=clin(ncou)
        ncou=ncou+1
       enddo
      enddo

c checking for RO case. Wrong P provided in FChk in some Gaussian
c versions. Building P and PS from MOs just in case
      irohf=int_locate(15,"IROHF",ilog)
      if(irohf.eq.1) then
       kop=1
       do i=1,igr
        do j=1,igr
         cb(i,j)=c(i,j)
        enddo
       enddo
      end if
      if(kop.eq.1) then
       do i=1,igr
        do j=1,igr
         pa(i,j)=0.0d0
         pb(i,j)=0.0d0
         do ij=1,nalf
          pa(i,j)=pa(i,j)+c(i,ij)*c(j,ij)
         end do
         do ij=1,nb
           pb(i,j)=pb(i,j)+cb(i,ij)*cb(j,ij)
         end do
         p(i,j)=pa(i,j)+pb(i,j)
         ps(i,j)=pa(i,j)-pb(i,j)
        end do
       end do
      else
       do i=1,igr
        do j=1,igr
         pa(i,j)=p(i,j)/2.0d0 
         pb(i,j)=p(i,j)/2.0d0 
        enddo
       enddo
      end if
c
c checking if spin density is available in fchk
c
      ntriang=int_locate(15,"Spin SCF ",ilog)
      if(ilog) then
       read(15,*)(clin(i),i=1,ntriang)
       ncou=1
       do i=1,igr
        do j=1,i
         ps(i,j)=clin(ncou)
         ps(j,i)=clin(ncou)
         ncou=ncou+1
        enddo
       enddo
      else if(icas.ne.1.or.icisd.ne.1) then
       print *, 'Spin density not found in FChk'
       print *, 'Reconstructing spin density (if any) from MOs'
      end if
      
c In case CASSCF
      ncasel=int_locate(15,"Number of CAS El",ilog)
      if(ilog) then
       icas=1
       ncasorb=int_locate(15,"Number of CAS Or",ilog)
       print *,'CAS specification is read:'
       print *,'Active electrons : ',ncasel 
       print *,'Active orbitals  : ',ncasorb
      endif
      
      nocc=nelectr/2
      if(icas.eq.1) then
      ncore=(nelectr-ncasel)/2
      norb=ncore+ncasorb
      nspinorb=norb*2

c       nx=ncasel/2
c       if(2*nx.ne.ncasel)nx=nx+1
c       nspinorb=2*(nocc-nx+ncasorb)
c       if(2*nocc.ne.nelectr)nspinorb=nspinorb+2
       print *,'Number of CAS spin-orbitals', nspinorb
       icass=0
       if(abs(ss2).gt.1.0d-2) then
        print *,'CAS WF other than pure singlet state'
        icass=1
       end if
      else if(icisd.eq.1) then 
       nspinorb=igr0*2
       print *,' Number of CI spin-orbitals', nspinorb
      endif
c      norb=nspinorb/2

      if(icas.eq.0)then
       print *,' '
       print *,'nocc',nocc,'    nalpha',nalf,'    nbeta',nb
      else
       print *,'nalpha',nalf,'    nbeta',nb
      endif
      print *,' '
      print *,'E(SCF/DFT)=',escf
      print *,' '

ccccccccccccccccccccccccccc
c end processing fchk file
ccccccccccccccccccccccccccc
c      allocate(vv(nbasis,nbasis))
c      call do_potential(zn,vv)
c      call printmat(nbasis,vv)

       end

      subroutine dm1input(dm1)
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      character*80 line
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      dimension  dm1(nspinorb,nspinorb)

      allocatable :: scr(:,:)

      iorca    = Iopt(43)
      ipyscf   = Iopt(87)
      thres   = 1.0d-8

      allocate(scr(igr,igr))
       
!! Nbasis HERE IS ncore + nact SPIN-ORBITALS (2x ncore + nact SPATIAL MOs)  

      do i=1,nspinorb
        do j=1,nspinorb
          dm1(i,j)=ZERO
        end do
      end do

c reading dm1 as provided by DMN node. full space
      if(iorca.ne.1.and.ipyscf.ne.1) then
          rewind(11)
          do while(.true.)
            read(11,end=99) i,j,vvv 
            if(ABS(vvv).gt.thres) then
              dm1(i,j)=dm1(i,j)+vvv
              dm1(j,i)=dm1(i,j)
            end if
          end do
      else 

c reading dm1 as provided by PySCF code. Only active space
c correcting for active orbital index

        ncore=(nspinorb-ncasorb*2)
        write(*,*) "Nspinorb, Ncasorb, ncoreorb : ",nspinorb,ncasorb,ncore 
        write(*,*) 'Generating DM1 for inactive spin orbitals, ',ncore
        write(*,*) 'Assuming diagonal dm1 '

        do ii=1,ncore  
          dm1(ii,ii)=ONE
        end do
        rewind(11)
        read(11,*) line
        do while(.true.)
          read(11,*,end=99) ii,jj,vvv
          if(ABS(vvv).gt.thres) then
            i=ii+ncore
            j=jj+ncore
            dm1(i,j)=vvv
            dm1(j,i)=dm1(i,j)
          end if
        end do
      end if
99    print *, "End DM1 reading "

      difi=ZERO
      do i=1,nspinorb
        do j=1,nspinorb
          x=ZERO
          do k=1,nspinorb
            x=x+dm1(i,k)*dm1(k,j)
          end do
          difi=dmax1(difi,dabs(x-dm1(i,j)))
        end do
      end do
      print *, " "
      print *, " Max. deviation from idempotency = ",difi
      print *, " "

      do i=1,norb
        do k=1,norb
          occ_no(i,k)=dm1((i-1)*2+1,(k-1)*2+1)+dm1((i-1)*2+2,(k-1)*2+2)
        end do
      end do

!! DIAGONALIZING TO GET NOs AND TRANSFORMATION TO THE AO BASIS !!

      call diagonalize(igr,norb,occ_no,scr,0)
      do i=1,igr          
        do j=1,norb         
          x=ZERO
          do k=1,norb         
            x=x+c(i,k)*scr(k,j) 
          end do
          c_no(i,j)=x                          
        end do
      end do

c Building Pa and Pb matrix from dm1
      write(*,*) 'Reconstructing P and Ps matrices from dm1 file'
      do nu=1,igr 
       do mu=1,igr  
        x=0.d0
        do i=1,nspinorb,2 
         ii=i/2+1
         do j=1,nspinorb,2 
          jj=j/2+1
          x=x+c(nu,ii)*dm1(i,j)*c(mu,jj)
         enddo 
        enddo
        pa(mu,nu)=x    
       enddo
      enddo

      do nu=1,igr 
       do mu=1,igr  
        x=0.d0
        do i=2,nspinorb,2 
         ii=i/2
         do j=2,nspinorb,2 
          jj=j/2
          x=x+c(nu,ii)*dm1(i,j)*c(mu,jj)
         enddo 
        enddo
        pb(mu,nu)=x    
       enddo
      enddo

      do mu=1,igr 
       do nu=1,igr 
        ps(mu,nu)=pa(mu,nu)-pb(mu,nu)
       enddo
      enddo

      dif=0.0d0
      do i=1,igr 
       do j=1,igr  
        dif=amax1(dif,dabs((pa(i,j)+pb(i,j))-p(i,j)))
        p(i,j)=pa(i,j)+pb(i,j)
       end do 
      end do 
      write(*,*) 'P (and Ps) in fchk overriden'    
      write(*,*) 'Max dif. P from dm1 and from FChk: ',dif        


      deallocate(scr)
      end

! *****

      subroutine mga_misc(iecp,eecp)
      use basis_set
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /achi/achi(maxat,maxat),ibcp
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/energ0/ekin0,eelnuc0,evee0,eone0
      dimension eecp(maxat)
      character*80 line

      allocatable :: xecpv(:),xecpm(:,:),xprod(:,:)


!! Reading KE, PE and EE from fchk. Printed to the fchk file using the following programs: !!
!! From Gaussian09: /users/mgimferrer/FORTRAN/APOST3D.3.1-devel/utils/get_energy xxx.log >> xxx.fchk !!
!! From Gaussian16: /users/mgimferrer/FORTRAN/APOST3D.3.1-devel/utils/get_energy_g16 xxx.log >> xxx.fchk !!
c also included in FChk with orca2fchk program
      ekin0=ZERO
      eelnuc0=ZERO
      ecoul0=ZERO
      eexch0=ZERO
      ecorr0=ZERO
      evee0=ZERO
      eone0=ZERO

      rewind(15)
991   read(15,'(a80)',end=891) line
      if(index(line,"Kinetic Energy").ne.0) then
       read(line(45:71),'(es27.15)')  ekin0
      else
       go to 991
      end if
891   continue
      rewind(15)
992   read(15,'(a80)',end=892) line
      if(index(line,"Electron-Nuclei Energy").ne.0) then
       read(line(45:71),'(es27.15)')  eelnuc0
      else
       go to 992
      end if
892   continue
      rewind(15)
993   read(15,'(a80)',end=893) line
      if(index(line,"Electron-Electron Energy").ne.0) then
       read(line(45:71),'(es27.15)')  evee0
      else
       go to 993
      end if

893   continue
c ECP?
      iecp=0
      rewind(15)
998   read(15,'(a80)',end=897) line
      if(index(line,"ECP Mat").ne.0) then
       iecp=1                                 
       write(*,*) ' Pseudopotential matrix found in fchk file'
      else
       go to 998
      end if
897   continue
c Reading ECP matrix
      if(iecp.eq.1) then
999   read(15,'(a80)',end=1000) line
      if(index(line,"ECP-Mat").ne.0) then !! MG: Bug corrected. The blank was giving problems !!
       read(line(51:66),'(i17)') nn
       ALLOCATE(xecpv(nn))
       read(15,*) (xecpv(i),i=1,nn)
      else
       go to 999
      end if
c Transform to matrix
       ALLOCATE(xecpm(igr,igr))
       ALLOCATE(xprod(igr,igr))
       k=1
       do i=1,igr
        do j=1,i
         xecpm(i,j)=xecpv(k)   
         xecpm(j,i)=xecpm(i,j) 
         k=k+1
        end do 
       end do 
       DEALLOCATE(xecpv)

c ECP energy values for each atom Mulliken-type
       xecp=0.0d0
       do iiat=1,nat
        xxx=0.0d0
        do i=llim(iiat),iulim(iiat)
         do k=1,igr
          xxx=xxx+p(i,k)*xecpm(k,i)
         end do
        end do
        eecp(iiat)=xxx
        xecp=xecp+xxx
       end do
       write(*,'(a27,es27.15)') ' Pseudopotential Energy :',xecp
       write(*,'(5ES16.5)') (eecp(i),i=1,nat)
      end if

      return
1000  stop ' ECP Matrix not found'
      end 

c ****
      subroutine field_misc(ifield)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common/efield/field(4),edipole
      character*80 line


      ifield=0
C External electric field?
C REading and taking into acocunt ony x.y and z dipole components
      rewind(15)
1     read(15,'(a80)',end=2) line
      if(index(line,"External E-field").ne.0) then
       read(15,*) (field(i),i=1,4)
       xx=abs(field(2))+abs(field(3))+abs(field(4))
       if(xx.gt.1.0d-4) ifield=1
      else
       go to 1
      end if
2     continue
      end


c*******************************
c PROCESSING INP FILE   ********
c*******************************
      subroutine readchar(section,keyword,ival)
      character section*(*), keyword*(*)
      character linea*80
      integer ipos,ii

        ival=0
      call locate(16,section,ii)
      ii=0
      do while(ii.eq.0)
       read(16,"(a80)") linea
         ipos=index(linea,keyword)
         if(ipos.ne.0) then
        ival=1
        ii=1
         end if
       if(index(linea,"#").ne.0) ii=2 
      end do
c        if(ii.eq.1) write(*,"(a8,a2,i3)") keyword,'= ',ival
      return
        end
C*****************************************************************
      subroutine readreal(section,keyword,intv,intdef,ilog)
      character section*(*), keyword*(*)
      character linea*80
      real*8 intv,intdef
      integer ilog

      call locate(16,section,ii)
       if(ii.eq.0) go to 10
      ii=0
      do while(ii.eq.0)
      read(16,"(a80)") linea
         ipos=index(linea,keyword)
         if(ipos.ne.0) then
          ipos=ipos+1+len(keyword)
          read(linea(ipos:),*,err=10,end=10) intv
      ii=1
         end if
      if(index(linea,"#").ne.0) ii=2 
      end do
      if(ilog.eq.0.and.ii.ne.1) then
      write(*,*) keyword,'is required.'
      stop
      else if(ii.eq.2) then
10       intv=intdef
      end if
c        write(*,"(a8,a2,2e10.4)") keyword,'= ',intv,intdef
      return
        end
C*****************************************************************
        subroutine readint(section,keyword,intv,intdef,ilog)
        character section*(*), keyword*(*)
        character linea*80
        integer ipos,ii,intv,intdef,ilog

        call locate(16,section,ii)
        if(ii.eq.0) goto 40
        ii=0
        do while(ii.eq.0)
         read(16,"(a80)") linea
         ipos=index(linea,keyword)
         if(ipos.ne.0) then
          ipos=ipos+1+len(keyword)
          read(linea(ipos:),*,end=40) intv
          ii=1
         end if
         if(index(linea,"#").ne.0) ii=2
        end do
        if(ilog.eq.0.and.ii.ne.1) then
         write(*,*) keyword,"is required."
         stop
        else if(ii.eq.2) then
40       intv=intdef
        end if
        write(*,"(a8,a2,i3)") keyword,'= ',intv
        return
        end

C*****************************************************************
      subroutine locate(iunit,string,ii)
      integer iunit
      character string*(*)
      character*80 linia

      rewind(iunit)

      ii=0
      do while(ii.eq.0)
       read(iunit,"(a80)",end=10)linia
       if(index(linia,string).ne.0) then
          ii=1
          return
         end if
      end do
10      write(*,*) string, 'section not found '
      return
      end

      function int_locate(iunit,text,ilog)
      implicit double precision(a-h,o-z)
      character*(*)text
      character*80 line
      integer int_locate
      logical ilog

      ilog=.false.
      rewind(iunit)
1     read(iunit,'(a80)',end=99) line
      if(index(line,text).eq.0) then
       goto 1
      else
       read(line(50:61),*,err=99) int_locate
       ilog=.true.
       go to 2
      end if
99    int_locate=0
2     continue
      end

      function real_locate(iunit,text,ilog)
      implicit double precision(a-h,o-z)
      character*(*)text
      character*80 line
      real*8  real_locate
      logical ilog

      ilog=.false.
      rewind(iunit)
1     read(iunit,'(a80)',end=99) line
      if(index(line,text).eq.0) then
       goto 1
      else
       read(line(50:71),*,err=99) real_locate
       ilog=.true.
       go to 2
      end if
99    real_locate=0
2     continue
      end


      subroutine dm2input(idmrg,dm1,dm2)
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      dimension dm2(norb,norb,norb,norb)
      dimension dm1(nspinorb,nspinorb)
      character*80 line
      logical ilog

      iorca    = Iopt(43)
      ipyscf   = Iopt(87)
      ispinsep = Iopt(88)

      thres   = 1.0d-8

!! PREPARATING FOR DM1 AND DM2 READING !!

      do ii=1,norb
        do jj=1,norb
          do kk=1,norb
            do ll=1,norb
              dm2(ii,jj,kk,ll)=ZERO
            end do 
          end do 
        end do 
      end do


!! READING THE DM2(1,2,1,2) MATRIX BUT SAVING AS DM2(1,1,2,2), IMPORTANT FOR LATER ON !!

c reading full rdm2 spin-separated as given by  DMN code
c reading only i<=j, k<=l, i<=k elements and reconstructing full rdm2
c reading from formatted dm2 file

       if(ipyscf.ne.1.and.iorca.ne.1) then
          rewind(12)
          do while(.true.)
            read(12,end=999) i0,k0,j0,l0,vvv
            if(abs(vvv).gt.thres) then
              ii=(i0-1)/2+1
              kk=(k0-1)/2+1
              jj=(j0-1)/2+1
              ll=(l0-1)/2+1
              ilog=.true.
              if((i0.eq.j0.and.k0.eq.l0).or.(k0.eq.j0.and.i0.eq.l0)) ilog=.false.

!! IJ0 FOR THE AAAA, ABAB,BABA AND BBBB TERMS !!
!! KJ0 FOR THE ABBA AND BAAB TERMS !!

              ij0=mod(i0+j0,2)+mod(k0+l0,2)
              kj0=mod(k0+j0,2)+mod(i0+l0,2)
              if(ij0.eq.0)then
                dm2(ii,jj,kk,ll)=dm2(ii,jj,kk,ll)+vvv
                if(ilog) dm2(jj,ii,ll,kk)=dm2(jj,ii,ll,kk)+vvv
                if(k0.ne.i0.and.l0.ne.j0) then
                  dm2(kk,ll,ii,jj)=dm2(kk,ll,ii,jj)+vvv
                  if(ilog) dm2(ll,kk,jj,ii)=dm2(ll,kk,jj,ii)+vvv
                end if
              end if
              if(kj0.eq.0) then
                if(k0.ne.i0) then
                  dm2(kk,jj,ii,ll)=dm2(kk,jj,ii,ll)-vvv
                  if(ilog) dm2(jj,kk,ll,ii)=dm2(jj,kk,ll,ii)-vvv
                end if
                if(l0.ne.j0) then
                  dm2(ii,ll,kk,jj)=dm2(ii,ll,kk,jj)-vvv
                  if(ilog) dm2(ll,ii,jj,kk)=dm2(ll,ii,jj,kk)-vvv
                end if
              end if
            end if
          end do

        else

c reading active-space rdm2 spin-separated as given by PySCF/ORCA

        ninact=(nspinorb-ncasorb*2)
        ncore=ninact/2
c        iorb=ncore+ncasorb

c reconstructing core and core-active blocks

!! ALPHA-ALPHA AND BETA-BETA CASES !!

        do ispin=1,2
          do ii=1,ncore
            do jj=1,ncore
              dm2(ii,ii,jj,jj)=dm2(ii,ii,jj,jj)+ONE
              dm2(ii,jj,jj,ii)=dm2(ii,jj,jj,ii)-ONE
            end do 
            do kk=ncore+1,norb
              do ll=ncore+1,norb
               dm2(ii,ii,kk,ll)=dm2(ii,ii,kk,ll)+dm1(2*(kk-1)+ispin,2*(ll-1)+ispin)
               dm2(kk,ll,ii,ii)=dm2(kk,ll,ii,ii)+dm1(2*(kk-1)+ispin,2*(ll-1)+ispin)
               dm2(ii,ll,kk,ii)=dm2(ii,ll,kk,ii)-dm1(2*(kk-1)+ispin,2*(ll-1)+ispin)
               dm2(kk,ii,ii,ll)=dm2(kk,ii,ii,ll)-dm1(2*(kk-1)+ispin,2*(ll-1)+ispin)
              end do 
            end do 
          end do 
        end do 

!! ALPHA-BETA CASE !!

        do ii=1,ncore
          dm2(ii,ii,ii,ii)=dm2(ii,ii,ii,ii)+TWO
          do jj=1,ncore
            if(jj.ne.ii) dm2(ii,ii,jj,jj)=dm2(ii,ii,jj,jj)+TWO
          end do
          do kk=ncore+1,norb
            do ll=ncore+1,norb
              dm2(ii,ii,kk,ll)=dm2(ii,ii,kk,ll)+dm1(2*(kk-1)+1,2*(ll-1)+1)
              dm2(kk,ll,ii,ii)=dm2(kk,ll,ii,ii)+dm1(2*(kk-1)+2,2*(ll-1)+2)
              dm2(ii,ii,ll,kk)=dm2(ii,ii,ll,kk)+dm1(2*(ll-1)+2,2*(kk-1)+2)
              dm2(ll,kk,ii,ii)=dm2(ll,kk,ii,ii)+dm1(2*(ll-1)+1,2*(kk-1)+1)
            end do
          end do
        end do

!! NOW ACTIVE-ACTIVE !!

c reading only i<=j, k<=l, i<=k elements and reconstructing full rdm2
c reading from unformatted dm2 file
c correcting for shift in the index of active orbs

        rewind(12)
        read(12,*) line
        do while(.true.)
         read(12,*,end=999) i0,k0,j0,l0,vvv
         if(abs(vvv).gt.thres) then
          if(idmrg.eq.0) then
            ii=(i0-1)/2+1+ncore
            jj=(j0-1)/2+1+ncore
            kk=(k0-1)/2+1+ncore
            ll=(l0-1)/2+1+ncore
            ilog=.true.
            if((i0.eq.j0.and.k0.eq.l0).or.(k0.eq.j0.and.i0.eq.l0)) ilog=.false.
            ij0=mod(i0+j0,2)+mod(k0+l0,2)
            kj0=mod(k0+j0,2)+mod(i0+l0,2)
            if(ij0.eq.0)then
              dm2(ii,jj,kk,ll)=dm2(ii,jj,kk,ll)+vvv
              if(ilog) dm2(jj,ii,ll,kk)=dm2(jj,ii,ll,kk)+vvv
              if(k0.ne.i0.and.l0.ne.j0) then
                dm2(kk,ll,ii,jj)=dm2(kk,ll,ii,jj)+vvv
                if(ilog) dm2(ll,kk,jj,ii)=dm2(ll,kk,jj,ii)+vvv
              end if
            end if
            if(kj0.eq.0) then
              if(k0.ne.i0) then
                dm2(kk,jj,ii,ll)=dm2(kk,jj,ii,ll)-vvv
                if(ilog) dm2(jj,kk,ll,ii)=dm2(jj,kk,ll,ii)-vvv
              end if
              if(l0.ne.j0) then
                dm2(ii,ll,kk,jj)=dm2(ii,ll,kk,jj)-vvv
                if(ilog) dm2(ll,ii,jj,kk)=dm2(ll,ii,jj,kk)-vvv
              end if
            end if
           else ! reading spinless rdm2 all elements
             ii=i0+ncore
             jj=j0+ncore
             kk=k0+ncore
             ll=l0+ncore
             dm2(ii,jj,kk,ll)=vvv
           end if
         end if
       end do
999   end if

!! CHECKING DM2 !!

      xx2=ZERO
      do ii=1,norb
        do jj=1,norb
          xx2=xx2+dm2(ii,ii,jj,jj)
        end do
      end do
      write(*,*) " TRACE OF THE DM2 (NORMALIZED TO N(N-1) : ",xx2
      write(*,*) " "
      call flush
      end

      subroutine do_potential(zn,vv)
      use basis_set
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'parameter.h'
      PARAMETER(TOL=1.0d-8)
      dimension AminusB(3)
      dimension zn(maxat)
      real*8, allocatable ::sp(:,:)
      real*8 vv(nbasis,nbasis)
     
      DIMENSION E(0:10,0:10,0:10)
      DIMENSION E_coulomb(3,0:10)
      DIMENSION max_ang(3),xnew_coord(3)
      DIMENSION Boys(0:20)
      DIMENSION xR(0:20,0:10,0:10,0:10)
     
      allocate(sp(numprim,numprim))
      do ia=1,numprim
       do ib=1,ia
       sp(ia,ib)=0.0d0
       do ixyz=1,3 ! screening of primitives
        AminusB(ixyz)=coord(ixyz,iptoat(ia))-coord(ixyz,iptoat(ib))
!        ii=mod(nlm(ia,ixyz)+nlm(ib,ixyz),2)
!        if(abs(AminusB(ixyz)).lt.TOL.and.ii.ne.0) go to 111
       end do
       gamma_p=expp(ia)+expp(ib)
       eta_p=(expp(ia)*expp(ib))/gamma_p

       do idir=1,3
         if (ABS(gauss_dist).lt.1.0d-15) gauss_dist=0.0d0
         xnew_coord(idir)=(expp(ia)*coord(idir,iptoat(ia))+expp(ib)*coord(idir,iptoat(ib)))/gamma_p
         xdist_a=xnew_coord(idir)-coord(idir,iptoat(ia))
         xdist_b=xnew_coord(idir)-coord(idir,iptoat(ib))
         imax=nlm(ia,idir) 
         jmax=nlm(ib,idir) 
         max_ang(idir)=imax+jmax

         E(0,0,0)=exp(-eta_p*gauss_dist**2.0d0)
         if(imax.eq.0.and.jmax.eq.0) E_coulomb(idir,0)=E(0,0,0)
         do ii=0,imax-1
           do it=0,ii+1
             if ((it-1).lt.0.or.(it-1).gt.ii) then
               E1=0.0d0
             else
               E1=E(ii,0,it-1)/(2.0d0*gamma_p)
             end if
             if (it.gt.ii) then
               E2=0.0d0
             else
               E2=E(ii,0,it)*xdist_a
             end if
             if ((it+1).gt.ii) then
               E3=0.0d0
             else
               E3=real((it+1))*E(ii,0,it+1)
             end if
             E(ii+1,0,it)=E1+E2+E3
             if(ii.eq.imax-1.and.jmax.eq.0) E_coulomb(idir,it)=E(ii+1,0,it)
           end do
         end do

         do jj=0,jmax-1
           do it=0,imax+jj+1
             if (((it-1).lt.0).or.((it-1).gt.(jj+imax))) then
               E1=0
             else
               E1=E(imax,jj,it-1)/(2.0d0*gamma_p)
             end if
             if (it.gt.(jj+imax)) then
               E2=0
             else
               E2=E(imax,jj,it)*xdist_b
             end if
             if ((it+1).gt.(jj+imax)) then
               E3=0
             else
               E3=real((it+1))*E(imax,jj,it+1)
             end if
             E(imax,jj+1,it)=E1+E2+E3
             if(jj+1.eq.jmax) E_coulomb(idir,it)=E(imax,jmax,it)
           end do
         end do
         E_ij=E(imax,jmax,0)
       end do !end of idir loop

!!! RECURRENCE FOR COULOMB BEGINS, USING THE SAVED E VALUES (ALL VALUES OF T FOR IMAX,JMAX) ACROSS ALL 3 DIRECTIONS !!
!!! CALCULATING BOYS FUNCTION OF ORDER n_max (lower n has to be obtained through recurrence) !!!

       n_max=max_ang(1)+max_ang(2)+max_ang(3)
      do icenter=1,natoms
       dist_x=xnew_coord(1)-coord(1,icenter) !assuming hydrogen is atom 1
       dist_y=xnew_coord(2)-coord(2,icenter) !assuming hydrogen is atom 1
       dist_z=xnew_coord(3)-coord(3,icenter) !assuming hydrogen is atom 1
       if(dist_x.lt.1.0d-12) dist_x=0.0d0
       if(dist_y.lt.1.0d-12) dist_y=0.0d0
       if(dist_z.lt.1.0d-12) dist_z=0.0d0
       Rmod_CN2=dist_x**2.0d0+dist_y**2.0d0+dist_z**2.0d0 
       if(Rmod_CN2.lt.1.0d-12) then
         Boys(n_max)=1.0d0/(2.0d0*real(n_max)+1.0d0)
         xT=0.0d0
       else
         xT=gamma_p*Rmod_CN2
         if(n_max.eq.0) then
           Boys(n_max)=DSQRT(PI/(4.0d0*xT))*erf(DSQRT(xT))
         else if(xT.gt.8.0d0.and.xT.lt.15.0d0) then !Taylor expansion
           CALL Boys_expansionn(n_max,xT,Boys(n_max))
         else if(xT.gt.15.0d0) then !asymptotic value
           xvalue=fact(2*n_max-1)/2.0d0**(real(n_max)+1.0d0)
           Boys(n_max)=xvalue*DSQRT(PI/(xT**(2.0d0*real(n_max)+1.0d0)))
         else !infinite series
           xa=real(n_max)+0.50d0
           gammaincomplete=0.0d0
           do i=0,30
             xnewsum=xT**(xa+real(i))/((xa+real(i))*FACT(i))
             if(MOD(i,2).ne.0) xnewsum=-xnewsum
             gammaincomplete=gammaincomplete+xnewsum
           end do
           Boys(n_max)=0.50d0*gammaincomplete*xT**(-xa)
         end if
       end if
       if(Boys(n_max).lt.1.0d-12) Boys(n_max)=0.0d0

!!!  RECURRENCE FOR BOYS FUNCTION OF ALL ORDERS !!!
       n=1
       do while(n.le.n_max)
         Boys(n_max-n)=(2.0d0*xT*Boys(n_max-n+1)+exp(-xT))/(2.0d0*(real(n_max-n+1))-1.0d0)
         n=n+1
       end do

!!! RATING R VALUES FOR t=0,u=0,v=0 FOR N=[0,n_max] !!!
       do n=0,n_max
         xR(n,0,0,0)=(-2.0d0*gamma_p)**n*Boys(n)
       end do

!!!  RECURRENCE LOOP, TO OBTAIN R VALUES FOR ALL VALUES OF n,t,u,v !!!
!!!  INCLUDES FINAL CALCULATIONS, USING ALL THE OBTAINED VALUES OF xR(0,t,u,v) AND E(imax,jmax,t/u/v) STORED AS E_coulomb(idir,it) !!!

       Vab=0.0d0
       do it=0,max_ang(1)
         do iu=0,max_ang(2)
           do iv=0,max_ang(3)
             do n=0,n_max-it-iu-iv
               if(it.eq.0) then
                 xRa=0.0d0
               else
                 xRa=real(it)*xR(n+1,it-1,iu,iv)
               end if
               xR(n,it+1,iu,iv)=xRa+dist_x*xR(n+1,it,iu,iv)

               if(iu.eq.0) then
                 xRa=0.0d0
               else
                 xRa=real(iu)*xR(n+1,it,iu-1,iv)
               end if
               xR(n,it,iu+1,iv)=xRa+dist_y*xR(n+1,it,iu,iv)

               if(iv.eq.0) then
                 xRa=0.0d0
               else
                 xRa=real(iv)*xR(n+1,it,iu,iv-1)
               end if
               xR(n,it,iu,iv+1)=xRa+dist_z*xR(n+1,it,iu,iv)
             end do
             Vab=Vab+E_coulomb(1,it)*E_coulomb(2,iu)*E_coulomb(3,iv)*xR(0,it,iu,iv)
           end do
         end do
       end do
       sp(ia,ib)=sp(ia,ib)+zn(icenter)*Vab*2.0d0*PI/gamma_p
       if(ia.ne.ib) sp(ib,ia)=sp(ia,ib)
      end do 
      end do
      end do
! basis functions potential
      do i=1,nbasis
      do j=1,i
      vv(i,j)=0.0d0
      k=1
      do while(nprimbas(k,i).ne.0) 
      l=1
      do while(nprimbas(l,j).ne.0) 
       vv(i,j)=vv(i,j)+coefpb(nprimbas(k,i),i)*coefpb(nprimbas(l,j),j)*sp(nprimbas(k,i),nprimbas(l,j)) 
       l=l+1
      end do
       k=k+1
      end do
      if(i.ne.j) vv(j,i)=vv(i,j)
      end do
      end do
      deallocate(sp)
      end subroutine do_potential

      SUBROUTINE Boys_expansionn(n_max,xT,Boys)
      use basis_set
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION WA_Boys(11,29) 
      DATA (WA_Boys(n,1),n=1,11)/0.0195608d0,0.00364669d0,0.00111862d0,0.000468431d0,0.000242526d0,
     + 0.00014577d0,0.000097472d0,0.0000704136d0,0.0000538481d0,0.0000429782d0,0.0000354424d0/,
     +(WA_Boys(n,2),n=1,11)/0.0186829d0,0.00338106d0,0.00100873d0,0.000412112d0,0.000208955d0,
     + 0.000123469d0,0.000081445d0,0.0000582071d0,0.0000441371d0,0.0000349906d0,0.0000286997d0/,
     + (WA_Boys(n,3),n=1,11)/0.0178681d0,0.00314123d0,0.000911924d0,0.000363529d0,0.000180488d0,
     + 0.000104818d0,0.0000681859d0,0.0000481953d0,0.0000362266d0,0.0000285198d0,0.0000232616d0/,
     + (WA_Boys(n,4),n=1,11)/0.0171104d0,0.00292416d0,0.000826419d0,0.000321513d0,0.000156294d0, 
     + 0.0000891873d0,0.0000571985d0,0.0000399724d0,0.0000297754d0,0.0000232726d0,0.0000188722d0/,
     + (WA_Boys(n,5),n=1,11)/0.0164044d0,0.00272721d0,0.000750702d0,0.000285084d0,0.000135686d0,
     + 0.0000760629d0,0.0000480782d0,0.0000332091d0,0.000024508d0,0.0000190135d0,0.0000153263d0/,
     + (WA_Boys(n,6),n=1,11)/0.0157453d0,0.0025481d0,0.000683481d0,0.000253419d0,0.00011809d0,
     + 0.0000650203d0,0.0000404947d0,0.0000276383d0,0.0000202022d0,0.0000155529d0,0.0000124595d0/,
     + (WA_Boys(n,7),n=1,11)/0.015129d0,0.00238485d0,0.000623653d0,0.000225827d0,0.000103031d0,
     + 0.0000557101d0,0.0000341779d0,0.000023043d0,0.0000166778d0,0.0000127383d0,0.0000101396d0/,
     + (WA_Boys(n,8),n=1,11)/0.0145517d0,0.00223574d0,0.000570276d0,0.000201725d0,0.0000901145d0,
     + 0.0000478443d0,0.0000289067d0,0.0000192465d0,0.0000137895d0,0.0000104465d0,8.260577971693436E-6/,
     + (WA_Boys(n,9),n=1,11)/0.0140101d0,0.00209924d0,0.000522541d0,0.000180619d0,0.0000790087d0,
     + 0.0000411848d0,0.0000245001d0,0.0000161051d0,0.0000114193d0,8.578378262164428E-6,6.737300687156685E-6/,
     + (WA_Boys(n,10),n=1,11)/0.0135012d0,0.00197405d0,0.000479752d0,0.000162093d0,0.000069438d0,
     + 0.0000355347d0,0.0000208094d0,0.0000135017d0,9.471736878316576E-6,7.053926821346582E-6,5.501217677944795E-6/,
     + (WA_Boys(n,11),n=1,11)/0.0130222d0,0.00185901d0,0.000441309d0,0.000145792d0,0.000061171d0,
     + 0.0000307307d0,0.0000177125d0,0.0000113405d0,7.869163916680592E-6,5.80846024129448E-6,4.497200748449379E-6/,
     + (WA_Boys(n,12),n=1,11)/0.0125709d0,0.00175308d0,0.000406696d0,0.000131415d0,0.0000540135d0,
     + 0.0000266374d0,0.0000151089d0,9.543599724113801E-6,6.5486412555044375E-6,4.789710490139303E-6,3.680861022155174E-6/,
     + (WA_Boys(n,13),n=1,11)/0.012145d0,0.00165538d0,0.000375463d0,0.000118706d0,0.0000478025d0,
     + 0.0000231421d0,0.0000129157d0,8.04699180325652E-6,5.458961812050689E-6,3.955389710843769E-6,3.016431051707373E-6/,
     + (WA_Boys(n,14),n=1,11)/0.0117426d0,0.0015651d0,0.000347222d0,0.000107447d0,0.0000424005d0,
     + 0.000020151d0,0.0000110647d0,6.798376276547872E-6,4.558448846544275E-6,3.2712546857899363E-6,2.475068922112048E-6/,
     + (WA_Boys(n,15),n=1,11)/0.0113619d0,0.00148155d0,0.000321635d0,0.0000974484d0,0.0000376915d0,
     + 0.0000175859d0,9.499443638459046E-6,5.75485047731543E-6,3.81314628329268E-6,2.7095515558230534E-6,
     + 2.033499524941453E-6/,
     + (WA_Boys(n,16),n=1,11)/0.0110013d0,0.00140409d0,0.000298406d0,0.0000885512d0,0.0000335775d0,0.0000153814d0,
     + 8.173154150633697E-6,4.8811909545661805E-6,3.1953583574648872E-6,2.247765275090749E-6,1.672925359561938E-6/,
     + (WA_Boys(n,17),n=1,11)/0.0106594d0,0.00133217d0,0.000277279d0,0.000080617d0,0.0000299754d0,0.0000134827d0,
     + 7.04711984415125E-6,4.148441054539189E-6,2.6824702322432485E-6,1.8676134191372295E-6,1.378152893686431E-6/,
     + (WA_Boys(n,18),n=1,11)/0.0103348d0,0.00126529d0,0.000258027d0,0.0000735268d0,0.0000268145d0,0.0000118439d0,
     + 6.089189123641323E-6,3.532764059693503E-6,2.255994760108594E-6,1.5542360428544605E-6,1.1368914084822314E-6/,
     + (WA_Boys(n,19),n=1,11)/0.0100264d0,0.00120301d0,0.000240454d0,0.000067178d0,0.000024035d0,0.0000104263d0,
     + 5.272628628800348E-6,3.014511050397064E-6,1.9008013873868572E-6,1.2955429275308632E-6,9.391899322427788E-7/,
     + (WA_Boys(n,20),n=1,11)/0.00973295d0,0.00114494d0,0.000224384d0,0.0000614818d0,0.0000215856d0,9.197635833852241E-6,
     + 4.5751743306442665E-6,2.577462531412298E-6,1.6044918676611233E-6,1.0816872579180755E-6,7.76984784609612E-7/,
     + (WA_Boys(n,21),n=1,11)/0.00945357d0,0.00109071d0,0.000209665d0,0.0000563613d0,0.0000194227d0,8.130379517397552E-6,
     + 3.9782540122764285E-6,2.208210799121746E-6,1.3568943914649463E-6,9.046409242635738E-7,6.437357693257796E-7/,
     + (WA_Boys(n,22),n=1,11)/0.0091873d0,0.00104001d0,0.000196161d0,0.0000517497d0,0.0000175089d0,7.201431182817885E-6,
     + 3.4663493986594863E-6,1.8956564025560807E-6,1.1496533030678184E-6,7.578515640049953E-7,5.341334540326893E-7/,
     + (WA_Boys(n,23),n=1,11)/0.0089333d0,0.000992538d0,0.000183753d0,0.0000475888d0,0.0000158121d0,6.391208339399082E-6,
     + 3.026472197252E-6,1.6305971804591078E-6,9.758960363489164E-7,6.359653927498272E-7,4.438634874578624E-7/,
     + (WA_Boys(n,24),n=1,11)/0.00869078d0,0.000948047d0,0.000172333d0,0.0000438278d0,0.0000143048d0,5.683105061880854E-6,
     + 2.6477331561605937E-6,1.405392484813605E-6,8.299624811562363E-7,5.346030229794071E-7,3.694167080892861E-7/,
     + (WA_Boys(n,25),n=1,11)/0.00845904d0,0.000906297d0,0.000161809d0,0.0000404225d0,0.0000129633d0,5.063013022410489E-6,
     + 2.320987163294029E-6,1.2136885260823884E-6,7.071848651534656E-7,4.5017798995758116E-7,3.0793603821185376E-7/,
     + (WA_Boys(n,26),n=1,11)/0.00823742d0,0.000867075d0,0.000152096d0,0.0000373341d0,0.000011767d0,4.5189232096371906E-6,
     + 2.0385405792175185E-6,1.0501934551115288E-6,6.037085445372549E-7,3.7974972381135496E-7,2.5709294675278E-7/,
     + (WA_Boys(n,27),n=1,11)/0.00802531d0,0.000830187d0,0.000143118d0,0.0000345284d0,0.0000106983d0,4.0405942296332715E-6,
     + 1.7939095628505054E-6,9.104929579375757E-7,5.163459524955484E-7,3.2090432533957015E-7,2.1498769550217542E-7/,
     + (WA_Boys(n,28),n=1,11)/0.00782215d0,0.00079546d0,0.00013481d0,0.0000319756d0,9.741947025074953E-6,
     + 3.6192756244500433E-6,1.5816202292644905E-6,7.908988842176024E-7,4.424574466830895E-7,2.716578010993421E-7,
     + 1.8006872774705342E-7/,
     + (WA_Boys(n,29),n=1,11)/0.00762742d0,0.000762731d0,0.000127112d0,0.0000296492d0,8.884563981972003E-6,
     + 3.24747671603967E-6,1.3970431662671294E-6,6.883248391168378E-7,3.798539981494806E-7,2.303774548112768E-7,
     + 1.510674743511661E-7/
 
      xcheck_min=0.25d0
        do iT=1,29 !find which pretabulated xT is closest to the xT of Boys being calculated
         xpoint=(8.d0+0.25d0*real(iT-1))
         xcheck=ABS(xpoint-xT)
         if(xcheck.lt.xcheck_min) then
           xcheck_min=xcheck
           iT_min=iT
         end if
       end do

       Boys=WA_Boys(n_max,iT_min)
       xdeltaT=xT-(8.d0+0.25d0*real(iT_min-1))
       do i=1,5 !Taylor expand around the found closest value to obtain the Boys function of interest
         Taylor_term=WA_Boys(n_max+i,iT_min)*xdeltaT**real(i)/FACT(i)
         if(MOD(i,2).ne.0) Taylor_term=-Taylor_term
         Boys=Boys+Taylor_term
       end do

      END SUBROUTINE
