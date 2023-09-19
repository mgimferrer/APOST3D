       program testing
       implicit double precision (a-h,o-z)
       common/units/iunit
       character*80 line,name,dummy,comment
       integer :: fcof
      
       double precision :: energy  ,norm 
       double precision, allocatable:: coord(:) 
       integer         , allocatable:: priat(:) 
       integer         , allocatable:: prity(:) 
       double precision, allocatable:: expo(:) 
       double precision, allocatable:: coef(:) 
       double precision, allocatable:: occup(:) 
       double precision, allocatable:: coefmo(:,:) 
       double precision, allocatable:: orben(:) 
       double precision, allocatable:: znuc(:) 
       integer, allocatable:: iznuc(:) 
       integer, allocatable:: ishtyp(:) 
       integer, allocatable:: jshtyp(:) 
       integer, allocatable:: nprimsh(:) 
       integer, allocatable:: shelltoat(:) 
       double precision, allocatable:: exppr(:) 
       double precision, allocatable:: coefpr(:) 
       double precision, allocatable:: coordsh(:) 
       double precision, allocatable:: c(:,:) 
       double precision, allocatable:: c0(:,:) 
       double precision, allocatable:: p(:) 
       integer, allocatable:: atlist(:) 
       
       dimension expo_byat(40,100),coef_byat(40,100),nprimsh_byat(40)
       dimension nprimsh_bysh(40,100)
       Dimension  nconsh_byat(40) 
       dimension ishtyp_byat(40,60)
       character*2  atsymb0,atsymb(40) 


!
! ido5d=0  assume cartesian basis
       ido5d=0


! READ WFN
        rewind(1)
        open(unit=1,file='WFN',status='OLD') 

        read(1,'(a80)') comment
! nmol, nprim, natoms
        read(1,'(a80)') line
        read(line(18:23),*) nmol 
        read(line(38:43),*) npri 
        read(line(58:63),*) nat  

        allocate( coord(3*nat), priat(npri), prity(npri),znuc(nat))
        allocate ( expo(npri) ,coefmo(nmol,npri), orben(nmol),occup(nmol) )  
        allocate(iznuc(nat),jshtyp(npri),coef(npri))

       
       ipos=0
       do i=1,nat
        read(1,'(a24,3f12.8,a10,f5.1)')dummy ,(coord(ipos+ii),ii=1,3),dummy,znuc(i)
        ipos=ipos+3
       end do
       icharge=0
       do i=1,nat
        iznuc(i)=int(znuc(i))
        icharge=icharge+iznuc(i)
       end do
       
       ii=npri/20
       ipos=0
       do i=1,ii  
        read(1,'(a20,20i3)')dummy ,(priat(ipos+ii),ii=1,20)
        ipos=ipos+20
       end do
       ii=mod(npri,20)
       if(ii.ne.0) read(1,'(a20,20i3)')dummy ,(priat(ipos+i),i=1,ii)

       ii=npri/20
       ipos=0
       do i=1,ii  
        read(1,'(a20,20i3)')dummy ,(prity(ipos+ii),ii=1,20)
        ipos=ipos+20
       end do
       ii=mod(npri,20)
       if(ii.ne.0) read(1,'(a20,20i3)')dummy ,(prity(ipos+i),i=1,ii)

       ii=npri/5
       ipos=0
       do i=1,ii  
        read(1,'(a10,5e14.7)')dummy ,(expo(ipos+ii),ii=1,5)
        ipos=ipos+5
       end do
       ii=mod(npri,5)
       if(ii.ne.0) read(1,'(a10,5E14.7)')dummy ,(expo(ipos+i),i=1,ii)

       icorr=0 
       elec=0.0d0
       do i=1,nmol
        read(1,'(a80)')line                             
        read(line(63:74),*) orben(i) 
        read(line(35:47),*) occup(i) 
        elec=elec+occup(i)
        if(occup(i).ne.2.0d0.and.occup(i).ne.1.0d0) icorr=1
        read(1,*)(coefmo(i,j),j=1,npri)                             
       end do
       nelec =anint(elec)
       icharge=icharge-nelec


       read(1,'(a80)')line    
       if (index(line,"END DATA").eq.0) stop 'Problem reading WFN'
       read(1,'(a80)')line   
       read(line(18:37),*) energy

       close(1)


       WRITe(*,*) '     '
       WRITe(*,*) 'Total number of primitives: ',npri
       WRITe(*,*) 'Number of atoms: ',nat  
       WRITe(*,*) 'Number of molecular orbitals: ',nmol
       WRITe(*,*) 'Number of electrons (from occupations): ',nelec
       WRITe(*,*) 'Molecular charge (from occupations): ',icharge
       WRITe(*,*) 'Total Energy: ',energy 
       WRITe(*,*) '     '
       WRITe(*,*) '-----------------------'
       WRITe(*,*) 'DATA READ from WFN FILE'
       WRITe(*,*) '-----------------------'


! READ PNOF/NWCHEM OUTPUT
       inwchem=0
       call getarg (1,name)
       write(*,*) 'Output file: ',name
       name=trim(name)
       open(unit=1,file=name,status='OLD') 
       call getarg (2,name)
       if(index(name,"gbs").ne.0) then
        write(*,*) 'Reading gbs file'
        stop
       else if(index(name,"nwchem").ne.0) then
        inwchem=1
        write(*,*) 'Reading NWCHEM file'
       else
        write(*,*) 'Reading PNOF file'
       end if

       if (inwchem.eq.1) then

       nattyp=0
       rewind(1)
51    read(1,'(a80)')line                             
       if(index(line,"NWChem Property Module").eq.0) go to 51
52    read(1,'(a80)')line                             
       if(index(line,"Functions and Types").eq.0) go to 52
      read(1,'(a80)')line                             

53    read(1,'(a80)')line                             
      if(index(line(2:2)," ").eq.0) then  
       nattyp=nattyp+1
       read(line(51:54),*) nconsh_byat(nattyp) 
       go to 53
      end if

       rewind(1)
54    read(1,'(a80)')line                             
       if(index(line,"Wavefunction type:  closed shell").eq.0) go to 54
       icalctype=0

       rewind(1)
55     read(1,'(a80)')line                             
       if(index(line,'Total DFT energy').eq.0) go to 55
       read(line(28:49),*) energy

       rewind(1)
56       read(1,'(a80)')line                             
       if(index(line,"AO basis - number of functions:").eq.0) go to 56
       read(line(42:47),*) nbas
!       nmol=nbas

       rewind(1)
57     read(1,'(a80)')line                             
       if(index(line,"           number of shells:").eq.0) go to 57
       read(line(39:44),*) nconsh 

       rewind(1)
58     read(1,'(a80)')line                             
       if(index(line,"No. of atoms  ").eq.0) go to 58
       read(line(29:34),*) natoms 
       read(1,'(a28,i6)') dummy,nelec  
       read(1,'(a28,i6)') dummy,nal    
       read(1,'(a28,i6)')dummy, nb     
       read(1,'(a28,i6)') dummy,icharge
       read(1,'(a28,i6)') dummy,multip 
       nat=natoms

       rewind(1)
59     read(1,'(a80)')line                             
       if(index(line,'Basis "ao basis" -> "ao basis"').eq.0) go to 59
        read(1,'(a80)')line                             



       do iattyp=1,nattyp 
        nprimsh0=0
         if(iattyp.eq.1) read(1,'(a80)')line                             
         read(line(3:4),'(a2)')  atsymb(iattyp)
         read(1,'(a80)')line                             
         read(1,'(a80)')line                             
         read(1,'(a80)')line                             

         do ii=1,nconsh_byat(iattyp)+1
          nprish=0
60        read(1,'(a80)')line         ! first line after atomic symbol                    
          if(index(line,".").ne.0) then  !not empty line
           nprish=nprish+1
           nprimsh0=nprimsh0+1
           read(line(6:21),*) expo_byat(iattyp,nprimsh0)
           read(line(22:31),*) coef_byat(iattyp,nprimsh0)    
           read(line(5:5),*) dummy    
           go to 60
          else
           nprimsh_bysh(iattyp,ii)=nprish
           if(index(dummy,"S").ne.0) then
            ishtyp_byat(iattyp,ii)=0
           else if(index(dummy,"P").ne.0) then
            ishtyp_byat(iattyp,ii)=1
           else if(index(dummy,"D").ne.0) then
            ishtyp_byat(iattyp,ii)=2
            if(ido5d.eq.1) ishtyp_byat(iattyp,ii)=-2
           else if(index(dummy,"F").ne.0) then
            ishtyp_byat(iattyp,ii)=3
            if(ido5d.eq.1) ishtyp_byat(iattyp,ii)=-3
           else if(index(dummy,"G").ne.0) then
            ishtyp_byat(iattyp,ii)=0
           else
            stop 'unknown shell type'
           end if
          end if
         end do
         nprimsh_byat(iattyp)=nprimsh0
         write(*,'(a5,2x,a2)') 'ATOM ',atsymb(iattyp)
         write(*,*) 'nprimsh by at',nprimsh_byat(iattyp)
         write(*,*) 'nconsh_byat',nconsh_byat(iattyp)
         write(*,*) 'nprimsh by sh',(nprimsh_bysh(iattyp,j),j=1,nconsh_byat(iattyp))
         write(*,*) 'ishtyp_byat',(ishtyp_byat(iattyp,j),j=1,nconsh_byat(iattyp))
         write(*,*) 'coeff_byat',(coef_byat(iattyp,j),j=1,nprimsh0)
         write(*,*) 'expo_byat',(expo_byat(iattyp,j),j=1,nprimsh0)
         write(*,*) '  '
       end do


! atomic coordinates
!       allocate( coord(3*nat),znuc(nat))
!       allocate ( orben(nbas),occup(nbas) )  
!       allocate(iznuc(nat))
       allocate(atlist(natoms))

       rewind(1)
61     read(1,'(a80)')line                             
       if(index(line,'Output coordinates in angstroms').eq.0) go to 61
       read(1,'(a80)')line                             
       read(1,'(a80)')line                             
       read(1,'(a80)')line                             

       nprish=0
       ipos=0
       do i=1,natoms
        read(1,'(a80)')line                             
        read(line(7:8),'(a2)') atsymb0
        write(*,*) 'atom read ',atsymb0
        write(*,'(3a2)') atsymb0,atsymb(1),atsymb(2)
        read(line(24:33),*) znuc(i)
        read(line(34:48),*) coord(ipos+1)
        read(line(49:63),*) coord(ipos+2)
        read(line(64:78),*) coord(ipos+3)
        ipos=ipos+3
        icontrol=0
        do j=1,nattyp
         if(index(atsymb(j),atsymb0).ne.0) then
          nprish=nprish+ nprimsh_byat(j)
          atlist(i)=j
          icontrol=1
         end if
        end do
        if(icontrol.eq.0) stop 'Problem getting atoms'
       end do
       write(*,*) (atlist(i),i=1,natoms)
       do i=1,3*natoms
        coord(i)=coord(i)*1.889725989
       end do
       icharge=0
        do i=1,nat
         iznuc(i)=int(znuc(i))
         icharge=icharge+iznuc(i)
        end do
! RHF
       do i=1,nmol
        orben(i)=-99.0d0
        if(i.le.nal) then 
         occup(i)=2.0d0
        else
         occup(i)=0.0d0
        end if
       end do

       WRITe(*,*) '     '
       WRITe(*,*) 'Number of contracted shells: ',nconsh
       WRITe(*,*) 'Number of primitive shells: ',nprish
       WRITe(*,*) 'Number of basis functions : ',nbas 
       WRITe(*,*) 'Multiplicity : ',multip
       WRITe(*,*) 'Total Energy: ',energy 
       WRITe(*,*) '     '

       WRITe(*,*) '--------------------------'
       WRITe(*,*) 'DATA READ from NWCHEM FILE'
       WRITe(*,*) '--------------------------'

       else if(inwchem.eq.0) then

       rewind(1)
10     read(1,'(a80)')line                             
       if(index(line,"TOTAL NUMBER OF PRIMITIVE EXPON").eq.0) go to 10
       read(line(48:52),*) nprish

       rewind(1)
11      read(1,'(a80)')line                             
       if(index(line,"TOTAL NUMBER OF BASIS SET").eq.0) go to 11
       read(line(48:52),*) nconsh

       rewind(1)
15       read(1,'(a80)')line                             
       if(index(line,"STATE MULTIPLICITY").eq.0) go to 15
       read(line(48:52),*) multip

       rewind(1)
20       read(1,'(a80)')line                             
       if(index(line,"NUMBER OF CARTESIAN GAUSSIAN").eq.0) go to 20
       read(line(48:52),*) nbas

       rewind(1)
25       read(1,'(a80)')line                            
       if(index(line,"ATOMIC BASIS SET").eq.0) go to 25
       do i=1,8
        read(1,'(a80)')line                          
       end do

       WRITe(*,*) '     '
       WRITe(*,*) 'Number of contracted shells: ',nconsh
       WRITe(*,*) 'Number of primitive shells: ',nprish
       WRITe(*,*) 'Number of basis functions : ',nbas 
       WRITe(*,*) 'Multiplicity : ',multip
       WRITe(*,*) 'Total Energy: ',energy 
       write(*,*)"Number of CAS Electrons (assuming PNOF5): ",nelec    
       write(*,*)"Number of CAS Orbitals (assuming PNOF5): ",nelec    
       WRITe(*,*) '     '
       
       WRITe(*,*) '--------------------------'
       WRITe(*,*) 'DATA READ from PNOF FILE'
       WRITe(*,*) '--------------------------'

       end if


!
! EVERYTHING IS READ. NOW BUILDING ARRAYS AND MATRICES FOR FCHK OUTPUT
!

       allocate(ishtyp(nconsh),nprimsh(nconsh),shelltoat(nconsh), exppr(nprish), coefpr(nprish)) 
       allocate(coordsh(3*nconsh))
       allocate(c(nmol,nbas))
       allocate(c0(nmol,nbas))
!       allocate(c(nbas,nbas))
!       allocate(c0(nbas,nbas))
       allocate(p(nbas*(nbas+1)/2))


       if(inwchem.eq.1) then

        ii=0
        jj=0
        do iat=1,natoms
         do i=1,nconsh_byat(atlist(iat))
          ii=ii+1
          shelltoat(ii)=iat
          nprimsh(ii)=nprimsh_bysh(atlist(iat),i)
          ishtyp(ii)=ishtyp_byat(atlist(iat),i)
         end do
         do i=1,nprimsh_byat(atlist(iat))
          jj=jj+1
          coefpr(jj) = coef_byat(atlist(iat),i)
          exppr(jj) = expo_byat(atlist(iat),i)
         end do
        end do
        maxang=0
        maxcon=1
        do i=1,nconsh
         if(abs(ishtyp(i)).gt.maxang) maxang=abs(ishtyp(i))
         if(nprimsh(i).gt.maxcon) maxcon=nprimsh(i)
        end do

       else

       iprsh0=0
       iconsh0=1
       iat0=1
       maxang=0
       maxcon=0

30      read(1,'(a80)')line                          
       if(line(6:7).eq." ") then
        read(1,'(a80)')line                          
        if(index(line,"-----").ne.0) go to 33
        read(1,'(a80)')line                          
        shelltoat(ii)=iat0 
        nprimsh(ii)=jj-iprsh0
        if(nprimsh(ii).gt.maxcon) maxcon=nprimsh(ii) 
        iprsh0=jj
        iconsh0=ii+1
        iat0=iat0+1
        go to 30
  !     else if(line(6:7).ne."ri") then
       else if(jj.ne.nprish) then
        read(line(1:7),*) ii 
        read(line(8:11),*) dummy 
        if(index(dummy,"S").ne.0) then
         ishtyp(ii)=0
        else if(index(dummy,"P").ne.0) then
         ishtyp(ii)=1
        else if(index(dummy,"D").ne.0) then
         ishtyp(ii)=2
        else if(index(dummy,"F").ne.0) then
         ishtyp(ii)=3
        else if(index(dummy,"G").ne.0) then
         ishtyp(ii)=0
        else
         stop 'unknown shell type'
        end if
        if(ishtyp(ii).gt.maxang) maxang=ishtyp(ii)
         
        read(line(12:19),*) jj 
        read(line(42:59),*) coefpr(jj) 
        read(line(20:41),*) exppr(jj) 

        if(ii.ne.iconsh0) then
         shelltoat(ii-1)=iat0 
         nprimsh(ii-1)=jj-iprsh0-1
         if(nprimsh(ii-1).gt.maxcon) maxcon=nprimsh(ii-1) 
         iprsh0=jj-1
         iconsh0=ii
        end if
        go to 30 
       else
33      shelltoat(ii)=iat0 
        nprimsh(ii)=jj-iprsh0
        iprsh0=jj
       end if

       if(iat0.ne.nat) stop 'problem reading output file.1'
       if(iprsh0.ne.nprish) stop 'problem reading output file.2'
       if(ii.ne.nconsh) stop 'problem reading output file.3'

       end if

! primitive type for normalization
       ipos=0
       jpos=0
       do k=1,nconsh
        jj=abs(ishtyp(k))
        do j=1,fcof(jj,ido5d) 
         do i=1,nprimsh(k)
          ipos=ipos+1
          jpos=jpos+1
          if(j.ne.1.and.i.eq.1) jpos=jpos-nprimsh(k)
          coef(ipos)=coefpr(jpos)
          jshtyp(ipos)=0
          if(jj.eq.2) then
           if(j.le.3) jshtyp(ipos)=1
          else if(jj.eq.3) then
           if(j.le.3) then 
            jshtyp(ipos)=2
           else if(j.ne.10)then 
            jshtyp(ipos)=1
           end if
          else if(abs(ishtyp(k)).eq.4) then
           if(j.le.3) then 
            jshtyp(ipos)=3
           else if(j.le.9)then 
            jshtyp(ipos)=2
           else if(j.le.12)then 
            jshtyp(ipos)=4
           else  
            jshtyp(ipos)=1
           end if
          end if
         end do
        end do
       end do
           
! renomrlaization of primitive coeeficients
      jpos=0
      do k=1,nconsh
       jj=abs(ishtyp(k))
       do j=1,fcof(jj,ido5d) 
        do i=1,nprimsh(k)
         jpos=jpos+1
         coef(jpos)=coef(jpos)*norm(jj,jshtyp(jpos),expo(jpos))
      write(*,*) 'title',jpos,coef(jpos)
        end do
       end do
      end do

! renormalization of basis functions
      ipos=0
      jpos1=0
      jpos=0
      do k=1,nconsh
       jj=abs(ishtyp(k))
       do j=1,fcof(jj,ido5d) 
        xx=0.0d0
        do i1=1,nprimsh(k)
         jpos1=jpos1+1
         do i2=1,nprimsh(k)
          jpos2=jpos1-i1+i2
          xx=xx+coef(jpos1)*coef(jpos2)*(norm(jj,jshtyp(jpos1),(expo(jpos1)+expo(jpos2))/2.0d0))**(-2.0d0)
         end do
        end do
!
        do i=1,nprimsh(k)
         jpos=jpos+1
         if(i.eq.1) then
          ipos=ipos+1
          do iorb=1,nmol
           c(iorb,ipos)=coefmo(iorb,jpos)/coef(jpos)*dsqrt(xx)
          end do
         end if
        end do
!
       end do
      end do
     
! rearranging columns if F orbitals are present
       if(maxang.ge.3) then
       ipos=0
       do k=1,nconsh
        jj=abs(ishtyp(k))
        do j=1,fcof(jj,ido5d) 
         ipos=ipos+1
         kk=0
         if(jj.eq.3.and.j.eq.4) kk=1
         if(jj.eq.3.and.j.eq.5) kk=1
         if(jj.eq.3.and.j.eq.6) kk=3
         if(jj.eq.3.and.j.eq.7) kk=-3
         if(jj.eq.3.and.j.eq.8) kk=-1 
         if(jj.eq.3.and.j.eq.9) kk=-1 
         do iorb=1,nmol
          c0(iorb,ipos+kk)=c(iorb,ipos)
         end do
        end do
       end do
       c=c0
       end if


       ipos=1
       do i=1,nconsh
        iat=shelltoat(i)
        coordsh(ipos)=coord(3*(iat-1)+1)
        coordsh(ipos+1)=coord(3*(iat-1)+2)
        coordsh(ipos+2)=coord(3*(iat-1)+3)
        ipos=ipos+3
       end do

!
       ipos=0
        do i=1,nbas
         do j=1,i    
          ipos=ipos+1
          p(ipos)=0.0d0
          do k=1,nmol
           p(ipos)=p(ipos)+occup(k)*c(k,i)*c(k,j)
          end do
         end do
        end do


       close(1)

!       do i=1,nbas
!        write(*,'(i5,6f12.5)') i,(c(k,i),k=1,6)
!       end do
 
! WRITTING FCHK FILE

       iunit=1
       open(unit=iunit,file='FCHK',status='UNKNOWN') 

       write(1,'(a80)') comment
       if(icorr.eq.1) then
!        write(1,'(a,a)') "SP        ","CASSCF"
        write(1,'(a)') "SP        "
       else
        write(1,'(a,a)') "SP        ","RHF"
       end if
       call ival("Number of atoms",nat)
       call ival("Charge",icharge)
       call ival("Multiplicity",multip)
       call ival("Number of electrons",nelec)
       nal=(multip-1+nelec)/2
       nbe=nelec-nal
       call ival("Number of alpha electrons",nal)
       call ival("Number of beta electrons",nbe)
       call ival("Number of basis functions",nbas)
       call ival("Number of independant functions",nmol)
       call ival("Number of contracted shells",nconsh)
       call ival("Highest angular momentum",maxang)
       call ival("Largest degree of contraction",maxcon)
       call ival("Number of primitive shells",nprish)

       call rval("Viriral Ratio",0.0d0)
       call rval("SCF Energy",energy)
       call rval("Total Energy",energy)

       call iarr("Atomic numbers",nat,iznuc)
       call rarr("Nuclear charges",nat,znuc)
       call rarr("Current cartesian coordinates",3*nat,coord)
       call iarr("Shell types",nconsh,ishtyp)
       call iarr("Number of primitives per shell",nconsh,nprimsh)
       call iarr("Shell to atom map",nconsh,shelltoat)
       call rarr("Coordinates of each shell",3*nconsh,coordsh)
       call rarr("Primitive exponents",nprish,exppr)
       call rarr("Contraction coefficients",nprish,coefpr)
       call rarr("Atomic Orbital Energies",nmol,orben)
       call rmat("Alpha MO coefficients",nmol,nbas,c)
       call rarr("Total SCF Density",nbas*(nbas+1)/2,p)
! for PNOF5 always equiv to CASSCF(n,n)
       if(inwchem.eq.0) then
       call ival("Number of CAS Electrons",nelec)     
       call ival("Number of CAS Orbitals",nelec)    
! printing nat orbital occupancies
       call rarr("Natural orbital occupancies",nmol,occup)
       end if

       WRITe(*,*) '  '
       WRITe(*,*) '---------------------------'
       WRITe(*,*) 'FCHK file has been created'
       WRITe(*,*) '---------------------------'
       end
       


!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine ival(key,ivalue)
       implicit double precision (a-h,o-z)
       common/units/iunit
       character(len=*) ::  key 
       integer ivalue
       character(len=43) ::  title
       title=adjustl(key)
       write(iunit,'(A43,A,I17)') title,"I",ivalue
       end

       subroutine rval(key,rvalue)
       implicit double precision (a-h,o-z)
       common/units/iunit
       character(len=*) ::  key 
       double precision :: rvalue
       character(len=43) ::  title
       title=adjustl(key)
       write(iunit,'(A43,A,ES27.15)') title,"R",rvalue
       end

       subroutine iarr(key,ival,iarray)
       implicit double precision (a-h,o-z)
       character(len=*) ::  key 
       common/units/iunit
       integer ival
       integer, dimension(ival) :: iarray
       character(len=43) ::  title
       title=adjustl(key)
       write(iunit,'(A43,A6,I12)') title,"I   N=",ival
       write(iunit,'(6I12)') (iarray(i),i=1,ival)
       end

       subroutine rarr(key,ival,rarray)
       implicit double precision (a-h,o-z)
       character(len=*) ::  key 
       integer ival
       common/units/iunit
       double precision, dimension(ival) :: rarray
       character(len=43) ::  title
       title=adjustl(key)
       write(iunit,'(A43,A6,I12)') title,"R   N=",ival
       write(iunit,'(5ES16.8)') (rarray(i),i=1,ival)
       end

       integer function fcof(ival,ido5d)
       implicit double precision (a-h,o-z)
       integer ival,icof

       if (ido5d.eq.1) then
         fcof=2*(ival-1)+1
       else
        ii=ival+2
        icof=1
        do i=ii,ii-1,-1
         icof=icof*i
        end do
        fcof=icof/2
       end if
       end
        
       subroutine rmat(key,ival,jval,rmatrix)
       implicit double precision (a-h,o-z)
       character(len=*) ::  key 
       integer ival,jval
       common/units/iunit
       double precision, dimension(ival,jval) :: rmatrix   
       character(len=43) ::  title
       title=adjustl(key)
       write(iunit,'(A43,A6,I17)') title,"R   N=",ival*jval
       write(iunit,'(5ES16.8)') ((rmatrix(i,j),j=1,jval),i=1,ival)
       end


       double precision  function norm(nn,itype,alpha)
       implicit double precision (a-h,o-z)
       integer nn,itype
       double precision alpha

       xn=real(nn)
       xval=2.0d0**xn*alpha**((2.0d0*xn+3.0d0)/4.0d0)*(2.0d0/dacos(-1.0d0))**(3.0d0/4.0d0)
       norm=1.0d0
! 200 210 211
       if(itype.eq.1) norm=1.0d0/dsqrt(3.0d0)
! 300 310 311
       if(itype.eq.2) norm=1.0d0/dsqrt(15.0d0)
! 400 410 411
       if(itype.eq.3) norm=1.0d0/dsqrt(105.0d0)
! 220 221    
       if(itype.eq.4) norm=1.0d0/dsqrt(9.0d0)

       norm=xval*norm
       end


