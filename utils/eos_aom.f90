! This program performs basic population and EOS analysis from
! AOM generated either by Multiwfn or Aimall. In addition, the 
! either the corresponding wfn or fchk files must be available

! Use keyword MULTIWFN to read a single $mol.aom file
! Use keyword AIMALL to read a multiple *.int files
! Aimall's int files must have first capital letter for atomname
! The following can be useful:
! for i in *.int; do mv "$i" "`echo $i | sed "s/^./\U&/"`"; done
!
! Use -wfn as second argument to indicate reading data from wfn
! instead of from fchk file (default)
!
! Not tested beyond RHF/RKS/UHF/UKS
! Not tested with effective core potentials
!

       program eos
       implicit double precision (a-h,o-z)
       parameter(maxat=500,maxfrag=50)
       character*80 line,name0,name1,name
       character*2 symb(maxat)
! fragments
       dimension ifrlist(maxat,maxfrag),nfrlist(maxfrag),jfrlist(maxat)
       dimension navect(maxat),qdev(maxat),ifchar(maxfrag)
! atoms
       dimension zn(maxat),iznuc(maxat)
!atomic symbols
      character*2 mend(92)
      data mend/"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg",&
      "Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn",&
      "Fe", "Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb",&
      "Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",&
      "Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sn",&
      "Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re",&
      "Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr",&
      "Ra","Ac","Th","Pa","U" /


! eos
       allocatable ip0(:),oxi(:),elec(:),occupg(:),occup2(:)
       allocatable p0net(:,:),occup(:,:),iorbat(:,:)
! pop
       allocatable ss(:,:,:),sa(:,:,:),sb(:,:,:),bo(:,:)
       allocatable q(:),qa(:),qb(:),qs(:)
       allocatable qfrag(:),qsfrag(:),bofrag(:,:)
       allocatable scr(:,:),c(:,:),occ(:)

       write(*,*)'This program performs basic population and EOS analysis from'
       write(*,*)'AOM generated either by Multiwfn or Aimall. In addition, the' 
       write(*,*)'either the corresponding wfn or fchk files must be available'
       write(*,*)''
       write(*,*)'Use keyword MULTIWFN to read a single $mol.aom file'
       write(*,*)'Use keyword AIMALL to read a multiple *.int files'
       write(*,*)"Aimall's int files must have first capital letter for "
       write(*,*)'atomname. The following can be useful:'
       write(*,*)'for i in *.int; do mv "$i" "`echo $i | sed "s/^./\U&/"`"; done'
       write(*,*)''
       write(*,*)'Use -wfn as second argument to indicate reading data from wfn'
       write(*,*)'instead of from fchk file (default)'
       write(*,*)''
       write(*,*)'Not tested beyond RHF/RKS/UHF/UKS'
       write(*,*)'Not tested with effective core potentials'
         
       call getarg(1,name0)
       call getarg(2,name1)

! default reading fchk file
       ifchk=1
       iwfn=0
       if(index(name1,"-wfn").ne.0) then
         iwfn=1
         ifchk=0
       end if
       iuhf=0
       icas=0

! read nmos, nprim, natoms from wfn
! unfortunately, no info about RHF vs UHF is provided
       if(iwfn.eq.1) then
        name=trim(name0)//".wfn"
        open(1, file=name,status="OLD")
        read(1,'(a80)')  line
        read(1,'(a80)')  line
        read(line(19:23),*) nmo
        read(line(55:63),*) nat  
        write(*,*) ' NUMBER OF MOs:   ',nmo
        write(*,*) ' NUMBER OF ATOMS: ',nat 
        do iat=1,nat
         read(1,'(a80)')  line
         k=50 
         do while(index(line(k:k),"=").eq.0)
            k=k+1
         end do
         read(line(k+1:),*) zn(iat)  
         read(line(3:4),*) symb(iat)  
          iznuc(iat)=int(zn(iat))
        end do
        close(1)

! read info from fchk
       else if(ifchk.eq.1) then

        name=trim(name0)//".fchk"
        open(1, file=name,status="OLD")
! number of atoms
        iflag=1
        do while (iflag.eq.1)
         read(1,'(a80)')  line
         if (index(line,"Number of atoms").ne.0) iflag=0
        end do
        read(line(45:),*) nat    
! read number electrons
        iflag=1
        do while (iflag.eq.1)
         read(1,'(a80)')  line
         if (index(line,"Number of alpha electrons").ne.0) iflag=0
        end do
        read(line(45:),*) nalpha 
        nmo=nalpha
        read(1,'(a80)')  line
        read(line(45:),*) nbeta  
        if(nalpha.ne.nbeta) then
         iuhf=1  ! warning the possibility of a RO WF
         nmo=nalpha+nbeta
        end if
! atomic numbers
        iflag=1
        do while (iflag.eq.1)
         read(1,'(a80)')  line
         if (index(line,"Nuclear charges").ne.0) iflag=0
        end do
        read(1,*) (zn(i),i=1,nat)
        do iat=1,nat
          iznuc(iat)=int(zn(iat))
          symb(iat)=mend(iznuc(iat))
        end do
        close(1)
       end if
       write(*,*) ' NUMBER OF ALPHA ELECTRONS:',nalpha
       write(*,*) ' NUMBER OF BETA ELECTRONS:',nbeta
       write(*,*) ' NUMBER OF ATOMS: ',nat 


! read  apost input file
       name=trim(name0)//".inp"
       open(1, file=name,status="OLD")
       call locate(1,"# METODE",ii)
       imultiwfn=0
       iaimall=0
       idofr=0
       idoat=0
       icube=0
       ieffthr=10
       thres=2.5d-3
       iflag=1
       do while (iflag.eq.1)
        read(1,'(a80)')  line
        if (index(line,"#").ne.0) iflag=0
        if (index(line,"AIMALL").ne.0) iaimall=1 
        if (index(line,"MULTIWFN").ne.0) imultiwfn=1   
!
        if (index(line,"DOFRAGS").ne.0) idofr=1
        if (index(line,"DOATOMS").ne.0) idoat=1
        if (index(line,"CUBE").ne.0) icube=1
        if (index(line,"EOS_THRESH").ne.0) then
         ii=index(line,"=") 
         read(line(ii+1:),*) thres 
        end if
        if (index(line,"EFF_THRESH").ne.0) then 
         ii=index(line,"=")
         read(line(ii+1:),*) ieffthr
        end if
       end do
!
       if(iaimall+imultiwfn.ne.1)stop 'indicate Multiwfn or Aimall data'
!
       if(icube.eq.1) then
        call locate(16,"# CUBE",ii)
        jcubthr=1000
        kcubthr=0
        if(ii.eq.1) then 
         read(1,'(a80)')  line
         ii=index(line,"MAX_OCC=")
         if(ii.ne.0) read(line(ii+1:),*) jcubthr
         read(1,'(a80)')  line
         ii=index(line,"MIN_OCC=")
         if(ii.ne.0) read(line(ii+1:),*) kcubthr
        end if
       end if
       xmaxocc=ieffthr*1.d-3

! reading fragment definition
       if(idofr.eq.1) then
       call locate(1,"# FRAGMENTS",ii)
       read(1,*) icufr
       do i=1,icufr
          read(1,*) nfrlist(i)
          if(i.eq.icufr.and.nfrlist(i).eq.-1) then 
            do l=1,nat
             navect(l)=0
            end do 
            do l=1,(icufr-1)
               do k=1,nfrlist(l)    
                 navect(ifrlist(k,l))=1
               end do 
            end do
            k=0
            do l=1,nat                     
              if(navect(l).eq.0) then
                k=k+1
                ifrlist(k,icufr)=l
              end if 
            end do  
            nfrlist(icufr)=k       
          else     
           read(1,*) (ifrlist(k,i),k=1,nfrlist(i))
          end if
        end do
        ixx=0
        do i=1,icufr
         write(*,*)'Fragment: ',i 
         write(*,'(20i4)') (ifrlist(k,i),k=1,nfrlist(i))
         ixx=ixx+nfrlist(i)
        end do
        if(ixx.ne.nat) stop 'Missing/Add atoms in fragment definition'

!  jfrlist tells which fragment a given atom belongs to
        do i=1,icufr
         do k=1,nfrlist(i)
          jfrlist(ifrlist(k,i))=i
         end do
        end do
        do i=1,nat
        if(jfrlist(i).eq.0) then
         write(*,*) 'Unassigned atom to fragment:',i
         stop
        end if
      end do
       
      allocate (qfrag(icufr),bofrag(icufr,icufr))

      else
! atom by atom
       icufr=nat
       do i=1,icufr
        nfrlist(i)=1
        ifrlist(1,i)=i
        jfrlist(i)=i
       end do

      end if
      close(1)
      
! end reading input file

! read Aimall's  int files. New first capital letter for atomname
! the following can be usefull:
! for i in *.int; do mv "$i" "`echo $i | sed "s/^./\U&/"`"; done

       if(iaimall.eq.1) then
! loop over atoms
        do iat=1,nat
! trimming blanks
         write(name,'(A2,i3,A4)') symb(iat),iat,".int"
         i=1
         do while (i.lt.len(trim(name)))
          if(index(name(i:i)," ").ne.0) then
           do j=i+1,len(name)
             name(j-1:j-1)=name(j:j)
           end do
          else
           i=i+1
          end if
         end do

         open(1, file=trim(name),status="OLD")
         write(*,*) "Reading AOM file ",trim(name)
! processing first atom...assuming params of remaining ones
        if(iat.eq.1) then
! WF type. Overrinding inf from fchk file
         iflag=1
         do while (iflag.eq.1)
          read(1,'(a80)')  line
          if(index(line,"Model: ").ne.0) iflag=0
         end do
         if(index(line,"Restricted CASSCF").ne.0) then
          write(*,*) 'DATA FROM RESTRICTED CASSCF ' 
          iuhf=0
          icas=1
         else if(index(line,"Restricted").ne.0) then
          write(*,*) 'DATA FROM RESTRICTED SCF CALCULATION ' 
          iuhf=0
         else if (index(line,"Unrestricted").ne.0) then 
          write(*,*) 'DATA FROM UNRESTRICTED SCF CALCULATION ' 
          iuhf=1
         else
          stop 'UNKNOWN WAVEFUNCTION ' 
         end if
! Number of electrons
         iflag=1
         do while (iflag.eq.1)
          read(1,'(a80)')  line
          if(index(line,"Number of Alpha electrons").ne.0) iflag=0
         end do
         read(line(58:),*) xalpha
         nalpha=nint(xalpha)
         read(1,'(a80)') line
         read(line(58:),*) xbeta
         nbeta=nint(xbeta)
         nmo=nalpha
         if(iuhf.eq.1.or.icas.eq.1) nmo=nalpha+nbeta !warning nonsinglet CASSCF
         allocate (ss(nmo,nmo,nat))

! reading NO occupations if CASSCF
         if(icas.eq.1) then
          iflag=1
          rewind(1)
          do while (iflag.eq.1)
           read(1,'(a80)')  line
           if(index(line,"Molecular Orbital (MO) Data").ne.0) iflag=0
          end do
          do i=1,11
          read(1,'(a80)') line
          end do
          allocate (occ(nmo))
          write(*,*) 'Number of Natural Orbitals:',nmo
           do inmo=1,nmo
            read(1,*)idummy,ino,occ(ino)
           end do
         end if
        end if

! reading atomic overlap matrices

         iflag=1
         do while (iflag.eq.1)
          read(1,'(a80)')  line
          if(index(line,"The Atomic Overlap Matrix: ").ne.0) iflag=0
         end do
         read(1,'(a80)') line
         read(1,'(a80)') line
         read(1,'(a80)') line
         read(1,*)((ss(i,j,iat),j=1,i),i=1,nmo)
         close(1)
        end do

! building s matrices
! alpha-beta block ignored at present
        allocate (sa(nmo,nmo,nat),q(nat),qa(nat),bo(nat,nat))
        allocate (sb(nmo,nmo,nat),qb(nat),qs(nat))
        sa=0.0d0
        sb=0.0d0

        do iat=1,nat
         do i=1,nmo
          do j=1,i
           ss(j,i,iat)=ss(i,j,iat)
          end do
         end do
         do i=1,nalpha
          do j=1,nalpha
           sa(i,j,iat)=ss(i,j,iat)
          end do
         end do
         if(iuhf.eq.1) then
          do i=1,nbeta 
           do j=1,nbeta
            sb(i,j,iat)=ss(i+nalpha,j+nalpha,iat)
           end do
          end do
         end if
         if(icas.eq.1) then
          do i=1,nmo   
           do j=1,nmo
            sa(i,j,iat)=ss(i,j,iat)
           end do
          end do
         end if
        end do

        deallocate (ss)
       end if
!end reading int files
        write(*,*) 'done with aimall input'

! read  Multiwfn AOM file
       if(imultiwfn.eq.1) then
       name=trim(name0)//".aom"
       open(1, file=name,status="OLD")
       read(1,'(a80)')  line
       if (index(line,"Atomic").ne.0) then
        write(*,*) 'RESTRICTED WFN '
        iuhf=0
        nalpha=nmo
        nbeta=nmo
       else if (index(line,"Alpha part of").ne.0) then
        iuhf=1
        write(*,*) 'UNRESTRICTED WFN '
        allocate (sb(nmo,nmo,nat),qb(nat),qs(nat))
! figuring out the spin state
        iflag=1
        read(1,'(a80)')  line
        do while (iflag.eq.1)
         read(1,'(a80)')  line
         if (index(line(1:7),"       ").ne.0) then
          iflag=0
          nalpha=idum
          nbeta=nmo-nalpha
! setting nmo to nalpha
          nmo=nalpha
          rewind(1)
          read(1,'(a80)')  line
         else
          read(line(1:7),*) idum
         end if
        end do
       else
        write(*,*) 'UNKNOWN WFN -- stop'
        stop
       end if

       write(*,*) 'NUMBER OF ALPHA ELECTRONS :',nalpha
       write(*,*) 'NUMBER OF BETA  ELECTRONS :',nbeta 
       do iat=1,nat
        icharge=icharge+zn(iat)
       end do
       icharge=icharge-nalpha-nbeta
       write(*,*) 'TOTAL CHARGE              :',icharge

! reading AOM atom by atom
       nbloc=nmo/5
       if(mod(nmo,5).ne.0) nbloc=nbloc+1
         
       do iat=1,nat
        read(1,'(a80)')  line
        do ibloc=1,nbloc
         j1=1+5*(ibloc-1)
         do imo=j1,nmo
          j2=min(imo,5*ibloc)
          read(1,'(I6,5(3x,f11.8))') idum,(sa(imo,j,iat),j=j1,j2)
         end do
         read(1,'(a80)')  line
        end do
!
        if (iuhf.eq.1) then
         read(1,'(a80)')  line
         if (index(line,"Beta").ne.0) stop 'error reading AOM file'
         do ibloc=1,nbloc
          j1=1+5*(ibloc-1)
          do imo=j1,nbeta
           j2=min(imo,5*ibloc)
           read(1,'(I6,5(3x,f11.8))') idum,(sb(imo,j,iat),j=j1,j2)
          end do
          read(1,'(a80)')  line
         end do
        end if

        if(iat.ne.nat) read(1,'(a80)')  line
       end do
! end reading AOM
       close(1)

! building s matrices
! unfortunately no alpha-beta bloc, so no local spin is possible
        do iat=1,nat
         do i=1,nmo
          do j=1,i-1
           sa(j,i,iat)=sa(i,j,iat)
           if(iuhf.eq.1) sb(j,i,iat)=sb(i,j,iat)
          end do
         end do
        end do

       end if
! printint occupatins with full accuracy for eos_edf
        name=trim(name0)//".occ"
        open(3, file=name)
        write(3,*) icufr,nalpha,nbeta,iuhf 
        do i=1,icufr
         ix=0
         do k=1,nfrlist(i)
          ix=ix+zn(ifrlist(k,i))
         end do
         ifchar(i)=ix
        end do
        write(3,*) (ifchar(ii),ii=1,icufr)

! CALCULATING DESCRIPTORS:
      print *,' '
      print *,' ---------------------------'
      print *,'  DOING POPULATION ANALYSIS '
      print *,'   Partial atomic charges'
      print *,'    Atomic spin densities '
      print *,'    Localization and  '
      print *,'   Delocalization indices '
      print *,' ---------------------------'
      print *,' '

! ATOMIC POPULATIONS
       xx0=0.0d0
       xx1=0.0d0
       xx2=0.0d0
       do iat=1,nat
        qa(iat)=0.0d0
        do i=1,nmo
         if(icas.eq.1) then
          qa(iat)=qa(iat)+occ(i)*sa(i,i,iat)*0.5d0  ! only for singlet CASCF
         else
          qa(iat)=qa(iat)+sa(i,i,iat)
         end if
        end do
        q(iat)=qa(iat)
        if(iuhf.eq.1) then
         qb(iat)=0.0d0
         do i=1,nbeta
          qb(iat)=qb(iat)+sb(i,i,iat)
         end do
         do i=1,nat
          q(iat)=qa(iat)+qb(iat)
          qs(iat)=qa(iat)-qb(iat)
         end do
        else
         q(iat)=2.0d0*qa(iat) 
        end if
        xx0=xx0+q(iat)
        xx1=xx1+qa(iat)
        if(iuhf.eq.1) xx2=xx2+qb(iat)
       end do
       
       write(*,*) 
       write(*,*) '**********************'  
       write(*,*) '* ATOMIC POPULATIONS *'  
       write(*,*) '**********************'  
       write(*,*) 

       if(iuhf.eq.1) then
        write(*,*) ' ATOM      NA        NB        N         Q      NA-NB   '  
        write(*,*) '--------------------------------------------------------'  
        do iat=1,nat
         write(*,'(i3,x,A2,x,5f10.6)') iat,symb(iat),qa(iat),qb(iat),q(iat),zn(iat)-q(iat),qs(iat)
        end do
        write(*,*) '--------------------------------------------------------' 
         write(*,'(6x,5f10.5)') xx1,xx2,xx1+xx2,icharge+nalpha+nbeta-xx1-xx2,xx1-xx2
       else
        write(*,*) ' ATOM       N       Q     '  
        write(*,*) '--------------------------'  
        do iat=1,nat
         write(*,'(i3,x,A2,x,2f10.6)') iat,symb(iat), q(iat),zn(iat)-q(iat)
        end do
        write(*,*) '--------------------------'  
         write(*,'(6x,2f10.5)') xx0, icharge+nalpha+nbeta-xx0
       end if

       if (idofr.eq.1) then
        write(*,*) '  '
        write(*,*) '   FRAGMENT ANALYSIS : Electron populations'
        xx=0.0d0
        qfrag=0.0d0
        do i=1,nat
         if(jfrlist(i).ne.0) qfrag(jfrlist(i))=qfrag(jfrlist(i))+q(i)
         xx=xx+q(i)
        end do
        write(*,*) '  '
        write(*,*) '  FRAG     Population'
        write(*,*) '------------------------'  
        do ifrag=1,icufr
         write(*,'(3x,i2,5x,f10.6)') ifrag,qfrag(ifrag)
        end do
        write(*,*) '------------------------'  
        write(*,'(a10,f10.6)') '   Sum :  ',xx
!
        write(*,*) '  '
        write(*,*) '   FRAGMENT ANALYSIS : Atomic charges'
        do i=1,nat
         if(jfrlist(i).ne.0) qfrag(jfrlist(i))=qfrag(jfrlist(i))-zn(i)
         xx=xx-zn(i)
        end do
        write(*,*) '  '
        write(*,*) '  FRAG        Charge'
        write(*,*) '------------------------'  
        do ifrag=1,icufr
         write(*,'(3x,i2,5x,f10.6)') ifrag,-qfrag(ifrag)
        end do
        write(*,*) '------------------------'  
        write(*,'(a10,f10.6)') '   Sum :  ',-xx
        
!
        if(iuhf.eq.1) then
        write(*,*) '  '
        write(*,*) '   FRAGMENT ANALYSIS : Spin Density'
        xx=0.0d0
        qfrag=0-0d0
        do i=1,nat
         if(jfrlist(i).ne.0) qfrag(jfrlist(i))=qfrag(jfrlist(i))+qs(i)
         xx=xx+qs(i)
        end do
        write(*,*) '  '
        write(*,*) '  FRAG    Spin density'
        write(*,*) '------------------------'  
        do ifrag=1,icufr
         write(*,'(3x,i2,5x,f10.6)') ifrag,qfrag(ifrag)
        end do
        write(*,*) '------------------------'  
        write(*,'(a10,f10.6)') '   Sum :  ',xx
       end if

       end if

        
!        if(jfrlist(i).ne.0.and.jfrlist(j).ne.0) B(jfrlist(i),jfrlist(j))=B(jfrlist(i),jfrlist(j))+a(i,j)
! LI/DI  
       write(*,*) 
       write(*,*) '**************************'  
       write(*,*) '*     LOCALIZATION /     *'
       write(*,*) '* DELOCALIZATION INDICES *' 
       write(*,*) '**************************'  
       write(*,*) 
!
       do iat=1,nat
        do jat=1,iat
         xx=0.0d0
         do i=1,nmo
          do j=1,nmo
           if(icas.eq.1) then
           xx=xx+sa(i,j,iat)*sa(j,i,jat)*occ(i)*occ(j)
           else
           xx=xx+sa(i,j,iat)*sa(j,i,jat)
           end if
          end do
         end do
         if(icas.eq.0) then
          bo(iat,jat)=2.0d0*xx
         else
          bo(iat,jat)=xx
         end if

         if(iuhf.eq.1) then
          xx=0.0d0
          do i=1,nbeta
           do j=1,nbeta
            xx=xx+sb(i,j,iat)*sb(j,i,jat)
           end do
          end do
          bo(iat,jat)=bo(iat,jat)+2.0d0*xx        
         else
          bo(iat,jat)=2.0d0*bo(iat,jat)
         end if
        end do
       end do
       do iat=1,nat
        do jat=1,iat
          bo(jat,iat)=bo(iat,jat) 
          if(iat.eq.jat) bo(iat,iat)=bo(iat,iat)*0.5d0
        end do
       end do
!
       call mprint(bo,nat,nat,symb)

      iprint=0
      if(iprint.eq.1) then
       write(*,*) ' ATOM PAIR       DI/LI     '  
       write(*,*) '---------------------------'  
       do iat=1,nat
        do jat=iat,nat 
         write(*,'(i4,a4,i4,f12.6)') iat,'   -',jat,bo(iat,jat)
        end do
       end do
       write(*,*) '---------------------------'  
       end if

  
!  EOS
       allocate (scr(nmo,nmo),c(nmo,nmo))
       allocate (p0net(nmo,icufr),ip0(icufr))
       allocate (oxi(icufr),elec(icufr))

        if(icas.eq.1) then
         do iat=1,nat
         do i=1,nmo
          do j=1,nmo
           if(occ(i)*occ(j).le.0.0d0) then
            sa(i,j,iat)=0.0d0
           else
            sa(i,j,iat)=dsqrt(occ(i)*occ(j))*sa(i,j,iat)*0.5d0 !warning only for singlet CASSCF
           end if
          end do
         end do
         end do
        end if


       print *,' '
       print *,' --------------------------------------'
       print *,'   DOING EFFAO-3D GENERAL FORMULATION'
       print *,' --------------------------------------'
       print *,' '
       print *,'     ALPHA  PART  '

! loop over fragments
       do ifrag=1,icufr

        scr=0.0d0
        do iat=1,nfrlist(ifrag)
         do i=1,nmo
          do j=1,nmo
           scr(i,j)=scr(i,j)+sa(i,j,ifrlist(iat,ifrag))
          end do
         end do
        end do
        call diagonalize(nmo,nmo,scr,c,0)

        imaxo=0
        do i=1,nmo
         if(scr(i,i).gt.xmaxocc) imaxo=imaxo+1
        end do

        xmaxo=0.0d0
        do i=1,imaxo
         xmaxo=xmaxo+scr(i,i)
        end do

        xx1=0.0d0
        do icenter=1,nfrlist(ifrag)
         xx1=xx1+qa(ifrlist(icenter,ifrag))
        end do
        write(*,*) '                  '
        write(*,'(a,i4,a)') ' ** FRAGMENT ',ifrag,' ** '
        write(*,*) '                  '

        write(*,'(a34,i3,f10.5)') ' Sum of occupations for fragment: ',ifrag, xmaxo
        if(iuhf.eq.0.and.icas.eq.0) then
         write(*,'(a34,f8.4)') ' Deviation from total population: ', 2.0d0*(xmaxo-xx1)
        else
         write(*,'(a34,f8.4)') ' Deviation from alpha population: ', xmaxo-xx1
        end if
        write(*,'(a27,f6.4)') ' EFOs occupations using >  ',xmaxocc
        write(*,60) (scr(mu,mu),mu=1,imaxo)
        write(3,*) imaxo,(scr(mu,mu),mu=1,imaxo)

! saving efo info for fragment
        do k=1,imaxo           
!         do mu=1,igr                       
!           p0(mu,k)=c(mu,k)
!         end do
         p0net(k,ifrag)=scr(k,k)
        end do
        ip0(ifrag)=imaxo

!  write cube file
!        if(icube.eq.1) call cubegen4(ifrag,icase)
       end do

! Doing EOS analysis of the ALPHA PART
      icase=1
      iimax=0
      do i=1,icufr
       iimax=iimax+ip0(i)
      end do
      allocate (occup(iimax,2),iorbat(iimax,2))
      allocate (occup2(iimax), occupg(iimax))

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
      write(*,*) ' Total number of EFOs for analysis:  ',iorb
      write(*,*) ' Total number of electrons to assign:',nalpha

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
      nnn=nalpha
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
       write(*,30) kk+1,kk+k+1
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

      write(*,*)'  '
      write(*,*)' EOS ANALYSIS  FOR ALPHA ELECTRONS'
      write(*,*)'  '

      print *,'  Fragm   Elect  Last occ. First unocc '
      print *,' -------------------------------------'

       xlast=1.0d0
       ilast=0
       do i=1,icufr
        nn=int(elec(i))
        if(elec(i)-nn.gt.thres) nn=nn+1
        if(nn+1.gt.ip0(i)) then 
         write(*,100) i,elec(i), p0net(nn,i)
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
       end do 

       xfirst=0.0d0
       do i=1,icufr
        if(i.ne.ilast) then
        nn=int(elec(i))
        if(elec(i)-nn.gt.thres) nn=nn+1
        if(nn+1.ne.0.and.p0net(nn+1,i).gt.xfirst) xfirst= p0net(nn+1,i)
        end if
       end do 

      print *,' -------------------------------------'
      confi=100.0*min(1.0d0,xlast-xfirst+0.5d0)
      write(*,'(a30,f8.3)') 'RELIABILITY INDEX R(%) =',confi

! PSS quadratic deviation of the fractional and integer occupations as new fragment indicator
      print *
      print *,'(warning in case of fractional occupations)'
      print *,'  Fragm    Quad Dev '
      print *,' -------------------'
      do i=1,icufr
       nn=int(elec(i))
       if(elec(i)-nn.gt.thres) nn=nn+1
       qdev(i)=0.0d0
       do k=1,ip0(i)
        xx=p0net(k,i)
        if(k.le.nn) xx=1.0d0-p0net(k,i)
        qdev(i)=qdev(i)+xx*xx
!        write(*,*) 'pollas',i,ip0(i),p0net(k,i),xx
       end do
       qdev(i)=sqrt(qdev(i))
       write(*,'(I4,5x,f8.4)') i, qdev(i)
      end do
      print *,' -------------------'
! PSS 

      confi0=confi

       zztot=0.0d0
       do i=1,icufr
        zzn=0.0d0
        do j=1,nfrlist(i)
         zzn=zzn+zn(ifrlist(j,i))
        end do
        if(iuhf.eq.0) then
         oxi(i)=zzn-2.0d0*elec(i)
         zztot=zztot+oxi(i)
        else
         oxi(i)=zzn-elec(i)
        end if
       end do

! output for eos_edf
        do ii=1,icufr
         write(3,*) (1,i=1,elec(ii)),(0,i=elec(ii)+1,imaxo)
        end do

! now beta part if unrestricted
! nmo > nbeta but it's padded with zeroes
       if(iuhf.eq.0) then
        write(*,*)'  '
        write(*,*)' SKIPPING EFOs FOR BETA ELECTRONS'
        write(*,*)'  '
       else
        print *,' '
        print *,'     BETA  PART  '
        if(nbeta.eq.0) then
         print *,' '
         print *,' Calculation has no beta electrons'
         go to 999
        end if
        icase=2

        write(3,*) 'fragment occupations for beta  part' 
! loop over fragments
       do ifrag=1,icufr

        scr=0.0d0
        do iat=1,nfrlist(ifrag)
         do i=1,nmo
          do j=1,nmo
           scr(i,j)=scr(i,j)+sb(i,j,ifrlist(iat,ifrag))
          end do
         end do
        end do
        call diagonalize(nmo,nmo,scr,c,0)

        imaxo=0
        do i=1,nmo
         if(scr(i,i).gt.xmaxocc) imaxo=imaxo+1
        end do

        xmaxo=0.0d0
        do i=1,imaxo
         xmaxo=xmaxo+scr(i,i)
        end do

        xx1=0.0d0
        do icenter=1,nfrlist(ifrag)
         xx1=xx1+qb(ifrlist(icenter,ifrag))
        end do
        write(*,*) '                  '
        write(*,'(a,i4,a)') ' ** FRAGMENT ',ifrag,' ** '
        write(*,*) '                  '

        write(*,'(a34,i3,f10.5)') ' Sum of occupations for fragment: ',ifrag, xmaxo
        write(*,'(a34,f8.4)') ' Deviation from beta population:  ', xmaxo-xx1
        write(*,'(a26,f6.4)') ' EFOs occupations using > ',xmaxocc
        write(*,60) (scr(mu,mu),mu=1,imaxo)
        write(3,*) imaxo,(scr(mu,mu),mu=1,imaxo)

! saving efo info for fragment
        do k=1,imaxo           
!         do mu=1,igr                       
!           p0(mu,k)=c(mu,k)
!         end do
         p0net(k,ifrag)=scr(k,k)
        end do
        ip0(ifrag)=imaxo
!  write cube file
!        if(icube.eq.1) call cubegen4(ifrag,icase)

       end do

      deallocate (occup,iorbat,occup2,occupg)
! Doing EOS analysis of the BETA PART
      icase=2
      iimax=0
      do i=1,icufr
       iimax=iimax+ip0(i)
      end do
      allocate (occup(iimax,2),iorbat(iimax,2))
      allocate (occup2(iimax), occupg(iimax))

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
      write(*,*) ' Total number of EFOs for analysis:  ',iorb
      write(*,*) ' Total number of electrons to assign:',nbeta 

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
      nnn=nbeta
335   k=k+1
      if(dabs(occup(nnn,icase)-occup(nnn+k,icase)).lt.thres) go to 335
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
336    kk=kk+1
       if(dabs(occup(nnn,icase)-occup(nnn-kk,icase)).lt.thres) go to 336
       kk=kk-1
       frac=float(kk+1)/float(kk+k+1)
       write(*,30) kk+1,kk+k+1
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

      write(*,*)'  '
      write(*,*)' EOS ANALYSIS  FOR BETA ELECTRONS'
      write(*,*)'  '

      print *,'  Fragm   Elect  Last occ. First unocc '
      print *,' -------------------------------------'

       xlast=1.0d0
       ilast=0
       do i=1,icufr
        nn=int(elec(i))
        if(elec(i)-nn.gt.thres) nn=nn+1
        if(nn+1.gt.ip0(i)) then 
         write(*,100) i,elec(i), p0net(nn,i)
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
       end do 

       xfirst=0.0d0
       do i=1,icufr
        if(i.ne.ilast) then
        nn=int(elec(i))
        if(elec(i)-nn.gt.thres) nn=nn+1
        if(nn+1.ne.0.and.p0net(nn+1,i).gt.xfirst) xfirst= p0net(nn+1,i)
        end if
       end do 

      print *,' -------------------------------------'
      confi=100.0*min(1.0d0,xlast-xfirst+0.5d0)
      write(*,'(a30,f8.3)') 'RELIABILITY INDEX R(%) =',confi

! PSS quadratic deviation of the fractional and integer occupations as new fragment indicator
      print *
      print *,'(warning in case of fractional occupations)'
      print *,'  Fragm    Quad Dev '
      print *,' -------------------'
      do i=1,icufr
       nn=int(elec(i))
       if(elec(i)-nn.gt.thres) nn=nn+1
       qdev(i)=0.0d0
       do k=1,ip0(i)
        xx=p0net(k,i)
        if(k.le.nn) xx=1.0d0-p0net(k,i)
        qdev(i)=qdev(i)+xx*xx
!        write(*,*) 'pollas',i,ip0(i),p0net(k,i),xx
       end do
       qdev(i)=sqrt(qdev(i))
       write(*,'(I4,5x,f8.4)') i, qdev(i)
      end do
      print *,' -------------------'
! PSS 
       zztot=0.0d0
       do i=1,icufr
         oxi(i)=oxi(i)-elec(i)
         zztot=zztot+oxi(i)
       end do

      end if
! output for eos_edf
!        do ii=1,icufr
!         write(3,*) (1,i=1,elec(ii)),(0,i=elec(ii)+1,imaxo)
!        end do

999   continue
! FINAL EOS RESULTS
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


100   format(i4,2x,f6.2,f12.3,'    < thresh',2f12.3)
15    format(i4,2x,f6.2,2f12.3)
20    format(i4,2x,f10.2)
25    format(i4,2x,f10.4,2x,f10.4)
30    format("Distributing ",i4," electrons over",i4,"pseudodegenerate EFOs")
60    format(7h OCCUP.   ,8F9.4)
      end 













!!!!!
!!!!!
        Subroutine diagonalize(M,N,A0,X,ival)
        implicit double precision(a-h,o-z)
        integer, intent(in) :: M,N
        dimension x(M,M),a0(M,M),d(M)


        call sdiag2(a0,m,n,d,ival)
        do i=1,n
         do j=1,n 
          x(j,i)=a0(j,i)
          a0(j,i)=0.d0
         end do
         a0(i,i)=d(i)
        end do
        return
        end

! DIAGONALIZATION OF THE REAL SYMMETRIC MATRIX X. IN D THE EIGENVALUES.
      Subroutine SDIAG2(X,M,N,D,inosort)
      IMPLICIT REAL*8(A-H,O-Z)
      integer, intent(in) :: M,N,inosort
      DIMENSION D(M)
      DIMENSION E(M)
      DIMENSION X(M,M)
      EPS=5.D-15
      TOL=5.D-40
      IF(N.EQ.1) GOTO 400
    5 DO 150 NI=2,N
      II=N+2-NI
      DO 150 I=II,II
      L=I-2
      H=0.0D0
      G=X(I,I-1)
      IF(L)140,140,20
   20 DO 30 K=1,L
   30 H=H+X(I,K)*X(I,K)
      S=H+G*G
      IF(S.GE.TOL) GOTO 50
   40 H=0.0D0
      GO TO 140
   50 IF(H)140,140,60
   60 L=L+1
      F=G
      G=DSQRT(S)
      IF(F)75,75,70
   70 G=-G
   75 H=S-F*G
      X(I,I-1)=F-G
      F=0.0D0
      DO 110 J=1,L
      X(J,I)=X(I,J)/H
      S=0.0D0
      DO 80 K=1,J
   80 S=S+X(J,K)*X(I,K)
      J1=J+1
      IF(J1.GT.L) GOTO 100
   85 DO 90 K=J1,L
   90 S=S+X(K,J)*X(I,K)
  100 E(J)=S/H
  110 F=F+S*X(J,I)
      F=F/(H+H)
      DO 120 J=1,L
  120 E(J)=E(J)-F*X(I,J)
      DO 130 J=1,L
      F=X(I,J)
      S=E(J)
      DO 130 K=1,J
  130 X(J,K)=X(J,K)-F*E(K)-X(I,K)*S
  140 D(I)=H
  150 E(I-1)=G
  160 D(1)=X(1,1)
      X(1,1)=1.0D0
      DO 220 I=2,N
      L=I-1
      IF(D(I))200,200,170
  170 DO 190 J=1,L
      S=0.0D0
      DO 180 K=1,L
  180 S=S+X(I,K)*X(K,J)
      DO 190 K=1,L
  190 X(K,J)=X(K,J)-S*X(K,I)
  200 D(I)=X(I,I)
      X(I,I)=1.0D0
  210 DO 220 J=1,L
      X(I,J)=0.0D0
  220 X(J,I)=0.0D0
      B=0.0D0
      F=0.0D0
      E(N)=0.0D0
      DO 340 L=1,N
      H=EPS*(DABS(D(L))+DABS(E(L)))
      IF(H.GT.B) B=H
  235 DO 240 J=L,N
      IF(DABS(E(J)).LE.B) GO TO 250
  240 CONTINUE
  250 IF(J.EQ.L) GOTO 340
  260 P=(D(L+1)-D(L))*.50D0/E(L)
      R=DSQRT(P*P+1.0D0)
      IF(P)270,280,280
  270 P=P-R
      GO TO 290
  280 P=P+R
  290 H=D(L)-E(L)/P
      DO 300 I=L,N
  300 D(I)=D(I)-H
      F=F+H
      P=D(J)
      C=1.0D0
      S=0.0D0
      J1=J-1
      DO 330 NI=L,J1
      II=L+J1-NI
      DO 330 I=II,II
      G=C*E(I)
      G=C*E(I)
      H=C*P
      IF(DABS(P).LT.DABS(E(I))) GO TO 310
  305 C=E(I)/P
      R=DSQRT(C*C+1.0D0)
      E(I+1)=S*P*R
      S=C/R
      C=1.0D0/R
      GO TO 320
  310 C=P/E(I)
      R=DSQRT(C*C+1.0D0)
      E(I+1)=S*E(I)*R
      S=1.0D0/R
      C=C/R
  320 P=C*D(I)-S*G
      D(I+1)=H+S*(C*G+S*D(I))
      DO 330 K=1,N
      H=X(K,I+1)
      X(K,I+1)=X(K,I)*S+H*C
  330 X(K,I)=X(K,I)*C-H*S
      E(L)=S*P
      D(L)=C*P
      IF(DABS(E(L)).GT.B) GO TO 260
  340 D(L)=D(L)+F
      NI=N-1
      if(inosort.eq.1) goto 400
  350 DO 380 I=1,NI
      K=I
      P=D(I)
      J1=I+1
      DO 360 J=J1,N
      IF(D(J).LE.P) GOTO 360
  355 K=J
      P=D(J)
  360 CONTINUE
      IF(K.EQ.I) GOTO 380
  365 D(K)=D(I)
      D(I)=P
      DO 370 J=1,N
      P=X(J,I)
      X(J,I)=X(J,K)
  370 X(J,K)=P
  380 CONTINUE
  390 GO TO 410
  400 D(1)=X(1,1)
      X(1,1)=1.0D0
  410 RETURN
      END
!

!*******************************
! PROCESSING INP FILE   ********
!*******************************
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

      SUBROUTINE MPRINT(H,N,ndim,atsymb)
      parameter(maxat=500,maxfrag=50)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer, intent(in) :: ndim
      character*2 atsymb(maxat)
      DIMENSION H(NDIM,NDIM)

      K=6
      NMIN=1
      NNMAX=MIN0(N,K)
   62 FORMAT(1X,I3,A4,6F12.6)
   1  PRINT 60, (I, atsymb(i),I=NMIN,NNMAX)
   60 FORMAT(10X,6(2X,I3,A4,3X))
      PRINT 64
   64 FORMAT(1X)
      DO 2 I=1,N
      PRINT 62,I,atsymb(i),(H(I,J),J=NMIN,NNMAX)
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

