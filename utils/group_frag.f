      program Agrupfrag

c     The program group by fragments the different energy contributions

      implicit double precision (a-h,o-z)
      character*80 cfile,ofile,ifile,tfile,line
      character*7 cthrash,flag1 
      parameter (maxat=200)
      dimension navect(maxat),nfrlist(maxat),ifrlist(maxat,maxat) 
      dimension enxc(maxat,maxat),jfrlist(maxat),enkin(maxat,maxat) 
      dimension fragxc(maxat,maxat),fragelen(maxat,maxat) 
      dimension fragkin(maxat,maxat),enelen(maxat,maxat) 
      dimension fragcoul(maxat,maxat),encoul(maxat,maxat) 
      dimension fragtot(maxat,maxat),entot(maxat,maxat) 
      dimension fragbond(maxat,maxat),enbond(maxat,maxat) 
      dimension fragnuc(maxat,maxat),fragstat(maxat,maxat) 

      call getarg(1,cfile)
      call getarg(2,flag1)
      iall=0
      if(index(flag1,"-all").eq.1) then
       iall=1
      end if

      i=1
      l=0
      do while(cfile(i:i).ne.' ')
       i=i+1
      end do
      l=i-1

      tfile=cfile(1:l)//".inp"
      ifile=cfile(1:l)//".tfvc"
      ofile=cfile(1:l)//".OUTEN"

      open(unit=1,file=tfile)
      open(unit=2,file=ifile)
      open(unit=3,file=ofile)

      rewind(1)
      rewind(2)
      rewind(3)

c     Reading input
      
      iecp=0
      read(2,'(a80)') line 
      do while(index(line,"Number of atoms").eq.0)
       read(2,'(a80)') line 
      end do 
      read(line(20:80),*) natoms
      rewind(2)
      read(2,'(a80)') line 
      do while(index(line,"Normal Termina").eq.0)
       read(2,'(a80)') line 
       if(index(line,"ECP ").ne.0) then
        iecp=1
       end if
      end do 
      rewind(2)

      write(3,*) "Starting to read the input file "
      idofr=0
      call readchar ("# METODE","DO FRAGS",idofr)
      if(idofr.eq.1) then
       call locate(1,"# FRAGMENTS",ii)
       if(ii.eq.0) then
        write(3,*) "DOFRAGS section missing in the .inp file! "
        stop
       end if
       read(1,*) icufr
       do i=1,icufr
        read(1,*) nfrlist(i)
        if(i.eq.icufr.and.nfrlist(i).eq.-1) then 
         do l=1,natoms
          navect(l)=0
         end do 
         do l=1,(icufr-1)
          do k=1,nfrlist(l)    
           navect(ifrlist(k,l))=1
          end do 
         end do
         k=0
         do l=1,natoms                     
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
      end if

      write(3,*) "Input file succesfully read "

c  jfrlist tells which fragment a given atom belongs to

      do i=1,icufr
       do k=1,nfrlist(i)
        jfrlist(ifrlist(k,i))=i
       end do
      end do

c     Printing fragment information

      write(3,*) " "
      write(3,*) "### Fragment definition ### "
      write(3,*) " "
      do i=1,icufr
       write(3,*) "Fragment:",i
       write(3,*) "Fragment number of atoms:",nfrlist(i)
       write(3,*) "List of atoms:",(ifrlist(k,i),k=1,nfrlist(i))
       write(3,*) " "
      end do
      write(3,*) " "

c     Reading TFVC 

      write(3,*) "Starting to read the TFVC file "

c     Electron-Nuclei energy contributions

c BOND ORDER AND TOTAL ENERGY ALWAYS. -all FLAG REQUIRED AS ARGUMENT 2
c IF ALL THE CONTRIBUTIONS ARE DESIRED.
      
c falta si no hi ha pseudopotencial
c refer el codi en funcio del nou output. (total i bond order no)

      ll=natoms/6
      if(ll*6.ne.natoms) ll=ll+1

      if(iall.eq.1) then

       rewind(2)
       if(iecp.eq.1) then
        read(2,'(a80)') line 
        do while(index(line,"Adding ECP").eq.0)
         read(2,'(a80)') line 
        end do 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        ii2=1
        do ii=1,ll 
         if(ii.ne.ll) then
          do jj=1,natoms
           read(2,*) ithrash,cthrash,(enelen(jj,kk),kk=ii2,ii2+5)
          end do 
         else 
          do jj=1,natoms
           read(2,*) ithrash,cthrash,(enelen(jj,kk),kk=ii2,natoms)
          end do 
         end if
         ii2=ii2+6
         read(2,'(a80)') line 
         read(2,'(a80)') line 
         read(2,'(a80)') line 
         read(2,'(a80)') line 
         read(2,'(a80)') line 
         read(2,'(a80)') line 
        end do 
       else 
c nopseudo part...
        continue
       end if 

c     Kinetic energy contributions

       rewind(2)
       read(2,'(a80)') line 
       do while(index(line,"tron-nuclei").eq.0)
        read(2,'(a80)') line 
       end do 
       do kk=1,17
        read(2,'(a80)') line 
       end do 
       ii2=1
       do ii=1,ll 
        if(ii.ne.ll) then
         do jj=1,natoms
          read(2,*) ithrash,cthrash,(enkin(jj,kk),kk=ii2,ii2+5)
         end do 
        else 
         do jj=1,natoms
          read(2,*) ithrash,cthrash,(enkin(jj,kk),kk=ii2,natoms)
         end do 
        end if
        ii2=ii2+6
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
       end do 
 
c     XC energy contributions
  
       rewind(2)
       read(2,'(a80)') line 
       do while(index(line,"Final ").eq.0)
        read(2,'(a80)') line 
       end do 
       read(2,'(a80)') line 
       read(2,'(a80)') line 
       ii2=1
       do ii=1,ll 
        if(ii.ne.ll) then
         do jj=1,natoms
          read(2,*) ithrash,cthrash,(enxc(jj,kk),kk=ii2,ii2+5)
         end do 
        else 
         do jj=1,natoms
          read(2,*) ithrash,cthrash,(enxc(jj,kk),kk=ii2,natoms)
         end do 
        end if
        ii2=ii2+6
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
       end do 
 
c     Coulomb energy contributions
  
       rewind(2)
       read(2,'(a80)') line 
       do while(index(line,"INTERPOLATED TW").eq.0)
        read(2,'(a80)') line 
       end do 
       read(2,'(a80)') line 
       read(2,'(a80)') line 
       ii2=1
       do ii=1,ll 
        if(ii.ne.ll) then
         do jj=1,natoms
          read(2,*) ithrash,cthrash,(encoul(jj,kk),kk=ii2,ii2+5)
         end do 
        else 
         do jj=1,natoms
          read(2,*) ithrash,cthrash,(encoul(jj,kk),kk=ii2,natoms)
         end do 
        end if
        ii2=ii2+6
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
        read(2,'(a80)') line 
       end do 

       end if

c     Total energy contributions
 
      rewind(2)
      read(2,'(a80)') line 
      do while(index(line,"DFT ENER").eq.0.and.index(line,"ock ENER").eq.0)
       read(2,'(a80)') line 
      end do 
      read(2,'(a80)') line 
      read(2,'(a80)') line 
      read(2,'(a80)') line 
      read(2,'(a80)') line 
      ii2=1
      do ii=1,ll 
       if(ii.ne.ll) then
        do jj=1,natoms
         read(2,*) ithrash,cthrash,(entot(jj,kk),kk=ii2,ii2+5)
        end do 
       else 
        do jj=1,natoms
         read(2,*) ithrash,cthrash,(entot(jj,kk),kk=ii2,natoms)
        end do 
       end if
       ii2=ii2+6
       read(2,'(a80)') line 
       read(2,'(a80)') line 
       read(2,'(a80)') line 
       read(2,'(a80)') line 
       read(2,'(a80)') line 
       read(2,'(a80)') line 
      end do 

c     Bond orders
 
      rewind(2)
      read(2,'(a80)') line 
      do while(index(line,"BOND ORDER ").eq.0)
       read(2,'(a80)') line 
      end do 
      read(2,'(a80)') line 
      read(2,'(a80)') line 
      read(2,'(a80)') line 
      read(2,'(a80)') line 
      ii2=1
      do ii=1,ll 
       if(ii.ne.ll) then
        do jj=1,natoms
         read(2,*) ithrash,cthrash,(enbond(jj,kk),kk=ii2,ii2+5)
        end do 
       else 
        do jj=1,natoms
         read(2,*) ithrash,cthrash,(enbond(jj,kk),kk=ii2,natoms)
        end do 
       end if
       ii2=ii2+6
       read(2,'(a80)') line 
       read(2,'(a80)') line 
       read(2,'(a80)') line 
       read(2,'(a80)') line 
       read(2,'(a80)') line 
       read(2,'(a80)') line 
      end do 

      write(3,*) "TFVC file succesfully read "
      write(3,*) " "

c     -----
c     Grouping by frags + calculating the nuclear repulsion

      if(iall.eq.1) then

       write(3,*) " Electron-Nuclei Energy by fragment analysis "
       call group_by_frag(1,natoms,icufr,jfrlist,enelen,fragelen)
       write(3,*) " Kinetic Energy by fragment analysis "
       call group_by_frag(1,natoms,icufr,jfrlist,enkin,fragkin)
       write(3,*) " XC Energy by fragment analysis "
       call group_by_frag(1,natoms,icufr,jfrlist,enxc,fragxc)
       write(3,*) " Coulomb Energy by fragment analysis "
       call group_by_frag(1,natoms,icufr,jfrlist,encoul,fragcoul)

      end if 

      write(3,*) " Total Calculated Energy by fragment analysis "
      call group_by_frag(1,natoms,icufr,jfrlist,entot,fragtot)
      write(3,*) " Total Calculated Bond Orders by fragment analysis "
      call group_by_frag(1,natoms,icufr,jfrlist,enbond,fragbond)

      if(iall.eq.1) then

       do ii=1,icufr
        do jj=1,icufr
         x1=fragelen(ii,jj)+fragkin(ii,jj)+fragxc(ii,jj)+fragcoul(ii,jj)
         fragnuc(ii,jj)=fragtot(ii,jj)-x1
        end do 
       end do      
       write(3,*) " Nuclear-Nuclear Energy by fragment analysis "
       write(3,*) " "
       do ii=1,icufr
        write(3,69) (fragnuc(ii,jj),jj=1,icufr)
       end do 

c     -----
c     Electrostatic energy by frags

       do ii=1,icufr
        do jj=1,icufr
         x1=fragelen(ii,jj)+fragcoul(ii,jj)+fragnuc(ii,jj)
         fragstat(ii,jj)=x1
        end do
       end do
       write(3,*) " "
       write(3,*) " Electrostatic Energy by fragment analysis "
       write(3,*) " "
       do ii=1,icufr
        write(3,69) (fragstat(ii,jj),jj=1,icufr)
       end do
       write(3,*) " "

      end if 

c     -----

      write(3,*) " End of Agrupfrag "

   69 FORMAT(6f12.6)
      close(1)
      close(2)
      close(3)
      end

c     *****

      subroutine group_by_frag(ilog,natoms,icufr,jfrlist,a,b)
      implicit double precision(a-h,o-z)
      parameter (maxat=200)
      dimension a(maxat,maxat),b(maxat,maxat)
      dimension jfrlist(maxat)

      do i=1,natoms
       do j=1,natoms
        b(i,j)=0.0d0
       end do
      end do

      nnat=natoms
      do i=1,natoms
       if(ilog.eq.1) nnat=i
       do j=1,nnat
        b(jfrlist(i),jfrlist(j))=b(jfrlist(i),jfrlist(j))+a(i,j)
       end do
      end do

      x=0.0d0
      if(ilog.eq.1) then
       do i=1,icufr
        x=x+b(i,i)
        do j=1,i-1
         b(i,j)=b(i,j)+b(j,i)
         x=x+b(i,j)
         b(j,i)=b(i,j)
        end do
       end do
      else
       do i=1,icufr
        do j=1,icufr
         x=x+b(i,j)
        end do
       end do
      end if

      write(3,*) ' '
      do i=1,icufr
       write(3,70) (b(i,j),j=1,icufr)
      end do 
      write(3,*) ' '
      write(3,71) ' Total contribution : ',x
      write(3,*) ' '

   70 FORMAT(6f19.12)
   71 FORMAT(a41,6f19.12)
      end

c     *****

      subroutine locate(iunit,string,ii)
      integer iunit
      character string*(*)
      character*80 linia

      rewind(iunit)
      ii=0
      do while(ii.eq.0)
       read(iunit,"(a80)")linia
       if(index(linia,string).ne.0) then
        ii=1
        return
       end if
      end do
      end

c     *****

      subroutine readchar(section,keyword,ival)
      character section*(*), keyword*(*)
      character linea*80
      integer ii,ilen,ilog

      ival=0
      call locate(1,section,ii)
      ii=0
      do while(ii.eq.0)
       read(1,"(a80)") linea
       if(index(linea,keyword).ne.0) then
        ii=1
        ival=1
       end if
       if(index(linea,"#").ne.0) then
        ii=2
       end if
      end do
      ival=1
      if(ii.eq.1) then
       write(3,"(a8,a2,i3)") keyword,'= ',ival
      end if
      end

c     *****

