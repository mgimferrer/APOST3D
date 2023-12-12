      program XXX

c     Program to know if the sum of the X contributions gives 
c     the total X energy. #P, iop(3/33=3) and pop=full required 
c     in the g09 input.

C reads .log file as argument and prints on screen (append to .fchk file)
 
      implicit double precision (a-h,o-z)
      integer b,d,i,j,k,bf
      parameter (maxdim=2000)
      dimension xkem(maxdim,maxdim),xpem(maxdim,maxdim)
      dimension xpm(maxdim,maxdim),xchm(maxdim,maxdim)
      dimension xecpm(maxdim,maxdim)
      character*80 line,name0,cc1,cc
           
      call getarg(1,name0)    
      open(unit=1,file=name0)

!! TO SEE IF THERE IS PSEUDOPOTENTIAL !!
 
      nn=0
      rewind(1)
      read(1,'(a80)') line
      do while(index(line,"Normal term").eq.0)
        read(1,'(a80)') line
        if(index(line,"ECP Int").ne.0)then
          nn=1
          go to 101
        end if
      end do
101   continue

!! EXTRACTING THE ENERGIES FROM log FILE !!

      rewind(1)
      read(1,'(a80)') line
      do while(index(line,"Normal term").eq.0)
        read(1,'(a80)') line
        if(index(line,"SCF Done").ne.0) read(line,*) cc,cc,cc,cc,xscf,cc,cc,icyc,cc
      end do 

      rewind(1)
      read(1,'(a80)') line
      do while(index(line,"Normal term").eq.0)
        read(1,'(a80)') line
        if(index(line," N-N=").eq.1) then
          read(line(6:25),*) xnn
          read(line(30:50),*) xen
          read(line(54:72),*) xkin
        end if
      end do

      xee=xscf-(xnn+xen+xkin)
      cc1="Kinetic Energy "
      write(*,'(A43,A6,ES22.15)') cc1,"R     ",xkin
      cc1="Electron-Nuclei Energy "
      write(*,'(A43,A6,ES22.15)') cc1,"R     ",xen
      cc1="Electron-Electron Energy "
      write(*,'(A43,A6,ES22.15)') cc1,"R     ",xee 

c     Reading matrices!

      rewind(1)
      read(1,'(a80)') line 
      do while(index(line,"NBasis").eq.0)
       read(1,'(a80)') line 
      end do 
      read(line(13:17),*) bf

c     P matrix
      
c      call readDmatrix(xpm,bf,"Density ")

c     Kinetic Energy
      
c     call readXmatrix(xkem,bf,"Kinetic ")
c      call rmat("Kinetic Energy Matrix ",bf,bf,xkem)
c      call calcenergy(xpm,xkem,bf,xke)
c      cc1="Kinetic Energy (Using KE Matrix)"
c      write(*,'(A43,A6,I12)') cc1,"R   N=",ii 
c      write(*,'(5ES16.8)') xke

c     Potential Energy

c     call readXmatrix(xpem,bf,"Potentia")
c      call rmat("Potential Energy Matrix ",bf,bf,xpem)
c      call calcenergy(xpm,xpem,bf,xpe)
c      cc1="Potential Energy (Using PE Matrix)"
c      write(*,'(A43,A6,I12)') cc1,"R   N=",ii 
c      write(*,'(5ES16.8)') xpe

c     Electron Core Potential

      if(nn.eq.1) then
       call readXmatrix(xecpm,bf,"ECP Inte")
       call rmat("ECP Matrix ",bf,bf,xecpm)
c       call calcenergy(xpm,xecpm,bf,xecp)
c       cc1="ECP Energy (Using ECP Integrals)"
c       write(*,'(A43,A6,I12)') cc1,"R   N=",ii 
c       write(*,'(5ES16.8)') xecp
      end if

c     Core Hamiltonian

c      call readXmatrix(xchm,bf,"Core Ham")
c      call rmat("Core Hamiltonian Matrix ",bf,bf,xchm)
c      call calcenergy(xpm,xchm,bf,xch)
c      write(*,*) "A saber... (Core Hamiltonian) = ",xch
c      write(*,*)" "                

      close(1)
      end

c ****       

      subroutine rmat(key,ival,jval,rmatrix)
      implicit double precision (a-h,o-z)
      character(len=*) ::  key
      integer ival,jval
      parameter (maxdim=10000)
      dimension rmatrix(maxdim,maxdim)
      character(len=43) ::  title

      title=adjustl(key)
      write(*,'(A43,A6,I12)') title,"R   N=",ival*(ival+1)/2
      write(*,'(5ES16.8)') ((rmatrix(i,j),j=1,i),i=1,ival)
      end

c ****

      subroutine readXmatrix(x,bf,str)
      implicit double precision (a-h,o-z)
      integer b,d,i,j,k,bf
      parameter (maxdim=10000)
      dimension x(maxdim,maxdim)
      character*80 line
      character*8 str

      rewind(1)

      read(1,'(a80)') line 
      do while(index(line,str).eq.0)
       read(1,'(a80)',end=100) line 
       if (index(line,"Normal Termination").eq.1) go to 100
      end do
                   
      b=bf/5
      if(B*5.ne.bf) then
       b=(bf/5)+1
      end if  
   
      d=1
      do k=1,b
       if(k.eq.b) then 
             
        read(1,'(a80)') line
        do  i=d,bf
         read(1,'(a80)') line
         read(line(9:80),*) (x(i,j),j=d,i)
        end do
       else  

c  De la fila 1 fins la 4 
            
        read(1,'(a80)') line
        do i=d,d+3
         read(1,'(a80)') line
         read(line(9:80),*) (x(i,j),j=d,i)
        end do
   
c  a partir de la fila 5 a la M 
              
        do i=d+4,bf         
         read(1,'(a80)') line
         read(line(9:80),*) (x(i,j),j=d,d+4)
        end do 
        d=d+5
       end if 
      end do 

100   continue

      end 

c ****

      subroutine calcenergy(xx,yy,bf,xenergy)
      implicit double precision (a-h,o-z)
      integer bf,n,m
      parameter (maxdim=10000)
      dimension xx(maxdim,maxdim),yy(maxdim,maxdim),zz(maxdim,maxdim)

      n=bf
      m=bf

c     Doing the full matrix 

      do i=1,n
       do j=1,m
        xx(i,j)=xx(j,i)
        yy(i,j)=yy(j,i)
       end do 
      end do 

      do i=1,n
       do j=1,m
        zz(i,j)=0.0d0
        do k=1,m
         zz(i,j)=zz(i,j)+(xx(i,k)*yy(k,j))
        end do
       end do
      end do
 
      xenergy=0.0d0
      do i=1,m
       xenergy=xenergy+zz(i,i)
      end do 
 
      end

c ****

      subroutine readDmatrix(x,bf,str)
      implicit double precision (a-h,o-z)
      integer b,d,i,j,k,bf
      parameter (maxdim=10000)
      dimension x(maxdim,maxdim)
      character*80 line
      character*8 str

      rewind(1)

      read(1,'(a80)') line 
      do while(index(line,str).eq.0)
       read(1,'(a80)') line 
       if (index(line,"Normal Termination").eq.1) go to 100
      end do
                   
      b=bf/5
      if(B*5.ne.bf) then
       b=(bf/5)+1
      end if  

      d=1
      do k=1,b
       if(k.eq.b) then 
             
        read(1,'(a80)') line
        do  i=d,bf
         read(1,'(a80)') line
         read(line(23:80),*) (x(i,j),j=d,i)
        end do
       else  

c  De la fila 1 fins la 4 
            
        read(1,'(a80)') line
        do i=d,d+3
         read(1,'(a80)') line
         read(line(23:80),*) (x(i,j),j=d,i)
        end do
   
c  a partir de la fila 5 a la M 
              
        do i=d+4,bf         
         read(1,'(a80)') line
         read(line(23:80),*) (x(i,j),j=d,d+4)
        end do 
        d=d+5
       end if 
      end do 

100   continue

      end 

c ****
