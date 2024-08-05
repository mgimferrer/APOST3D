      PROGRAM gen_hirs

      IMPLICIT DOUBLE PRECISION (a-h,o-z)

      INTEGER, PARAMETER :: maxat=20
      INTEGER, PARAMETER :: maxbasisf=200
      INTEGER, PARAMETER :: chargestates=5
      INTEGER, PARAMETER :: maxshells=50
      INTEGER, PARAMETER :: maxprimicoeff=20
      INTEGER, PARAMETER :: NMAX=100
      

      INTEGER :: status
      CHARACTER (len=80) :: line      
      CHARACTER (len=50) :: fileout   
      CHARACTER (len=80), DIMENSION(10)  :: gaussinp  
      CHARACTER (len=80), DIMENSION(chargestates)  :: atomsinp  
      CHARACTER (len=30) :: filename
      CHARACTER (len=30) :: keyword
      CHARACTER (len=20) :: gautitle   ! letters atoom en multipliciteit
      CHARACTER (len=8) :: lading
      CHARACTER (len=80) :: glocal,title
      CHARACTER (len=2), DIMENSION (maxat) :: atoms
     
      INTEGER, DIMENSION (maxat,chargestates) :: charge
      INTEGER, DIMENSION (maxat,chargestates) :: multipl

      DIMENSION  xlin((maxbasisf*(maxbasisf+1))/2)
      DIMENSION  Ca(maxbasisf,maxbasisf)
      DOUBLE PRECISION, DIMENSION (maxbasisf,maxbasisf) :: densmatrix
      INTEGER, DIMENSION (maxshells) :: shelltypes
      INTEGER, DIMENSION (maxshells) :: nprimi
      DOUBLE PRECISION, DIMENSION (maxshells,maxprimicoeff) ::
     * primexp
      DOUBLE PRECISION, DIMENSION (maxshells,maxprimicoeff) ::
     * contracoeff
      DOUBLE PRECISION, DIMENSION (maxshells,maxprimicoeff) ::
     * contracoeffSP

      DOUBLE PRECISION, DIMENSION (NMAX):: y2
      
      DOUBLE PRECISION :: xabs, yabs, zabs
      
      DOUBLE PRECISION :: x 
      DOUBLE PRECISION, DIMENSION (NMAX) :: radial     
       common /quadrat/th(1000),ph(1000),w(1000),wr(500),Xr(500)
       common/pha2/pha,phb,rr00

c      common /density/nsh,nprimi,primexp(maxshells,maxpmaxshells,
c     *maxprimicoeffrimicoeff),
c     *contracoeff(maxshells,maxprimicoeff),
c     *contracoeffSP(maxshells,maxprimicoeff),densmatrix
c     *(maxbasisf,maxbasisf),nbasf,shelltypes(maxshells)

! getting filename as argument

      CALL getarg(1,filename)
      filename=trim(filename)//".inp"
!      WRITE (*,*) 'de filename is: ', filename


! file opening and getting info

      OPEN (UNIT=12, FILE=filename, STATUS='OLD', ACTION='READ', 
     *IOSTAT=status)

      call reading(12,"$Gaussian",nlines, gaussinp)
      

      call locate2(12,"$Atoms")
      read(12,*) natoms
      do i=1,natoms
       read(12,*) atoms(i),(charge(i,j),multipl(i,j),j=1,chargestates)
      end do

      call locate2(12,"$Title")
      read(12,'(a80)') title 

!
! inlezen gegevens voor grid
      call locate2(12,"$Grid")
      read(12,*) nrad,nang

! get G03 location
      call locate2(12,"$GaussLocal")
      read(12,'(a80)') glocal

      idens=0
      irohf=0
      icalc=0
      igauss=0
      call locate2(12,"$Options")
      read(12,'(a80)') line  
      if(index(line,"NOCALC").ne.0) icalc=1
      if(index(line,"DENS2").ne.0) idens=1
      if(index(line,"ROHF").ne.0) irohf=1
      if(index(line,"g09").ne.0) igauss=1
      if(idens.eq.1.and.irohf.eq.1) then
       write(*,*) 'RO only for SCF density'
       idens=0
      end if

      CLOSE (UNIT=12)
      

c writing
!      write(*,*) '########'
!      do i=1,nlines,1
!       write(*,'(a80)') gaussinp(i)
!      end do
!      write(*,*) '########'
!
!
!      write(*,*) '########'
!      do i=2,nlines2,1
!
!       write(*,'(a80)') atomsinp(i)
!      end do

       nlines2=natoms+1


       if(icalc.eq.0) then

c write to gaussian com file

      OPEN (UNIT=99, FILE="runatoms.com")
      write(99,'(a11)') '#!/bin/bash'
      write(99,'(a43)') '# File generated automatically by gen_hirsh'
      write(99,*) 
      write(99,'(a80)') glocal  
      write(99,*) 


      DO i=1, nlines2-1, 1

       DO k=1, chargestates, 1

      IF (trim(atoms(i))/="H".OR.charge(i,k).lt.1) then


      SELECT CASE (charge(i,k))
      CASE (-1)
      lading='min'
      CASE (0)
      lading='nul'
      CASE (1)
      lading='plus'
      CASE (-2)
      lading= 'mintwee'
      CASE (+2)
      lading= 'plustwee'
      END SELECT


      gautitle= trim(atoms(i))//'_'//lading
      l = 13+i+k

      OPEN (UNIT=l, FILE="gauss_"//trim(gautitle)//".com")
      
      WRITE (l,'(a13)') '%chk=atom.chk'
      DO j=1, nlines, 1 
      WRITE (l,'(a80)') gaussinp (j)
      END DO

      WRITE (l, *) ''

      WRITE (l, *) gautitle,': this title was automatically generated'
      WRITE (l, *) ''       

      WRITE (l,*) charge(i,k), multipl(i,k)
      WRITE (l,*) atoms(i)
      WRITE (l,*) ''
      
      
      CLOSE (UNIT=l)

! aanmaak van de file die gaussian aanroept 
      write(99,*) "echo Running for ","gauss_"//trim(gautitle)//".com"
      if(igauss.eq.1) then
       write(99,*) "g09 ","gauss_"//trim(gautitle)//".com"
      else
       write(99,*) "g16 ","gauss_"//trim(gautitle)//".com"
      end if
      write(99,*) "formchk atom.chk ","gauss_"//trim(gautitle)//".fchk"
      
      END IF

      END DO
      END DO
      CLOSE (UNIT=99)

! aanroepen van gaussian
      call system( 'chmod +x runatoms.com')
      call system( './runatoms.com' )

      else
      write(*,*) 'NOCALC has been set. Atomic FChk files read from disk'

      end if
 
! bepalen vd grid
      rr00=0.50d0
      call quad(nrad, nang)
!omdraaien volgorde x waarden om spline en splint te doen werken
!voorlopige opslag in xlin
 
      n=1
      DO ik=0, nrad-1, 1
      xlin(n)=xr(nrad-ik)
      n=n+1
      END DO
 
      DO n=1, nrad,1
      xr(n)=xlin(n)
      END DO
      n=1
      DO ik=0, nrad-1, 1
      xlin(n)=wr(nrad-ik)
      n=n+1
      END DO
 
      DO n=1, nrad,1
      wr(n)=xlin(n)
      END DO
      
      
!      write(*,*) 'RAdial points: ',nrad
!      write(*,*) (xr(irad),irad=1,nrad)
!      write(*,*) 'Theta  points: ',nang
!      write(*,*) (th(irad),irad=1,nang)
!      write(*,*) 'Phi    points: ',nang
!      write(*,*) (ph(irad),irad=1,nang)

      OPEN (UNIT=98, FILE="densoutput")

! using different .fchk files

      WRITE (98,'(a80)') title           
! Total different atoms and radial points 
      WRITE (98,*) nlines2-1,nrad
!The gridpoints are 
      WRITE (98,*) (xr(irad),irad=1,nrad)

      DO i=1, nlines2-1, 1

! atomic symbol and number of states of the atom
      WRITE (98,'(a2)') atoms(i)
      WRITE (98,*)chargestates

       DO k=1, chargestates, 1

       IF ((trim(atoms(i)).EQ."H".AND.charge(i,k).GE.1).OR.
     & (trim(atoms(i)).EQ."He".AND.charge(i,k).GE.2)) THEN   
!einde ve lange if om Hplus te voorzien
       WRITE (98,*) charge(i,k), multipl(i,k)
       WRITE (98,*) (0.0d0,irad=1,nrad)

      ELSE


       SELECT CASE (charge(i,k))
       CASE (-1)
       lading='min'
       CASE (0)
       lading='nul'
       CASE (1)
       lading='plus'
       CASE (-2)
       lading= 'mintwee'
       CASE (+2)
       lading= 'plustwee'
       END SELECT

! building .fchk name

       gautitle= trim(atoms(i))//'_'//lading
       l = 43+i+k
       fileout="gauss_"//trim(gautitle)//".fchk"
 
! reading .fchk file

       OPEN (UNIT=l, FILE=fileout,IOSTAT=log, ACTION='READ')
       if(log/=0) then
        write(*,*) 'A problem occurred processing file ',fileout
        stop
       end if 

      CALL sameline("Number of basis functions",l,57,nbasf)

      rewind(l)
      read(l,'(a80)') line
      read(l,'(a80)') line
      if(irohf.eq.1) then
       CALL sameline("Number of alpha",l,57,nalf)
       CALL sameline("Number of beta",l,57,nbet )
       CALL sameline("Alpha MO coeff",l,57,ndim)  
       if(ndim.ne.nbasf*nbasf) stop 'inconsistency problem'
       READ(l,'(5e16.8)') ((ca(ii,jj),ii=1,nbasf),jj=1,nbasf)
! build P matrix from orbitals
       do ii=1,nbasf
        do jj=1,nbasf
         densmatrix(ii,jj)=0.0d0
         do kk=1,nbet
          densmatrix(ii,jj)=densmatrix(ii,jj)+2.0d0*ca(ii,kk)*ca(jj,kk)
         end do
         do kk=nbet+1,nalf
          densmatrix(ii,jj)=densmatrix(ii,jj)+ca(ii,kk)*ca(jj,kk)
         end do
        end do
       end do
       write(*,*) 'Atomic ROHF WaveFunctions. Building P from Orbitals' 
       write(*,*) 'Number of Alpha and Beta electrons: ',nalf,nbet  

      else 


      CALL sameline("Total SCF Density",l,57,ndim)  !Opletten met 
      !volgorde voor read
      if(idens.eq.1) then
       icont=0
       rewind(l)
13     read(l,'(a80)') line
       if(index(line,'Total ').ne.0) icont=icont+1
       if(icont.ne.3) go to 13
      end if
 
      READ(l,*) (xlin(io),io=1,ndim) !readen met do lus
                            !deze matrix w telkens overschreven
      
!      write(*,*) 'ndim: ',ndim
!      write (*,*) fileout
!      write (*,*) 'basisfunctis: ',nbasf

      m=1
      n=1
      io=1

      DO m=1, nbasf, 1
        DO n=1, m, 1   ! halve matrix inlezen, variabel #keren
        densmatrix(m,n)= xlin(io)   
        io=io+1
        END DO
      END DO


!symmetrisch maken v matrix
      DO m=1, nbasf, 1
        DO n=1, m, 1   
        densmatrix(n,m)=densmatrix(m,n)
        END DO
      END DO

      END IF


!lezen van de basisfuncties en primitieven


      CALL sameline("Shell types",l,57,nsh)        
      write (*,*) 'shelltypes: ',nsh
      READ (l,*) shelltypes(1:nsh)
      write(*,*) shelltypes(1:nsh)




      CALL sameline("Number of primitives per shell",l,57,npr)  
      READ (l,*) nprimi(1:npr)
       

      CALL sameline("Primitive exponents",l,57,nprco)
      write (*,*) 'Prim exp: ',nprco
      READ (l,*) xlin(1:nprco)
      write(*,*) xlin(1:nprco)

      n=1
      nn=1
      ii=1
      DO nn=1, npr,1 
       DO n=1, nprimi(nn), 1
       primexp(nn,n)=xlin(ii) !deze matrix is dezelfde voor atomen
                              !ongeacht de lading, k heeft gn belang
       
       ii=ii+1
       
       END DO
      END DO

    
!zelfde procedure voor contractiecoefficienten
      CALL sameline("Contraction coefficients",l,57,ncoco)
      READ (l,*) xlin(1:ncoco)
      write (*,*) 'Ncoco: ',ncoco
      write(*,*) xlin(1:ncoco)

      n=1
      nn=1
      ii=1
      DO nn=1, npr,1 
       DO n=1, nprimi(nn), 1
       contracoeff(nn,n)=xlin(ii) !deze matrix is dezelfde voor atomen
                                  !ongeacht de lading, k heeft gn belang

       ii=ii+1
       END DO
      END DO


! H atomen hebben geen SP orb
!        IF (trim(atoms(i))/="H") THEN


! contractiecoeff voor SP shells
       CALL sameline("P(S=P) Contraction coefficients",l,57,ncoco)
      if(ncoco.ne.0) then
      write (*,*) 'SP co coeff: ',ncoco
       READ (l,*) xlin(1:ncoco)
       write(*,*) xlin(1:ncoco)

       n=1
       nn=1 
       ii=1
      DO nn=1, npr , 1
       DO n=1, nprimi(nn), 1

       contracoeffSP(nn,n)=xlin(ii) !deze matrix is dezelfde voor atomen
       ii=ii+1
       END DO
      END DO

      END IF

      CLOSE (UNIT=l)




! densiteit rond een atoom berekenen
! werken met radiale en angulaire punten
! voor bep rad afstand, rho uitmiddelen over angulaire punten

! writen ie vooraf geopende file
      WRITE (98,*) charge (i,k) , multipl (i,k)

      xpedro=0.0d0

! via een dubbele loop de punten uit de common halen
! id commen geraakt door sub quad

      DO irad=1, nrad
! herinitialisatie bij versch radiale punten
       angav=0.0d0
       angmax=-1.0d0
       angmin= 10000.0d0
       DO iang=1, nang
        rr=xr(irad)
        th1=th(iang)
        ph1=ph(iang)

! omzetting van sferische coordinaten naar carthesiaanse
        xabs=rr*sin(th1)*cos(ph1) 
        yabs=rr*sin(th1)*sin(ph1) 
        zabs=rr*cos(th1) 

! uitrekenen dens op een punt en opstapelen in angav
      x = rho1(xabs,yabs,zabs,nsh,nprimi,primexp,contracoeff,
     *contracoeffSP,densmatrix,nbasf,shelltypes)
        IF (x>angmax) then
        angmax=x
        END IF
        IF (x<angmin) then
        angmin=x
        END IF
        angav=angav+x     

       xpedro=xpedro+x*wr(irad)*w(iang)*rr**2
       END DO

! uitmiddelen van de geintegreerd dens voor alle ang punten 
! met bepaalde radiale afstand
       radial(irad)=angav/nang
       angmax=max(angmax-radial(irad),-angmin+radial(irad))
      WRITE (*,'(a20,4f8.3)') 'Dist rho Maxdev: ', xr(irad),
     * radial(irad),angmax 
      IF (angmax>1.0) then
      WRITE (*,*) 'ANGULAR SYMMETRY IS LOW'
      END IF

      END DO


 

!density on radial points
      WRITE (98,*) (radial(irad),irad=1,nrad)

! checken of integratie over alle ang en rad punten, effectief 
! totaal aantal elektronen geeft
! 4pi r**2 dr

       write(*,*) 'Integration :',atoms(i),charge(i,k),multipl(i,k)
       write(*,*) 'Integration :',xpedro*4.0d0*dacos(-1.0d0)



      END IF

        END DO
      END DO


    

      call spline(xr, radial, nrad, -10.0d0, 0.0d0, y2)

  

      OPEN (UNIT=55, FILE="fitt.dat") 

      x1=0.10d0
          
      DO n=1, 10, 1
        call ssplint(xr, radial, y2, nrad , x1 ,y1)
        WRITE (55,'(a20,4f8.3)') 'De interpol geeft ',x1, y1
        x1=x1+0.1d0
      END DO



      END 









      subroutine reading(iunit,keyword,iline,gauinp)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      CHARACTER (len=80) :: line     
      CHARACTER*80, DIMENSION(10)  :: gauinp  
      CHARACTER*(*)  keyword 

      rewind(iunit)
c looking for gaussian
      iline=0
1     READ (iunit,'(a80)') line
      if(index(line,trim(keyword))/= 0) then
2      READ (iunit,'(a80)') line
       if(index(line,"$end")== 0) then
        iline=iline+1
        gauinp(iline)=line
        go to 2
       end if
      else
       go to 1
      end if


      end



      subroutine locate2(iunit,keyword)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      CHARACTER (len=80) :: line     
      CHARACTER*(*)  keyword 

      rewind(iunit)
c looking for gaussian
      iline=0
1     READ (iunit,'(a80)') line
      if(index(line,trim(keyword))== 0) then
       go to 1
      end if


      end



      subroutine sameline(keyword,iunit,iplaats,nd)
      
      CHARACTER*(*)  keyword 
      INTEGER :: outgetal
      INTEGER :: lg
      CHARACTER (len=80) :: line

      nd=0
      rewind(iunit)
      lg=1
      do while (lg==1)
       READ(iunit,'(a80)',err=10,end=10) line
       if(index(line,keyword)/=0) lg=0
      end do
      read(line(iplaats:),'(i5)') nd  !zoeken ie string, geg op zlfd lijn
10    continue

      end
 
   

      function rho1(xabs,yabs,zabs,nsh,nprimi,primexp,contracoeff,
     *contracoeffSP,densmatrix,
     *nbasf,shelltypes)

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      INTEGER, PARAMETER :: maxbasisf=200
      INTEGER, PARAMETER :: maxshells=50
      INTEGER, PARAMETER :: maxprimicoeff=20

       INTEGER :: expx, expy, expz
       DOUBLE PRECISION :: xabs, yabs, zabs
       DOUBLE PRECISION :: alfa  !voordurend overschreven contractiecoeff
       DOUBLE PRECISION :: rho
       DOUBLE PRECISION, DIMENSION (maxbasisf) :: chi
       INTEGER, DIMENSION (maxshells) :: nprimi 
      DOUBLE PRECISION, DIMENSION (maxshells,maxprimicoeff) :: 
     * primexp
      DOUBLE PRECISION, DIMENSION (maxshells,maxprimicoeff) :: 
     * contracoeff
      DOUBLE PRECISION, DIMENSION (maxshells,maxprimicoeff) ::  
     * contracoeffSP
      INTEGER, DIMENSION (maxshells) :: shelltypes
      DOUBLE PRECISION, DIMENSION (maxbasisf,
     * maxbasisf) :: densmatrix
      DOUBLE PRECISION :: x

         iorb=1           !absolute orbitaalnummerteller

       chi=0.0d0
       rho=0.0d0

        DO nn=1,nsh,1
         DO n=1, nprimi(nn), 1

          
         ishelltypes= shelltypes (nn)
C PSS
           alfa = primexp (nn,n)
          SELECT CASE (ishelltypes)
            
c           alfa = primexp (nn,n)

           CASE (0)
            expx=0
            expy=0
            expz=0
            chi(iorb)= (contracoeff(nn,n)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *  sqrt(primovlp(0,0,0,alfa,alfa)) +chi(iorb)
            
            IF (n == nprimi(nn)) then
!          WRITE (*,*) 'x', iorb,  chi(iorb)
            iorb = iorb + 1
            ELSE 
            END IF


           CASE (-1)
            expx=0
            expy=0
            expz=0
            chi(iorb+0)=(contracoeff(nn,n)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,0,0,alfa,alfa)) +chi(iorb+0)
            expx=1
            expy=0
            expz=0
       !opletten, contracoeff in andere matrix
            chi(iorb+1)=(contracoeffSP(nn,n)*(xabs**expx)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(1,0,0,alfa,alfa))+chi(iorb+1)
            expx=0
            expy=1
            expz=0
            chi(iorb+2)=(contracoeffSP(nn,n)*(yabs**expy)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,1,0,alfa,alfa))+chi(iorb+2)
            expx=0
            expy=0
            expz=1
            chi(iorb+3)=(contracoeffSP(nn,n)*(zabs**expz)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,0,1,alfa,alfa))+chi(iorb+3)

            IF (n == nprimi(nn)) then
c        WRITE (*,*) 'x', iorb,  chi(iorb)
c        WRITE (*,*) 'x', iorb+1,  chi(iorb+1)
c        WRITE (*,*) 'x', iorb+2,  chi(iorb+2)
c        WRITE (*,*) 'x', iorb+3,  chi(iorb+3)
            iorb = iorb + 4
            ELSE
            END IF

           CASE (1)
            expx=1
            expy=0
            expz=0
            chi(iorb)=(contracoeff(nn,n)*(xabs**expx)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(1,0,0,alfa,alfa))+chi(iorb)
            expx=0
            expy=1
            expz=0
            chi(iorb+1)=(contracoeff(nn,n)*(yabs**expy)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,1,0,alfa,alfa))+chi(iorb+1)
            expx=0
            expy=0
            expz=1
            chi(iorb+2)=(contracoeff(nn,n)*(zabs**expz)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,0,1,alfa,alfa))+chi(iorb+2)
            IF (n == nprimi(nn)) then
c         WRITE (*,*) 'x', iorb,  chi(iorb)
c         WRITE (*,*) 'x', iorb+1,  chi(iorb+1)
c         WRITE (*,*) 'x', iorb+2,  chi(iorb+2)
            iorb=iorb+3
            ELSE
            END IF



!opletten met orbvolgorde, moet overeenstmmn met densmatr
           CASE (2)
            expx=2
            expy=0
            expz=0
            chi(iorb+0)=(contracoeff(nn,n)*(xabs**expx)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(2,0,0,alfa,alfa))+chi(iorb)
            expx=0
            expy=2
            expz=0
            chi(iorb+1)=(contracoeff(nn,n)*(yabs**expy)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,2,0,alfa,alfa))+chi(iorb+1)
            expx=0
            expy=0
            expz=2
            chi(iorb+2)= (contracoeff(nn,n)*(zabs**expz)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,0,2,alfa,alfa))+chi(iorb+2)
            expx=1
            expy=1
            expz=0
            chi(iorb+3)=(contra
     *coeff(nn,n)*(xabs**expx)*(yabs**expy)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(1,1,0,alfa,alfa))+chi(iorb+3)
            expx=1
            expy=0
            expz=1
            chi(iorb+4)=(contra
     *coeff(nn,n)*(xabs**expx)*(zabs**expz)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(1,0,1,alfa,alfa))+chi(iorb+4)
            expx=0
            expy=1
            expz=1
            chi(iorb+5)=(contra
     *coeff(nn,n)*(yabs**expy)*(zabs**expz)*
     *exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,1,1,alfa,alfa))+chi(iorb+5)

            IF (n == nprimi(nn)) then
c         WRITE (*,*) 'x', iorb,  chi(iorb)
c         WRITE (*,*) 'x', iorb+1,  chi(iorb+1)
c         WRITE (*,*) 'x', iorb+2,  chi(iorb+2)
c         WRITE (*,*) 'x', iorb+3,  chi(iorb+3)
c         WRITE (*,*) 'x', iorb+4,  chi(iorb+4)
c         WRITE (*,*) 'x', iorb+5,  chi(iorb+5)
            iorb=iorb+6
            ELSE 
            END IF

           CASE (3)
            expx=3
            expy=0
            expz=0
            chi(iorb+0)=(contracoeff(nn,n)*(xabs**expx)*
     1exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     1sqrt(primovlp(3,0,0,alfa,alfa))+chi(iorb)
            expx=0
            expy=3
            expz=0
            chi(iorb+1)=(contracoeff(nn,n)*(yabs**expy)*
     1exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,3,0,alfa,alfa))+chi(iorb+1)
            expx=0
            expy=0
            expz=3
            chi(iorb+2)=(contracoeff(nn,n)*
     1(zabs**expz)*exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,0,3,alfa,alfa))+chi(iorb+2)
            expx=1
            expy=2
            expz=0
            chi(iorb+3)=(contracoeff(nn,n)*(xabs**expx)*(yabs**expy)*
     1exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(1,2,0,alfa,alfa))+chi(iorb+3)
            expx=2
            expy=1
            expz=0
            chi(iorb+4)=(contracoeff(nn,n)*(xabs**expx)*(yabs**expy)*
     1exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(2,1,0,alfa,alfa))+chi(iorb+4)
            expx=2
            expy=0
            expz=1
            chi(iorb+5)=(contracoeff(nn,n)*(xabs**expx)*
     1(zabs**expz)*exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(2,0,1,alfa,alfa))+chi(iorb+5)
            expx=1
            expy=0
            expz=2
            chi(iorb+6)=(contracoeff(nn,n)*(xabs**expx)*
     1(zabs**expz)*exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(1,0,2,alfa,alfa))+chi(iorb+6)
            expx=0
            expy=1
            expz=2
            chi(iorb+7)=(contracoeff(nn,n)*(yabs**expy)*
     1(zabs**expz)*exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,1,2,alfa,alfa))+chi(iorb+7)
            expx=0
            expy=2
            expz=1
            chi(iorb+8)=(contracoeff(nn,n)*(yabs**expy)*
     1(zabs**expz)*exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(0,2,1,alfa,alfa))+chi(iorb+8)
            expx=1
            expy=1
            expz=1
            chi(iorb+9)=(contracoeff(nn,n)*(xabs**expx)*(yabs**expy)*
     1(zabs**expz)*exp(-alfa*(xabs**2+yabs**2+zabs**2)))/
     *sqrt(primovlp(1,1,1,alfa,alfa))+chi(iorb+9)

            IF (n == nprimi(nn)) then
c         WRITE (*,*) 'x', iorb,  chi(iorb)
c         WRITE (*,*) 'x', iorb+1,  chi(iorb+1)
c         WRITE (*,*) 'x', iorb+2,  chi(iorb+2)
c         WRITE (*,*) 'x', iorb+3,  chi(iorb+3)
c         WRITE (*,*) 'x', iorb+4,  chi(iorb+4)
c         WRITE (*,*) 'x', iorb+5,  chi(iorb+5)
c         WRITE (*,*) 'x', iorb+6,  chi(iorb+6)
c         WRITE (*,*) 'x', iorb+7,  chi(iorb+7)
c         WRITE (*,*) 'x', iorb+8,  chi(iorb+8)
c         WRITE (*,*) 'x', iorb+9,  chi(iorb+9)
            iorb=iorb+10
            ELSE 
            END IF


          END SELECT


         END DO
        END DO
     
        rho1=0.0d0
        DO mu=1, nbasf, 1
          DO nu=1, nbasf, 1

          rho1 = densmatrix(mu,nu)*chi(mu)*chi(nu)+rho1

          END DO
        END DO

c       WRITE (*,*) 'rho1', rho1

      END




      function gamma(getal)
      implicit double precision(a-h,o-z)

      DOUBLE PRECISION :: pi
      DOUBLE PRECISION :: getal
      DOUBLE PRECISION :: x

      pi=dacos(-1.0d0)

      IF (getal==0.5) then
      x=SQRT(pi)
      END IF

      IF (getal==1.5) then
      x=SQRT(pi)/2
      END IF

      IF (getal==2.5) then
      x=(3*SQRT(pi))/4
      END IF

      IF (getal==3.5) then
      x=(15*SQRT(pi))/8
      END IF

      gamma=x
c       WRITE (*,*) 'gamma',x

      end



      function fac(getal)

      INTEGER :: getal
      DOUBLE PRECISION :: x=1.0d0
      DOUBLE PRECISION :: fac
      i=0

      DO WHILE (getal-i>=2) !maal 1 geen zin
       x = x*(getal-i)
       i=i+1
      END DO

      fac=x
       WRITE(*,*) 'fac', fac

      end

      function dfac(getal)

      INTEGER :: getal
      DOUBLE PRECISION :: u=1.0d0
      DOUBLE PRECISION :: dfac
      i=0

      DO WHILE (getal-i>=2)
       u = u*(getal-i)
       i=i+2
      END DO

      dfac=u
c       WRITE(*,*) 'dfac', dfac

      end

 
      function primovlp(xexp,yexp,zexp,alfaa,alfab)
      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      
 
      INTEGER :: lambda 
      DOUBLE PRECISION :: alfaunic
      INTEGER :: xexp, yexp, zexp
      DOUBLE PRECISION :: NO
      INTEGER :: a, b
      DOUBLE PRECISION :: xdoub, ydoub, zdoub
      DOUBLE PRECISION :: primovlp             

      pi=dacos(-1.0d0) 

      NO = 0.0d0

      xdoub=xexp+0.5
      ydoub=yexp+0.5
      zdoub=zexp+0.5

      overlap= (alfaa+alfab)**(-1.0-xexp-yexp-zexp)*
     *gamma(xdoub)*gamma(ydoub)*gamma(zdoub)
      overlap= overlap / ( SQRT(alfaa+alfab))
      
      primovlp = overlap
     
      end 

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)

      INTEGER n,NMAX
      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=100)
      INTEGER i,k
      DOUBLE PRECISION p,qn,sig,un,u(NMAX)

      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
      end do

      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      END DO

      END

      SUBROUTINE ssplint(xa,ya,y2a,n,x,y)

      INTEGER n
      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      DOUBLE PRECISION a,b,h
      klo=1
      khi=n

1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif

      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return

      END

