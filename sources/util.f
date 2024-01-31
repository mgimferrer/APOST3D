      subroutine kiir
      write(*,'(a80)') '------------------------------------------------------------------------------'
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '                    Program APOST-3D, Version 4                               '
      write(*,'(a80)') '                             04-03-2021                                       '
      write(*,'(a80)') '                   --------------------------------                           '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Real-space and Hilbert-space tools for wave function analysis             '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Available atomic definitions:                                             '
      write(*,'(a80)') '    ----------------------------                                              '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Real space:                                                               '
      write(*,'(a80)') '      Becke, J. Chem. Phys. 88 2547 1988                                      '
      write(*,'(a80)') '      Hirshfeld, Theor. Chim. Acta 44  129 1977                               '
      write(*,'(a80)') '      Hirshfeld-Iterative, J Chem Phys 126 144111 2007                        '
      write(*,'(a80)') '      Topological fuzzy Voronoi cells (TFVC), J Chem Phys, 139 071103 2013    '
      write(*,'(a80)') '      QTAIM, J. Comput. Chem 30 1082 2009                                     '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Hilbert-space : Mulliken, Lowdin, Davidson-Lowdin                         '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Calculating:                                                              '
      write(*,'(a80)') '    ------------                                                              '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      A)  Atomic and overlap populations, bond orders and valences            '
      write(*,'(a80)') '         I. Mayer and P. Salvador, Chem. Phys. Lett. 383 368-375 2004	      '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      B) Hartree-Fock molecular energy decomposition                          '
      write(*,'(a80)') '         P. Salvador, M. Duran, I.Mayer, J. Chem. Phys. 115 1153-1157 2001    '
      write(*,'(a80)') '         P. Salvador and I. Mayer, J. Chem. Phys. 120 5046-5052 2004          '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      C) KS-DFT molecular energy decomposition                                '
      write(*,'(a80)') '         P. Salvador, I. Mayer, J. Chem. Phys. 126 234113 2007	              '
      write(*,'(a80)') '         M. Gimferrer, P. Salvador, J. Chem. Phys. 158 234105 2023            '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      D) Molecular energy decomposition for CAS/DMRG wavefunctions            '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      E) Effective atomic orbitals:                                           '
      write(*,'(a80)') '         I. Mayer, J. Phys. Chem. 100 6249 1996                               '
      write(*,'(a80)') '         I. Mayer and P. Salvador, J. Chem. Phys. 130 234106 2009             '
      write(*,'(a80)') '         E. Ramos-Cordoba et.al., J. Chem. Phys. 138 214107 2013              '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      F) Local spin analysis                                                  '
      write(*,'(a80)') '         E. Ramos-Cordoba et.al., J. Chem. Theor. Comput. 8, 1270-1279 2012   '
      write(*,'(a80)') '         E. Ramos-Cordoba et.al., Phys. Chem. Chem. Phys. 14 15291-15298 2012 '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      G) Effective Oxidation states analysis                                  '
      write(*,'(a80)') '         E. Ramos-Cordoba et.al., J. Chem. Theor. Comput. 11 1501-1508 2015   '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      H) Oxidation states from localized orbitals                             '
      write(*,'(a80)') '         M. Gimferrer et al., J. Chem. Theor. Comput. 18 309-322 2022         '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '      I) Decomposition of EDA quantities into one- and two-center IQA terms   '
      write(*,'(a80)') '         M. Gimferrer et al., J. Chem. Theory Comput. 19 3469-3485 2023       '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '    Cite this program as:                                                     '
      write(*,'(a80)') '    ---------------------                                                     '
      write(*,'(a80)') '            P. Salvador, E. Ramos-Cordoba, M. Gimferrer, M. Montilla          '
      write(*,'(a80)') '            Program APOST-3D, Version 4, Girona, 2020                       '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '  e-mail: psalse@gmail.com, eloy.raco@gmail.com, mgimferrer18@gmail.com       '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '------------------------------------------------------------------------------'
      write(*,'(a80)') '  The program has been written by using parts of the program APOST by         '
      write(*,'(a80)') '  I. Mayer and A. Hamza, Budapest, 2000-2003.                                 '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '  The numerical integration utilizes the subroutines for Lebedev              '
      write(*,'(a80)') '  quadrature downloaded from CCL. The appropriate reference is:               '
      write(*,'(a80)') '  V.I. Lebedev, and D.N. Laikov "A quadrature formula for the sphere of the   '
      write(*,'(a80)') '  131st algebraic order of accuracy" Doklady Mathematics, 59 477-481 1999.    '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '  The program makes use of libxc library when necessary, using the F90        '
      write(*,'(a80)') '  interfaces provided by the authors.                                         '
      write(*,'(a80)') '                 (see http://www.tddft.org/programs/libxc)                    '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') '  We are extremely grateful for the possibility of using these routines!      '
      write(*,'(a80)') '                                                                              '
      write(*,'(a80)') ' -----------------------------------------------------------------------------'
      end

! *****

      subroutine gennatural()
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension dm1(igr,igr)

      allocatable :: s0(:,:),sm(:,:),splus(:,:)

      nbasis=igr

      allocate(s0(igr,igr),sm(igr,igr),splus(igr,igr))

c make S1/2
       s0=s
       call build_Smp(igr,s0,sm,splus,0)

c tranform P with S1/2
       do i=1,igr
        do j=1,igr
        occ_no(i,j)=p(i,j)
        end do
       end do
      call to_lowdin_basis(igr,splus,occ_no)
      call diagonalize(igr,igr,occ_no,c_no,0)
      call to_AO_basis(igr,igr,sm,c_no)

      i=1
      do while(abs(occ_no(i,i)).gt.1.0d-6) 
       i=i+1
      end do
      write(*,*) 'Natural orbital occupations (>1.0e-6): '
      write(*,'(8f10.5)') (occ_no(ii,ii),ii=1,i-1)

       write(*,*) 'Checking NO overlap matrix '
! reusing s0
        do i=1,igr
         do j=1,igr
          xx=0.0d0
          do k=1,igr
           xx=xx+c_no(k,i)*s(k,j)
          end do
          s0(i,j)=xx
         end do
        end do
        do i=1,igr
         do j=1,igr
          xx=0.0d0
          do k=1,igr
           xx=xx+s0(i,k)*c_no(k,j)
          end do
         if(i.eq.j.and.abs(xx-1.0d0).gt.1.0d-5) write(*,*) i,j,xx
         if(i.ne.j.and.abs(xx).gt.1.0d-5) write(*,*) i,j,xx
         end do
        end do

      deallocate(s0,sm,splus)
      return
      end

      subroutine gennatural_old()
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      dimension dm1(igr,igr)

      allocatable :: scr1(:,:),scr2(:,:)

      nbasis=igr

!! COMPUTING DM1, scr1 DIMENSIONATED TO nmax TO MATCH occ_no DIMENSIONS !!

      ALLOCATE(scr2(igr,igr),scr1(nmax,nmax))

!! Ct*S !!

      x1=ZERO
       do i=1,nbasis
        do nu=1,igr
         x=0.0d0
         x1=x1+p(i,nu)*s(nu,i)
          do mu=1,igr
           x=x+c(mu,i)*s(mu,nu)
          end do
          scr1(i,nu)=x
        end do
       end do
c P*SC=P*(Ct*S)t
       do i=1,nbasis
        do nu=1,igr
         x=0.0d0
         do mu=1,igr
          x=x+scr1(i,mu)*p(mu,nu)
         end do
         scr2(i,nu)=x
        end do
       end do
c (Ct*S)*(P*SC)
       do i=1,nbasis
        do j=1,nbasis
         x=0.0d0
         do nu=1,igr
          x=x+scr2(i,nu)*scr1(j,nu)
         end do
         dm1(i,j)=x
        end do
       end do

      do i=1,nbasis
       do k=1,nbasis
        occ_no(i,k)=dm1(i,k)
       end do
      end do
C diagonalize to get NOs
      call diagonalize(nbasis,nbasis,occ_no,scr1,0)

      i=1
      do while(abs(occ_no(i,i)).gt.1.0d-6) 
       i=i+1
      end do
      write(*,*) 'Natural orbital occupations (>1.0e-6): '
      write(*,'(8f10.5)') (occ_no(ii,ii),ii=1,i-1)
C trasnform to the AO basis
      do i=1,igr
       do j=1,nbasis
        x=0.0d0
        do k=1,nbasis
         x=x+c(i,k)*scr1(k,j)
        end do
        c_no(i,j)=x
       end do
      end do

       write(*,*) 'Checking NO overlap matrix '
        do i=1,igr
         do j=1,igr
          xx=0.0d0
          do k=1,igr
           xx=xx+c_no(k,i)*s(k,j)
          end do
          scr1(i,j)=xx
         end do
        end do
        do i=1,igr
         do j=1,igr
          xx=0.0d0
          do k=1,igr
           xx=xx+scr1(i,k)*c_no(k,j)
          end do
         if(i.eq.j.and.abs(xx-1.0d0).gt.1.0d-5) write(*,*) i,j,xx
         if(i.ne.j.and.abs(xx).gt.1.0d-5) write(*,*) i,j,xx
         end do
        end do


      DEALLOCATE(scr1,scr2)
      return
      end

! *****s
      subroutine readintfiles(sat)
      use basis_set
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      character*80 linia      
      dimension sat(nbasis,nbasis,natoms)

c loop over atoms
       call locate(16,"# INTFILES",ii)
       if(ii.eq.0)  stop 'Required section not found in input file'
        do ia=1,nat
         read(16,'(a80)') linia   
c         k=1
c         do while(index(linia(k:k)," ").ne.0)
c            k=k+1
c         end do
c         ka=k 
c         do while(index(linia(k:k)," ").eq.0)
c            k=k+1
c         end do
c         kb=k 
c         open(unit=22,file=linia(ka:kb))
         open(unit=22,file=trim(linia))
C ONLY FOR RESTRICTED
         do while (index(linia,"The Atomic Overlap Matrix").eq.0)
          read(22,'(a80)') linia
         end do 
         read(22,'(a80)') linia   
         read(22,'(a80)') linia   
         if(index(linia,"Restricted").eq.0) stop 'only resctricted'
         do i=1,nocc
           read(22,*) (sat(i,j,ia),j=1,i)
           do j=1,i
            sat(j,i,ia)=sat(i,j,ia)
           end do
         end do 
        end do 
       end
       
 
! *****s
        Subroutine build_Smp(n,s0,sm,sp,ip)
        implicit double precision(a-h,o-z)
        include 'parameter.h'
        integer, intent(in) :: n
        dimension s0(n,n),sm(n,n),sp(n,n)
        allocatable :: x(:,:),s(:,:)

        allocate(x(n,n),s(n,n))
 
C       Builds Sm(p) matrix from eigenvectors and eigenvalues of
C       overlap matrix S.IP controls whether to calculate Sp.
        call move_to(N,s0,s)
        call diagonalize(n,n,s,x,0)
 
        do i=1,n
         do j=i,n
           Sm(j,i)=0.0d0
           if (ip.ne.1) Sp(j,i)=0.0d0
           do k=1,n
            if(s(k,k).gt.1.0d-8) then
             xx=x(i,k)*x(j,k)
             xs=sqrt(s(k,k))
             Sm(j,i)=Sm(j,i)+xx/xs
             if(ip.ne.1) Sp(j,i)=Sp(j,i)+xx*xs
            end if
           end do
           Sm(i,j)=Sm(j,i)
           if (ip.ne.1) Sp(i,j)=Sp(j,i)
         end do
        end do
        deallocate(x,s)
        return
        end

*********************************************************************
        Subroutine to_lowdin_basis(N,X,F)
        implicit double precision(a-h,o-z)
        include 'parameter.h'
        integer,  intent(in) :: n
        dimension f(n,n), X(n,n)
        allocatable :: f0(:,:)
 
        allocate (f0(n,n))
C       Similarity transformation of Fock matrix into lowdin-basis
 
        do j=1,n
         do i=1,n
          xx=0.0d0
          do k=1,n
           xx=xx+f(i,k)*x(k,j)
          end do
           f0(i,j)=xx
         end do
        end do
        do j=1,n
          do i=1,n
           f(i,j)=f0(i,j)
          end do
        end do
        do j=1,n
         do i=1,n
          xx=0.0d0
          do k=1,n
           xx=xx+x(k,i)*f(k,j)
          end do
           f0(i,j)=xx
         end do
        end do
        do j=1,n
          do i=1,n
           f(i,j)=f0(i,j)
          end do
        end do
        deallocate (f0)
        return
        end

        Subroutine trans_rect_mat(N,M,X,F)
        implicit double precision(a-h,o-z)
        include 'parameter.h'
        common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
        dimension x(igr,igr),f(igr,igr)
        allocatable :: f0(:,:)

        allocate (f0(igr,igr))
c X(N,M) 
        do j=1,m
         do i=1,n
          xx=0.0d0
          do k=1,n
           xx=xx+f(i,k)*x(k,j)
          end do
           f0(i,j)=xx
         end do
        end do
        do j=1,m
          do i=1,n
           f(i,j)=f0(i,j)
          end do
        end do
        do j=1,m
         do i=1,m
          xx=0.0d0
          do k=1,n
           xx=xx+x(k,i)*f(k,j)
          end do
           f0(i,j)=xx
         end do
        end do
        do j=1,m
          do i=1,m
           f(i,j)=f0(i,j)
          end do
        end do
        deallocate(f0)
        return
        end


**********************************7*********************************
        Subroutine diagonalize(M,N,A0,X,ival)
        implicit double precision(a-h,o-z)
        include 'parameter.h'
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

C DIAGONALIZATION OF THE REAL SYMMETRIC MATRIX X. IN D THE EIGENVALUES.
      Subroutine SDIAG2(X,M,N,D,inosort)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
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
C*********************************************************************
        Subroutine move_to(n,A,B)
        implicit double precision(a-h,o-z)
        include 'parameter.h'
        integer, intent(in) :: n
        dimension A(n,n),B(n,n)
 
        do j=1,N
         do i=1,N
          b(i,j)=a(i,j)
         end do
        end do
        return
        end
**********************************************************************
        Subroutine to_AO_basis(n,na,X,C) 
        implicit double precision(a-h,o-z)
        include 'parameter.h'
        integer, intent(in) :: n
        dimension x(n,n),c(n,n)
        allocatable :: cx(:,:)

        allocate (cx(n,n))
         
C       Transforms first na  orbitals from lowdin basis to AO basis
 
        do j=1,na
         do i=1,n
           cx(i,j)=0.0d0
            do k=1,n
             cx(i,j)=cx(i,j)+x(i,k)*c(k,j)
            end do
         end do
        end do
 
        do j=1,na
         do i=1,n
          c(i,j)=cx(i,j)
         end do
        end do
 
        deallocate (cx)
        return
       end


        subroutine invert(n,a) 
        implicit double precision (a-h,o-z) 
        include 'parameter.h'
        integer, intent(in) :: n
        dimension a(n,n) 
        dimension x(n,n),b(n)  
        real*8 d,e,an,t,s,ta,xh,vxkp,vxkq,vpk,vkp,vkq,vqk 
        logical cont 
 
 
C       DEFINIM LA TOLERANCIA DEL METODE         
        tol=1.0d-12
 
C       INICIALITZEM LA MATRIU X COM A MATRIU IDENTITAT 
         
        do i=1,n 
         do j=1,n 
          x(i,j)=0.0d0 
         end do 
         x(i,i)=1.0d0 
        end do 
         
        ia=0.0 
 
C       ALGORITME PRINCIPAL         
 
        cont=.true. 
        do while (cont) 
         ia =ia+1 
         xmax=0.0 
          do i=1,n-1 
           do j=i+1,n 
            if (abs(a(i,j)).gt.xmax) then 
             xmax=abs(a(i,j)) 
             ip=i 
             iq=j 
            end if 
           end do 
          end do 
         
         if (xmax.gt.tol) then 
          d=a(iq,iq)-a(ip,ip) 
          if (d.ne.0.0)then 
           e=a(ip,iq)/d 
            if (abs(e).gt.tol) then 
             an=1/(2*e) 
             t=1/(abs(an)+sqrt(1+(an)**2)) 
             if (an.lt.0.0) t=-t 
            else 
             t=e 
            endif 
          else 
           t=1.0 
          end if 
           
          c=1/(sqrt(1+t**2)) 
          s=t*c 
          ta=s/(1+c) 
          xh=t*(a(ip,iq)) 
          
          a(ip,iq)=0.0 
          a(iq,ip)=0.0 
          a(ip,ip)=(a(ip,ip)-(xh)) 
          a(iq,iq)=(a(iq,iq)+(xh)) 
 
           do k=1,n 
            vxkq=(X(k,ip)-(ta*(x(k,iq)))) 
            vxkp=(X(k,iq)+(ta*(x(k,ip)))) 
            
            x(k,ip)=x(k,ip)-(s*(vxkp)) 
            x(k,iq)=x(k,iq)+(s*(vxkq)) 
 
             if (k.ne.ip.and.k.ne.iq) then 
            
             vpk=(a(iq,k)+(ta*(a(ip,k)))) 
             vqk=(a(ip,k)-(ta*(a(iq,k)))) 
             vkp=(a(k,iq)+(ta*(a(k,ip)))) 
             vkq=(a(k,ip)-(ta*(a(k,iq)))) 
 
             a(ip,k)=a(ip,k)-(s*(vpk)) 
             a(iq,k)=a(iq,k)+(s*(vqk)) 
             a(k,ip)=a(k,ip)-(s*(vkp)) 
             a(k,iq)=a(k,iq)+(s*(vkq)) 
            end if 
           end do 
        
         else  
          cont=.false. 
         end if 
         end do 
 
 
c       Obtencion de matriz inversa 
 
             do k=1,n 
               b(k)=1.0d0/a(k,k) 
             end do 
         
        do i=1,n 
           do j=1,n 
             a(i,j)=0.0d0 
             do k=1,n 
               a(i,j)=a(i,j) + x(i,k)*x(j,k)*b(k) 
             end do 
           end do 
         end do 
 
         return 
         end 
 
c************************************************************************ 
c*  SVD SUBROUTINES ***************************************************** 
c************************************************************************ 

      subroutine svd(m,n,a,w,matu,u,matv,v,ierr)
      integer i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      real*8 a(m,n),w(n),u(m,n),v(m,n),rv1(n)
      real*8 c,f,g,h,s,x,y,z,tst1,tst2,scale,pythag
      logical matu,matv
c
c     this subroutine is a translation of the algol procedure svd,
c     num. math. 14, 403-420(1970) by golub and reinsch.
c     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).
c
c     this subroutine determines the singular value decomposition
c          t
c     a=usv  of a real m by n rectangular matrix.  householder
c     bidiagonalization and a variant of the qr algorithm are used.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.  note that nm must be at least
c          as large as the maximum of m and n.
c
c        m is the number of rows of a (and u).
c
c        n is the number of columns of a (and u) and the order of v.
c
c        a contains the rectangular input matrix to be decomposed.
c
c        matu should be set to .true. if the u matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c        matv should be set to .true. if the v matrix in the
c          decomposition is desired, and to .false. otherwise.
c
c     on output
c
c        a is unaltered (unless overwritten by u or v).
c
c        w contains the n (non-negative) singular values of a (the
c          diagonal elements of s).  they are unordered.  if an
c          error exit is made, the singular values should be correct
c          for indices ierr+1,ierr+2,...,n.
c
c        u contains the matrix u (orthogonal column vectors) of the
c          decomposition if matu has been set to .true.  otherwise
c          u is used as a temporary array.  u may coincide with a.
c          if an error exit is made, the columns of u corresponding
c          to indices of correct singular values should be correct.
c
c        v contains the matrix v (orthogonal) of the decomposition if
c          matv has been set to .true.  otherwise v is not referenced.
c          v may also coincide with a if u is not needed.  if an error
c          exit is made, the columns of v corresponding to indices of
c          correct singular values should be correct.
c
c        ierr is set to
c          zero       for normal return,
c          k          if the k-th singular value has not been
c                     determined after 30 iterations.
c
c        rv1 is a temporary storage array.
c
c     calls pythag for  sqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      pythag(x,y)= dsqrt(x*x + y*y) 
      ierr = 0
c
      do 100 i = 1, m
c
         do 100 j = 1, n
            u(i,j) = a(i,j)
  100 continue
c     .......... householder reduction to bidiagonal form ..........
      g = 0.0d0
      scale = 0.0d0
      x = 0.0d0
c
      do 300 i = 1, n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m) go to 210
c
         do 120 k = i, m
  120    scale = scale +dabs(u(k,i))
c
         if (scale .eq. 0.0d0) go to 210
c
         do 130 k = i, m
            u(k,i) = u(k,i) / scale
            s = s + u(k,i)**2
  130    continue
c
         f = u(i,i)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,i) = f - g
         if (i .eq. n) go to 190
c
         do 150 j = l, n
            s = 0.0d0
c
            do 140 k = i, m
  140       s = s + u(k,i) * u(k,j)
c
            f = s / h
c
            do 150 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  150    continue
c
  190    do 200 k = i, m
  200    u(k,i) = scale * u(k,i)
c
  210    w(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i .gt. m .or. i .eq. n) go to 290
c
         do 220 k = l, n
  220    scale = scale +dabs(u(i,k))
c
         if (scale .eq. 0.0d0) go to 290
c
         do 230 k = l, n
            u(i,k) = u(i,k) / scale
            s = s + u(i,k)**2
  230    continue
c
         f = u(i,l)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,l) = f - g
c
         do 240 k = l, n
  240    rv1(k) = u(i,k) / h
c
         if (i .eq. m) go to 270
c
         do 260 j = l, m
            s = 0.0d0
c
            do 250 k = l, n
  250       s = s + u(j,k) * u(i,k)
c
            do 260 k = l, n
               u(j,k) = u(j,k) + s * rv1(k)
  260    continue
c
  270    do 280 k = l, n
  280    u(i,k) = scale * u(i,k)
c
  290    x = dmax1(x,dabs(w(i))+dabs(rv1(i)))
  300 continue
c     .......... accumulation of right-hand transformations ..........
      if (.not. matv) go to 410
c     .......... for i=n step -1 until 1 do -- ..........
      do 400 ii = 1, n
         i = n + 1 - ii
         if (i .eq. n) go to 390
         if (g .eq. 0.0d0) go to 360
c
         do 320 j = l, n
c     .......... double division avoids possible underflow ..........
  320    v(j,i) = (u(i,j) / u(i,l)) / g
c
         do 350 j = l, n
            s = 0.0d0
c
            do 340 k = l, n
  340       s = s + u(i,k) * v(k,j)
c
            do 350 k = l, n
               v(k,j) = v(k,j) + s * v(k,i)
  350    continue
c
  360    do 380 j = l, n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
  380    continue
c
  390    v(i,i) = 1.0d0
         g = rv1(i)
         l = i
  400 continue
c     .......... accumulation of left-hand transformations ..........
  410 if (.not. matu) go to 510
c     ..........for i=min(m,n) step -1 until 1 do -- ..........
      mn = n
      if (m .lt. n) mn = m
c
      do 500 ii = 1, mn
         i = mn + 1 - ii
         l = i + 1
         g = w(i)
         if (i .eq. n) go to 430
c
         do 420 j = l, n
  420    u(i,j) = 0.0d0
c
  430    if (g .eq. 0.0d0) go to 475
         if (i .eq. mn) go to 460
c
         do 450 j = l, n
            s = 0.0d0
c
            do 440 k = l, m
  440       s = s + u(k,i) * u(k,j)
c     .......... double division avoids possible underflow ..........
            f = (s / u(i,i)) / g
c
            do 450 k = i, m
               u(k,j) = u(k,j) + f * u(k,i)
  450    continue
c
  460    do 470 j = i, m
  470    u(j,i) = u(j,i) / g
c
         go to 490
c
  475    do 480 j = i, m
  480    u(j,i) = 0.0d0
c
  490    u(i,i) = u(i,i) + 1.0d0
  500 continue
c     .......... diagonalization of the bidiagonal form ..........
  510 tst1 = x
c     .......... for k=n step -1 until 1 do -- ..........
      do 700 kk = 1, n
         k1 = n - kk
         k = k1 + 1
         its = 0
c     .......... test for splitting.
c                for l=k step -1 until 1 do -- ..........
  520    do 530 ll = 1, k
            l1 = k - ll
            l = l1 + 1
            tst2 = tst1 +dabs(rv1(l))
            if (tst2 .eq. tst1) go to 565
c     .......... rv1(1) is always zero, so there is no exit
c                through the bottom of the loop ..........
            tst2 = tst1 +dabs(w(l1))
            if (tst2 .eq. tst1) go to 540
  530    continue
c     .......... cancellation of rv1(l) if l greater than 1 ..........
  540    c = 0.0d0
         s = 1.0d0
c
         do 560 i = l, k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            tst2 = tst1 +dabs(f)
            if (tst2 .eq. tst1) go to 565
            g = w(i)
            h = pythag(f,g)
            w(i) = h
            c = g / h
            s = -f / h
            if (.not. matu) go to 560
c
            do 550 j = 1, m
               y = u(j,l1)
               z = u(j,i)
               u(j,l1) = y * c + z * s
               u(j,i) = -y * s + z * c
  550       continue
c
  560    continue
c     .......... test for convergence ..........
  565    z = w(k)
         if (l .eq. k) go to 650
c     .......... shift from bottom 2 by 2 minor ..........
         if (its .eq. 30) go to 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = 0.5d0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
         g = pythag(f,1.0d0)
         f = x - (z / x) * z + (h / x) * (y / (f + dsign(g,f)) - h)
c     .......... next qr transformation ..........
         c = 1.0d0
         s = 1.0d0
c
         do 600 i1 = l, k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = pythag(f,h)
            rv1(i1) = z
            c = f / z
            s = h / z
            f = x * c + g * s
            g = -x * s + g * c
            h = y * s
            y = y * c
            if (.not. matv) go to 575
c
            do 570 j = 1, n
               x = v(j,i1)
               z = v(j,i)
               v(j,i1) = x * c + z * s
               v(j,i) = -x * s + z * c
  570       continue
c
  575       z = pythag(f,h)
            w(i1) = z
c     .......... rotation can be arbitrary if z is zero ..........
            if (z .eq. 0.0d0) go to 580
            c = f / z
            s = h / z
  580       f = c * g + s * y
            x = -s * g + c * y
            if (.not. matu) go to 600
c
            do 590 j = 1, m
               y = u(j,i1)
               z = u(j,i)
               u(j,i1) = y * c + z * s
               u(j,i) = -y * s + z * c
  590       continue
c
  600    continue
c
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
         go to 520
c     .......... convergence ..........
  650    if (z .ge. 0.0d0) go to 700
c     .......... w(k) is made non-negative ..........
         w(k) = -z
         if (.not. matv) go to 700
c
         do 690 j = 1, n
  690    v(j,k) = -v(j,k)
c
  700 continue
c
      go to 1001
c     .......... set error -- no convergence to a
c                singular value after 30 iterations ..........
 1000 ierr = k
 1001 return
      end
      subroutine uvprint(a,ndec,n,m)
      implicit real*8(a-h,o-z)
      dimension a(ndec,m)
      kmax=min0(8,n)
      write(*,62) (i,i=1,kmax)
  62  format(3x,8i9)
      do i=1,n
      write(*,61)i,(a(i,j),j=1,kmax)
      enddo
      print *,' '
  61  format(1x,i4,8f9.5)
   1  if(kmax.ge.m) return
      k1=kmax+1
      kmax=kmax+8
      kmax=min0(kmax,n)
      write(*,62) (i,i=k1,kmax)
      do i=1,n
      write(*,61)i,(a(i,j),j=k1,kmax)
      enddo
      print *,' '
      goto 1
      end
      subroutine inc(C,m,n)
c      subroutine inc(C,m,n,n1)
      implicit real*8(a-h,o-z)
      parameter (Large=1000)
      dimension C(Large,Large),ifal(2),ibe(20)
      dimension mol(2)
      data mol /4H Mol,4Hecul/
c      read(23,*)mm,nn,n1
c      if(mm.ne.m) then
c      print *,' Wrong m'
c      stop 11
c      endif
c      if(nn.ne.n) then
c      print *,' Wrong n'
c      stop 12
c     print *,'  Entering inc'
       print *,' M,N',m,n

      do i=1,m
      do j=1,m
      c(i,j)=0.d0
      enddo
      enddo
      rewind 15
 11    read(15,150,end=220,err=220)ibe
 150  format(20a4)
      if(ibe(2).eq.mol(1).and.ibe(3).eq.mol(2))go to 10
      go to 11
c 10  print *,' Found'
  10  read(15,61)ifal
      read(15,61)ifal
      read(15,61)ifal
 61   format(2a4)
c     print *,' Reading orbitals:'
      print *,' '

      kmax=min0(5,n)         
      do i=1,m
      read (15,60,end=220,err=220) ii,(c(i,j),j=1,kmax)
c     write (*,60) ii,(c(i,j),j=1,kmax)

  60  format(i4,17x,5f10.5)
      if(i.ne.ii) stop 111
      enddo
       print *,' '
  1   if(kmax.ge.n) return
      read(15,61,end=220,err=220)ifal
      read(15,61,end=220,err=220)ifal
      read(15,61,end=220,err=220)ifal
      k1=kmax+1
      kmax=kmax+5
      kmax=min0(kmax,n)
      do i=1,m
      read (15,60,end=220,err=220) ii,(c(i,j),j=k1,kmax)
c     write (*,60) ii,(c(i,j),j=k1,kmax)

c 62  format(1x,i4,5f9.5)
      if(i.ne.ii) stop 112
      enddo
c      print *,' '
      goto 1
  220  return
c     x=0.d0
c     do i=1,m
c     do j=1,m
c     x=x+c(i,j)**2
c     enddo
c     enddo
c     print *,'x',x
      end 
      subroutine outc(c,m,n,sigma,m0)
       implicit real*8(a-h,o-z)
      parameter (Large=1000)
      dimension C(Large,Large),sigma(Large)
      dimension iki(10)
      data icom/4H(A)-/
      data iocc/4H-O  /
      data ivirt/4H-V  /
      print *,' Transformed orbitals:'
      print *,' '
      print *,'       Occupied'
      print *,' '

      kmax=min0(5,n)
      write(*,61)(j,j=1,kmax)
      l=1
      do j=1,kmax
      iki(l)=icom
      iki(l+1)=iocc
      if(j.ge.m0)iki(l+1)=ivirt
      l=l+2
      enddo
      write(*,66)iki
  66  format(24x,5(2a4,2x))
c     print *,' '
      write(*,64)(sigma(j),j=1,kmax)
  64  format('   LAMBDA VALUES --  ',5f10.5)
      do i=1,m
      write (*,60) i,(c(i,j),j=1,kmax)
  60  format(i4,17x,5f10.5)
      enddo
  61  format(18x,5i10)
  1   if(kmax.ge.n) return
      print *,' '
      print *,'        Virtual'
      print *,' '
      k1=m0+1
      kmax=m0+5
      kmax=min0(kmax,n)
      write(*,61)(j,j=k1-m0,kmax-m0)
      l=1
      do j=k1,kmax
      iki(l)=icom
      iki(l+1)=iocc
      if(j.ge.m0)iki(l+1)=ivirt
      l=l+2
      enddo
      write(*,66)iki
      write(*,64)(sigma(j),j=k1,kmax)
c     print *,'    EIGENVALUES --     '
      do i=1,m
      write (*,60) i,(c(i,j),j=k1,kmax)
      enddo
c     goto 1
      return
      end
      subroutine outc23(c,m,n,sigma,m0)
       implicit real*8(a-h,o-z)
      parameter (Large=1000)
      dimension C(Large,Large),sigma(Large)
      dimension iki(10)
      data icom/4H(A)-/
      data iocc/4H-O  /
      data ivirt/4H-V  /
c     print *,' Entering outc23'

      kmax=min0(5,n)
      write(23,61)(j,j=1,kmax)
      l=1
      do j=1,kmax
      iki(l)=icom
      iki(l+1)=iocc
      if(j.ge.m0)iki(l+1)=ivirt
      l=l+2
      enddo
      write(23,66)iki
  66  format(24x,5(2a4,2x))
c     print *,' '
      write(23,65)(200*sigma(j)**2,j=1,kmax)
  64  format('     EIGENVALUES --  ',5f10.5)
  65  format('     EIGENVALUES --  ',5(f7.2,'   '))
      do i=1,m
      write (23,60) i,(c(i,j),j=1,kmax)
  60  format(i4,17x,5f10.5)
      enddo
  61  format(18x,5i10)
  1   if(kmax.ge.n) return
      k1=m0+1
      kmax=m0+5
      kmax=min0(kmax,n)
      write(23,61)(j,j=k1-m0,kmax-m0)
      l=1
      do j=k1,kmax
      iki(l)=icom
      iki(l+1)=iocc
      if(j.ge.m0)iki(l+1)=ivirt
      l=l+2
      enddo
      write(23,66)iki
      write(23,65)(200*sigma(j)**2,j=k1,kmax)
c     print *,'    EIGENVALUES --     '
      do i=1,m
      write (23,60) i,(c(i,j),j=k1,kmax)
      enddo
c     goto 1
      return
      end
      subroutine outc2(c,nbas,m2,sigma,m0)
      implicit real*8(a-h,o-z)
      parameter (Large=1000)
      dimension C(Large,Large),ifal(2),ibe(20),sigma(Large)
      dimension mol(2)
      data mol /4H Mol,4Hecul/
c      read(23,*)mm,nn,n1
c      if(mm.ne.m) then
c      print *,' Wrong m'
c      stop 11
c      endif
c      if(nn.ne.n) then
c      print *,' Wrong n'
c      stop 12
c      print *,' M, N, N1',m,n,n1
c     print *,'  Entering outc2'

      rewind 23
      rewind 15
c     print *,'  Entering outc2'

c     iii=1
 11    read(15,150,end=220)ibe
c     iii=iii+1
c     write(*,150)ibe

      write(23,150)ibe

 150  format(20a4)
      if(ibe(2).eq.mol(1).and.ibe(3).eq.mol(2)) then
      print *,' Calling OUTC23'
      call outc23(c,nbas,m2,sigma,m0) 
      endif
      go to 11
c 220 print *,' III=',iii 
  220 return
      end 

       subroutine solvesystem(n,ndim,A,C,ilog)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
       integer, intent(in) :: n,ndim
       DIMENSION A(ndim,ndim+1),C(ndim)
       logical ilog
 
       ilog=.true.

      do i=1,N
c pivoting
       do l=i,N
        do m=l,N
         if (abs(A(m,i)).gt.abs(A(l,i))) then
          do kk=1,N+1
           f=A(m,kk)
           A(m,kk)=A(l,kk)
           A(l,kk)=f
          end do
         end if
        end do
       end do
c  final pivotin
c compatibilidad del sistema
       if (Abs(A(i,i)).lt.1e-15)then
        write(*,*) 'sistema no comp. determinat'
        ilog=.false.
        return
       end if
       
       do k=1,N
        if (k.ne.i)then
         alfa=A(k,i)/A(i,i)
         do j=1,N+1
          A(k,j)=A(k,j)-alfa*A(i,j)
         end do
        end if
       end do
      end do
   
      do i=1,N
       P=A(i,i)
       do j=1,N+1
        A(i,j)=A(i,j)/P
       end do
      end do

      do i=1,n   
       c(i)= A(i,n+1)
      end do
      end
