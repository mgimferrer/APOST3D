      subroutine top_3d(iorb,itype,ifunc,iatpairs)
      use ao_matrices
      implicit double precision(a-h,o-z)
      include 'parameter.h'
      logical ilog
      character*80 namecube
      character*60 name0
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /atlist/iatlist(maxat),icuat
      common /coord/ coord(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common /filename/name0

      dimension :: iatpairs(2,maxat)
      dimension :: igrid(3),xgrid(3,3),xpt(3)

      allocatable :: pgrid(:,:),pcoord(:,:),wp(:),omp2(:,:),omp(:)
      allocatable :: xao(:),xmo(:),xmo2(:,:),gxchp(:),gychp(:),gzchp(:)
      allocatable :: ftot(:,:,:),fab(:,:,:),faa(:,:,:)

      allocatable :: dm1aa(:,:),dm1bb(:,:),dm2(:,:,:,:),cumul(:,:,:,:)
      allocatable :: dm2c(:,:),eigenmat(:,:) 
      allocatable :: f1(:),f2(:,:),ff1(:),ff2(:,:)

!! OPTION ifunc CONTROLS IF WE INTRODUCE THE OPERATOR OR NOT !!

      ipairs   = Iopt(85)
      ietop    = Iopt(86)
      iorca    = Iopt(43)
      ipyscf   = Iopt(87)
      ispinsep = Iopt(88)
      if(ipairs.eq.0) iallmo=1

      iorb2=iorb*(iorb+1)/2
      ALLOCATE(dm2c(iorb2,iorb2))

!! itype = 2 FOR PHF EXCHANGE OR CORRELATION TOPOLOGY !!

      if(itype.eq.2) then 

!! READING DM1 !!

        ALLOCATE(dm1aa(iorb,iorb),dm1bb(iorb,iorb))
        thres=1.0d-6
        if(ispinsep.eq.0) then
          rewind(11)
          do while(.true.)
            read(11,end=997) i0,j0,vvv
            if(abs(vvv).gt.thres) then
              if(mod(i0,2).eq.0) then
                ii=i0/2
                jj=j0/2
                dm1bb(ii,jj)=dm1bb(ii,jj)+vvv
                dm1bb(jj,ii)=dm1bb(ii,jj)
              end if
              if(mod(i0,2).eq.1) then
                ii=(i0+1)/2
                jj=(j0+1)/2
                dm1aa(ii,jj)=dm1aa(ii,jj)+vvv
                dm1aa(jj,ii)=dm1aa(ii,jj)
              end if
            end if
          end do
        else
          rewind(30)
          do while(.true.)
            read(30,end=995) ii,jj,vvv
            if(abs(vvv).gt.thres) then
              dm1aa(ii,jj)=dm1aa(ii,jj)+vvv
              dm1aa(jj,ii)=dm1aa(ii,jj)
            end if
          end do
995       rewind(31)
          do while(.true.)
            read(31,end=997) ii,jj,vvv
            if(abs(vvv).gt.thres) then
              dm1bb(ii,jj)=dm1bb(ii,jj)+vvv
              dm1bb(jj,ii)=dm1bb(ii,jj)
            end if
          end do
        end if
997     continue

!! READING THE DM2(1,2,1,2) MATRIX BUT SAVING AS DM2(1,1,2,2), IMPORTANT FOR LATER ON !!

        ALLOCATE(dm2(iorb,iorb,iorb,iorb))
        icont=0
        rewind(12)
        do while(.true.)
          read(12,end=999) i0,k0,j0,l0,vvv
          if(abs(vvv).gt.thres) then
            ii=(i0-1)/2+1
            kk=(k0-1)/2+1
            jj=(j0-1)/2+1
            ll=(l0-1)/2+1
            if((ii.le.0.or.jj.le.0.or.kk.le.0.or.ll.le.0)) then
              write(*,*) " WARNING: IGNORED DM2 ELEMENT ",ii,kk,jj,ll,vvv
              write(*,*) " PROBLEM WITH THE CODE "
              stop
            else
              ilog=.true.
              if((i0.eq.j0.and.k0.eq.l0).or.(k0.eq.j0.and.i0.eq.l0)) ilog=.false.

!! IJ0 FOR THE AAAA, ABAB,BABA AND BBBB TERMS !!
!! KJ0 FOR THE ABBA AND BAAB TERMS !!

              ij0=mod(i0+j0,2)+mod(k0+l0,2)
              kj0=mod(k0+j0,2)+mod(i0+l0,2)
              if(ij0.eq.0)then
                icont=icont+1
                dm2(ii,jj,kk,ll)=dm2(ii,jj,kk,ll)+vvv
                if(ilog) dm2(jj,ii,ll,kk)=dm2(jj,ii,ll,kk)+vvv
                if(k0.ne.i0.and.l0.ne.j0) then
                  dm2(kk,ll,ii,jj)=dm2(kk,ll,ii,jj)+vvv
                  if(ilog) dm2(ll,kk,jj,ii)=dm2(ll,kk,jj,ii)+vvv
                end if
              end if
              if(kj0.eq.0) then
                icont=icont+1
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
          end if
        end do
999     continue

!! COMPUTING CUMULANT TENSOR (EXCHANGE FOR ietop = 1, CORRELATION FOR ietop = 2 AND XC FOR ietop = 3) !!

        ALLOCATE(cumul(iorb,iorb,iorb,iorb))
        do ii=1,iorb
          do jj=1,iorb
            do kk=1,iorb
              do ll=1,iorb
                if(ietop.eq.1) then
                  xx2=(dm1aa(ii,kk)+dm1bb(ii,kk))*(dm1aa(jj,ll)+dm1bb(jj,ll))/TWO
                  xx3=(dm1aa(ii,kk)-dm1bb(ii,kk))*(dm1aa(jj,ll)-dm1bb(jj,ll))/TWO 
                  cumul(ii,jj,kk,ll)=-xx2-xx3
                else if(ietop.eq.2) then
                  xx1=(dm1aa(ii,jj)+dm1bb(ii,jj))*(dm1aa(kk,ll)+dm1bb(kk,ll))
                  xx2=(dm1aa(ii,kk)+dm1bb(ii,kk))*(dm1aa(jj,ll)+dm1bb(jj,ll))/TWO
                  xx3=(dm1aa(ii,kk)-dm1bb(ii,kk))*(dm1aa(jj,ll)-dm1bb(jj,ll))/TWO
                  cumul(ii,jj,kk,ll)=dm2(ii,jj,kk,ll)-xx1+xx2+xx3
                else if(ietop.eq.3) then
                  xx1=(dm1aa(ii,jj)+dm1bb(ii,jj))*(dm1aa(kk,ll)+dm1bb(kk,ll))
                  cumul(ii,jj,kk,ll)=dm2(ii,jj,kk,ll)-xx1
                end if 
              end do
            end do
          end do
        end do
        DEALLOCATE(dm1aa,dm1bb,dm2)

!! COMPACTING INDEXES (I,J --> IJ AND K,L --> KL) !!

        ij=1
        do ii=1,iorb
          do jj=ii,iorb
            kl=1
            do kk=1,iorb
              do ll=kk,iorb
                if(jj.ne.ii) cumul(ii,jj,kk,ll)=cumul(ii,jj,kk,ll)+cumul(jj,ii,kk,ll)
                if(ll.ne.kk) cumul(ii,jj,kk,ll)=cumul(ii,jj,kk,ll)+cumul(ii,jj,ll,kk)
                if(jj.ne.ii.and.ll.ne.kk) cumul(ii,jj,kk,ll)=cumul(ii,jj,kk,ll)+cumul(jj,ii,ll,kk)
                dm2c(ij,kl)=cumul(ii,jj,kk,ll)
                kl=kl+1
              end do
            end do
            ij=ij+1
          end do
        end do
        DEALLOCATE(cumul)

!! DIAGONALIZING THE CORRELATION COMPACTED TENSOR !!

        ALLOCATE(eigenmat(iorb2,iorb2))
        call diagonalize(iorb2,iorb2,dm2c,eigenmat,0)

!! RHF CASE !!

      else 
        ij=0
        do ii=1,iorb
          do jj=1,iorb
            ij=ij+1
            if(ii.eq.jj) dm2c(ij,ij)=FOUR 
            if(ii.ne.jj) dm2c(ij,ij)=ZERO 
          end do 
        end do  

!! END IF WAVEFUNCTION TYPE !!

      end if 

!! SETTING GRID FOR NUMERICAL INTEGRATION !!
!! FIRST GRID WILL BE CUBE-LIKE AND ASSUMING RECTANGULAR GRID !! 

      npt1=0
      ngridp=200
      do ii=1,3
        igrid(ii)=ngridp
        npt1=npt1+igrid(ii)
        do jj=1,3
          xgrid(ii,jj)=ZERO
        end do
      end do
      ixx=igrid(1)*igrid(2)*igrid(3)
      write(*,*) " Total number of grid points for CUBE : ",ixx
      write(*,*) " "

!! SECOND ONE (AS USED TO INTEGRATE) WILL BE SPHERICAL !! 

      nrad=20
      nang=74
      pha=ZERO
      phb=ZERO
      rr00=0.1d0
      itotps=nrad*nang
      ALLOCATE(pcoord(itotps,3),wp(itotps))
      call quad_iso(nrad,nang,rr00,pha,phb,pcoord,wp)

!! ALLOCATING EVERYTHING HERE !!
!! gxchp, gychp AND gzchp NOT USED BUT NECESSARY !!

      ALLOCATE(pgrid(maxgrid,3),xao(igr),xmo(iorb),xmo2(itotps,iorb))
      ALLOCATE(gxchp(igr),gychp(igr),gzchp(igr),omp2(itotps,nat),omp(nat))
      ALLOCATE(ftot(igrid(1),igrid(2),igrid(3)),faa(igrid(1),igrid(2),igrid(3)))
      ALLOCATE(fab(igrid(1),igrid(2),igrid(3)))

!! ENTIRE MOLECULE CUBE CASE !!

      if(iallmo.eq.1) then

!! ADDING EXTRA SIZE TO THE CAGE !!

        extra=THREE
        do ii=1,3
          xmax=-1.0d8
          xmin=1.0d8
          do iatom=1,nat
            if(coord(ii,iatom).lt.xmin) xmin=coord(ii,iatom)
            if(coord(ii,iatom).gt.xmax) xmax=coord(ii,iatom)
          end do
          xmin=xmin-extra
          xmax=xmax+extra
          dist=xmax-xmin
          xgrid(ii,ii)=dist/(real(igrid(ii))-ONE)
          do jj=1,igrid(ii)
           pgrid(jj,ii)=xmin+(jj-1)*xgrid(ii,ii)
          end do
        end do

!! TOPOLOGY CALCULATION !!

        irun=0
        diff=ZERO
        do i1=1,igrid(1)
          do j1=1,igrid(2)
            do k1=1,igrid(3)
              irun=irun+1
              xabs=pgrid(i1,1)
              yabs=pgrid(j1,2)
              zabs=pgrid(k1,3)

!! EVALUATING AOs, CONVERSION TO MOs AND ATOMIC WEIGHT COMPUTING FOR 1st GRID !!

              call gpoints(xabs,yabs,zabs,gxchp,gychp,gzchp,xao)
              do kk=1,iorb
                xx=ZERO
                do ll=1,igr
                  xx=xx+c(ll,kk)*xao(ll)
                end do
                xmo(kk)=xx
              end do
              do iat=1,nat
                omp(iat)=wat(iat,xabs,yabs,zabs)
              end do

!! NOW FOR THE SECOND ELECTRON GRID POINTS !!

              do jj=1,itotps
                xabs2=pcoord(jj,1)+xabs
                yabs2=pcoord(jj,2)+yabs
                zabs2=pcoord(jj,3)+zabs
                call gpoints(xabs2,yabs2,zabs2,gxchp,gychp,gzchp,xao)
                do kk=1,iorb
                  xx=ZERO
                  do ll=1,igr
                    xx=xx+c(ll,kk)*xao(ll)
                  end do
                  xmo2(jj,kk)=xx
                end do
                do iat=1,nat
                  omp2(jj,iat)=wat(iat,xabs2,yabs2,zabs2)
                end do
              end do

!! FINISHING INDEX COMPACTATION, PRODUCT OF ORBITALS FROM FIRST GRID IN f1, FROM SECOND IN f2 !!

              ALLOCATE(f1(iorb2),f2(itotps,iorb2),ff1(iorb2),ff2(itotps,iorb2))
              ij=0
              do kk=1,iorb
                do jj=kk,iorb
                  ij=ij+1
                  f1(ij)=xmo(kk)*xmo(jj)
                end do
              end do
              do ifut=1,itotps
                ij=0
                do kk=1,iorb
                  do jj=kk,iorb
                    ij=ij+1
                    f2(ifut,ij)=xmo2(ifut,kk)*xmo2(ifut,jj)
                  end do
                end do
              end do

!! TRANSFORMATION OF f1 AND f2 FROM DIAGONALIZATION OF COMPACTED CUMULANT !!

              if(itype.eq.2) then 
                do jj=1,iorb2
                  xx=ZERO
                  do kk=1,iorb2
                    xx=xx+eigenmat(kk,jj)*f1(kk)
                  end do
                  ff1(jj)=xx
                end do
                do ifut=1,itotps
                  do jj=1,iorb2
                    xxb=ZERO
                    do kk=1,iorb2
                      xxb=xxb+eigenmat(kk,jj)*f2(ifut,kk)
                    end do
                    ff2(ifut,jj)=xxb
                  end do
                end do
              else 

!! RHF CASE !!

                do ifut=1,itotps
                  do jj=1,iorb2
                    if(ifut.eq.1) ff1(jj)=f1(jj)
                    ff2(ifut,jj)=f2(ifut,jj)
                  end do
                end do
              end if 
              DEALLOCATE(f1,f2)

!! INTEGRATING THE CUMMULANT !!

              xtot=ZERO
              do jj=1,itotps
                dist=dsqrt(pcoord(jj,1)**TWO+pcoord(jj,2)**TWO+pcoord(jj,3)**TWO)
                if(ifunc.eq.1) dist=ONE
                if(dist.gt.1.0d-12) then
                  do ij=1,iorb2
                    xtot=xtot+dm2c(ij,ij)*wp(jj)*ff1(ij)*ff2(jj,ij)/dist
                  end do
                end if 
              end do 
              ftot(i1,j1,k1)=xtot

!! ONLY ATOMIC AND INTERATOMIC CONTRIBUTIONS ARE FEASIBLE (aa = ATOMIC AND ab = INTERATOMIC) !!

              xaa=ZERO
              xab=ZERO
              xba=ZERO
              do jj=1,itotps
                do iat=1,nat
                  do jat=1,nat
!                 do jat=iat,nat
                    if(iat.eq.jat) xwaa=wp(jj)*omp(iat)*omp2(jj,iat)
                    if(iat.ne.jat) then
                      xwab=wp(jj)*omp(iat)*omp2(jj,jat)
!                     xwba=wp(jj)*omp(jat)*omp2(jj,iat)
                    end if 
                    dist=dsqrt(pcoord(jj,1)**TWO+pcoord(jj,2)**TWO+pcoord(jj,3)**TWO)
                    if(ifunc.eq.1) dist=ONE
                    if(dist.gt.1.0d-12) then
                      do ij=1,iorb2
                        if(iat.eq.jat) xaa=xaa+dm2c(ij,ij)*xwaa*ff1(ij)*ff2(jj,ij)/dist
                        if(iat.ne.jat) then
                          xab=xab+dm2c(ij,ij)*xwab*ff1(ij)*ff2(jj,ij)/dist
!                         xba=xba+dm2c(ij,ij)*xwba*ff1(ij)*ff2(jj,ij)/dist
                        end if 
                      end do
                    end if 
                  end do 
                end do 
              end do 
              faa(i1,j1,k1)=faa(i1,j1,k1)+xaa
              fab(i1,j1,k1)=fab(i1,j1,k1)+xab+xba
              diff=diff+xtot-(xaa+xab+xba)

!! END DOs OF GRID !!

              DEALLOCATE(ff1,ff2)
            end do
          end do
          write(*,*) " Topology Process (%) : ",REAL(irun)*100.0d0/REAL(ixx)
        end do
        write(*,*) " Difference between Total and Atomic + Diatomic : ",diff


!! PRINTING ftot !!

        namecube=trim(name0)//"-ftot.cube"
        open(44,file=namecube,status="unknown")
        rewind(44)
        write(44,*)'Cube generated with APOST-3D code '
        if(iallmo.eq.1) then
          if(ietop.eq.1) write(44,*) " Molecular Exchange energy "
          if(ietop.eq.2) write(44,*) " Molecular Correlation energy "
        else 
!! (TO DO)
          if(ietop.eq.1) write(44,*) " Exchange energy atom pair : "
          if(ietop.eq.2) write(44,*) " Correlation energy atom pair : "
        end if 
        write(44,42) -nat,(pgrid(1,j),j=1,3)
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
            write(44,40)(ftot(i,j,k),k=1,igrid(3))
          end do
        end do
        close(44)

!! PRINTING faa !!

        namecube=trim(name0)//"-faa.cube"
        open(44,file=namecube,status="unknown")
        rewind(44)
        write(44,*)'Cube generated with APOST-3D code '
        if(iallmo.eq.1) then
          if(ietop.eq.1) write(44,*) " Molecular Exchange energy "
          if(ietop.eq.2) write(44,*) " Molecular Correlation energy "
        else 
!! (TO DO)
          if(ietop.eq.1) write(44,*) " Exchange energy atom pair : "
          if(ietop.eq.2) write(44,*) " Correlation energy atom pair : "
        end if 
        write(44,42) -nat,(pgrid(1,j),j=1,3)
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
            write(44,40)(faa(i,j,k),k=1,igrid(3))
          end do
        end do
        close(44)

!! PRINTING fab !!

        namecube=trim(name0)//"-fab.cube"
        open(44,file=namecube,status="unknown")
        rewind(44)
        write(44,*)'Cube generated with APOST-3D code '
        if(iallmo.eq.1) then
          if(ietop.eq.1) write(44,*) " Molecular Exchange energy "
          if(ietop.eq.2) write(44,*) " Molecular Correlation energy "
        else 
!! (TO DO)
          if(ietop.eq.1) write(44,*) " Exchange energy atom pair : "
          if(ietop.eq.2) write(44,*) " Correlation energy atom pair : "
        end if 
        write(44,42) -nat,(pgrid(1,j),j=1,3)
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
            write(44,40)(fab(i,j,k),k=1,igrid(3))
          end do
        end do
        close(44)
     
        DEALLOCATE(ftot,fab,faa,eigenmat,dm2c)


!**************
!! CUBE BETWEEN PAIR OF ATOMS CASE !!

      else 
!       do iiatp=1,ipairs

!! ADDING EXTRA SIZE TO THE CAGE !!

!         extra=TWO
!         do ii=1,3
!           xmax=-1.0d8
!           xmin=1.0d8
!           do jatom=1,2
!             iatom=iatpairs(jatom,iiatp)
!             if(coord(ii,iatom).lt.xmin) xmin=coord(ii,iatom)
!             if(coord(ii,iatom).gt.xmax) xmax=coord(ii,iatom)
!           end do
!           xmin=xmin-extra
!           xmax=xmax+extra
!           dist=xmax-xmin
!           xgrid(ii,ii)=dist/(real(igrid(ii))-ONE)
!           do jj=1,igrid(ii)
!            pgrid(jj,ii)=xmin+(jj-1)*xgrid(ii,ii)
!           end do
!         end do
!
!! MAKING ZERO THE VECTORS !!

!         do ii=1,npt1
!           totcorr(ii)=ZERO
!           aacorr(ii)=ZERO
!           abcorr(ii)=ZERO
!           bacorr(ii)=ZERO
!           bbcorr(ii)=ZERO
!         end do

!         icenter=iatpairs(1,iiatp)
!         jcenter=iatpairs(2,iiatp)
!         do ii=1,igrid(1)
!           do jj=1,igrid(2)
!             do kk=1,igrid(3)
!
!               xabs=pgrid(ii,1)
!               yabs=pgrid(jj,2)
!               zabs=pgrid(kk,3)
!
!
!
!
!             end do
!           end do
!         end do

!       end do 
      end if 

!! DEALLOCATING !!
      DEALLOCATE(pgrid,xao,xmo,xmo2,gxchp,gychp,gzchp,omp2,omp)
      DEALLOCATE(pcoord,wp)

!! PRINTING FORMATS !!
40    FORMAT(6e13.5)
41    FORMAT(a8,x,a8,a7,i3,a9,a2,a11,f8.4,a9,f8.4)
42    FORMAT(i5,3f12.6)
43    FORMAT(i5,4f12.6)
44    FORMAT(a8,x,a8,a7,i3,a9,i2,a11,f8.4,a9,f8.4)

      end

!*****

      subroutine quad_iso(nrad,nang,rr00,pha,phb,pcoord,wp)
      implicit double precision(a-h,o-z)
      include 'parameter.h'

      dimension :: wp(nrad*nang),pcoord(nrad*nang,3)

      allocatable :: wr(:),xr(:)
      allocatable :: X(:),Y(:),Z(:),W(:)
      allocatable :: th1(:),ph1(:)

      ALLOCATE(wr(nrad),xr(nrad))
      call LEGZO(nrad,xr,wr)

!! RADIAL PART, GAUSS-LEGENDRE QUADRATURE MODIFIED RMAX TRUNCATING PARAMETER !!
!! IF APHA = -1, THEN NO TRUNCATION !!

!      alpha=(2.0d0*rr00/rmax)-1.0d0
      alpha=-ONE
      do i=1,nrad
        wr(i)=rr00*wr(i)*(ONE-alpha)/((ONE+alpha*xr(i))**TWO)
        xr(i)=rr00*((ONE+xr(i))/(ONE+(alpha*xr(i))))
      end do

!! ANGULAR PART, LEBEDEV-LAIKOV GRIDS !!

      ALLOCATE(X(nang),Y(nang),Z(nang),W(nang))
      if(nang.eq.6)     CALL LD0006(X,Y,Z,W,N)
      if(nang.eq.14)    CALL LD0014(X,Y,Z,W,N)
      if(nang.eq.26)    CALL LD0026(X,Y,Z,W,N)
      if(nang.eq.38)    CALL LD0038(X,Y,Z,W,N)
      if(nang.eq.50)    CALL LD0050(X,Y,Z,W,N)
      if(nang.eq.74)    CALL LD0074(X,Y,Z,W,N)
      if(nang.eq.86)    CALL LD0086(X,Y,Z,W,N)
      if(nang.eq.110)   CALL LD0110(X,Y,Z,W,N)
      if(nang.eq.146)   CALL LD0146(X,Y,Z,W,N)
      if(nang.eq.170)   CALL LD0170(X,Y,Z,W,N)
      if(nang.eq.194)   CALL LD0194(X,Y,Z,W,N)
      if(nang.eq.230)   CALL LD0230(X,Y,Z,W,N)
      if(nang.eq.266)   CALL LD0266(X,Y,Z,W,N)
      if(nang.eq.302)   CALL LD0302(X,Y,Z,W,N)
      if(nang.eq.350)   CALL LD0350(X,Y,Z,W,N)
      if(nang.eq.434)   CALL LD0434(X,Y,Z,W,N)
      if(nang.eq.590)   CALL LD0590(X,Y,Z,W,N)
      if(nang.eq.770)   CALL LD0770(X,Y,Z,W,N)
      if(nang.eq.974)   CALL LD0974(X,Y,Z,W,N)

!! TO SPHERICAL COORDINATES !!

      ALLOCATE(th1(nang),ph1(nang))
      do iang=1,nang
        xx=X(iang)
        yy=Y(iang)
        zz=Z(iang)
        th1(iang)=dacos(zz)
        if(dsin(th1(iang)).ne.ZERO) then
          xs=xx/dsin(th1(iang))
          if(xs.gt.ONE) xs=ONE
          if(xs.lt.-ONE) xs=-ONE
          ph1(iang)=dacos(xs)
          if(yy.lt.ZERO) ph1(iang)=-ph1(iang)
        else
          ph1(iang)=ZERO
        end if
      end do

!! BACK TO CARTESIAN COORDINATES WITH ROTATION INCLUDED !!

      irun=1
      do k=1,nrad
        do i=1,nang
          thx=th1(i)
          fix=ph1(i)
          rr=xr(k)
          xabs0=rr*dsin(thx)*dcos(fix)
          yabs0=rr*dsin(thx)*dsin(fix)
          zabs0=rr*dcos(thx)
          xabs=dcos(phb)*xabs0+dsin(phb)*yabs0
          yabs=-dcos(pha)*dsin(phb)*xabs0
          yabs=yabs+dcos(pha)*dcos(phb)*yabs0+dsin(pha)*zabs0
          zabs=dsin(pha)*dsin(phb)*xabs0
          zabs=zabs-dsin(pha)*dcos(phb)*yabs0+dcos(pha)*zabs0

!! SAVING COORDINATES INTO VECTOR PCOORD !!

          pcoord(irun,1)=xabs
          pcoord(irun,2)=yabs
          pcoord(irun,3)=zabs
          irun=irun+1
        end do
      end do

!! FINISHING TO CREATE THE INTEGRATION WEIGHTS !!

      irun=1
      do k=1,nrad
        xxr=wr(k)*xr(k)*xr(k)
        do i=1,nang
          wp(irun)=W(i)*xxr*FOUR*pi
          irun=irun+1
        end do
      end do

      DEALLOCATE(wr,xr)
      DEALLOCATE(X,Y,Z,W)
      DEALLOCATE(th1,ph1)

      end
      
!*****

