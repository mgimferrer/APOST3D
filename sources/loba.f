      subroutine eos_loba(sat)

      use basis_set
      use ao_matrices
      use integration_grid
      implicit double precision(a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord2(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(100)
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      allocatable :: orbpop(:)
      allocatable :: flo(:,:),clo(:,:)

      dimension sat(igr,igr,nat)

      !! VALUES FOR THE CI ACCORDING TO OUR INORG. CHEM PAPER !!
      PPP=0.200d0
      WWW=0.100d0

      idofr=iopt(40)
      if(idofr.eq.0) then
        write(*,*) " FRAGMENT DEFINITION REQUIRED FOR LOBA, CHECK .inp FILE "
        stop
      end if

      !! FRAGMENT CHARGE EXTRACTED FROM zn (AVOIDS PROBLEMS WHEN PSEUDOPOTENTIALS ARE USED) !!
      iznrest=0
      do ifrg=1,icufr
        izn=0
        do icenter=1,nfrlist(ifrg)
          iiat=ifrlist(icenter,ifrg)
          izn=izn+INT(zn(iiat))
        end do

        !! ONLY 2 FRAGMENTS CONSIDERED FOR LOBA, FRG1 VS THE REST !!
        if(ifrg.eq.1) then
          iznfrg1=izn
        else
          iznrest=iznrest+izn
        end if
      end do

      !! STARTING LOBA !!
      if(kop.ne.1) then

        !! CALCULATING FIRST FRAGMENT POPULATIONS FOR EACH LO !!
        ALLOCATE(orbpop(nocc))
        ALLOCATE(flo(nocc,icufr),clo(igr,igr))

        flo=ZERO
        clo=c

        !! LOOP OVER FRAGMENTS !!
        do jfrg=1,icufr
          orbpop=ZERO
          call rwf_uwf_frg_pop(jfrg,sat,nocc,clo,orbpop)

          !! SAVING FRAGMENT POPULATIONS !!
          do ii=1,nocc
            flo(ii,jfrg)=orbpop(ii)
          end do
        end do

        write(*,*) " ----------------------------------------- "
        write(*,*) "  Fragment population analysis of each LO  "
        write(*,*) " ----------------------------------------- "
        write(*,*) " "
        call rwf_uwf_print_LO_pop(nocc,flo)

        !! ASSIGNING ELECTRONS AND CI VALUES !!
        write(*,*) " -------------------------------- "
        write(*,*) "  Assigning electrons of each LO  "
        write(*,*) " -------------------------------- "
        write(*,*) " "

        !! EVALUATING CLARITY INDEX (CI) FOR EACH LOCALIZED ORBITAL !!
        write(*,*) " INFO: Considering Frg. 1 for OS assignment "
        write(*,*) " "
        write(*,*) "  Orbital     CI     Frg. 1 Pop    Assignment   "
        write(*,*) " ---------------------------------------------- "

        iel1=0
        ielrest=0
        do ii=1,nocc
          xfrg1=ZERO
          xfrgrest=ZERO
          do jfrg=1,icufr
            if(jfrg.eq.1) xfrg1=flo(ii,jfrg)
            if(jfrg.gt.1) xfrgrest=xfrgrest+flo(ii,jfrg)
          end do

          !! EVALUATING x FIRST, AND THEN THE CI VALUE ASSOCIATED TO IT !!
          xx=ABS((xfrg1-xfrgrest))/(xfrg1+xfrgrest)

          !! CI EVALUATION IN THREE RANGES !!
          if(xx.ge.ZERO.and.xx.le.PPP-WWW) then !! x = [0, P-W] !!
            xCI=100.0d0
            write(*,110) ii,xCI,xfrg1,"Covalent"
            iel1=iel1+1
            ielrest=ielrest+1
          else if(xx.ge.PPP+WWW.and.xx.le.ONE) then !! x = [P+W, 1] !!
            xCI=100.0d0

            !! ELECTRONS TO FRAGMENT 1 !!
            if(xfrg1-xfrgrest.ge.ZERO) then
              write(*,111) ii,xCI,xfrg1,"Ionic(Frg.1)"
              iel1=iel1+2

            !! ELECTRONS TO THE REST !!
            else if(xfrg1-xfrgrest.lt.ZERO) then
              write(*,112) ii,xCI,xfrg1,"Ionic(Rest)"
              ielrest=ielrest+2
            end if

!            iio=1
          else !! x = (P-W, P+W) !!
            xCI=100.0d0*DCOS(DCOS(pi*(PPP+WWW-xx)/TWO*WWW))
            if(xx.le.PPP) then
              write(*,110) ii,xCI,xfrg1,"Covalent"
              iel1=iel1+1
              ielrest=ielrest+1  
            else if(xx.gt.PPP) then

              !! ELECTRONS TO FRAGMENT 1 !!
              if(xfrg1-xfrgrest.ge.ZERO) then
                write(*,111) ii,xCI,xfrg1,"Ionic(Frg.1)"
                iel1=iel1+2

              !! ELECTRONS TO THE REST !!
              else if(xfrg1-xfrgrest.lt.ZERO) then
                write(*,112) ii,xCI,xfrg1,"Ionic(Rest)"
                ielrest=ielrest+2
              end if
            end if
          end if
        end do

        !! CLOSING THE TABLE !!
        write(*,*) " ---------------------------------------------- "
        write(*,*) " "

        !! FINAL PRINTING !!
        write(*,*) " --------------------------- "
        write(*,*) "  FRAGMENT OXIDATION STATES  "
        write(*,*) " --------------------------- "
        write(*,*) " "
        write(*,*) "  Frag.  Oxidation State  "
        write(*,*) " ------------------------ "
        write(*,'(3x,a3,6x,f8.2)') "  1",REAL(iznfrg1-iel1)
        write(*,'(3x,a5,4x,f8.2)') " Rest",REAL(iznrest-ielrest) 
        write(*,*) " ------------------------ "
        write(*,*) " "

        !! DEALLOCATING MATRICES !!
        DEALLOCATE(orbpop,flo)

      !! FOR THE MOMENT ONLY CLOSED-SHELL !!
      else
        write(*,*) " LOBA IMPLEMENTED ONLY FOR CLOSED-SHELL "
        stop
      end if

      !! PRINTING FORMATS !!
110   FORMAT(5x,i3,3x,f7.2,5x,f7.3,5x,a8)
111   FORMAT(5x,i3,3x,f7.2,5x,f7.3,5x,a12)
112   FORMAT(5x,i3,3x,f7.2,5x,f7.3,5x,a11)

      end

!! ****** !!

      subroutine rwf_uwf_print_LO_pop(norb,frgpop)

!! ROUTINE FOR PRINTING THE FRAGMENT POPULATIONS OF EACH ORBITAL (INTENDED FOR LOCALIZED) !!

      implicit double precision (a-h,o-z)
      include 'parameter.h'

      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /frlist/ifrlist(maxat,maxfrag),nfrlist(maxfrag),icufr,jfrlist(maxat)

      dimension frgpop(norb,icufr)

!! TRICK OF THE INTEGER ROUNDING FOR NUMBER OF COLUMNS !!
      b=norb/5
      if(b*5.ne.norb) b=(norb/5)+1

!! PRINTING !!
      dd=1
      do k=1,b

!! FOR THE LAST PACK OF COLUMNS !!
        if(k.eq.b) then
          write(*,'(2x,a15,5(i7,3x))') "Orbital Number:",(jj,jj=dd,norb)
          do jfrg=1,icufr
            write(*,'(2x,a9,i3,a1,5f10.5)') "Frg. Pop.",jfrg,":",(frgpop(jj,jfrg),jj=dd,norb)
          end do
          write(*,*) " "

!! FOR PACKS OF 5 COLUMNS !!
        else
          write(*,'(2x,a15,5(i7,3x))') "Orbital Number:",(jj,jj=dd,dd+4)
          do jfrg=1,icufr
            write(*,'(2x,a9,i3,a1,5f10.5)') "Frg. Pop.",jfrg,":",(frgpop(jj,jfrg),jj=dd,dd+4)
          end do
          write(*,*) " "
          dd=dd+5
        end if
      end do

      end

!! ****** !!