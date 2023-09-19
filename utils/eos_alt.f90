        program eos
        implicit double precision(a-h,o-z)
        parameter(maxfrag=10,maxeff=100,maxval=120,maxconf=300000,maxdist=1000)
        integer nat,nelec(2)
        character*80 line
        character*80 name,val1,val2,val3
        dimension effao(maxfrag*maxeff,2),ieffao(maxfrag*maxeff,2),neffao(2)
        dimension ncore(2),nval(2),neval(2)
        dimension kpop(maxfrag,2)
        dimension ivect(maxval),jvect(maxval),svect(maxval),fvect(maxval)
        dimension lconf(maxval,maxconf),www(maxconf)
        dimension ileos(maxdist),weos(maxdist),ieos(maxdist,2*maxfrag), weos0(maxdist)
        dimension ileos_t(maxdist),weos_t(maxdist),ieos_t(maxdist,maxfrag), nntot_t(maxfrag)
        dimension ipop(maxfrag,2)
        dimension iznuc(maxfrag)
        dimension isort(maxdist)
        common /range/xmax,xmin
!
        dimension iprim(25)
        data iprim /2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97/

        call getarg(1,name)
        call getarg(2,val1)
        call getarg(3,val2)
        call getarg(4,val3)
        read(val1,*) xmax
        read(val2,*) xmin        
! pseudo T for Boltzmann-like distribution
        tt=70.0d0
        read(val3,*) tt          
        write(*,*) 'Temperature: ',tt

        open(1,file=name)
        rewind(1)

1       read(1,'(A80)') line
        ii=index(line,"nalpha")
        if(ii.eq.0) go to 1
        read(line(ii:),*)name,nelec(1),name,nelec(2)
        write(*,*) 'Number of alpha electrons:',nelec(1)
        write(*,*) 'Number of alpha electrons:',nelec(2)
        if(index(line,"nocc").eq.0) then 
         write(*,*) 'Multireference calculation'
        else
         write(*,*) 'Single-determinant calculation'
        end if

! compatbility with different outputs of apost3d
        icorr=1

11       read(1,'(A80)') line
        if(index(line,"Number of atoms").eq.0) go to 11
        read(line(24:),*)nat 
        write(*,*) 'Number of atoms:',nat

!  in case no fragments
        nfrag=nat
        nfrag0=0
2       read(1,'(A80)',end=20) line
        if(index(line,"Fragment").ne.0) then
         nfrag0=nfrag0+1
        else if (nfrag0.ne.0.and.index(line,"DOING POPULATION ").ne.0) then
         go to 20
        end if
        go to 2
20      continue
        if(nfrag0.ne.0)nfrag=nfrag0
       
21     continue

       write(*,*) 'Number of fragments :',nfrag

! inferring Znuc from population analysis
       call getz(nfrag,iznuc)
       ix=0
       do i=1,nfrag
        ix=ix+iznuc(i)
       end do
       if(ix.ne.nelec(1)+nelec(2)) then
        write(*,*) 'Systems appears to be charged ',ix-nelec(1)-nelec(2)
        write(*,*)
       end if

       write(*,*) 'Processing alpha eff-AOs...'
       call readeffao(1,icorr,nfrag,effao,ieffao,neffao)

! check if beta
       ibeta=0 
5     read(1,'(A80)',end=55) line
      if(index(line,"SKIPPING").ne.0) then 
       ibeta=0
       go to 55
      else if(index(line,"BETA   PART").ne.0) then 
       ibeta=1
       go to 55
      else
       go to 5
      end if
55    continue

      if(ibeta.eq.1) then
       write(*,*) 'Processing beta eff-AOs...'
       call readeffao(2,icorr,nfrag,effao,ieffao,neffao)
      else
       write(*,*) 'Copying alpha to beta eff-AOs...'
       neffao(2)=neffao(1)
       do i=1,neffao(1)
        effao(i,2)=effao(i,1)
        ieffao(i,2)=ieffao(i,1)
       end do
      end if

! start process
      do i=1,nfrag
       kpop(i,1)=0
       kpop(i,2)=0
      end do

      do icase=1,2
       ncore(icase)=0
       nval(icase)=0
       do i=1,neffao(icase)
        if(effao(i,icase).ge.xmax) then
         ncore(icase)=ncore(icase)+1
         kpop(ieffao(i,icase),icase)= kpop(ieffao(i,icase),icase)+1
        else if(effao(i,icase).ge.xmin) then
         nval(icase)=nval(icase)+1
        end if
       end do
       neval(icase)=nelec(icase)-ncore(icase)
      end do

      write(*,'(a19,2i4)') 'Number of eff-AOs :',(neffao(icase),icase=1,2)
      write(*,'(a19,2i4)') 'Total electrons   :',(nelec(icase),icase=1,2)
      write(*,'(a19,2i4)') 'Core electrons    :',(ncore(icase),icase=1,2)
      write(*,*) ' Core populations  :'
      write(*,*) '                    FRAG  ALPHA  BETA '
      write(*,*) '                    ------------------'
      do ifrag=1,nfrag
       write(*,'(20x,i3,2x,i4,x,i4)') ifrag,(kpop(ifrag,icase),icase=1,2)
      end do
      write(*,*) '                    ------------------'

      nevaltot=neval(1)+neval(2)
      norbtot=nval(1)+nval(2)

      write(*,'(a19,3i4)') 'Valence electrons :',(neval(icase),icase=1,2),nevaltot
      write(*,'(a19,3i4)') 'Valence orbitals  :',(nval(icase),icase=1,2),norbtot

      if(nval(1)+nval(2).eq.0) stop 'NO VALENCE eff-AOs IN THE RANGE SELECTED'
      write(*,*)  
      write(*,*) ' ********************************************************'
      write(*,123) 'eff-AO range of (',xmax,',',xmin,') leads to a eff-CAS (',nevaltot,',',norbtot,')'
      write(*,*) ' ******************************************************** '
      write(*,*)  
123   FORMAT(a19,f5.3,a,f5.3,a22,i3,a,i3,a)



! NESTED NSS 
       n=norbtot
       nconf=0

       ipos=0
       do k=1,n
          ivect(k)=0   
          fvect(k)=1
          svect(k)=1   
          jvect(k)=ivect(k)
       end do
       k=n        
       do while (k.gt.0)
        if ((jvect(k)-fvect(k))*svect(k).gt.0) then 
           jvect(k)=ivect(k)
           k=k-1
        else 
!
! checksum for alpha and beta electrons
         ii=0
         ilog=0
         do i=1,nval(1)
          ii=ii+jvect(i)
         end do
         if(ii.eq.neval(1)) then
          ii=0
          do i=nval(1)+1,n
           ii=ii+jvect(i)
          end do
          if(ii.eq.neval(2)) ilog=1 
         end if
         if(ilog.eq.1) then 
          nconf=nconf+1
          do i=1,n
           lconf(i,nconf)=jvect(i)
          end do
         end if
!
         k=n
        end if
        if (k.gt.0) jvect(k)=jvect(k)+svect(k)  
       end do
       
!
! processing configurations. Calculating weights
!
       wwtot=0.0d0
       do ic=1,nconf
        ww=1.0d0
        do i=1,nval(1)
         if(lconf(i,ic).eq.1) ww=ww*effao(ncore(1)+i,1)
        end do
        do i=nval(1)+1,n
         if(lconf(i,ic).eq.1) ww=ww*effao(ncore(2)+i-nval(1),2)
        end do
! 
        www(ic)=exp(ww*tt)
!
        wwtot=wwtot+www(ic)
       end do
       do i=1,nconf
        www(i)=www(i)/wwtot
       end do

       write(*,*) 'Total number of configurations :',nconf
       do ic=1,nconf
         write(*,'(a2,f8.5,a2,40(x,i2))') 'W=',www(ic),'::',(lconf(i,ic),i=1,n)
       end do
       write(*,*) 

       write(*,*) ' Transforming configurations to distinct spin-resolved EOS distributions...'
       write(*,*) 

      neos=0
      do ic=1,nconf

       do i=1,nfrag   
        ipop(i,1)=0
        ipop(i,2)=0
       end do
       do i=1,nval(1)
        if(lconf(i,ic).eq.1) ipop(ieffao(ncore(1)+i,1),1)= ipop(ieffao(ncore(1)+i,1),1)+1
       end do
       do i=nval(1)+1,n
       if(lconf(i,ic).eq.1) ipop(ieffao(ncore(2)+i-nval(1),2),2)=ipop(ieffao(ncore(2)+i-nval(1),2),2)+ 1
       end do

       index_eos=1  ! assign unique index to each eos distribution
       do i=1,nfrag
        index_eos=index_eos*iprim(2*i-1)**ipop(i,1)*iprim(2*i)**ipop(i,2)
       end do
       if(index_eos.lt.0) stop 'integer overflow in index_eos'

       ilog=0
       do i=1,neos
        if(index_eos.eq.ileos(i)) then
         ilog=1
         j=i
        end if
       end do
       if (ilog.eq.1) then    !existing eos distribution, j
         weos(j)=weos(j)+www(ic)
       else  !new eos distribution
        neos=neos+1
        ileos(neos)=index_eos
        weos(neos)=www(ic)
        do i=1,nfrag
         ieos(neos,2*i-1)=ipop(i,1)
         ieos(neos,2*i)=ipop(i,2)
        end do
        write(*,'(a6,i5,x,a2,i8,x,a2,30(x,i2))') 'Dist. ',ic,'ID',  index_eos,'::',(ipop(i,1),ipop(i,2),i=1,nfrag) 
       end if
      end do

! sort by weight
       do i=1,neos
        weos0(i)=weos(i)
        isort(i)=i
       end do
       do i=1,neos
        wmax=weos0(i)
        imax=i
        do j=i+1,neos
         if(weos0(j).gt.wmax) then
           wmax=weos0(j)
           imax=j
         end if
        end do
        weos0(imax)=weos0(i)
        weos0(i)=wmax 
        ii=isort(imax)
        isort(imax)=isort(i)
        isort(i)=ii         
       end do
!       write(*,*) 'sort array', (isort(i),i=1,neos)


! neos : number of eos distributions
! weos(:) : weight of each eos distributions
! ileos(:) : unique index of each distribution
! ieos(i,j) : eos distribution i in alpha,beta alternate for each fragment j

      write(*,*)
      write(*,*) '*********************************************'
      write(*,*) '         SPIN-RESOLVED RESULTS '
      write(*,'(a40,i4)') 'Number of distinct EOS distributions:',neos
      write(*,*) '*********************************************'
      do i=1,neos
       write(*,*) 
       write(*,'(i3,x,a18,i3,a9,f6.4)') i,'EOS distribution :',isort(i), ' Weight :',weos(isort(i))
       write(*,*) 
       write(*,*) '  FRAG    N_alpha   N_beta   N_tot   EOS   '
       write(*,*) '------------------------------------------'
       do ifrag=1,nfrag
        nna=kpop(ifrag,1)+ieos(isort(i),2*ifrag-1)
        nnb=kpop(ifrag,2)+ieos(isort(i),2*ifrag)
        nntot=nna+nnb
        write(*,'(5(i6,3x))') ifrag, nna, nnb,nntot, -nntot+iznuc(ifrag)
       end do
       write(*,*) '------------------------------------------'
      end do
      write(*,*) 

      write(*,*) 'Getting unique distributions of EOS (spinless)...'
      write(*,*) 


! neos_t : number of eos distributions
! weos_t(:) : weight of each eos distributions
! ileos_t(:) : unique index of each distribution
! ieos_t(i,j) : eos distribution i  for each fragment j


      neos_t=0
      do ic=1,neos 

       do ifrag=1,nfrag
        nntot_t(ifrag)=ieos(ic,2*ifrag)+ieos(ic,2*ifrag-1)
       end do

       index_eos_t=1  ! assign unique index to each eos distribution
       do i=1,nfrag
        index_eos_t=index_eos_t*iprim(i)**nntot_t(i)
       end do
       if(index_eos_t.lt.0) stop 'integer overflow in index_eos'

       ilog=0
       do i=1,neos_t
        if(index_eos_t.eq.ileos_t(i)) then
         ilog=1
         j=i
        end if
       end do
       if (ilog.eq.1) then    !existing eos_t distribution, j
         weos_t(j)=weos_t(j)+weos(ic)
       else  !new eos_t distribution
        neos_t=neos_t+1
        ileos_t(neos_t)=index_eos_t
        weos_t(neos_t)=weos(ic)
        do i=1,nfrag
         ieos_t(neos_t,i)=nntot_t(i)
        end do
        write(*,'(a6,i5,x,a2,i8,x,a2,30(x,i2))') 'Dist. ',ic,'ID',  index_eos_t,'::',(nntot_t(i),i=1,nfrag) 
       end if
      end do

! sort by weight
       do i=1,neos_t
        weos0(i)=weos_t(i)
        isort(i)=i
       end do
       do i=1,neos_t
        wmax=weos0(i)
        imax=i
        do j=i+1,neos_t
         if(weos0(j).gt.wmax) then
           wmax=weos0(j)
           imax=j
         end if
        end do
        weos0(imax)=weos0(i)
        weos0(i)=wmax 
        ii=isort(imax)
        isort(imax)=isort(i)
        isort(i)=ii         
       end do
!       write(*,*) 'sort array', (isort(i),i=1,neos_t)


      write(*,*) 
      write(*,*) '******************************************'
      write(*,*) '               EOS  RESULTS '
      write(*,'(a38,i4)') 'Number of unique EOS distributions:',neos_t
      write(*,*) '******************************************'
      do i=1,neos_t
       write(*,*) 
       write(*,'(i3,x,a18,i3,a9,f6.4)') i,'EOS distribution :',isort(i), ' Weight :',weos_t(isort(i))
       write(*,*) 
       write(*,*) '  FRAG    N_tot      EOS   '
       write(*,*) '------------------------------------------'
       do ifrag=1,nfrag
        nntot=ieos_t(isort(i),ifrag)+kpop(ifrag,1)+kpop(ifrag,2)
        write(*,'(3(i6,3x))') ifrag, nntot, -nntot+iznuc(ifrag)
       end do
       write(*,*) '------------------------------------------'
      end do

! Average configuration..if it makes any sense
      write(*,*) 
      write(*,*) '******************************************'
      write(*,*) '               EOS  RESULTS '
      write(*,*) '           AVERAGE CONFIGURATION'
      write(*,*) '******************************************'
      write(*,*) 
      write(*,*) '  FRAG   N_tot    EOS   '
      write(*,*) '------------------------'
      do ifrag=1,nfrag
       xx=0.0d0
       do i=1,neos_t
        xx=xx+ieos_t(i,ifrag)*weos_t(i)
       end do
       xx=xx+kpop(ifrag,1)+kpop(ifrag,2)
       write(*,'(2x,i3,2x,2f8.3)') ifrag, xx, -xx+iznuc(ifrag)
      end do
      write(*,*) '------------------------'

      end
        




        subroutine readeffao(icase,icorr,nfrag,effao,ieffao,neffao)
        implicit double precision(a-h,o-z)
        parameter(maxfrag=10,maxeff=100,maxval=120,maxconf=300000,maxdist=1000)
        integer nat,nalpha, nbeta,nocc
        character*80 line
        character*80 name
        character*16 string
        logical ilog
        dimension effao(maxfrag*maxeff,2),ieffao(maxfrag*maxeff,2),neffao(2)
        dimension nocc(maxfrag),occ(maxeff,maxfrag)
        dimension icount(maxfrag*maxeff)
        common /range/xmax,xmin

       string="Net occupations "
       if(icorr.eq.1) string="et occupation us"

       rewind(1)
1      read(1,'(A80)',err=1000) line
       if(index(line,"UEFFAO").eq.0) go to 1

       if(icase.eq.2) then
2      read(1,'(A80)',err=1000) line
       if(index(line,"BETA   PART").eq.0) go to 2
       end if

       do ifrag=1,nfrag

3      read(1,'(A80)') line
       if(index(line,(string)).eq.0) go to 3
       j0=1
4      read(1,'(A80)') line
       if(index(line,"OCCUP").ne.0) then
        read(line(8:),*,end=44)(occ(j,ifrag),j=j0,j0+8)
44      j0=j0+8
        go to 4
       else
        nocc(ifrag)=j0-1
       end if
       write(*,*) ' '    
       write(*,*) 'eff-AOs processed for fragment ', ifrag
       write(*,*) 'nocc(ifrag)',nocc(ifrag)
       write(*,'(10f8.5)')( occ(i,ifrag),i=1,nocc(ifrag))
      end do

!! Sorting valence effaos in decreasing order
!! Core eff-Aos are assumed for occ >xmax
       neffao(icase)=0
       do ifrag=1,nfrag
        icount(ifrag)=1
       end do

! copying core/lone pairs
       do ifrag=1,nfrag
        do j=1,nocc(ifrag)
         if(occ(j,ifrag).gt.xmax) then
          neffao(icase)=neffao(icase)+1
          effao(neffao(icase),icase)=occ(j,ifrag)
          ieffao(neffao(icase),icase)=ifrag
          icount(ifrag)=icount(ifrag)+1
         end if
        end do
       end do

! dealing with valence effaos
      xminocc=xmin
      ilog=.true.
      do while(ilog)
       occmax=xminocc
       iifrag=0
       iieff=0
       do ifrag=1,nfrag
        do jeff=icount(ifrag),nocc(ifrag)
         if(occ(jeff,ifrag).gt.occmax) then
          occmax=occ(jeff,ifrag)
          iifrag=ifrag
          iieff=jeff
         end if
        end do
       end do
       if(iifrag.eq.0) then 
        ilog=.false.
       else
       neffao(icase)=neffao(icase)+1
       effao(neffao(icase),icase)=occmax
       ieffao(neffao(icase),icase)=iifrag
       icount(iifrag)=icount(iifrag)+1
      end if
      end do

      write(*,*) ' '
      write(*,*) 'eff-AOs'
      do i=1,neffao(icase)
        write(*,'(f8.5,i3)') effao(i,icase),ieffao(i,icase)
      end do
      write(*,*) ' '
       


1000    continue
        end 

        subroutine getz(nfrag,iznuc)
        implicit double precision(a-h,o-z)
        parameter(maxfrag=10,maxeff=100,maxval=120,maxconf=300000,maxdist=1000)
        character*80 line
        character*4 name  
        dimension iznuc(maxfrag)
        dimension pop(maxfrag)

1       read(1,'(A80)',err=100) line
        ii=index(line," FRAGMENT ANALYSIS : Electron populations")
        if(ii.eq.0) go to 1
        read(1,'(A80)') line
        read(1,'(A80)') line
        read(1,'(A80)') line
   
        do i=1,nfrag
         read(1,*)j, xx,pop(i)
        end do

2       read(1,'(A80)') line
        ii=index(line," FRAGMENT ANALYSIS : Atomic Charges")
        if(ii.eq.0) go to 2
        read(1,'(A80)') line
        read(1,'(A80)') line
        read(1,'(A80)') line
        do i=1,nfrag
         read(1,*)j, xy,xx
         pop(i)=pop(i)+xx
         iznuc(i)=nint(pop(i))
        end do
        do i=1,nfrag
         write(*,*) pop(i),iznuc(i)
        end do
        return
100     continue
! Fragments not defined
        rewind(1)
3       read(1,'(A80)') line
        ii=index(line,"ELECTRON POPULATIONS")
        if(ii.eq.0) go to 3
        read(1,'(A80)') line
        read(1,'(A80)') line
        read(1,'(A80)') line
        write(*,'(A80)') line
        do i=1,nfrag
         read(1,*)j, name,xx,pop(i)
        end do

4       read(1,'(A80)') line
        ii=index(line,"ATOMIC CHARGES")
        if(ii.eq.0) go to 4
        read(1,'(A80)') line
        read(1,'(A80)') line
        read(1,'(A80)') line
        do i=1,nfrag
         read(1,*)j, name,xy,xx
         pop(i)=pop(i)+xx
         iznuc(i)=nint(pop(i))
        end do
        do i=1,nfrag
         write(*,*) pop(i),iznuc(i)
        end do

        end 
