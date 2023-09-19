      implicit real*8(a-h,o-z)
      include 'parameter.h'
      character*60 name,name2,namepat,name3,name0
      character*30 integ1,integ2
      character*80 linea ,line        
      common /filename/name0
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(100)
      dimension ilog(40)

CCCCCCCCCCCCCC
C PROCESS ARGUMENTS
CCCCCCCCCCCCCC

      CALL GETARG(1,name0)
      if(name0.ne."".and.name0.ne."%") then
       j=len(name0)     
       do i=1,j
        if(name0(i:i).eq.' ') then
         l=i-1
         go to 10
        end if
       end do
  10   continue
       name=name0(1:l)//".fchk"
       name3=name0(1:l)//".inp"
       namepat=name0(1:l)
      else
       name='Test.FChk'
      end if
      open (15,file=name)

CCCCCCCCCCCCCC
C PROCESS INP FILE
CCCCCCCCCCCCCC

      open (16,file=name3,err=2210)

      idoint=0
      iwfn=0  
      ndens0=1

c look for options      
      call readchar("# METODE","ORCA",iorca)
      call readchar("# METODE","MOBAS",imobas)
      if(imobas.eq.0) then
       call readint("# METODE","DENS",ndens0,1,1)
       iopt(9) = ndens0
      end if
      call readchar("# METODE","WFN",iwfn)
      call readchar("# METODE","ALLPOINTS",iallpo)
      call readchar("# METODE","NOCALCS",inocalcs)
      call readchar("# METODE","FULLPRECISION",iaccur)

C Atoms in molecules
      call readchar("# METODE","MULLI",imulli)
      call readchar("# METODE","LOWDIN",ilow)
      if(ilow.eq.1) imulli=2
      call readchar("# METODE","LOWDIN-DAVIDSON",ilow)
      if(ilow.eq.1) imulli=3
      call readchar("# METODE","HIRSH",ihirsh)
      call readchar("# METODE","HIRSH-IT",ihirsh0)
      if(ihirsh0.eq.1) ihirsh=2
      call readchar("# METODE","BECKE-RHO",ibcp)
      call readchar("# METODE","NEWBEC",inewbec)
      call readint("# METODE","STIFFNESS",istiff,3,1)
      call readchar("# METODE","TFVC",itfvc)
      if(itfvc.eq.1) then
        ibcp=1
        inewbec=1
        istiff=4
      end if
      call readchar("# METODE","WMATRAD",iradmat)
      call readchar("# METODE","RMATRAD",iradmat0)
      if(iradmat0.eq.1) iradmat=2

      call readchar("# METODE","ERF_PROF",ierf)
      call readreal("# METODE","ERF_PROF",aerf,6.266d0,1)
      if(ierf.eq.1)  istiff=0

c ERC QTAIM input module
      call readchar("# METODE","QTAIM",iqtaim)
      call readchar("# METODE","READINT",ireadint)
      if(ireadint.eq.1) iqtaim=2
      if(iqtaim.eq.1) then
       call readint("# QTAIM","STEP",istep,300,1)
       call readint("# QTAIM","NNA",inna,0,1)
       call readint("# QTAIM","MAXDIST",imaxdist,12,1)
       call readint("# QTAIM","SCREENING",iscreening,1000,1)
       call readint("# QTAIM","PATH",ipath,0,1)
      end if 

c Miscellaneous options
      call readchar("# METODE","OPOP",iopop)
      call readchar("# METODE","SHANNON",isha)
      call readchar("# METODE","DOINT",idoint)
      call readchar("# METODE","PCA",ipca)
      call readchar("# METODE","LAPLACIAN",ilaplacian)
      call readchar("# METODE","FINEGRID",ifinegrid)
 
C eff-AO-s and EOS
      call readchar("# METODE","EFFAO",ieffao)
      call readchar("# METODE","UEFFAO",idummy)
      if(idummy.eq.1) ieffao=2
      call readint("# METODE","EFF_THRESH",ieffthr,1,1)

      call readchar("# METODE","CUBE",icube)
      if(icube.eq.1) then
       call locate(16,"# CUBE",ii)
       if(ii.eq.0) stop'Required section # CUBE not found in input file'
       call readint("# CUBE","MAX_OCC",jcubthr,1000,1)
       call readint("# CUBE","MIN_OCC",kcubthr,0,1)
      end if

      call readchar("# METODE","EOS",ieos)
      if(ieos.eq.1) then 
       iopop=1
       call readreal("# METODE","EOS_THRESH",xthresh,2.5d-3,1)
       ieffao=2
      end if

c Local spin and methods for correlated WFs
      call readchar("# METODE","SPIN",ispin)
      call readint("# METODE","DM",icorr,0,1)
      call readchar("# METODE","PNOF",ipnof)
      call readchar("# METODE","DAFH",idafh)


C Energy partitioning options
      call readchar("# METODE","ENPART",ienpart )
      if(ienpart.eq.1) then
       xmix=0.0d0
       call readchar("# ENPART","HF ",ihf)
       if(ihf.eq.1) then
        imeth=-1
        xmix=1.0d0
        goto 233
       end if
       call readchar("# ENPART","LDA",ival )
       if(ival.eq.1) then 
        imeth=0 
        go to 233
       end if
       call readchar("# ENPART","B88",ival)
       if(ival.eq.1) then 
        imeth=2
        go to 233
       end if
c warning, order is important
       call readchar("# ENPART","B3LYP",ival)
       if(ival.eq.1) then
        imeth=5
        xmix=0.20d0
        go to 233
       end if
       call readchar("# ENPART","BLYP",ival)
       if(ival.eq.1) then
        imeth=1
        goto 233
       end if
       call readchar("# ENPART","LYP",ival)
       if(ival.eq.1) then    
        imeth=3
        go to 233
       end if
c warning, order is important
       call readchar("# ENPART","PBE0",ival)
       if(ival.eq.1) then
        imeth=6
        xmix=0.25d0
        go to 233
       end if
       call readchar("# ENPART","PBE",ival)
       if(ival.eq.1) then
        imeth=4
        goto 233
       end if
       call readchar("# ENPART","B3ONLY",ival)
       if(ival.eq.1) then
        imeth=7
        xmix=0.20d0
        goto 233
       end if
       stop ' No DFT/HF functional for energy decomposition '
233    continue

       call readint("# ENPART","THREBOD",ithrebod,100,1)
       call readchar("# ENPART","EXACT",iexact)
       call readchar("# ENPART","HOMO",ihomo)
       call readchar("# ENPART","DEKIN",idek)
       call readchar("# ENPART","IONIC",iionic)
       call readreal("# ENPART","TWOELTOLER",twoeltoler,0.25d0,1)
      end if
C NLOPs                       
      call readchar("# METODE","POLAR",ipolar )
      if(ipolar.eq.1) iaccur=1

CCCCCCCCCCCCCCCCCCCCCCC
C DEPENDENCIES & TO DO
CCCCCCCCCCCCCCCCCCCCCCC
      iposthf=0
      idono=0
      if(icas.eq.1.or.icisd.eq.1) iposthf=1
      if(iposthf.eq.1.or.kop.eq.1) idono=1

      if(iposthf.eq.1) then
       if(ispin.eq.1.and.icorr.lt.2) stop ' Local Spin needs dm1 and dm2 for correlated WFs'
       if(ieffao.eq.2.and.(icorr.lt.1.and.nalf.ne.nb)) stop ' EOS/ueff-AOs need dm1 for correlated non-singlet WFs'
      end if

      if (ipca.eq.1.and.iqtaim.ne.1) iopop=1
      if(imulli.gt.1.or.iqtaim.eq.1) iopop=0
      if(i5d.eq.1.and.imobas.eq.1)stop'MOBAS can not be used with 5d/7f'
      if(imulli.ne.0.and.inocalcs.eq.1) stop'MULLI cant do with NOCALCS'
      if(ihirsh.ne.0.and.idoatoms.eq.1) stop'HIRSH cant do with DOATOMS'
      if(ispin.eq.1.and.idono.eq.0)  then
       write(*,*) 'No Local Spin Analysis needed for Restricted SD WFs'
       ispin=0
      end if
      if(imulli.ne.0.and.ienpart.ne.0)  then
       write(*,*) 'Can not do ENPART with Hilbert-space analysis'
       ienpart=0
      end if
      if(imulli.ne.0.and.icube.ne.0.and.idofr.ne.0)  then
       write(*,*) 'Can not do Hilbert-space cubes using fragments'
       icube=0
      end if

CCCCCCCCCCCCCC
C END PROCESS INP FILE
CCCCCCCCCCCCCC


CCCCCCCCCCCCC
c OPTIONS LIST
CCCCCCCCCCCCC
      iopt(1) = idoint
      iopt(2) = imobas
      iopt(3) = iwfn
      iopt(4) = idono
      iopt(5) = imulli
      iopt(6) = ihirsh
      iopt(7) = iallpo
      iopt(8) = iopop
      iopt(9) = ndens0
      iopt(10) = inocalcs
c     iopt(11) = i5d
      iopt(12) = ieffao
      iopt(13) = icube 
      iopt(14) = ibcp  
      iopt(15) = ispin 
      iopt(16) = iqtaim 
      iopt(17) = ienpart
      iopt(18) = ihf
      iopt(19) = imeth
      iopt(20) = iexact
      iopt(21) = ihomo
      iopt(22) = idek
      iopt(23) = iionic
      iopt(24) = ieffthr
      iopt(25) = istiff 
      iopt(26) = icorr  
      iopt(27) = isha   
      iopt(28) = ipca   
      iopt(29) = ipnof  
      iopt(30) = idafh  
      iopt(31) = inewbec
c      iopt(32) =   
c      iopt(33) =  
      iopt(34) = ilaplacian
c ERC qtaim 
      iopt(35) = istep 
      iopt(36) = inna
      iopt(37) = imaxdist
      iopt(38) = iscreening
      iopt(39) = ipath
c
      iopt(40) = idofr
      iopt(41) = jcubthr
      iopt(42) = kcubthr
      iopt(43) = iorca   
      iopt(44) = ithrebod
      iopt(45) = ifinegrid
c
      iopt(46) = ipolar
      iopt(47) = ifield
      iopt(48) = iradmat

CCCCCCCCCCCCCCCCCCC
C PREPARE FOR NUMERICAL INTEGRATIONS
CCCCCCCCCCCCCCCCCCC
       call getarg(2,integ1)
       call getarg(3,integ2)
       if(integ1.ne.' '.and.integ2.ne.' ') then
        read(integ1,'(i4)') Nrad        
        read(integ2,'(i4)') Nang        
        if(Nrad.gt.500) stop 'Max number of radial points  is 500 '
        ilog(1)=6
        ilog(2)= 14 
        ilog(3)= 26
        ilog(4)= 38 
        ilog(5)= 50
        ilog(6)= 74 
        ilog(7)= 86 
        ilog(8)= 110
        ilog(9)= 146
        ilog(10)= 170
        ilog(11)= 194
        ilog(12)= 230
        ilog(13)= 266
        ilog(14)= 302
        ilog(15)= 350
        ilog(16)= 434
        ilog(17)= 590
        ilog(18)= 770
        ilog(19)= 974
        do 111 i=1,18
        npoints=ilog(i)
        if(nang.lt.ilog(i+1)) goto 211
  111 continue
        npoints=ilog(19)
  211 continue
          print *,' Angular points:',npoints
          nang=npoints
          rr00=0.500d0
       else if(ienpart.eq.1.or.ipolar.eq.1) then
C Enpart defaults for high-accuracyone-el integrations 
        nrad=150
        nang=590
        if(ifinegrid.eq.1) nang=974
        rr00=0.500d0
       else
C APOST old defaults for one-el integrations 
        nrad=40 
        nang=146
        rr00=0.5d0
       end if

CCCCCCCCCCCCC
C PROCESS FCHK FILE
CCCCCCCCCCCCC
      rewind(15)
22    read(15,'(a80)')line
      if (index(line,"Number of atoms").ne.0) then
       read(line(45:),*) nat
      else
       go to 22
      end if
26    read(15,'(a80)')line
      if (index(line,"Multiplicity").ne.0) then
       read(line(45:),*) kop0
      else
       go to 26
      end if
      kop=0
      if(kop0.ne.1) kop=1
23    read(15,'(a80)')line
      if (index(line,"Number of alpha").ne.0) then
       read(line(45:),*) nalf 
      else
       go to 23
      end if
24    read(15,'(a80)')line
      if (index(line,"Number of beta").ne.0) then
       read(line(45:),*) nb   
      else
       go to 24
      end if
25    read(15,'(a80)')line
      if (index(line,"Number of basis").ne.0) then
       read(line(45:),*) igr  
      else
       go to 25
      end if
      close(15)

CCCCCCCCCCCCC
C END PROCESS FCHK FILE
CCCCCCCCCCCCC

CCCC
C ESTIMATED MEMORY USAGE
CCCC
      iatps=nrad*nang
      itotps=iatps*nat
      write(*,*) kop,igr,nat,nalf,nb,iatps
      xmem0=(9.3d7+(12*nmax*nmax+29*maxat*maxat)*8)/(1024.0**3)
      write(*,*)
      write(*,'(a26,f8.2)') ' Static Memory aprox. (GB)',xmem0
c dynamic memory..
      xmem=0.0d0
      if(imulli.lt.1) xmem=xmem+igr*igr*nat*2+iatps*nat*(6+igr+nat)
      if(ispin.ne.0) xmem=xmem+igr*igr*(nat+1)
      if(ieffao.ne.0) xmem=xmem+igr*(igr*5+2)+iatps*nat
      ymem=0.0d0
      if(ienpart.ne.0) then
        ymem=-igr*igr*nat*2
        ymem=ymem+(iatps)*nat*(nalf+1+igr) +40*146*nat*(2*nalf+igr+nat+8)
        if(kop.eq.1) ymem=ymem+ iatps*nalf+40*146*nat*(2*nalf+igr+nat+8)
        if(imeth.ne.-1) ymem=ymem+iatps*nat*(6+4*nalf+2*igr)+nalf*nalf*nat
        if(kop.eq.1) ymem=ymem+iatps*nat*(3+4*nb)+nb*nb*nat
      end if
      xmem=(xmem+ymem)*8.0/(1024.0**3)
      write(*,'(a33,f8.2)') ' Max allocated memory aprox. (GB)',xmem
      write(*,'(a31,f8.2)') ' Total memory usage aprox. (GB)',xmem+xmem0
      write(*,*)' '
CCCC
C ESTIMATED MEMORY USAGE
CCCC
2210      end
