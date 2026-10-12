!! ********************************************************************** !!
!! .fchk READING, DM1/DM2 INPUT, AND .inp KEYWORD PARSING                 !!
!! .fchk reader (native Gaussian/ORCA wavefunction file, unit 15):        !!
!!   input           -- MO formalism/type, basis size, P/PS/spin density, !!
!!                      CAS/CI spin-orbital counts                        !!
!!   mga_misc        -- reference KE/E-N/E-E energies (project's own      !!
!!                      post-processing fields) + ECP matrix if present   !!
!!   field_misc      -- external static electric field, if present        !!
!! post-HF 1-/2-RDM input (# DM keyword, PySCF/DMN codes:                 !!
!!   dm1input        -- spin-orbital 1-RDM + natural orbitals + P/Ps      !!
!!                      reconstruction from it                            !!
!!   dm2input_pyscf  -- active-space 2-RDM from PySCF's write_dm12        !!
!!   dm2input_dmn    -- full spin-separated 2-RDM from the DMN code       !!
!! .inp keyword parsers (unit 16, rewind+rescan per call and their        !!
!! shared file-scanning                                                   !!
!! primitives:                                                            !!
!!   readchar        -- bare-keyword presence check                       !!
!!   readreal        -- keyword + real value, with default                !!
!!   readint         -- keyword + integer value, with default             !!
!!   locate          -- rewind + scan for a matching line (used by the    !!
!!                      three above and by other files' .inp readers)     !!
!!   int_locate      -- rewind + scan + fixed-column integer read (.fchk) !!
!!   real_locate     -- rewind + scan + fixed-column real*8 read (.fchk)  !!
!! ********************************************************************** !!

!! ***** !!

!! ********************************************************************* !!
!! subroutine: input                                                     !!
!! purpose: reads the native Gaussian .fchk (unit 15) into the           !!
!! wavefunction's global state -- MO formalism/type, basis size (via     !!
!! build_basis/build_ao_matrices), P/PS density matrices, spin density,  !!
!! CAS/CI spin-orbital counts. Prints the WAVEFUNCTION SUMMARY block as  !!
!! each value is determined (see the print_box call below).              !!
!! arguments: none (unit 15 already open; output via COMMON /nat/,       !!
!!   /coord/, /energ/, /cas/ and ao_matrices' P/PS/C/CB module arrays)   !!
!! author:                                                               !!
!! ********************************************************************* !!
      subroutine input
      use basis_set
      use ao_matrices
      use effao_mod, only: allocate_effao
      use nao_mod, only: allocate_nao
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /coord/ coord0(3,maxat),zn(maxat),iznuc(maxat)
      common /iops/iopt(200)
      common /energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      dimension clin(nmax**2)
      character*80 line,line2
      character*20 cform,ctype
      logical ilog

      ndens0=iopt(9)

      irhf=0
      iuhf=0
      icas=0
      icisd=0
      irohf=0

!! wavefunction summary -- each line below prints as its value is determined !!
      call print_box('WAVEFUNCTION SUMMARY')

!! start processing fchk file !!
      rewind 15
      read(15,'(a80)')line
      read(15,'(a80)')line
      if(index(line(11:11),"R").ne.0) then
        irhf=1
        cform='Restricted'
      else if(index(line(11:11),"U").ne.0) then
        cform='Unrestricted'
        iuhf=1
      else

!! no R/U in the header (pySCF CASSCF, FCI, CCSD): unrestricted only if  !!
!! the file has Beta MO coefficients, as for kop below                   !!
        cform='Restricted'
        do while(.true.)
          read(15,'(a80)',end=5) line2
          if(index(line2,"Beta MO co").ne.0) then
            cform='Unrestricted'
            exit
          end if
        end do
5       rewind 15
        read(15,'(a80)') line2
        read(15,'(a80)') line2
      end if
      write(*,'(2x,a,1x,a)') 'MO formalism                     :',
     +  trim(cform)

      if(index(line(11:25),"CASSCF").ne.0) then
        icas=1
        ctype='CASSCF'
      else if(index(line(11:25),"HF").ne.0) then
        ctype='Hartree-Fock'
      else if(index(line(11:25),"MP").ne.0) then
        ctype='Moller-Plesset PT'
      else if(index(line(11:25),"CI").ne.0) then
        icisd=1
        ctype='CI'
      else if(index(line(11:25),"CCSD").ne.0) then
        ctype='Coupled-cluster'
        icisd=1
      else
        ctype='KS-DFT'
      end if
      write(*,'(2x,a,1x,a)') 'Calculation type                 :',
     +  trim(ctype)
      if(ctype.eq.'KS-DFT') write(*,'(4x,a)')
     +  '(single det. built from the Kohn-Sham orbitals)'

      call build_basis()

      igr=nbasis
      nat=natoms

      call build_ao_matrices(igr) !MMO- WIP!

!! EFO and natural-AO arrays, sized now that igr and nat are known    !!
      call allocate_effao(igr,nat)
      call allocate_nao(igr)

      igr0=int_locate(15,"Number of independ",ilog)
      if(igr.ne.igr0) write(*,'(2x,a)')
     +  'WARNING: some basis functions were removed (linear dependency)'

      nelectr=int_locate(15,"Number of electr",ilog)

      nalf=int_locate(15,"Number of alpha electr",ilog)

      nb=int_locate(15,"Number of beta electr",ilog)

      dummy=int_locate(15,"Atomic numbers",ilog)
      read(15,*)(iznuc(i),i=1,nat)   

      dummy=int_locate(15,"Nuclear charg",ilog)
      read(15,*)(zn(i),i=1,nat)

      dummy=int_locate(15,"Current cartesian",ilog)
      read(15,*)(coord0(1,i),coord0(2,i),coord0(3,i),i=1,nat)

      escf=real_locate(15,"SCF Energ",ilog)
      ss2=real_locate(15,"S**2    ",ilog)

      dummy=int_locate(15,"Alpha MO co",ilog)
      read(15,*)((c(i,j),i=1,igr),j=1,igr0)
      if(igr0.ne.igr) then
        do i=igr0+1,igr
          do j=1,igr
            c(j,i)=0.0d0
          end do
        end do
      end if

      dummy=int_locate(15,"Beta MO co",ilog)
      if(ilog) then
        read(15,*)((cb(i,j),i=1,igr),j=1,igr0)
        if(igr0.ne.igr) then
          do i=igr0+1,igr
            do j=1,igr
              cb(j,i)=0.0d0
            end do
          end do
        end if

!! just in case not detected above !!
        iuhf=1
        kop=1
      end if

!! checking for RO case -- wrong P provided in FChk in some Gaussian    !!
!! versions, building P and PS from MOs just in case                    !!
      irohf=int_locate(15,"IROHF",ilog)
      if(irohf.eq.1) then
        kop=1
        do i=1,igr
          do j=1,igr
            cb(i,j)=c(i,j)
          enddo
        enddo
      end if

!! reading P-matrix from fchk !!
      ntriang=int_locate(15,"otal SCF D",ilog)
      read(15,*)(clin(i),i=1,ntriang)

!! DENS n: the n-th density of the .fchk, counting the blocks whose      !!
!! header starts with "Total" and names a "Density" (the same rule for   !!
!! Gaussian, Q-Chem and pySCF files)                                    !!
      if(ndens0.gt.1) then
        ndens=0
        rewind 15
 200    read(15,'(a80)',end=211) line
        if(line(1:5).eq.'Total'.and.index(line,'Density').ne.0) then
          ndens=ndens+1
          write(*,'(2x,a,i0,a,a)') 'Density ',ndens,' in the .fchk: ',
     +      trim(line(1:40))
        end if
        goto 200
 211    continue
!! the density of an unrestricted or ROHF wavefunction is rebuilt from  !!
!! its SCF orbitals below, so another one cannot be used                !!
        if(kop.eq.1) then
          write(*,'(2x,a)') 'DENS above 1 is only available for '//
     +      'restricted wavefunctions'
          call apost_stop(' DENS: unrestricted wavefunction')
        end if
        if(ndens0.gt.ndens) then
          write(*,'(2x,a,i0,a,i0,a)') 'DENS ',ndens0,
     +      ' asked for, the number of densities in the .fchk is ',ndens
          call apost_stop(' DENS: density not found')
        end if
        ncou=0
        rewind 15
        do while(ncou.lt.ndens0)
          read(15,'(a80)') line
          if(line(1:5).eq.'Total'.and.index(line,'Density').ne.0)
     +      ncou=ncou+1
        end do
        write(*,'(2x,a,i0,a,a)') 'Using density ',ndens0,': ',
     +    trim(line(1:40))
        read(15,*)(clin(i),i=1,ntriang)
      end if

      ncou=1
      do i=1,igr
        do j=1,i
          p(i,j)=clin(ncou)
          p(j,i)=clin(ncou)
          ncou=ncou+1
        enddo
      enddo

      if(kop.eq.1) then
        do i=1,igr
          do j=1,igr
            pa(i,j)=0.0d0
            pb(i,j)=0.0d0
            do ij=1,nalf
              pa(i,j)=pa(i,j)+c(i,ij)*c(j,ij)
            end do
            do ij=1,nb
              pb(i,j)=pb(i,j)+cb(i,ij)*cb(j,ij)
            end do
            p(i,j)=pa(i,j)+pb(i,j)
            ps(i,j)=pa(i,j)-pb(i,j)
          end do
        end do
      else
        do i=1,igr
          do j=1,igr
            pa(i,j)=p(i,j)/2.0d0 
            pb(i,j)=p(i,j)/2.0d0 
          enddo
        enddo
      end if

!! checking if spin density is available in fchk !!
      ntriang=int_locate(15,"Spin SCF ",ilog)
      if(ilog) then
        read(15,*)(clin(i),i=1,ntriang)
        ncou=1
        do i=1,igr
          do j=1,i
            ps(i,j)=clin(ncou)
            ps(j,i)=clin(ncou)
            ncou=ncou+1
          enddo
        enddo
        write(*,'(2x,a,1x,a)') 'Spin density (FChk)              :',
     +    'found'
      else if(icas.ne.1.and.icisd.ne.1) then
        write(*,'(2x,a,1x,a)') 'Spin density (FChk)              :',
     +    'not found, reconstructed from MOs'
      end if

!! in case CASSCF !!
      ncasel=int_locate(15,"Number of CAS El",ilog)
      if(ilog) then
        icas=1
        ncasorb=int_locate(15,"Number of CAS Or",ilog)
        write(*,'(2x,a,1x,i0)') 'CAS active electrons             :',
     +    ncasel
        write(*,'(2x,a,1x,i0)') 'CAS active orbitals              :',
     +    ncasorb
      endif

      nocc=nelectr/2
      if(icas.eq.1) then
        ncore=(nelectr-ncasel)/2
        norb=ncore+ncasorb
        nspinorb=norb*2

        write(*,'(2x,a,1x,i0)') 'CAS spin-orbitals                :',
     +    nspinorb
        icass=0
        if(abs(ss2).gt.1.0d-2) then
          write(*,'(2x,a,1x,a)') 'CAS spin state                   :',
     +      'non-singlet'
          icass=1
        end if
      else if(icisd.eq.1) then
        nspinorb=igr0*2
        write(*,'(2x,a,1x,i0)') 'CI spin-orbitals                 :',
     +    nspinorb
      endif

      if(icas.eq.0)then
        write(*,'(2x,a,1x,i0,a,i0,a,i0)')
     +    'Occupied MOs (nocc/alpha/beta)   :',nocc,' / ',nalf,' / ',nb
      else
        write(*,'(2x,a,1x,i0,a,i0)')
     +    'Occupied MOs (alpha/beta)        :',nalf,' / ',nb
      endif
      write(*,'(2x,a,1x,f20.10)') 'SCF/DFT energy (au)              :',
     +  escf

      end

!! ***** !!

!! ********************************************************************** !!
!! subroutine: dm1input                                                   !!
!! purpose: reads the spin-orbital 1-RDM (dm1) for a post-HF wavefunction !!
!! -- either the full-space matrix from a DMN-code file (unit 11,         !!
!! unformatted), or the active-space-only matrix from PySCF (unit 11,     !!
!! formatted, inactive block assumed diagonally doubly occupied) --       !!
!! then reconstructs the natural orbitals (occ_no/c_no) and overwrites    !!
!! the AO-basis P/Pa/Pb/Ps density matrices from it, replacing whatever   !!
!! the .fchk itself provided. The test LiH-35-CAS22 covers the pySCF      !!
!! path; the ORCA path has no test.                                       !!
!! arguments:                                                             !!
!!   dm1 (out) -- spin-orbital 1-RDM, (nspinorb,nspinorb)                 !!
!! author:                                                                !!
!! ********************************************************************** !!
      subroutine dm1input(dm1)
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      character*80 line
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      dimension  dm1(nspinorb,nspinorb)

      allocatable :: scr(:,:)

      iorca    = Iopt(43)
      ipyscf   = Iopt(87)
      thres    = 1.0d-8

      allocate(scr(igr,igr))

!! nspinorb here is ncore + nact spin-orbitals (2x ncore + nact spatial MOs). !!

      do i=1,nspinorb
        do j=1,nspinorb
          dm1(i,j)=ZERO
        end do
      end do

!! full-space dm1, as provided by the DMN code. !!
      if(iorca.ne.1.and.ipyscf.ne.1) then
        rewind(11)
        do while(.true.)
          read(11,end=99) i,j,vvv
          if(ABS(vvv).gt.thres) then
            dm1(i,j)=dm1(i,j)+vvv
            dm1(j,i)=dm1(i,j)
          end if
        end do
      else

!! active-space-only dm1, as provided by PySCF -- inactive spin-orbitals !!
!! assumed doubly occupied/diagonal, active-orbital indices offset by    !!
!! ncore to land in the full spin-orbital space.                        !!
        ncore=(nspinorb-ncasorb*2)
        write(*,'(2x,a,1x,i0,a,i0,a,i0)') 'Spin-orbitals (total/active/core):',nspinorb,' / ',ncasorb,' / ',ncore
        write(*,'(2x,a,1x,i0)') 'Generating dm1 for inactive spin-orbitals:',ncore
        write(*,*) 'Assuming diagonal dm1'

        do ii=1,ncore
          dm1(ii,ii)=ONE
        end do
        rewind(11)
        read(11,*) line
        do while(.true.)
          read(11,*,end=99) ii,jj,vvv
          if(ABS(vvv).gt.thres) then
            i=ii+ncore
            j=jj+ncore
            dm1(i,j)=vvv
            dm1(j,i)=dm1(i,j)
          end if
        end do
      end if
99    write(*,*) "End DM1 reading "

!! sanity check: dm1 should be idempotent (dm1^2 = dm1) for a single- !!
!! determinant-equivalent density; reported, not enforced. !!
      difi=ZERO
      do i=1,nspinorb
        do j=1,nspinorb
          x=ZERO
          do k=1,nspinorb
            x=x+dm1(i,k)*dm1(k,j)
          end do
          difi=dmax1(difi,dabs(x-dm1(i,j)))
        end do
      end do
      write(*,*) " "
      write(*,'(2x,a,1x,f12.8)') 'Max. deviation from idempotency:',difi
      write(*,*) " "

      do i=1,norb
        do k=1,norb
          occ_no(i,k)=dm1((i-1)*2+1,(k-1)*2+1)+dm1((i-1)*2+2,(k-1)*2+2)
        end do
      end do

!! diagonalizing to get NOs and transform to the AO basis. !!
      call diagonalize(igr,norb,occ_no,scr,0)
      do i=1,igr
        do j=1,norb
          x=ZERO
          do k=1,norb
            x=x+c(i,k)*scr(k,j)
          end do
          c_no(i,j)=x
        end do
      end do

!! building Pa and Pb from dm1. !!
      write(*,*) 'Reconstructing P and Ps matrices from dm1 file'
      do nu=1,igr
        do mu=1,igr
          x=0.d0
          do i=1,nspinorb,2
            ii=i/2+1
            do j=1,nspinorb,2
              jj=j/2+1
              x=x+c(nu,ii)*dm1(i,j)*c(mu,jj)
            end do
          end do
          pa(mu,nu)=x
        end do
      end do

      do nu=1,igr
        do mu=1,igr
          x=0.d0
          do i=2,nspinorb,2
            ii=i/2
            do j=2,nspinorb,2
              jj=j/2
              x=x+c(nu,ii)*dm1(i,j)*c(mu,jj)
            end do
          end do
          pb(mu,nu)=x
        end do
      end do

      do mu=1,igr
        do nu=1,igr
          ps(mu,nu)=pa(mu,nu)-pb(mu,nu)
        end do
      end do

      dif=0.0d0
      do i=1,igr
        do j=1,igr
          dif=dmax1(dif,dabs((pa(i,j)+pb(i,j))-p(i,j)))
          p(i,j)=pa(i,j)+pb(i,j)
        end do
      end do
      write(*,*) 'P (and Ps) in fchk overridden'
      write(*,'(2x,a,1x,f12.8)') 'Max. deviation, P from dm1 vs. from .fchk:',dif

      deallocate(scr)
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: dm2input_pyscf                                            !!
!! purpose: reads the active-space-only spinless 2-RDM (dm2) produced by !!
!! PySCF's write_dm12 (unit 12, formatted *.dm2), then reconstructs the  !!
!! core and core-active blocks analytically from the spin-resolved 1-RDM !!
!! (dm1, already read by dm1input) assuming a closed-shell inactive      !!
!! space. Covered by the test LiH-35-CAS22.                              !!
!! arguments:                                                             !!
!!   dm1 (in)  -- spin-orbital 1-RDM, (nspinorb,nspinorb)                 !!
!!   dm2 (out) -- spinless active+core 2-RDM, (norb,norb,norb,norb)       !!
!! author:                                                                 !!
!! ********************************************************************* !!
      subroutine dm2input_pyscf(dm1,dm2)
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      dimension dm2(norb,norb,norb,norb)
      dimension dm1(nspinorb,nspinorb)
      character*80 line
      real*8 vvv
      logical ilog

      dm2=0.0d0

      ninact=(nspinorb-ncasorb*2)
      ncore=ninact/2
      iorb=ncore+ncasorb
      if(iorb.ne.norb) call apost_stop('inconsistency in norb')
      write(*,'(2x,a,1x,i0)') 'Reconstructing rdm2 for inactive orbitals:',ncore
      write(*,'(2x,a,1x,i0)') 'Active orbitals:',ncasorb
      write(*,'(2x,a,1x,i0)') 'Total spinless rdm2 dimension:',norb

!! active space, from *.dm2. !!
      rewind(12)
      read(12,*)
      icount=0
      do while(.true.)
        read(12,*,end=999) ii,jj,kk,ll,vvv
        dm2(ii+ncore,kk+ncore,jj+ncore,ll+ncore)=vvv
        icount=icount+1
      end do
999   continue
      write(*,'(2x,a,1x,i0)') 'rdm2 elements read:',icount

!! core and core-active blocks, reconstructed from the spin-resolved dm1. !!

!! alpha-alpha and beta-beta cases. !!
      do ispin=1,2
        do ii=1,ncore
          do jj=1,ncore
            dm2(ii,ii,jj,jj)=dm2(ii,ii,jj,jj)+ONE
            dm2(ii,jj,jj,ii)=dm2(ii,jj,jj,ii)-ONE
          end do
          do kk=ncore+1,iorb
            do ll=ncore+1,iorb
              dm2(ii,ii,kk,ll)=dm2(ii,ii,kk,ll)+dm1(2*(kk-1)+ispin,2*(ll-1)+ispin)
              dm2(kk,ll,ii,ii)=dm2(kk,ll,ii,ii)+dm1(2*(kk-1)+ispin,2*(ll-1)+ispin)
              dm2(ii,ll,kk,ii)=dm2(ii,ll,kk,ii)-dm1(2*(kk-1)+ispin,2*(ll-1)+ispin)
              dm2(kk,ii,ii,ll)=dm2(kk,ii,ii,ll)-dm1(2*(kk-1)+ispin,2*(ll-1)+ispin)
            end do
          end do
        end do
      end do

!! alpha-beta case. !!
      do ii=1,ncore
        dm2(ii,ii,ii,ii)=dm2(ii,ii,ii,ii)+TWO
        do jj=1,ncore
          if(jj.ne.ii) dm2(ii,ii,jj,jj)=dm2(ii,ii,jj,jj)+TWO
        end do
        do kk=ncore+1,iorb
          do ll=ncore+1,iorb
            dm2(ii,ii,kk,ll)=dm2(ii,ii,kk,ll)+dm1(2*(kk-1)+1,2*(ll-1)+1)
            dm2(kk,ll,ii,ii)=dm2(kk,ll,ii,ii)+dm1(2*(kk-1)+2,2*(ll-1)+2)
            dm2(ii,ii,ll,kk)=dm2(ii,ii,ll,kk)+dm1(2*(ll-1)+2,2*(kk-1)+2)
            dm2(ll,kk,ii,ii)=dm2(ll,kk,ii,ii)+dm1(2*(ll-1)+1,2*(kk-1)+1)
          end do
        end do
      end do

!! checking dm2: traces should reproduce N(N-1) when normalized. !!
      xx2=ZERO
      do ii=1,ncore
        do jj=1,ncore
          xx2=xx2+dm2(ii,ii,jj,jj)
        end do
      end do
      write(*,'(2x,a,1x,f14.8)') 'Trace of inactive orbitals:',xx2
      xx2=ZERO
      do ii=ncore+1,norb
        do jj=ncore+1,norb
          xx2=xx2+dm2(ii,ii,jj,jj)
        end do
      end do
      write(*,'(2x,a,1x,f14.8)') 'Trace of active orbitals:',xx2
      xx2=ZERO
      do ii=1,norb
        do jj=1,norb
          xx2=xx2+dm2(ii,ii,jj,jj)
        end do
      end do
      nelect=nalf+nb
      write(*,*) " "
      write(*,'(2x,a,1x,f14.8,a,1x,i0)') 'Trace of the DM2, normalized to N(N-1):',xx2,' / ',nelect*(nelect-1)
      write(*,*) " "
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: dm2input_dmn                                              !!
!! purpose: reads the full spin-separated 2-RDM produced by the DMN code !!
!! (unit 12, unformatted) -- only the i<=j, k<=l, i<=k elements are      !!
!! stored on file, reconstructed here into the full spinless dm2(norb^4) !!
!! via the AAAA/ABAB/BABA/BBBB and ABBA/BAAB symmetry relations.         !!
!! arguments:                                                            !!
!!   dm1 (in)  -- spin-orbital 1-RDM, unused here (kept for a uniform    !!
!!               call signature with dm2input_pyscf)                     !!
!!   dm2 (out) -- spinless 2-RDM, (norb,norb,norb,norb)                  !!
!! author:                                                               !!
!! ********************************************************************* !!
      subroutine dm2input_dmn(dm1,dm2)
      use ao_matrices
      implicit real*8(a-h,o-z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /cas/icas,ncasel,ncasorb,nspinorb,norb,icisd,icass
      dimension dm2(norb,norb,norb,norb)
      dimension dm1(nspinorb,nspinorb)
      character*80 line
      logical ilog

      thres = 1.0d-8

      dm2=ZERO

!! reading the dm2(1,2,1,2) matrix but saving as dm2(1,1,2,2) (important !!
!! for downstream consumers). Reads only i<=j, k<=l, i<=k elements from  !!
!! the formatted (unlabeled) *.dm2 file and reconstructs the full rdm2   !!
!! via symmetry -- ij0 for the AAAA/ABAB/BABA/BBBB terms, kj0 for the    !!
!! ABBA/BAAB terms.                                                      !!
      rewind(12)
      do while(.true.)
        read(12,end=999) i0,k0,j0,l0,vvv
        if(abs(vvv).gt.thres) then
          ii=(i0-1)/2+1
          kk=(k0-1)/2+1
          jj=(j0-1)/2+1
          ll=(l0-1)/2+1
          ilog=.true.
          if((i0.eq.j0.and.k0.eq.l0).or.(k0.eq.j0.and.i0.eq.l0)) ilog=.false.

          ij0=mod(i0+j0,2)+mod(k0+l0,2)
          kj0=mod(k0+j0,2)+mod(i0+l0,2)
          if(ij0.eq.0)then
            dm2(ii,jj,kk,ll)=dm2(ii,jj,kk,ll)+vvv
            if(ilog) dm2(jj,ii,ll,kk)=dm2(jj,ii,ll,kk)+vvv
            if(k0.ne.i0.and.l0.ne.j0) then
              dm2(kk,ll,ii,jj)=dm2(kk,ll,ii,jj)+vvv
              if(ilog) dm2(ll,kk,jj,ii)=dm2(ll,kk,jj,ii)+vvv
            end if
          end if
          if(kj0.eq.0) then
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
      end do

999   continue

!! checking dm2: trace should reproduce N(N-1) when normalized. !!
      xx2=ZERO
      do ii=1,norb
        do jj=1,norb
          xx2=xx2+dm2(ii,ii,jj,jj)
        end do
      end do
      write(*,'(2x,a,1x,f14.8)') 'Trace of the DM2, normalized to N(N-1):',xx2
      write(*,*) " "
      end

!! ***** !!

!! ********************************************************************** !!
!! subroutine: mga_misc                                                   !!
!! purpose: reads the reference one-/two-electron energy components       !!
!! (kinetic, electron-nuclear, electron-electron) that ENPART's           !!
!! integration-error checks compare against, from custom fields appended  !!
!! to the .fchk by the project's own post-processing scripts (Gaussian:   !!
!! utils/get_energy(_g16); ORCA: orca2fchk) -- not a stock Gaussian/ORCA  !!
!! field; pySCF: utils/apost3d.py. Also reads an ECP matrix if present,   !!
!! returning per-atom (Mulliken-type) ECP energies; the reference E-N     !!
!! then includes the ECP energy (ENPART adds eecp to the E-N terms).      !!
!! arguments:                                                             !!
!!   iecp (out) -- 1 if an ECP matrix was found in the .fchk, 0 otherwise !!
!!   eecp (out) -- per-atom ECP energy (only filled if iecp=1)            !!
!! author: MGimf                                                          !!
!! ********************************************************************** !!
      subroutine mga_misc(iecp,eecp)
      use basis_set
      use ao_matrices
      use integration_grid
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common /nat/ nat,igr,ifg,nocc,nalf,nb,kop
      common /iops/iopt(200)
      common /achi/achi(maxat,maxat),ibcp
      common/energ/escf,eelnuc,ekinen,erep,coulen,exchen,exchen_hf,etot
      common/energ0/ekin0,eelnuc0,evee0,eone0
      dimension eecp(maxat)
      character*80 line

      allocatable :: xecpv(:),xecpm(:,:)

      ekin0=ZERO
      eelnuc0=ZERO
      evee0=ZERO

      rewind(15)
991   read(15,'(a80)',end=891) line
      if(index(line,"Kinetic Energy").ne.0) then
        read(line(45:71),'(es27.15)')  ekin0
      else
        go to 991
      end if
891   continue
      rewind(15)
992   read(15,'(a80)',end=892) line
      if(index(line,"Electron-Nuclei Energy").ne.0) then
        read(line(45:71),'(es27.15)')  eelnuc0
      else
        go to 992
      end if
892   continue
      rewind(15)
993   read(15,'(a80)',end=893) line
      if(index(line,"Electron-Electron Energy").ne.0) then
        read(line(45:71),'(es27.15)')  evee0
      else
        go to 993
      end if
893   continue

!! pseudopotential (ECP) matrix, if present: header "ECP Matrix ... N=" !!
!! then the lower triangle of the AO matrix, as written by the utilities.  !!
      iecp=0
      rewind(15)
998   read(15,'(a80)',end=897) line
      if(index(line,"ECP Mat").eq.0) go to 998
      iecp=1
897   continue

      if(iecp.eq.1) then
        write(*,'(2x,a)') 'Pseudopotential matrix found in .fchk file'
        ntri=igr*(igr+1)/2
        nn=-1
        i=index(line,"N=")
        if(i.gt.0) read(line(i+2:),*,iostat=ios) nn
        if(nn.ne.ntri) then
          write(*,'(2x,a,i0,a,i0,a)') 'ECP Matrix: ',nn,
     +      ' values announced, ',ntri,' expected (lower triangle)'
          call apost_stop(' Wrong ECP Matrix size in the .fchk file')
        end if
        ALLOCATE(xecpv(ntri))
        read(15,*,iostat=ios) (xecpv(i),i=1,ntri)
        if(ios.ne.0) call apost_stop(' ECP Matrix in the .fchk file is incomplete')

!! lower-triangular packed -> full symmetric matrix. !!
        ALLOCATE(xecpm(igr,igr))
        k=1
        do i=1,igr
          do j=1,i
            xecpm(i,j)=xecpv(k)
            xecpm(j,i)=xecpm(i,j)
            k=k+1
          end do
        end do
        DEALLOCATE(xecpv)

!! per-atom ECP energy, Mulliken-type (P contracted against the ECP matrix). !!
        xecp=0.0d0
        do iiat=1,nat
          xxx=0.0d0
          do i=llim(iiat),iulim(iiat)
            do k=1,igr
              xxx=xxx+p(i,k)*xecpm(k,i)
            end do
          end do
          eecp(iiat)=xxx
          xecp=xecp+xxx
        end do
        DEALLOCATE(xecpm)
        write(*,'(2x,a,1x,f16.8)') 'Pseudopotential energy (au)      :',
     +    xecp
      end if

      return
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: field_misc                                                !!
!! purpose: scans the .fchk for an external static electric field entry  !!
!! and reads its 4 components; ifield=1 if the x/y/z components are      !!
!! non-negligible (the 4th component, magnitude, is read but not used    !!
!! for this check).                                                      !!
!! arguments:                                                            !!
!!   ifield (out) -- 1 if a non-negligible field is present, 0 otherwise !!
!! author:                                                               !!
!! ********************************************************************* !!
      subroutine field_misc(ifield)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'parameter.h'
      common/efield/field(4),edipole
      character*80 line

      ifield=0

      rewind(15)
1     read(15,'(a80)',end=2) line
      if(index(line,"External E-field").ne.0) then
        read(15,*) (field(i),i=1,4)
        xx=abs(field(2))+abs(field(3))+abs(field(4))
        if(xx.gt.1.0d-4) ifield=1
      else
        go to 1
      end if
2     continue
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: readchar                                                  !!
!! purpose: scans an .inp section (unit 16, already open) for a bare     !!
!! keyword; ival=1 if found before the section closes, 0 otherwise.      !!
!! arguments:                                                            !!
!!   section (in)  -- .inp section header, e.g. "# METHOD"               !!
!!   keyword (in)  -- keyword text to search for                         !!
!!   ival    (out) -- 1 if found, 0 otherwise                            !!
!! author:                                                               !!
!! ********************************************************************* !!
      subroutine readchar(section,keyword,ival)
      character section*(*), keyword*(*)
      character linea*80
      integer ipos,ii

      ival=0
      call locate_block(16,section,ii)

!! a missing block: the keyword is absent (the caller stops if the     !!
!! block is required); the end of the file also ends the block         !!
      if(ii.eq.0) return
      ii=0
      do while(ii.eq.0)
        read(16,"(a80)",end=20) linea
        call keyword_value(linea,keyword,ipos)
        if(ipos.gt.0) then
          ival=1
          ii=1
        end if
        if(index(linea,"#").ne.0) ii=2
      end do
20    return
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: readreal                                                  !!
!! purpose: scans an .inp section (unit 16, already open) for a keyword  !!
!! followed by a real value; falls back to intdef if not found.          !!
!! arguments:                                                            !!
!!   section (in)  -- .inp section header, e.g. "# METHOD"               !!
!!   keyword (in)  -- keyword text to search for                         !!
!!   intv    (out) -- parsed value, or intdef if not found               !!
!!   intdef  (in)  -- default value                                      !!
!!   ilog    (in)  -- unused, kept for the call sites                    !!
!! author:                                                               !!
!! ********************************************************************* !!
      subroutine readreal(section,keyword,intv,intdef,ilog)
      character section*(*), keyword*(*)
      character linea*80
      real*8 intv,intdef,xval
      integer ilog,ios

      intv=intdef
      call locate_block(16,section,ii)
      if(ii.eq.0) return
      ii=0
      do while(ii.eq.0)
        read(16,"(a80)",end=10) linea
        call keyword_value(linea,keyword,ipos)
        if(ipos.gt.0) then
          if(linea(ipos:).ne.' ') then
            read(linea(ipos:),*,iostat=ios) xval
            if(ios.ne.0) call bad_value(section,keyword,linea)
            ii=1
          end if
        end if
        if(index(linea,"#").ne.0) ii=2
      end do
      if(ii.eq.1) intv=xval
10    return
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: readint                                                   !!
!! purpose: scans an .inp section (unit 16, already open) for a keyword  !!
!! followed by an integer value; falls back to intdef if not found.      !!
!! arguments:                                                            !!
!!   section (in)  -- .inp section header, e.g. "# METHOD"               !!
!!   keyword (in)  -- keyword text to search for                         !!
!!   intv    (out) -- parsed value, or intdef if not found               !!
!!   intdef  (in)  -- default value                                      !!
!!   ilog    (in)  -- unused, kept for the call sites                    !!
!! author:                                                               !!
!! ********************************************************************* !!
      subroutine readint(section,keyword,intv,intdef,ilog)
      character section*(*), keyword*(*)
      character linea*80
      integer ipos,ii,intv,intdef,ilog,ival,ios

      intv=intdef
      call locate_block(16,section,ii)
      if(ii.eq.0) return
      ii=0
      do while(ii.eq.0)
        read(16,"(a80)",end=40) linea
        call keyword_value(linea,keyword,ipos)
        if(ipos.gt.0) then
          if(linea(ipos:).ne.' ') then
            read(linea(ipos:),*,iostat=ios) ival
            if(ios.ne.0) call bad_value(section,keyword,linea)
            ii=1
          end if
        end if
        if(index(linea,"#").ne.0) ii=2
      end do
      if(ii.eq.1) intv=ival
40    return
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: readintlist                                               !!
!! purpose: like readint, for a keyword followed by one or more integers !!
!!   on its line; nval=0 if the keyword is absent or has no value.       !!
!! arguments:                                                            !!
!!   section (in)  -- .inp section header, e.g. "# OSLO"                 !!
!!   keyword (in)  -- keyword text to search for                         !!
!!   nmax    (in)  -- largest number of values accepted                  !!
!!   nval    (out) -- number of values read                              !!
!!   ivals   (out) -- (nmax) the values                                  !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine readintlist(section,keyword,nmax,nval,ivals)
      character section*(*), keyword*(*)
      character linea*80
      integer nmax,nval,ivals(nmax),ipos,ii,i,j,ios

      nval=0
      call locate_block(16,section,ii)
      if(ii.eq.0) return
      ii=0
      do while(ii.eq.0)
        read(16,"(a80)",end=40) linea
        call keyword_value(linea,keyword,ipos)
        if(ipos.gt.0) then
          i=ipos
          do while(i.le.len(linea))
            if(linea(i:i).eq.' ') then
              i=i+1
            else
              j=index(linea(i:),' ')
              if(j.eq.0) then
                j=len(linea)
              else
                j=i+j-2
              end if
              if(nval.eq.nmax) then
                write(*,'(2x,a,i0,a,a)') 'At most ',nmax,
     +            ' values are accepted for ',trim(keyword)
                call apost_stop(' Too many values for an input keyword')
              end if
              nval=nval+1
              read(linea(i:j),*,iostat=ios) ivals(nval)
              if(ios.ne.0) call bad_value(section,keyword,linea)
              i=j+1
            end if
          end do
          ii=1
        end if
        if(index(linea,"#").ne.0) ii=2
      end do
40    return
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: keyword_value                                             !!
!! purpose: finds keyword in an .inp line as a whole word -- at the      !!
!!   start of the line or after a blank, and followed by a blank, "=" or !!
!!   the end of the line ("-" and "_" belong to the word) -- and returns !!
!!   where its value starts, after blanks and one optional "=".          !!
!! arguments:                                                            !!
!!   linea   (in)  -- the .inp line                                      !!
!!   keyword (in)  -- keyword text to search for                         !!
!!   ipos    (out) -- start of the value (beyond the line if none), or   !!
!!                    0 if the keyword is not there                      !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine keyword_value(linea,keyword,ipos)
      character linea*(*), keyword*(*)
      integer ipos,i,j,n,lk
      logical lword

      n=len(linea)
      lk=len_trim(keyword)
      ipos=0
      if(lk.eq.0) return
      j=0
      lword=.false.
      do while(.not.lword)
        i=index(linea(j+1:),keyword(1:lk))
        if(i.eq.0) return
        i=i+j
        j=i
        lword=.true.
        if(i.gt.1) then
          if(linea(i-1:i-1).ne.' ') lword=.false.
        end if
        i=i+lk
        if(i.le.n) then
          if(linea(i:i).ne.' '.and.linea(i:i).ne.'=') lword=.false.
        end if
      end do
      do while(i.le.n)
        if(linea(i:i).ne.' ') exit
        i=i+1
      end do
      if(i.le.n) then
        if(linea(i:i).eq.'=') i=i+1
      end if
      do while(i.le.n)
        if(linea(i:i).ne.' ') exit
        i=i+1
      end do
      ipos=i

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: bad_value                                                 !!
!! purpose: stops the run when a keyword's value cannot be read, showing !!
!!   the block and the line.                                             !!
!! arguments:                                                            !!
!!   section, keyword, linea (in) -- block header, keyword, .inp line    !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine bad_value(section,keyword,linea)
      character section*(*), keyword*(*), linea*(*)

      write(*,'(2x,a,a,a,a,a)') 'Cannot read the value of ',
     +  trim(keyword),' in ',trim(section),':'
      write(*,'(4x,a)') trim(linea)
      call apost_stop(' Wrong keyword value in the input')

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: locate                                                    !!
!! purpose: rewinds iunit and scans it for a line containing string,     !!
!! leaving the file positioned just past that line for the caller to     !!
!! keep reading from. Used throughout input()/read_input() to find       !!
!! .fchk/.inp section markers before parsing what follows.               !!
!! arguments:                                                            !!
!!   iunit  (in)  -- file unit to scan (already open)                    !!
!!   string (in)  -- text to search for                                  !!
!!   ii     (out) -- 1 if found, 0 otherwise                             !!
!! author:                                                               !!
!! ********************************************************************* !!
      subroutine locate(iunit,string,ii)
      integer iunit
      character string*(*)
      character*80 linia

      rewind(iunit)

      ii=0
      do while(ii.eq.0)
        read(iunit,'(a80)',end=10) linia
        if(index(linia,string).ne.0) then
          ii=1
          return
        end if
      end do

10    write(*,'(2x,a,1x,a)') trim(string),'section not found'
      return
      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: locate_block                                              !!
!! purpose: like locate, for an .inp block header: the line must be the  !!
!!   header itself (blanks ignored), so "# DFTDM1" does not find         !!
!!   "# DFTDM1_FUNCTIONAL". Leaves the file just past the header, and    !!
!!   reports a missing block.                                            !!
!! arguments:                                                            !!
!!   iunit  (in)  -- file unit to scan (already open)                    !!
!!   string (in)  -- block header, e.g. "# METHOD"                       !!
!!   ii     (out) -- 1 if found, 0 otherwise                             !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine locate_block(iunit,string,ii)
      integer iunit,ii
      character string*(*)

      call find_block(iunit,string,ii)
      if(ii.eq.0) write(*,'(2x,a,1x,a)') trim(string),'section not found'

      end

!! ***** !!

!! ********************************************************************* !!
!! subroutine: find_block                                                !!
!! purpose: silent core of locate_block (see there).                     !!
!! arguments: as locate_block                                            !!
!! author: MGimf                                                         !!
!! ********************************************************************* !!
      subroutine find_block(iunit,string,ii)
      integer iunit,ii,i,k
      character string*(*)
      character*80 linia,a,b

      b=' '
      k=0
      do i=1,len(string)
        if(string(i:i).ne.' ') then
          k=k+1
          b(k:k)=string(i:i)
        end if
      end do

      rewind(iunit)
      ii=0
      do while(ii.eq.0)
        read(iunit,'(a80)',end=10) linia
        a=' '
        k=0
        do i=1,80
          if(linia(i:i).ne.' '.and.linia(i:i).ne.char(9)) then
            k=k+1
            a(k:k)=linia(i:i)
          end if
        end do
        if(a.eq.b) ii=1
      end do
10    return
      end

!! ***** !!

!! ********************************************************************* !!
!! function: int_locate                                                  !!
!! purpose: rewinds iunit and scans it for a line containing text, then  !!
!! reads an integer value from its fixed .fchk column position. Returns  !!
!! 0 and ilog=.false. if text isn't found or the value can't be parsed.  !!
!! See real_locate for the real*8 twin.                                  !!
!! arguments:                                                            !!
!!   iunit (in)  -- file unit to scan (already open)                     !!
!!   text  (in)  -- text to search for                                   !!
!!   ilog  (out) -- .true. if found and parsed, .false. otherwise        !!
!! author:                                                               !!
!! ********************************************************************* !!
      function int_locate(iunit,text,ilog)
      implicit double precision(a-h,o-z)
      character*(*) text
      character*80 line
      integer int_locate
      logical ilog

      ilog=.false.
      rewind(iunit)
1     read(iunit,'(a80)',end=99) line
      if(index(line,text).eq.0) then
        goto 1
      else
        read(line(50:61),*,err=99) int_locate
        ilog=.true.
        go to 2
      end if
99    int_locate=0
2     continue
      end

!! ***** !!

!! ********************************************************************* !!
!! function: real_locate                                                 !!
!! purpose: rewinds iunit and scans it for a line containing text, then  !!
!! reads a real*8 value from its fixed .fchk column position. Returns 0  !!
!! and ilog=.false. if text isn't found or the value can't be parsed.    !!
!! See int_locate for the integer twin.                                  !!
!! arguments:                                                            !!
!!   iunit (in)  -- file unit to scan (already open)                     !!
!!   text  (in)  -- text to search for                                   !!
!!   ilog  (out) -- .true. if found and parsed, .false. otherwise        !!
!! author:                                                               !!
!! ********************************************************************* !!
      function real_locate(iunit,text,ilog)
      implicit double precision(a-h,o-z)
      character*(*) text
      character*80 line
      real*8  real_locate
      logical ilog

      ilog=.false.
      rewind(iunit)
1     read(iunit,'(a80)',end=99) line
      if(index(line,text).eq.0) then
        goto 1
      else
        read(line(50:71),*,err=99) real_locate
        ilog=.true.
        go to 2
      end if
99    real_locate=0
2     continue
      end
