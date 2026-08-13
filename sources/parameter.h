
c main limits
c nmax was 8000 through the ifort/PGO era. On gfortran/macOS (Apple Silicon,
c confirmed on macOS 26.5.1) the resulting ~2.6GB of static BSS from the
c nmax*nmax COMMON arrays in effao.f/mulliken.f/qtaim.f (/effao/, /nao/,
c /stv/) makes apost3d/apost3d-eos fail to even launch: dyld cannot map its
c shared cache alongside such a large fixed data segment ("dyld cache '(null)'
c not loaded: syscall to map cache into shared region failed"). Confirmed by
c bisection: nmax=8000 fails to launch, nmax=2000 launches fine.
c
c Lowered to 3000 as an interim default (still ~9x more basis functions than
c any current test case needs) until those COMMON blocks are converted to
c ALLOCATABLE -- the proper fix, which removes this compile-time cap
c entirely rather than just picking a smaller magic number. Until then, a
c system requiring more than nmax basis functions will silently overrun
c these fixed-size arrays (undefined behavior, not a clean bounds-checked
c error) rather than failing loudly.
      parameter (nmax=3000)
      parameter (maxat=350) 
c other limits and equivalences 
      parameter (maxp=10000,maxg=nmax,maxc=36) 
      parameter (maxnna=20)
      parameter (maxfrag=maxat)
c subroutine cubegen3
      parameter (maxgrid=10000)
      parameter (thresh=1.0d-8)

c numerical constants
      parameter (pi=3.141592653589793d0)
      parameter (pi3=31.0062766802998d0)      
      parameter (pi52=34.98683665524973d0)   

!! Numerical parameters !!
!! Double-precision !!
      parameter (ZERO=0.0d0)
      parameter (HALF=0.5d0)
      parameter (ONE=1.0d0)
      parameter (TWO=2.0d0)
      parameter (THREE=3.0d0)
      parameter (FOUR=4.0d0)
      parameter (FIVE=5.0d0)

!! Integer !!
      parameter (IZERO=0)

!! Others !!
      parameter(angtoau=0.52917721067d0) !! From Angstroms to a.u. !!
      parameter(tokcal=627.5096d0) !! From a.u. to kcal/mol !!
