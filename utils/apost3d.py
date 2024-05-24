#####################
## init write_fchk ##
#####################
from pyscf import gto, scf, tools, fci,dft,cc
import numpy
def write_fchk(mol, mf, titol,matriu_overlap,unrest=None,myhf=None):

    if myhf is None:
     mo_coeff=mf.mo_coeff
    else:
     mo_coeff=myhf.mo_coeff
    if unrest is None:
     unrest=0

# Variables #
    label=str(type(mf))
    rohf=cas=ksdft=fullci=ccsd=0
    if label.find("UHF") != -1 or label.find("UKS") != -1: unrest=1
    if label.find("ROHF") != -1 or label.find("ROKS") != -1: rohf=1
    if label.find("CASSCF") != -1: cas=1  
    if label.find("KS") != -1  : ksdft=1
    if label.find("dftd3") != -1  : ksdft=1
    if label.find("FCI") != -1  : fullci=1
    if label.find("CCSD") != -1  : ccsd=1

#Getting number of alpha and beta electrons
    tot_elect = mol.tot_electrons()
    nalpha = (tot_elect + mol.spin) // 2
    nbeta = nalpha - mol.spin
    natoms = len(mol._atom)
    spinqn =  mol.ms 

# Getting P matrix. Spin resolved if available
    if cas == 1 : 
      dm = mf.make_rdm1s() 
    elif fullci == 1: # in this cas it is generated in MO basis, trasnforming
      dma = mo_coeff.dot(mf.make_rdm1s(mf.ci,mf.norb,mf.nelec)[0]).dot(mo_coeff.T)
      dmb = mo_coeff.dot(mf.make_rdm1s(mf.ci,mf.norb,mf.nelec)[1]).dot(mo_coeff.T)
      dm=[dma,dmb]
    elif ccsd == 1: # in this cas it is generated in MO basis, trasnforming
      dm1 = mf.make_rdm1() 
      dm1_ao = numpy.einsum('pi,ij,qj->pq', mf.mo_coeff, dm1, mf.mo_coeff.conj())
      dm=[dm1_ao/2,dm1_ao/2]
    else: 
      dm = mf.make_rdm1() 

# Getting energies and other stuff
    energy = mf.e_tot
    repul = mol.energy_nuc()	
    if cas == 1:
      nelecas=mf.nelecas #  tuple of alpha and beta electrons
      ncore=mf.ncore #  number of inactive orbitals
      ncas=mf.ncas #  number of orbitals in active space
      kin=numpy.trace(numpy.matmul(mol.intor_symmetric('int1e_kin'),dm[0]+dm[1]))
      elnuc=numpy.trace(numpy.matmul(mol.intor_symmetric('int1e_nuc'),dm[0]+dm[1]))
      coul=0.5*numpy.trace(numpy.matmul(dm[0]+dm[1],scf.hf.SCF.get_j(mf,mol,dm[0]+dm[1])))
      exch_hf=-0.5*numpy.trace(numpy.matmul(dm[0],scf.hf.SCF.get_k(mf,mol,dm[0])))
      exch_hf+=-0.5*numpy.trace(numpy.matmul(dm[1],scf.hf.SCF.get_k(mf,mol,dm[1])))
      if (mf.ci.all() == 0):
         s=(nalpha - nbeta) * .5
         ss=(s,s * (s+1))
      else:
         ss=fci.spin_op.spin_square0(mf.ci,ncas,nelecas) # two floats, <S^2> and multiplicity
    elif fullci == 1 :
      kin=numpy.trace(numpy.matmul(mol.intor_symmetric('int1e_kin'),dm[0]+dm[1]))
      elnuc=numpy.trace(numpy.matmul(mol.intor_symmetric('int1e_nuc'),dm[0]+dm[1]))
      coul=0.5*numpy.trace(numpy.matmul(dm[0]+dm[1],scf.hf.SCF.get_j(myhf,mol,dm[0]+dm[1])))
      exch_hf=-0.5*numpy.trace(numpy.matmul(dm[0],scf.hf.SCF.get_k(myhf,mol,dm[0])))
      exch_hf+=-0.5*numpy.trace(numpy.matmul(dm[1],scf.hf.SCF.get_k(myhf,mol,dm[1])))
      ss=fci.spin_op.spin_square0(mf.ci,mf.norb,mf.nelec) # two floats, <S^2> and multiplicity
    elif ccsd == 1:
#      kin=elnuc=coul=exch_hf=0      
      kin=numpy.trace(numpy.matmul(mol.intor_symmetric('int1e_kin'),dm[0]+dm[1]))
      elnuc=numpy.trace(numpy.matmul(mol.intor_symmetric('int1e_nuc'),dm[0]+dm[1]))
#      coul=0.5*numpy.trace(numpy.matmul(dm[0]+dm[1],scf.hf.SCF.get_j(myhf,mol,dm[0]+dm[1])))
#      exch_hf=-0.5*numpy.trace(numpy.matmul(dm[0],scf.hf.SCF.get_k(myhf,mol,dm[0])))
#      exch_hf+=-0.5*numpy.trace(numpy.matmul(dm[1],scf.hf.SCF.get_k(myhf,mol,dm[1])))
      coul=exch_hf=0      
      ss=[0.0,1]
    else: # ok for KS-DFT and single determinant
      if (unrest == 0 and rohf == 0):
        kin=numpy.trace(numpy.matmul(mol.intor_symmetric('int1e_kin'),dm))
        elnuc=numpy.trace(numpy.matmul(mol.intor_symmetric('int1e_nuc'),dm))
        coul=0.5*numpy.trace(numpy.matmul(dm,mf.get_j()))
        exch_hf=-0.25*numpy.trace(numpy.matmul(dm,mf.get_k()))
      else:
        kin=numpy.trace(numpy.matmul(mol.intor_symmetric('int1e_kin'),dm[0]+dm[1]))
        elnuc=numpy.trace(numpy.matmul(mol.intor_symmetric('int1e_nuc'),dm[0]+dm[1]))
        coul=0.5*numpy.trace(numpy.matmul(dm[0]+dm[1],mf.get_j()[0]+mf.get_j()[1]))
        exch_hf=-0.5*numpy.trace(numpy.matmul(dm[0],mf.get_k()[0]))
        exch_hf+=-0.5*numpy.trace(numpy.matmul(dm[1],mf.get_k()[1]))
      ss=mf.spin_square()
    vee=energy -repul - kin - elnuc
    if ksdft:
      xmix=dft.libxc.hybrid_coeff(mf.xc)
      vhf = mf.get_veff(mol, dm)
      exch_dft=vhf.exc-xmix*exch_hf

# if basis functions have been removed should be somewhere. assuming none
    basis_func=mol.nao_nr()
    indep_bf=basis_func

#atomic numbers and atomic charges
    llista_num_atom = mol.atom_charges()
    llista_sym_atom=numpy.zeros(mol.natm)
    for ia in range(mol.natm):
      llista_sym_atom[ia]=mol.atom_nelec_core(ia)+mol.atom_charge(ia)
#    nuclear_charges = []
#    for x in range(0,len(mol.atom)):
#        llista_num_atom.append(mol.atom[x][0])
#    nuclear_charges = [str(("%10.8E"%(i))) for i in llista_num_atom]
#    llista_num_atom_str = [str(i) for i in llista_num_atom]


#shells, primitives per shell, highest angular momentum
    num_ctr_primitives = 0
    num_ctr_shells = 0
    for i in mol._bas:
        num_ctr_primitives += i[2]*i[3]
        num_ctr_shells += i[3]
    num_coord_shells = (num_ctr_shells*3)
    angular_momentum = []
    for i in mol._bas:
        for k in range(0,i[3]):
            angular_momentum.append(i[1])

# highest angular momentum
    highest_angular_mom = 0
    for i in angular_momentum:
        if i > highest_angular_mom:
            highest_angular_mom = i
    if highest_angular_mom >3 :
      print('CAN NOT HANDLE G-TYPE FUNCTIONS')
      exit()
#shells per atom
    ctr_total = mol.bas_nprim(range(0,mol.nbas))       
    highest_contraction = 0    
    for i in ctr_total:
        if i > highest_contraction:
            highest_contraction = i

    if mol.cart == False:
        for i in range(0,len(angular_momentum)):
            if angular_momentum[i] > 1:
                angular_momentum[i] = angular_momentum[i]*-1 
           

#Creating map array for reordering basis functions 
# For Cartesian basis functions
    if mol.cart == True:
        imap= list(range(int(basis_func)+1))
        for i in range(0,int(basis_func)):
            imap[i]=i

        contador = 0
        for i in mol.ao_labels():
            if i.find('dxy') != -1:
                imap[contador] += 2
            elif i.find('dxz') != -1:
                imap[contador] +=3 
            elif i.find('dyy') != -1:
                imap[contador] -= 2
            elif i.find('dyz') != -1:
                imap[contador] -= 2
            elif i.find('dzz') != -1:
                imap[contador] -= 1 
      
            elif i.find('fxxx') != -1:
                imap[contador] -= 0  
            elif i.find('fxxy') != -1:
                imap[contador] += 5
            elif i.find('fxxz') != -1:
                imap[contador] += 7
            elif i.find('fxyy') != -1:
                imap[contador] += 0
            elif i.find('fxyz') != -1:
                imap[contador] -= 3
            elif i.find('fxzz') != -1:
                imap[contador] -= 3
            elif i.find('fyyy') != -1:
                imap[contador] -= 1
            elif i.find('fyyz') != -1:
                imap[contador] += 1
            elif i.find('fyzz') != -1:
                imap[contador] -= 1
            elif i.find('fzzz') != -1:
                imap[contador] -= 5
# ara les G
            contador += 1

    else:

# For Pure basis functions
        imap= list(range(int(basis_func)+1))
        for i in range(0,int(basis_func)):
            imap[i]=i

        contador = 0
        for i in mol.ao_labels():
            if i.find('dxy') != -1:
                imap[contador] += 2
            elif i.find('dyz') != -1:
                imap[contador] += 2
            elif i.find('dz^2') != -1:
                imap[contador] -= 1
            elif i.find('dxz') != -1:
                imap[contador] += 1
            elif i.find('dx2-y2') != -1:
                imap[contador] -= 4

            elif i.find('f-3') != -1:
                imap[contador] += 3
            elif i.find('f-2') != -1:
                imap[contador] += 3
            elif  i.find('f-1') != -1:
                imap[contador] += 0
            elif i.find('f+0') != -1 or i.find('f 0') != -1:
                imap[contador] += 2
            elif i.find('f+1') != -1 or i.find('f 1') != -1 :
                imap[contador] -= 3
            elif i.find('f+2') != -1 or i.find('f 2') != -1:
                imap[contador] += 1
            elif i.find('f+3') != -1 or i.find('f 3') != -1:
                imap[contador] -= 6

# ara les G
            contador +=1

## End basis functions section

####################
## WRITTING FCHK  ##
####################

    nameHandle = open(titol+'.fchk', 'w')

    nameHandle.write('Automatically generated file for job '+titol + '\n')
    if cas==1:
      nameHandle.write('SP        CASSCF '   + mol.basis+ '\n')
    elif unrest==1:
      nameHandle.write('SP        Unrestricted'+ mol.basis+ '\n')
    elif rohf==1 :
      nameHandle.write('SP        RO calculation  '   + mol.basis+ '\n')
    elif ccsd==1 :
      nameHandle.write('SP        CCSD   '   + mol.basis+ '\n')
    else:
      nameHandle.write('SP        Restricted  '   + mol.basis+ '\n')

    print('Number of atoms'.ljust(43)+'I     ',"%11i"% (natoms),file=nameHandle)

#    nameHandle.write('Info1-9'+'\n')
    print('Charge'.ljust(43)+'I     ',"%11i"% (mol.charge),file=nameHandle)

    print('Number of electrons'.ljust(43)+'I     ',"%11i"% (tot_elect),file=nameHandle)
    print('Number of alpha electrons'.ljust(43)+'I     ',"%11i"% (nalpha),file=nameHandle)
    print('Number of beta electrons'.ljust(43)+'I     ',"%11i"% (nbeta),file=nameHandle)
    print('Number of basis functions'.ljust(43)+'I     ',"%11i"% (basis_func),file=nameHandle)
    print('Number of independent functions'.ljust(43)+'I     ',"%11i"% (indep_bf),file=nameHandle)
    
#    nameHandle.write ('Number of point charges in /Mol/' + '           ' + 'I' + '\n')
#    nameHandle.write('Number of translation vectors' + '              ' + 'I' + '\n')

    print('Atomic numbers'.ljust(43)+'I   N=',"%11i"%(natoms),file=nameHandle)
    nums=["{:12.0f}".format(i) for i in llista_sym_atom]
    print('\n'.join(''.join(nums[i:i+6]) for i in range(0, len(nums), 6)),file=nameHandle)

    print('Nuclear charges'.ljust(43)+'R   N=',"%11i"%(natoms),file=nameHandle)
    nums=["{:16.8E}".format(i) for i in llista_num_atom]
    print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle)

    print('Current cartesian coordinates'.ljust(43)+'R   N=',"%11i"%(3*natoms),file=nameHandle)
    coords=mol.atom_coords().ravel()
    nums=["{:16.8E}".format(i) for i in coords]
    print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle)

# irrelevant stuff
    nameHandle.write('Force Field' + '\n')
    nameHandle.write('Int Atom Types' + '\n')
    nameHandle.write('MM charges' + '\n')
    nameHandle.write('Atom fragment info' + '\n')
    nameHandle.write('Atom residue num' + '\n')
    nameHandle.write('Nuclear spins' + '\n')
    nameHandle.write('Nuclear QMom' + '\n')
    nameHandle.write('Nuclear GFac' + '\n')
    nameHandle.write('MicOpt' + '\n')
    nameHandle.write('SymInf integers' + '\n')
    nameHandle.write('RotTr to input orientation'+'\n')
# possibly relevant stuff
    nameHandle.write('Nuclear ZEff' + '\n')
    nameHandle.write('Nuclear ZNuc' + '\n')
    nameHandle.write('Integer atomic weights' + '\n')
    nameHandle.write('Real atomic weights' + '\n')

##
## Basis set info
##
    print('Number of contracted shells'.ljust(43)+'I     ',"%11i"% (num_ctr_shells),file=nameHandle)
    print('Number of primitive shells'.ljust(43)+'I     ',"%11i"% (num_ctr_primitives),file=nameHandle)
    val=0
    if  mol.cart == True: val=1
    print('Pure/Cartesian d shells'.ljust(43)+'I     ',"%11i"% (1),file=nameHandle)
    print('Pure/Cartesian f shells'.ljust(43)+'I     ',"%11i"% (1),file=nameHandle)
    print('Highest angular momentum'.ljust(43)+'I     ',"%11i"% (highest_angular_mom),file=nameHandle)
    print('Largest degree of contraction'.ljust(43)+'I     ',"%11i"% (highest_contraction),file=nameHandle)

    print('Shell types'.ljust(43)+'I   N=',"%11i"%(num_ctr_shells),file=nameHandle)
    nums=["{:12.0f}".format(i) for i in angular_momentum]
    print('\n'.join(''.join(nums[i:i+6]) for i in range(0, len(nums), 6)),file=nameHandle)

    val=[]
    for j in range(0,mol.nbas):
        for x in range(0,mol._bas[j][3]):
            val.append(mol._bas[j][2])
    print('Number of primitives per shell'.ljust(43)+'I   N=',"%11i"%(num_ctr_shells),file=nameHandle)
    nums=["{:12.0f}".format(i) for i in val]
    print('\n'.join(''.join(nums[i:i+6]) for i in range(0, len(nums), 6)),file=nameHandle)

    val=[]
    for j in range(0,mol.nbas):
        for x in range(0,mol._bas[j][3]):
            val.append(mol._bas[j][0]+1)
    print('Shell to atom map'.ljust(43)+'I   N=',"%11i"%(num_ctr_shells),file=nameHandle)
    nums=["{:12.0f}".format(i) for i in val]
    print('\n'.join(''.join(nums[i:i+6]) for i in range(0, len(nums), 6)),file=nameHandle)
    
    val=[]
    for j in range(0,mol.nbas):
     for x in range(0,mol._bas[j][3]):
      for i in mol.bas_exp(j):
       val.append(i)
    print('Primitive exponents'.ljust(43)+'R   N=',"%11i"%(num_ctr_primitives),file=nameHandle)
    nums=["{:16.8E}".format(i) for i in (val)]
    print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle) 

    val=[]
    for j in range(0,mol.nbas):
     for x in range(0,mol._bas[j][3]):
      for i in mol.bas_ctr_coeff(j):
       val.append(i[x])
    print('Contraction coefficients'.ljust(43)+'R   N=',"%11i"%(num_ctr_primitives),file=nameHandle)
    nums=["{:16.8E}".format(i) for i in (val)]
    print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle) 

    val=[]
    for j in range(0,mol.nbas):
     for x in range(0,mol._bas[j][3]):
      for i in mol.bas_coord(j):
       val.append(i)
    print('Coordinates of each shell'.ljust(43)+'R   N=',"%11i"%(num_coord_shells),file=nameHandle)
    nums=["{:16.8E}".format(i) for i in (val)]
    print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle) 

# some irrelevant stuff
#    print('Constraint Structure'.ljust(43)+'R   N=',"%11i"%(3*natoms),file=nameHandle)
    print('Num ILSW'.ljust(43)+'I     ',"%11i"% (0),file=nameHandle)
    print('ILSW'.ljust(43)+'I   N=',"%11i"%(0),file=nameHandle)
    print('Num RLSW'.ljust(43)+'I     ',"%11i"% (0),file=nameHandle)
    print('RLSW'.ljust(43)+'R   N=',"%11i"%(0),file=nameHandle)
    print('MxBond'.ljust(43)+'I     ',"%11i"% (0),file=nameHandle)
    print('NBond'.ljust(43)+'I   N=',"%11i"%(0),file=nameHandle)
    print('IBond'.ljust(43)+'I   N=',"%11i"%(0),file=nameHandle)
    print('RBond'.ljust(43)+'R   N=',"%11i"%(0),file=nameHandle)

# energetics and calculation data
    print('SCF Energy'.ljust(43)+'R    ',"%22.15E"% (energy),file=nameHandle)
    print('Total Energy'.ljust(43)+'R    ',"%22.15E"% (energy),file=nameHandle)
    print('S**2'.ljust(43)+'R    ',"%22.15E"% (ss[0]),file=nameHandle)
    print('Virial Ratio'.ljust(43)+'R    ',"%22.15E"% ((kin - energy)/kin),file=nameHandle)
#    print('S**2 after annihilation'.ljust(43)+'R     ',"%22.15E"% (ss[0]),file=nameHandle)
    print('RMS Density'.ljust(43)+'R    ',"%22.15E"% (0),file=nameHandle)
    print('Job Status'.ljust(43)+'I     ',"%11i"% (1),file=nameHandle)
#    print('Nuclear derivative order'.ljust(43)+'I     ',"%11i"% (0),file=nameHandle)

# to be fixed
    print('External E-field'.ljust(43)+'R   N=',"%11i"%(35),file=nameHandle)
    nums=["{:16.8E}".format(0) for i in range(35)]
    print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle) 

    print('IOpCl'.ljust(43)+'I     ',"%11i"% (0),file=nameHandle)
    print('IROHF'.ljust(43)+'I     ',"%11i"% (rohf),file=nameHandle)

# MO  and density matrix information
# Printing Spin density when available
    if cas == 1:
      print('Number of CAS Electrons'.ljust(43)+'I     ',"%11i"% (nelecas[0]+nelecas[1]),file=nameHandle)
      print('Number of CAS Orbitals'.ljust(43)+'I     ',"%11i"% (ncas),file=nameHandle)
      val=numpy.zeros(indep_bf)
    elif fullci == 1:
      print('Number of CAS Electrons'.ljust(43)+'I     ',"%11i"% (mf.nelec[0]+mf.nelec[1]),file=nameHandle)
      print('Number of CAS Orbitals'.ljust(43)+'I     ',"%11i"% (mf.norb),file=nameHandle)
      val=numpy.zeros(indep_bf)
    elif ccsd == 1:
      val=numpy.zeros(indep_bf)
    else:
       if unrest :
         val=mf.mo_energy[0]
       else:
         val=mf.mo_energy
    print('Alpha Orbital Energies'.ljust(43)+'R   N=',"%11i"%(indep_bf),file=nameHandle)
    nums=["{:16.8E}".format(val[i]) for i in range(indep_bf)]
    print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle)

    if unrest: 
      print('Beta Orbital Energies'.ljust(43)+'R   N=',"%11i"%(indep_bf),file=nameHandle)
      nums=["{:16.8E}".format(i) for i in mf.mo_energy[1]]
      print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle)

    if cas == 1 or fullci== 1 or unrest == 0:
      val=[]
      for j in range(indep_bf):
       for i in range(basis_func):
           val.append(mo_coeff[imap[i]][j]*numpy.sqrt(matriu_overlap[imap[i]][imap[i]]))
      print('Alpha MO coefficients'.ljust(43)+'R   N=',"%11i"%(basis_func*indep_bf),file=nameHandle)
      if len(val) != basis_func*indep_bf: 
       print('Wrong number of items')
       return
      nums=["{:16.8E}".format(i) for i in (val)]
      print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle) 

    if unrest:
      val=[]
      valb=[]
      for j in range(indep_bf):
       for i in range(basis_func):
           val.append(mo_coeff[0][imap[i]][j]*numpy.sqrt(matriu_overlap[imap[i]][imap[i]]))
           valb.append(mo_coeff[1][imap[i]][j]*numpy.sqrt(matriu_overlap[imap[i]][imap[i]]))
      if len(val) != basis_func*indep_bf: 
       print('Wrong number of items')
       return
      print('Alpha MO coefficients'.ljust(43)+'R   N=',"%11i"%(basis_func*indep_bf),file=nameHandle)
      nums=["{:16.8E}".format(i) for i in (val)]
      print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle) 
      print('Beta MO coefficients'.ljust(43)+'R   N=',"%11i"%(basis_func*indep_bf),file=nameHandle)
      nums=["{:16.8E}".format(i) for i in (valb)]
      print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle) 
     
# Total SCF density 
    n2=int(indep_bf*(indep_bf+1)/2)
    if (unrest == 0 and rohf == 0 and cas == 0 and fullci == 0 and ccsd == 0):
       val=[]
       for j in range(indep_bf):
         for i in range(j+1):
             overlap=numpy.sqrt(matriu_overlap[imap[i]][imap[i]])*numpy.sqrt(matriu_overlap[imap[j]][imap[j]])
             val.append(dm[imap[i]][imap[j]]*overlap)
       if len(val) != n2: 
        print('Wrong number of items')
        return
       print('Total SCF Density'.ljust(43)+'R   N=',"%11i"%(n2),file=nameHandle)
       nums=["{:16.8E}".format(i) for i in (val)]
       print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle) 
#   aaa=dm[numpy.triu_indices(dm.shape[0])]
    else:
       val=[]
       valb=[]
       for j in range(indep_bf):
         for i in range(j+1):
             overlap=numpy.sqrt(matriu_overlap[imap[i]][imap[i]])*numpy.sqrt(matriu_overlap[imap[j]][imap[j]])
             val.append(dm[0][imap[i]][imap[j]]*overlap)
             valb.append(dm[1][imap[i]][imap[j]]*overlap)
       if len(val) != n2: 
        print('Wrong number of items')
        return
       print('Total SCF Density'.ljust(43)+'R   N=',"%11i"%(n2),file=nameHandle)
       nums=["{:16.8E}".format(val[i]+valb[i]) for i in range(n2)]
       print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle) 
       print('Spin SCF Density'.ljust(43)+'R   N=',"%11i"%(n2),file=nameHandle)
       nums=["{:16.8E}".format(val[i]-valb[i]) for i in range(n2)]
       print('\n'.join(''.join(nums[i:i+5]) for i in range(0, len(nums), 5)),file=nameHandle) 

##
## Additional infor for APOST3D
##
    print('Kinetic Energy'.ljust(43)+'R    ',"%22.15E"% (kin),file=nameHandle)
    print('Electron-Nuclei Energy'.ljust(43)+'R    ',"%22.15E"% (elnuc),file=nameHandle)
    print('Electron-Electron Energy'.ljust(43)+'R    ',"%22.15E"% (vee),file=nameHandle)
#    if cas != 1 and fullci != 1:
    print('Coulomb Energy'.ljust(43)+'R    ',"%22.15E"% (coul),file=nameHandle)
    print('Exact-exchange Energy'.ljust(43)+'R    ',"%22.15E"% (exch_hf),file=nameHandle)
    if ksdft == 1:
      print('DFT-exchange Energy'.ljust(43)+'R    ',"%22.15E"% (exch_dft),file=nameHandle)
      print('Total Exchange Energy'.ljust(43)+'R    ',"%22.15E"% (vhf.exc),file=nameHandle)
      if xmix != 0 : print('% of Exact-exchange'.ljust(43)+'R    ',"%22.15E"% (xmix),file=nameHandle)

####################
## end write_fchk ##
####################

###################
## init write_dm12 ##
###################
from pyscf import mcscf,fci
#import numpy
def write_dm12(mol, mycas, titol):

# Variables #
    label=str(type(mycas))
    cas=ccsd=0
    if label.find("CASSCF") != -1: cas=1  
    if label.find("CCSD") != -1  : ccsd=1
    toler=1.0e-10

    if (cas == 1):
     ncore = mycas.ncore
     ncas = int(mycas.ncas)
     nelecas = mycas.nelecas
# Spin resolved rdm1 
     casdm1s = mycas.fcisolver.make_rdm1s(mycas.ci,mycas.ncas, mycas.nelecas)
# Spinless rdm2 of the active space
     dm2  = mycas.fcisolver.make_rdm2(mycas.ci,mycas.ncas, mycas.nelecas)
    else: 
# trying for CCSD. only rdm1 for closed-shell..should call uccsd
     ncore = 0 
     ncas = mycas.get_nmo()
     dmrg = False
     casdm1s = [0,0]
     casdm1s[0] = casdm1s[1] = mycas.make_rdm1()/2
#     casdm2s = [0,0,0]
#     casdm2s[0] = casdm2s[1] = casdm2s[2] = mycas.make_rdm2()/2

# rdm1s
    nameHandle =open(titol+'.dm1', 'w')
    if( cas == 1) :
      print('mcscf spin resolved rdm1 for ',nelecas,'electrons in',ncas,'orbitals',file=nameHandle)
    else:
      print('ccsd spin resolved rdm1 for ',ncas,'orbitals',file=nameHandle)

    for case in (0,1):   
      for i in range(ncas):
        for j in range(i,ncas):
          if abs(casdm1s[case][i,j]) >= (toler): 
            print(2*i+1+case,2*j+1+case,'{:22.15E}'.format(casdm1s[case][i,j]),file=nameHandle)
    nameHandle.close()

#
# rdm2s; printing all elements in [1212] convention
#
    if(cas ==1) :
     nameHandle =open(titol+'.dm2', 'w')
     print('spinless rdm2 for ',nelecas,'electrons in',ncas,'orbitals',file=nameHandle)
     for i in range(ncas):
       for j in range(ncas):
         for k in range(ncas):
           for l in range(ncas):
             if abs(dm2[i,k,j,l]) >= toler:
               print(i+1,j+1,k+1,l+1,'{:22.15E}'.format(dm2[i,k,j,l]),file=nameHandle)
     nameHandle.close()
    else:
      print('not yet ready for CCSD or other WF types')



###########################
## end write_dm12
###########################

