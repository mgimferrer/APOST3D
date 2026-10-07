"""
apost3d.py -- APOST-3D input files from a pySCF calculation.

    import apost3d
    apost3d.write_fchk(mol, obj, 'name')          # name.fchk
    apost3d.write_dm12(mol, mycas, 'name')        # name.dm1, name.dm2

write_fchk(mol, obj, name, overlap=None, myhf=None)
    Writes name.fchk (Gaussian formatted-checkpoint layout) with the basis,
    orbitals, densities and the reference energies ENPART needs. obj is the
    pySCF object of the calculation, after its kernel():
      - SCF: RHF, UHF, ROHF, RKS, UKS, ROKS (density fitting and other
        wrappers included);
      - CASSCF/CASCI on RHF/ROHF orbitals;
      - FCI (pass the mean-field object it started from as myhf);
      - RCCSD, closed shells (its unrelaxed 1-RDM).
    Correlated wavefunctions on UHF orbitals (UHF-based CASSCF or FCI,
    UCCSD) are refused: APOST-3D takes the alpha and beta densities of a
    correlated wavefunction from the RDMs in one set of orbitals.
    overlap (the AO overlap matrix) is optional; it is computed when not
    given. Pseudopotentials are supported. An empirical dispersion energy
    (D3/D4) is removed from the reference electron-electron energy, with a
    printed note: APOST-3D decomposes the electronic energy only.

write_dm12(mol, mycas, name)
    Writes the spin-resolved 1-RDM (name.dm1) and the spinless 2-RDM
    (name.dm2) of a CASSCF/CASCI or FCI calculation on RHF orbitals, for
    SPIN and ENPART (# DM block, pySCF). For RCCSD only name.dm1 is written.

Errors raise Apost3dError with the reason; no partial file is left behind.
Tested with pySCF 2.7 and 2.14.
"""

import numpy as np
from pyscf import scf, dft, mcscf, fci, cc


class Apost3dError(Exception):
    """A pySCF object or basis that APOST-3D cannot take."""


#####################################################################
## basis-function order: pySCF -> Gaussian/APOST-3D                 ##
#####################################################################

# Cartesian components in the .fchk (Gaussian) order, as APOST-3D reads
# them (sources/modules.f90, build_basis), as (lx, ly, lz)
_CART_ORDER = {
    0: [(0, 0, 0)],
    1: [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    2: [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)],
    3: [(3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 0), (2, 1, 0), (2, 0, 1),
        (1, 0, 2), (0, 1, 2), (0, 2, 1), (1, 1, 1)],
    4: [(0, 0, 4), (0, 1, 3), (0, 2, 2), (0, 3, 1), (0, 4, 0), (1, 0, 3),
        (1, 1, 2), (1, 2, 1), (1, 3, 0), (2, 0, 2), (2, 1, 1), (2, 2, 0),
        (3, 0, 1), (3, 1, 0), (4, 0, 0)],
}


def _pyscf_cart(l):
    """pySCF's Cartesian order: x power descending, then y descending."""
    return [(lx, ly, l - lx - ly)
            for lx in range(l, -1, -1) for ly in range(l - lx, -1, -1)]


def _shell_perm(l, cart):
    """Positions, within one pySCF shell, of the .fchk functions in order."""
    if l <= 1:
        return list(range(2 * l + 1))
    if cart:
        py = _pyscf_cart(l)
        return [py.index(c) for c in _CART_ORDER[l]]
    # pure: pySCF m = -l..l; .fchk m = 0, +1, -1, +2, -2, ...
    order = [0]
    for m in range(1, l + 1):
        order += [m, -m]
    return [m + l for m in order]


def _basis(mol):
    """Shells in .fchk form and the pySCF AO index of each .fchk function."""
    shells = []          # (type, nprim, atom, exps, coeffs, coord)
    perm = []
    ao0 = 0
    for j in range(mol.nbas):
        l = mol.bas_angular(j)
        if l > 4:
            raise Apost3dError('basis functions with l > 4 (h and higher) '
                               'are not supported by APOST-3D')
        nctr = mol.bas_nctr(j)
        ncomp = (l + 1) * (l + 2) // 2 if mol.cart else 2 * l + 1
        exps = mol.bas_exp(j)
        ctr = mol.bas_ctr_coeff(j)
        stype = l if (mol.cart or l < 2) else -l
        local = _shell_perm(l, mol.cart)
        for x in range(nctr):
            shells.append((stype, len(exps), mol.bas_atom(j) + 1, exps,
                           ctr[:, x], mol.bas_coord(j)))
            perm += [ao0 + x * ncomp + k for k in local]
        ao0 += nctr * ncomp
    if ao0 != mol.nao_nr():
        raise Apost3dError('basis bookkeeping mismatch (%d vs %d functions)'
                           % (ao0, mol.nao_nr()))
    return shells, np.array(perm)


#####################################################################
## wavefunction                                                     ##
#####################################################################

_UNRESTRICTED_CORRELATED = (
    'APOST-3D takes the alpha and beta densities of a correlated '
    'wavefunction from the # DM files, which are written in one set of '
    'orbitals; a UHF-based one has two sets, so its spin-resolved analyses '
    'would be wrong. Use an RHF/ROHF-based CASSCF or FCI (with write_dm12) '
    'for open shells')


def _ao(c, d):
    return c @ d @ c.T


def _ss_scf(mf):
    try:
        return float(mf.spin_square()[0])
    except Exception:
        s = 0.5 * mf.mol.spin
        return s * (s + 1)


def _wavefunction(mol, obj, myhf):
    """Kind, reference SCF, AO alpha/beta densities, orbitals, energies."""
    w = {'cas': None, 'rohf': 0}

    if isinstance(obj, mcscf.casci.CASBase):
        if obj.ci is None:
            raise Apost3dError('run the CASSCF/CASCI calculation first')
        if isinstance(obj, mcscf.ucasci.UCASBase):
            raise Apost3dError('UHF-based CASSCF: ' + _UNRESTRICTED_CORRELATED)
        w['kind'] = 'cas'
        w['ref'] = obj._scf
        da, db = obj.make_rdm1s()
        w['mo'] = [obj.mo_coeff]
        w['moe'] = None
        w['header'] = 'CASSCF'
        w['cas'] = (int(sum(obj.nelecas)), int(obj.ncas))
        w['ss'] = float(mcscf.addons.spin_square(obj)[0])
        w['e_tot'] = obj.e_tot

    elif isinstance(obj, fci.direct_spin1.FCIBase):
        if myhf is None:
            raise Apost3dError('FCI: pass the mean-field object it started '
                               'from, write_fchk(mol, cisolver, name, '
                               'myhf=mf)')
        if getattr(obj, 'ci', None) is None:
            raise Apost3dError('run the FCI calculation first')
        if isinstance(myhf, scf.uhf.UHF):
            raise Apost3dError('FCI on UHF orbitals: ' +
                               _UNRESTRICTED_CORRELATED)
        w['ref'] = myhf
        norb = myhf.mo_coeff.shape[1]
        nelec = obj.nelec if getattr(obj, 'nelec', None) is not None \
            else mol.nelec
        dma, dmb = obj.make_rdm1s(obj.ci, norb, nelec)
        w['kind'] = 'fci'
        w['mo'] = [myhf.mo_coeff]
        da, db = _ao(myhf.mo_coeff, dma), _ao(myhf.mo_coeff, dmb)
        w['moe'] = None
        w['header'] = 'CASSCF'
        w['cas'] = (int(sum(nelec)) if not np.isscalar(nelec)
                    else int(nelec), int(norb))
        w['ss'] = float(fci.spin_op.spin_square0(obj.ci, norb, nelec)[0])
        w['e_tot'] = obj.e_tot

    elif isinstance(obj, cc.ccsd.CCSDBase):
        if getattr(obj, 't2', None) is None:
            raise Apost3dError('run the CCSD calculation first')
        w['ref'] = obj._scf
        if isinstance(obj, cc.uccsd.UCCSD):
            raise Apost3dError('UCCSD: ' + _UNRESTRICTED_CORRELATED)
        rdm = obj.make_rdm1()
        w['kind'] = 'ccsd'
        w['mo'] = [obj.mo_coeff]
        da = db = 0.5 * _ao(obj.mo_coeff, rdm)
        w['moe'] = None
        w['header'] = 'CCSD'
        w['ss'] = _ss_scf(obj._scf)
        w['e_tot'] = obj.e_tot

    elif isinstance(obj, scf.hf.SCF):
        if isinstance(obj, scf.ghf.GHF) or not isinstance(
                obj, (scf.hf.RHF, scf.uhf.UHF)):
            raise Apost3dError('only RHF/UHF/ROHF/RKS/UKS/ROKS references '
                               'are supported (not %s)' % type(obj).__name__)
        if obj.mo_coeff is None:
            raise Apost3dError('run the SCF calculation first')
        w['ref'] = obj
        ks = isinstance(obj, dft.rks.KohnShamDFT)
        dm = obj.make_rdm1()
        if isinstance(obj, scf.uhf.UHF):
            w['kind'] = 'uhf'
            da, db = dm
            w['mo'] = list(obj.mo_coeff)
            w['moe'] = list(obj.mo_energy)
            w['header'] = 'UKS' if ks else 'UHF'
        elif isinstance(obj, scf.rohf.ROHF):
            w['kind'] = 'rohf'
            w['rohf'] = 1
            da, db = dm
            w['mo'] = [obj.mo_coeff]
            w['moe'] = [obj.mo_energy]
            w['header'] = 'ROKS' if ks else 'ROHF'
        else:
            w['kind'] = 'rhf'
            da = db = 0.5 * dm
            w['mo'] = [obj.mo_coeff]
            w['moe'] = [obj.mo_energy]
            w['header'] = 'RKS' if ks else 'RHF'
        w['ss'] = _ss_scf(obj)
        w['e_tot'] = obj.e_tot

    else:
        raise Apost3dError('unsupported pySCF object %s' % type(obj).__name__)

    w['da'], w['db'] = np.asarray(da), np.asarray(db)
    nmo = {m.shape[1] for m in w['mo']}
    if len(nmo) != 1:
        raise Apost3dError('alpha and beta orbitals of different number')
    w['nmo'] = nmo.pop()
    if w['moe'] is None:
        w['moe'] = [np.zeros(w['nmo']) for _ in w['mo']]
    return w


def _dispersion(ref):
    """Empirical dispersion energy included in the reference SCF energy."""
    e = getattr(ref, 'scf_summary', {}).get('dispersion')
    if e:
        return float(e)
    try:                   # older wrappers add it to energy_nuc
        d = float(ref.energy_nuc() - ref.mol.energy_nuc())
        if abs(d) > 1e-12:
            return d
    except Exception:
        pass
    return 0.0


def _energies(mol, w):
    """Reference energy components (ENPART) and the exchange breakdown."""
    tr = lambda a, b: float(np.einsum('ij,ji->', a, b))
    ref = w['ref']
    da, db = w['da'], w['db']
    p = da + db
    e = {}
    e['kin'] = tr(mol.intor_symmetric('int1e_kin'), p)
    e['en'] = tr(mol.intor_symmetric('int1e_nuc'), p)
    e['vecp'] = None
    if mol.has_ecp():
        # the reference E-N includes the ECP energy; APOST-3D splits the
        # ECP matrix among the atoms (ENPART)
        e['vecp'] = mol.intor_symmetric('ECPscalar')
        e['en'] += tr(e['vecp'], p)
    e['disp'] = _dispersion(ref)
    e['vee'] = (w['e_tot'] - e['disp'] - mol.energy_nuc()
                - e['kin'] - e['en'])
    e['coul'] = 0.5 * tr(p, ref.get_j(mol, p))
    e['exch'] = -0.5 * (tr(da, ref.get_k(mol, da)) +
                        tr(db, ref.get_k(mol, db)))
    e['xmix'] = None
    if w['kind'] in ('rhf', 'uhf', 'rohf') and \
            isinstance(ref, dft.rks.KohnShamDFT):
        ni = ref._numint
        e['xmix'] = float(ni.hybrid_coeff(ref.xc)) if hasattr(
            ni, 'hybrid_coeff') else float(dft.libxc.hybrid_coeff(ref.xc))
        dm = p if w['kind'] == 'rhf' else np.array([da, db])
        e['exc'] = float(ref.get_veff(mol, dm).exc)
        e['exch_dft'] = e['exc'] - e['xmix'] * e['exch']
    return e


#####################################################################
## .fchk writing                                                    ##
#####################################################################

def _i(f, label, value):
    f.write('%-43sI     %12d\n' % (label, value))


def _r(f, label, value):
    f.write('%-43sR     %22.15E\n' % (label, value))


def _iarr(f, label, values):
    values = list(values)
    f.write('%-43sI   N=%12d\n' % (label, len(values)))
    for k in range(0, len(values), 6):
        f.write(''.join('%12d' % v for v in values[k:k + 6]) + '\n')


def _rarr(f, label, values):
    values = np.ravel(values)
    f.write('%-43sR   N=%12d\n' % (label, len(values)))
    for k in range(0, len(values), 5):
        f.write(''.join('%16.8E' % v for v in values[k:k + 5]) + '\n')


def write_fchk(mol, obj, name, overlap=None, unrest=None, myhf=None):
    """Writes name.fchk for APOST-3D (see the module docstring).
    unrest is accepted for old scripts and ignored."""
    w = _wavefunction(mol, obj, myhf)
    shells, perm = _basis(mol)
    s = mol.intor_symmetric('int1e_ovlp') if overlap is None \
        else np.asarray(overlap)
    # .fchk functions are normalized; pySCF's Cartesian ones are not all
    norm = np.sqrt(np.diag(s))[perm]
    e = _energies(mol, w)

    nao = mol.nao_nr()
    tril = np.tril_indices(nao)
    def lower(m, scale):            # row-wise lower triangle, .fchk order
        mm = m[np.ix_(perm, perm)] * scale
        return mm[tril]
    pp = np.outer(norm, norm)
    p = w['da'] + w['db']
    ps = w['da'] - w['db']
    mos = [(c[perm, :] * norm[:, None]).T for c in w['mo']]

    nalpha, nbeta = mol.nelec
    ztrue = [mol.atom_charge(i) + mol.atom_nelec_core(i)
             for i in range(mol.natm)]
    lmax = max(abs(sh[0]) for sh in shells)
    basis_name = mol.basis if isinstance(mol.basis, str) else 'Gen'

    if w['kind'] in ('cas', 'fci') and w['ss'] > 0.01:
        print('apost3d.py: open-shell correlated wavefunction: APOST-3D takes '
              'its alpha and beta densities from the # DM files '
              '(write_dm12), give them for spin-resolved analyses')
    if e['disp']:
        print('apost3d.py: dispersion energy %.10f au removed from the '
              'reference electron-electron energy (APOST-3D decomposes the '
              'electronic energy only)' % e['disp'])

    with open(name + '.fchk', 'w') as f:
        f.write('Automatically generated by apost3d.py for job %s\n' % name)
        # method in columns 11-25 (APOST-3D reads its flags there), basis after
        f.write('SP        %-15s %s\n' % (w['header'], basis_name))
        _i(f, 'Number of atoms', mol.natm)
        _i(f, 'Charge', mol.charge)
        _i(f, 'Multiplicity', mol.spin + 1)
        _i(f, 'Number of electrons', nalpha + nbeta)
        _i(f, 'Number of alpha electrons', nalpha)
        _i(f, 'Number of beta electrons', nbeta)
        _i(f, 'Number of basis functions', nao)
        _i(f, 'Number of independent functions', w['nmo'])
        _iarr(f, 'Atomic numbers', ztrue)
        _rarr(f, 'Nuclear charges', mol.atom_charges())
        _rarr(f, 'Current cartesian coordinates', mol.atom_coords())
        _i(f, 'Number of contracted shells', len(shells))
        _i(f, 'Number of primitive shells', sum(sh[1] for sh in shells))
        _i(f, 'Pure/Cartesian d shells', 1 if mol.cart else 0)
        _i(f, 'Pure/Cartesian f shells', 1 if mol.cart else 0)
        _i(f, 'Highest angular momentum', lmax)
        _i(f, 'Largest degree of contraction', max(sh[1] for sh in shells))
        _iarr(f, 'Shell types', [sh[0] for sh in shells])
        _iarr(f, 'Number of primitives per shell', [sh[1] for sh in shells])
        _iarr(f, 'Shell to atom map', [sh[2] for sh in shells])
        _rarr(f, 'Primitive exponents', np.concatenate([sh[3] for sh in shells]))
        _rarr(f, 'Contraction coefficients',
              np.concatenate([sh[4] for sh in shells]))
        _rarr(f, 'Coordinates of each shell', [sh[5] for sh in shells])
        _r(f, 'SCF Energy', w['e_tot'])
        _r(f, 'Total Energy', w['e_tot'])
        _r(f, 'S**2', w['ss'])
        _r(f, 'Virial Ratio', (e['kin'] - w['e_tot']) / e['kin'])
        _rarr(f, 'External E-field', np.zeros(35))
        _i(f, 'IOpCl', 0)
        _i(f, 'IROHF', w['rohf'])
        if w['cas'] is not None:
            _i(f, 'Number of CAS Electrons', w['cas'][0])
            _i(f, 'Number of CAS Orbitals', w['cas'][1])
        _rarr(f, 'Alpha Orbital Energies', w['moe'][0])
        if len(mos) == 2:
            _rarr(f, 'Beta Orbital Energies', w['moe'][1])
        _rarr(f, 'Alpha MO coefficients', mos[0])
        if len(mos) == 2:
            _rarr(f, 'Beta MO coefficients', mos[1])
        _rarr(f, 'Total SCF Density', lower(p, pp))
        if w['kind'] != 'rhf':
            _rarr(f, 'Spin SCF Density', lower(ps, pp))
        # reference energies for ENPART (not stock .fchk fields)
        _r(f, 'Kinetic Energy', e['kin'])
        _r(f, 'Electron-Nuclei Energy', e['en'])
        _r(f, 'Electron-Electron Energy', e['vee'])
        _r(f, 'Coulomb Energy', e['coul'])
        _r(f, 'Exact-exchange Energy', e['exch'])
        if e['xmix'] is not None:
            _r(f, 'DFT-exchange Energy', e['exch_dft'])
            _r(f, 'Total Exchange Energy', e['exc'])
            if e['xmix'] != 0:
                _r(f, '% of Exact-exchange', e['xmix'])
        if e['vecp'] is not None:
            _rarr(f, 'ECP Matrix', lower(e['vecp'], 1.0 / pp))


#####################################################################
## RDM files                                                        ##
#####################################################################

def _nel(nelec):
    return str(tuple(int(n) for n in nelec))


def write_dm12(mol, mycas, name):
    """Writes name.dm1 (and name.dm2) for the # DM block (pySCF)."""
    toler = 1.0e-10
    if isinstance(mycas, mcscf.ucasci.UCASBase):
        raise Apost3dError('write_dm12: UHF-based CASSCF: ' +
                           _UNRESTRICTED_CORRELATED)
    if isinstance(mycas, mcscf.casci.CASBase):
        ncas, nelecas = int(mycas.ncas), mycas.nelecas
        dm1s = mycas.fcisolver.make_rdm1s(mycas.ci, ncas, nelecas)
        dm2 = mycas.fcisolver.make_rdm2(mycas.ci, ncas, nelecas)
        head = 'mcscf spin resolved rdm1 for  %s electrons in %d orbitals' \
            % (_nel(nelecas), ncas)
    elif isinstance(mycas, fci.direct_spin1.FCIBase):
        ncas = int(mycas.norb)
        nelecas = mycas.nelec
        if isinstance(ncas, tuple) or mycas.ci is None:
            raise Apost3dError('write_dm12: run the FCI calculation first')
        dm1s = mycas.make_rdm1s(mycas.ci, ncas, nelecas)
        dm2 = mycas.make_rdm2(mycas.ci, ncas, nelecas)
        head = 'fci spin resolved rdm1 for  %s electrons in %d orbitals' \
            % (_nel(nelecas), ncas)
    elif isinstance(mycas, cc.ccsd.CCSDBase):
        if isinstance(mycas, cc.uccsd.UCCSD):
            raise Apost3dError('write_dm12: UCCSD has different alpha and '
                               'beta orbitals; not available')
        ncas = int(mycas.nmo)
        d = 0.5 * mycas.make_rdm1()
        dm1s, dm2, nelecas = (d, d), None, None
        head = 'ccsd spin resolved rdm1 for  %d orbitals' % ncas
    else:
        raise Apost3dError('write_dm12: CASSCF/CASCI, FCI or CCSD objects '
                           'only (not %s)' % type(mycas).__name__)

    with open(name + '.dm1', 'w') as f:
        f.write(head + '\n')
        for case in (0, 1):
            for i in range(ncas):
                for j in range(i, ncas):
                    if abs(dm1s[case][i, j]) >= toler:
                        f.write('%d %d %22.15E\n' % (2 * i + 1 + case,
                                                     2 * j + 1 + case,
                                                     dm1s[case][i, j]))
    if dm2 is None:
        print('apost3d.py: no 2-RDM for CCSD (SPIN and ENPART need it); '
              'only %s.dm1 written' % name)
        return
    # spinless 2-RDM, every element in the [1212] convention
    with open(name + '.dm2', 'w') as f:
        f.write('spinless rdm2 for  %s electrons in %d orbitals\n'
                % (_nel(nelecas), ncas))
        for i in range(ncas):
            for j in range(ncas):
                for k in range(ncas):
                    for l in range(ncas):
                        if abs(dm2[i, k, j, l]) >= toler:
                            f.write('%d %d %d %d %22.15E\n'
                                    % (i + 1, j + 1, k + 1, l + 1,
                                       dm2[i, k, j, l]))
