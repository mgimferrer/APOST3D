#!/usr/bin/env python3
"""
orca2fchk.py -- APOST-3D .fchk file from an ORCA calculation.

    orca2fchk.py name                    # name.gbw, name.property.txt -> name.fchk
    orca2fchk.py name -o other.fchk      # another output name
    orca2fchk.py name --out job.log      # ORCA output file other than name.out

Writes name.fchk (Gaussian formatted-checkpoint layout) with the basis,
orbitals, densities and the reference energies ENPART needs, from the files
an ORCA single point leaves: name.gbw and name.property.txt, converted with
ORCA's orca_2json (ORCA_PATH or orca_2json in PATH: load the ORCA module
first). Without orca_2json, name.json and name.property.json written by it
beforehand are read instead (name.json.conf: MOCoefficients, Basisset and
1elIntegrals S, T, V).

Wavefunctions: RHF, UHF, ROHF, RKS, UKS, ROKS, and single-root CASSCF.
For RKS/UKS/ROKS the method field of the .fchk carries the R/U/RO prefix
and the functional keyword of the ORCA input (from name.out), as Gaussian
writes it: RB3LYP, UPBE0, RB3LYP/G; without one known keyword (custom
%method functional, no output file) it is RKS/UKS/ROKS, with a note.
APOST-3D does not read the name: give the functional in # ENPART; for the
ORCA functionals checked against ORCA's exchange-correlation energy a
note gives their libxc ids (or says that libxc has no equivalent).
Pseudopotentials are supported (the ECP matrix is written for ENPART). The
energies written are those of the vacuum Hamiltonian: a CPCM dielectric
energy is left out, with a printed note; ORCA's SCF energy excludes an
empirical dispersion (D3/D4) already. The ORCA output (name.out; optional
for SCF) gives the basis name, a cross-check of the reference energies,
and stops the conversion of a scalar-relativistic (ZORA, DKH, X2C)
calculation or one with an external electric field or point charges.

CASSCF: name.dm1 and name.dm2 (active-space spin-resolved 1-RDM and
spinless 2-RDM, the # DM files of APOST-3D) are written too, built from the
CI vector that ORCA prints in determinants: the %casscf block needs
"PrintWF det" and "ci TPrintWF 1e-14 end", and name.out is then required.
The RDMs are checked against ORCA's CASSCF density.

ORCA's pure functions are reordered and their signs changed to the .fchk
convention (p: z,x,y -> x,y,z; f(+-3), g(+-3), g(+-4): sign flipped). The
overlap matrix of the basis written is recomputed and compared with ORCA's,
so a change of these conventions in a new ORCA version stops the
conversion. Tested with ORCA 6.1.1.

Errors are printed with the reason (exit code 1); no partial file is left
behind. Needs Python 3 (3.6 or newer) and numpy (1.x or 2.x).
"""

import argparse
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np


class Orca2fchkError(Exception):
    """An ORCA calculation that APOST-3D cannot take."""


BOHR = 0.529177210903       # Angstrom; ORCA 6's value (CODATA 2018)
TESTED_MAJOR = 6

_SHELL_L = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4, 'h': 5, 'i': 6}

# density-functional keywords of ORCA's simple input, upper case, each one
# run with ORCA 6.1.1 (accepted, DFT settings printed); the -3c composite
# methods are left out on purpose
_ORCA_FUNCTIONALS = frozenset('''
    HFS LDA LSD VWN VWN3 VWN5 PWLDA XHF HFLDA LB94
    BNULL BVWN BP BP86 BPW91 BLYP GP GLYP PW91 PBE BPBE REVPBE RPBE PWP
    OLYP OPBE XLYP MPWPW MPWLYP PW86PBE RPW86PBE B97-D B97-D3 B97-D4
    TPSS REVTPSS M06L SCANFUNC RSCAN R2SCAN B97M-V B97M-D3BJ B97M-D4
    B1LYP B3LYP B3LYP/G B3LYP_G O3LYP X3LYP G1LYP G3LYP B1P B3P B3P86 G1P
    G3P B1PW B1PW91 B3PW B3PW91 PW1PW MPW1PW MPW1LYP B1PBE B3PBE PWP1
    PBE0 PW91_0 REVPBE0 REVPBE38 BHANDHLYP
    TPSSH TPSS0 PW6B95 M06 M062X R2SCANH R2SCAN0 R2SCAN50
    WB97 WB97X WB97X-D3 WB97X-D3BJ WB97X-V WB97X-D4 WB97X-D4REV WB97M-V
    WB97M-D3BJ WB97M-D4 WB97M-D4REV CAM-B3LYP LC-BLYP LC-PBE WR2SCAN
    B2PLYP B2GP-PLYP B2K-PLYP B2T-PLYP B2NC-PLYP MPW2PLYP PWPB95 DSD-BLYP
    DSD-PBEP86 DSD-PBEB95 REVDSD-PBEP86-D4/2021 REVDOD-PBEP86-D4/2021
    PBE0-DH PBE0-2 PBE-QIDH R2SCAN0-DH R2SCAN-CIDH R2SCAN-QIDH R2SCAN0-2
    PR2SCAN50 PR2SCAN69 KPR2SCAN50 WPR2SCAN50 WB2PLYP WB2GP-PLYP WB97X-2
'''.split())

# libxc ids of # ENPART (LIBRARY) for ORCA functionals: the exchange-
# correlation energy of ENPART equals ORCA's E(XC) within 0.02 kcal/mol
# (water, def2-TZVP, ORCA 6.1.1; B3LYP, PBE0, TPSSh also UKS)
_ENPART_LIBXC = {
    'LSD': 'EX_FUNCTIONAL 1 + EC_FUNCTIONAL 7',
    'BLYP': 'EX_FUNCTIONAL 106 + EC_FUNCTIONAL 131',
    'PBE': 'EX_FUNCTIONAL 101 + EC_FUNCTIONAL 130',
    'REVPBE': 'EX_FUNCTIONAL 102 + EC_FUNCTIONAL 130',
    'TPSS': 'EX_FUNCTIONAL 202 + EC_FUNCTIONAL 231',
    'R2SCAN': 'EX_FUNCTIONAL 497 + EC_FUNCTIONAL 498',
    'M06L': 'EX_FUNCTIONAL 203 + EC_FUNCTIONAL 233',
    'B3LYP': 'EXC_FUNCTIONAL 475',
    'B3LYP/G': 'EXC_FUNCTIONAL 402',
    'B3LYP_G': 'EXC_FUNCTIONAL 402',
    'PBE0': 'EXC_FUNCTIONAL 406',
    'B3PW91': 'EXC_FUNCTIONAL 401',
    'BHANDHLYP': 'EXC_FUNCTIONAL 436',
    'TPSSH': 'EXC_FUNCTIONAL 457',
    'TPSS0': 'EXC_FUNCTIONAL 396',
    'M06': 'EX_FUNCTIONAL 449 + EC_FUNCTIONAL 235',
    'M062X': 'EX_FUNCTIONAL 450 + EC_FUNCTIONAL 236',
    'PW6B95': 'EXC_FUNCTIONAL 451',
    'R2SCANH': 'EXC_FUNCTIONAL 659',
    'R2SCAN0': 'EXC_FUNCTIONAL 660',
}
# ORCA functionals that no libxc id reproduces (same check)
_ENPART_NO_LIBXC = {
    'BP86': 'its P86 correlation is built on PW92',
    'B3P86': 'its P86 correlation is built on PW92',
    'X3LYP': 'its local correlation is VWN5, libxc\'s X3LYP has VWN-RPA',
}

# wavefunction methods of APOST-3D's method-field rule (sources/input2.f,
# subroutine input): every other method word reads as KS-DFT
_WFN_METHODS = {'CASSCF', 'HF', 'CID', 'CIS', 'CISD', 'QCISD', 'QCISD(T)',
                'QCISD(TQ)', 'CCD', 'CCSD', 'CCSD(T)', 'BD', 'BD(T)'}


def _note(text):
    print('orca2fchk.py: ' + text)


#####################################################################
## ORCA data (orca_2json)                                           ##
#####################################################################

_CONF = {
    'MOCoefficients': True,
    'Basisset': True,
    '1elIntegrals': ['S', 'T', 'V'],
    'JSONFormats': ['json'],
}


def _orca_2json():
    path = os.environ.get('ORCA_PATH')
    if path and os.path.isfile(os.path.join(path, 'orca_2json')):
        return os.path.join(path, 'orca_2json')
    return shutil.which('orca_2json')


def _run(cmd, cwd):
    r = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True)
    if r.returncode != 0:
        raise Orca2fchkError('%s failed:\n%s' % (' '.join(cmd),
                                                 r.stdout[-2000:]))


def _load_orca(name):
    """The orca_2json data of the wavefunction and of the property file.
    orca_2json runs in a temporary directory: nothing is added next to the
    ORCA files."""
    gbw, prop = name + '.gbw', name + '.property.txt'
    tool = _orca_2json()
    if tool is None:
        if os.path.isfile(name + '.json') and \
                os.path.isfile(name + '.property.json'):
            with open(name + '.json') as f:
                mol = json.load(f)
            with open(name + '.property.json') as f:
                props = json.load(f)
            return mol, props
        raise Orca2fchkError('orca_2json not found: load the ORCA module '
                             '(ORCA_PATH) or put orca_2json in PATH')
    for fil in (gbw, prop):
        if not os.path.isfile(fil):
            raise Orca2fchkError('%s not found' % fil)
    conf = dict(_CONF)
    try:            # the density file names its densities: ask for scfp
        with open(name + '.densitiesinfo', 'rb') as f:
            if b'scfp' in f.read():
                conf['Densities'] = ['scfp']
    except OSError:
        pass
    # the job name is kept: the density files refer to it
    base = os.path.basename(name)
    with tempfile.TemporaryDirectory(prefix='orca2fchk-') as tmp:
        for ext in ('.gbw', '.property.txt', '.densities', '.densitiesinfo'):
            if os.path.isfile(name + ext):
                os.symlink(os.path.abspath(name + ext),
                           os.path.join(tmp, base + ext))
        with open(os.path.join(tmp, base + '.json.conf'), 'w') as f:
            json.dump(conf, f)
        try:
            _run([tool, base + '.gbw', '-json'], tmp)
        except Orca2fchkError:
            if 'Densities' not in conf:
                raise
            # density files of another job name: go without them
            del conf['Densities']
            with open(os.path.join(tmp, base + '.json.conf'), 'w') as f:
                json.dump(conf, f)
            _run([tool, base + '.gbw', '-json'], tmp)
        _run([tool, base, '-property'], tmp)
        with open(os.path.join(tmp, base + '.json')) as f:
            mol = json.load(f)
        with open(os.path.join(tmp, base + '.property.json')) as f:
            props = json.load(f)
    return mol, props


def _orca_output(fil):
    """Basis name and energy components from the ORCA output; None
    without the file."""
    out = {}
    if fil is None or not os.path.isfile(fil):
        return None
    pats = {
        'basis': r'Your calculation utilizes the basis:\s*(\S+)',
        'enuc': r'^Nuclear Repulsion\s*:\s*(\S+)',
        'eone': r'^One Electron Energy:\s*(\S+)',
        'etwo': r'^Two Electron Energy:\s*(\S+)',
        'ekin': r'^Kinetic Energy\s*:\s*(\S+)',
        'xfun': r'Exchange Functional\s+Exchange\s+\.+\s*(\S+)',
        'cfun': r'Correlation Functional\s+Correlation\s+\.+\s*(\S+)',
        'ldac': r'LDA part of GGA corr\.\s+LDAOpt\s+\.+\s*(\S+)',
        'hfx': r'Fraction HF Exchange ScalHFX\s+\.+\s*(\S+)',
        # COSX: exchange recomputed on a larger grid after the SCF, in the
        # total energy but not in the components printed with it
        'dxfinal': r'Exchange energy change after final integration\s*:'
                   r'\s*(\S+)',
    }
    with open(fil, errors='replace') as f:
        text = f.read()
    for key, pat in pats.items():
        found = re.findall(pat, text, re.MULTILINE)
        if found:
            out[key] = found[-1]          # the last SCF of the file
    # Hamiltonians whose one-electron energies APOST-3D cannot reproduce
    echo = '\n'.join(re.findall(r'^\|\s*\d+>(.*)$', text, re.MULTILINE))
    if 'ORCA RELATIVISTIC HAMILTONIAN' in text:
        m = re.search(r'Relativistic Method\s+\.+\s*(\S+)', text)
        raise Orca2fchkError('scalar-relativistic Hamiltonian (%s): not '
                             'supported (APOST-3D is nonrelativistic)'
                             % (m.group(1) if m else '?'))
    if 'An Electric Field will be added' in text or \
            re.search(r'\befield\b', echo, re.IGNORECASE):
        raise Orca2fchkError('external electric field: not supported')
    if re.search(r'pointcharges|\bqmmm\b', echo, re.IGNORECASE):
        raise Orca2fchkError('external point charges (QM/MM): not '
                             'supported')
    out['functional'], out['functional_why'] = _functional(echo)
    return out


def _functional(echo):
    """The density-functional keyword of the ORCA input (its ! lines and
    %method Functional), or None and the reason."""
    text = '\n'.join(line.split('#')[0] for line in echo.splitlines())
    if re.search(r'\$new_job|%compound', text, re.IGNORECASE):
        return None, 'input with several jobs'
    words = []
    for line in text.splitlines():
        if line.strip().startswith('!'):
            words += line.strip()[1:].upper().split()
    for block in re.findall(r'%method\b(.*?)\bend\b', text,
                            re.IGNORECASE | re.DOTALL):
        if re.search(r'\b(exchange|correlation|ldaopt|scalhfx|scaldfx|'
                     r'scaldfc|scalldac)\b', block, re.IGNORECASE):
            return None, 'functional set in %method'
        for value in re.findall(r'\bfunctional\s+"?([^\s"]+)', block,
                                re.IGNORECASE):
            if value.upper() not in _ORCA_FUNCTIONALS:
                return None, 'functional %s set in %%method' % value
            words.append(value.upper())
    names = sorted({w for w in words if w in _ORCA_FUNCTIONALS})
    if len(names) != 1:
        return None, ('no known functional keyword in the input' if not names
                      else 'several functional keywords (%s)'
                      % ', '.join(names))
    return names[0], None


def _reads_as_wavefunction(word):
    """True when APOST-3D reads the method word as a wavefunction method:
    cut at '-', without the R/U/RO prefix, compared whole."""
    w = word.split('-')[0] if word.find('-') > 0 else word
    if w.startswith('RO'):
        w = w[2:]
    elif w[:1] in ('R', 'U'):
        w = w[1:]
    return w in _WFN_METHODS or (w[:2] == 'MP' and w[2:3].isdigit())


#####################################################################
## basis: ORCA -> .fchk order and signs                             ##
#####################################################################

def _shell_map(l):
    """Positions, within one ORCA shell, of the .fchk functions in order,
    and the sign of each. ORCA's pure order is the .fchk one (m = 0, +1,
    -1, +2, -2, ...) with the sign of |m| >= 3 flipped; its p order is
    z, x, y."""
    if l == 0:
        return [0], [1.0]
    if l == 1:
        return [1, 2, 0], [1.0, 1.0, 1.0]
    return list(range(2 * l + 1)), [1.0 if k < 5 else -1.0
                                     for k in range(2 * l + 1)]


def _renormalized(l, exps, coefs):
    """Contraction coefficients of normalized primitives, rescaled to a
    normalized contraction."""
    a = np.asarray(exps)
    c = np.asarray(coefs, dtype=float)
    s = (2.0 * np.sqrt(np.outer(a, a)) / np.add.outer(a, a)) ** (l + 1.5)
    return c / math.sqrt(c @ s @ c)


def _basis(mol, coords):
    """Shells in .fchk form, the ORCA AO index and sign of each .fchk
    function."""
    shells = []          # (l, exps, coeffs, atom, coord)
    perm, sign = [], []
    ao0 = 0
    for iat, atom in enumerate(mol['Atoms']):
        for sh in atom['Basis']:
            l = _SHELL_L[sh['Shell']]
            if l > 4:
                raise Orca2fchkError('basis functions with l > 4 (h and '
                                     'higher) are not supported by APOST-3D')
            exps = np.asarray(sh['Exponents'], dtype=float)
            shells.append((l, exps, _renormalized(l, exps,
                                                  sh['Coefficients']),
                           iat + 1, coords[iat]))
            local, sgn = _shell_map(l)
            perm += [ao0 + k for k in local]
            sign += sgn
            ao0 += 2 * l + 1
    return shells, np.array(perm), np.array(sign)


#####################################################################
## one-electron integrals (overlap check, ECP matrix)               ##
#####################################################################
# McMurchie-Davidson over Cartesian primitives, then real solid harmonics
# in the .fchk order and sign (the Gaussian/pySCF convention).

def _binom(n, k):
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))


def _cart(l):
    return [(lx, ly, l - lx - ly)
            for lx in range(l, -1, -1) for ly in range(l - lx, -1, -1)]


def _solid_harmonic(l, m):
    """Real solid harmonic S_lm as {(lx, ly, lz): coefficient}
    (Helgaker, Jorgensen, Olsen, eq. 6.4.48; overall scale irrelevant)."""
    am = abs(m)
    vm = 0.0 if m >= 0 else 0.5
    terms = {}
    for t in range((l - am) // 2 + 1):
        for u in range(t + 1):
            v = vm
            while v <= math.floor(am / 2.0 - vm) + vm + 1e-9:
                c = ((-1) ** int(round(t + v - vm)) * 0.25 ** t
                     * _binom(l, t) * _binom(l - t, am + t)
                     * _binom(t, u) * _binom(am, int(round(2 * v))))
                key = (int(round(2 * t + am - 2 * (u + v))),
                       int(round(2 * (u + v))), l - 2 * t - am)
                terms[key] = terms.get(key, 0.0) + c
                v += 1.0
    return terms


_SPH = {}


def _sph(l):
    """Cartesian -> pure transformation, rows in the .fchk order."""
    if l not in _SPH:
        cart = _cart(l)
        if l == 1:
            ms = [1, -1, 0]                   # x, y, z
        else:
            ms = [0] + [s * m for m in range(1, l + 1) for s in (1, -1)]
        t = np.zeros((len(ms), len(cart)))
        for i, m in enumerate(ms):
            for key, c in _solid_harmonic(l, m).items():
                t[i, cart.index(key)] = c
        _SPH[l] = t
    return _SPH[l]


def _dfact(n):
    return 1.0 if n <= 0 else n * _dfact(n - 2)


# Boys function: 6-term Taylor expansion around a tabulated grid,
# asymptotic form beyond it
_BOYS_NMAX = 8
_BOYS_DT = 0.05
_BOYS_TMAX = 40.0
_BOYS_GRID = np.arange(0.0, _BOYS_TMAX + 2 * _BOYS_DT, _BOYS_DT)


def _boys_series(n, t):
    term = np.full_like(t, 1.0 / (2 * n + 1))
    total = term.copy()
    for k in range(1, 300):
        term = term * 2.0 * t / (2 * n + 2 * k + 1)
        total += term
    return np.exp(-t) * total


_BOYS_TABLE = np.array([_boys_series(n, _BOYS_GRID)
                        for n in range(_BOYS_NMAX + 7)])


def _boys(nmax, t):
    """F_n(t), n = 0..nmax, as an array (nmax+1,) + t.shape."""
    f = np.empty((nmax + 1,) + t.shape)
    small = t < _BOYS_TMAX
    ts = t[small]
    ig = np.rint(ts / _BOYS_DT).astype(int)
    dt = _BOYS_GRID[ig] - ts
    tl = t[~small]
    for n in range(nmax + 1):
        acc = np.zeros_like(ts)
        fac = 1.0
        for k in range(6):
            acc += _BOYS_TABLE[n + k, ig] * fac
            fac *= dt / (k + 1)
        f[n][small] = acc
        f[n][~small] = (_dfact(2 * n - 1) / 2.0 ** (n + 1)
                        * np.sqrt(math.pi / tl ** (2 * n + 1)))
    return f


def _hermite(la, lb, p, xpa, xpb, e00):
    """Hermite expansion coefficients E[i, j, t] of one Cartesian
    direction, over the primitive pairs."""
    e = np.zeros((la + 1, lb + 1, la + lb + 1) + p.shape)
    e[0, 0, 0] = e00
    h = 0.5 / p
    for i in range(la + 1):
        for j in range(lb + 1):
            if i == 0 and j == 0:
                continue
            if i > 0:
                prev, x = e[i - 1, j], xpa
            else:
                prev, x = e[i, j - 1], xpb
            cur = e[i, j]
            cur += x * prev
            cur[1:] += h * prev[:-1]
            tt = np.arange(1, la + lb + 1).reshape((-1,) + (1,) * p.ndim)
            cur[:-1] += tt * prev[1:]
    return e


def _hermite_r(lmax, p, pc):
    """R[t, u, v] of the nuclear attraction, t+u+v <= lmax, over
    (primitive pairs, nuclei); pc = P - C, shape (..., 3)."""
    t2 = p[..., None] * np.sum(pc * pc, axis=-1)
    f = _boys(lmax, t2)
    m2p = -2.0 * p[..., None]
    prev = None
    for n in range(lmax, -1, -1):
        cur = {(0, 0, 0): m2p ** n * f[n]}
        for tot in range(1, lmax - n + 1):
            for t in range(tot, -1, -1):
                for u in range(tot - t, -1, -1):
                    v = tot - t - u
                    if t > 0:
                        val = pc[..., 0] * prev[(t - 1, u, v)]
                        if t > 1:
                            val = val + (t - 1) * prev[(t - 2, u, v)]
                    elif u > 0:
                        val = pc[..., 1] * prev[(t, u - 1, v)]
                        if u > 1:
                            val = val + (u - 1) * prev[(t, u - 2, v)]
                    else:
                        val = pc[..., 2] * prev[(t, u, v - 1)]
                        if v > 1:
                            val = val + (v - 1) * prev[(t, u, v - 2)]
                    cur[(t, u, v)] = val
        prev = cur
    r = np.zeros((lmax + 1,) * 3 + t2.shape)
    for key, val in prev.items():
        r[key] = val
    return r


def _prim_norm(l, a):
    return ((2.0 * a / math.pi) ** 0.75 * (4.0 * a) ** (0.5 * l)
            / math.sqrt(_dfact(2 * l - 1)))


def _pair(sha, shb, nuclei=None):
    """Pure overlap (and nuclear attraction of the point charges nuclei =
    (Z, coords)) blocks of two shells, before the normalization of the
    pure contracted functions."""
    la, ea, ca, _, ra = sha
    lb, eb, cb, _, rb = shb
    a, b = ea[:, None], eb[None, :]
    p = a + b
    ab = np.asarray(ra) - np.asarray(rb)
    w = (ca * _prim_norm(la, ea))[:, None] * (cb * _prim_norm(lb, eb))[None, :]
    pcent = (a[..., None] * np.asarray(ra) + b[..., None] * np.asarray(rb)) \
        / p[..., None]
    mu = a * b / p
    es = []
    for k in range(3):
        xpa = pcent[..., k] - ra[k]
        xpb = pcent[..., k] - rb[k]
        es.append(_hermite(la, lb, p, xpa, xpb, np.exp(-mu * ab[k] ** 2)))
    cart_a, cart_b = _cart(la), _cart(lb)
    idx = [(np.array([c[k] for c in cart_a])[:, None],
            np.array([c[k] for c in cart_b])[None, :]) for k in range(3)]
    ex, ey, ez = (es[k][idx[k][0], idx[k][1]] for k in range(3))
    # ex: (ncart_a, ncart_b, t, prim_a, prim_b)
    s = np.einsum('abpq,abpq,abpq,pq->ab', ex[:, :, 0], ey[:, :, 0],
                  ez[:, :, 0], w * (math.pi / p) ** 1.5)
    s = _sph(la) @ s @ _sph(lb).T
    if nuclei is None:
        return s, None
    z, rc = nuclei
    pc = pcent[:, :, None, :] - rc[None, None, :, :]
    r = _hermite_r(la + lb, p, pc)
    rz = np.einsum('tuvpqc,c->tuvpq', r, z)
    x1 = np.einsum('abvpq,tuvpq->abtupq', ez, rz)
    x2 = np.einsum('abupq,abtupq->abtpq', ey, x1)
    v = -np.einsum('abtpq,abtpq,pq->ab', ex, x2, w * 2.0 * math.pi / p)
    return s, _sph(la) @ v @ _sph(lb).T


def _integrals(shells, nuclei=None):
    """Overlap (and nuclear attraction) matrices of the .fchk basis."""
    nsh = len(shells)
    off = np.cumsum([0] + [2 * sh[0] + 1 for sh in shells])
    nao = off[-1]
    s = np.zeros((nao, nao))
    v = np.zeros((nao, nao)) if nuclei is not None else None
    for i in range(nsh):
        for j in range(i + 1):
            sij, vij = _pair(shells[i], shells[j], nuclei)
            s[off[i]:off[i + 1], off[j]:off[j + 1]] = sij
            s[off[j]:off[j + 1], off[i]:off[i + 1]] = sij.T
            if v is not None:
                v[off[i]:off[i + 1], off[j]:off[j + 1]] = vij
                v[off[j]:off[j + 1], off[i]:off[i + 1]] = vij.T
    norm = 1.0 / np.sqrt(np.diag(s))
    s *= np.outer(norm, norm)
    if v is not None:
        v *= np.outer(norm, norm)
    return s, v


#####################################################################
## wavefunction                                                     ##
#####################################################################

def _orbitals(mol, nao):
    """Orbital sets of the calculation, occupied orbitals first (CASSCF:
    as ORCA orders them, inactive, active, external)."""
    mos = mol['MolecularOrbitals']['MOs']
    hftyp = mol['HFTyp']
    nset = 2 if hftyp == 'UHF' else 1
    if len(mos) % nset:
        raise Orca2fchkError('odd number of orbitals for a UHF wavefunction')
    nmo = len(mos) // nset
    sets = []
    for k in range(nset):
        part = mos[k * nmo:(k + 1) * nmo]
        occ = np.array([m['Occupancy'] for m in part])
        c = np.array([m['MOCoefficients'] for m in part]).T
        e = np.array([m['OrbitalEnergy'] for m in part])
        if c.shape[0] != nao:
            raise Orca2fchkError('orbital coefficients do not match the '
                                 'basis (%d vs %d)' % (c.shape[0], nao))
        if hftyp == 'CASSCF':
            sets.append((c, e, occ))
            continue
        if np.any(abs(occ - np.rint(occ)) > 1e-6):
            raise Orca2fchkError('fractional orbital occupations (smearing, '
                                 'natural orbitals): not supported')
        order = np.argsort(-np.rint(occ), kind='stable')
        if np.any(order != np.arange(nmo)):
            _note('orbitals reordered: occupied orbitals first '
                  '(non-aufbau occupation)')
        sets.append((c[:, order], e[order], np.rint(occ[order])))
    return sets


def _energy(geom, method):
    """Energy entry of the method (SCF: no dispersion correction) of the
    last geometry."""
    for item in reversed(geom.get('Energy', [])):
        if item.get('Method') == method:
            return item
    raise Orca2fchkError('no %s energy in the property file' % method)


#####################################################################
## CASSCF: active-space RDMs from the printed CI vector             ##
#####################################################################

_CAS_INPUT = ('add "PrintWF det" and "ci TPrintWF 1e-14 end" to the '
              '%casscf block and rerun')


def _ci_vector(fil, nact):
    """Determinants (active occupation strings of 2, u, d, 0) and
    coefficients of the root, and ORCA's active density matrix, from the
    last CI printing of the ORCA output."""
    if fil is None or not os.path.isfile(fil):
        raise Orca2fchkError('CASSCF: the ORCA output %s (CI vector) is '
                             'needed' % fil)
    with open(fil, errors='replace') as f:
        text = f.read()
    k = text.rfind('Spin-Determinant CI Printing')
    if k < 0:
        raise Orca2fchkError('CASSCF: no CI vector in determinants in %s: '
                             % fil + _CAS_INPUT)
    dets = []
    lines = text[k:].split('\n')
    i = next(n for n, l in enumerate(lines) if l.startswith('ROOT')) + 1
    for line in lines[i:]:
        m = re.match(r'\s*\[([2ud0]+)\]\s+(\S+)\s*$', line)
        if m:
            dets.append((m.group(1), float(m.group(2))))
        elif dets or line.strip():
            break
    if not dets or any(len(s) != nact for s, _ in dets):
        raise Orca2fchkError('CASSCF: CI vector of %d active orbitals not '
                             'found in %s' % (nact, fil))
    k = text.rfind('\nDENSITY MATRIX\n')
    rows = text[k:].split('\n')[5:5 + nact] if k >= 0 else []
    try:
        dens = np.array([[float(x) for x in r.split()[1:]] for r in rows])
    except ValueError:
        dens = None
    if dens is None or dens.shape != (nact, nact):
        dens = None
    return dets, dens


def _parity(det, i):
    """Sign of moving an operator past the occupied orbitals below i."""
    return -1 if bin(det & ((1 << i) - 1)).count('1') % 2 else 1


def _strings(n, nel):
    """Occupation strings (bit p = orbital p) of nel electrons in n
    orbitals, and the excitations a+_r a_s of each: {(r, s): (target
    index, source index, sign) arrays}."""
    strs = [sum(1 << p for p in occ)
            for occ in itertools.combinations(range(n), nel)]
    index = {s: i for i, s in enumerate(strs)}
    exc = {}
    for j, s0 in enumerate(strs):
        for s in range(n):
            if not s0 >> s & 1:
                continue
            k = s0 ^ (1 << s)
            for r in range(n):
                if k >> r & 1:
                    continue
                sign = _parity(s0, s) * _parity(k, r)
                exc.setdefault((r, s), []).append((index[k | (1 << r)], j,
                                                   sign))
    exc = {key: tuple(np.array(v).T) for key, v in exc.items()}
    return index, exc


def _cas_rdms(dets, nact):
    """Spin-resolved 1-RDM (2, n, n) and spinless 2-RDM in pySCF's layout,
    dm2[p,q,r,s] = sum_st <p+_s r+_t s_t q_s>, of the active space.
    ORCA's determinants order the spin-orbitals orbital by orbital (1a 1b
    2a 2b ...); the CI vector is turned into a matrix C[alpha string, beta
    string] (alpha operators first), the excitations E_rs = E^a_rs + E^b_rs
    are applied to it, and dm2[p,q,r,s] = <E_qp Psi|E_rs Psi> -
    delta_qr gamma_ps."""
    occa = [sum(1 << p for p, ch in enumerate(s) if ch in '2u')
            for s, _ in dets]
    occb = [sum(1 << p for p, ch in enumerate(s) if ch in '2d')
            for s, _ in dets]
    na, nb = bin(occa[0]).count('1'), bin(occb[0]).count('1')
    ia, exa = _strings(nact, na)
    ib, exb = _strings(nact, nb)
    ci = np.zeros((len(ia), len(ib)))
    for (s, c), a, b in zip(dets, occa, occb):
        # interleaved -> alpha first: each beta electron moves past the
        # alpha electrons of the orbitals above it
        swaps = sum(bin(a >> (q + 1)).count('1')
                    for q in range(nact) if b >> q & 1)
        ci[ia[a], ib[b]] += -c if swaps % 2 else c
    n = nact
    t = np.zeros((n, n) + ci.shape)       # E_rs Psi
    g1s = np.zeros((2, n, n))
    for r in range(n):
        for s in range(n):
            ta = np.zeros_like(ci)
            if (r, s) in exa:
                tgt, src, sgn = exa[(r, s)]
                ta[tgt] = sgn[:, None] * ci[src]
            tb = np.zeros_like(ci)
            if (r, s) in exb:
                tgt, src, sgn = exb[(r, s)]
                tb[:, tgt] = ci[:, src] * sgn[None, :]
            g1s[0, r, s] = np.sum(ci * ta)
            g1s[1, r, s] = np.sum(ci * tb)
            t[r, s] = ta + tb
    tm = t.reshape(n * n, -1)
    dm2 = (tm @ tm.T).reshape((n,) * 4).transpose(1, 0, 2, 3).copy()
    gamma = g1s[0] + g1s[1]
    for q in range(n):
        dm2[:, q, q, :] -= gamma
    return g1s, dm2


def _write_rdms(base, g1s, dm2, nel):
    """name.dm1 and name.dm2 in the layout of utils/apost3d.py (# DM)."""
    toler = 1.0e-10
    nact = g1s.shape[1]
    tag = '(%d, %d) electrons in %d orbitals' % (nel[0], nel[1], nact)
    with open(base + '.dm1', 'w') as f:
        f.write('orca casscf spin resolved rdm1 for  %s\n' % tag)
        for case in (0, 1):
            for i in range(nact):
                for j in range(i, nact):
                    if abs(g1s[case][i, j]) >= toler:
                        f.write('%d %d %22.15E\n' % (2 * i + 1 + case,
                                                     2 * j + 1 + case,
                                                     g1s[case][i, j]))
    # spinless 2-RDM, every element in the [1212] convention
    with open(base + '.dm2', 'w') as f:
        f.write('spinless rdm2 for  %s\n' % tag)
        for i in range(nact):
            for j in range(nact):
                for k in range(nact):
                    for l in range(nact):
                        if abs(dm2[i, k, j, l]) >= toler:
                            f.write('%d %d %d %d %22.15E\n'
                                    % (i + 1, j + 1, k + 1, l + 1,
                                       dm2[i, k, j, l]))


def _coords(mol, geom):
    """Coordinates in bohr (property file), checked against the .gbw."""
    cart = geom['Geometry']['Coordinates']
    xyz = np.array([row[1:4] for row in cart['Cartesians']], dtype=float)
    if not cart.get('Units', 'a.u.').startswith('a.u'):
        xyz /= BOHR
    ang = np.array([a['Coords'] for a in mol['Atoms']], dtype=float)
    if not mol.get('CoordinateUnits', 'Angs').startswith('Ang'):
        ang *= BOHR
    if xyz.shape != ang.shape or abs(xyz - ang / BOHR).max() > 1e-5:
        raise Orca2fchkError('the geometry of the property file is not the '
                             'one of the .gbw file')
    return xyz


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


def convert(name, fchk=None, outfile=None):
    """Writes the .fchk of the ORCA calculation name (see the module
    docstring)."""
    if name.endswith('.gbw'):
        name = name[:-4]
    fchk = fchk or name + '.fchk'
    if outfile is None:
        outfile = name + '.out'
    # the output first: it stops the calculations APOST-3D cannot take
    out = _orca_output(outfile)
    if out is None:
        _note('ORCA output %s not found: Hamiltonian and energy checks '
              'skipped' % outfile)
        out = {}
    mol, props = _load_orca(name)

    version = mol.get('ORCA Header', {}).get('Version', '?')
    if not version.startswith('%d.' % TESTED_MAJOR):
        _note('ORCA %s: only ORCA %d.x is tested' % (version, TESTED_MAJOR))
    mol = mol['Molecule']
    geom = props['Geometries'][-1]

    coords = _coords(mol, geom)
    shells, perm, sign = _basis(mol, coords)
    nao = len(perm)
    hftyp = mol['HFTyp']
    if hftyp not in ('RHF', 'UHF', 'ROHF', 'CASSCF'):
        raise Orca2fchkError('wavefunction type %s not supported' % hftyp)
    cas = hftyp == 'CASSCF'
    ks = 'DFT_Energy' in geom and not cas
    main = 'CASSCF' if cas else 'SCF'
    others = sorted({e.get('Method') for e in geom.get('Energy', [])}
                    - {main})
    if others:
        _note('only the %s wavefunction is written (%s not converted)'
              % (main, ', '.join(others)))
    e_item = _energy(geom, main)
    e_scf = float(np.ravel(e_item['totalEnergy'])[-1])

    # the basis written must reproduce ORCA's overlap matrix
    zeff = np.array([a['NuclearCharge'] for a in mol['Atoms']], dtype=float)
    zatom = [int(a['ElementNumber']) for a in mol['Atoms']]
    ecp = any(abs(z - q) > 1e-8 for z, q in zip(zatom, zeff))
    s_orca = np.array(mol['S-Matrix'])
    s_fchk, v_nuc = _integrals(shells, (zeff, coords) if ecp else None)
    sgn = np.outer(sign, sign)
    err = abs(s_orca[np.ix_(perm, perm)] * sgn - s_fchk).max()
    if err > 1e-6:
        raise Orca2fchkError('the basis written does not reproduce ORCA\'s '
                             'overlap matrix (max. error %.1e): ORCA\'s '
                             'function order or sign convention changed?'
                             % err)

    sets = _orbitals(mol, nao)
    c_orca = [c for c, _, _ in sets]
    if cas:
        if e_item.get('NBlocks', 1) != 1 or e_item.get('NTotalRoots', 1) != 1:
            raise Orca2fchkError('CASSCF with several roots or blocks '
                                 '(state-averaged): not supported')
        nelcas = int(e_item['numOfActiveEl'])
        ncas = int(e_item['numOfActiveOrbs'])
        ncore = (int(e_item['numOfElectrons']) - nelcas) // 2
        dets, dens_orca = _ci_vector(outfile, ncas)
        norm = sum(c * c for _, c in dets)
        if abs(norm - 1.0) > 1e-6:
            raise Orca2fchkError('CASSCF: the CI vector printed has norm '
                                 '%.8f (truncated?): ' % norm + _CAS_INPUT)
        # the 9 printed decimals leave the norm off by ~1e-9
        dets = [(s, c / math.sqrt(norm)) for s, c in dets]
        g1s, dm2 = _cas_rdms(dets, ncas)
        nelact = (int(round(np.trace(g1s[0]))), int(round(np.trace(g1s[1]))))
        if dens_orca is not None and \
                abs(g1s[0] + g1s[1] - dens_orca).max() > 1e-5:
            raise Orca2fchkError('CASSCF: the 1-RDM of the CI vector is not '
                                 'ORCA\'s active density matrix')
        nalpha, nbeta = ncore + nelact[0], ncore + nelact[1]
        core = c_orca[0][:, :ncore]
        act = c_orca[0][:, ncore:ncore + ncas]
        da = core @ core.T + act @ g1s[0] @ act.T
        db = core @ core.T + act @ g1s[1] @ act.T
        if 'scfp' in mol.get('Densities', {}):
            dev = abs(da + db - np.array(mol['Densities']['scfp'])).max()
            if dev > 1e-6:
                raise Orca2fchkError('CASSCF: the density of the orbitals '
                                     'and CI vector differs from ORCA\'s by '
                                     '%.1e' % dev)
        else:
            _note('CASSCF: ORCA\'s density file not found, only the active '
                  'density matrix of the output checked')
        ca = cb = None
    elif hftyp == 'UHF':
        nalpha = int(sets[0][2].sum())
        nbeta = int(sets[1][2].sum())
        ca, cb = c_orca[0][:, :nalpha], c_orca[1][:, :nbeta]
    else:
        occ = sets[0][2]
        nbeta = int(np.sum(occ == 2))
        nalpha = nbeta + int(np.sum(occ == 1))
        ca, cb = c_orca[0][:, :nalpha], c_orca[0][:, :nbeta]
        if hftyp == 'RHF' and nalpha != nbeta:
            raise Orca2fchkError('RHF wavefunction with singly occupied '
                                 'orbitals')
    if not cas:
        da, db = ca @ ca.T, cb @ cb.T
    p = da + db

    # reference energies (ENPART), traces in ORCA's basis
    tr = lambda x, y: float(np.sum(x * y))
    t_orca, v_orca = np.array(mol['T-Matrix']), np.array(mol['V-Matrix'])
    e_kin, e_en = tr(p, t_orca), tr(p, v_orca)    # ORCA's V includes the ECP
    rab = np.sqrt(((coords[:, None] - coords[None]) ** 2).sum(-1))
    iu = np.triu_indices(len(zeff), 1)
    e_nuc = float(np.sum(np.outer(zeff, zeff)[iu] / rab[iu]))
    e_cpcm = float(geom.get('Solvation_Details', {}).get('CPCMDielEnergy',
                                                         0.0))
    # the energies written are those of the vacuum Hamiltonian (with the
    # density polarized by the solvent): ENPART checks their sum
    e_scf -= e_cpcm
    e_vee = e_scf - e_nuc - e_kin - e_en
    if e_cpcm:
        _note('CPCM dielectric energy %.10f au left out of the energies '
              'written (APOST-3D decomposes the energy in vacuo)' % e_cpcm)
    disp = geom.get('VdW_Correction', {}).get('vdW')
    if disp:
        _note('dispersion correction %.10f au: not in the energies written '
              '(ORCA\'s SCF energy excludes it)' % disp)
    if 'etwo' in out:
        out['etwo'] = float(out['etwo']) + float(out.get('dxfinal', 0.0))
    if cas:             # the components printed are the SCF guess's
        out = {k: v for k, v in out.items() if k not in ('ekin', 'etwo')}
    for key, value, label in (('enuc', e_nuc, 'nuclear repulsion'),
                              ('ekin', e_kin, 'kinetic energy'),
                              ('etwo', e_vee, 'two-electron energy')):
        if key in out and abs(float(out[key]) - value) > 1e-6:
            _note('WARNING: %s %.10f differs from the ORCA output\'s %s'
                  % (label, value, out[key]))

    # spin (CASSCF: ORCA's CI is spin-adapted)
    if cas:
        sz = 0.5 * (mol['Multiplicity'] - 1)
        ss = sz * (sz + 1)
    elif hftyp == 'UHF':
        sz = 0.5 * (nalpha - nbeta)
        ss = sz * (sz + 1) + nbeta - float(np.sum((ca.T @ s_orca @ cb) ** 2))
    else:
        sz = 0.5 * (nalpha - nbeta)
        ss = sz * (sz + 1)

    # ECP matrix: ORCA's V minus the nuclear attraction of the core charges
    v_ecp = None
    if ecp:
        v_ecp = v_orca[np.ix_(perm, perm)] * sgn - v_nuc

    # to the .fchk order and signs
    mos = [(c[perm, :] * sign[:, None]).T for c in c_orca]
    tril = np.tril_indices(nao)
    lower = lambda m: (m[np.ix_(perm, perm)] * sgn)[tril]

    if ks:
        # the functional's name as Gaussian writes it (RB3LYP); RKS, UKS,
        # ROKS when the input does not give one known keyword
        prefix = {'RHF': 'R', 'UHF': 'U', 'ROHF': 'RO'}[hftyp]
        header = prefix + 'KS'
        func = out.get('functional')
        if func is None:
            _note('functional not identified (%s): method field %s in the '
                  '.fchk' % (out.get('functional_why', 'no ORCA output'),
                             header))
        elif _reads_as_wavefunction(prefix + func):
            _note('functional %s would read as a wavefunction method in '
                  'APOST-3D: method field %s in the .fchk' % (func, header))
        else:
            header = prefix + func
        if func in _ENPART_LIBXC:
            _note('# ENPART for ORCA\'s %s: LIBRARY + %s (its exchange-'
                  'correlation energy checked against ORCA\'s)%s'
                  % (func, _ENPART_LIBXC[func], ', not the B3LYP keyword'
                     if func == 'B3LYP' else ''))
        elif func in _ENPART_NO_LIBXC:
            _note('ORCA\'s %s has no libxc equivalent (%s): for ENPART use a '
                  'Gaussian 16 or pySCF wavefunction'
                  % (func, _ENPART_NO_LIBXC[func]))
        elif out.get('xfun') == 'B88' and out.get('cfun') == 'LYP' and \
                out.get('ldac', '').upper().startswith('VWN-5'):
            _note('ORCA\'s B3LYP/BLYP-type functional uses VWN5: in # ENPART '
                  'use LIBRARY + EXC_FUNCTIONAL (B3LYP: 475), not B3LYP')
    else:
        header = hftyp
    basis_name = out.get('basis', 'Gen')
    lmax = max(sh[0] for sh in shells)
    nprim = [len(sh[1]) for sh in shells]

    tmp = fchk + '.tmp'
    try:
        with open(tmp, 'w') as f:
            f.write('Automatically generated by orca2fchk.py for job %s\n'
                    % os.path.basename(name))
            # method from column 11 (KS: the functional's name), basis after
            f.write('SP        %-15s %s\n' % (header, basis_name))
            _i(f, 'Number of atoms', len(zatom))
            _i(f, 'Charge', mol['Charge'])
            _i(f, 'Multiplicity', mol['Multiplicity'])
            _i(f, 'Number of electrons', nalpha + nbeta)
            _i(f, 'Number of alpha electrons', nalpha)
            _i(f, 'Number of beta electrons', nbeta)
            _i(f, 'Number of basis functions', nao)
            _i(f, 'Number of independent functions', mos[0].shape[0])
            _iarr(f, 'Atomic numbers', zatom)
            _rarr(f, 'Nuclear charges', zeff)
            _rarr(f, 'Current cartesian coordinates', coords)
            _i(f, 'Number of contracted shells', len(shells))
            _i(f, 'Number of primitive shells', sum(nprim))
            _i(f, 'Pure/Cartesian d shells', 0)
            _i(f, 'Pure/Cartesian f shells', 0)
            _i(f, 'Highest angular momentum', lmax)
            _i(f, 'Largest degree of contraction', max(nprim))
            _iarr(f, 'Shell types', [sh[0] if sh[0] < 2 else -sh[0]
                                     for sh in shells])
            _iarr(f, 'Number of primitives per shell', nprim)
            _iarr(f, 'Shell to atom map', [sh[3] for sh in shells])
            _rarr(f, 'Primitive exponents',
                  np.concatenate([sh[1] for sh in shells]))
            _rarr(f, 'Contraction coefficients',
                  np.concatenate([sh[2] for sh in shells]))
            _rarr(f, 'Coordinates of each shell', [sh[4] for sh in shells])
            _r(f, 'SCF Energy', e_scf)
            _r(f, 'Total Energy', e_scf)
            _r(f, 'S**2', ss)
            _r(f, 'Virial Ratio', (e_kin - e_scf) / e_kin)
            _rarr(f, 'External E-field', np.zeros(35))
            _i(f, 'IOpCl', 0)
            _i(f, 'IROHF', 1 if hftyp == 'ROHF' else 0)
            if cas:
                _i(f, 'Number of CAS Electrons', nelcas)
                _i(f, 'Number of CAS Orbitals', ncas)
            _rarr(f, 'Alpha Orbital Energies', sets[0][1])
            if len(mos) == 2:
                _rarr(f, 'Beta Orbital Energies', sets[1][1])
            _rarr(f, 'Alpha MO coefficients', mos[0])
            if len(mos) == 2:
                _rarr(f, 'Beta MO coefficients', mos[1])
            _rarr(f, 'Total SCF Density', lower(p))
            if hftyp != 'RHF':
                _rarr(f, 'Spin SCF Density', lower(da - db))
            # reference energies for ENPART (not stock .fchk fields)
            _r(f, 'Kinetic Energy', e_kin)
            _r(f, 'Electron-Nuclei Energy', e_en)
            _r(f, 'Electron-Electron Energy', e_vee)
            if v_ecp is not None:
                _rarr(f, 'ECP Matrix', v_ecp[tril])
        os.replace(tmp, fchk)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    if v_ecp is not None:
        _note('pseudopotential energy %.10f au (ECP matrix written)'
              % tr(p[np.ix_(perm, perm)] * sgn, v_ecp))
    _note('%s written (%s, %d basis functions, overlap check %.1e)'
          % (fchk, header, nao, err))
    if cas:
        base = fchk[:-5] if fchk.endswith('.fchk') else fchk
        _write_rdms(base, g1s, dm2, nelact)
        _note('%s.dm1 and %s.dm2 written: CAS(%d,%d), %d determinants; '
              'give them in the # DM block, with ORCA'
              % (base, base, nelcas, ncas, len(dets)))


def main():
    ap = argparse.ArgumentParser(
        description='APOST-3D .fchk file from an ORCA calculation '
                    '(name.gbw + name.property.txt).')
    ap.add_argument('name', help='ORCA job name (name.gbw)')
    ap.add_argument('-o', dest='fchk', help='.fchk file (default name.fchk)')
    ap.add_argument('--out', dest='outfile',
                    help='ORCA output file (default name.out; optional)')
    args = ap.parse_args()
    try:
        convert(args.name, args.fchk, args.outfile)
    except Orca2fchkError as err:
        sys.stderr.write('orca2fchk.py: %s\n' % err)
        sys.exit(1)


if __name__ == '__main__':
    main()
