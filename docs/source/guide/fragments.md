# Fragments

Chemical questions are often about groups of atoms rather than single
atoms: the oxidation state of a metal and its ligands, the charge
transferred between two molecules, the bond between a carbene and a
metal. With `DOFRAGS` in `# METHOD`, atoms are grouped into
**fragments**, defined in a [`# FRAGMENTS` block](../input/fragments.md),
and the analyses work with them:

```text
# METHOD
TFVC
EOS
DOFRAGS
#
# FRAGMENTS
2
1
1
-1
#
```

(fragment 1 is atom 1, e.g. the metal; fragment 2 all the other atoms).

| Analysis | Effect of fragments |
|---|---|
| [Population analysis](../methods/population.md) | Every table is repeated with the populations, charges and valences summed over each fragment, and the bond orders between fragments. |
| [EFFAO](../methods/effao.md) | The effective orbitals belong to the fragment as a whole: for a ligand they resemble the orbitals of the free ligand. |
| [EOS](../methods/eos.md), [GEOS](../methods/geos.md) | Oxidation states of the fragments. |
| [OSLO](../methods/oslo.md) | Required: orbitals are localized onto fragments. |
| [SPIN](../methods/spin.md) | Local spins and couplings summed over each fragment. |
| [ENPART](../methods/enpart.md) | Every energy matrix is also condensed to fragments: self-energies and interaction energies of the fragments. |

Without `DOFRAGS`, every atom is its own fragment (EOS then gives atomic
oxidation states).

**Choosing fragments.** For a transition-metal complex, the natural choice
is one fragment per metal center and one per ligand, as a chemist would
count electrons. Ligands that are chemically one unit (a cyanide, a
carbonyl, a cyclopentadienyl ring) should not be split, and every atom
must belong to exactly one fragment.
