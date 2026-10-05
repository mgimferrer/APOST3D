# Block section # FRAGMENTS

Defines the fragments used with `DOFRAGS`; see [Fragments](../guide/fragments.md)
for what they change in the results. Atoms are numbered as in the `.fchk`
file (the order of the geometry of the calculation).

The block holds the number of fragments, then, for each fragment, a line
with its number of atoms and a line with their numbers:

```text
# FRAGMENTS
3
1
1
2
3 5
3
2 4 6
#
```

Here fragment 1 is atom 1, fragment 2 atoms 3 and 5, and fragment 3 atoms
2, 4 and 6.

**All remaining atoms.** Writing `-1` as the number of atoms of the
**last** fragment puts every atom not listed before in it, with no atom
line:

```text
# FRAGMENTS
2
1
1
-1
#
```

(fragment 1 is atom 1, fragment 2 all the others).

- The run stops, with a message naming the problem, for an atom number
  that does not exist, an atom listed twice, a fragment with no atoms,
  more fragments than atoms, or `-1` on a fragment that is not the last.
- `EOS`, `GEOS`, `OSLO` and `ENPART` need every atom in a fragment and stop
  otherwise (`atoms in no fragment: ...`). The other analyses run with the
  listed fragments only and print a warning naming the atoms left out.
- `EOS`, `GEOS` and `OSLO` need fragments: without `DOFRAGS` they stop. For
  one result per atom, make every atom a fragment of its own.
- A long atom list may continue on the next line.
- The program echoes the fragments in the `INPUT SUMMARY` of the output
  (`Fragment  2 :   3   5`): check them there.
