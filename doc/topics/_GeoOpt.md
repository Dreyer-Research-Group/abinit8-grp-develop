---
description: How to perform a geometry optimization
authors: GG
---
<!--- This is the source file for this topics. Can be edited. -->

This page gives hints on how to perform a geometry optimization with the ABINIT package.

## Introduction

Different algorithms (Broyden; modified Broyden; Verlet with sudden stop of
atoms) allows to find the equilibrium configuration of the nuclei, for which
the forces vanish, see [[ionmov]]. The cell parameters can also be optimized
concurently with the atomic positions [[optcell]], possibly with a state of
given stress as a target, [[strtarget]].

Specified lattice parameters, or angles, or atomic positions, can be kept
fixed if needed, see [[topic:GeoConstraints]].

A genetic algorithm has been coded, for global optimisation. Not in production
yet. See [[ga_rules]].



## Related Input Variables

{{ related_variables }}

## Selected Input Files

{{ selected_input_files }}

