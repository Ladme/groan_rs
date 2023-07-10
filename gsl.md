# Groan Selection Language

Groan selection language (GSL) is a query language for selecting atoms in the `groan_rs` library and associated programs. Inside the `groan_rs` library, it is mainly used in the `System::group_create()` function.

In general, GSL is similar to the selection language used by VMD.

*This document outlines the GSL v2.0 syntax. Most queries from GSL v1.0 are compatible with GSL v2.0, but the reverse is not true.*

## Basic queries
You can select atoms based on:
- Their **residue names** using `resname XYZ`. For instance, `resname POPE` will select all atoms of the system corresponding to residues named POPE.
- Their **residue numbers** using `resid XYZ` or `resnum XYZ`. For instance, `resid 17` will select all atoms of the system corresponding to residues with number 17.
- Their **atom names** using `name XYZ` or `atomname XYZ`. For instance, `name P` will select all atoms of the system which name is P.
- Their **real atom numbers** using `serial XYZ`. For instance, `serial 256` will select the atom with the number 256 (this is guaranteed to be a single atom). Note that in this case, the atoms are selected based on their "real" atom numbers, as understood by gromacs, **not** by the atom numbers specified in the `gro file` the system originates from.
- Their **gro atom numbers** using `atomid XYZ` or `atomnum XYZ`. For instance, `atomid 124` will select all atoms which have the atom number 124 in the `gro file` the system originates from. This can be multiple atoms.

## Multiple identifiers
You can specify multiple identifiers in the query. For instance, by using `resname POPE POPG`, you will select all atoms of the system corresponding to residues named POPE as well as atoms corresponding to residues named POPG. See examples of similar queries below:

- `resid 13 15 16 17` will select all atoms corresponding to residues with numbers 13, 15, 16, or 17.
- `name P CA HA` will select all atoms with atom names P, CA, or HA.
- `serial 245 267 269 271` will select atoms numbered 245, 267, 269, or 271.

## Selecting atoms using groups
You can also select atoms using previously created groups of atoms. For instance, if you have previously created groups named "Protein" and "Membrane", you can use query `group Protein Membrane` or just `Protein Membrane` to select atoms of these two groups. In case any of the groups does not exist, an error will be raised.

If you load an `ndx file` for your system using `System::read_ndx()`, you can also use groups from the `ndx file`.

In case your group consists of multiple words, you have to enclose it into quotes (' or "). For instance `Protein 'My Group'` will select all atoms of the group "Protein" as well as all atoms of the group "My Group".

## Selecting atoms by autodetection
You can select atoms using internally defined macros. Currently, `groan_rs` library provides four of such macros:
- `@membrane` will select all atoms corresponding to membrane lipids (supports over 200 membrane lipid types).
- `@protein` will select all atoms corresponding to either of the 20 standard amino acids.
- `@water` will select all atoms of water. Using this macro is not recommended as it only supports some water models.
- `@ions` will select all atoms of ions. Using this macro is not recommended as it only supports small number of ion types.

## Selecting all atoms
You can select all atoms of the system by using `all`.

## Ranges
Instead of writing residue or atom numbers explicitly, you can use keyword `to` or `-` to specify a range. For example, instead of writing `resid 14 15 16 17 18 19 20`, you can use `resid 14 to 20` or `resid 14-20`. This will select all atoms corresponding to residues with residue numbers 14, 15, 16, 17, 18, 19, or 20.

You can also specify multiple ranges in a single query or combine ranges with explicitly provided numbers. For example, `serial 1 3 to 6 10 12 - 14 17` will expand to `serial 1 3 4 5 6 10 12 13 14 17`.

## Negations
Using keyword `not` or `!` in front of a query will negate it. For example, the query `not name CA` or `! name CA` will select all atoms which name does **not** correspond to CA. Similarly, `not resname POPE POPG` will select all atoms that correspond to residues with names other than POPE or POPG. `!Protein` will then select all atoms that are not part of the group named `Protein`.

## Binary operations
You can combine basic queries by using `and` (`&&`) and `or` (`||`) operators.

Joining two queries by `and` will select only atoms that were selected by **both** of the queries. For example, `resname POPE and name P` will select all atoms that belong to residues named POPE and that have the name P. Similarly, `resid 17 18 && serial 256 to 271` will select only atoms corresponding to residue 17 or 18 and with atom numbers between 256 and 271 (including 271).

Joining two queries by `or` will select atoms that were selected by **either** of the queries (at least one of them). For example, `resname POPE or name P` will select all atoms that belong to residues named POPE as well as all atoms with the name P. Similarly, `resid 17 18 || serial 256 to 271` will select all atoms corresponding to residue 17 or 18 as well as all atoms with atom numbers between 256 and 271.

In case multiple `and` and/or `or` operators are used in a single query, they are evaluated from left to right. For example, `resname POPE or name CA and not Protein` will select all atoms belonging to residues named POPE or having the atom name CA but all these atoms can not belong to the group named `Protein`.

Autodetection macros can also be combined with other sub-queries using operators, i.e. `@membrane or group 'My Lipids'` will select all autodetected membrane lipids and all atoms of the group "My Lipids". 

## Parentheses
You can change the order in which the individual sub-queries and operations are evaluated by using parentheses `(` and `)`. Expressions enclosed in parentheses are evaluated first (think math). For example, `resname POPE or (name CA and not resid 18 to 21)` will select all atoms belonging to residues named POPE along with all atoms that
- have the atom name P **and**
- do not correspond to residues numbered 18 to 21. 

Meanwhile `(resname POPE or name CA) and not resid 18 to 21` is equivalent to `resname POPE or name CA and not resid 18 to 21`. This will select all atoms belonging to residues named POPE or having the atom name CA but all of these atoms can not belong to residues 18 to 21.

You can place parenthetical expressions into other parenthetical expressions. For example `serial 1 to 6 or (name CA and resname POPE || (resid 1 to 7 or serial 123 to 128)) and Protein` is a valid query, albeit possibly way too convoluted.

You can also place `not` (`!`) operator in front of a parenthetical expression. For example, `!(serial 1 to 6 && name P)` will select all atoms that do **not** have atom number between 1 and 6 while also having the atom name P.

## Note about whitespace
Operators and parentheses do not have to be separated by whitespace from the rest of the query, unless the meaning of the query would become unclear. For instance, `not(name CA)or(serial 1to45||Protein)` is a valid query, while `not(name CA)or(serial 1to45orProtein)` is **not** valid, as `orProtein` becomes uninterpretable. However, enclosing the `Protein` in parentheses, i.e. `not(name CA)or(serial 1to45or(Protein))`, turns the query into a valid one again.