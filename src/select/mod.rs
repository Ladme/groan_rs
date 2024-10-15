// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of the Groan Selection Language for selecting groups of atoms.

use fancy_regex::Regex;
use std::collections::HashMap;
use std::fmt::{self, Write};

use crate::errors::SelectError;
use crate::system::System;

use self::name::Name;

mod name;
mod numbers;

#[derive(Debug, PartialEq, Clone)]
pub enum Select {
    ResidueName(Vec<Name>),
    AtomName(Vec<Name>),
    ResidueNumber(Vec<(usize, usize)>),
    GmxAtomNumber(Vec<(usize, usize)>),
    AtomNumber(Vec<(usize, usize)>),
    Chain(Vec<char>),
    GroupName(Vec<Name>),
    LabeledAtom(Vec<Name>),
    ElementName(Vec<Name>),
    ElementSymbol(Vec<Name>),
    And(Box<Select>, Box<Select>),
    Or(Box<Select>, Box<Select>),
    Not(Box<Select>),
    Molecule(Box<Select>),
}

#[derive(Debug, PartialEq)]
enum Operator {
    And,
    Or,
    Not,
    Molecule,
}

impl Select {
    /// Construct a Selection tree (`Select` structure) from the given Groan Selection language query.
    pub fn parse_query(query: &str) -> Result<Box<Select>, SelectError> {
        // check that the expression is not empty
        if query.trim().is_empty() {
            return Err(SelectError::EmptyQuery);
        }

        // check that number of '(' is balanced with the number of ')'
        if !par_balanced(query) {
            return Err(SelectError::InvalidParentheses(query.to_string()));
        }

        // check that the number of ' and " quotes is even, i.e. all quote-blocks are closed
        if !quotes_balanced(query) {
            return Err(SelectError::InvalidQuotes(query.to_string()));
        }

        // expand macros
        let mut expression = query.to_string();
        if expression.contains('@') {
            let macros = get_macros();
            expand_macros(&mut expression, macros);
        }

        // replace `mol with` and `molecule with` with `@@`
        let molwith_pattern = Regex::new(
            r"(molecule\s*with|mol\s*with)(?=(?:[^']*'[^']*')*[^']*$)",
        )
        .expect("FATAL GROAN ERROR | select::parse_query | Could not construct regex pattern.");
        expression = molwith_pattern.replace_all(&expression, "@@").to_string();

        // replace word operators with their symbolic equivalents
        expression = replace_keywords(&expression);

        match parse_subquery(&expression, 0, expression.chars().count()) {
            Ok(x) => Ok(x),
            Err(SelectError::InvalidOperator(_)) => {
                Err(SelectError::InvalidOperator(query.to_string()))
            }
            Err(SelectError::MissingArgument(_)) => {
                Err(SelectError::MissingArgument(query.to_string()))
            }
            Err(SelectError::EmptyArgument(_)) => {
                Err(SelectError::EmptyArgument(query.to_string()))
            }
            Err(SelectError::InvalidParentheses(_)) => {
                Err(SelectError::InvalidParentheses(query.to_string()))
            }
            Err(SelectError::InvalidNumber(_)) => {
                Err(SelectError::InvalidNumber(query.to_string()))
            }
            Err(SelectError::InvalidChainId(_)) => {
                Err(SelectError::InvalidChainId(query.to_string()))
            }
            Err(SelectError::InvalidRegex(e)) => Err(SelectError::InvalidRegex(e)),
            Err(SelectError::InvalidTokenParentheses(_)) => {
                Err(SelectError::InvalidTokenParentheses(query.to_string()))
            }
            Err(_) => Err(SelectError::UnknownError(query.to_string())),
        }
    }

    /// Expand regular expressions for `GroupName` and `LabeledAtom` into actuall group names
    /// and labeled atoms matching the regex.
    ///
    /// Checks that all `String` groups actually exist and ensures that at least one group is present in each `Select::GroupName`.
    /// Checks that all `String` labeled atoms actually exist and ensures that at least one label is present in each `Select::LabeledAtom`.
    ///
    /// Performing this expansion once before applying the `Select` operation
    /// is more efficient than performing similar expansion for each individual atom.
    pub(crate) fn expand_regex_group_label(self, system: &System) -> Result<Self, SelectError> {
        match self {
            Select::GroupName(vector) => {
                let new_vector = Self::expand_vector(
                    &vector,
                    system,
                    |s| system.group_exists(s),
                    |system| system.get_groups_as_ref().keys(),
                    SelectError::GroupNotFound,
                )?;
                Ok(Select::GroupName(new_vector))
            }
            Select::LabeledAtom(vector) => {
                let new_vector = Self::expand_vector(
                    &vector,
                    system,
                    |s| system.label_exists(s),
                    |system| system.get_labeled_atoms().keys(),
                    SelectError::LabelNotFound,
                )?;
                Ok(Select::LabeledAtom(new_vector))
            }
            Select::And(left, right) => Ok(Select::And(
                Box::from(left.expand_regex_group_label(system)?),
                Box::from(right.expand_regex_group_label(system)?),
            )),
            Select::Or(left, right) => Ok(Select::Or(
                Box::from(left.expand_regex_group_label(system)?),
                Box::from(right.expand_regex_group_label(system)?),
            )),
            Select::Not(op) => Ok(Select::Not(Box::from(op.expand_regex_group_label(system)?))),
            Select::Molecule(op) => Ok(Select::Molecule(Box::from(
                op.expand_regex_group_label(system)?,
            ))),
            other => Ok(other),
        }
    }

    /// Helper function for expanding regular expressions for groups and labeled atoms.
    fn expand_vector<'a, F, I, E>(
        vector: &'a [Name],
        system: &'a System,
        exists_fn: F,
        keys_fn: fn(&'a System) -> I,
        error_fn: fn(String) -> E,
    ) -> Result<Vec<Name>, SelectError>
    where
        F: Fn(&str) -> bool,
        I: Iterator<Item = &'a String>,
        E: Into<SelectError>,
    {
        let mut new_vector = Vec::new();

        for name in vector {
            match name {
                Name::String(s) => {
                    if !exists_fn(s) {
                        return Err(error_fn(s.clone()).into());
                    }
                    new_vector.push(Name::String(s.to_string()));
                }
                Name::Regex(r) => {
                    for key in keys_fn(system) {
                        if r.is_match(key) {
                            new_vector.push(Name::String(key.to_owned()))
                        }
                    }
                }
            }
        }

        if new_vector.is_empty() {
            return Err(SelectError::NoRegexMatch(vector[0].to_string()));
        }

        Ok(new_vector)
    }

    /// Validate that all groups specified in the `Select` structure actually exist in the system.
    /// Returns `Ok` if they all exist, otherwise returns a `SelectError`
    pub(crate) fn validate_groups(&self, system: &System) -> Result<(), SelectError> {
        match self {
            Select::GroupName(names) => {
                for name in names.iter() {
                    // we can use the index 0 even if the system contains no atoms
                    // due to the way `AtomContainer::isin` works
                    match name.match_groups(system, 0) {
                        Ok(_) => (),
                        Err(e) => return Err(e),
                    }
                }
            }

            Select::And(left, right) | Select::Or(left, right) => {
                left.validate_groups(system)?;
                right.validate_groups(system)?;
            }

            Select::Not(operand) | Select::Molecule(operand) => operand.validate_groups(system)?,

            _ => (),
        }

        Ok(())
    }

    fn write_range(string: &mut String, start: usize, end: usize) {
        if start == end {
            write!(string, "{} ", start).unwrap();
        } else if start == 1 {
            write!(string, "<= {} ", end).unwrap();
        } else if end == usize::MAX {
            write!(string, ">= {} ", start).unwrap();
        } else {
            write!(string, "{} to {} ", start, end).unwrap();
        }
    }

    fn write_name(string: &mut String, name: &Name) {
        match name {
            Name::String(x) => {
                if x.contains(char::is_whitespace) {
                    write!(string, "'{}' ", x)
                } else {
                    write!(string, "{} ", x)
                }
            }
            Name::Regex(x) => write!(string, "r'{}' ", x),
        }
        .unwrap()
    }

    /// Check whether the target Selection contains an empty vector.
    fn is_empty(&self) -> bool {
        match self {
            Select::ResidueName(v)
            | Select::AtomName(v)
            | Select::GroupName(v)
            | Select::LabeledAtom(v)
            | Select::ElementName(v)
            | Select::ElementSymbol(v) => v.is_empty(),

            Select::ResidueNumber(v) | Select::GmxAtomNumber(v) | Select::AtomNumber(v) => {
                v.is_empty()
            }

            Select::Chain(v) => v.is_empty(),

            Select::And(..) | Select::Or(..) | Select::Not(_) | Select::Molecule(_) => false,
        }
    }
}

impl fmt::Display for Select {
    /// Convert `Select` structure into a valid Groan Selection Language query.
    ///
    /// ## Warning
    /// - The GSL query returned by this function may be different from the query
    ///   used to construct the Select structure in the first place.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut string = String::new();

        // check for empty vectors
        // this is necessary since the select can actually be empty
        // (for instance, the query `resid < 0` will lead to empty vector)
        // (similarly, expanding REGEX query to no matching groups will lead to an empty vector)
        // we can't just ignore Selections with empty vectors
        // we use 'none' written as '!(all)' to represent empty Selections
        if self.is_empty() {
            write!(&mut string, "!(all)").unwrap();
            return write!(f, "{}", string);
        }

        match self {
            Self::ResidueName(vector) => {
                write!(&mut string, "resname ").unwrap();
                for name in vector {
                    Select::write_name(&mut string, name);
                }
            }

            Self::AtomName(vector) => {
                write!(&mut string, "name ").unwrap();
                for name in vector {
                    Select::write_name(&mut string, name);
                }
            }

            Self::ResidueNumber(vector) => {
                write!(&mut string, "resid ").unwrap();
                for (start, end) in vector {
                    Select::write_range(&mut string, *start, *end);
                }
            }

            Self::GmxAtomNumber(vector) => {
                write!(&mut string, "serial ").unwrap();
                for (start, end) in vector {
                    Select::write_range(&mut string, *start, *end);
                }
            }

            Self::AtomNumber(vector) => {
                write!(&mut string, "atomid ").unwrap();
                for (start, end) in vector {
                    Select::write_range(&mut string, *start, *end);
                }
            }

            Self::Chain(vector) => {
                write!(&mut string, "chain ").unwrap();
                for chain in vector {
                    write!(&mut string, "{} ", chain).unwrap();
                }
            }

            Self::GroupName(vector) => {
                write!(&mut string, "group ").unwrap();
                for name in vector {
                    Select::write_name(&mut string, name);
                }
            }

            Self::LabeledAtom(vector) => {
                write!(&mut string, "label ").unwrap();
                for name in vector {
                    Select::write_name(&mut string, name);
                }
            }

            Self::ElementName(vector) => {
                write!(&mut string, "element name ").unwrap();
                for name in vector {
                    Select::write_name(&mut string, name);
                }
            }

            Self::ElementSymbol(vector) => {
                write!(&mut string, "element symbol ").unwrap();
                for name in vector {
                    Select::write_name(&mut string, name);
                }
            }

            Self::And(left, right) => {
                write!(
                    &mut string,
                    "({} and {})",
                    Select::to_string(left),
                    Select::to_string(right),
                )
                .unwrap();
            }

            Select::Or(left, right) => {
                write!(
                    &mut string,
                    "({} or {})",
                    Select::to_string(left),
                    Select::to_string(right),
                )
                .unwrap();
            }

            Select::Not(op) => {
                write!(&mut string, "not {}", Select::to_string(op),).unwrap();
            }

            Select::Molecule(op) => {
                write!(&mut string, "molecule with {}", Select::to_string(op),).unwrap();
            }
        }

        write!(f, "{}", string.trim())
    }
}

fn parse_subquery(expression: &str, start: usize, end: usize) -> Result<Box<Select>, SelectError> {
    let mut tree: Option<Box<Select>> = None;

    let mut i = start;

    let mut token = String::new();
    let mut unary_operators: Vec<Operator> = Vec::new();
    let mut binary_operator: Option<Operator> = None;

    let mut inside_regex = false;

    while i < end {
        let c = expression.chars().nth(i).unwrap();

        // when inside of a regex block, ignore all the operators
        if inside_regex {
            if c == '\'' {
                inside_regex = false;
            }

            token.push(c);
            i += 1;
            continue;
        }

        match c {
            '(' => {
                let new_end = match find_parenthesis(expression, i, end) {
                    Some(x) => x,
                    None => {
                        panic!("FATAL GROAN ERROR | select::parse_subquery | Despite all efforts, parentheses are not balanced.");
                    }
                };

                let parsed = parse_subquery(expression, i + 1, new_end)?;
                tree = process_operation(tree, parsed, &mut unary_operators, &binary_operator)?;

                binary_operator = None;

                i = new_end + 1;
            }

            ')' => i += 1,

            // binary operators
            '&' | '|' => {
                let operator = find_operator(expression, c, i);
                match operator {
                    Some(Operator::Or) | Some(Operator::And) => (),
                    _ => return Err(SelectError::InvalidOperator("".to_string())),
                }

                if !token.trim().is_empty() {
                    // parse the token and process the queued operations
                    let parsed = Box::from(parse_token(&token)?);
                    tree = process_operation(tree, parsed, &mut unary_operators, &binary_operator)?;
                    token.clear();
                }

                // set the new binary operator
                binary_operator = operator;
                i += 2;
            }

            // NOT operator
            '!' => {
                unary_operators.push(Operator::Not);
                i += 1;
            }

            // `molecule` operator (represented as `@@`)
            '@' => {
                let operator = find_operator(expression, c, i);
                match operator {
                    Some(Operator::Molecule) => (),
                    _ => return Err(SelectError::InvalidOperator("".to_string())),
                }

                unary_operators.push(Operator::Molecule);
                i += 2;
            }

            // regex
            'r' => {
                if expression.get(i + 1..i + 2) == Some("'") {
                    token.push('r');
                    token.push('\'');
                    i += 2;
                    inside_regex = true;
                } else {
                    token.push(c);
                    i += 1;
                }
            }

            _ => {
                token.push(c);
                i += 1;
            }
        }
    }

    // process the last operation
    if !token.trim().is_empty() {
        let parsed = Box::from(parse_token(&token)?);
        tree = process_operation(tree, parsed, &mut unary_operators, &binary_operator)?;
    } else if binary_operator.is_some() {
        return Err(SelectError::MissingArgument("".to_string()));
    }

    match tree {
        Some(x) => Ok(x),
        None => Err(SelectError::UnknownError("".to_string())),
    }
}

fn process_operation(
    tree: Option<Box<Select>>,
    mut parsed: Box<Select>,
    unary: &mut Vec<Operator>,
    binary: &Option<Operator>,
) -> Result<Option<Box<Select>>, SelectError> {
    // modify the parsed token using unary operators
    for operator in unary.iter() {
        parsed = match operator {
            Operator::Not => Box::from(Select::Not(parsed)),
            Operator::Molecule => Box::from(Select::Molecule(parsed)),
            Operator::And => panic!(
                "FATAL GROAN ERROR | select::process_operation | AND operator is being treated as an unary operator."
            ),
            Operator::Or => panic!(
                "FATAL GROAN ERROR | select::process_operation | OR operator is being treated as an unary operator."
            ),
        };
    }
    unary.clear();

    // apply the previous binary operator
    if let Some(op) = binary {
        if let Some(t) = tree {
            match op {
                Operator::And => Ok(Some(Box::from(Select::And(t, parsed)))),
                Operator::Or => Ok(Some(Box::from(Select::Or(t, parsed)))),
                Operator::Molecule => panic!(
                    "FATAL GROAN ERROR | select::process_operation | Molecule operator is being treated as a binary operator."
                ),
                Operator::Not => panic!(
                    "FATAL GROAN ERROR | select::process_operation | NOT operator is being treated as a binary operator."
                ),
            }
        } else {
            Err(SelectError::MissingArgument("".to_string()))
        }
    // or create a new tree
    } else {
        if tree.is_some() {
            return Err(SelectError::InvalidTokenParentheses("".to_string()));
        }
        Ok(Some(parsed))
    }
}

fn find_operator(string: &str, op_symbol: char, start: usize) -> Option<Operator> {
    if string.get(start + 1..start + 2) == Some(op_symbol.to_string().as_str()) {
        match op_symbol {
            '&' => Some(Operator::And),
            '|' => Some(Operator::Or),
            '@' => Some(Operator::Molecule),
            _ => None,
        }
    } else {
        None
    }
}

/// Check whether the number of '(' and ')' matches each other.
fn par_balanced(string: &str) -> bool {
    string.chars().fold(0, |acc, c| {
        if c == '(' {
            acc + 1
        } else if c == ')' {
            acc - 1
        } else {
            acc
        }
    }) == 0
}

/// Check whether the number of ' and " is even.
fn quotes_balanced(string: &str) -> bool {
    let single = string.chars().filter(|&c| c == '\'').count();
    let double = string.chars().filter(|&c| c == '"').count();

    single % 2 == 0 && double % 2 == 0
}

/// Macros @protein, @water, @ion, @dna, @rna are partly based on https://github.com/gromacs/gromacs/blob/main/share/top/residuetypes.dat
fn get_macros() -> HashMap<&'static str, &'static str> {
    let mut macros = HashMap::new();

    macros.insert(
        "@membrane",
        "(resname DAPC DBPC DFPC DGPC DIPC DLPC DNPC DOPC DPPC DRPC DTPC DVPC DXPC DYPC LPPC PAPC PEPC PGPC PIPC POPC PRPC PUPC DAPE DBPE DFPE DGPE DIPE DLPE DNPE DOPE DPPE DRPE DTPE DUPE DVPE DXPE DYPE LPPE PAPE PGPE PIPE POPE PQPE PRPE PUPE DAPS DBPS DFPS DGPS DIPS DLPS DNPS DOPS DPPS DRPS DTPS DUPS DVPS DXPS DYPS LPPS PAPS PGPS PIPS POPS PQPS PRPS PUPS DAPG DBPG DFPG DGPG DIPG DLPG DNPG DOPG DPPG DRPG DTPG DVPG DXPG DYPG JFPG JPPG LPPG OPPG PAPG PGPG PIPG POPG PRPG DAPA DBPA DFPA DGPA DIPA DLPA DNPA DOPA DPPA DRPA DTPA DVPA DXPA DYPA LPPA PAPA PGPA PIPA POPA PRPA PUPA DPP1 DPP2 DPPI PAPI PIPI POP1 POP2 POP3 POPI PUPI PVP1 PVP2 PVP3 PVPI PADG PIDG PODG PUDG PVDG TOG APC CPC IPC LPC OPC PPC TPC UPC VPC BNSM DBSM DPSM DXSM PGSM PNSM POSM PVSM XNSM DPCE DXCE PNCE XNCE DBG1 DPG1 DPG3 DPGS DXG1 DXG3 PNG1 PNG3 XNG1 XNG3 DFGG DFMG DPGG DPMG DPSG FPGG FPMG FPSG OPGG OPMG OPSG CHOA CHOL CHYO BOG DDM DPC EO5 SDS BOLA BOLB CDL0 CDL1 CDL2 CDL DBG3 ERGO HBHT HDPT HHOP HOPR ACA ACN BCA BCN LCA LCN PCA PCN UCA UCN XCA XCN RAMP REMP OANT)");
    macros.insert(
        "@protein",
        "(resname ABU ACE AIB ALA ARG ARGN ASN ASN1 ASP ASP1 ASPH ASPP ASH CT3 CYS CYS1 CYS2 CYSH DALA GLN GLU GLUH GLUP GLH GLY HIS HIS1 HISA HISB HISH HISD HISE HISP HSD HSE HSP HYP ILE LEU LSN LYS LYSN LYSH MELEU MET MEVAL NAC NME NHE NH2 PHE PHEH PHEU PHL PRO SER THR TRP TRPH TRPU TYR TYRH TYRU VAL PGLU HID HIE HIP LYP LYN CYN CYM CYX DAB ORN HYP NALA NGLY NSER NTHR NLEU NILE NVAL NASN NGLN NARG NHID NHIE NHIP NHISD NHISE NHISH NTRP NPHE NTYR NGLU NASP NLYS NORN NDAB NLYSN NPRO NHYP NCYS NCYS2 NMET NASPH NGLUH CALA CGLY CSER CTHR CLEU CILE CVAL CASN CGLN CARG CHID CHIE CHIP CHISD CHISE CHISH CTRP CPHE CTYR CGLU CASP CLYS CORN CDAB CLYSN CPRO CHYP CCYS CCYS2 CMET CASPH CGLUH)",
    );
    macros.insert(
        "@water",
        "(name W OW HW1 HW2 OH2 H1 H2 and resname SOL WAT HOH OHH TIP T3P T4P T5P T3H W TIP3 TIP4 SPC SPCE)",
    );
    macros.insert(
        "@ion",
        "(name NA NA+ CL CL- K K+ SOD CLA CA CA2+ MG ZN CU1 CU LI RB CS F BR I OH Cal CAL IB+ and resname ION NA NA+ CL CL- K K+ SOD CLA CA CA2+ MG ZN CU1 CU LI RB CS F BR I OH Cal CAL IB+)",
    );
    macros.insert(
        "@dna",
        "(resname DA DG DC DT DA5 DG5 DC5 DT5 DA3 DG3 DC3 DT3 DAN DGN DCN DTN)",
    );
    macros.insert(
        "@rna",
        "(resname A U C G RA RU RC RG RA5 RT5 RU5 RC5 RG5 RA3 RT3 RU3 RC3 RG3 RAN RTN RUN RCN RGN)",
    );
    // deprecated since v0.6.0
    //macros.insert("@hydrogen", "(name r'^[1-9]?[Hh].*')");

    macros
}

fn expand_macros(string: &mut String, macros: HashMap<&str, &str>) {
    for (m, expand) in macros {
        *string = string.replace(m, expand);
    }
}

fn find_parenthesis(query: &str, start: usize, end: usize) -> Option<usize> {
    let mut open = 0;
    let mut closed = 0;

    for (index, c) in query.char_indices().skip(start) {
        if c == '(' {
            open += 1;
        } else if c == ')' {
            closed += 1;
            if open == closed {
                return Some(index);
            }
        }

        if index > end {
            return None;
        }
    }

    None
}

/// Replace alphabetical keywords with their symbolic representations.
/// Ignores quote blocks.
fn replace_keywords(input: &str) -> String {
    let mut result = String::new();
    let mut input_chars = input.chars().peekable();
    let mut inside_quotes = false;

    while let Some(c) = input_chars.next() {
        if c == '\'' || c == '"' {
            inside_quotes = !inside_quotes;
            result.push(c);
            continue;
        }

        if inside_quotes {
            result.push(c);
            continue;
        }

        if c.is_alphabetic() {
            let keyword = get_keyword(&mut input_chars, c);
            let replaced_keyword = match keyword.as_str() {
                "and" => "&&",
                "or" => "||",
                "not" => "!",
                "to" => "-",
                _ => keyword.as_str(),
            };
            result.push_str(replaced_keyword);
        } else {
            result.push(c);
        }
    }

    result
}

fn get_keyword<I: Iterator<Item = char>>(
    iter: &mut std::iter::Peekable<I>,
    first_char: char,
) -> String {
    let mut keyword = String::new();
    keyword.push(first_char);

    while let Some(&c) = iter.peek() {
        if c.is_alphabetic() {
            keyword.push(iter.next().unwrap());
        } else {
            break;
        }
    }

    keyword
}

/// Split a string by whitespace while keeping the items enclosed in ' or " together.
fn split_with_quotes(string: &str) -> Vec<String> {
    let mut result = vec![String::new()];
    let mut inside = false;
    let mut block = 0;
    let mut regex = false;

    let mut iterator = string.chars().peekable();

    while let Some(c) = iterator.next() {
        if c == 'r' && !inside {
            if let Some('\'') = iterator.peek() {
                regex = true;
                inside = !inside;
                result[block].push('r');
                result[block].push('\'');
                iterator.next();
                continue;
            }
        }

        if c == '\'' || c == '"' {
            inside = !inside;
            if regex {
                result[block].push(c);
                regex = false;
            }
            continue;
        }

        if c.is_whitespace() && !inside {
            result.push(String::new());
            block += 1;
            continue;
        }

        result[block].push(c);
    }

    // filter out empty strings and trim all strings
    result
        .into_iter()
        .filter(|s| !s.trim().is_empty())
        .collect()
}

/// Collect words from the query and conver them to the `Name` enum.
fn collect_words(token: &[String]) -> Result<Vec<Name>, SelectError> {
    let names: Result<Vec<Name>, SelectError> = token.iter().map(|s| Name::new(s)).collect();

    names
}

fn parse_token(string: &str) -> Result<Select, SelectError> {
    if string.trim().is_empty() {
        return Err(SelectError::MissingArgument("".to_string()));
    }

    let token = split_with_quotes(string);

    match token[0].as_str() {
        "resname" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            Ok(Select::ResidueName(collect_words(&token[1..])?))
        }
        "name" | "atomname" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            Ok(Select::AtomName(collect_words(&token[1..])?))
        }
        "resid" | "resnum" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            let range = numbers::parse_numbers(&token[1..])?;
            Ok(Select::ResidueNumber(fix_ranges(range)))
        }
        "serial" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            let range = numbers::parse_numbers(&token[1..])?;
            Ok(Select::GmxAtomNumber(fix_ranges(range)))
        }

        "atomid" | "atomnum" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            let range = numbers::parse_numbers(&token[1..])?;
            Ok(Select::AtomNumber(fix_ranges(range)))
        }

        "chain" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            // collect tokens into the Select::Chain enum and check that each token is only one character long
            token
                .iter()
                .skip(1)
                .try_fold(Vec::new(), |mut acc, t| {
                    if t.len() == 1 {
                        acc.push(t.chars().next().unwrap());
                        Ok(acc)
                    } else {
                        Err(SelectError::InvalidChainId("".to_string()))
                    }
                })
                .map(Select::Chain)
        }

        "group" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            Ok(Select::GroupName(collect_words(&token[1..])?))
        }

        "label" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            Ok(Select::LabeledAtom(collect_words(&token[1..])?))
        }

        "element" if token.len() >= 2 && token[1] == "name" => {
            if token.len() <= 2 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            Ok(Select::ElementName(collect_words(&token[2..])?))
        }

        "elname" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            Ok(Select::ElementName(collect_words(&token[1..])?))
        }

        "element" if token.len() >= 2 && token[1] == "symbol" => {
            if token.len() <= 2 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            Ok(Select::ElementSymbol(collect_words(&token[2..])?))
        }

        "elsymbol" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            Ok(Select::ElementSymbol(collect_words(&token[1..])?))
        }

        // it is not necessary to provide group identifier for groups
        _ => Ok(Select::GroupName(collect_words(&token)?)),
    }
}

fn fix_ranges(mut ranges: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    if ranges.is_empty() {
        return ranges;
    }

    // sort the ranges in ascending order
    ranges.sort_unstable();

    let mut merged_indices = Vec::new();
    let mut current_start = usize::MAX;
    let mut current_end = 0usize;

    for (start, end) in &ranges {
        // start must not not larger than end
        if *start > *end {
            continue;
        }

        // current range does not overlap with the previous one nor is adjacent to it
        if *start > current_end + 1 || (current_end == 0usize && current_start != 0usize) {
            if current_start != usize::MAX {
                merged_indices.push((current_start, current_end));
            }
            current_start = *start;
            current_end = *end;
        // current range overlaps
        } else if *end > current_end {
            current_end = *end;
        }
    }

    if current_start != usize::MAX {
        // add the last merged range to the result if it exists
        merged_indices.push((current_start, current_end));
    }

    merged_indices
}

/******************************/
/*       FEATURE: SERDE       */
/******************************/

#[cfg(feature = "serde")]
mod serde {
    use serde::{
        de::{self, Visitor},
        ser::{Serialize, Serializer},
        Deserialize, Deserializer,
    };
    use std::fmt;

    use crate::select::Select;

    impl Serialize for Select {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer,
        {
            let query_str = self.to_string();
            serializer.serialize_str(&query_str)
        }
    }

    impl<'de> Deserialize<'de> for Select {
        fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: Deserializer<'de>,
        {
            struct SelectVisitor;

            impl<'de> Visitor<'de> for SelectVisitor {
                type Value = Select;

                fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                    formatter.write_str("a string representing a Groan Selection Language query")
                }

                fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
                where
                    E: de::Error,
                {
                    match Select::parse_query(value) {
                        Ok(select_box) => Ok(*select_box),
                        Err(err) => Err(E::custom(err.to_string())),
                    }
                }
            }

            deserializer.deserialize_str(SelectVisitor)
        }
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod pass_tests {
    use super::*;

    macro_rules! parsing_success {
        ($name:ident, $expression:expr, $expected:expr) => {
            #[test]
            fn $name() {
                let query = $expression;

                match Select::parse_query(query) {
                    Ok(x) => assert_eq!(*x, $expected),
                    Err(e) => panic!("Parsing failed, returning {:?}", e),
                }
            }
        };
    }

    macro_rules! parsing_success_string {
        ($name:ident, $expression:expr, $expected:expr) => {
            #[test]
            fn $name() {
                let query = $expression;

                match Select::parse_query(query) {
                    Ok(x) => {
                        let string1 = format!("{:?}", *x);
                        assert_eq!(string1, $expected)
                    }
                    Err(e) => panic!("Parsing failed, returning {:?}", e),
                }
            }
        };
    }

    parsing_success!(
        simple_resname,
        "resname LYS",
        Select::ResidueName(vec![Name::new("LYS").unwrap()])
    );
    parsing_success!(
        multiple_resname,
        "resname  LYS   CYS   LEU  POPE ",
        Select::ResidueName(vec![
            Name::new("LYS").unwrap(),
            Name::new("CYS").unwrap(),
            Name::new("LEU").unwrap(),
            Name::new("POPE").unwrap()
        ])
    );
    parsing_success!(
        simple_atomname,
        "name BB",
        Select::AtomName(vec![Name::new("BB").unwrap()])
    );
    parsing_success!(
        multiple_atomname,
        "  name   BB SC1   PO4   C1 C2B",
        Select::AtomName(vec![
            Name::new("BB").unwrap(),
            Name::new("SC1").unwrap(),
            Name::new("PO4").unwrap(),
            Name::new("C1").unwrap(),
            Name::new("C2B").unwrap()
        ])
    );

    parsing_success!(
        simple_resid,
        "resid 15",
        Select::ResidueNumber(vec![(15, 15)])
    );
    parsing_success!(
        multiple_resid,
        "resid 15 16 22 24 48   ",
        Select::ResidueNumber(vec![(15, 16), (22, 22), (24, 24), (48, 48)])
    );
    parsing_success!(
        range_resid_1,
        " resid 4 to 11",
        Select::ResidueNumber(vec![(4, 11)])
    );
    parsing_success!(
        range_resid_2,
        "resid 4 - 11",
        Select::ResidueNumber(vec![(4, 11)])
    );
    parsing_success!(
        range_resid_3,
        "resid 4-11",
        Select::ResidueNumber(vec![(4, 11)])
    );
    parsing_success!(
        range_resid_4,
        "  resnum 4- 11",
        Select::ResidueNumber(vec![(4, 11)])
    );
    parsing_success!(
        range_resid_5,
        "resid 4 -11",
        Select::ResidueNumber(vec![(4, 11)])
    );
    parsing_success!(
        complex_resid,
        " resnum 4 8   to 11 to 14   20  21",
        Select::ResidueNumber(vec![(4, 4), (8, 14), (20, 21)])
    );

    parsing_success!(
        simple_serial,
        "serial 15",
        Select::GmxAtomNumber(vec![(15, 15)])
    );
    parsing_success!(
        multiple_serial,
        "serial 15 16 22 24 48",
        Select::GmxAtomNumber(vec![(15, 16), (22, 22), (24, 24), (48, 48)])
    );
    parsing_success!(
        range_serial_1,
        "serial 4-11   ",
        Select::GmxAtomNumber(vec![(4, 11)])
    );
    parsing_success!(
        range_serial_2,
        "  serial 4 to11   ",
        Select::GmxAtomNumber(vec![(4, 11)])
    );
    parsing_success!(
        range_atomid_1,
        "  atomid 4 to 11   ",
        Select::AtomNumber(vec![(4, 11)])
    );
    parsing_success!(
        range_atomid_2,
        "atomnum 4 - 11",
        Select::AtomNumber(vec![(4, 11)])
    );

    parsing_success!(
        range_atomid_3,
        "   atomid 4- 11",
        Select::AtomNumber(vec![(4, 11)])
    );
    parsing_success!(
        range_atomid_5,
        "atomnum 4 -11",
        Select::AtomNumber(vec![(4, 11)])
    );

    parsing_success!(
        range_atomid_6,
        "  atomid 4to 11   ",
        Select::AtomNumber(vec![(4, 11)])
    );

    parsing_success!(
        open_ended_range_1,
        "serial  <   44",
        Select::GmxAtomNumber(vec![(1, 43)])
    );

    parsing_success!(
        open_ended_range_2,
        "atomid <44",
        Select::AtomNumber(vec![(1, 43)])
    );

    parsing_success!(
        open_ended_range_3,
        "  resid    <    30   ",
        Select::ResidueNumber(vec![(1, 29)])
    );

    parsing_success!(
        open_ended_range_4,
        "serial  <=44",
        Select::GmxAtomNumber(vec![(1, 44)])
    );

    parsing_success!(
        open_ended_range_5,
        "atomid  <=  44",
        Select::AtomNumber(vec![(1, 44)])
    );

    parsing_success!(
        open_ended_range_6,
        " resnum  >44",
        Select::ResidueNumber(vec![(45, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_7,
        "serial  >  44",
        Select::GmxAtomNumber(vec![(45, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_8,
        "atomnum  >=44",
        Select::AtomNumber(vec![(44, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_9,
        "resid  >=  44   ",
        Select::ResidueNumber(vec![(44, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_10,
        "serial 1 4 > 50 7-11",
        Select::GmxAtomNumber(vec![(1, 1), (4, 4), (7, 11), (51, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_11,
        "resid   1   4 >50 7 - 11",
        Select::ResidueNumber(vec![(1, 1), (4, 4), (7, 11), (51, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_12,
        "atomid 1 4 >= 50 7-11",
        Select::AtomNumber(vec![(1, 1), (4, 4), (7, 11), (50, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_13,
        "serial 1 4 >=50 7-11",
        Select::GmxAtomNumber(vec![(1, 1), (4, 4), (7, 11), (50, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_14,
        "serial 1 4>=50 7-11",
        Select::GmxAtomNumber(vec![(1, 1), (4, 4), (7, 11), (50, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_15,
        "resid 4>50",
        Select::ResidueNumber(vec![(4, 4), (51, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_16,
        "atomid 50<20",
        Select::AtomNumber(vec![(1, 19), (50, 50)])
    );

    parsing_success!(
        open_ended_range_17,
        "serial 50<=20",
        Select::GmxAtomNumber(vec![(1, 20), (50, 50)])
    );

    parsing_success!(
        open_ended_range_18,
        "serial 50<=20 - 25",
        Select::GmxAtomNumber(vec![(1, 25), (50, 50)])
    );

    parsing_success!(
        open_ended_range_19,
        "serial <= 20 -40",
        Select::GmxAtomNumber(vec![(1, 40)])
    );

    parsing_success!(
        open_ended_range_20,
        "resnum 5> 20",
        Select::ResidueNumber(vec![(5, 5), (21, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_21,
        "serial 5>= 20",
        Select::GmxAtomNumber(vec![(5, 5), (20, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_22,
        "serial 24-37<=10 55",
        Select::GmxAtomNumber(vec![(1, 10), (24, 37), (55, 55)])
    );

    parsing_success!(
        open_ended_range_23,
        "serial 24-37<=10-13to17>88",
        Select::GmxAtomNumber(vec![(1, 17), (24, 37), (89, usize::MAX)])
    );

    parsing_success!(
        open_ended_range_24,
        "serial < 0",
        Select::GmxAtomNumber(vec![])
    );

    parsing_success!(
        open_ended_range_25,
        "serial < 1",
        Select::GmxAtomNumber(vec![])
    );

    parsing_success!(
        open_ended_range_26,
        "serial <= 1",
        Select::GmxAtomNumber(vec![(1, 1)])
    );

    parsing_success!(
        open_ended_range_negation_1,
        "not resid >20",
        Select::Not(Box::from(Select::ResidueNumber(vec![(21, usize::MAX)])))
    );

    parsing_success!(
        complex_serial,
        "serial 4 8   to 11 to 14   20  21",
        Select::GmxAtomNumber(vec![(4, 4), (8, 14), (20, 21)])
    );

    parsing_success!(
        simple_not_resname,
        "!resname LYS   SER",
        Select::Not(Box::from(Select::ResidueName(vec![
            Name::new("LYS").unwrap(),
            Name::new("SER").unwrap()
        ])))
    );
    parsing_success!(
        simple_not_name,
        " not name   BB SC1",
        Select::Not(Box::from(Select::AtomName(vec![
            Name::new("BB").unwrap(),
            Name::new("SC1").unwrap()
        ])))
    );
    parsing_success!(
        simple_not_resid,
        "not   resid 5 6 to 9",
        Select::Not(Box::from(Select::ResidueNumber(vec![(5, 9)])))
    );
    parsing_success!(
        simple_not_serial,
        "!   serial 15-18",
        Select::Not(Box::from(Select::GmxAtomNumber(vec![(15, 18)])))
    );
    parsing_success!(
        simple_not_explicit_group,
        "! group Protein Membrane",
        Select::Not(Box::from(Select::GroupName(vec![
            Name::new("Protein").unwrap(),
            Name::new("Membrane").unwrap()
        ])))
    );
    parsing_success!(
        simple_not_explicit,
        "! Protein Membrane",
        Select::Not(Box::from(Select::GroupName(vec![
            Name::new("Protein").unwrap(),
            Name::new("Membrane").unwrap()
        ])))
    );
    parsing_success!(
        simple_not_not,
        " not ! name BB SC1",
        Select::Not(Box::from(Select::Not(Box::from(Select::AtomName(vec![
            Name::new("BB").unwrap(),
            Name::new("SC1").unwrap()
        ])))))
    );
    parsing_success!(
        simple_not_not_not,
        "!!! name BB SC1",
        Select::Not(Box::from(Select::Not(Box::from(Select::Not(Box::from(
            Select::AtomName(vec![Name::new("BB").unwrap(), Name::new("SC1").unwrap()])
        ))))))
    );

    parsing_success!(
        not_parentheses_1,
        "! (name BB or resname LYS)",
        Select::Not(Box::from(Select::Or(
            Box::from(Select::AtomName(vec![Name::new("BB").unwrap()])),
            Box::from(Select::ResidueName(vec![Name::new("LYS").unwrap()]))
        )))
    );

    parsing_success!(
        not_parentheses_2,
        "not(name BB or resname LYS) ||serial 1-5",
        Select::Or(
            Box::from(Select::Not(Box::from(Select::Or(
                Box::from(Select::AtomName(vec![Name::new("BB").unwrap()])),
                Box::from(Select::ResidueName(vec![Name::new("LYS").unwrap()]))
            )))),
            Box::from(Select::GmxAtomNumber(vec![(1, 5)]))
        )
    );

    parsing_success!(
        not_parentheses_3,
        "(name BB or resname LYS) || ! serial 1-5",
        Select::Or(
            Box::from(Select::Or(
                Box::from(Select::AtomName(vec![Name::new("BB").unwrap()])),
                Box::from(Select::ResidueName(vec![Name::new("LYS").unwrap()]))
            )),
            Box::from(Select::Not(Box::from(Select::GmxAtomNumber(vec![(1, 5)]))))
        )
    );

    parsing_success!(
        not_parentheses_4,
        "not(name BB or not resname LYS) ||serial 1-5",
        Select::Or(
            Box::from(Select::Not(Box::from(Select::Or(
                Box::from(Select::AtomName(vec![Name::new("BB").unwrap()])),
                Box::from(Select::Not(Box::from(Select::ResidueName(vec![
                    Name::new("LYS").unwrap()
                ]))))
            )))),
            Box::from(Select::GmxAtomNumber(vec![(1, 5)]))
        )
    );

    parsing_success!(
        advanced_ranges_1,
        "atomid 4 to 6 to 12",
        Select::AtomNumber(vec![(4, 12)])
    );
    parsing_success!(
        advanced_ranges_2,
        "serial 1 4 to 6 to 12",
        Select::GmxAtomNumber(vec![(1, 1), (4, 12)])
    );
    parsing_success!(
        advanced_ranges_3,
        "atomnum 4 to 6 to 12 15 1",
        Select::AtomNumber(vec![(1, 1), (4, 12), (15, 15)])
    );
    parsing_success!(
        advanced_ranges_4,
        "serial 1 4 - 6 - 12 15",
        Select::GmxAtomNumber(vec![(1, 1), (4, 12), (15, 15)])
    );
    parsing_success!(
        advanced_ranges_5,
        "atomid 4-6 -12 15 1",
        Select::AtomNumber(vec![(1, 1), (4, 12), (15, 15)])
    );
    parsing_success!(
        advanced_ranges_6,
        "serial 15 1 4- 6 -12",
        Select::GmxAtomNumber(vec![(1, 1), (4, 12), (15, 15)])
    );
    parsing_success!(
        advanced_ranges_7,
        "atomnum 1 4 - 6-12 15",
        Select::AtomNumber(vec![(1, 1), (4, 12), (15, 15)])
    );
    parsing_success!(
        advanced_ranges_8,
        "atomnum 1 15 4- 6-12",
        Select::AtomNumber(vec![(1, 1), (4, 12), (15, 15)])
    );
    parsing_success!(
        advanced_ranges_9,
        "atomnum 1 15 4-6-12-18",
        Select::AtomNumber(vec![(1, 1), (4, 18)])
    );
    parsing_success!(
        advanced_ranges_10,
        "resid 1 15 4-6-12-18 23",
        Select::ResidueNumber(vec![(1, 1), (4, 18), (23, 23)])
    );

    parsing_success!(chain_simple, "chain B", Select::Chain(vec!['B']));

    parsing_success!(
        chain_multi,
        "chain B C E F",
        Select::Chain(vec!['B', 'C', 'E', 'F'])
    );

    parsing_success!(
        chain_combined,
        "(serial 1 to 3 and chain A) or chain C D",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::GmxAtomNumber(vec![(1, 3)])),
                Box::from(Select::Chain(vec!['A']))
            )),
            Box::from(Select::Chain(vec!['C', 'D']))
        )
    );

    parsing_success!(
        simple_group,
        "Protein",
        Select::GroupName(vec![Name::new("Protein").unwrap()])
    );
    parsing_success!(
        simple_explicit_group,
        "group Protein",
        Select::GroupName(vec![Name::new("Protein").unwrap()])
    );
    parsing_success!(
        multiword_group_1,
        "  'Protein Membrane ION'  ",
        Select::GroupName(vec![Name::new("Protein Membrane ION").unwrap()])
    );
    parsing_success!(
        multiword_group_2,
        "\"Protein   Membrane  ION  \"",
        Select::GroupName(vec![Name::new("Protein   Membrane  ION  ").unwrap()])
    );
    parsing_success!(
        multiword_explicit_group,
        "  group  '  Protein Membrane ION'  ",
        Select::GroupName(vec![Name::new("  Protein Membrane ION").unwrap()])
    );
    parsing_success!(
        multiple_groups,
        "Protein Membrane ION",
        Select::GroupName(vec![
            Name::new("Protein").unwrap(),
            Name::new("Membrane").unwrap(),
            Name::new("ION").unwrap()
        ])
    );
    parsing_success!(
        multiple_groups_binary,
        "Protein and Membrane or ION",
        Select::Or(
            Box::new(Select::And(
                Box::new(Select::GroupName(vec![Name::new("Protein").unwrap()])),
                Box::new(Select::GroupName(vec![Name::new("Membrane").unwrap()]))
            )),
            Box::new(Select::GroupName(vec![Name::new("ION").unwrap()]))
        )
    );
    parsing_success!(
        multiple_explicit_groups,
        "group Protein Membrane ION",
        Select::GroupName(vec![
            Name::new("Protein").unwrap(),
            Name::new("Membrane").unwrap(),
            Name::new("ION").unwrap()
        ])
    );
    parsing_success!(
        multiple_multiword_groups,
        "Protein 'Membrane Only POPC  ' ' ION with  Water'",
        Select::GroupName(vec![
            Name::new("Protein").unwrap(),
            Name::new("Membrane Only POPC  ").unwrap(),
            Name::new(" ION with  Water").unwrap()
        ])
    );
    parsing_success!(
        multiple_multiword_explicit_groups,
        "group Protein 'Membrane Only POPC  ' ' ION with  Water'",
        Select::GroupName(vec![
            Name::new("Protein").unwrap(),
            Name::new("Membrane Only POPC  ").unwrap(),
            Name::new(" ION with  Water").unwrap()
        ])
    );

    parsing_success!(
        simple_label,
        "label  atom1   atom2   atomX",
        Select::LabeledAtom(vec![
            Name::new("atom1").unwrap(),
            Name::new("atom2").unwrap(),
            Name::new("atomX").unwrap()
        ])
    );
    parsing_success!(
        multiword_label_1,
        "label   'this is a labeled Atom'",
        Select::LabeledAtom(vec![Name::new("this is a labeled Atom").unwrap()])
    );
    parsing_success!(
        multiword_label_2,
        "label \"this is a labeled Atom\"",
        Select::LabeledAtom(vec![Name::new("this is a labeled Atom").unwrap()])
    );
    parsing_success!(
        multiple_multiword_labels,
        "label 'Atom 1  ' \"  another atom\"   'atom X'",
        Select::LabeledAtom(vec![
            Name::new("Atom 1  ").unwrap(),
            Name::new("  another atom").unwrap(),
            Name::new("atom X").unwrap()
        ])
    );
    parsing_success!(
        complex_label,
        "label 'Atom 1 ' MyAtom   &&!label 'atom X'",
        Select::And(
            Box::new(Select::LabeledAtom(vec![
                Name::new("Atom 1 ").unwrap(),
                Name::new("MyAtom").unwrap()
            ])),
            Box::new(Select::Not(Box::new(Select::LabeledAtom(vec![Name::new(
                "atom X"
            )
            .unwrap()]))))
        )
    );
    parsing_success!(
        label_regex,
        "label r'^A' other_atom",
        Select::LabeledAtom(vec![
            Name::new("r'^A'").unwrap(),
            Name::new("other_atom").unwrap()
        ])
    );

    parsing_success!(
        element_symbol_1,
        "element symbol C",
        Select::ElementSymbol(vec![Name::new("C").unwrap()])
    );
    parsing_success!(
        element_symbol_2,
        "elsymbol C",
        Select::ElementSymbol(vec![Name::new("C").unwrap()])
    );
    parsing_success!(
        element_symbol_3,
        "element symbol C Na K Br",
        Select::ElementSymbol(vec![
            Name::new("C").unwrap(),
            Name::new("Na").unwrap(),
            Name::new("K").unwrap(),
            Name::new("Br").unwrap()
        ])
    );
    parsing_success!(
        element_symbol_4,
        "elsymbol C Na K Br",
        Select::ElementSymbol(vec![
            Name::new("C").unwrap(),
            Name::new("Na").unwrap(),
            Name::new("K").unwrap(),
            Name::new("Br").unwrap()
        ])
    );
    parsing_success!(
        element_symbol_regex,
        "elsymbol r'^N' C K",
        Select::ElementSymbol(vec![
            Name::new("r'^N'").unwrap(),
            Name::new("C").unwrap(),
            Name::new("K").unwrap()
        ])
    );
    parsing_success!(
        element_symbol_complex_1,
        "element symbol C K and serial 4 to 8",
        Select::And(
            Box::new(Select::ElementSymbol(vec![
                Name::new("C").unwrap(),
                Name::new("K").unwrap()
            ])),
            Box::new(Select::GmxAtomNumber(vec![(4, 8)])),
        )
    );
    parsing_success!(
        element_symbol_complex_2,
        "resname LYS SER || elsymbol C K",
        Select::Or(
            Box::new(Select::ResidueName(vec![
                Name::new("LYS").unwrap(),
                Name::new("SER").unwrap()
            ])),
            Box::new(Select::ElementSymbol(vec![
                Name::new("C").unwrap(),
                Name::new("K").unwrap()
            ])),
        )
    );
    parsing_success!(
        element_name_1,
        "element name carbon",
        Select::ElementName(vec![Name::new("carbon").unwrap()])
    );
    parsing_success!(
        element_name_2,
        "elname carbon",
        Select::ElementName(vec![Name::new("carbon").unwrap()])
    );
    parsing_success!(
        element_name_multiword,
        "element name 'carbon number 2' hydrogen",
        Select::ElementName(vec![
            Name::new("carbon number 2").unwrap(),
            Name::new("hydrogen").unwrap()
        ])
    );
    parsing_success!(
        element_name_3,
        "element name carbon sodium potassium bromine",
        Select::ElementName(vec![
            Name::new("carbon").unwrap(),
            Name::new("sodium").unwrap(),
            Name::new("potassium").unwrap(),
            Name::new("bromine").unwrap()
        ])
    );
    parsing_success!(
        element_name_4,
        "elname carbon sodium potassium bromine",
        Select::ElementName(vec![
            Name::new("carbon").unwrap(),
            Name::new("sodium").unwrap(),
            Name::new("potassium").unwrap(),
            Name::new("bromine").unwrap()
        ])
    );
    parsing_success!(
        element_name_regex,
        "elname r'^ni' carbon potassium",
        Select::ElementName(vec![
            Name::new("r'^ni'").unwrap(),
            Name::new("carbon").unwrap(),
            Name::new("potassium").unwrap()
        ])
    );
    parsing_success!(
        element_name_complex_1,
        "element name carbon potassium and serial 4 to 8",
        Select::And(
            Box::new(Select::ElementName(vec![
                Name::new("carbon").unwrap(),
                Name::new("potassium").unwrap()
            ])),
            Box::new(Select::GmxAtomNumber(vec![(4, 8)])),
        )
    );
    parsing_success!(
        element_name_complex_2,
        "resname LYS SER || elname carbon potassium",
        Select::Or(
            Box::new(Select::ResidueName(vec![
                Name::new("LYS").unwrap(),
                Name::new("SER").unwrap()
            ])),
            Box::new(Select::ElementName(vec![
                Name::new("carbon").unwrap(),
                Name::new("potassium").unwrap()
            ])),
        )
    );
    parsing_success!(
        not_element,
        "element C",
        Select::GroupName(vec![Name::new("element").unwrap(), Name::new("C").unwrap()])
    );

    parsing_success!(
        molwith_1,
        "molecule with serial 1 to 12 17",
        Select::Molecule(Box::new(Select::GmxAtomNumber(vec![(1, 12), (17, 17)])))
    );
    parsing_success!(
        molwith_2,
        "mol with serial 1 to 12 17",
        Select::Molecule(Box::new(Select::GmxAtomNumber(vec![(1, 12), (17, 17)])))
    );
    parsing_success!(
        molwith_3,
        "molecule with resname LYS && resid 1 to 3",
        Select::And(
            Box::new(Select::Molecule(Box::new(Select::ResidueName(vec![
                Name::new("LYS").unwrap()
            ])))),
            Box::new(Select::ResidueNumber(vec![(1, 3)])),
        )
    );
    parsing_success!(
        molwith_4,
        "molecule with (resname LYS && resid 1 to 3)",
        Select::Molecule(Box::new(Select::And(
            Box::new(Select::ResidueName(vec![Name::new("LYS").unwrap()])),
            Box::new(Select::ResidueNumber(vec![(1, 3)])),
        )))
    );
    parsing_success!(
        molwith_5,
        "molecule with (mol with name P and resid 50 to 76)",
        Select::Molecule(Box::new(Select::And(
            Box::new(Select::Molecule(Box::new(Select::AtomName(vec![
                Name::new("P").unwrap()
            ])))),
            Box::new(Select::ResidueNumber(vec![(50, 76)])),
        )))
    );
    parsing_success!(
        molwith_6,
        "molecule      with    serial 1 to 12 17",
        Select::Molecule(Box::new(Select::GmxAtomNumber(vec![(1, 12), (17, 17)])))
    );
    parsing_success!(
        molwith_7,
        "moleculewith serial 1 to 12 17",
        Select::Molecule(Box::new(Select::GmxAtomNumber(vec![(1, 12), (17, 17)])))
    );
    parsing_success!(
        molwith_8,
        "mol          with serial 1 to 12 17",
        Select::Molecule(Box::new(Select::GmxAtomNumber(vec![(1, 12), (17, 17)])))
    );
    parsing_success!(
        molwith_9,
        "molwith serial 1 to 12 17",
        Select::Molecule(Box::new(Select::GmxAtomNumber(vec![(1, 12), (17, 17)])))
    );
    parsing_success!(
        molwith_10,
        "molecule with 'molecule with Group' Membrane",
        Select::Molecule(Box::new(Select::GroupName(vec![
            Name::new("molecule with Group").unwrap(),
            Name::new("Membrane").unwrap()
        ])))
    );
    parsing_success!(
        not_molwith,
        "'molecule with serial 17' Membrane",
        Select::GroupName(vec![
            Name::new("molecule with serial 17").unwrap(),
            Name::new("Membrane").unwrap()
        ])
    );
    parsing_success!(
        molwith_label_1,
        "molecule with label MyAtom",
        Select::Molecule(Box::new(Select::LabeledAtom(vec![
            Name::new("MyAtom").unwrap()
        ])))
    );

    parsing_success!(
        molwith_label_2,
        "molwith  label MyAtom  MyAtom2 'Interesting atom '",
        Select::Molecule(Box::new(Select::LabeledAtom(vec![
            Name::new("MyAtom").unwrap(),
            Name::new("MyAtom2").unwrap(),
            Name::new("Interesting atom ").unwrap(),
        ])))
    );

    parsing_success!(
        hyphen_group,
        "Protein-No1",
        Select::GroupName(vec![Name::new("Protein-No1").unwrap()])
    );
    parsing_success!(
        hyphen_multiword_group,
        "'Protein - No1'",
        Select::GroupName(vec![Name::new("Protein - No1").unwrap()])
    );
    parsing_success!(
        group_with_unfortunate_name,
        "group resname",
        Select::GroupName(vec![Name::new("resname").unwrap()])
    );
    parsing_success!(
        group_with_very_unfortunate_name,
        "group group",
        Select::GroupName(vec![Name::new("group").unwrap()])
    );

    parsing_success!(
        simple_parentheses,
        "(resname LYS)",
        Select::ResidueName(vec![Name::new("LYS").unwrap()])
    );

    // residue names
    parsing_success!(
        resnames_or_1,
        "resname POPE  'LYS' LEU or resname POPG",
        Select::Or(
            Box::from(Select::ResidueName(vec![
                Name::new("POPE").unwrap(),
                Name::new("LYS").unwrap(),
                Name::new("LEU").unwrap(),
            ])),
            Box::from(Select::ResidueName(vec![Name::new("POPG").unwrap()]))
        )
    );

    parsing_success!(
        resnames_or_2,
        "resname POPE   LYS LEU ||resname POPG",
        Select::Or(
            Box::from(Select::ResidueName(vec![
                Name::new("POPE").unwrap(),
                Name::new("LYS").unwrap(),
                Name::new("LEU").unwrap(),
            ])),
            Box::from(Select::ResidueName(vec![Name::new("POPG").unwrap()]))
        )
    );

    parsing_success!(
        resnames_and_1,
        "resname 'POPE' LYS LEU and resname 'POPG'",
        Select::And(
            Box::from(Select::ResidueName(vec![
                Name::new("POPE").unwrap(),
                Name::new("LYS").unwrap(),
                Name::new("LEU").unwrap(),
            ])),
            Box::from(Select::ResidueName(vec![Name::new("POPG").unwrap()]))
        )
    );

    parsing_success!(
        resnames_and_2,
        "resname POPE LYS LEU&& resname POPG",
        Select::And(
            Box::from(Select::ResidueName(vec![
                Name::new("POPE").unwrap(),
                Name::new("LYS").unwrap(),
                Name::new("LEU").unwrap(),
            ])),
            Box::from(Select::ResidueName(vec![Name::new("POPG").unwrap()]))
        )
    );

    parsing_success!(
        resnames_complex,
        "resname POPE 'LYS'   LEU && resname POPG ||resname LYS",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![
                    Name::new("POPE").unwrap(),
                    Name::new("LYS").unwrap(),
                    Name::new("LEU").unwrap(),
                ])),
                Box::from(Select::ResidueName(vec![Name::new("POPG").unwrap()]))
            )),
            Box::from(Select::ResidueName(vec![Name::new("LYS").unwrap()]))
        )
    );

    parsing_success!(
        resnames_complex_par_1,
        "((resname POPE   LYS LEU&& resname POPG)) or resname LYS",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![
                    Name::new("POPE").unwrap(),
                    Name::new("LYS").unwrap(),
                    Name::new("LEU").unwrap(),
                ])),
                Box::from(Select::ResidueName(vec![Name::new("POPG").unwrap()]))
            )),
            Box::from(Select::ResidueName(vec![Name::new("LYS").unwrap()]))
        )
    );

    parsing_success!(
        resnames_complex_par_2,
        "resname 'POPE' LYS LEU&& (resname POPG   ||   resname LYS)",
        Select::And(
            Box::from(Select::ResidueName(vec![
                Name::new("POPE").unwrap(),
                Name::new("LYS").unwrap(),
                Name::new("LEU").unwrap(),
            ])),
            Box::from(Select::Or(
                Box::from(Select::ResidueName(vec![Name::new("POPG").unwrap()])),
                Box::from(Select::ResidueName(vec![Name::new("LYS").unwrap()]))
            )),
        )
    );

    // atom names
    parsing_success!(
        names_or_1,
        "name BB  'PO4' D2A or name W",
        Select::Or(
            Box::from(Select::AtomName(vec![
                Name::new("BB").unwrap(),
                Name::new("PO4").unwrap(),
                Name::new("D2A").unwrap(),
            ])),
            Box::from(Select::AtomName(vec![Name::new("W").unwrap()]))
        )
    );

    parsing_success!(
        names_or_2,
        "name BB   PO4 D2A ||atomname W",
        Select::Or(
            Box::from(Select::AtomName(vec![
                Name::new("BB").unwrap(),
                Name::new("PO4").unwrap(),
                Name::new("D2A").unwrap(),
            ])),
            Box::from(Select::AtomName(vec![Name::new("W").unwrap()]))
        )
    );

    parsing_success!(
        names_and_1,
        "atomname 'BB' PO4 D2A and name 'W'",
        Select::And(
            Box::from(Select::AtomName(vec![
                Name::new("BB").unwrap(),
                Name::new("PO4").unwrap(),
                Name::new("D2A").unwrap(),
            ])),
            Box::from(Select::AtomName(vec![Name::new("W").unwrap()]))
        )
    );

    parsing_success!(
        names_and_2,
        "atomname BB PO4 D2A&& atomname W",
        Select::And(
            Box::from(Select::AtomName(vec![
                Name::new("BB").unwrap(),
                Name::new("PO4").unwrap(),
                Name::new("D2A").unwrap(),
            ])),
            Box::from(Select::AtomName(vec![Name::new("W").unwrap()]))
        )
    );

    parsing_success!(
        names_complex,
        "atomname BB 'PO4'   D2A && name W ||name PO4",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::AtomName(vec![
                    Name::new("BB").unwrap(),
                    Name::new("PO4").unwrap(),
                    Name::new("D2A").unwrap(),
                ])),
                Box::from(Select::AtomName(vec![Name::new("W").unwrap()]))
            )),
            Box::from(Select::AtomName(vec![Name::new("PO4").unwrap()]))
        )
    );

    parsing_success!(
        names_complex_par_1,
        "((name BB   PO4 D2A&& name W)) or name PO4",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::AtomName(vec![
                    Name::new("BB").unwrap(),
                    Name::new("PO4").unwrap(),
                    Name::new("D2A").unwrap(),
                ])),
                Box::from(Select::AtomName(vec![Name::new("W").unwrap()]))
            )),
            Box::from(Select::AtomName(vec![Name::new("PO4").unwrap()]))
        )
    );

    parsing_success!(
        names_complex_par_2,
        "atomname 'BB' PO4 D2A&& (name W   ||   atomname PO4)",
        Select::And(
            Box::from(Select::AtomName(vec![
                Name::new("BB").unwrap(),
                Name::new("PO4").unwrap(),
                Name::new("D2A").unwrap(),
            ])),
            Box::from(Select::Or(
                Box::from(Select::AtomName(vec![Name::new("W").unwrap()])),
                Box::from(Select::AtomName(vec![Name::new("PO4").unwrap()]))
            )),
        )
    );

    // residue numbers
    parsing_success!(
        resid_or_1,
        "resid  1  or   resnum 9   5 -7 9",
        Select::Or(
            Box::from(Select::ResidueNumber(vec![(1, 1)])),
            Box::from(Select::ResidueNumber(vec![(5, 7), (9, 9)]))
        )
    );

    parsing_success!(
        resid_or_2,
        "resnum 1   ||resid   5 - 7 9  ",
        Select::Or(
            Box::from(Select::ResidueNumber(vec![(1, 1)])),
            Box::from(Select::ResidueNumber(vec![(5, 7), (9, 9)]))
        )
    );

    parsing_success!(
        resid_and_1,
        "resid 1 and resid 5- 7 9",
        Select::And(
            Box::from(Select::ResidueNumber(vec![(1, 1)])),
            Box::from(Select::ResidueNumber(vec![(5, 7), (9, 9)]))
        )
    );

    parsing_success!(
        resid_and_2,
        "  resnum    1  && resnum   9  5-7   9",
        Select::And(
            Box::from(Select::ResidueNumber(vec![(1, 1)])),
            Box::from(Select::ResidueNumber(vec![(5, 7), (9, 9)]))
        )
    );

    parsing_success!(
        resid_complex,
        "resnum 1 and resnum 5-7 9 || resid 11 12to15",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueNumber(vec![(1, 1)])),
                Box::from(Select::ResidueNumber(vec![(5, 7), (9, 9)]))
            )),
            Box::from(Select::ResidueNumber(vec![(11, 15)]))
        )
    );

    parsing_success!(
        resid_complex_par_1,
        "(resid 1 &&resid 9  5-7 9) or resnum 11 12 to 15",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueNumber(vec![(1, 1)])),
                Box::from(Select::ResidueNumber(vec![(5, 7), (9, 9)]))
            )),
            Box::from(Select::ResidueNumber(vec![(11, 15)]))
        )
    );

    parsing_success!(
        resid_complex_par_2,
        "resid 1 &&(resnum 9   5to7 9 or resid 11 12 to15)",
        Select::And(
            Box::from(Select::ResidueNumber(vec![(1, 1)])),
            Box::from(Select::Or(
                Box::from(Select::ResidueNumber(vec![(5, 7), (9, 9)])),
                Box::from(Select::ResidueNumber(vec![(11, 15)]))
            )),
        )
    );

    // atom numbers
    parsing_success!(
        serial_or_1,
        "serial  1  or   atomid 9  5 -7 9",
        Select::Or(
            Box::from(Select::GmxAtomNumber(vec![(1, 1)])),
            Box::from(Select::AtomNumber(vec![(5, 7), (9, 9)]))
        )
    );

    parsing_success!(
        serial_or_2,
        "atomid 1   ||serial   5 - 7 9  ",
        Select::Or(
            Box::from(Select::AtomNumber(vec![(1, 1)])),
            Box::from(Select::GmxAtomNumber(vec![(5, 7), (9, 9)]))
        )
    );

    parsing_success!(
        serial_and_1,
        "serial 1 and serial 5- 7 9",
        Select::And(
            Box::from(Select::GmxAtomNumber(vec![(1, 1)])),
            Box::from(Select::GmxAtomNumber(vec![(5, 7), (9, 9)]))
        )
    );

    parsing_success!(
        serial_and_2,
        "  atomnum    1  && atomid   9  5-7   9",
        Select::And(
            Box::from(Select::AtomNumber(vec![(1, 1)])),
            Box::from(Select::AtomNumber(vec![(5, 7), (9, 9)]))
        )
    );

    parsing_success!(
        serial_complex,
        "serial 1 and atomnum 9  5-7 9 || serial 11 12-15",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::GmxAtomNumber(vec![(1, 1)])),
                Box::from(Select::AtomNumber(vec![(5, 7), (9, 9)]))
            )),
            Box::from(Select::GmxAtomNumber(vec![(11, 15)]))
        )
    );

    parsing_success!(
        serial_complex_par_1,
        "(serial 1 &&atomnum 9  5-7 9) or atomid 11 12 to 15",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::GmxAtomNumber(vec![(1, 1)])),
                Box::from(Select::AtomNumber(vec![(5, 7), (9, 9)]))
            )),
            Box::from(Select::AtomNumber(vec![(11, 15)]))
        )
    );

    parsing_success!(
        serial_complex_par_2,
        "serial 1 &&(serial 5to 7 9 or serial 11 12 - 15)",
        Select::And(
            Box::from(Select::GmxAtomNumber(vec![(1, 1)])),
            Box::from(Select::Or(
                Box::from(Select::GmxAtomNumber(vec![(5, 7), (9, 9)])),
                Box::from(Select::GmxAtomNumber(vec![(11, 15)]))
            )),
        )
    );

    // residue names with atom names
    parsing_success!(complex_resnames_names, "(resname 'POPE'  LYS LEU and(name   BB PO4   D2A) || name C1A ) ||(resname ION&& name CL) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![Name::new("POPE").unwrap(), Name::new("LYS").unwrap(), Name::new("LEU").unwrap()])),
                Box::from(Select::AtomName(vec![Name::new("BB").unwrap(), Name::new("PO4").unwrap(), Name::new("D2A").unwrap()]))
            )),
            Box::from(Select::AtomName(vec![Name::new("C1A").unwrap()]))
        )),
        Box::from(Select::And(
            Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
            Box::from(Select::AtomName(vec![Name::new("CL").unwrap()]))
        ))
    ));

    // residue names with labeled atom regexes
    parsing_success!(complex_resnames_label_regex, "(resname 'POPE'  LYS LEU and(label   MyAtom r'^A.*m$') || label Atom ) ||(resname ION&& label ' another atom') ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![Name::new("POPE").unwrap(), Name::new("LYS").unwrap(), Name::new("LEU").unwrap()])),
                Box::from(Select::LabeledAtom(vec![Name::new("MyAtom").unwrap(), Name::new("r'^A.*m$'").unwrap()]))
            )),
            Box::from(Select::LabeledAtom(vec![Name::new("Atom").unwrap()]))
        )),
        Box::from(Select::And(
            Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
            Box::from(Select::LabeledAtom(vec![Name::new(" another atom").unwrap()]))
        ))
    ));

    // residue names with residue numbers
    parsing_success!(complex_resnames_resid, "(resname 'POPE'  LYS LEU && (resid   15 22-25 33)or(resid 5 to 10) ) or(resname ION&& resnum 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![Name::new("POPE").unwrap(), Name::new("LYS").unwrap(), Name::new("LEU").unwrap()])),
                Box::from(Select::ResidueNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::ResidueNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
            Box::from(Select::ResidueNumber(vec![(6, 6)]))
        ))
    ));

    // residue names with atom numbers
    parsing_success!(complex_resnames_serial, "(resname 'POPE'  LYS LEU && (serial   33 22 -25 15) || atomid 5  to  10 ) ||(resname ION&& atomnum 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![Name::new("POPE").unwrap(), Name::new("LYS").unwrap(), Name::new("LEU").unwrap()])),
                Box::from(Select::GmxAtomNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::AtomNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
            Box::from(Select::AtomNumber(vec![(6, 6)]))
        ))
    ));

    // residue names with group names
    parsing_success!(complex_resnames_group, "(resname 'POPE'  LYS LEU && (Protein)or'Charged   Membrane' ) ||(resname ION&& group Membrane ION W) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![Name::new("POPE").unwrap(), Name::new("LYS").unwrap(), Name::new("LEU").unwrap()])),
                Box::from(Select::GroupName(vec![Name::new("Protein").unwrap()]))
            )),
            Box::from(Select::GroupName(vec![Name::new("Charged   Membrane").unwrap()]))
        )),
        Box::from(Select::And(
            Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
            Box::from(Select::GroupName(vec![Name::new("Membrane").unwrap(), Name::new("ION").unwrap(), Name::new("W").unwrap()]))
        ))
    ));

    // atom names with residue numbers
    parsing_success!(complex_names_resid, "(name 'BB'  PO4 D2A && (resid   33 22-25 15) || resid 5 to 10 ) or(atomname NA&& resnum 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::AtomName(vec![Name::new("BB").unwrap(), Name::new("PO4").unwrap(), Name::new("D2A").unwrap()])),
                Box::from(Select::ResidueNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::ResidueNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::AtomName(vec![Name::new("NA").unwrap()])),
            Box::from(Select::ResidueNumber(vec![(6, 6)]))
        ))
    ));

    // atom names with atom numbers
    parsing_success!(complex_names_serial, "(name 'BB'  PO4 D2A && (atomnum   15 22- 25 33) or serial 5 to 10 ) ||(atomname NA&& atomid 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::AtomName(vec![Name::new("BB").unwrap(), Name::new("PO4").unwrap(), Name::new("D2A").unwrap()])),
                Box::from(Select::AtomNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::GmxAtomNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::AtomName(vec![Name::new("NA").unwrap()])),
            Box::from(Select::AtomNumber(vec![(6, 6)]))
        ))
    ));

    // atom names with group names
    parsing_success!(complex_names_group, "(name 'BB'  PO4 D2A and(group Protein)|| 'Charged   Membrane'   W   )or(atomname NA&& Membrane) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::AtomName(vec![Name::new("BB").unwrap(), Name::new("PO4").unwrap(), Name::new("D2A").unwrap()])),
                Box::from(Select::GroupName(vec![Name::new("Protein").unwrap()]))
            )),
            Box::from(Select::GroupName(vec![Name::new("Charged   Membrane").unwrap(), Name::new("W").unwrap()]))
        )),
        Box::from(Select::And(
            Box::from(Select::AtomName(vec![Name::new("NA").unwrap()])),
            Box::from(Select::GroupName(vec![Name::new("Membrane").unwrap()]))
        ))
    ));

    // atom names with labeled atoms
    parsing_success!(complex_names_label, "(name 'BB'  PO4 D2A and(label MyAtom)|| label 'Atom with    long space'   This   )or(atomname NA&& label MyAtom) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::AtomName(vec![Name::new("BB").unwrap(), Name::new("PO4").unwrap(), Name::new("D2A").unwrap()])),
                Box::from(Select::LabeledAtom(vec![Name::new("MyAtom").unwrap()]))
            )),
            Box::from(Select::LabeledAtom(vec![Name::new("Atom with    long space").unwrap(), Name::new("This").unwrap()]))
        )),
        Box::from(Select::And(
            Box::from(Select::AtomName(vec![Name::new("NA").unwrap()])),
            Box::from(Select::LabeledAtom(vec![Name::new("MyAtom").unwrap()]))
        ))
    ));

    // residue numbers with atom numbers
    parsing_success!(complex_resid_serial, "(serial 4 to 8- 12 && (resid   15 22-25 33) || resid 5 to 10 ) ||(atomid 9 10 11 12&& resnum 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::GmxAtomNumber(vec![(4, 12)])),
                Box::from(Select::ResidueNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::ResidueNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::AtomNumber(vec![(9, 12)])),
            Box::from(Select::ResidueNumber(vec![(6, 6)]))
        ))
    ));

    // residue numbers with group names
    parsing_success!(complex_resid_group, "(Protein Membrane && (resid   15 22-25 33) || resid 5 to 10 ) ||('Water with Ions' && resnum 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::GroupName(vec![Name::new("Protein").unwrap(), Name::new("Membrane").unwrap()])),
                Box::from(Select::ResidueNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::ResidueNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::GroupName(vec![Name::new("Water with Ions").unwrap()])),
            Box::from(Select::ResidueNumber(vec![(6, 6)]))
        ))
    ));

    // residue numbers with labeled atoms
    parsing_success!(complex_resid_label, "(label MyAtom 'selected atom' && (resid   15 22-25 33) || resid 5 to 10 ) ||(label 'My favorite atom' && resnum 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::LabeledAtom(vec![Name::new("MyAtom").unwrap(), Name::new("selected atom").unwrap()])),
                Box::from(Select::ResidueNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::ResidueNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::LabeledAtom(vec![Name::new("My favorite atom").unwrap()])),
            Box::from(Select::ResidueNumber(vec![(6, 6)]))
        ))
    ));

    parsing_success!(complex_serial_group, "(Protein Membrane && (serial   15 22-25 33) || atomid 5 to 10 ) ||('Water with Ions' && atomnum 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::GroupName(vec![Name::new("Protein").unwrap(), Name::new("Membrane").unwrap()])),
                Box::from(Select::GmxAtomNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::AtomNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::GroupName(vec![Name::new("Water with Ions").unwrap()])),
            Box::from(Select::AtomNumber(vec![(6, 6)]))
        ))
    ));

    parsing_success!(complex_not_1, "! (name 'BB'  PO4 D2A && (atomnum   15 22- 25 33) or serial 5 to 10 ) ||(atomname NA&& atomid 6) ", 
    Select::Or(
        Box::from(Select::Not(
            Box::from(Select::Or(
                Box::from(Select::And(
                    Box::from(Select::AtomName(vec![Name::new("BB").unwrap(), Name::new("PO4").unwrap(), Name::new("D2A").unwrap()])),
                    Box::from(Select::AtomNumber(vec![(15, 15), (22, 25), (33, 33)]))
                )),
                Box::from(Select::GmxAtomNumber(vec![(5, 10)]))
            ))
        )),
        Box::from(Select::And(
            Box::from(Select::AtomName(vec![Name::new("NA").unwrap()])),
            Box::from(Select::AtomNumber(vec![(6, 6)]))
        ))
    ));

    parsing_success!(complex_not_2, "(resname 'POPE'  LYS LEU and(name   BB PO4   D2A) || name C1A ) ||not(resname ION&& name CL) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![Name::new("POPE").unwrap(), Name::new("LYS").unwrap(), Name::new("LEU").unwrap()])),
                Box::from(Select::AtomName(vec![Name::new("BB").unwrap(), Name::new("PO4").unwrap(), Name::new("D2A").unwrap()]))
            )),
            Box::from(Select::AtomName(vec![Name::new("C1A").unwrap()]))
        )),
        Box::from(Select::Not(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
                Box::from(Select::AtomName(vec![Name::new("CL").unwrap()]))
            ))
        ))
    ));

    parsing_success!(complex_not_3, "  not(!(resname 'POPE' LYS LEU and(!name   BB PO4   D2A) || not name C1A ) ||(resname ION&& name CL) )", 
    Select::Not(
        Box::from(Select::Or(
            Box::from(Select::Not(
                Box::from(Select::Or(
                    Box::from(Select::And(
                        Box::from(Select::ResidueName(vec![Name::new("POPE").unwrap(), Name::new("LYS").unwrap(), Name::new("LEU").unwrap()])),
                        Box::from(Select::Not(
                            Box::from(Select::AtomName(vec![Name::new("BB").unwrap(), Name::new("PO4").unwrap(), Name::new("D2A").unwrap()]))
                        ))
                    )),
                    Box::from(Select::Not(Box::from(Select::AtomName(vec![Name::new("C1A").unwrap()]))))
                ))
            )),
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
                Box::from(Select::AtomName(vec![Name::new("CL").unwrap()]))
            ))
        ))
    ));

    parsing_success!(
        complex_not_3_multiline,
        "  not(!(resname 'POPE' LYS LEU and
    (!name   BB 
        PO4   D2A) || not name C1A ) ||(resname ION
            && name CL) )",
        Select::Not(Box::from(Select::Or(
            Box::from(Select::Not(Box::from(Select::Or(
                Box::from(Select::And(
                    Box::from(Select::ResidueName(vec![
                        Name::new("POPE").unwrap(),
                        Name::new("LYS").unwrap(),
                        Name::new("LEU").unwrap(),
                    ])),
                    Box::from(Select::Not(Box::from(Select::AtomName(vec![
                        Name::new("BB").unwrap(),
                        Name::new("PO4").unwrap(),
                        Name::new("D2A").unwrap(),
                    ]))))
                )),
                Box::from(Select::Not(Box::from(Select::AtomName(vec![Name::new(
                    "C1A"
                )
                .unwrap()]))))
            )))),
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
                Box::from(Select::AtomName(vec![Name::new("CL").unwrap()]))
            ))
        )))
    );

    parsing_success_string!(complex_parentheses_1,
        "!(!(name BB and resid 15 to 18) || ((resname   POPE POPG &&name PO4  )or not(name C1A||(serial 5 to 12 or group 'Protein 2' Membrane and !resid 1 2 3) )))",
        "Not(Or(Not(And(AtomName([String(\"BB\")]), ResidueNumber([(15, 18)]))), Or(And(ResidueName([String(\"POPE\"), String(\"POPG\")]), AtomName([String(\"PO4\")])), Not(Or(AtomName([String(\"C1A\")]), And(Or(GmxAtomNumber([(5, 12)]), GroupName([String(\"Protein 2\"), String(\"Membrane\")])), Not(ResidueNumber([(1, 3)]))))))))"
    );

    parsing_success!(
        keywords_in_quotes,
        "group 'Protein and  Membrane' or Membrane ' W or not ION' \"Lipids   to Count\"",
        Select::Or(
            Box::from(Select::GroupName(vec![
                Name::new("Protein and  Membrane").unwrap()
            ])),
            Box::from(Select::GroupName(vec![
                Name::new("Membrane").unwrap(),
                Name::new(" W or not ION").unwrap(),
                Name::new("Lipids   to Count").unwrap()
            ]))
        )
    );

    parsing_success!(
        coalesced_keywords_1,
        "resname BB andserial 1 2 3",
        Select::ResidueName(vec![
            Name::new("BB").unwrap(),
            Name::new("andserial").unwrap(),
            Name::new("1").unwrap(),
            Name::new("2").unwrap(),
            Name::new("3").unwrap(),
        ])
    );

    parsing_success!(
        coalesced_keywords_2,
        "resname BB or notserial 1 2 3",
        Select::Or(
            Box::from(Select::ResidueName(vec![Name::new("BB").unwrap()])),
            Box::from(Select::GroupName(vec![
                Name::new("notserial").unwrap(),
                Name::new("1").unwrap(),
                Name::new("2").unwrap(),
                Name::new("3").unwrap(),
            ]))
        )
    );

    parsing_success!(
        macro_dna,
        "@dna",
        Select::ResidueName(
            vec![
                "DA", "DG", "DC", "DT", "DA5", "DG5", "DC5", "DT5", "DA3", "DG3", "DC3", "DT3",
                "DAN", "DGN", "DCN", "DTN"
            ]
            .iter()
            .map(|&s| Name::new(s).unwrap())
            .collect()
        )
    );

    parsing_success!(
        regex_1,
        "resname r'^.*L.*' and name CA r'^[0-9]*H.*'",
        Select::And(
            Box::from(Select::ResidueName(vec![Name::new("r'^.*L.*'").unwrap()])),
            Box::from(Select::AtomName(vec![
                Name::new("CA").unwrap(),
                Name::new("r'^[0-9]*H.*'").unwrap()
            ]))
        )
    );

    parsing_success!(
        regex_2,
        "name r'^C.*' r'^[0-9]*H.*' r'.*'",
        Select::AtomName(vec![
            Name::new("r'^C.*'").unwrap(),
            Name::new("r'^[0-9]*H.*'").unwrap(),
            Name::new("r'.*'").unwrap()
        ])
    );

    parsing_success!(
        nonascii_group_1,
        "'velmi česká řřř skupina, že jo?'",
        Select::GroupName(vec![Name::new("velmi česká řřř skupina, že jo?").unwrap()])
    );

    parsing_success!(
        nonascii_group_2,
        "'這是一個中文名字，希望是對的'",
        Select::GroupName(vec![Name::new("這是一個中文名字，希望是對的").unwrap()])
    );

    parsing_success!(
        nonascii_atomname,
        "name underhålls uppföljnings",
        Select::AtomName(vec![
            Name::new("underhålls").unwrap(),
            Name::new("uppföljnings").unwrap()
        ])
    );

    parsing_success!(
        operator_in_regex,
        "name r'C3[2-9]|C3[1][0-6]|C2[2-9]|C2[1][0-8]' P",
        Select::AtomName(vec![
            Name::new("r'C3[2-9]|C3[1][0-6]|C2[2-9]|C2[1][0-8]'").unwrap(),
            Name::new("P").unwrap()
        ])
    );

    parsing_success!(
        operator_in_regex_2,
        "name r'C3[2-9]|C3[1][0-6]|C2[2-9]|C2[1][0-8]' and resname POPC",
        Select::And(
            Box::from(Select::AtomName(vec![Name::new(
                "r'C3[2-9]|C3[1][0-6]|C2[2-9]|C2[1][0-8]'"
            )
            .unwrap()])),
            Box::from(Select::ResidueName(vec![Name::new("POPC").unwrap()]))
        )
    );

    parsing_success!(
        operator_in_regex_3,
        "name r'C3[2-9]&C3[1][0-6]&C2[2-9]&C2[1][0-8]' and resname POPC",
        Select::And(
            Box::from(Select::AtomName(vec![Name::new(
                "r'C3[2-9]&C3[1][0-6]&C2[2-9]&C2[1][0-8]'"
            )
            .unwrap()])),
            Box::from(Select::ResidueName(vec![Name::new("POPC").unwrap()]))
        )
    );

    /*
    // deprecated since v0.6.0
    parsing_success!(
        hydrogen_macro,
        "@hydrogen",
        Select::AtomName(vec![Name::new("r'^[1-9]?[Hh].*'").unwrap()])
    );*/
}

#[cfg(test)]
mod fail_tests {
    use super::*;

    macro_rules! parsing_fails {
        ($name:ident, $expression:expr, $variant:path) => {
            #[test]
            fn $name() {
                let query = $expression;

                match Select::parse_query(query) {
                    Err($variant(e)) => assert_eq!(e, query),
                    Ok(_) => panic!("Parsing should have failed, but it succeeded."),
                    Err(e) => panic!("Parsing successfully failed but incorrect error type `{:?}` was returned.", e),
                }
            }
        };
    }

    #[test]
    fn empty_query() {
        let query = "";

        match Select::parse_query(query) {
            Err(SelectError::EmptyQuery) => (),
            Ok(_) => panic!("Parsing should have failed, but it succeeded."),
            Err(e) => panic!(
                "Parsing successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    parsing_fails!(
        empty_resname,
        "resname LYS and resname",
        SelectError::EmptyArgument
    );
    parsing_fails!(empty_name, "name BB or name", SelectError::EmptyArgument);
    parsing_fails!(
        empty_resid,
        "resid 1-3 or resid",
        SelectError::EmptyArgument
    );
    parsing_fails!(
        empty_serial,
        "atomnum 65 66 69 and serial",
        SelectError::EmptyArgument
    );
    parsing_fails!(
        empty_group,
        "Protein Membrane and group",
        SelectError::EmptyArgument
    );

    parsing_fails!(invalid_number_1, "resid 1 to x", SelectError::InvalidNumber);
    parsing_fails!(
        invalid_number_2,
        "serial 25 24 23 22 21 2O 19 18 17",
        SelectError::InvalidNumber
    );
    parsing_fails!(
        invalid_number_3,
        "serial 19 20 21 21.5 22 22.3",
        SelectError::InvalidNumber
    );

    parsing_fails!(invalid_range_1, "resid 25-20", SelectError::InvalidNumber);
    parsing_fails!(invalid_range_2, "resid 25 -20", SelectError::InvalidNumber);
    parsing_fails!(invalid_range_3, "resid 25- 20", SelectError::InvalidNumber);
    parsing_fails!(invalid_range_4, "resid 25 - 20", SelectError::InvalidNumber);

    parsing_fails!(
        invalid_open_ended_range_1,
        "serial <",
        SelectError::InvalidNumber
    );
    parsing_fails!(
        invalid_open_ended_range_2,
        "serial <==7",
        SelectError::InvalidNumber
    );
    parsing_fails!(
        invalid_open_ended_range_3,
        "resid <<=  7",
        SelectError::InvalidNumber
    );
    parsing_fails!(
        invalid_open_ended_range_4,
        "atomnum <<=7",
        SelectError::InvalidNumber
    );
    parsing_fails!(
        invalid_open_ended_range_5,
        "atomnum 1 5 >",
        SelectError::InvalidNumber
    );
    parsing_fails!(
        invalid_open_ended_range_6,
        "atomnum 1 -> 8",
        SelectError::InvalidNumber
    );

    parsing_fails!(chain_multichar_1, "chain AB", SelectError::InvalidChainId);
    parsing_fails!(
        chain_multichar_2,
        "chain myidentifier",
        SelectError::InvalidChainId
    );
    parsing_fails!(
        chain_multichar_3,
        "chain A B C chain",
        SelectError::InvalidChainId
    );

    parsing_fails!(
        fail_parentheses_1,
        "(resname LYS and name SC1",
        SelectError::InvalidParentheses
    );
    parsing_fails!(
        fail_parentheses_2,
        "resname LYS and name SC1)",
        SelectError::InvalidParentheses
    );
    parsing_fails!(
        fail_parentheses_3,
        "((resname LYS and name SC1)",
        SelectError::InvalidParentheses
    );
    parsing_fails!(
        fail_parentheses_4,
        "(((resname LYS and name SC1))))",
        SelectError::InvalidParentheses
    );
    parsing_fails!(
        fail_parentheses_5,
        "(resname LYS) and (name SC1))",
        SelectError::InvalidParentheses
    );

    parsing_fails!(
        missing_argument_1,
        "resname LYS and",
        SelectError::MissingArgument
    );
    parsing_fails!(
        missing_argument_2,
        "or serial 2-154",
        SelectError::MissingArgument
    );
    parsing_fails!(
        missing_argument_3,
        "or resid 15 and",
        SelectError::MissingArgument
    );
    parsing_fails!(
        missing_argument_4,
        "(name BB) and ",
        SelectError::MissingArgument
    );
    parsing_fails!(
        missing_argument_5,
        "resname POPG && (&& serial 1 to 44 || resname POPE)",
        SelectError::MissingArgument
    );

    parsing_fails!(
        messed_up_1,
        "(resname LYS) (and (name SC1))",
        SelectError::MissingArgument
    );

    parsing_fails!(
        fail_quotes_1,
        "group 'Group number 1' && resname LYS and (Protein 'Membrane With Ions)",
        SelectError::InvalidQuotes
    );
    parsing_fails!(
        fail_quotes_2,
        "group \"Group number 1\" resname LYS and (Protein \"Membrane With Ions)",
        SelectError::InvalidQuotes
    );

    parsing_fails!(
        invalid_operator_1,
        "resname LYS & name SC1",
        SelectError::InvalidOperator
    );
    parsing_fails!(
        invalid_operator_2,
        "resname LYS | name SC1",
        SelectError::InvalidOperator
    );
    parsing_fails!(
        invalid_operator_3,
        "resname LYS &&& name SC1",
        SelectError::InvalidOperator
    );
    parsing_fails!(
        invalid_operator_4,
        "resname LYS &@ name SC1",
        SelectError::InvalidOperator
    );
    parsing_fails!(
        invalid_operator_5,
        "resname LYS @& name SC1",
        SelectError::InvalidOperator
    );
    parsing_fails!(
        invalid_token_after_parentheses_1,
        "(name CA CB) resname LYS",
        SelectError::InvalidTokenParentheses
    );
    parsing_fails!(
        invalid_token_after_parentheses_2,
        "(name BB or group Protein) (resid 1 to 54)",
        SelectError::InvalidTokenParentheses
    );
    parsing_fails!(
        invalid_token_after_parentheses_3,
        "(name BB or group Protein)!serial 7",
        SelectError::InvalidTokenParentheses
    );
    parsing_fails!(
        invalid_token_before_parentheses_1,
        "element (name CA)",
        SelectError::InvalidTokenParentheses
    );
    parsing_fails!(
        invalid_token_before_parentheses_2,
        "Membrane (serial 1 to 3)",
        SelectError::InvalidTokenParentheses
    );

    #[test]
    fn invalid_regex() {
        let query = "name r'*L*'";

        match Select::parse_query(query) {
            Err(SelectError::InvalidRegex(e)) => assert_eq!(e, "r'*L*'"),
            Ok(_) => panic!("Parsing should have failed, but it succeeded."),
            Err(e) => panic!(
                "Parsing successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    parsing_fails!(deep_error_1, "!(!(name BB and resid 15 to 18) || ((resname   POPE POPG &&name PO4  )or not(name C1A||(serial x to 12 or group 'Protein 2' Membrane and !resid 1 2 3) )))", 
    SelectError::InvalidNumber);
}

#[cfg(test)]
mod select_impl {
    use super::*;

    #[test]
    fn expand_regex_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();

        let selection =
            Select::parse_query("(Protein r'membrane$' or r'^P' ION r'-') and !r'C'").unwrap();

        let selection = selection.expand_regex_group_label(&system).unwrap();

        let string = format!("{:?}", selection);
        assert_eq!(string, "And(Or(GroupName([String(\"Protein\"), String(\"Transmembrane\")]), GroupName([String(\"Protein\"), String(\"Protein-H\"), String(\"Prot-Masses\"), String(\"POPC\"), String(\"Protein_Membrane\"), String(\"ION\"), String(\"Protein-H\"), String(\"C-alpha\"), String(\"SideChain-H\"), String(\"Prot-Masses\"), String(\"non-Protein\")])), Not(GroupName([String(\"C-alpha\"), String(\"MainChain\"), String(\"MainChain+Cb\"), String(\"MainChain+H\"), String(\"SideChain\"), String(\"SideChain-H\"), String(\"POPC\")])))");
    }

    #[test]
    fn expand_regex_label() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.label_atom("atom", 174).unwrap();
        system.label_atom("atom2", 184).unwrap();
        system.label_atom("atom_new", 1923).unwrap();
        system.label_atom("different", 438).unwrap();

        let selection = Select::parse_query("resname ION and label different r'a'").unwrap();
        let selection = selection.expand_regex_group_label(&system).unwrap();

        let string = format!("{:?}", selection);
        assert_eq!(string, "And(ResidueName([String(\"ION\")]), LabeledAtom([String(\"different\"), String(\"atom_new\"), String(\"atom\"), String(\"atom2\")]))");
    }

    #[test]
    fn expand_regex_group_label() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();
        system.label_atom("atom", 174).unwrap();
        system.label_atom("atom2", 184).unwrap();
        system.label_atom("atom_new", 1923).unwrap();
        system.label_atom("different", 438).unwrap();

        let selection = Select::parse_query(
            "((Protein r'membrane$' or r'^P' ION r'-') and !r'C') or label r'a.*_new' different",
        )
        .unwrap();

        let selection = selection.expand_regex_group_label(&system).unwrap();
        let string = format!("{:?}", selection);
        assert_eq!(string, "Or(And(Or(GroupName([String(\"Protein\"), String(\"Transmembrane\")]), GroupName([String(\"Protein\"), String(\"Protein-H\"), String(\"Prot-Masses\"), String(\"POPC\"), String(\"Protein_Membrane\"), String(\"ION\"), String(\"Protein-H\"), String(\"C-alpha\"), String(\"SideChain-H\"), String(\"Prot-Masses\"), String(\"non-Protein\")])), Not(GroupName([String(\"C-alpha\"), String(\"MainChain\"), String(\"MainChain+Cb\"), String(\"MainChain+H\"), String(\"SideChain\"), String(\"SideChain-H\"), String(\"POPC\")]))), LabeledAtom([String(\"atom_new\"), String(\"different\")]))");
    }

    #[test]
    fn expand_regex_group_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();

        let selection =
            Select::parse_query("(Protein r'membrane$' or r'^P' Nonexistent r'-') and !r'C'")
                .unwrap();

        match selection.expand_regex_group_label(&system) {
            Ok(_) => panic!("Expansion should have failed."),
            Err(SelectError::GroupNotFound(e)) => assert_eq!(e, "Nonexistent"),
            Err(e) => panic!("Incorrect error '{}' returned.", e),
        }
    }

    #[test]
    fn expand_regex_label_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.label_atom("atom", 174).unwrap();
        system.label_atom("atom2", 184).unwrap();
        system.label_atom("atom_new", 1923).unwrap();

        let selection = Select::parse_query("resname ION and label different r'a'").unwrap();
        match selection.expand_regex_group_label(&system) {
            Ok(_) => panic!("Expansion should have failed."),
            Err(SelectError::LabelNotFound(e)) => assert_eq!(e, "different"),
            Err(e) => panic!("Incorrect error '{}' returned.", e),
        }
    }

    #[test]
    fn expand_regex_group_nomatch() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();

        let selection = Select::parse_query("r'^x'").unwrap();

        match selection.expand_regex_group_label(&system) {
            Ok(_) => panic!("Expansion should have failed."),
            Err(SelectError::NoRegexMatch(e)) => assert_eq!(e, "^x"),
            Err(e) => panic!("Incorrect error '{}' returned.", e),
        }
    }

    #[test]
    fn expand_regex_label_nomatch() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.label_atom("atom", 174).unwrap();
        system.label_atom("atom2", 184).unwrap();
        system.label_atom("atom_new", 1923).unwrap();

        let selection = Select::parse_query("resname ION and label r'atx'").unwrap();
        match selection.expand_regex_group_label(&system) {
            Ok(_) => panic!("Expansion should have failed."),
            Err(SelectError::NoRegexMatch(e)) => assert_eq!(e, "atx"),
            Err(e) => panic!("Incorrect error '{}' returned.", e),
        }
    }

    #[test]
    fn expand_regex_group_match() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();

        let selection = Select::parse_query("r'^x' r'^P'").unwrap();

        let selection = selection.expand_regex_group_label(&system).unwrap();

        let string = format!("{:?}", selection);
        assert_eq!(string, "GroupName([String(\"Protein\"), String(\"Protein-H\"), String(\"Prot-Masses\"), String(\"POPC\"), String(\"Protein_Membrane\")])")
    }

    macro_rules! convert_to_string {
        ($name:ident, $expression:expr, $expected:expr) => {
            #[test]
            fn $name() {
                let select = $expression;

                assert_eq!(select.to_string(), $expected);
            }
        };
    }

    convert_to_string!(
        empty_selections_1,
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::Chain(vec![])),
                Box::from(Select::GmxAtomNumber(vec![]))
            )),
            Box::from(Select::GroupName(vec![]))
        ),
        "((!(all) and !(all)) or !(all))"
    );

    convert_to_string!(
        empty_selections_2,
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![])),
                Box::from(Select::ElementSymbol(vec![]))
            )),
            Box::from(Select::AtomNumber(vec![]))
        ),
        "((!(all) and !(all)) or !(all))"
    );

    convert_to_string!(
        resname_to_string,
        Select::ResidueName(vec![Name::new("LYS").unwrap()]),
        "resname LYS"
    );

    convert_to_string!(
        resnames_to_string,
        Select::ResidueName(vec![
            Name::new("LYS").unwrap(),
            Name::new("CYS").unwrap(),
            Name::new("LEU").unwrap(),
            Name::new("POPE").unwrap()
        ]),
        "resname LYS CYS LEU POPE"
    );

    convert_to_string!(
        atomnames_to_string,
        Select::AtomName(vec![
            Name::new("BB").unwrap(),
            Name::new("SC1").unwrap(),
            Name::new("PO4").unwrap(),
            Name::new("C1").unwrap(),
            Name::new("C2B").unwrap()
        ]),
        "name BB SC1 PO4 C1 C2B"
    );

    convert_to_string!(
        resids_to_string,
        Select::ResidueNumber(vec![(4, 4), (8, 14), (20, 21)]),
        "resid 4 8 to 14 20 to 21"
    );

    convert_to_string!(
        serial_to_string,
        Select::GmxAtomNumber(vec![(15, 15)]),
        "serial 15"
    );

    convert_to_string!(
        serials_to_string,
        Select::GmxAtomNumber(vec![(1, 17), (24, 37), (89, usize::MAX)]),
        "serial <= 17 24 to 37 >= 89"
    );

    convert_to_string!(
        atomids_to_string,
        Select::AtomNumber(vec![(1, 1), (4, 4), (7, 11), (50, usize::MAX)]),
        "atomid 1 4 7 to 11 >= 50"
    );

    convert_to_string!(
        groups_to_string,
        Select::GroupName(vec![
            Name::new("Protein").unwrap(),
            Name::new("Membrane Only POPC  ").unwrap(),
            Name::new(" ION with  Water").unwrap()
        ]),
        "group Protein 'Membrane Only POPC  ' ' ION with  Water'"
    );

    convert_to_string!(
        chains_to_string,
        Select::Chain(vec!['B', 'C', 'E', 'F']),
        "chain B C E F"
    );

    convert_to_string!(
        elsymbols_to_string,
        Select::ElementSymbol(vec![
            Name::new("C").unwrap(),
            Name::new("Na").unwrap(),
            Name::new("K").unwrap(),
            Name::new("Br").unwrap()
        ]),
        "element symbol C Na K Br"
    );

    convert_to_string!(
        elnames_to_string,
        Select::ElementName(vec![
            Name::new("carbon number 2  ").unwrap(),
            Name::new("hydrogen").unwrap()
        ]),
        "element name 'carbon number 2  ' hydrogen"
    );

    convert_to_string!(
        regex_to_string,
        Select::And(
            Box::from(Select::ResidueName(vec![Name::new("r'^.*L.*'").unwrap()])),
            Box::from(Select::AtomName(vec![
                Name::new("CA").unwrap(),
                Name::new("r'^[0-9]*H.*'").unwrap()
            ]))
        ),
        "(resname r'^.*L.*' and name CA r'^[0-9]*H.*')"
    );

    convert_to_string!(
        complex_to_string_1,
        Select::Not(
            Box::from(Select::Or(
                Box::from(Select::Not(
                    Box::from(Select::Or(
                        Box::from(Select::And(
                            Box::from(Select::ResidueName(vec![Name::new("POPE").unwrap(), Name::new("LYS").unwrap(), Name::new("LEU").unwrap()])),
                            Box::from(Select::Not(
                                Box::from(Select::AtomName(vec![Name::new("BB").unwrap(), Name::new("PO4").unwrap(), Name::new("D2A").unwrap()]))
                            ))
                        )),
                        Box::from(Select::Not(Box::from(Select::AtomName(vec![Name::new("C1A").unwrap()]))))
                    ))
                )),
                Box::from(Select::And(
                    Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
                    Box::from(Select::AtomName(vec![Name::new("CL").unwrap()]))
                ))
            ))
        ),
        "not (not ((resname POPE LYS LEU and not name BB PO4 D2A) or not name C1A) or (resname ION and name CL))"
    );

    convert_to_string!(
        complex_to_string_2,
        Select::Or(
            Box::from(Select::Not(
                Box::from(Select::Or(
                    Box::from(Select::And(
                        Box::from(Select::AtomName(vec![Name::new("BB").unwrap(), Name::new("PO4").unwrap(), Name::new("D2A").unwrap()])),
                        Box::from(Select::AtomNumber(vec![(15, 15), (22, 25), (33, 33)]))
                    )),
                    Box::from(Select::GmxAtomNumber(vec![(5, 10)]))
                ))
            )),
            Box::from(Select::And(
                Box::from(Select::AtomName(vec![Name::new("NA").unwrap()])),
                Box::from(Select::AtomNumber(vec![(6, 6)]))
            ))
        ),
        "(not ((name BB PO4 D2A and atomid 15 22 to 25 33) or serial 5 to 10) or (name NA and atomid 6))"
    );

    convert_to_string!(
        molecule_with_to_string_1,
        Select::Molecule(Box::new(Select::GmxAtomNumber(vec![(1, 12), (17, 17)]))),
        "molecule with serial <= 12 17"
    );

    convert_to_string!(
        molecule_with_to_string_2,
        Select::Molecule(Box::new(Select::And(
            Box::new(Select::Molecule(Box::new(Select::AtomName(vec![
                Name::new("P").unwrap()
            ])))),
            Box::new(Select::ResidueNumber(vec![(50, 76)])),
        ))),
        "molecule with (molecule with name P and resid 50 to 76)"
    );

    macro_rules! convert_forward_back {
        ($name:ident, $expression:expr) => {
            #[test]
            fn $name() {
                let select = $expression;

                assert_eq!(
                    Select::parse_query(&select.to_string()).unwrap(),
                    Box::from($expression)
                );
            }
        };
    }

    convert_forward_back!(
        forward_back_1,
        Select::Not(Box::from(Select::Or(
            Box::from(Select::Not(Box::from(Select::Or(
                Box::from(Select::And(
                    Box::from(Select::ResidueName(vec![
                        Name::new("POPE").unwrap(),
                        Name::new("LYS").unwrap(),
                        Name::new("LEU").unwrap()
                    ])),
                    Box::from(Select::Not(Box::from(Select::AtomName(vec![
                        Name::new("BB").unwrap(),
                        Name::new("PO4").unwrap(),
                        Name::new("D2A").unwrap()
                    ]))))
                )),
                Box::from(Select::Not(Box::from(Select::AtomName(vec![Name::new(
                    "C1A"
                )
                .unwrap()]))))
            )))),
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
                Box::from(Select::AtomName(vec![Name::new("CL").unwrap()]))
            ))
        )))
    );

    convert_forward_back!(
        forward_back_2,
        Select::Or(
            Box::from(Select::Not(Box::from(Select::Or(
                Box::from(Select::And(
                    Box::from(Select::AtomName(vec![
                        Name::new("BB").unwrap(),
                        Name::new("PO4").unwrap(),
                        Name::new("D2A").unwrap()
                    ])),
                    Box::from(Select::AtomNumber(vec![(15, 15), (22, 25), (33, 33)]))
                )),
                Box::from(Select::GmxAtomNumber(vec![(5, 10)]))
            )))),
            Box::from(Select::And(
                Box::from(Select::AtomName(vec![Name::new("NA").unwrap()])),
                Box::from(Select::AtomNumber(vec![(6, 6)]))
            ))
        )
    );

    convert_forward_back!(
        forward_back_3,
        Select::Molecule(Box::new(Select::And(
            Box::new(Select::Molecule(Box::new(Select::AtomName(vec![
                Name::new("P").unwrap()
            ])))),
            Box::new(Select::ResidueNumber(vec![(50, 76)])),
        )))
    );

    #[test]
    fn convert_forward_back_fail() {
        let select = Select::Or(
            Box::from(Select::And(
                Box::from(Select::Chain(vec![])),
                Box::from(Select::GmxAtomNumber(vec![])),
            )),
            Box::from(Select::GroupName(vec![])),
        );

        assert_ne!(
            Select::parse_query(&select.to_string()).unwrap(),
            Box::from(select)
        );
    }
}

#[cfg(test)]
#[cfg(feature = "serde")]
mod serde_tests {
    use super::*;

    #[test]
    fn simple_selection_to_yaml() {
        let select = Select::ResidueNumber(vec![(4, 11)]);
        let string = serde_yaml::to_string(&select).unwrap();

        assert_eq!(string, "resid 4 to 11\n");
    }

    #[test]
    fn simple_query_from_yaml() {
        let string = "resid 4 to 11";
        let select: Select = serde_yaml::from_str(&string).unwrap();

        assert_eq!(select, Select::ResidueNumber(vec![(4, 11)]));
    }

    #[test]
    fn complex_selection_to_yaml() {
        let select = Select::Not(Box::from(Select::Or(
            Box::from(Select::Not(Box::from(Select::Or(
                Box::from(Select::And(
                    Box::from(Select::ResidueName(vec![
                        Name::new("POPE").unwrap(),
                        Name::new("LYS").unwrap(),
                        Name::new("LEU").unwrap(),
                    ])),
                    Box::from(Select::Not(Box::from(Select::AtomName(vec![
                        Name::new("BB").unwrap(),
                        Name::new("PO4").unwrap(),
                        Name::new("D2A").unwrap(),
                    ])))),
                )),
                Box::from(Select::Not(Box::from(Select::AtomName(vec![Name::new(
                    "C1A",
                )
                .unwrap()])))),
            )))),
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
                Box::from(Select::AtomName(vec![Name::new("CL").unwrap()])),
            )),
        )));

        let string = serde_yaml::to_string(&select).unwrap();
        assert_eq!(string, "not (not ((resname POPE LYS LEU and not name BB PO4 D2A) or not name C1A) or (resname ION and name CL))\n");
    }

    #[test]
    fn complex_query_from_yaml() {
        let string = "  not(!(resname 'POPE' LYS LEU and(!name   BB PO4   D2A) || not name C1A ) ||(resname ION&& name CL) )";
        let select: Select = serde_yaml::from_str(&string).unwrap();

        assert_eq!(
            select,
            Select::Not(Box::from(Select::Or(
                Box::from(Select::Not(Box::from(Select::Or(
                    Box::from(Select::And(
                        Box::from(Select::ResidueName(vec![
                            Name::new("POPE").unwrap(),
                            Name::new("LYS").unwrap(),
                            Name::new("LEU").unwrap()
                        ])),
                        Box::from(Select::Not(Box::from(Select::AtomName(vec![
                            Name::new("BB").unwrap(),
                            Name::new("PO4").unwrap(),
                            Name::new("D2A").unwrap()
                        ]))))
                    )),
                    Box::from(Select::Not(Box::from(Select::AtomName(vec![Name::new(
                        "C1A"
                    )
                    .unwrap()]))))
                )))),
                Box::from(Select::And(
                    Box::from(Select::ResidueName(vec![Name::new("ION").unwrap()])),
                    Box::from(Select::AtomName(vec![Name::new("CL").unwrap()]))
                ))
            )))
        );
    }

    #[test]
    fn vector_multiple_selects_to_yaml() {
        let selection1 = Select::ElementName(vec![
            Name::new("r'^ni'").unwrap(),
            Name::new("carbon").unwrap(),
            Name::new("potassium").unwrap(),
        ]);

        let selection2 = Select::Or(
            Box::from(Select::And(
                Box::from(Select::GmxAtomNumber(vec![(1, 3)])),
                Box::from(Select::Chain(vec!['A'])),
            )),
            Box::from(Select::Chain(vec!['C', 'D'])),
        );

        let selection3 = Select::AtomName(vec![
            Name::new("BB").unwrap(),
            Name::new("SC1").unwrap(),
            Name::new("PO4").unwrap(),
            Name::new("C1").unwrap(),
            Name::new("C2B").unwrap(),
        ]);

        let selections = vec![selection1, selection2, selection3];
        let string = serde_yaml::to_string(&selections).unwrap();
        let expected = "- element name r'^ni' carbon potassium
- ((serial <= 3 and chain A) or chain C D)
- name BB SC1 PO4 C1 C2B
";

        assert_eq!(string, expected);
    }

    #[test]
    fn vector_multiple_queries_from_yaml() {
        let string = "- elname r'^ni' carbon potassium
- (serial 1 to 3 and chain A) or chain C D
-   name   BB SC1   PO4   C1 C2B
";

        let selections: Vec<Select> = serde_yaml::from_str(&string).unwrap();

        assert_eq!(
            selections[0],
            Select::ElementName(vec![
                Name::new("r'^ni'").unwrap(),
                Name::new("carbon").unwrap(),
                Name::new("potassium").unwrap()
            ])
        );

        assert_eq!(
            selections[1],
            Select::Or(
                Box::from(Select::And(
                    Box::from(Select::GmxAtomNumber(vec![(1, 3)])),
                    Box::from(Select::Chain(vec!['A']))
                )),
                Box::from(Select::Chain(vec!['C', 'D']))
            )
        );

        assert_eq!(
            selections[2],
            Select::AtomName(vec![
                Name::new("BB").unwrap(),
                Name::new("SC1").unwrap(),
                Name::new("PO4").unwrap(),
                Name::new("C1").unwrap(),
                Name::new("C2B").unwrap()
            ])
        );
    }

    #[test]
    fn parsing_fail_yaml() {
        let string = "((resname LYS and name SC1)";

        let serde_error = match serde_yaml::from_str::<Select>(&string) {
            Ok(_) => panic!("Parsing should have failed. (YAML)"),
            Err(e) => e,
        };

        let select_error = match Select::parse_query(&string) {
            Ok(_) => panic!("Parsing should have failed. (QUERY)"),
            Err(e) => e,
        };

        let serde_formatted = format!("{}", serde_error);
        let select_formatted = format!("{}", select_error);

        assert_eq!(serde_formatted, select_formatted);
    }
}
