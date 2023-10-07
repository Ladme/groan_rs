// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the Groan selection language.

use std::collections::HashMap;

use crate::errors::SelectError;
use crate::selections::numbers;
use crate::structures::group::Group;

#[derive(Debug, PartialEq)]
pub enum Select {
    ResidueName(Vec<String>),
    AtomName(Vec<String>),
    ResidueNumber(Vec<(usize, usize)>),
    GmxAtomNumber(Vec<(usize, usize)>),
    AtomNumber(Vec<(usize, usize)>),
    Chain(Vec<char>),
    GroupName(Vec<String>),
    And(Box<Select>, Box<Select>),
    Or(Box<Select>, Box<Select>),
    Not(Box<Select>),
}

#[derive(Debug, PartialEq)]
enum Operator {
    And,
    Or,
    Not,
}

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

    // replace word operators with their symbolic equivalents
    expression = replace_keywords(&expression);

    match parse_subquery(&expression, 0, expression.len()) {
        Ok(x) => Ok(x),
        Err(SelectError::InvalidOperator(_)) => {
            Err(SelectError::InvalidOperator(query.to_string()))
        }
        Err(SelectError::MissingArgument(_)) => {
            Err(SelectError::MissingArgument(query.to_string()))
        }
        Err(SelectError::EmptyArgument(_)) => Err(SelectError::EmptyArgument(query.to_string())),
        Err(SelectError::InvalidParentheses(_)) => {
            Err(SelectError::InvalidParentheses(query.to_string()))
        }
        Err(SelectError::InvalidNumber(_)) => Err(SelectError::InvalidNumber(query.to_string())),
        Err(SelectError::InvalidChainId(_)) => Err(SelectError::InvalidChainId(query.to_string())),
        Err(_) => Err(SelectError::UnknownError(query.to_string())),
    }
}

fn parse_subquery(expression: &str, start: usize, end: usize) -> Result<Box<Select>, SelectError> {
    let mut tree: Option<Box<Select>> = None;

    let mut i = start;

    let mut token = String::new();
    let mut unary_operators: Vec<Operator> = Vec::new();
    let mut binary_operator: Option<Operator> = None;

    while i < end {
        let c = expression.chars().nth(i).unwrap();

        match c {
            '(' => {
                let new_end = match find_parenthesis(expression, i, end) {
                    Some(x) => x,
                    None => {
                        panic!("Groan error. Parentheses should be balanced but they are not.")
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
                // unknown operator
                if operator.is_none() {
                    return Err(SelectError::InvalidOperator("".to_string()));
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
    for _ in unary.iter() {
        parsed = Box::from(Select::Not(parsed));
    }
    unary.clear();

    // apply the previous binary operator
    if let Some(op) = binary {
        if let Some(t) = tree {
            match op {
                Operator::And => Ok(Some(Box::from(Select::And(t, parsed)))),
                Operator::Or => Ok(Some(Box::from(Select::Or(t, parsed)))),
                Operator::Not => panic!(
                    "Groan error. Somehow, NOT operator is being treated as binary operator."
                ),
            }
        } else {
            Err(SelectError::MissingArgument("".to_string()))
        }
    // or create a new tree
    } else {
        if tree.is_some() {
            panic!("Groan error. No binary operator detected but the tree already exists.")
        }
        Ok(Some(parsed))
    }
}

fn find_operator(string: &str, op_symbol: char, start: usize) -> Option<Operator> {
    if string.get(start + 1..start + 2) == Some(op_symbol.to_string().as_str()) {
        match op_symbol {
            '&' => Some(Operator::And),
            '|' => Some(Operator::Or),
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
        "(name W OW HW1 HW2 OH2 H1 H2 or resname SOL WAT HOH OHH TIP T3P T4P T5P T3H W TIP3 TIP4 SPC SPCE)",
    );
    macros.insert(
        "@ion",
        "(name NA NA+ CL CL- K K+ SOD CLA CA CA2+ MG ZN CU1 CU LI RB CS F BR I OH Cal IB+ or resname ION NA NA+ CL CL- K K+ SOD CLA CA CA2+ MG ZN CU1 CU LI RB CS F BR I OH Cal IB+)",
    );
    macros.insert(
        "@dna",
        "(resname DA DG DC DT DA5 DG5 DC5 DT5 DA3 DG3 DC3 DT3 DAN DGN DCN DTN)",
    );
    macros.insert(
        "@rna",
        "(resname A U C G RA RU RC RG RA5 RT5 RU5 RC5 RG5 RA3 RT3 RU3 RC3 RG3 RAN RTN RUN RCN RGN)",
    );

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

    for c in string.chars() {
        if c == '\'' || c == '"' {
            inside = !inside;
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

            Ok(Select::ResidueName(token[1..].to_vec()))
        }
        "name" | "atomname" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            Ok(Select::AtomName(token[1..].to_vec()))
        }
        "resid" | "resnum" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            let range = numbers::parse_numbers(&token[1..])?;
            Ok(Select::ResidueNumber(Group::fix_atom_ranges(
                range,
                usize::MAX,
            )))
        }
        "serial" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            let range = numbers::parse_numbers(&token[1..])?;
            Ok(Select::GmxAtomNumber(Group::fix_atom_ranges(
                range,
                usize::MAX,
            )))
        }

        "atomid" | "atomnum" => {
            if token.len() <= 1 {
                return Err(SelectError::EmptyArgument("".to_string()));
            }

            let range = numbers::parse_numbers(&token[1..])?;
            Ok(Select::AtomNumber(Group::fix_atom_ranges(
                range,
                usize::MAX,
            )))
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

            Ok(Select::GroupName(token[1..].to_vec()))
        }
        // it is not necessary to provide group identifier for groups
        _ => Ok(Select::GroupName(token)),
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

                match parse_query(query) {
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

                match parse_query(query) {
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
        Select::ResidueName(vec!["LYS".to_string()])
    );
    parsing_success!(
        multiple_resname,
        "resname  LYS   CYS   LEU  POPE ",
        Select::ResidueName(vec![
            "LYS".to_string(),
            "CYS".to_string(),
            "LEU".to_string(),
            "POPE".to_string()
        ])
    );
    parsing_success!(
        simple_atomname,
        "name BB",
        Select::AtomName(vec!["BB".to_string()])
    );
    parsing_success!(
        multiple_atomname,
        "  name   BB SC1   PO4   C1 C2B",
        Select::AtomName(vec![
            "BB".to_string(),
            "SC1".to_string(),
            "PO4".to_string(),
            "C1".to_string(),
            "C2B".to_string()
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
        Select::ResidueNumber(vec![(45, usize::MAX - 1)])
    );

    parsing_success!(
        open_ended_range_7,
        "serial  >  44",
        Select::GmxAtomNumber(vec![(45, usize::MAX - 1)])
    );

    parsing_success!(
        open_ended_range_8,
        "atomnum  >=44",
        Select::AtomNumber(vec![(44, usize::MAX - 1)])
    );

    parsing_success!(
        open_ended_range_9,
        "resid  >=  44   ",
        Select::ResidueNumber(vec![(44, usize::MAX - 1)])
    );

    parsing_success!(
        open_ended_range_10,
        "serial 1 4 > 50 7-11",
        Select::GmxAtomNumber(vec![(1, 1), (4, 4), (7, 11), (51, usize::MAX - 1)])
    );

    parsing_success!(
        open_ended_range_11,
        "resid   1   4 >50 7 - 11",
        Select::ResidueNumber(vec![(1, 1), (4, 4), (7, 11), (51, usize::MAX - 1)])
    );

    parsing_success!(
        open_ended_range_12,
        "atomid 1 4 >= 50 7-11",
        Select::AtomNumber(vec![(1, 1), (4, 4), (7, 11), (50, usize::MAX - 1)])
    );

    parsing_success!(
        open_ended_range_13,
        "serial 1 4 >=50 7-11",
        Select::GmxAtomNumber(vec![(1, 1), (4, 4), (7, 11), (50, usize::MAX - 1)])
    );

    parsing_success!(
        open_ended_range_14,
        "serial 1 4>=50 7-11",
        Select::GmxAtomNumber(vec![(1, 1), (4, 4), (7, 11), (50, usize::MAX - 1)])
    );

    parsing_success!(
        open_ended_range_15,
        "resid 4>50",
        Select::ResidueNumber(vec![(4, 4), (51, usize::MAX - 1)])
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
        Select::ResidueNumber(vec![(5, 5), (21, usize::MAX - 1)])
    );

    parsing_success!(
        open_ended_range_21,
        "serial 5>= 20",
        Select::GmxAtomNumber(vec![(5, 5), (20, usize::MAX - 1)])
    );

    parsing_success!(
        open_ended_range_22,
        "serial 24-37<=10 55",
        Select::GmxAtomNumber(vec![(1, 10), (24, 37), (55, 55)])
    );

    parsing_success!(
        open_ended_range_23,
        "serial 24-37<=10-13to17>88",
        Select::GmxAtomNumber(vec![(1, 17), (24, 37), (89, usize::MAX - 1)])
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
        Select::Not(Box::from(Select::ResidueNumber(vec![(21, usize::MAX - 1)])))
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
            "LYS".to_string(),
            "SER".to_string()
        ])))
    );
    parsing_success!(
        simple_not_name,
        " not name   BB SC1",
        Select::Not(Box::from(Select::AtomName(vec![
            "BB".to_string(),
            "SC1".to_string()
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
            "Protein".to_string(),
            "Membrane".to_string()
        ])))
    );
    parsing_success!(
        simple_not_explicit,
        "! Protein Membrane",
        Select::Not(Box::from(Select::GroupName(vec![
            "Protein".to_string(),
            "Membrane".to_string()
        ])))
    );
    parsing_success!(
        simple_not_not,
        " not ! name BB SC1",
        Select::Not(Box::from(Select::Not(Box::from(Select::AtomName(vec![
            "BB".to_string(),
            "SC1".to_string()
        ])))))
    );
    parsing_success!(
        simple_not_not_not,
        "!!! name BB SC1",
        Select::Not(Box::from(Select::Not(Box::from(Select::Not(Box::from(
            Select::AtomName(vec!["BB".to_string(), "SC1".to_string()])
        ))))))
    );

    parsing_success!(
        not_parentheses_1,
        "! (name BB or resname LYS)",
        Select::Not(Box::from(Select::Or(
            Box::from(Select::AtomName(vec!["BB".to_string()])),
            Box::from(Select::ResidueName(vec!["LYS".to_string()]))
        )))
    );

    parsing_success!(
        not_parentheses_2,
        "not(name BB or resname LYS) ||serial 1-5",
        Select::Or(
            Box::from(Select::Not(Box::from(Select::Or(
                Box::from(Select::AtomName(vec!["BB".to_string()])),
                Box::from(Select::ResidueName(vec!["LYS".to_string()]))
            )))),
            Box::from(Select::GmxAtomNumber(vec![(1, 5)]))
        )
    );

    parsing_success!(
        not_parentheses_3,
        "(name BB or resname LYS) || ! serial 1-5",
        Select::Or(
            Box::from(Select::Or(
                Box::from(Select::AtomName(vec!["BB".to_string()])),
                Box::from(Select::ResidueName(vec!["LYS".to_string()]))
            )),
            Box::from(Select::Not(Box::from(Select::GmxAtomNumber(vec![(1, 5)]))))
        )
    );

    parsing_success!(
        not_parentheses_4,
        "not(name BB or not resname LYS) ||serial 1-5",
        Select::Or(
            Box::from(Select::Not(Box::from(Select::Or(
                Box::from(Select::AtomName(vec!["BB".to_string()])),
                Box::from(Select::Not(Box::from(Select::ResidueName(vec![
                    "LYS".to_string()
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
        Select::GroupName(vec!["Protein".to_string()])
    );
    parsing_success!(
        simple_explicit_group,
        "group Protein",
        Select::GroupName(vec!["Protein".to_string()])
    );
    parsing_success!(
        multiword_group_1,
        "  'Protein Membrane ION'  ",
        Select::GroupName(vec!["Protein Membrane ION".to_string()])
    );
    parsing_success!(
        multiword_group_2,
        "\"Protein   Membrane  ION  \"",
        Select::GroupName(vec!["Protein   Membrane  ION  ".to_string()])
    );
    parsing_success!(
        multiword_explicit_group,
        "  group  '  Protein Membrane ION'  ",
        Select::GroupName(vec!["  Protein Membrane ION".to_string()])
    );
    parsing_success!(
        multiple_groups,
        "Protein Membrane ION",
        Select::GroupName(vec![
            "Protein".to_string(),
            "Membrane".to_string(),
            "ION".to_string()
        ])
    );
    parsing_success!(
        multiple_explicit_groups,
        "group Protein Membrane ION",
        Select::GroupName(vec![
            "Protein".to_string(),
            "Membrane".to_string(),
            "ION".to_string()
        ])
    );
    parsing_success!(
        multiple_multiword_groups,
        "Protein 'Membrane Only POPC  ' ' ION with  Water'",
        Select::GroupName(vec![
            "Protein".to_string(),
            "Membrane Only POPC  ".to_string(),
            " ION with  Water".to_string()
        ])
    );
    parsing_success!(
        multiple_multiword_explicit_groups,
        "group Protein 'Membrane Only POPC  ' ' ION with  Water'",
        Select::GroupName(vec![
            "Protein".to_string(),
            "Membrane Only POPC  ".to_string(),
            " ION with  Water".to_string()
        ])
    );
    parsing_success!(
        hyphen_group,
        "Protein-No1",
        Select::GroupName(vec!["Protein-No1".to_string()])
    );
    parsing_success!(
        hyphen_multiword_group,
        "'Protein - No1'",
        Select::GroupName(vec!["Protein - No1".to_string()])
    );
    parsing_success!(
        group_with_unfortunate_name,
        "group resname",
        Select::GroupName(vec!["resname".to_string()])
    );
    parsing_success!(
        group_with_very_unfortunate_name,
        "group group",
        Select::GroupName(vec!["group".to_string()])
    );

    parsing_success!(
        simple_parentheses,
        "(resname LYS)",
        Select::ResidueName(vec!["LYS".to_string()])
    );

    // residue names
    parsing_success!(
        resnames_or_1,
        "resname POPE  'LYS' LEU or resname POPG",
        Select::Or(
            Box::from(Select::ResidueName(vec![
                "POPE".to_string(),
                "LYS".to_string(),
                "LEU".to_string()
            ])),
            Box::from(Select::ResidueName(vec!["POPG".to_string()]))
        )
    );

    parsing_success!(
        resnames_or_2,
        "resname POPE   LYS LEU ||resname POPG",
        Select::Or(
            Box::from(Select::ResidueName(vec![
                "POPE".to_string(),
                "LYS".to_string(),
                "LEU".to_string()
            ])),
            Box::from(Select::ResidueName(vec!["POPG".to_string()]))
        )
    );

    parsing_success!(
        resnames_and_1,
        "resname 'POPE' LYS LEU and resname 'POPG'",
        Select::And(
            Box::from(Select::ResidueName(vec![
                "POPE".to_string(),
                "LYS".to_string(),
                "LEU".to_string()
            ])),
            Box::from(Select::ResidueName(vec!["POPG".to_string()]))
        )
    );

    parsing_success!(
        resnames_and_2,
        "resname POPE LYS LEU&& resname POPG",
        Select::And(
            Box::from(Select::ResidueName(vec![
                "POPE".to_string(),
                "LYS".to_string(),
                "LEU".to_string()
            ])),
            Box::from(Select::ResidueName(vec!["POPG".to_string()]))
        )
    );

    parsing_success!(
        resnames_complex,
        "resname POPE 'LYS'   LEU && resname POPG ||resname LYS",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![
                    "POPE".to_string(),
                    "LYS".to_string(),
                    "LEU".to_string()
                ])),
                Box::from(Select::ResidueName(vec!["POPG".to_string()]))
            )),
            Box::from(Select::ResidueName(vec!["LYS".to_string()]))
        )
    );

    parsing_success!(
        resnames_complex_par_1,
        "((resname POPE   LYS LEU&& resname POPG)) or resname LYS",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec![
                    "POPE".to_string(),
                    "LYS".to_string(),
                    "LEU".to_string()
                ])),
                Box::from(Select::ResidueName(vec!["POPG".to_string()]))
            )),
            Box::from(Select::ResidueName(vec!["LYS".to_string()]))
        )
    );

    parsing_success!(
        resnames_complex_par_2,
        "resname 'POPE' LYS LEU&& (resname POPG   ||   resname LYS)",
        Select::And(
            Box::from(Select::ResidueName(vec![
                "POPE".to_string(),
                "LYS".to_string(),
                "LEU".to_string()
            ])),
            Box::from(Select::Or(
                Box::from(Select::ResidueName(vec!["POPG".to_string()])),
                Box::from(Select::ResidueName(vec!["LYS".to_string()]))
            )),
        )
    );

    // atom names
    parsing_success!(
        names_or_1,
        "name BB  'PO4' D2A or name W",
        Select::Or(
            Box::from(Select::AtomName(vec![
                "BB".to_string(),
                "PO4".to_string(),
                "D2A".to_string()
            ])),
            Box::from(Select::AtomName(vec!["W".to_string()]))
        )
    );

    parsing_success!(
        names_or_2,
        "name BB   PO4 D2A ||atomname W",
        Select::Or(
            Box::from(Select::AtomName(vec![
                "BB".to_string(),
                "PO4".to_string(),
                "D2A".to_string()
            ])),
            Box::from(Select::AtomName(vec!["W".to_string()]))
        )
    );

    parsing_success!(
        names_and_1,
        "atomname 'BB' PO4 D2A and name 'W'",
        Select::And(
            Box::from(Select::AtomName(vec![
                "BB".to_string(),
                "PO4".to_string(),
                "D2A".to_string()
            ])),
            Box::from(Select::AtomName(vec!["W".to_string()]))
        )
    );

    parsing_success!(
        names_and_2,
        "atomname BB PO4 D2A&& atomname W",
        Select::And(
            Box::from(Select::AtomName(vec![
                "BB".to_string(),
                "PO4".to_string(),
                "D2A".to_string()
            ])),
            Box::from(Select::AtomName(vec!["W".to_string()]))
        )
    );

    parsing_success!(
        names_complex,
        "atomname BB 'PO4'   D2A && name W ||name PO4",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::AtomName(vec![
                    "BB".to_string(),
                    "PO4".to_string(),
                    "D2A".to_string()
                ])),
                Box::from(Select::AtomName(vec!["W".to_string()]))
            )),
            Box::from(Select::AtomName(vec!["PO4".to_string()]))
        )
    );

    parsing_success!(
        names_complex_par_1,
        "((name BB   PO4 D2A&& name W)) or name PO4",
        Select::Or(
            Box::from(Select::And(
                Box::from(Select::AtomName(vec![
                    "BB".to_string(),
                    "PO4".to_string(),
                    "D2A".to_string()
                ])),
                Box::from(Select::AtomName(vec!["W".to_string()]))
            )),
            Box::from(Select::AtomName(vec!["PO4".to_string()]))
        )
    );

    parsing_success!(
        names_complex_par_2,
        "atomname 'BB' PO4 D2A&& (name W   ||   atomname PO4)",
        Select::And(
            Box::from(Select::AtomName(vec![
                "BB".to_string(),
                "PO4".to_string(),
                "D2A".to_string()
            ])),
            Box::from(Select::Or(
                Box::from(Select::AtomName(vec!["W".to_string()])),
                Box::from(Select::AtomName(vec!["PO4".to_string()]))
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
                Box::from(Select::ResidueName(vec!["POPE".to_string(), "LYS".to_string(), "LEU".to_string()])),
                Box::from(Select::AtomName(vec!["BB".to_string(), "PO4".to_string(), "D2A".to_string()]))
            )),
            Box::from(Select::AtomName(vec!["C1A".to_string()]))
        )),
        Box::from(Select::And(
            Box::from(Select::ResidueName(vec!["ION".to_string()])),
            Box::from(Select::AtomName(vec!["CL".to_string()]))
        ))
    ));

    // residue names with residue numbers
    parsing_success!(complex_resnames_resid, "(resname 'POPE'  LYS LEU && (resid   15 22-25 33)or(resid 5 to 10) ) or(resname ION&& resnum 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec!["POPE".to_string(), "LYS".to_string(), "LEU".to_string()])),
                Box::from(Select::ResidueNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::ResidueNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::ResidueName(vec!["ION".to_string()])),
            Box::from(Select::ResidueNumber(vec![(6, 6)]))
        ))
    ));

    // residue names with atom numbers
    parsing_success!(complex_resnames_serial, "(resname 'POPE'  LYS LEU && (serial   33 22 -25 15) || atomid 5  to  10 ) ||(resname ION&& atomnum 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec!["POPE".to_string(), "LYS".to_string(), "LEU".to_string()])),
                Box::from(Select::GmxAtomNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::AtomNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::ResidueName(vec!["ION".to_string()])),
            Box::from(Select::AtomNumber(vec![(6, 6)]))
        ))
    ));

    // residue names with group names
    parsing_success!(complex_resnames_group, "(resname 'POPE'  LYS LEU && (Protein)or'Charged   Membrane' ) ||(resname ION&& group Membrane ION W) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec!["POPE".to_string(), "LYS".to_string(), "LEU".to_string()])),
                Box::from(Select::GroupName(vec!["Protein".to_string()]))
            )),
            Box::from(Select::GroupName(vec!["Charged   Membrane".to_string()]))
        )),
        Box::from(Select::And(
            Box::from(Select::ResidueName(vec!["ION".to_string()])),
            Box::from(Select::GroupName(vec!["Membrane".to_string(), "ION".to_string(), "W".to_string()]))
        ))
    ));

    // atom names with residue numbers
    parsing_success!(complex_names_resid, "(name 'BB'  PO4 D2A && (resid   33 22-25 15) || resid 5 to 10 ) or(atomname NA&& resnum 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::AtomName(vec!["BB".to_string(), "PO4".to_string(), "D2A".to_string()])),
                Box::from(Select::ResidueNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::ResidueNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::AtomName(vec!["NA".to_string()])),
            Box::from(Select::ResidueNumber(vec![(6, 6)]))
        ))
    ));

    // atom names with atom numbers
    parsing_success!(complex_names_serial, "(name 'BB'  PO4 D2A && (atomnum   15 22- 25 33) or serial 5 to 10 ) ||(atomname NA&& atomid 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::AtomName(vec!["BB".to_string(), "PO4".to_string(), "D2A".to_string()])),
                Box::from(Select::AtomNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::GmxAtomNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::AtomName(vec!["NA".to_string()])),
            Box::from(Select::AtomNumber(vec![(6, 6)]))
        ))
    ));

    // atom names with group names
    parsing_success!(complex_names_group, "(name 'BB'  PO4 D2A and(group Protein)|| 'Charged   Membrane'   W   )or(atomname NA&& Membrane) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::AtomName(vec!["BB".to_string(), "PO4".to_string(), "D2A".to_string()])),
                Box::from(Select::GroupName(vec!["Protein".to_string()]))
            )),
            Box::from(Select::GroupName(vec!["Charged   Membrane".to_string(), "W".to_string()]))
        )),
        Box::from(Select::And(
            Box::from(Select::AtomName(vec!["NA".to_string()])),
            Box::from(Select::GroupName(vec!["Membrane".to_string()]))
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
                Box::from(Select::GroupName(vec!["Protein".to_string(), "Membrane".to_string()])),
                Box::from(Select::ResidueNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::ResidueNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::GroupName(vec!["Water with Ions".to_string()])),
            Box::from(Select::ResidueNumber(vec![(6, 6)]))
        ))
    ));

    parsing_success!(complex_serial_group, "(Protein Membrane && (serial   15 22-25 33) || atomid 5 to 10 ) ||('Water with Ions' && atomnum 6) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::GroupName(vec!["Protein".to_string(), "Membrane".to_string()])),
                Box::from(Select::GmxAtomNumber(vec![(15, 15), (22, 25), (33, 33)]))
            )),
            Box::from(Select::AtomNumber(vec![(5, 10)]))
        )),
        Box::from(Select::And(
            Box::from(Select::GroupName(vec!["Water with Ions".to_string()])),
            Box::from(Select::AtomNumber(vec![(6, 6)]))
        ))
    ));

    parsing_success!(complex_not_1, "! (name 'BB'  PO4 D2A && (atomnum   15 22- 25 33) or serial 5 to 10 ) ||(atomname NA&& atomid 6) ", 
    Select::Or(
        Box::from(Select::Not(
            Box::from(Select::Or(
                Box::from(Select::And(
                    Box::from(Select::AtomName(vec!["BB".to_string(), "PO4".to_string(), "D2A".to_string()])),
                    Box::from(Select::AtomNumber(vec![(15, 15), (22, 25), (33, 33)]))
                )),
                Box::from(Select::GmxAtomNumber(vec![(5, 10)]))
            ))
        )),
        Box::from(Select::And(
            Box::from(Select::AtomName(vec!["NA".to_string()])),
            Box::from(Select::AtomNumber(vec![(6, 6)]))
        ))
    ));

    parsing_success!(complex_not_2, "(resname 'POPE'  LYS LEU and(name   BB PO4   D2A) || name C1A ) ||not(resname ION&& name CL) ", 
    Select::Or(
        Box::from(Select::Or(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec!["POPE".to_string(), "LYS".to_string(), "LEU".to_string()])),
                Box::from(Select::AtomName(vec!["BB".to_string(), "PO4".to_string(), "D2A".to_string()]))
            )),
            Box::from(Select::AtomName(vec!["C1A".to_string()]))
        )),
        Box::from(Select::Not(
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec!["ION".to_string()])),
                Box::from(Select::AtomName(vec!["CL".to_string()]))
            ))
        ))
    ));

    parsing_success!(complex_not_3, "  not(!(resname 'POPE' LYS LEU and(!name   BB PO4   D2A) || not name C1A ) ||(resname ION&& name CL) )", 
    Select::Not(
        Box::from(Select::Or(
            Box::from(Select::Not(
                Box::from(Select::Or(
                    Box::from(Select::And(
                        Box::from(Select::ResidueName(vec!["POPE".to_string(), "LYS".to_string(), "LEU".to_string()])),
                        Box::from(Select::Not(
                            Box::from(Select::AtomName(vec!["BB".to_string(), "PO4".to_string(), "D2A".to_string()]))
                        ))
                    )),
                    Box::from(Select::Not(Box::from(Select::AtomName(vec!["C1A".to_string()]))))
                ))
            )),
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec!["ION".to_string()])),
                Box::from(Select::AtomName(vec!["CL".to_string()]))
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
                        "POPE".to_string(),
                        "LYS".to_string(),
                        "LEU".to_string()
                    ])),
                    Box::from(Select::Not(Box::from(Select::AtomName(vec![
                        "BB".to_string(),
                        "PO4".to_string(),
                        "D2A".to_string()
                    ]))))
                )),
                Box::from(Select::Not(Box::from(Select::AtomName(vec![
                    "C1A".to_string()
                ]))))
            )))),
            Box::from(Select::And(
                Box::from(Select::ResidueName(vec!["ION".to_string()])),
                Box::from(Select::AtomName(vec!["CL".to_string()]))
            ))
        )))
    );

    parsing_success_string!(complex_parentheses_1,
        "!(!(name BB and resid 15 to 18) || ((resname   POPE POPG &&name PO4  )or not(name C1A||(serial 5 to 12 or group 'Protein 2' Membrane and !resid 1 2 3) )))",
        "Not(Or(Not(And(AtomName([\"BB\"]), ResidueNumber([(15, 18)]))), Or(And(ResidueName([\"POPE\", \"POPG\"]), AtomName([\"PO4\"])), Not(Or(AtomName([\"C1A\"]), And(Or(GmxAtomNumber([(5, 12)]), GroupName([\"Protein 2\", \"Membrane\"])), Not(ResidueNumber([(1, 3)]))))))))"
    );

    parsing_success!(
        keywords_in_quotes,
        "group 'Protein and  Membrane' or Membrane ' W or not ION' \"Lipids   to Count\"",
        Select::Or(
            Box::from(Select::GroupName(vec!["Protein and  Membrane".to_string()])),
            Box::from(Select::GroupName(vec![
                "Membrane".to_string(),
                " W or not ION".to_string(),
                "Lipids   to Count".to_string()
            ]))
        )
    );

    parsing_success!(
        coalesced_keywords_1,
        "resname BB andserial 1 2 3",
        Select::ResidueName(vec![
            "BB".to_string(),
            "andserial".to_string(),
            "1".to_string(),
            "2".to_string(),
            "3".to_string()
        ])
    );

    parsing_success!(
        coalesced_keywords_2,
        "resname BB or notserial 1 2 3",
        Select::Or(
            Box::from(Select::ResidueName(vec!["BB".to_string()])),
            Box::from(Select::GroupName(vec![
                "notserial".to_string(),
                "1".to_string(),
                "2".to_string(),
                "3".to_string()
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
            .map(|&s| s.to_string())
            .collect()
        )
    );
}

#[cfg(test)]
mod fail_tests {
    use super::*;

    macro_rules! parsing_fails {
        ($name:ident, $expression:expr, $variant:path) => {
            #[test]
            fn $name() {
                let query = $expression;

                match parse_query(query) {
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

        match parse_query(query) {
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

    parsing_fails!(deep_error_1, "!(!(name BB and resid 15 to 18) || ((resname   POPE POPG &&name PO4  )or not(name C1A||(serial x to 12 or group 'Protein 2' Membrane and !resid 1 2 3) )))", 
    SelectError::InvalidNumber);
}
