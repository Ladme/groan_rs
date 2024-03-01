// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of parser for atom and residue numbers.

use crate::errors::SelectError;

#[derive(Debug, Clone, Copy, PartialEq)]
enum NumberToken {
    Number(usize),
    Range,
    Lower,
    LowerOrEqual,
    Greater,
    GreaterOrEqual,
}

impl NumberToken {
    fn str2number(string: &str) -> Result<Self, SelectError> {
        match string.parse::<usize>() {
            Ok(number) => Ok(NumberToken::Number(number)),
            Err(_) => Err(SelectError::InvalidNumber("".to_string())),
        }
    }

    fn extract_number(self) -> Result<usize, SelectError> {
        match self {
            NumberToken::Number(n) => Ok(n),
            _ => Err(SelectError::InvalidNumber("".to_string())),
        }
    }
}

fn tokenize_numbers(token: &[String]) -> Result<Vec<NumberToken>, SelectError> {
    let mut tokens: Vec<NumberToken> = Vec::new();

    let mut current = String::new();
    let joined = token.join(" ");

    for char in joined.chars() {
        // if '-' is reached, immediately tokenize it as Range
        if char == '-' {
            if !current.is_empty() {
                tokens.push(NumberToken::str2number(&current)?);
                current.clear();
            }

            tokens.push(NumberToken::Range);
        }
        // if '<' or '>' is reached, parse the current fragment
        // store the current character
        else if char == '<' || char == '>' {
            if !current.is_empty() {
                tokens.push(NumberToken::str2number(&current)?);
                current.clear();
            }

            current.push(char);
        }
        // if '=' is reached, either '>' or '<' must have preceeded it
        else if char == '=' {
            if current.is_empty() {
                return Err(SelectError::InvalidNumber("".to_string()));
            }

            current.push(char);
            match current.as_str() {
                ">=" => tokens.push(NumberToken::GreaterOrEqual),
                "<=" => tokens.push(NumberToken::LowerOrEqual),
                _ => return Err(SelectError::InvalidNumber("".to_string())),
            }

            current.clear();
        }
        // if whitespace is reached, parse the current fragment
        else if char.is_whitespace() {
            if !current.is_empty() {
                match current.as_str() {
                    ">" => tokens.push(NumberToken::Greater),
                    "<" => tokens.push(NumberToken::Lower),
                    anything => tokens.push(NumberToken::str2number(anything)?),
                }

                current.clear();
            }
        }
        // if digit is reached, end loading of any ">" or "<" operator
        // save the digit
        else if char.is_ascii_digit() {
            match current.as_str() {
                ">" => {
                    tokens.push(NumberToken::Greater);
                    current.clear();
                }
                "<" => {
                    tokens.push(NumberToken::Lower);
                    current.clear();
                }
                _ => (),
            }

            current.push(char);
        } else {
            return Err(SelectError::InvalidNumber("".to_string()));
        }
    }

    // must be a number
    if !current.is_empty() {
        tokens.push(NumberToken::str2number(&current)?);
    }

    Ok(tokens)
}

pub fn parse_numbers(token: &[String]) -> Result<Vec<(usize, usize)>, SelectError> {
    let tokens = tokenize_numbers(token)?;

    let mut numbers: Vec<(usize, usize)> = Vec::new();
    let mut t = 0;
    while t < tokens.len() {
        let token = tokens[t];

        match token {
            NumberToken::Number(n) => {
                if t + 1 < tokens.len() && tokens[t + 1] == NumberToken::Range {
                    t += 1;
                    continue;
                }

                numbers.push((n, n));
                t += 1;
            }

            NumberToken::Range => {
                if t == 0 || t + 1 == tokens.len() {
                    return Err(SelectError::InvalidNumber("".to_string()));
                }

                let previous = tokens[t - 1].extract_number()?;
                let next = tokens[t + 1].extract_number()?;

                if previous > next {
                    return Err(SelectError::InvalidNumber("".to_string()));
                }

                numbers.push((previous, next));
                t += 2;
            }

            _ => {
                if t + 1 == tokens.len() {
                    return Err(SelectError::InvalidNumber("".to_string()));
                }

                let next = tokens[t + 1].extract_number()?;

                match token {
                    NumberToken::Greater => numbers.push((next + 1, usize::MAX)),
                    NumberToken::GreaterOrEqual => numbers.push((next, usize::MAX)),
                    NumberToken::Lower if next > 1 => numbers.push((1, next - 1)),
                    NumberToken::Lower if next <= 1 => (),
                    NumberToken::LowerOrEqual => numbers.push((1, next)),
                    _ => panic!("FATAL GROAN ERROR | numbers::parse_numbers | Impossible match condition reached."),
                }

                t += 2;
            }
        }
    }

    Ok(numbers)
}
