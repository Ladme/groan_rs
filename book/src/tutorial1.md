# Tutorial 1: Analyzing a protein

In the first tutorial of The Groan Book, we will be performing the simplest possible analyses I could come up with. We will load a molecular structure of a protein from a PDB file and calculate distances between various atoms and groups of atoms. We will also look at properties of atoms that are available in the groan crate.

## 0. What will you learn?
- Load a molecular structure.
- Obtain basic information about the system.
- Iterate over atoms.
- Work with atom properties.
- Select atoms and create groups.
- Calculate distances between atoms and groups.
- The basic philosophy behind the groan crate.

## 1. Starting the project
We will start the project as any other Rust project:

```bash
cargo new groan_tutorial_1
cd groan_tutorial_1
```

We then include the groan crate:

```bash
cargo add groan_rs
```


## 2. Getting the input structure
Then, we will be needing some molecular structure to work with. I will be using a structure of human hemoglobin solved by X-ray diffraction and deposited to the RCSB Protein Data Bank under the code 1A3N ([here](https://www.rcsb.org/structure/1a3n)). You can use the same structure or, if you are feeling brave, you can use any other protein or any other molecular system whatsoever.

If you want to work with the same structure, run:
```bash
wget https://files.rcsb.org/download/1A3N.pdb
```

No matter what structure you decide to use, you should place the file containing the structure into the `groan_tutorial_1` directory containing the Rust project.

## 3. Loading the structure using groan
Now we can finally open our favorite IDE or editor and start writing some code. Loading the molecular structure is simple:

```rust
use groan_rs::prelude::*;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
    let system = System::from_file("1A3N.pdb")?;
    Ok(())
}
```

The `System::from_file` will automatically recognize the file format of the input file based on the file extension. In my case, the method will recognize that the structure is a PDB file and will read it as such. In case your file is missing an extension or uses some non-standard extension, you can enforce a specific file format using:

```rust
// this will read the "1A3N.pdb" file as a Gromacs GRO file
// no matter the file extension
let system = System::from_file_with_format("1A3N.pdb", FileType::GRO)?;
```

Note that both `System::from_file` and `System::from_file_with_format` are fallible methods. They return an error if the file could not be found or its content could not be understood. We propagate the error outside the main function using the `?` operator.

## 4. Obtaining information about the structure

## 5. Calculating distances between atoms

## 6. Creating groups

## 7. Simpler ways to calculate distances