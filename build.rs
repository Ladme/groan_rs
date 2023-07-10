extern crate cc;

use std::io::Result;

fn main() -> Result<()> {
    // compile the xdrfile library
    let source_files = vec![
        "external/xdrfile/xdrfile.c",
        "external/xdrfile/xdrfile_xtc.c",
        "external/xdrfile/xdrfile_trr.c",
    ];
    cc::Build::new()
        .files(source_files)
        .include("external/xdrfile")
        .warnings(false)
        .flag("-O3")
        .compile("libxdrfile.a");

    println!("cargo:rerun-if-changed=src/xdrfile.c");
    println!("cargo:rerun-if-changed=src/xdrfile.h");
    println!("cargo:rerun-if-changed=src/xdrfile_xtc.c");
    println!("cargo:rerun-if-changed=src/xdrfile_xtc.h");
    println!("cargo:rerun-if-changed=src/xdrfile_trr.c");
    println!("cargo:rerun-if-changed=src/xdrfile_trr.h");

    Ok(())
}
