// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

extern crate cc;

fn main() {
    if std::env::var("CARGO_FEATURE_NO_XDRFILE").is_ok() {
        // do nothing if `no-xdrfile` feature is enabled
        return;
    } else {
        // compile the xdrfile library
        let source_files = vec![
            "external/xdrfile/xdrfile.c",
            "external/xdrfile/xdrfile_xtc.c",
            "external/xdrfile/xdrfile_trr.c",
            "external/xdrfile/xdrfile_jump.c",
        ];
        cc::Build::new()
            .files(source_files)
            .include("external/xdrfile")
            .warnings(false)
            .flag("-O3")
            .compile("libxdrfile.a");

        println!("cargo:rerun-if-changed=external/xdrfile/xdrfile.c");
        println!("cargo:rerun-if-changed=external/xdrfile/xdrfile.h");
        println!("cargo:rerun-if-changed=external/xdrfile/xdrfile_xtc.c");
        println!("cargo:rerun-if-changed=external/xdrfile/xdrfile_xtc.h");
        println!("cargo:rerun-if-changed=external/xdrfile/xdrfile_trr.c");
        println!("cargo:rerun-if-changed=external/xdrfile/xdrfile_trr.h");
        println!("cargo:rerun-if-changed=external/xdrfile/xdrfile_jump.c");
        println!("cargo:rerun-if-changed=external/xdrfile/xdrfile_jump.h");
    }
}
