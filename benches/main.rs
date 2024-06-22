// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

use criterion::{criterion_group, criterion_main, Criterion};
use groan_rs::{
    prelude::TrajMasterRead,
    progress::ProgressPrinter,
    structures::{dimension::Dimension, vector3d::Vector3D},
    system::System,
};

fn benchmark(c: &mut Criterion) {
    let mut system = System::from_file("test_files/example.gro").unwrap();
    system.read_ndx("test_files/index.ndx").unwrap();

    c.bench_function("System::atoms_iter", |b| {
        b.iter(|| {
            std::hint::black_box(
                system
                    .atoms_iter()
                    .fold(0, |sum, atom| sum + atom.get_atom_number()),
            );
        })
    });

    c.bench_function("System::group_iter (Membrane)", |b| {
        b.iter(|| {
            std::hint::black_box(
                system
                    .group_iter("Membrane")
                    .unwrap()
                    .fold(0, |sum, atom| sum + atom.get_atom_number()),
            );
        })
    });

    c.bench_function("System::atoms_iter (creation)", |b| {
        b.iter(|| {
            std::hint::black_box(system.atoms_iter());
        })
    });

    c.bench_function("System::group_iter (Membrane, creation)", |b| {
        b.iter(|| {
            std::hint::black_box(system.group_iter("Membrane").unwrap());
        })
    });

    c.bench_function("System::get_atoms_as_ref and iter", |b| {
        b.iter(|| {
            std::hint::black_box(
                system
                    .get_atoms_as_ref()
                    .iter()
                    .fold(0, |sum, atom| sum + atom.get_atom_number()),
            );
        })
    });

    c.bench_function("System::group_get_center (Protein)", |b| {
        b.iter(|| {
            std::hint::black_box(system.group_get_center("Protein").unwrap());
        })
    });

    c.bench_function("System::group_get_center (Membrane)", |b| {
        b.iter(|| {
            std::hint::black_box(system.group_get_center("Membrane").unwrap());
        })
    });

    c.bench_function("System::atoms_center (Protein, xyz)", |b| {
        b.iter(|| {
            std::hint::black_box(system.atoms_center("Protein", Dimension::XYZ).unwrap());
        })
    });

    c.bench_function("System::atoms_center (Membrane, xyz)", |b| {
        b.iter(|| {
            std::hint::black_box(system.atoms_center("Membrane", Dimension::XYZ).unwrap());
        })
    });

    let translate = Vector3D::new(0.432, 0.721, -0.453);

    c.bench_function("System::atoms_translate (short)", |b| {
        b.iter(|| {
            std::hint::black_box(system.atoms_translate(&translate).unwrap());
        })
    });

    let translate = Vector3D::new(2.432, 1.721, -1.453);

    c.bench_function("System::atoms_translate (medium)", |b| {
        b.iter(|| {
            std::hint::black_box(system.atoms_translate(&translate).unwrap());
        })
    });

    let translate = Vector3D::new(12.432, 8.721, -9.453);

    c.bench_function("System::atoms_translate (long)", |b| {
        b.iter(|| {
            std::hint::black_box(system.atoms_translate(&translate).unwrap());
        })
    });

    c.bench_function("System::xtc_iter (no progress printing)", |b| {
        b.iter(|| {
            std::hint::black_box(system.xtc_iter("test_files/short_trajectory.xtc").unwrap());
        })
    });

    c.bench_function("System::xtc_iter (with progress printing)", |b| {
        b.iter(|| {
            std::hint::black_box(
                system
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .print_progress(ProgressPrinter::new().with_print_freq(1)),
            );
        })
    });
}

criterion_group!(benches, benchmark);
criterion_main!(benches);
