# Getting started

To begin using the groan crate in your project, you need to add it as a dependency. This is done the same way you would add any other crate in Rust.

You can add the groan crate to your project by using the command line:

```bash
cargo add groan_rs
```

Alternatively, you can manually add the crate to your Cargo.toml file. Open your Cargo.toml and add the following line under [dependencies]:

```toml
[dependencies]
groan_rs = "0.9"
```

Once you have added the groan crate to your project, the simplest way to import it into your Rust code is by including a *prelude*:

```rust
use groan_rs::prelude::*;
```

*(The prelude of the groan crate includes most of its public structures, methods, and traits, making it very convenient to use. However, using the prelude may reduce some level of control over your namespace, as future versions of the crate might introduce new items that could potentially clash with names in your code.)*

With these steps, you are now ready to start using the groan crate in your project. The next sections will guide you through some basic functionalities and common tasks you can perform with this crate.