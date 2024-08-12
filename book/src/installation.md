# Installation

You can add the `groan_rs` crate into your project in the same way you do that for all other crates.

You can either use the command line:

```bash
cargo add groan_rs
```

or add the crate into your dependencies in the `.toml` file:

```
[dependencies]
groan_rs = "0.9"
```

Then you can import the crate in your Rust code as:

```rust
use groan_rs::prelude::*;
```