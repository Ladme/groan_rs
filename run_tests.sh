#!/bin/bash

# Runs tests for all possible combinations of feature flags.

cargo test && \
cargo test --features molly && \
cargo test --features parallel && \
cargo test --features "molly parallel" && \
cargo test --features serde && \
cargo test --features "serde molly" && \
cargo test --features "serde parallel" && \
cargo test --features "serde parallel molly"