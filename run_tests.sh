#!/bin/bash

# Runs cargo test for various combinations of feature flags.

set -e

# list of features
features=("molly" "parallel" "chemfiles" "no-xdrfile") # serde only tested once with all features

# function to generate all combinations of features (powerset)
powerset() {
  local array=("$@")
  local result=("")
  for element in "${array[@]}"; do
    local temp=("${result[@]}")
    for subset in "${temp[@]}"; do
      if [[ -z "$subset" ]]; then
        result+=("$element")
      else
        result+=("$subset,$element")
      fi
    done
  done
  echo "${result[@]}"
}

# generate all combinations of features
combinations=($(powerset "${features[@]}"))

echo "Running: cargo test --no-default-features"
cargo test --no-default-features

# run `cargo check` for each combination
for combo in "${combinations[@]}"; do
  if [[ -n "$combo" ]]; then
    echo "Running: cargo test --no-default-features --features $combo"
    cargo test --no-default-features --features "$combo"
  fi
done

# testing serde
echo "Running: cargo test --all-features"
cargo test --all-features

# basic default test again to be super sure
echo "Running: cargo test"
cargo test