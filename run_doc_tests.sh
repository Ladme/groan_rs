#!/bin/bash

# Runs cargo test --doc for all possible combinations of feature flags.

set -e

# list of features
features=("molly" "parallel" "serde" "no-xdrfile")

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

echo "Running: cargo test --doc --no-default-features"
cargo test --doc --no-default-features

# run `cargo check` for each combination
for combo in "${combinations[@]}"; do
  if [[ -n "$combo" ]]; then
    echo "Running: cargo test --doc --no-default-features --features $combo"
    cargo test --doc --no-default-features --features "$combo"
  fi
done