#!/bin/bash

# Runs cargo check for all possible combinations of feature flags.

set -e

# list of features
features=("molly" "parallel" "chemfiles" "serde" "no-xdrfile")

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

echo "Running: cargo check --no-default-features"
cargo check --no-default-features

# run `cargo check` for each combination
for combo in "${combinations[@]}"; do
  if [[ -n "$combo" ]]; then
    echo "Running: cargo check --no-default-features --features $combo"
    cargo check --no-default-features --features "$combo"
  fi
done