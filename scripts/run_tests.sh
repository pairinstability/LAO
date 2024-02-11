#!/bin/bash

set -e

ROOT_SRC_DIR="$(cd -P "$(dirname "${BASH_SOURCE[0]}")/.." && pwd )"
BUILD_DIR="$ROOT_SRC_DIR/build"

for target in $(grep -oP 'create_test\((\w+)\)' "$ROOT_SRC_DIR/tests/CMakeLists.txt" | grep -oP '\(\K\w+'); do
    "$BUILD_DIR/$target"
done