#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p bin
flags=(-std=c++17 -pthread -O2 -Wno-comment -I./src/OpenCL/include -L./src/OpenCL/lib)
build() {
	local name=$1; shift
	g++ src/*.cpp -o "bin/$name" "${flags[@]}" "$@" -lOpenCL
}
case "${1:-all}" in
	FluidX3D) build FluidX3D ;;
	FluidX3D_force) build FluidX3D_force -DFORCE_FIELD ;;
	all)
		build FluidX3D &
		a=$!
		build FluidX3D_force -DFORCE_FIELD &
		b=$!
		wait "$a"
		wait "$b"
		;;
	*) echo "unknown target $1" >&2; exit 2 ;;
esac
