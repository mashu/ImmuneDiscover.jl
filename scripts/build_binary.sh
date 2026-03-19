#!/usr/bin/env sh
# Build immunediscover standalone binary for the current platform.
# Requires: Julia 1.9+
# Usage: ./scripts/build_binary.sh [output_dir]
#   If output_dir is omitted, uses build/immunediscover_app/
#   Do NOT pass "build" as output_dir (that folder holds the build environment).

set -eu
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

BUILD_ENV="build"
if [ ! -f "$BUILD_ENV/Project.toml" ]; then
  echo "Creating $BUILD_ENV/Project.toml with PackageCompiler..."
  mkdir -p "$BUILD_ENV"
  cat > "$BUILD_ENV/Project.toml" << 'EOF'
name = "immunediscover-build"

[deps]
PackageCompiler = "9b87118b-4619-50d2-8e1e-99f35a4d4d9d"

[compat]
julia = "1.9"
EOF
fi

julia --project="$BUILD_ENV" -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
julia --project="$BUILD_ENV" scripts/build_binary.jl "$@"
