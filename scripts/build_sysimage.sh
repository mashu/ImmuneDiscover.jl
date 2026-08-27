#!/usr/bin/env sh
# Build a local sysimage (immunediscover + dependencies as native code).
# Output: build/immunediscover.so — scripts/run.sh loads it when present.
set -eu
SCRIPT_DIR="$(CDPATH= cd -- "$(dirname "$0")" && pwd)"
REPO_ROOT="$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)"
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
julia --project="$BUILD_ENV" scripts/build_sysimage.jl "$@"
