#!/bin/bash
# mkvenv.sh <lib> <py> <req...>: fresh venv under the repository's tmp/newlib/<lib>/venv, wheels only
L=$1; P=$2; shift 2; B=$(cd "$(dirname "$0")/../../../.." && pwd)/tmp/newlib; export UV_CACHE_DIR=$B/_uvcache
uv venv --python $P $B/$L/venv -q && uv pip install -q --only-binary :all: --python $B/$L/venv/Scripts/python.exe "$@" && uv pip list --python $B/$L/venv/Scripts/python.exe 2>/dev/null | grep -iE "^($(echo "$@" | sed -E 's/==[^ ]*//g; s/ /|/g')|numpy|scipy|pandas) "; du -sh $B/$L/venv
