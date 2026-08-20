#!/usr/bin/env bash
# One AlphaGenome panel per gene, on that gene's biggest DEEP-INTRONIC density peak.
#
# --limit 1 per gene on purpose: N panels covering N genes answers a different question from N
# panels covering three peaks in one gene. The peak must sit at least --peakmin bp from any splice
# boundary, so nothing here is a donor/acceptor effect in disguise.
#
# Needs a python with numpy and the alphagenome package. Point PY at it, or let it try `python3`:
#   PY=/path/to/python scripts/run_peak_panels.sh
#   PY=/path/to/python GENES="F8 F9 VWF" scripts/run_peak_panels.sh
#
# The API key is read from .env in the repo root (which is gitignored -- never commit it).
set -u
cd "$(dirname "$0")/.."
PY="${PY:-python3}"
command -v "$PY" >/dev/null 2>&1 || { echo "no python at '$PY' -- set PY=/path/to/python"; exit 1; }
"$PY" -c "import numpy, alphagenome" 2>/dev/null || {
  echo "'$PY' is missing numpy or alphagenome -- set PY to an interpreter that has both"; exit 1; }

[ -f .env ] || { echo "no .env in the repo root; it must define ALPHAGENOME_API_KEY"; exit 1; }
ALPHAGENOME_API_KEY="$(grep -m1 '^ALPHAGENOME_API_KEY' .env | cut -d= -f2- | tr -d "\r\"'")"
export ALPHAGENOME_API_KEY
[ -z "$ALPHAGENOME_API_KEY" ] && { echo "ALPHAGENOME_API_KEY not set in .env"; exit 1; }

GENES="${GENES:-F8 F9 F10 F11 F7 F3 F13B VWF FGB SERPINC1 PROC PROCR SERPINA5 PROZ KLKB1 KNG1}"
PEAKMIN="${PEAKMIN:-300}"
for g in $GENES; do
  echo "############ $g ############"
  "$PY" scripts/make_panels.py --gene "$g" --bucket peak --limit 1 --peakmin "$PEAKMIN" 2>&1 \
    | grep -vE "^\s*$" | tail -20
done
echo "ALL DONE"
