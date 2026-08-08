#!/bin/bash
# Run CAVER 3.0 tunnel analysis for every CYP2B6 system.
#
# Requires: caver_analysis/snapshots/<SYS>/snapshots/*.pdb (prep_caver_snapshots.py)
#           caver_analysis/starting_points.tsv
# Reads:    env CAVER_HOME, CAVER_JAR, JAVA (defaults set below)
# Usage:    ./run_all_caver.sh [P]   (P = parallel jobs, default 3)

set -u

CAVER_HOME="${CAVER_HOME:-/tmp/caver/caver_3.0.3/caver}"
CAVER_JAR="${CAVER_JAR:-$CAVER_HOME/caver.jar}"
JAVA="${JAVA:-/opt/anaconda3/envs/cyp2b6/bin/java}"
P="${1:-3}"

BASE="$(cd "$(dirname "$0")/.." && pwd)"
CONF="$BASE/caver_analysis/config_caver_template.txt"
POINTS="$BASE/caver_analysis/starting_points.tsv"
OUTBASE="$BASE/caver_analysis/caver_results"
LOGBASE="$BASE/caver_analysis/logs"
mkdir -p "$OUTBASE" "$LOGBASE"

if [ ! -f "$POINTS" ]; then
    echo "ERROR: $POINTS missing (run prep_caver_snapshots.py first)" >&2
    exit 1
fi

run_one() {
    sys="$1"; x="$2"; y="$3"; z="$4"
    conf_tmp="$LOGBASE/config_${sys}.txt"
    snaps="$BASE/caver_analysis/snapshots/$sys/snapshots"
    outdir="$OUTBASE/$sys"
    if [ ! -d "$snaps" ]; then
        echo "SKIP $sys: no snapshots"
        return
    fi
    if [ -f "$outdir/summary.txt" ]; then
        echo "SKIP $sys: already done ($outdir/summary.txt)"
        return
    fi
    sed -e "s/__X__/$x/g" -e "s/__Y__/$y/g" -e "s/__Z__/$z/g" \
        "$CONF" > "$conf_tmp"
    mkdir -p "$outdir"
    echo "== $sys: starting point ($x $y $z) =="
    if timeout 5400 "$JAVA" -Xmx4000m -cp "$CAVER_HOME/lib" \
            -jar "$CAVER_JAR" \
            -home "$CAVER_HOME" \
            -pdb "$snaps" \
            -conf "$conf_tmp" \
            -out "$outdir" > "$LOGBASE/caver_$sys.log" 2>&1; then
        echo "DONE $sys"
    else
        echo "FAIL $sys (rc=$?)"
    fi
}

export -f run_one
export CAVER_HOME CAVER_JAR JAVA BASE CONF OUTBASE LOGBASE

# shellcheck disable=SC2162
tail -n +2 "$POINTS" | while read sys x y z; do
    echo "$sys $x $y $z"
done | xargs -P "$P" -n 4 bash -c 'run_one "$@"' _

echo "CAVER driver finished."
