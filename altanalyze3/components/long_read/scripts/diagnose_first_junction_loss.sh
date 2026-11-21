#!/usr/bin/env bash
# Read-only diagnostic for the per-library long-read junction matrix.
#
# Usage:
#   bash diagnose_first_junction_loss.sh \
#     <transcript_associations_raw.txt[.gz]> \
#     <transcript_associations.txt[.gz]> \
#     <junction-counts.txt> [isoform-id] [junction-prefix]
#
# The association tables are five-column, tab-delimited files:
# gene, strand, exon structure, molecule/isoform id, source.
# This script deliberately never opens a BAM or modifies any run output.

set -euo pipefail

if (( $# < 3 || $# > 5 )); then
    echo "Usage: $0 <raw-associations> <filtered-associations> <junction-counts.txt> [isoform-id] [junction-prefix]" >&2
    exit 2
fi

raw_associations=$1
filtered_associations=$2
junction_counts=$3
isoform_id=${4:-34754578}
junction_prefix=${5:-ENSG00000213658:E1.1-E2.9}

for path in "$raw_associations" "$filtered_associations" "$junction_counts"; do
    [[ -r "$path" ]] || { echo "Unreadable input: $path" >&2; exit 2; }
done

stream() {
    case $1 in
        *.gz) gzip -cd -- "$1" ;;
        *)    cat -- "$1" ;;
    esac
}

echo "== Raw association for molecule $isoform_id =="
stream "$raw_associations" | awk -F '\t' -v id="$isoform_id" '$4 == id { print }'

echo
echo "== Filtered association used by isoform_junctions_chunk.py =="
stream "$filtered_associations" | awk -F '\t' -v id="$isoform_id" '$4 == id { print }'

echo
echo "== First-junction comparison from filtered association =="
stream "$filtered_associations" | awk -F '\t' -v id="$isoform_id" '
    $4 == id {
        n = split($3, exon, "|")
        if (n < 2) {
            print "No junction can be formed (fewer than two exon tokens)."
            next
        }
        print "Expected first junction: " $1 ":" exon[1] "-" exon[2]
        if (n < 3)
            print "Current [1:] parser output: no junction (only two tokens remain after dropping exon[1])."
        else
            print "Current [1:] parser first junction: " $1 ":" exon[2] "-" exon[3]
    }
'

echo
echo "== Junction-count text search: $junction_prefix =="
stream "$junction_counts" | grep -F -- "$junction_prefix" || true

echo
echo "== Scope in filtered association table =="
stream "$filtered_associations" | awk -F '\t' '
    NF >= 4 && $1 !~ /UNK/ {
        n = split($3, exon, "|")
        if (n >= 2) { models++; suppressed_first_junctions++ }
    }
    END {
        printf("%d non-UNK molecule models contain >=2 exon tokens; [1:] suppresses one first junction from each.\n", models, suppressed_first_junctions)
    }
'

echo
echo "Interpretation: an empty search above is expected if the target is the first adjacent exon pair."
echo "The active exporter is components/long_read/isoform_junctions_chunk.py; it currently uses exon_structure.split('|')[1:]."
