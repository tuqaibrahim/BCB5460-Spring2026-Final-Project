#!/usr/bin/env bash
set -euo pipefail

mkdir -p reference/gencode_release_scan
mkdir -p metadata

IDS_REGEX='AC024651|AC244131|AL031767|AL732292|FP671120|FP236383|AC129926|AC026469'

out="metadata/gencode_release_missing_id_scan.tsv"
echo -e "release\tid_prefix\tcount" > "$out"

for rel in {19..46}; do
    echo "===== Checking GENCODE release $rel ====="

    url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${rel}/gencode.v${rel}.annotation.gtf.gz"
    gtf_gz="reference/gencode_release_scan/gencode.v${rel}.annotation.gtf.gz"

    if [ ! -f "$gtf_gz" ]; then
        echo "Downloading $url"
        wget -q -c -O "$gtf_gz" "$url" || {
            echo "Release $rel download failed or not available"
            rm -f "$gtf_gz"
            continue
        }
    fi

    for id in AC024651 AC244131 AL031767 AL732292 FP671120 FP236383 AC129926 AC026469; do
        count=$(zgrep -c "$id" "$gtf_gz" || true)
        echo -e "${rel}\t${id}\t${count}" >> "$out"
    done
done

cat "$out"
