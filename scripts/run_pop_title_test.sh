#!/usr/bin/env bash
set -u
OUT=pop_title_output
mkdir -p "$OUT" "$HOME/.publish_or_perish/cache"
if [ ! -x pop_bin/pop8query ]; then
  mkdir -p pop_bin
  curl -L --retry 3 --connect-timeout 30 -o pop_bin/pop8tools.tar.gz https://harzing.com/download/pop8tools_linux_ub2510.tar.gz
  (cd pop_bin && tar xzf pop8tools.tar.gz && chmod +x pop8query pop8metrics pop8error)
fi
TITLE='Rheological mechanism of polymer nanocomposites filled with spherical nanoparticles: Insight from molecular dynamics simulation'
echo "$TITLE" > "$OUT/query_title.txt"
# Print syntax and run direct search in multiple output formats to inspect all internal fields.
./pop_bin/pop8query --gscholar --title "$TITLE" --syntax > "$OUT/syntax.txt" 2> "$OUT/syntax_stderr.txt" || true
timeout 300 ./pop_bin/pop8query --direct --gscholar --title "$TITLE" --max 10 --wait 10 --format json "$OUT/title.json" > "$OUT/title_stdout.txt" 2> "$OUT/title_stderr.txt"; echo $? > "$OUT/title_exit.txt"
timeout 300 ./pop_bin/pop8query --direct --gscholar --title "$TITLE" --max 10 --wait 10 --format tdf "$OUT/title.tdf" > "$OUT/title_tdf_stdout.txt" 2> "$OUT/title_tdf_stderr.txt"; echo $? > "$OUT/title_tdf_exit.txt"
timeout 300 ./pop_bin/pop8query --direct --gscholar --title "$TITLE" --max 10 --wait 10 --format csv "$OUT/title.csv" > "$OUT/title_csv_stdout.txt" 2> "$OUT/title_csv_stderr.txt"; echo $? > "$OUT/title_csv_exit.txt"
find "$OUT" -type f -printf '%f %s bytes\n' | sort > "$OUT/manifest.txt"
cat "$OUT/syntax.txt" || true
cat "$OUT/title_stderr.txt" || true
