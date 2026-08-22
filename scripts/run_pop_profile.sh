#!/usr/bin/env bash
set -u
OUT=pop_profile_output
mkdir -p "$OUT" pop_bin pop_data
cd pop_bin
curl -L --retry 3 --connect-timeout 30 -o pop8tools.tar.gz https://harzing.com/download/pop8tools_linux_ub2510.tar.gz
file pop8tools.tar.gz > ../$OUT/download_file_type.txt 2>&1 || true
tar xzf pop8tools.tar.gz
chmod +x pop8query pop8metrics pop8error || true
./pop8query --version > ../$OUT/version.txt 2>&1 || true
./pop8query --help > ../$OUT/help.txt 2>&1 || true
./pop8query --info --all > ../$OUT/info_all.txt 2>&1 || true
strings ./pop8query | grep -Ei 'profile|profileid|citedid|gsciting|gsprofile|scholar' | sort -u > ../$OUT/strings_profile.txt 2>&1 || true
cd ..

try_query() {
  local name="$1"; shift
  echo "COMMAND: $*" > "$OUT/${name}_command.txt"
  timeout 240 ./pop_bin/pop8query --datadir "$PWD/pop_data" "$@" --format json "$OUT/${name}.json" > "$OUT/${name}_stdout.txt" 2> "$OUT/${name}_stderr.txt"
  local rc=$?
  echo "$rc" > "$OUT/${name}_exit_code.txt"
  return $rc
}

# The current PoP CLI documentation does not expose a dedicated --profileid switch;
# test the supported raw/profile field encodings, stopping after the first nonempty JSON result.
for spec in \
  'profile_raw|--gsprofile|--raw=NSpr644AAAAJ' \
  'profile_author|--gsprofile|--author=NSpr644AAAAJ' \
  'profile_title|--gsprofile|--title=NSpr644AAAAJ' \
  'profile_keywords|--gsprofile|--keywords=NSpr644AAAAJ' \
  'profile_raw_url|--gsprofile|--raw=https://scholar.google.com/citations?user=NSpr644AAAAJ'; do
  IFS='|' read -r name a b <<< "$spec"
  try_query "$name" "$a" "$b" || true
  if [ -s "$OUT/${name}.json" ] && grep -q 'NSpr644AAAAJ\|Multichannel Flexible\|Haoyu Wu' "$OUT/${name}.json" 2>/dev/null; then
    echo "$name" > "$OUT/successful_profile_method.txt"
    cp "$OUT/${name}.json" "$OUT/profile_success.json"
    break
  fi
  sleep 8
done

find "$OUT" -maxdepth 1 -type f -printf '%f %s bytes\n' | sort > "$OUT/manifest.txt"
cat "$OUT/version.txt" || true
cat "$OUT/successful_profile_method.txt" 2>/dev/null || true
