#!/bin/sh
# Emit make dependencies relative to src, including headers in subdirectories.
cd src || exit 1
headers=$(find . ../inst/include -type f \( -name '*.h' -o -name '*.hpp' \) -print) || exit 1
printf '%s\n' "$headers" | LC_ALL=C sort |
  sed 's|^\./||; s|^|$(OBJECTS): |'
