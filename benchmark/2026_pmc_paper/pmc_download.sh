#!/usr/bin/env bash

pmcids=(
  PMC13044143
  PMC12992144
  PMC12858432
)

for id in "${pmcids[@]}"; do
  echo "Checking $id ..."

  xml=$(curl -s "https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=${id}")

  # Extract tgz href from XML
  url=$(echo "$xml" | grep -o 'href="[^"]*\.tar\.gz"' | head -n 1 | cut -d'"' -f2)

  if [ -z "$url" ]; then
    echo "No .tar.gz href found for $id"
    continue
  fi

  # Insert /deprecated after /pub/pmc
  deprecated_url=$(echo "$url" | sed 's#ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/#ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/deprecated/#')

  echo "Original:   $url"
  echo "Download:   $deprecated_url"

  wget -O "${id}.tar.gz" "$deprecated_url"
done
