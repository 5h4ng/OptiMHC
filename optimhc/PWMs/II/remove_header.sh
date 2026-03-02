#!/bin/bash

find . -type f -name "*.txt" | while read -r file; do
    sed '1,2d' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"
    echo "Processed $file"
done
