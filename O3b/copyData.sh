#!/bin/bash

zenodo_get -d 10.5281/zenodo.5546663 -w files_tmp
grep nocosmo.h5 files_tmp > files
rm files_tmp
cat files | while read f; do
    filename=`echo $f | sed -e 's:.*/::g'`
    if [[ -f "$filename" ]]; then
        echo "$filename exists. Skipping..."
    else
        echo "Downloading $f"
        curl -O $f done
    fi
done
