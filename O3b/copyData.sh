#!/bin/bash

zenodo_get -d 10.5281/zenodo.5546663 -w files_tmp
grep nocosmo.h5 files_tmp > files
rm files_tmp
cat files | while read f; do curl -O $f; done;
