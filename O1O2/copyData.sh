#!/bin/bash

curl -O https://dcc.ligo.org/public/0171/P2000434/001/injections_O1O2an_spin.h5
curl -O https://dcc.ligo.org/public/0157/P1800370/005/GWTC-1_sample_release.tar.gz
tar -xvf GWTC-1_sample_release.tar.gz
rm GWTC-1_sample_release.tar.gz
