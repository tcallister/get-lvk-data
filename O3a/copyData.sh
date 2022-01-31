#!/bin/bash

curl -O https://dcc.ligo.org/public/0168/P2000217/002/o3a_bbhpop_inj_info.hdf
curl -O https://dcc.ligo.org/public/0169/P2000223/005/all_posterior_samples.tar
tar -xvf all_posterior_samples.tar
rm all_posterior_samples.tar
