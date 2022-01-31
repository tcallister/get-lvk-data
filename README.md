# get-lvk-data

This repo contains some tools and scripts for downloading and prepping O1, O2, and O3 parameter estimation data for population inference.
To start, make sure [conda](https://docs.conda.io/en/latest/miniconda.html) is installed, then create a python environment via

```bash
conda env create -f environment.yml
```

Activate your new environment with

```bash
conda activate get-lvk-data
```

This will guarantee that you have all the python libraries and command line tools needed.

Next, download PE samples:

```bash
cd O1O2/
. copyData.sh

cd ../O3a/
. copyData.sh

cd ../O3a-GWTC-2.1/
. copyData.sh

cd ../O3b/
. copyData.sh
```

**WARNING**: This will download approximately 25GB of data. Make sure you and your computer are ready.

Finally, start a jupyter session (using the command `jupyter notebook`) to open and run the notebook `preprocess_samples.ipynb`. The result will be  a single file, `sampleDict_FAR_1_in_1_yr.pickle`, that contains PE samples and redshift priors for every event passing a one-per-year false-alarm rate threshold.
