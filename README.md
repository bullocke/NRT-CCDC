#Near-Real Time CCDC (NRT-CCDC)
This repository contains code to perform forest monitoring in near-real time. Monitoring depends on existing time series models to perform predictions. Right now, the only supported model results are from the CCDC algorithm contained in the Yet Another Time Series Model (YATSM) package from Chris Holden (https://github.com/ceholden/yatsm). 

#Installation
It is recommended that NRT-CCDC be installed in a virtual Python environment using virtualenv. 

With virtualenv installed:

```
virtualenv venv
```

The dependencies can then be installed using pip:

```
pip install -r requirements/requirements.txt
```

Once the dependencies are installed, NRT-CCDC can be installed:

```
pip install .
``` 

#Model creation
NRT-CCDC relies on existing time series models to create predictions that are compared with near-real time data. At the moment, these models can be created using the YATSM package (https://github.com/ceholden/yatsm), or when using MODIS data, a specific version of YATSM altered to filter images by view angle (https://github.com/bullocke/yatsm/tree/NRT). 

