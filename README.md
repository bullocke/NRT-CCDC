#Continuous Forest Monitoring (CFM)
This repository contains code to perform forest monitoring in near-real time. Monitoring depends on existing time series models to perform predictions. Right now, the only supported model results are from the CCDC algorithm contained in the Yet Another Time Series Model (YATSM) package from Chris Holden (https://github.com/ceholden/yatsm). 

#Installation
It is recommended that CFM be installed in a virtual Python environment using virtualenv. 

With virtualenv installed:

```
virtualenv venv
```

The dependencies can then be installed using pip:

```
pip install -r requirements/requirements.txt
```

Once the dependencies are installed, CFM can be installed:

```
pip install .
``` 
