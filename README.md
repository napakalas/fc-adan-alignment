# fc-adan-alignment

This is a script to align FC map vascular to ADAN

Requirements:

- python = "^3.10"
- pandas = "^1.5.3"
- notebook = "^6.4.12"
- jupyterlab = "^3.5.0"
- anatomy-lookup = {url = "https://github.com/napakalas/anatomy-lookup/releases/download/v0.0.5a/anatomy_lookup-0.0.5-py3-none-any.whl"}

The easiest way to run this code:

```
poetry install

jupyter notebook
```

open `fc_adan_alignement.ipynb`, and make sure it uses the correct kernel.

If you don't like poetry, that's fine, create a new virtual environment, activate it, and then pip install all required packages.
