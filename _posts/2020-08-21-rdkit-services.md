---
layout: post
title: Creating an RDKit webservice
subtitle: for chemical structure image generation
tags: [rdkit, python, cherrypy]
comments: true
---

In my [previous post](2020-08-14-cherrypy.md) I outlined how to install and setup CherryPy for creating web services. In this post we will create a basic web service for generating chemical structure images from the given input.

## Creating the baseline

In this section we will setup the basic harness and check if the installations are working correctly. We are going to assume that we are running CherryPy in standalone mode so that it works for all readers also on Windows. 

In a directory of your choosing create a new file `rdkit_services.py`. The name is actually not important and you can change it if you want. In that file paste below code.

```python
import cherrypy
from rdkit import Chem


class RDKit_Services(object):

    @cherrypy.expose
    def smilesToInchi(self, smiles):
        m = Chem.MolFromSmiles(smiles)
        cherrypy.response.headers['Content-Type'] = 'text/plain'        
        return Chem.MolToInchi(m)   


if __name__ == '__main__':
    cherrypy.quickstart(RDKit_Services())
```

In a cli activate your cherrypy conda environment

```bash
conda activate cherrypy
```

Change to the directory you created the `rdkit_services.py` file in and run it.

```bash
python rdkit_services.py
```

If you now open a browser and enter below url

http://localhost:8080/smilesToInchi?smiles=c1ccccc1

you should see the InChI displayed for benzol. With that we are certain CherryPy and RDKit are installed properly and you can launch the application.

## Request Parameters

The smart reader obviously immediately realized that the python method arguments are automatically taken from the request parameters. If you omit a required parameter (method argument with no default), then CherryPy will throw an exception. If the request has additional unspecified parameters, CherryPy will also throw an exception. If you want to have optional parameters, the method argument must have a default value. 

