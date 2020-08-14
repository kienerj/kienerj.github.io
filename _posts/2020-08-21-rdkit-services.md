---
layout: post
title: Creating an RDKit webservice
subtitle: for chemical structure image generation
tags: [rdkit, python, cherrypy]
comments: true
---

In my [previous post](2020-08-14-cherrypy.md) I outlined how to install and setup CherryPy for creating web services. In this post we will create a basic web service for generating chemical structure images from the given input. 

## Testing the setup

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

you should see the InChI displayed for benzene. With that we are certain CherryPy and RDKit are installed properly and you can launch the application.

### Request Parameters

The smart reader obviously immediately realized that the python method arguments are automatically taken from the request parameters. If you omit a required parameter (method argument with no default), then CherryPy will throw an exception. If the request has additional unspecified parameters, CherryPy will also throw an exception. If you want to have optional parameters, the method argument must have a default value. This will become relevant once we start building our image generation web service.

## Building the web service

The goal is to generate an svg image of the passed-in chemical structure (or identifier!) and return the svg so that it can be shown in a html img tag.

### First basic implementation

To our `rdkit_services.py`file we need to add an additional import:

```
from rdkit.Chem.Draw import rdMolDraw2D
```

We create a very basic drawing method:

```python
@cherrypy.expose
def getStructureImage(self, smiles):    
    m = Chem.MolFromSmiles(smiles)
    mol = rdMolDraw2D.PrepareMolForDrawing(m) 
    drawer = rdMolDraw2D.MolDraw2DSVG(200, 200)         
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()  
    cherrypy.response.headers['Content-Type'] = 'image/svg+xml'        
    return svg.encode('utf8')
```
The method accepts a smiles as parameter and returns by default an 200x200 image. The most important aspect here are the last two lines of the method. The first one set the correct content type and that trigger a peculiarity of CherryPy as it disables automatic encoding of text into bytes. Therefore before returning the svg we must explicitly encode it. If you forget the encoding when setting a content type that isn't text you will see this error:

```
ValueError: Page handlers MUST return bytes. Use tools.encode if you wish to return unicode.
```

Note that the default content type is `text/html`and hence if you don't set it to `image/svg+xml`, the html `img` tag wouldn't display the svg data.

We can test that above method is working correctly with a trivial html test page:

```html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">  
</head>
<body>
<img src="http://localhost:8080/getStructureImage?smiles=c1ccccc1">
</body>
```

Opening this in a browser will display a huge benzene ring

![benzene_huge](/assets/img/benzene_huge.png)

With that we realize two things in our web service we need to improve. The caller should be able to define the size of the image. On top of that I personally really prefer to have all structures shown in approximately the same size. Of course large structures need to be shrinked to fit in the image but small ones should not become overly huge. Since RDKit 2020.x.x this is now possible by setting the `fixedBondLength` drawing options. It makes sense that the caller can also control this parameter.

Making these adjustments and setting a grey background on the test html page will lead to this output:

![benzene_white_background](/assets/img/benzene_white_background.png)

This is to show another problem. By default the svg image has a white background. Personally I would prefer that to be transparent by default. I don't really see the use case for actually having a background color. Hence as a next step the method is enhanced to generate a transparent background. Having this configurable isn't really needed for our use-case. If a background is needed, one can just set the `img` tags `background-color` via css which is much simpler.

```python
    @cherrypy.expose
    def getStructureImage(self, smiles, image_width = 200, image_height = 100, max_bond_length=15):
    
        m = Chem.MolFromSmiles(smiles)
        mol = rdMolDraw2D.PrepareMolForDrawing(m) 
        
        drawer = rdMolDraw2D.MolDraw2DSVG(image_width, image_height)         
        opts = drawer.drawOptions()
        opts.fixedBondLength = max_bond_length
        opts.bgColor = None
        opts.clearBackground = False
        
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        svg = drawer.GetDrawingText()  
        cherrypy.response.headers['Content-Type'] = 'image/svg+xml'        
        return svg.encode('utf8')
```

Setting the`img` tags `background-color` to red with the transparent image gives this result:

![benzene_red](/assets/img/benzene_red.png)