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

Note that the default content type is `text/html`and hence if you don't set it to `image/svg+xml`, the html `img` tag wouldn't display the svg data.

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

### Adding more options 

With that we realize two things in our web service we need to improve. The caller should be able to define the size of the image. On top of that I personally really prefer to have all structures shown in approximately the same size. Of course large structures need to shrink to fit in the image but small ones should not become overly huge. Since RDKit 2020.x.x this is now possible by setting the `fixedBondLength` drawing options. It makes sense that the caller can also control this parameter.

Making these adjustments and setting a grey background on the test html page will lead to this output:

![benzene_white_background](/assets/img/benzene_white_background.png)

This is to show another problem. By default the svg image has a white background. Personally I would prefer that to be transparent by default. I don't really see the use case for actually having a background color. Hence as a next step the method is enhanced to generate a transparent background. Having this configurable isn't really needed for our use-case. If a background is needed, one can just set the `img` tags `background-color` via css which is much simpler. We also change the font size to 0.8 as the default atom labels are far too small for my taste.

```python
    @cherrypy.expose
    def getStructureImage(self, smiles, image_width = 200, image_height = 100, max_bond_length=15):
    
        m = Chem.MolFromSmiles(smiles)
        mol = rdMolDraw2D.PrepareMolForDrawing(m) 
        
        drawer = rdMolDraw2D.MolDraw2DSVG(image_width, image_height)
        drawer.SetFontSize(0.8)
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

And done! 

## Not so fast there!

If you set a low bar, then yes you are done. But I can tell you, you don't want me as a customer for your software product because mostly likely, I will make it sound worse than it is. At least I try to apply the same high standards to my own creations.

So what is still wrong with the above implementation?

- It's limited to SMILES
- The original orientation of the structure drawing is lost (2D coordinates)
- If it works with 2D /molfile, we should add the bond scaling feature from a [previous post](2020-08-06-depiction-scaling-rdkit.md) 

- No Exception handling

#### Multiple chemical structure formats

We address the first 3 points by changing the first method argument from `smiles` to `identifier`. The goal is to accept multiple types of molecule data. For that we create a new method for converting the string data from the request into an RDKit molecule. I'm simply posting the method below and then explain it afterwards.

```python
def string_to_molecule(self, identifier) -> Mol:
    if identifier is None:
        raise ValueError("Identifier can not be None.")
    # value is a string of type smiles, molfile, inchi, base64-encoded binary RDKit molecule
    # convert it to a rdkit Mol object
    if isinstance(identifier, str):
        mol = (Chem.MolFromSmiles(identifier)
               or Chem.MolFromMolBlock(identifier)
               or Chem.inchi.MolFromInchi(identifier)
               #Converts base64 string created from a binary RDKit molecule representation
               #back into a rdkit molecule
               #base64.b64encode(molecule.ToBinary())
               or Chem.Mol(base64.b64decode(identifier)))
        if mol is None:
            raise ValueError("Input string can not be converted to valid molecule")
        else:
            return mol
    elif isinstance(identifier, Mol):
        return identifier
    else:
        raise TypeError("Expected str but got {}".format(type(identifier))) 
```

The method allows us to create a RDKit molecule from SMILES, InChI, molfile and more esoteric but in my opinion very interesting!!! from a base64-encoded RDKit molecule. This is actually far more efficient than using a molfile. Note that the "base64-decoding" must be tried last as it will throw an exception if the input is invalid base64 data. To make this more robust, there are some basic checks done and exceptions are thrown in case of invalid input. These exceptions are then caught in the calling code (shown later).

#### Bond scaling

We just reuse the code from my [previous post](2020-08-06-depiction-scaling-rdkit.md) on bond scaling and add it as a new method.

```python
def scaleBondLength(self, mol):
    
    default_bond_length = 1.5 #rdkit default bond length
    
    bonds = mol.GetBonds()
    
    if len(bonds) == 0:
        return mol
    
    total = 0
    for bond in bonds:
        
        ai = bond.GetBeginAtomIdx()
        aj = bond.GetEndAtomIdx()
        bl = AllChem.GetBondLength(mol.GetConformer(), ai, aj)
        if not math.isnan(bl):
            total += bl
                    
    avg_bl = (total / len(bonds))
    if avg_bl > 0.0:
        
        bl_ratio = default_bond_length / avg_bl

        tm = np.zeros((4,4),np.double)
        for i in range(3): 
            tm[i,i] = bl_ratio

        AllChem.TransformMol(mol,tm)
        return mol 
```

#### Putting it all together

After above changes the `getStructureImage`method now as below.

```python
@cherrypy.expose
def getStructureImage(self, identifier, image_width = 200, image_height = 100, max_bond_length=15):

    try:
        m = self.string_to_molecule(identifier)
    except (TypeError, ValueError, binascii.Error) as err:
        cherrypy.response.headers['Content-Type'] = 'text/plain'
        raise cherrypy.HTTPError(400, message="The submitted molecule string data is invalid.")
        
    mol = rdMolDraw2D.PrepareMolForDrawing(m)
    mol = self.scaleBondLength(mol)
    #...rest if method stays the same    
```

The input is checked if it is an acceptable format to create an RDKit molecule. If not an HTTP response with client error status 400 is returned. Then the molecules is prepared for drawing and the bond length is scaled. The rest of the method is unchanged.

We can look at the result of the service with below simple html page:

```html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">  
</head>
<body>
<img src="http://localhost:8080/getStructureImage?identifier=c1ccccc1" style="margin-top:10px;background-color: red">
<img src="http://localhost:8080/getStructureImage?identifier=%0A%20%20ChemDraw08142013292D%0A%0A%20%207%20%207%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200999%20V2000%0A%20%20%20-1.0717%20%20%20%200.4125%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20-1.0717%20%20%20-0.4125%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20-0.3572%20%20%20-0.8250%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%200.3572%20%20%20-0.4125%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%200.3572%20%20%20%200.4125%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20-0.3572%20%20%20%200.8250%20%20%20%200.0000%20C%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%20%20%201.0717%20%20%20%200.8250%20%20%20%200.0000%20O%20%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%20%200%0A%20%201%20%202%20%202%20%200%20%20%20%20%20%20%0A%20%202%20%203%20%201%20%200%20%20%20%20%20%20%0A%20%203%20%204%20%202%20%200%20%20%20%20%20%20%0A%20%204%20%205%20%201%20%200%20%20%20%20%20%20%0A%20%205%20%206%20%202%20%200%20%20%20%20%20%20%0A%20%206%20%201%20%201%20%200%20%20%20%20%20%20%0A%20%205%20%207%20%201%20%200%20%20%20%20%20%20%0AM%20%20END%0A">
<img src="http://localhost:8080/getStructureImage?identifier=776t3gAAAAALAAAAAAAAAAAAAAALAAAACwAAAIABBgBgAAAAAQMIACgAAAADAgZAKAAAAAMEBkAoAAAAAwQGQGgAAAADAwEGQGgAAAADAwEGQCgAAAADBAZAaAAAAAMDAQYAaAAAAAMDAQgAKAAAAAMCCABoAAAAAwEBCwABAAECIAIDaAwDBGgMBAVoDAUGaAwGB2gMBgggCAkoAgMKIAcCaAwUAQYCBwYFBAMXAAAAABY=">
<img src="http://localhost:8080/getStructureImage?identifier=ergerg">
</body>
```

The first image is generated from SMILES and has a red background color. The second image is generated from a molfile which must be url-escaped. The third image is a base64-encoded RDKit molecule (see how much less space it uses compared to molfile!) and the last one is simply non-sense to trigger an error response.

#### Additional Remarks

The method can be further extended to allow the caller to set even more options like settings for `PrepareMoleculeForDrawing` or the font size. I would also expect someone to say: 

> "Why do I need this. I can just use the new RDKit JS wrappers."

True. But if you go the web service route you will profit from the browsers built-in image caching. This can be a blessing or a burden. Burden because one needs to take care not to show outdated images. To have some control,  you can add `cache-control` header to the response:

```python
cherrypy.response.headers['Cache-Control'] = 'max-age=3600'
```

In this case the image will be cached in the browser for 1 hour. 

However what also speaks against building the images client-side? In a real application you would most likely push your structure id and not the structure itself and have the image generation service be able to fetch the structure from the database using this id. In that case you can also implement cache revalidation by checking directly in the database if the structure has changed and accordingly return early (304 Not Modified). This will require setting the <pre><a href="https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Last-Modified"> Last-Modified</a></pre> or <pre><a href="https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/ETag">ETag</a></pre> header on the image response. On top of that I'm not a fan of moving too much "compute" stuff onto the client if it can easily be done server-side. 



