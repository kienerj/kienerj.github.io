import cherrypy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdchem import Mol
import base64
import binascii
import numpy as np
import math

class RDKit_Services(object):
    
    @cherrypy.expose
    def smilesToInchi(self, smiles):
    
        m = Chem.MolFromSmiles(smiles)
        cherrypy.response.headers['Content-Type'] = 'text/plain'        
        return Chem.MolToInchi(m)
    
    @cherrypy.expose
    def getStructureImage(self, identifier, image_width = 200, image_height = 100, max_bond_length=15):
    
        try:
            m = self.string_to_molecule(identifier)
        except (TypeError, ValueError, binascii.Error) as err:
            cherrypy.response.headers['Content-Type'] = 'text/plain'
            raise cherrypy.HTTPError(400, message="The submitted molecule string data is invalid.")
            
        mol = rdMolDraw2D.PrepareMolForDrawing(m)
        mol = self.scaleBondLength(mol)    
        
        drawer = rdMolDraw2D.MolDraw2DSVG(image_width, image_height)         
        opts = drawer.drawOptions()
        opts.fixedBondLength = max_bond_length
        opts.bgColor = None
        opts.clearBackground = False
        
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        
        svg = drawer.GetDrawingText()  
        cherrypy.response.headers['Content-Type'] = 'image/svg+xml'
        cherrypy.response.headers['Cache-Control'] = 'max-age=3600'
        return svg.encode('utf8')    
        
    
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


if __name__ == '__main__':
    cherrypy.quickstart(RDKit_Services())