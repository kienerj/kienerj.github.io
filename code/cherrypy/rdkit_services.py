import cherrypy
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


class RDKit_Services(object):
    
    @cherrypy.expose
    def smilesToInchi(self, smiles):
    
        m = Chem.MolFromSmiles(smiles)
        cherrypy.response.headers['Content-Type'] = 'text/plain'        
        return Chem.MolToInchi(m)
    
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
    


if __name__ == '__main__':
    cherrypy.quickstart(RDKit_Services())