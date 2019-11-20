import sys
import os

import mbuild as mb
import foyer


class Propane(mb.Compound):
    def __init__(self):
        super(Propane, self).__init__()
        # PDB contains coordinates and some atom names
        # But we need to touch up the atom names 
        # and add bonds
        mb.load('TraPPE_UA_3_propane_monomer1.pdb', compound=self,
                relative_to_module=self.__module__)
        self.name = 'Pro'
        particles = [a for a in self.particles()]
        particles[0].name = particles[0].name[1:]
        particles[2].name = particles[2].name[1:]
        self.add_bond((particles[0], particles[1]))
        self.add_bond((particles[1], particles[2]))

        # For convenience, each mbcompound carries its xml file
        xml_path = os.path.realpath(
                sys.modules[self.__module__].__file__)
        file_dir = os.path.dirname(xml_path)
        self.xml = os.path.join(file_dir, "TraPPE_UA_3_fully_flexible_propane.xml" )
