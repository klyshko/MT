from string import rjust

class PDBAtom:
    """ Just a class representing one PDB entry for one atom
    By now has only two methods - constructor (__init__) and stringifier(__str__)
    for input and output respectively
    """
    def __init__ (self, line):
        """ PDBAtom(string)
        Builds new PDBAtom instance using string read from PDB file"""
        self.prec = '%.3f'
        if not PDBAtom.Atomizable(line):
            raise Exception('Trying to initialize PDBAtom from string not starting with "ATOM  "')
        self.record_name = line[0:6]
        self.number      = int(line[6:11])
        self.name        = line[12:16]
        self.altlocation = line[16:17]
        self.residue     = line[17:20]
        self.chainid     = line[21:22]
        self.residueseq  = int(line[22:26])
        self.insertcode  = line[26:27]
        self.x           = float(line[30:38])
        self.y           = float(line[38:46])
        self.z           = float(line[46:54])
        self.occupancy   = float(line[54:60])
        self.tempfactor  = float(line[60:66])
        self.segmentid   = line[72:76]
        self.element     = line[76:78]
#       self.charge      = line[78:80]
    def __str__(self):
        s = 'ATOM  '
        s += rjust(str(self.number)    , 5) + ' '
        s += rjust(self.name[:4]       , 4)
        s += rjust(self.altlocation[:1], 1)
        s += rjust(self.residue[:3]    , 3) + ' '
        s += rjust(self.chainid[:1]    , 1)
        s += rjust(str(self.residueseq), 4)
        s += rjust(self.insertcode[:1] , 1) + '   '
        s += rjust(self.prec%(self.x)     , 8)
        s += rjust(self.prec%(self.y)     , 8)
        s += rjust(self.prec%(self.z)     , 8)
        s += rjust('%.2f'%(self.occupancy) , 6)
        s += rjust('%.2f'%(self.tempfactor), 6) + '      '
        s += rjust(self.segmentid[:4]  , 4)
        s += rjust(self.element[:2]    , 2)
#        s += rjust(self.charge[:2]     , 2)
        return s
    @staticmethod
    def Atomizable(line):
        """Checks whether given line can be used to construct PDBAtom"""
        return line[0:6] == 'ATOM  '

def PDBReadFile(name):
    return [PDBAtom(i) for i in open(name,'r').readlines() if PDBAtom.Atomizable(i)]

