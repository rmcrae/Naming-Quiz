from numpy import *

# W22:  Added physical state
# change .elements from list of tuples to dictionary
# add charges

# define custom record
record = dtype([('symbol', str, 2), ('atWt', float)])
AWdict = dict(loadtxt('AWdata.txt', usecols = (1, 2), dtype = record))
states = ["(g)", "(s)", "(l)", "(aq)", "(graph)"]
class FormulaError(Exception):
    pass

class Formula:
    def __init__(self, formula):
        self.state_formula = formula.strip()
        self.formula = self.state_formula
        self.state = ""
        self.Chemdict = AWdict
        self.punc_dict = {"(":")", "[":"]"}
        for state in states:
            if state in formula:
                self.state = state
                idx = self.state_formula.find(state)
                self.formula = self.state_formula[:idx].strip()
                self.state_formula = self.formula + " " + state
        
        if  self.wrong_bracket_count():
            print( f"Invalid formula: Punctuation error for {formula}")
        elif self.wrong_caps():
            print (f"Invalid formula: Capitalization error for {formula}")
        else:
            self.parse()
    
    def wrong_bracket_count(self):          #ensures same number of [ and ], or ( and )
        return not prod([self.formula.count(key) == self.formula.count(value) for
                key, value in self.punc_dict.items()])
        
    def wrong_caps(self):         #ensures there are not two lower case letters in a row
        out = array([char.islower() for char in self.formula])
        return any(out[1:] * out[:-1])

    def __str__(self):
        return f"{self.state_formula}"
    
    def __repr__(self):
        return f"<Formula: {self.state_formula}>"
    
    def parse(self):
        self._collect_charge()
        self.collect_symbols()
        self.elim_brackets()
        self.elim_numbers()
        self.make_tuple()

    def _collect_charge(self):
        split_result = self.formula.split('^')
        if len(split_result) == 1:
            self.charge = 0
            self.elem_formula = self.formula
            return
        self.elem_formula, self.strCharge = split_result
        if '+' not in self.strCharge and '-' not in self.strCharge:
            print(f"Unrecognizable charge in {self.strCharge}")
            raise FormulaError
        factor = 1 if '+' in self.strCharge else -1
        for strSign in ['+', '-']:
            if self.strCharge.endswith(strSign):
                self.charge = 1 if len(self.strCharge)==1 else int(self.strCharge[:-1])
            elif self.strCharge.startswith('+'):
                self.charge = 1 if len(self.strCharge)==1 else int(self.strCharge[1:])

        # correct sign of charge
        self.charge *= factor
        self.strSign = '+' if factor == 1 else '-'
        self.strCharge = (f"{abs(self.charge)}{self.strSign}" if abs(self.charge) > 1
                          else f"{self.strSign}")
        self.formula = f"{self.elem_formula}^{self.strCharge}"
        self.state_formula = self.formula
        # correct for state if there is one
        if self.state:
            self.state_formula = self.formula + " " + self.state
                
    def collect_symbols(self):
        '''this method collects letters into one or two-letter symbols and
        collects consecutive digits into numbers.  Resulting list of tokens
        is stored in self.atoms'''
        
        x = list(self.elem_formula)
        for i, char in enumerate(x):
            if x[i].islower() or (x[i-1].isdigit() and x[i].isdigit()):
                x[i-1] = x[i-1] + x[i]
                del x[i]
        self.atoms = x
            
    def elim_brackets(self):                           
        '''This method eliminates all brackets and replaces the tokens inside
            brackets with a list multiplied by any number following right bracket'''
        
        for left, right in self.punc_dict.items():                     
            for i, char in enumerate(self.atoms):    
                try:                                
                    start = self.atoms.index(left)
                    stop = self.atoms.index(right)
                except ValueError:
                    break
                if (stop+1) < len(self.atoms) and self.atoms[stop+1].isdigit():
                    num = int(self.atoms[stop+1])
                    del self.atoms[stop+1]
                    self.atoms[start:stop+1] = self.atoms[start:stop+1]*num
                self.atoms.remove(left)
                self.atoms.remove(right)
        return self.atoms

    def elim_numbers(self):                              
        '''This method removes any number n and replaces any preceding token
            with a list of n such tokens'''
        
        for i, char in enumerate(self.atoms):    
            if char.isdigit():                  
                num = int(char)
                del self.atoms[i]
                self.atoms[i-1:i] = self.atoms[i-1:i]*num

    def make_tuple(self):
        '''create list of tuples of elements and numbers of each'''

        self.elements = dict([(elem, self.atoms.count(elem)) for
                         elem in set(self.atoms)])

    @property
    def FW(self):       #gives formula weight of compound
        weight = [self.Chemdict[char]*num for char,num in self.elements.items()]
        return sum(weight)
    
        
    def __add__(self, other):
        return Compound(self.formula + other.formula)

    def __mul__(self, N):
        return Compound(self.formula*N)

    def __rmul__(self, N):
        return Compound(N*self.formula)
    
if __name__ == '__main__':

    formTests = ["NaCl (s)","RuF4","NH4NO3","CH3CH3", "(NH4)3PO4", "Ca3(PO4)2",
             "Na3[Fe(CN)6]","CH3(CH2OCH2)4CHNH4(CH2)2CH3","C12H22O11",
             "Ca3[Fe(CN)6]2", "C108H141N39O31", "C (graph)", "CO2 (g)", "NO3^- (aq)",
                 "Ca^+2", "HPO4^2-"]
    ans = []
    for c in formTests:
        f = Formula(c)
        ans.append(f)
        print( c)
        print( f)
        print( f.FW)
        print( f.elements)
    
