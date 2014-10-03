# Created by Noud Aldenhoven (2014)

import sys
from numbers import Number

class NonComVar():
    def __init__(self, label):
        self.label = label
        self.var = var(label)
    
    def getLabel(self):
        return self.label

    def getVar(self):
        return self.var

    def printVar(self, newline=True):
        sys.stdout.write(self.getLabel())
        if newline: sys.stdout('\n')
        sys.stdout.flush()

    def __str__(self):
        return self.getLabel()

    def __repr__(self):
        return str(self)

class NonComMonomial():
    def __init__(self, coef, noncomvars):
        self.coef = coef
        self.noncomvars = noncomvars

    def __len__(self):
        return len(self.noncomvars)
    
    def __getitem__(self, key):
        return self.noncomvars[key]

    def __str__(self):
        sret = ''
        if self.getCoef() != 1:
            sret += '(' + str(self.getCoef()) + ')'
            for ncvar in self.getNonComVars():
                sret += '*' + str(ncvar)
        else:
            ncvars = self.getNonComVars()
            for ncvar in ncvars[:-1]:
                sret += str(ncvar) + '*'
            sret += str(ncvars[-1])
    
        return sret

    def __repr__(self):
        return str(self)

    def negatif(self):
        return NonComMonomial(-self.coef, self.noncomvars)

    def getCoef(self):
        return self.coef

    def getNonComVars(self):
        return self.noncomvars

    def multLeft(self, nmmonomial):
        newcoef = self.getCoef() * nmmonomial.getCoef()
        newcomvars = nmmonomial.getNonComVars() + self.getNonComVars()

        return NonComMonomial(newcoef, newcomvars)
    
    def multRight(self, nmmonomial):
        newcoef = self.coef * nmmonomial.getCoef()
        newcomvars = self.getNonComVars() + nmmonomial.getNonComVars()

        return NonComMonomial(newcoef, newcomvars)

    def printMonomial(self, mnewline=True):
        if self.getCoef() != 1:
            sys.stdout.write('(' + str(self.getCoef()) + ')')
            for ncvar in self.getNonComVars():
                sys.stdout.write('*')
                ncvar.printVar(newline=False)
        else:
            ncvars = self.getNonComVars()
            for ncvar in ncvars[:-1]:
                ncvar.printVar(newline=False)
                sys.stdout.write('*')
            ncvars[-1].printVar(newline=False)
        
        if mnewline: sys.stdout.write('\n')

        sys.stdout.flush()

class NonComPolynomial():
    def __init__(self, polynomial=[]):
        self.polynomial = polynomial

    def __str__(self):
        if len(self.getPolynomial()) == 0:
            return '0'
        if len(self.getPolynomial()) > 1:
            return ' + '.join(map(str, self.getPolynomial()))
        else:
            return str(self.getPolynomial()[0])

    def __repr__(self):
        return str(self)

    def __mul__(self, other):
        if isinstance(other, (Number, Integer, Rational, Expression)):
            return self.multCoef(other)
        return self.multRightPolynomial(other)

    def __rmul__(self, other):
        if isinstance(other, (Number, Integer, Rational, Expression)):
            return self.multCoef(other)
        return self.multLeftPolynomial(other)

    def __add__(self, other):
        return self.addPolynomial(other)

    def __sub__(self, other):
        return self.subPolynomial(other)

    def __pow__(self, power):
        ret = 1
        for n in range(power):
            ret *= self
        return ret

    def __getitem__(self, key):
        return self.getPolynomial()[key]

    def getPolynomial(self):
        return self.polynomial

    def addMonomial(self, nmmonomial):
        return NonComPolynomial(self.getPolynomial().append(nmmonomial))

    def addPolynomial(self, npol):
        return NonComPolynomial(self.getPolynomial() + npol.getPolynomial())

    def subPolynomial(self, npol):
        newmonomials = []
        for mon in npol.getPolynomial():
            newmonomials.append(mon.negatif())
        return NonComPolynomial(self.getPolynomial() + newmonomials)
    
    def multCoef(self, coef):
        newmons = []
        
        for mon in self.getPolynomial():
            newmons.append(NonComMonomial(coef*mon.getCoef(), mon.getNonComVars()))

        return NonComPolynomial(newmons)

    def multLeftMonomial(self, nmmonomial):
        newpolynomial = []
        for elm in self.getPolynomial():
            newpolynomial.append(elm.multLeft(nmmonomial))

        return NonComPolynomial(newpolynomial)

    def multRightMonomial(self, nmmonomial):
        newpolynomial = []
        for elm in self.getPolynomial():
            newpolynomial.append(elm.multRight(nmmonomial))

        return NonComPolynomial(newpolynomial)

    def multLeftPolynomial(self, nmpolynomial):
        newpolynomial = []
        for mon in nmpolynomial.getPolynomial():
            newpolynomial += self.multLeftMonomial(mon).getPolynomial()
        
        return NonComPolynomial(newpolynomial)

    def multRightPolynomial(self, nmpolynomial):
        newpolynomial = []
        for mon in nmpolynomial.getPolynomial():
            newpolynomial += self.multRightMonomial(mon).getPolynomial()
        
        return NonComPolynomial(newpolynomial)

    def printPolynomial(self, pnewline=True):
        monomials = self.getPolynomial()

        for mon in monomials[:-1]:
            mon.printMonomial(mnewline=False)
            sys.stdout.write(' + ')
        monomials[-1].printMonomial(mnewline=False)
        
        if pnewline: sys.stdout.write('\n')

        sys.stdout.flush()

    def simplify(self):
        oldpol =  self.getPolynomial()
        usedmonomials = []
        newpol = []

        for m1 in range(len(oldpol)):
            breakout = False
            for mon in usedmonomials:
                if oldpol[m1].getNonComVars() == mon:
                    breakout = True
                    break
            if breakout:
                continue
            newmon = oldpol[m1]
            ncoef = newmon.getCoef()
            for m2 in range(m1+1, len(oldpol)):
                if oldpol[m1].getNonComVars() == oldpol[m2].getNonComVars():
                    ncoef += oldpol[m2].getCoef()
            try:
                ncoef = ncoef.full_simplify()
            except:
                pass
            usedmonomials.append(newmon.getNonComVars())
            if ncoef != 0:
                newpol.append(NonComMonomial(ncoef, newmon.getNonComVars()))

        return NonComPolynomial(newpol)

    def split(self):
        ret1 = []

        for monomial in self.getPolynomial():
            ret2 = []
            for ncvar in monomial.getNonComVars():
                ret2.append(NonComPolynomial([NonComMonomial(1, [ncvar])]))
            ret1.append((monomial.getCoef(), ret2))

        return ret1


def ncvar(label):
    return NonComPolynomial([NonComMonomial(1, [NonComVar(label)])])
            

if __name__ == '__main__':
    print 'loading...'
    q = var('q')
    var1 = NonComVar('a')
    var2 = NonComVar('b')
    mon1 = NonComMonomial(1, [var1])
    mon2 = NonComMonomial((1 - q), [var2])
    pol1 = NonComPolynomial([mon1, mon2])
    pol2 = NonComPolynomial([mon1])
