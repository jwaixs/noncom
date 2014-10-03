load('noncom.sage')

q = var('q')
c1, c2 = var('c1 c2')
d1, d2 = var('d1 d2')
l1, l2 = var('l1, l2')

E1 = ncvar('E1')
E2 = ncvar('E2')
F1 = ncvar('F1')
F2 = ncvar('F2')
K1 = ncvar('K1')
K2 = ncvar('K2')
K1inv = ncvar('K1inv')
K2inv = ncvar('K2inv')
B1c = ncvar('B1c')
B2c = ncvar('B2c')
B1d = ncvar('B1d')
B2d = ncvar('B2d')
Kl = ncvar('Kl')

def part1(fterm):
    msplit = fterm.split()[0][1]

    if len(msplit) == 0:
        return Kl
    
    mainF = msplit[-1]

    if str(mainF) == str(F1):
        return Kl*prod(msplit[0:-1])*B1d + d1*part2(prod(msplit[0:-1]), E2*K1inv)
    else:
        return Kl*prod(msplit[0:-1])*B2d + d2*part2(prod(msplit[0:-1]), E1*K2inv)

def part2(fterm, eterm):
    '''l1 = (lambda, alpha_1), l2 = (lambda, alpha_2).'''
    if fterm == 1:
        if str(eterm.split()[0][1][0]) == str(E1):
            return 1/c2*q**(l1+l2)*Kl*F2 - 1/c2*q**l1*B2c*Kl
        else:
            return 1/c1*q**(l1+l2)*Kl*F1 - 1/c1*q**l2*B1c*Kl

    fsplit = fterm.split()[0][1]
    esplit = eterm.split()[0][1]
    mainE = esplit[0]
    mainF = fsplit[-1]

    if str(mainE) == str(E1) and str(mainF) == str(F1):
        return q*part2(prod(fsplit[0:-1]), E1*K2inv)*F1 \
            - prod(fsplit[0:-1])*(K1 - K1inv)*(1/(q-q**(-1)))*K2inv
    if str(mainE) == str(E1) and str(mainF) == str(F2):
        return q**(-2)*part2(prod(fsplit[0:-1]), E1*K2inv)*F2
    if str(mainE) == str(E2) and str(mainF) == str(F1):
        return q**(-2)*part2(prod(fsplit[0:-1]), E2*K1inv)*F1
    if str(mainE) == str(E2) and str(mainF) == str(F2):
        return q*part2(prod(fsplit[0:-1]), E2*K1inv)*F2 - (K2 - K2inv)*(1/(q-q**(-1)))*K1inv
                
