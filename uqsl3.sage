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

def join_non_com_term(lterm):
    if len(lterm) == 0:
        return 1
    
    coef, elms = lterm[0]
    ret = coef*prod(elms)

    for coef, elms in lterm[1:]:
        ret += coef*prod(elms)

    return ret

def part1(fterm):
    msplit = fterm.split()[0][1]

    if len(msplit) == 0:
        return Kl
    
    mainF = msplit[-1]

    rterm = 0
    if str(mainF) == str(F1):
        rterm = d1*part2(prod(msplit[0:-1]), E2*K1inv) + Kl*prod(msplit[0:-1])*B1d 
    else:
        rterm = d2*part2(prod(msplit[0:-1]), E1*K2inv) + Kl*prod(msplit[0:-1])*B2d 

    return rterm

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
            - Kl*prod(fsplit[0:-1])*(K1 - K1inv)*(1/(q-q**(-1)))*K2inv
    if str(mainE) == str(E1) and str(mainF) == str(F2):
        return q**(-2)*part2(prod(fsplit[0:-1]), E1*K2inv)*F2
    if str(mainE) == str(E2) and str(mainF) == str(F1):
        return q**(-2)*part2(prod(fsplit[0:-1]), E2*K1inv)*F1
    if str(mainE) == str(E2) and str(mainF) == str(F2):
        return q*part2(prod(fsplit[0:-1]), E2*K1inv)*F2 \
            - Kl*prod(fsplit[0:-1])*(K2 - K2inv)*(1/(q-q**(-1)))*K1inv

def part3(fterm):
    rterm = part1(fterm).split()
    final = (Kl*fterm).split()[0][1]

    while str(rterm[0][1]) != str(final):
        rterm = (rterm[0][0]*part1(prod(rterm[0][1][1:]))).split() + rterm[1:]

    return 1/(1 - rterm[0][0])*join_non_com_term(rterm[1:])
                
klf1 = 1/(1 - (d1*q^(l1 + l2)/c1)) * (Kl*B1d + (-d1*q^l2/c1)*B1c*Kl)
klf2 = 1/(1 - (d2*q^(l1 + l2)/c2)) * (Kl*B2d + (-d2*q^l1/c2)*B2c*Kl)
