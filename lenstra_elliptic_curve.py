import random
from number_theory import ModInt, gcd

def E_add(P, Q, a, N):
    (x1, y1) = P
    (x2, y2) = Q

    if P == None:
        return Q
    
    if Q == None:
        return P
    
    if P == Q:
        d = gcd(2*int(y1), N)
        if d != 1:
            return (False, None, d)

        l = (3*x1**2 + a)/(2*y1)
        
    else:
        d = gcd(int(x1)-int(x2), N)
        if d != 1:
            return (False, None, d)
        
        l = (y1-y2)/(x1-x2)
       
    v = y1-l*x1
    x_ = l**2-x1-x2
    R = (x_, -l*x_**3-v)
    return (True, R, None)

def lenstra_elliptic_curve(N, B, K):
    for _ in range(K):
        P = (ModInt(0,N), ModInt(1,N))
        while True:
            a = ModInt(random.randint(2,N), N)
            if gcd(int(4*a**3 + 27), N) == 1:
                break
        
        P_ = P
           
        for b in range(B):
            if P == None:
                    break
            for _ in range(b):
                if P == None:
                    break
                (work, P, d) = E_add(P, P_, a, N)
                if not work:
                    if 1 < d and d < N:
                        return d
    
    return None