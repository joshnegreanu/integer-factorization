class ModInt:
    def __init__(self, val, mod):
        self.mod = mod
        self.val = val % mod

    def __add__(self, other):
        if type(other) == ModInt:
            if self.mod != other.mod:
                raise TypeError("Modulos between ModInt must match.")
            return (self.val + other.val) % self.mod
        
        if type(other) == int:
            return ModInt(self.val + other, self.mod)

        else:
            raise TypeError("Type must be ModInt or int.")
        
    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if type(other) == ModInt:
            if self.mod != other.mod:
                raise TypeError("Modulos between ModInt must match.")
            return ModInt(self.val - other.val, self.mod)
        
        if type(other) == int:
            return ModInt(self.val - other, self.mod)

        else:
            raise TypeError("Type must be ModInt or int.")
    
    def __rsub__(self, other):
        return self.__sub__(other)

    def __mul__(self, other):
        if type(other) == ModInt:
            if self.mod != other.mod:
                raise TypeError("Modulos between ModInt must match.")
            return ModInt(self.val * other.val, self.mod)
        
        if type(other) == int:
            return ModInt(self.val * other, self.mod)

        else:
            raise TypeError("Type must be ModInt or int.")
    
    def __rmul__(self, other):
        return self.__mul__(other)
        
    def __truediv__(self, other):
        if type(other) == ModInt:
            if self.mod != other.mod:
                raise TypeError("Modulos between ModInt must match.")
        
        if type(other) == int:
            other = ModInt(other, self.mod)
        
        if gcd(other.val, self.mod) != 1:
            raise TypeError("Divisor is not invertible.")
        
        return ModInt(self.val * xgcd(other.val, self.mod)[0], self.mod)
    
    def __rtruediv__(self, other):
        if gcd(self.val, self.mod) != 1:
            raise TypeError("Divisor is not invertible.")
        
        return ModInt(other * xgcd(self.val, self.mod)[0], self.mod)
    
    def __neg__(self):
        return ModInt(-self.val, self.mod)
    
    def __pow__(self, x):
        r = 1
        k = self.val
        while x != 0:
            if x % 2 == 1:
                r = (r*k) % self.mod
            k = (k^2) % self.mod
            x >>= 1
        
        return ModInt(r, self.mod)

    def __lt__(self, other):
        if self.mod != other.mod:
            raise TypeError("Modulos between ModInt must match.")
        
        return True if self.val < other.val else False

    def __le__(self, other):
        if self.mod != other.mod:
            raise TypeError("Modulos between ModInt must match.")
        
        return True if self.val <= other.val else False
    
    def __gt__(self, other):
        if self.mod != other.mod:
            raise TypeError("Modulos between ModInt must match.")
        
        return True if self.val > other.val else False
    
    def __ge__(self, other):
        if self.mod != other.mod:
            raise TypeError("Modulos between ModInt must match.")
        
        return True if self.val >= other.val else False
    
    def __eq__(self, other):
        if self.mod != other.mod:
            raise TypeError("Modulos between ModInt must match.")
        
        return True if self.val == other.val else False
    
    def __ne__(self, other):
        if self.mod != other.mod:
            raise TypeError("Modulos between ModInt must match.")
        
        return True if self.val != other.val else False
    
    def __str__(self):
        return str(self.val)

    def __int__(self):
        return int(self.val)

def gcd(x, y):
    while x != 0:
        x_ = x
        x = y % x
        y = x_
    
    return y

def xgcd(x, y):
    a, a_ = 0, 1
    b, b_ = 1, 0
    
    while x != 0:
        q = y // x
        r = y % x
        
        m = a-a_*q
        n = b-b_*q

        y, x = x, r
        a, b = a_, b_
        a_, b_ = m, n
    
    return a, b

a = ModInt(3,10)
b = ModInt(4,10)
print(-a)