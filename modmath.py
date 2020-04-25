def num_sqrts_mod_p(p):
    """Return a table giving the number of squares mod p of each number
    from 0,...,p-1

    >>> num_sqrts_mod_p(5)
    [1, 2, 0, 0, 2]
    """
    table = [0] * p

    for i in range(p):
        table[pow(i, 2, p)] += 1

    return table

def is_square_mod_p(n, p):
    """Checks if n is a square (mod p) using Euler's criterion 

    >>> is_square_mod_p(0,5)
    True
    >>> is_square_mod_p(1,5)
    True
    >>> is_square_mod_p(2,5)
    False
    """
    if n == 0:
        return True
    else:
        return pow(n, int((p-1)/2), p) == 1

    
def mod_sqrt(a, p):
    """Returns a square root (mod p) of a using Tonelli-Shanks algorithm

    >>> mod_sqrt(2, 113) in [51, 62]
    True
    >>> mod_sqrt(5, 40961) in [19424, 21537]
    True
    """
    assert is_square_mod_p(a, p)

    if a == 0:
        return 0

    # write p-1 = s*2^e for an odd s
    s = p - 1
    e = 0
    while s % 2 == 0:
        s /= 2
        e += 1
    s = int(s)

    # find a number that is not square mod p
    n = 2
    while is_square_mod_p(n, p):
        n += 1

    # the below is mostly copied from
    # http://www.math.vt.edu/people/brown/doc/sqrts.pdf pg 91
    x = pow(a, int((s + 1) / 2), p) # guess of square root
    b = pow(a, s, p) # "fudge factor" - loop invariant is x^2 = ab (mod p)
    g = pow(n, s, p) # used to update x and b
    r = e # exponent - decreases with each update

    while True:
        t = b
        m = 0
        for m in range(r):
            if t == 1:
                break
            t = pow(t, 2, p)

        if m == 0:
            return x

        gs = pow(g, 2 ** (r - m - 1), p)
        g = (gs * gs) % p
        x = (x * gs) % p
        b = (b * g) % p
        r = m

def mod_div(q,d,p):
    """ 
    >>> mod_div(2,3,5)
    4
    """
    # use Fermat's Little Theorem (d^{p-1} = 1 mod p) to calculate 1/d
    d_inv = pow(d, p - 2, p)

    return (q*d_inv) % p
