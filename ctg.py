""" Counts the number of points on an elliptic curve over Z/p.

To run from terminal:

$ python3 ctg.py p a b

Where a, b parameterize the curve with Weierstrass equation Y^2 = X^3 + aX^2 + b

Must use Python3 or the math breaks!

Sources:

http://www.mat.uniroma2.it/~schoof/ctg.pdf
https://math.mit.edu/classes/18.783/2017/LectureNotes8.pdf
https://en.wikipedia.org/wiki/Counting_points_on_elliptic_curves
https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm
http://www.math.vt.edu/people/brown/doc/sqrts.pdf
"""

import random

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

def naive_count(p, a, b):
    """Return the number of points on the elliptic curve over F_p given by 
    Weierstrass equation Y^2 = X^3 + aX + b 

    >>> naive_count(5, 1, 1)
    9
    """
    # precompute the number of square roots of each field element
    num_sqrts = num_sqrts_mod_p(p)

    # count number of points with a given x coordinate for all x
    num_pts = 0 
    for x in range(p):
        y_squared = (pow(x,3,p) + ((a*x) % p) + (b % p)) % p
        num_pts += num_sqrts[y_squared]

    return 1 + num_pts # add 1 for the point at infinity


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

    # the below is essentially copied from
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
    d_inv = pow(d, p - 2, p) # uses FLT: a^{p-1} = 1 mod p
    return (q*d_inv) % p

def get_random_point(p,a,b):
    """Returns a random point on the curve Y^2 = X^3 + AX + B over F_p

    >>> possible_pts = [(0,1),(0,4),(2,1),(2,4),(3,1),(3,4),(4,2),(4,3)]
    >>> random_point = get_random_point(5,1,1) 
    >>> random_point in possible_pts
    True
    """
    # pick a random x value until one is found on the curve
    while True:
        x = random.randint(0, p-1)
        y_squared = (pow(x,3,p) + ((a*x) % p) + (b % p)) % p
        if is_square_mod_p(y_squared, p):
            break
        
    y = mod_sqrt(y_squared, p)

    return (x,y)

def add_points(p,q,a,pr):
    """Takes two points as tuples (x, y) or the string "inf"
    and returns their sum under the group law for elliptic curves

    Also needs A from the Weierstrass equation of the curve and the prime

    >>> add_points((0,1),(0,4),1,5)
    'inf'

    >> add_points((0,1),(0,1),1,5)
    (4,2)

    >> add_points((0,1),(4,3),1,5)
    (0,4)
    """
    if p == "inf":
        return q

    if q == "inf":
        return p

    if p == negate_point(q, pr):
        return "inf"

    # calculate slope (http://www.mat.uniroma2.it/~schoof/ctg.pdf pg 222)
    if p == q:
        m = mod_div((3 * pow(p[0],2,pr) + a) % pr, (2 * p[1]) % pr, pr)
    else:
        m = mod_div((q[1] - p[1]) % pr, (q[0] - p[0]) % pr, pr)

    x = (-p[0] - q[0]) % pr + pow(m, 2, pr)
    y = m * (p[0] - x) - p[1]
    
    return (x % pr, y % pr)

def negate_point(pt, p):
    """ 
    >>> negate_point((0,1),5)
    (0, 4)
    """
    if pt == 'inf':
        return 'inf'
    
    return (pt[0], -pt[1] % p)

def bsgs(p,a,b):
    """Return the number of points on the curve using the baby-step-giant-step
    algorithm.

    >>> bsgs(5,1,1) in [9, -1]
    True
    """
    
    P = get_random_point(p,a,b)

    # compute the values m in interval of possible #E(F_p) such that mP = 0
    # by taking baby steps then giant steps
    ms = set()

    # baby steps: calculate the first p^(-4) multiples of P
    # (these are stored in a dictionary for fastest lookup, with the index saved)
    baby_steps = {'inf': 0} 
    s = int(p**(1/4))
    next_multiple = 'inf'
    
    for i in range(1, s+1):
        next_multiple = add_points(next_multiple, P, a, p)
        baby_steps[next_multiple] = i
        baby_steps[negate_point(next_multiple, p)] = -i

    # the final value of next_multiple is sP
    two_sP = add_points(next_multiple, next_multiple, a, p)
    Q = add_points(two_sP, P, a, p) # Q = (2s + 1)P

    # compute R = (p+1)P by using the binary expansion of p+1
    bin_str = bin(p+1)[2:]
    R = 'inf'
    doubled = P
    for bit in reversed(bin_str):
        if bit == '1':
            R = add_points(R, doubled, a, p)
        doubled = add_points(doubled, doubled, a, p)

    # giant steps: compute R, R+/-Q, R+/-2Q,..., R+/-tQ to find i,j such that R + iQ = jP
    # we are guaranteed this happens by Theorem 2.1 in Schoof.
    
    t = int(round((2 * p**(1/2)) / (2*s + 1))) # num "giant steps" to take,approx p^(1/4)
    #assert(t*(2*s+1) + s >= 2 * int(p**(1/2)))
    
    next_Q_multiple = 'inf'

    for i in range(t+1):
        iQ = next_Q_multiple
        R_plus_iQ = add_points(R, iQ, a, p)
        if R_plus_iQ in baby_steps: # then R + iQ = jP
            j = baby_steps[R_plus_iQ]
            m = p + 1 + (2*s+1)*i - j # R + iQ - jP = mP = 0
            ms.add(m)

        R_minus_iQ = add_points(R, negate_point(iQ, p), a, p)
        if R_minus_iQ in baby_steps: # then R - iQ = jP
            j = baby_steps[R_minus_iQ]
            m = p + 1 + (2*s+1)*(-i) - j # R - iQ - jP = (p+1)P - i(2s+1)P - jP = mP = 0
            ms.add(m)

        next_Q_multiple = add_points(next_Q_multiple, Q, a, p)

    if len(ms) == 1:
        return ms.pop()
    else:
        return -1 # failure


def count_points(p,a,b):
    """Uses Big-Steps-Little-Steps and Mestre's algorithm to count points 

    >>> count_points(13,1,0)
    20

    >>> count_points(449,1,0)
    464

    >>> count_points(773,1,2)
    724
    """

    if p < 230:
        return naive_count(p, a, b)

    num_points = bsgs(p,a,b)
    twisted = False

    # find a non-square g
    for n in range(p):
        if not is_square_mod_p(n,p):
            g = n
            break
    
    while num_points == -1: #bsgs failed; try Mestre's algorithm
        # alternate between trying BSGS on E and on E's quadratic twist
        twisted = not twisted

        if twisted:
            # count #E'(F_p): E' is the quadratic twist with eqn Y^2 = X^3 + Ag^2X + Bg^3
            twist_points = bsgs(p,(a * pow(g,2,p)) % p,(b * pow(g,3,p)) % p) # try again

            if twist_points != -1:
                num_points = 2*(p+1) - twist_points # bc #E'(F_p) + #E(F_p) = 2(p+1)

        else:
            num_points = bsgs(p,a,b)

    return num_points

def scale_point(s, pt, a, p):
    """ Return sP -- not used in the above code but useful for testing / checking things

    >>> scale_point(1, (1,0), 1, 5)
    (1, 0)

    >>> scale_point(2, (0,1), 1, 5)
    (4, 2)
    """

    sum = 'inf'
    for i in range(s):
        sum = add_points(sum, pt, a, p)

    return sum

def oracle():
    """Tests that naive_count() and count_points() yield the same answer
    """

    primes = [233, 277, 331, 379, 431, 467, 523, 587, 631, 677, 739, 797, 853, 907, 967, 1019, 1063, 1117, 1187, 1237, 1297, 1367, 1433, 1483, 1543, 1597, 1637, 1709, 1777, 1847, 1901, 1979, 2027, 2087, 2141, 2221, 2281, 2341, 2389, 2447, 2539, 2609, 2671, 2711, 2767, 2833, 2897, 2963, 3037, 3109, 3187, 3253, 3319, 3371, 3457, 3517, 3559, 3623, 3691, 3761, 3823, 3889, 3943, 4019, 4091, 4153, 4229, 4273, 4357, 4441, 4507, 4567, 4643, 4703, 4787, 4861, 4933, 4987, 5039, 5107, 5189, 5273, 5347, 5417, 5477, 5527, 5623, 5669, 5741, 5813, 5861, 5927, 6037, 6091, 6163, 6229, 6299, 6353, 6421, 6521, 6577, 6661, 6719, 6793, 6863, 6947, 6991, 7057, 7151, 7219, 7307, 7393, 7481, 7537, 7583, 7649, 7717, 7793, 7877, 7937, 8039, 8101, 8179, 8243, 8311, 8389, 8467, 8563, 8629, 8693, 8747, 8821, 8887, 8969, 9041, 9127, 9187, 9257, 9337, 9403, 9461, 9521, 9619, 9679, 9749, 9817, 9883, 9949, 10067, 10133, 10181, 10267, 10331, 10427, 10487, 10589, 10651, 10723, 10799, 10883, 10957, 11057, 11117, 11177, 11273, 11351, 11437, 11497, 11593, 11689, 11779, 11831, 11909, 11969, 12043, 12113, 12203, 12269, 12347, 12421, 12491, 12547, 12613, 12689, 12763, 12841, 12919, 12983, 13049, 13147, 13217, 13297, 13381, 13457, 13537, 13627, 13693, 13757, 13831, 13903, 13997, 14071, 14159, 14251, 14347, 14423, 14503, 14561, 14639, 14723, 14771, 14843, 14923, 15013, 15091, 15161, 15241, 15299, 15361, 15439, 15511, 15601]

    curves = 0
    failed = False
    
    for p in primes:
        for a in range(2,5):
            for b in range(2,5):
                curves += 1
                if naive_count(p,a,b) != count_points(p,a,b):
                    failed = True
                    print("Algorithm failed for prime ", p, ", A = ", a, ", B = ", b)
                    print(curves, " curves checked so far")

    print("Checked ", curves, " elliptic curves")
    if not failed:
        print("Algorithm did not fail for any! :)")

if __name__ == "__main__":
    import doctest
    doctest.testmod()

    import sys
    if len(sys.argv) != 4:
        print("Usage: $ python3 ctg.py p A B")
    else:
        p = int(sys.argv[1])
        a = int(sys.argv[2])
        b = int(sys.argv[3])
        print(count_points(p,a,b))
