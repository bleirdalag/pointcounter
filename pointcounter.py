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
from modmath import num_sqrts_mod_p, is_square_mod_p, mod_sqrt, mod_div
from primes import primes # a large list of primes

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

    # calculate slope (see http://www.mat.uniroma2.it/~schoof/ctg.pdf pg 222)
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
    # (stored in a dictionary for fastest lookup, with the index saved)
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

    # giant steps: compute R, R+/-Q, R+/-2Q,..., R+/-tQ to find i,j
    # such that R + iQ = jP (Theorem 2.1 in Schoof guarantees this happens)
    
    t = int(round((2 * p**(1/2)) / (2*s + 1))) # num giant steps--approx p^(1/4)
    
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
            m = p + 1 + (2*s+1)*(-i) - j # R-iQ-jP = (p+1)P-i(2s+1)P-jP = mP = 0
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
            # count #E'(F_p): E' is the quadratic twist Y^2 = X^3 + Ag^2X + Bg^3
            twist_points = bsgs(p,(a * pow(g,2,p)) % p,(b * pow(g,3,p)) % p)

            if twist_points != -1:
                num_points = 2*(p+1) - twist_points # #E'(F_p)+#E(F_p) = 2(p+1)

        else:
            num_points = bsgs(p,a,b)

    return num_points

def scale_point(s, pt, a, p):
    """ Given point P and scalar s, return sP -- not used in the above
    code but a useful helper function for testing / checking things

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
    curves = 0
    failed = False
    
    for p in primes:
        for a in range(2,5):
            for b in range(2,5):
                curves += 1
                if naive_count(p,a,b) != count_points(p,a,b):
                    failed = True
                    print("Algorithm failed for prime ",
                          p, ", A = ", a, ", B = ", b)
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
