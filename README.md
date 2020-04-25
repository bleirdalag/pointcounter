## Pointcounter

Counts the number of points on an elliptic curve over Z/p 

[Demo](http://spynx.de/cgi-bin/ctg.py).

### Usage

To run from terminal:

```$ python3 pointcounter.py p a b```

Where `p`, `a`, `b` parameterize the curve over Z/p with Weierstrass
equation Y^2 = X^3 + aX^2 + b

Must use Python3 or the math breaks!

### Info

To write the code I read the first 8 pages of [this paper by
Schoof][1]. It describes how to use the big-steps-little-steps
algorithm to count points for most curves (section 2), and Mestre's
algorithm to get around cases where big-steps-little-steps fails
(section 3). Mestre's algorithm only works for p > 229, so I also
implemented a naive solution for small primes. I also referred to
[these lecture notes][2] which cover much of the same material as the
Schoof paper to help clarify things. I also wrote code for adding
points on elliptic curves, modular division, and [modular square
roots][4] using the [Tonelli-Shanks algorithm][5], as helper
functions.

To test the code, I wrote a simple oracle that tests the result of the
algorithm against a naively-calculated solution. This oracle goes
through a list of primes between 233 and 15601 - from these the oracle
generates about 2000 curves to compare answers for. I cannot test the
correctness of primes larger than 15601 because the naive solution
runs out of memory, but the BSGS/Mestre algorithm runs reasonably fast
for primes of up to 16 digits. Since it is correct for every curve the
oracle generates, I have reasonable faith that it has been correctly
implemented.

### Sources

[1] http://www.mat.uniroma2.it/~schoof/ctg.pdf
[2] https://math.mit.edu/classes/18.783/2017/LectureNotes8.pdf
[3] https://en.wikipedia.org/wiki/Counting_points_on_elliptic_curves
[4] http://www.math.vt.edu/people/brown/doc/sqrts.pdf
[5] https://en.wikipedia.org/wiki/Tonelli%E2%80%93Shanks_algorithm