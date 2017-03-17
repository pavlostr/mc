︠bab0f467-3cb2-46a4-9400-80253ac66b7a︠
# -*- coding: utf-8 -*-


# initialize our finite fields F = GF2 and Phi = GF(2^m)
F.<x> = GF(2);
Phi.<x> = GF(2^m);
codelocators = [x^i for i in range(N)];


# returns a Goppa polynomial in z over field F
def goppapolynomial(F,z):
X = PolynomialRing(F,repr(z)).gen();
return X^(N-K)+X+1;

# checks whether the Goppa polynomial is irreducible
g = goppapolynomial(Phi,z);
if g.is_irreducible():
print 'g(z) =',g,'is irreducible';

# checks that the code locators are not zeroes of g
for i in range(N):
if g(codelocators[i])==Phi(0):
print 'alarm: g(alpha_'+str(i)+')=0';

# set up the parity matrix H_gRS
H_gRS = matrix([[codelocators[j]^(i) for j in range(N)] for i in range(N-K)]);
H_gRS = H_gRS*diagonal_matrix([1/g(codelocators[i]) for i in range(N)]);

# From H_gRS we construct H_Goppa
H_Goppa = matrix(F,m*H_gRS.nrows(),H_gRS.ncols());
for i in range(H_gRS.nrows()):
for j in range(H_gRS.ncols()):
be = bin(eval(H_gRS[i,j].int_repr()))[2:];
be = '0'*(m-len(be))+be;
be = list(be);
H_Goppa[m*i:m*(i+1),j]=vector(map(int,be));

Krnl = H_Goppa.right_kernel();
G_Goppa = Krnl.basis_matrix();

# encodes information words u to code words c = uG_Goppa
def encode(u):
return u*G_Goppa;






#decoding
def split(p):
Phi = p.parent()
p0 = Phi([sqrt(c) for c in p.list()[0::2]]);
p1 = Phi([sqrt(c) for c in p.list()[1::2]]);
return (p0,p1);

# returns the g(z) inverse of polynomial p
def g_inverse(p):
(d,u,v) = xgcd(p,g);
return u.mod(g);

(g0,g1) = split(g);
w = g0*g_inverse(g1);
T = g_inverse(s);
(T0,T1) = split(T+z);
R = T0+w*T1;
(d,u,v) = xgcd(PR(1),PR(R.list()));
a = g*u;
b = g*v;
sigma = a^2+z*b^2;

u = vector(F,[randint(0,1) for _ in range(k)]);
c = encode(u);
e = vector(F,n); # e = zero vector
for trial in range(t):
j = randint(0,n-1);
e[j] += 1;
y = c+e; c_y = decode(y);

if (c_y == c):
print 'corrected received word == sent word';
else:
print 'alarm: c_y =',c_y,'<>',c,'= c';









