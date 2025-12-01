from .message import decode, encode
from .functions import sample_discrete_gaussian_ZN

def encrypt(plaintext, N, r, q, rng, R_p, moduli, Aq, p, pk):
    x = encode(plaintext, R_p, moduli, Aq)
    a, b = pk

    u_coeff, v_coeff, w_coeff = sample_discrete_gaussian_ZN(N, r, q, rng, num_samples=3)
    u = Aq(u_coeff.tolist())
    v = Aq(v_coeff.tolist())
    w = Aq(w_coeff.tolist())

    c_0 = (b*v) + (p*w) + x
    c_1 = (a*v) + (p*u)
    ciphertext = (c_0, c_1, 0)
    
    return ciphertext

def decrypt(ciphertext, R_p, moduli, sk):
    c0, c1, c2 = ciphertext
    
    t = c0 - (sk*c1) - (sk*sk*c2)
    
    return decode(t, R_p, moduli)
