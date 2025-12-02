from preprocessing.encryption import encrypt
from preprocessing.functions import sample_Fpk

def Pi_prep_init(s, N, r, q, rng, R_p, slots_moduli, Aq, p, pk, k, n):
    alpha_list = sample_Fpk(p, k, n)
    alpha = sum(alpha_list)
    e_a_list = []

    # MAC
    beta_list = sample_Fpk(p, k, n)
    e_b_list = []

    for ai, bi in zip(alpha_list, beta_list):
        e_a_list.append(encrypt([ai]*s, N, r, q, rng, R_p, slots_moduli, Aq, p, pk))
        e_b_list.append(encrypt([bi]*s, N, r, q, rng, R_p, slots_moduli, Aq, p, pk))
    
    return alpha, alpha_list, beta_list, e_a_list, e_b_list
