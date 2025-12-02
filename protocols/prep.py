import numpy as np

from preprocessing.encryption import decrypt, encrypt
from preprocessing.functions import sample_Fpk

# =============================================================================
# 1. FUNCIONES AUXILIARES DE OPERACIÓN HOMOMÓRFICA
# =============================================================================

def add_ciphertexts(ct1, ct2):
    """Suma dos cifrados BGV/BFV (c0, c1, c2) componente a componente."""
    # Normalizar a tuplas de 3 elementos
    c = (list(ct1) + [0]*3)[:3]
    d = (list(ct2) + [0]*3)[:3]
    return (c[0] + d[0], c[1] + d[1], c[2] + d[2])

def mult_ciphertexts(ct1, ct2):
    """
    Multiplicación tensorial de cifrados.
    Entrada: ct1=(c0, c1, 0), ct2=(d0, d1, 0)
    Salida: (c0*d0, c1*d0 + c0*d1, -c1*d1)
    """
    # Asumimos que los cifrados de entrada tienen c2=0 (frescos)
    c0, c1 = ct1[0], ct1[1]
    d0, d1 = ct2[0], ct2[1]
    
    # Fórmulas del paper para el producto tensorial
    new_c0 = c0 * d0
    new_c1 = (c1 * d0) + (c0 * d1)
    new_c2 = -(c1 * d1) 
    
    return (new_c0, new_c1, new_c2)

def sum_cipher_vec(cipher_list, Aq):
    """Suma una lista de cifrados."""
    res = (Aq(0), Aq(0), Aq(0))
    for ct in cipher_list:
        res = add_ciphertexts(res, ct)
    return res

def sum_plain_vec(vec1, vec2):
    """Suma vectores de elementos de cuerpo finito (listas)."""
    return [v1 + v2 for v1, v2 in zip(vec1, vec2)]

# =============================================================================
# 2. SUB-PROTOCOLOS (Reshare, PBracket, PAngle)
# =============================================================================

def reshare(e_m, s, N, k, r, q, rng, R_p, slots_moduli, Aq, p, pk, n, sk, enc=False):
    """
    Protocolo Reshare Corregido.
    """
    # 1. Cada jugador muestrea f_i y lo cifra
    f_list = []      
    e_fi_list = []   
    
    for _ in range(n):
        # f_i es un vector de s elementos
        fi = sample_Fpk(p, k, s) 
        f_list.append(fi)
        
        # Cifrar f_i.
        e_fi = encrypt(fi, N, r, q, rng, R_p, slots_moduli, Aq, p, pk, encode_text=True)
        e_fi_list.append(e_fi)
    
    # 2. Sumar cifrados: e_f = sum(e_fi)
    e_f = sum_cipher_vec(e_fi_list, Aq)
    
    # 3. e_mf = e_m + e_f
    e_mf = add_ciphertexts(e_m, e_f)
    
    # 4. Decriptación Distribuida
    # --- CORRECCIÓN AQUÍ: Pasamos R_p y slots_moduli ---
    mf_plain = decrypt(e_mf, R_p, slots_moduli, sk)
    
    # Nota: Como decrypt ya devuelve la lista de valores (slots), 
    # mf_plain ya es el vector [v1, v2, ...]

    # 5. Calcular shares m_i
    m_i_list = []
    
    # P1: m1 = (m+f) - f1
    m1 = []
    # Aseguramos que mf_plain y f_list[0] sean listas del mismo tamaño
    for val_mf, val_f1 in zip(mf_plain, f_list[0]):
        # Operación en F_pk (o enteros si k=1)
        m1.append(val_mf - val_f1) 
    m_i_list.append(m1)
    
    # Pi: mi = -fi
    for i in range(1, n):
        mi = [-x for x in f_list[i]]
        m_i_list.append(mi)
        
    # 6. NewCiphertext (si enc=True)
    e_m_prime = None
    if enc:
        # e'_m = Enc(m+f) - e_f
        e_mf_fresh = encrypt(mf_plain, N, r, q, rng, R_p, slots_moduli, Aq, p, pk, encode_text=True)
        
        # Restar es sumar el negativo
        e_f_neg = tuple(-x for x in e_f)
        e_m_prime = add_ciphertexts(e_mf_fresh, e_f_neg)
        
    return m_i_list, e_m_prime


def PBracket(shares_v, e_v, e_beta_list, s, N, k, r, q, rng, R_p, slots_moduli, Aq, p, pk, n, sk):
    """
    Protocolo PBracket (Fig 5). Genera la representación [v].
    Entrada: shares v_i (privados), cifrado público e_v, y lista de cifrados de claves MAC e_beta.
    Salida: shares gamma_i (MACs).
    """
    # 1. Calcular e_gamma_i = e_beta_i * e_v para cada jugador i
    # e_beta_list debe contener los cifrados de las claves beta de cada jugador.
    
    gamma_shares_matrix = [] # Matriz n x n: fila j tiene shares de gamma_j
    
    for i in range(n):
        e_beta_i = e_beta_list[i]
        
        # Multiplicación homomórfica
        e_gamma_i_ct = mult_ciphertexts(e_beta_i, e_v)
        
        # 2. Reshare(e_gamma_i) -> cada jugador j recibe un share gamma_i^j
        # Nota: No necesitamos NewCiphertext (enc=False)
        shares_gamma_i, _ = reshare(e_gamma_i_ct, s, N, k, r, q, rng, R_p, slots_moduli, Aq, p, pk, n, sk, enc=False)
        
        gamma_shares_matrix.append(shares_gamma_i)
    
    # Reorganizar para devolver lo que cada jugador almacena.
    # El jugador J necesita tener: (v_j, beta_j, {gamma_1^j, ..., gamma_n^j})
    # Aquí devolvemos la matriz de gammas para que el orquestador la distribuya.
    return gamma_shares_matrix


def PAngle(shares_v, e_v, e_alpha, s, N, k, r, q, rng, R_p, slots_moduli, Aq, p, pk, n, sk):
    """
    Protocolo PAngle (Fig 6). Genera la representación <v>.
    Entrada: shares v_i, cifrado público e_v, cifrado clave global e_alpha.
    """
    # 1. e_v_alpha = e_v * e_alpha
    e_valpha = mult_ciphertexts(e_v, e_alpha)
    
    # 2. Reshare(e_v_alpha) -> shares gamma_i
    shares_gamma, _ = reshare(e_valpha, s, N, k, r, q, rng, R_p, slots_moduli, Aq, p, pk, n, sk, enc=False)
    
    return shares_gamma # Cada elemento i es el share del jugador i del MAC global

# =============================================================================
# 3. PROTOCOLO PRINCIPAL (Initialize, Pair, Triple)
# =============================================================================

class PreprocessingProtocol:
    def __init__(self, context):
        self.ctx = context # Diccionario con N, r, q, etc.
        # Desempaquetar claves comunes
        self.N, self.r, self.q = context['N'], context['r'], context['q']
        self.rng, self.Aq = context['rng'], context['Aq']
        self.pk, self.sk = context['pk'], context['sk']
        self.R_p, self.slots_moduli = context['R_p'], context['slots_moduli']
        self.p, self.n, self.s = context['p'], context['n'], context['s']
        self.k = context['k']
        
        # Estado interno
        self.e_alpha = None
        self.e_beta_list = [] # Cifrados de las claves beta personales

    def run_initialize(self):
        """Fase Initialize (Fig 7)"""
        print("--- Running Initialize ---")
        
        # 1-3. Generar claves alpha y beta
        alpha_i_list = [sample_Fpk(self.p, self.k, self.n) for _ in range(self.n)] # alpha_i es escalar o vector pequeño? Asumimos escalar Fpk repetido s veces
        beta_i_list = [sample_Fpk(self.p, self.k, self.n) for _ in range(self.n)]
        
        # En el paper alpha y beta son escalares en F_pk, pero operamos SIMD (s slots).
        # Convertimos a vectores de tamaño s (diagonales)
        
        e_alpha_i_list = []
        self.e_beta_list = [] # Guardar para uso futuro
        
        for i in range(self.n):
            # Diag(alpha_i)
            msg_alpha = [alpha_i_list[i]] * self.s
            msg_beta = [beta_i_list[i]] * self.s
            
            # 4. Encrypt y Broadcast
            e_a = encrypt(msg_alpha, self.N, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, encode_text=True)
            e_b = encrypt(msg_beta, self.N, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, encode_text=True)
            
            e_alpha_i_list.append(e_a)
            self.e_beta_list.append(e_b)
            
            # 5. ZKPoPK (Omitido por brevedad, llamaríamos a zkpopk_protocol aquí)
            
        # 6. Calcular e_alpha global = sum(e_alpha_i)
        self.e_alpha = sum_cipher_vec(e_alpha_i_list, self.Aq)
        
        # Generar [alpha] usando PBracket sobre e_alpha
        # Nota: PBracket espera shares del valor que está cifrado.
        # Los shares de alpha son alpha_i_list (expandidos a vectores)
        shares_alpha_vec = [[alpha_i_list[i]]*self.s for i in range(self.n)]
        
        bracket_alpha = PBracket(shares_alpha_vec, self.e_alpha, self.e_beta_list, 
                                 self.s, self.N, self.k, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, self.n, self.sk)
        
        print("Initialize completado. Claves generadas.")
        return bracket_alpha

    def run_pair(self):
        """Fase Pair (Fig 7). Genera par [r], <r>."""
        print("--- Running Pair ---")
        
        # 1. Generar r_i
        r_i_list = [sample_Fpk(self.p, self.k, self.s) for _ in range(self.n)]
        
        # 2. Encrypt y Broadcast e_r_i
        e_ri_list = []
        for ri in r_i_list:
            e_ri = encrypt(ri, self.N, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, encode_text=True)
            e_ri_list.append(e_ri)
            
        # 3. ZKPoPK (Omitido)
        
        # e_r = sum(e_ri)
        e_r = sum_cipher_vec(e_ri_list, self.Aq)
        
        # 4. Generar [r] y <r>
        bracket_r = PBracket(r_i_list, e_r, self.e_beta_list, 
                             self.s, self.N, self.k, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, self.n, self.sk)
                             
        angle_r = PAngle(r_i_list, e_r, self.e_alpha, 
                         self.s, self.N, self.k, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, self.n, self.sk)
        
        return r_i_list, bracket_r, angle_r

    def run_triple(self):
        """Fase Triple (Fig 7). Genera <a>, <b>, <c>."""
        print("--- Running Triple ---")
        
        # 1. Generar a_i, b_i
        a_i_list = [sample_Fpk(self.p, self.k, self.s) for _ in range(self.n)]
        b_i_list = [sample_Fpk(self.p, self.k, self.s) for _ in range(self.n)]
        
        # 2. Encrypt y Broadcast
        e_ai_list = []
        e_bi_list = []
        for i in range(self.n):
            e_a = encrypt(a_i_list[i], self.N, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, encode_text=True)
            e_b = encrypt(b_i_list[i], self.N, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, encode_text=True)
            e_ai_list.append(e_a)
            e_bi_list.append(e_b)
            
        # 4. Sumar cifrados e_a, e_b
        e_a = sum_cipher_vec(e_ai_list, self.Aq)
        e_b = sum_cipher_vec(e_bi_list, self.Aq)
        
        # 5. Generar <a>, <b>
        angle_a = PAngle(a_i_list, e_a, self.e_alpha, 
                         self.s, self.N, self.k, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, self.n, self.sk)
        angle_b = PAngle(b_i_list, e_b, self.e_alpha, 
                         self.s, self.N, self.k, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, self.n, self.sk)
        
        # 6. e_c = e_a * e_b
        e_c = mult_ciphertexts(e_a, e_b)
        
        # 7. Reshare(e_c) -> shares c_i, e_c_prime (fresh)
        c_i_list, e_c_prime = reshare(e_c, self.s, self.N, self.k, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, self.n, self.sk, enc=True)
        
        # 8. Generar <c> usando e_c_prime (que es un cifrado fresco de c)
        angle_c = PAngle(c_i_list, e_c_prime, self.e_alpha, 
                         self.s, self.N, self.k, self.r, self.q, self.rng, self.R_p, self.slots_moduli, self.Aq, self.p, self.pk, self.n, self.sk)
                         
        return a_i_list, b_i_list, c_i_list, angle_a, angle_b, angle_c
