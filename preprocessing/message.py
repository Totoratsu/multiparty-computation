import gmpy2
import multiprocessing
from sage.all import GF, PolynomialRing, ZZ, CRT


# -----------------------------
# Preparar ranuras (usa tu prepare_slots si prefieres)
# -----------------------------
def prepare_slots(p, coeffs):
    """
    Calc factores irreducibles de Phi_m mod p y devuelve R_p y moduli.
    Asume coeffs como lista/tupla de coeficientes para Phi_m (orden estándar).
    """
    R_p = PolynomialRing(GF(p), 'x')
    Phi_p = R_p(coeffs)
    factors = Phi_p.factor()
    moduli = [f for f, exp in factors]
    return R_p, moduli

def poly_crt(residues, moduli, R_p):
    """
    Reconstruye un polinomio en R_p a partir de residuos (lista de polinomios)
    y los módulos (moduli), usando el CRT polinomial clásico.
    """
    M = R_p(1)
    for f in moduli:
        M *= f
    total = R_p(0)
    for r, f in zip(residues, moduli):
        Mi = M // f                    # cociente exacto M/f en R_p
        # Nota: La implementación nativa de Sage CRT(...) es mucho más rápida y recomendada
        inv = Mi.inverse_mod(f)        # inverso polinomial de Mi mod f
        total += r * Mi * inv
    return total % M

# -----------------------------
# Encode (robusto)
# -----------------------------
def center_worker(coeffs_chunk, p):
    """
    Función ejecutada por cada núcleo para centrar una porción de coeficientes.
    """
    centered_coeffs = []
    limit = p // 2
    p_int = int(p)
    
    for c in coeffs_chunk:
        val = int(c)
        if val > limit:
            centered_coeffs.append(val - p_int)
        else:
            centered_coeffs.append(val)
    return centered_coeffs

# -----------------------------
# Encode (PARALELIZADO y CORREGIDO)
# -----------------------------
def encode(messages, R_p, moduli, Aq, use_parallel=True):
    p = R_p.characteristic()
    if len(messages) > len(moduli):
        raise ValueError("Demasiados mensajes para el número de slots")

    # 1. Preparación de residuos (Sequential)
    poly_msgs = [R_p(m) for m in messages]
    while len(poly_msgs) < len(moduli):
        poly_msgs.append(R_p(0))

    # 2. CRT Polinomial (Usando la implementación nativa y optimizada de Sage)
    # ESTE ES EL PASO MÁS LENTO Y NO DEBE SER PARALELIZADO MANUALMENTE.
    plaintext_poly_mod_p = CRT(poly_msgs, moduli)

    # 3. Centrado exacto en (-p/2, ..., p/2] (PARALELIZABLE)
    coeffs_list = plaintext_poly_mod_p.list()

    if use_parallel:
        num_cores = multiprocessing.cpu_count()
        chunk_size = int(gmpy2.ceil(len(coeffs_list) / num_cores))
        
        # Dividir la lista de coeficientes en trozos
        chunks = [coeffs_list[i:i + chunk_size] 
                  for i in range(0, len(coeffs_list), chunk_size)]

        with multiprocessing.Pool(processes=num_cores) as pool:
            # pool.starmap aplica center_worker a cada trozo (chunk)
            # Pasamos p a cada worker
            results = pool.starmap(center_worker, [(chunk, p) for chunk in chunks])

        # Aplanar la lista de resultados de los chunks
        centered_coeffs = [item for sublist in results for item in sublist]
    else:
        # Versión secuencial de centrado (para referencia)
        centered_coeffs = []
        limit = p // 2
        p_int = int(p)
        for c in coeffs_list:
            val = int(c)
            if val > limit:
                centered_coeffs.append(val - p_int)
            else:
                centered_coeffs.append(val)

    # 4. Reconstrucción sobre ZZ y pasar a Aq (Sequential)
    R_ZZ = PolynomialRing(ZZ, 'x')
    poly_integers = R_ZZ(centered_coeffs)
    
    return Aq(poly_integers)

# -----------------------------
# FUNCIÓN AUXILIAR DE CENTRADO MODULAR
# -----------------------------
def center_mod_q(poly_Aq):
    """
    Toma un polinomio con coeficientes mod q y lo centra en Z
    en el rango [-q/2, q/2].
    """
    # Obtenemos el módulo q del anillo Aq
    q = poly_Aq.parent().characteristic()
    limit = q // 2
    
    # 1. Obtenemos los coeficientes como lista de enteros [0, q-1]
    try:
        coeffs_raw = poly_Aq.lift().list()
    except Exception:
        coeffs_raw = poly_Aq.list()

    # 2. Aplicamos el centrado
    centered_coeffs = []
    for c in coeffs_raw:
        val = int(c)
        if val > limit:
            centered_coeffs.append(val - int(q))
        else:
            centered_coeffs.append(val)
            
    # 3. Reconstruimos en el anillo Z[x]
    R_ZZ = PolynomialRing(ZZ, 'x')
    return R_ZZ(centered_coeffs)

# -----------------------------
# Decode (CORREGIDO PARA DOBLE CENTRADO)
# -----------------------------
def decode(element_Aq, R_p, moduli):
    """
    Recupera el vector de mensajes. Aplica:
    1. Centrado mod q (para obtener el polinomio más pequeño en Z[x]).
    2. Reducción mod p.
    3. Descenso final mod p (para obtener el mensaje pequeño).
    """
    p = R_p.characteristic() 
    
    # 1. Centrar el elemento en el anillo de cifrado (mod q)
    # Esto es crucial si el ruido es grande: obtenemos el polinomio más pequeño en Z[x].
    poly_zz_centered = center_mod_q(element_Aq)
        
    # 2. Reducción módulo p
    # Convertir coeficientes de ZZ a R_p (mod p).
    clean_coeffs = [int(c) for c in poly_zz_centered.list()]
    poly_mod_p = R_p(clean_coeffs)

    # 3. Recuperación de slots y Descentrado final
    limit = p // 2           
    decoded = []
    
    for f in moduli:
        res = poly_mod_p % f
        m_standard = int(res.constant_coefficient())
        
        # APLICAR DESCENTRADO: De [0, p-1] a [-p/2, p/2]
        if m_standard > limit:
            # Si el residuo es grande, es un número negativo.
            m_centered = m_standard - int(p)
        else:
            m_centered = m_standard

        decoded.append(m_centered) 
            
    return decoded
