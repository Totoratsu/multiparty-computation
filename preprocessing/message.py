from sage.all import GF, PolynomialRing, ZZ


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
def encode(messages, R_p, moduli, Aq):
    p = R_p.characteristic()
    if len(messages) > len(moduli):
        raise ValueError("Demasiados mensajes para el número de slots")

    # Convertir mensajes a residuos en R_p
    poly_msgs = [R_p(m) for m in messages]
    while len(poly_msgs) < len(moduli):
        poly_msgs.append(R_p(0))

    # CRT polinomial robusto
    plaintext_poly_mod_p = poly_crt(poly_msgs, moduli, R_p)

    # Centrado exacto en (-p/2, ..., p/2]
    coeffs_list = plaintext_poly_mod_p.list()
    centered_coeffs = []
    limit = p // 2
    for c in coeffs_list:
        val = int(c)
        if val > limit:
            centered_coeffs.append(val - int(p))
        else:
            centered_coeffs.append(val)

    # Reconstrucción sobre ZZ y pasar a Aq
    R_ZZ = PolynomialRing(ZZ, 'x')
    poly_integers = R_ZZ(centered_coeffs)
    try:
        return Aq(poly_integers)
    except Exception:
        # si Aq es un anillo cociente, usamos el constructor de la forma Aq(poly)
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
    # Si poly_Aq es un elemento de un cociente (ej. Aq), usamos lift()
    try:
        coeffs_raw = poly_Aq.lift().list()
    except Exception:
        # Si ya es un polinomio simple (ej. en Rq), usamos list()
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
