from functools import partial
from sage.all import Zmod, PolynomialRing, ZZ, GF
from multiprocessing import Pool, cpu_count


################################################################
#                        Aq Generation                         #
################################################################

def parameters_worker(m, p, k, s):
    from sage.all import euler_phi, gcd, Zmod, cyclotomic_polynomial, PolynomialRing

    # 1. Filtro rápido de GCD
    if gcd(m, p) != 1:
        return None

    # 2. Filtro rápido de phi(m)
    # Calculamos phi primero porque es barato y descarta muchos candidatos
    phi_m = int(euler_phi(m))
    if phi_m < 14300:
        return None

    # 3. EL GRAN CAMBIO: Calcular el orden multiplicativo
    # Esto reemplaza a la factorización polinómica. Es O(log m) vs O(poly_degree).
    try:
        Zn = Zmod(m)
        d = int(Zn(p).multiplicative_order())
    except ArithmeticError:
        return None

    # Esto garantiza que Phi_m mod p se rompa en factores lineales (p ≡ 1 (mod m)).
    if d != 1:
        return None

    # 4. Verificaciones aritméticas (sin polinomios aún)
    # Condición A: El grado d debe contener a la extensión k (d debe ser múltiplo de k)
    if d % k != 0:
        return None
        
    # Condición B: Debe haber suficientes slots
    # El número de factores es phi(m) / d
    num_slots = phi_m // d
    if num_slots < s:
        return None

    # --- SI LLEGAMOS AQUÍ, EL M ES CORRECTO ---
    # Recién ahora gastamos recursos en construir el objeto pesado para devolverlo.
    
    Phi_ZZ = cyclotomic_polynomial(m)
    R = PolynomialRing(Zmod(p), 'x')
    Phi = R(Phi_ZZ)
    
    # Ni siquiera necesitamos factorizar 'Phi' realmente para devolver la info,
    # porque ya sabemos que hay 'num_slots' factores de grado 'd'.
    # Pero si tu código posterior necesita los factores explícitos serializados:
    fac = Phi.factor()
    
    # Reconstruimos outputs para compatibilidad con tu código original
    coeffs = [int(c) for c in Phi.list()]
    fac_serial = [(str(f), int(e)) for (f, e) in fac]
    degrees = [d] * num_slots # Todos tienen el mismo grado d

    return (m, coeffs, fac_serial, degrees, phi_m)

def calculate_m(p, k, s, min_m, max_m):
    # Preconfiguramos el worker para que reciba sólo m
    worker_partial = partial(parameters_worker, p=p, k=k, s=s)

    found = None
    pool = Pool(processes=cpu_count())
    try:
        # Recorremos los valores de m en paralelo.
        # Mat. hablando: probamos distintos índices ciclotómicos
        # buscando uno cuya reducción mod p produzca los grados deseados.
        for res in pool.imap_unordered(worker_partial, range(min_m, max_m + 1), chunksize=16):
            if res is not None:
                found = res
                break
    finally:
        pool.terminate()
        pool.join()

    if not found:
        print("No encontrado en rango; aumenta max_m o relaja condiciones.")
        return None

    m, coeffs, fac_serial, degrees, phi_m = found

    print("m =", m, "phi(m) =", phi_m, "degrees:", degrees)
    print("coeff:", coeffs)
    
    return m, phi_m, coeffs

def generate_Aq(q, coeffs):
# 1. Definir el anillo base Z_q (enteros modulo q)
    # q suele ser una potencia de 2 (ej. 2^64 o 2^128)
    R_q = Zmod(q)
    
    # 2. Definir el anillo de polinomios sobre Z_q
    PolyRing_q = PolynomialRing(R_q, 'x')
    
    # 3. Reconstruir el polinomio ciclotómico Phi_m usando los coeficientes
    # Los coeficientes vienen de Phi_m calculado en Z (o modulo p), 
    # pero aquí los interpretamos como elementos de Z_q.
    Phi = PolyRing_q(coeffs)

    # 4. Construir el anillo cociente A_q = Z_q[x] / (Phi_m)
    # Este es el espacio donde viven los textos cifrados.
    Aq = PolyRing_q.quotient(Phi, 'a')
    
    print(f"Construido A_q = {Aq}")
    
    return Aq