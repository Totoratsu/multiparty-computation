import numpy as np
import hashlib
import random

from preprocessing import encrypt

# -------------------------------
# Util: construir M_e con NumPy
# -------------------------------
def build_Me_numpy(e_bits):
    """
    e_bits: iterable de 0/1 (length = sec)
    Returns: Me as numpy.ndarray shape (V, sec) dtype=int
    """
    sec = len(e_bits)
    V = 2 * sec - 1
    Me = np.zeros((V, sec), dtype=np.int8)
    for i in range(V):
        for k in range(sec):
            j = i - k
            if 0 <= j < sec and e_bits[j]:
                Me[i, k] = 1
    return Me

# -------------------------------
# Fiat-Shamir: derive challenge bits from a_list + pk
# -------------------------------
def fiat_shamir_challenge_from_serialized(a_serial_list, pk_bytes, sec, hash_fn='sha256'):
    h = hashlib.new(hash_fn)
    for a in a_serial_list:
        if isinstance(a, bytes):
            h.update(a)
        else:
            h.update(str(a).encode())
    h.update(pk_bytes or b'')
    digest = h.digest()
    bits = []
    while len(bits) < sec:
        digest = hashlib.new(hash_fn, digest).digest()
        for byte in digest:
            for b in range(8):
                bits.append((byte >> b) & 1)
                if len(bits) >= sec:
                    break
            if len(bits) >= sec:
                break
    return bits[:sec]

# -------------------------------
# Default comparator (componente a componente)
# -------------------------------
def default_compare_ciphertexts(c1, c2):
    """
    Compara ciphertexts componente-a-componente.
    Soporta que ciphertext sea lista/tuple/np.array de componentes.
    """
    try:
        for a, b in zip(c1, c2):
            if a != b:
                return False
        return True
    except Exception:
        # fallback: compare serialized forms
        return str(c1) == str(c2)

# -------------------------------
# Prover: Fiat-Shamir non-interactive (usa numpy)
# -------------------------------
def zk_pok_prover_fiat_shamir_numpy(pk,
                                    ciphertexts_c,      # list length = sec of c_k
                                    witness_xr,         # list length = sec of (x_k, r_k)
                                    sample_random_plaintext, # callable -> random plaintext (in message-space)
                                    sec,
                                    N, r, q, rng, R_p, slots_moduli, Aq, p,
                                    serialize_ciphertext=lambda C: str(C).encode()):
    """
    Devuelve un 'transcript' dict con keys: a_list (np.array dtype=object),
    a_serial (list bytes), e_bits (list int), z_list (np.array dtype=object), T_list (np.array dtype=object)
    """

    V = 2 * sec - 1

    # 1) Prover samples y_i and s_i, computes a_i = Encpk(y_i, s_i)
    y_list = [None] * V
    s_list = [None] * V
    a_list = [None] * V

    for i in range(V):
        mi = sample_random_plaintext()               # user-provided
        ai, yi, si = encrypt(mi, N, r, q, rng, R_p, slots_moduli, Aq, p, pk, verbose=True)

        y_list[i] = yi # encoded
        s_list[i] = si # random vector
        a_list[i] = ai # ciphertext

    # Serializaciones para Fiat-Shamir
    a_serial = [serialize_ciphertext(a) for a in a_list]
    pk_bytes = serialize_ciphertext(pk) if not isinstance(pk, bytes) else pk

    # 2) derive e via Fiat-Shamir
    e_bits = fiat_shamir_challenge_from_serialized(a_serial, pk_bytes, sec)

    # 3) Build Me matrix (numpy)
    Me = build_Me_numpy(e_bits)  # shape (V, sec), dtype=int8

    # 4) Build x_vec and r_vec from witness_xr
    x_vec = np.array([xr[0] for xr in witness_xr], dtype=object)   # length sec
    r_vec = np.array([xr[1] for xr in witness_xr], dtype=object)

    # convert y_list/s_list to numpy object arrays
    y_arr = np.array(y_list, dtype=object)   # length V
    s_arr = np.array(s_list, dtype=object)

    # 5) compute z_list and T_list using numpy for indexing (but object-sums)
    z_list = np.empty(V, dtype=object)
    T_list = np.empty(V, dtype=object)

    for i in range(V):
        # start from y_i and s_i
        zi = y_arr[i]
        Ti = s_arr[i]
        # add Me[i,k] * x_k or r_k for those k with Me[i,k]==1
        ks = np.nonzero(Me[i, :])[0]
        for k in ks:
            zi = zi + x_vec[k]
            Ti = Ti + r_vec[k]
        z_list[i] = zi
        T_list[i] = Ti

    # 6) Return transcript
    return {
        'a_list': np.array(a_list, dtype=object),
        'a_serial': a_serial,
        'e_bits': e_bits,
        'z_list': z_list,   # numpy array dtype object
        'T_list': T_list,   # numpy array dtype object
        'Me': Me
    }

# -------------------------------
# Verifier: checks given transcript (uses numpy)
# -------------------------------
def zk_pok_verify_fiat_shamir(pk,
                             ciphertexts_c,    # list of original c_k (len sec)
                             transcript,       # dict returned by prover
                             N, r, q, rng, R_p, slots_moduli, Aq, p,
                             compare_ciphertexts_fn=default_compare_ciphertexts,
                             serialize_ciphertext=lambda C: str(C).encode()):
    """
    Verifier: reconstruye e y chequea d_i = Encpk([z_i], T_i) == a_i + sum_k Me[i,k]*c_k
    Observación: encrypt espera el mensaje como lista -> usamos encrypt([z], ...)
    """
    sec = len(ciphertexts_c)
    V = 2 * sec - 1

    a_list = list(transcript['a_list'])
    a_serial = transcript['a_serial']
    e_bits = transcript['e_bits']
    z_list = list(transcript['z_list'])
    T_list = list(transcript['T_list'])
    Me = transcript.get('Me', build_Me_numpy(e_bits))

    # Recompute challenge para asegurar integridad del transcript
    pk_bytes = serialize_ciphertext(pk) if not isinstance(pk, bytes) else pk
    e_check = fiat_shamir_challenge_from_serialized(a_serial, pk_bytes, sec)
    if e_check != e_bits:
        print("[VERIFY] Fiat-Shamir challenge mismatch")
        return False

    # compute d_i = Encpk([z_i], T_i)   <-- observe la lista alrededor de z_i
    d_list = [
        encrypt(z_list[i], N, T_list[i], q, rng, R_p, slots_moduli, Aq, p, pk, custom_r=True, encode_text=False)
        for i in range(V)
    ]

    # compute rhs_i = a_i + sum_{k:Me[i,k]==1} c_k
    rhs_list = [None] * V
    for i in range(V):
        rhs = a_list[i]
        ks = np.nonzero(Me[i, :])[0]
        for k in ks:
            rhs = rhs + ciphertexts_c[k]
        rhs_list[i] = rhs

    # compare d_list y rhs_list componente a componente
    for i in range(V):
        if not compare_ciphertexts_fn(d_list[i], rhs_list[i]):
            print(f"[VERIFY FAIL] Mismatch en i={i}")
            # debug helpful info
            print("  e_bits (sample):", e_bits[:min(10, len(e_bits))])
            print("  tipo d_list[i]:", type(d_list[i]), " tipo rhs_list[i]:", type(rhs_list[i]))
            # opcional: mostrar repr parciales si es útil (cuidado con tamaño)
            # print("  repr(d_list[i])[:200]:", repr(d_list[i])[:200])
            # print("  repr(rhs_list[i])[:200]:", repr(rhs_list[i])[:200])
            return False

    # Si todo coincide
    return True
