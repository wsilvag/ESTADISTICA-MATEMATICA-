"""
Intervalo de confianza MEDIA
  1) Muestra chica  (n < 30)  -> t de Student
  2) Muestra grande (n >= 30) -> z

"""

# ====================== MATEMATICA BASE ======================
PI = 3.141592653589793
E = 2.718281828459045
LN2 = 0.6931471805599453

def raiz(x): return x ** 0.5
def exp(x): return E ** x

def ln(x):
    k = 0
    while x > 2: x /= 2; k += 1
    while x < 0.5: x *= 2; k -= 1
    s = (x - 1) / (x + 1); t = s; suma = 0.0; i = 1
    while abs(t) > 1e-17:
        suma += t / i; t *= s * s; i += 2
    return 2 * suma + k * LN2

_L = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
      771.32342877765313, -176.61502916214059, 12.507343278686905,
      -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
def lgamma(z):
    z -= 1; x = _L[0]
    for i in range(1, 9): x += _L[i] / (z + i)
    t = z + 7.5
    return 0.5 * ln(2 * PI) + (z + 0.5) * ln(t) - t + ln(x)

def gammp(a, x):             #gamma 
    if x <= 0: return 0.0
    if x < a + 1:
        ap, s = a, 1 / a; d = s
        for _ in range(1000):
            ap += 1; d *= x / ap; s += d
            if abs(d) < abs(s) * 1e-16: break
        return s * exp(-x + a * ln(x) - lgamma(a))
    b = x + 1 - a; c = 1e30; d = 1 / b; h = d
    for i in range(1, 1000):
        an = -i * (i - a); b += 2
        d = an * d + b; d = 1e-30 if abs(d) < 1e-30 else d
        c = b + an / c; c = 1e-30 if abs(c) < 1e-30 else c
        d = 1 / d; de = d * c; h *= de
        if abs(de - 1) < 1e-16: break
    return 1 - exp(-x + a * ln(x) - lgamma(a)) * h

def _betacf(a, b, x):
    qab, qap, qam = a + b, a + 1, a - 1
    c = 1.0; d = 1 - qab * x / qap
    d = 1e-30 if abs(d) < 1e-30 else d
    d = 1 / d; h = d
    for m in range(1, 1000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d; d = 1e-30 if abs(d) < 1e-30 else d
        c = 1 + aa / c; c = 1e-30 if abs(c) < 1e-30 else c
        d = 1 / d; h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d; d = 1e-30 if abs(d) < 1e-30 else d
        c = 1 + aa / c; c = 1e-30 if abs(c) < 1e-30 else c
        d = 1 / d; de = d * c; h *= de
        if abs(de - 1) < 1e-16: break
    return h

def betai(a, b, x):                   # beta
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * ln(x) + b * ln(1 - x))
    if x < (a + 1) / (a + b + 2): return bt * _betacf(a, b, x) / a
    return 1 - bt * _betacf(b, a, 1 - x) / b

# ---------------- distribuciones ----------------
def cdf_z(x):
    p = 0.5 * gammp(0.5, x * x / 2)
    return 0.5 + p if x >= 0 else 0.5 - p

def cdf_t(x, v):
    q = 0.5 * betai(v / 2, 0.5, v / (v + x * x))
    return 1 - q if x >= 0 else q

def inversa(cdf, p, lo, hi):          # biseccion
    while cdf(hi) < p: hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if cdf(mid) < p: lo = mid
        else: hi = mid
    return (lo + hi) / 2

def z_tabla(p):    return inversa(cdf_z, p, -10, 10)
def t_tabla(p, v): return inversa(lambda x: cdf_t(x, v), p, 0, 10)

# ====================== ENTRADA DE DATOS ======================
MODO_TABLA = False

def pedir_num(msg):
    while True:
        try: return float(input(msg).replace(",", "."))
        except ValueError: print("  Valor no valido.")

def pedir_conf():
    c = pedir_num("Nivel de confianza (ej. 95 o 0.95): ")
    return c / 100 if c > 1 else c

def critico(nombre, calculado):
    if not MODO_TABLA: return calculado
    r = input(f"  Valor de tabla para {nombre} [Enter = {calculado:.4f}]: ").strip()
    return float(r.replace(",", ".")) if r else calculado

def decision(lo, hi, valor):
    dentro = lo <= valor <= hi
    print(f"  {valor} {'SI' if dentro else 'NO'} esta en el IC -> " +
          ("NO se rechaza H0" if dentro else "se RECHAZA H0 (la media es distinta)"))

# ====================== PROBLEMAS ======================
def media_chica():
    print("\n[MEDIA - MUESTRA CHICA (t de Student)]")
    print("1) Tengo los datos   2) Tengo n, media y desviacion s")
    if input("Opcion: ").strip() == "1":
        datos = [float(x.replace(",", ".")) for x in input("Datos separados por espacio: ").split()]
        n = len(datos); xb = sum(datos) / n
        s = raiz(sum((d - xb) ** 2 for d in datos) / (n - 1))
    else:
        n = int(pedir_num("n: ")); xb = pedir_num("Media: "); s = pedir_num("Desviacion s: ")
    c = pedir_conf(); gl = n - 1
    t = critico(f"t({1-(1-c)/2:.4f}; gl={gl})", t_tabla(1 - (1 - c) / 2, gl))
    e = t * s / raiz(n)
    print(f"\n  x̄={xb:.4f}  s={s:.4f}  n={n}  gl={gl}  t={t:.4f}  error={e:.4f}")
    print(f"  IC {c*100:.1f}% para mu: ({xb-e:.4f}, {xb+e:.4f})")
    m0 = input("  Valor a contrastar mu0 (Enter para omitir): ").strip()
    if m0: decision(xb - e, xb + e, float(m0))

def media_grande():
    print("\n[MEDIA - MUESTRA GRANDE (z)]")
    n = int(pedir_num("n: ")); xb = pedir_num("Media muestral: ")
    s = pedir_num("Desviacion (sigma poblacional o s): "); c = pedir_conf()
    z = critico(f"z({1-(1-c)/2:.4f})", z_tabla(1 - (1 - c) / 2))
    e = z * s / raiz(n)
    print(f"\n  z={z:.4f}  error={e:.4f}")
    print(f"  IC {c*100:.1f}% para mu: ({xb-e:.4f}, {xb+e:.4f})")
    m0 = input("  Valor a contrastar mu0 (Enter para omitir): ").strip()
    if m0: decision(xb - e, xb + e, float(m0))

# ====================== MENU ======================
def main():
    global MODO_TABLA
    print("INTERVALOS DE CONFIANZA PARA LA MEDIA - Estadistica Matematica I")
    MODO_TABLA = input("¿Quieres escribir tu los valores de tabla? (s/n): ").strip().lower() == "s"
    while True:
        print("\n1) Media muestra chica (t)   2) Media muestra grande (z)   0) Salir")
        op = input("Opcion: ").strip()
        if op == "0": break
        elif op == "1": media_chica()
        elif op == "2": media_grande()
        else: print("Opcion no valida.")

if __name__ == "__main__":
    main()