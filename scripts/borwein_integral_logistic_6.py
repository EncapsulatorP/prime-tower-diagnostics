import mpmath as mp
import numpy as np


mp.mp.dps = 60


# -----------------------------
# Basic sinc and Borwein-like integral
# -----------------------------

def sinc(z):
    if abs(z) < mp.mpf("1e-40"):
        return mp.mpf(1)
    return mp.sin(z) / z


def borwein_like_integral(scales, cutoff=250, chunks=500):
    """
    Computes I = ∫_0^∞ Π sinc(x/a_j) dx approximately.

    This is a diagnostic integral, not a theorem-grade quadrature.
    Increase cutoff/chunks for stability checks.
    """
    scales = [mp.mpf(s) for s in scales]

    def f(x):
        y = mp.mpf(1)
        for a in scales:
            y *= sinc(x / a)
        return y

    total = mp.mpf(0)
    step = mp.mpf(cutoff) / chunks

    for k in range(chunks):
        a = k * step
        b = (k + 1) * step
        total += mp.quad(f, [a, b])

    return total


# -----------------------------
# Logistic Lyapunov exponent
# -----------------------------

def logistic_lyapunov(r, x0=mp.mpf("0.3141592653"), burn=5000, N=20000):
    r = mp.mpf(r)
    x = mp.mpf(x0)

    for _ in range(burn):
        x = r * x * (1 - x)

    s = mp.mpf(0)

    for _ in range(N):
        d = abs(r * (1 - 2 * x))
        if d <= 0:
            d = mp.mpf("1e-50")
        s += mp.log(d)
        x = r * x * (1 - x)

    return s / N


# -----------------------------
# c10 positional signature
# -----------------------------

def c10_digit_signature(value, digits=80):
    """
    Converts a high-precision number into a base-10 digit frequency
    and positional weighted code.

    This is intentionally lossy: it is your compressed c10 layer.
    """
    s = mp.nstr(abs(value), n=digits + 10)
    s = s.replace(".", "").replace("-", "")

    digs = [int(ch) for ch in s if ch.isdigit()]
    digs = digs[:digits]

    if len(digs) == 0:
        return np.zeros(20)

    hist = np.bincount(digs, minlength=10).astype(float)
    hist = hist / hist.sum()

    # positional moments
    pos = np.arange(1, len(digs) + 1, dtype=float)
    arr = np.array(digs, dtype=float)

    m1 = np.sum(arr / pos)
    m2 = np.sum((arr ** 2) / (pos ** 2))
    alt = np.sum(((-1) ** pos) * arr / pos)

    # digit transition counts mod 10
    trans = np.zeros(10)
    for a, b in zip(digs[:-1], digs[1:]):
        trans[(b - a) % 10] += 1

    if trans.sum() > 0:
        trans /= trans.sum()

    return np.concatenate([hist, trans, [m1, m2, alt]])


def cosine_distance(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    den = np.linalg.norm(a) * np.linalg.norm(b)
    if den == 0:
        return np.nan

    return 1 - float(np.dot(a, b) / den)


# -----------------------------
# Operation-tower construction
# -----------------------------

def tower_scales(level):
    """
    Operation introduction ladder:

    L0: {2,3}
    L1: add 6 = 2*3
    L2: add sqrt(6)
    L3: add 1+sqrt(6)
    L4: add reciprocal critical scale 1/(1+sqrt(6))
    """
    two = mp.mpf(2)
    three = mp.mpf(3)
    six = two * three
    root6 = mp.sqrt(six)
    crit = 1 + root6

    scales = [two, three]

    if level >= 1:
        scales.append(six)
    if level >= 2:
        scales.append(root6)
    if level >= 3:
        scales.append(crit)
    if level >= 4:
        scales.append(1 / crit)

    return scales


def run_test():
    print("=== L4 positional tower / Borwein-like / logistic test ===")
    print()

    pi = mp.pi
    rcrit = 1 + mp.sqrt(6)

    lyap_crit = logistic_lyapunov(rcrit)

    print(f"rcrit = 1 + sqrt(6) = {mp.nstr(rcrit, 30)}")
    print(f"logistic Lyapunov at rcrit ≈ {mp.nstr(lyap_crit, 20)}")
    print()

    rows = []

    for L in range(0, 5):
        scales = tower_scales(L)
        I = borwein_like_integral(scales, cutoff=220, chunks=350)

        delta_pi = abs((2 * I / pi) - 1)

        sig_I = c10_digit_signature(I)
        sig_pi = c10_digit_signature(pi / 2)
        sig_r = c10_digit_signature(rcrit)

        d_I_pi = cosine_distance(sig_I, sig_pi)
        d_I_r = cosine_distance(sig_I, sig_r)

        rows.append((L, scales, I, delta_pi, d_I_pi, d_I_r))

        print(f"L={L}")
        print("  scales =", [mp.nstr(s, 12) for s in scales])
        print(f"  B_L    = {mp.nstr(I, 30)}")
        print(f"  |2B_L/pi - 1| = {mp.nstr(delta_pi, 12)}")
        print(f"  c10 distance B_L vs pi/2       = {d_I_pi:.6f}")
        print(f"  c10 distance B_L vs 1+sqrt(6)  = {d_I_r:.6f}")
        print()

    print("=== Operation introduction jumps ===")
    for i in range(1, len(rows)):
        L0, _, I0, _, _, _ = rows[i - 1]
        L1, _, I1, _, _, _ = rows[i]

        sig0 = c10_digit_signature(I0)
        sig1 = c10_digit_signature(I1)

        jump = cosine_distance(sig0, sig1)

        print(f"L{L0} -> L{L1}: c10 jump = {jump:.6f}")


if __name__ == "__main__":
    run_test()