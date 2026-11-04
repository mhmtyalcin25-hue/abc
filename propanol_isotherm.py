# Saf 1-propanol, 350 C izotermi: PC-SAFT ile hacim/yogunluk - basinc
# Kullanim: python3 propanol_isotherm.py  -> propanol_350C.csv, propanol_350C.png
import numpy as np
import pcsaft_density as ps

T_C = 350.0
P_BAR = [1, 5, 10, 20, 30, 40, 50, 60, 80, 100, 150, 200, 300, 400, 500]

ps.x = np.array([0.0, 1.0])
ps.T = T_C + 273.15
T = ps.T

def molar_density(P):
    """Verilen P (Pa) icin yogunluk, 1/A^3.
    350 C propanol icin kritik ustu: P(rho) monoton, tek kok -> dogrudan bisection
    (6000 noktalik tarama yerine ~60 degerlendirme)."""
    d = ps.s * (1 - 0.12 * np.exp(-3 * ps.e / T))
    rmax = 0.74 / (np.pi / 6 * np.sum(ps.x * ps.m * d**3))
    return ps.brentq(lambda r: ps.Pcalc(r) - P, 1e-8, rmax * 0.999)

rows = []
for Pb in P_BAR:
    r = molar_density(Pb * 1e5)
    rm = r * 1e30 / ps.NA                # mol/m3
    rho = rm * ps.M[1] / 1e3             # kg/m3
    rows.append((Pb, 1e6 / rm, 1e3 / rho, rho, ps.Z(r)))

with open("propanol_350C.csv", "w") as fh:
    fh.write("P_bar,V_cm3_per_mol,V_L_per_kg,rho_kg_per_m3,Z\n")
    for row in rows:
        fh.write("{},{:.2f},{:.4f},{:.3f},{:.4f}\n".format(*row))
print("P(bar)  V(cm3/mol)  rho(kg/m3)  Z")
for Pb, V, _, rho, z in rows:
    print(f"{Pb:5d}  {V:10.1f}  {rho:9.2f}  {z:.4f}")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("matplotlib yok; grafik atlandi")
else:
    Pf = np.linspace(1, 500, 120)
    Vf, Zf = [], []
    for P in Pf:
        r = molar_density(P * 1e5)
        Vf.append(1e6 / (r * 1e30 / ps.NA)); Zf.append(ps.Z(r))
    P_ = [r[0] for r in rows]; V_ = [r[1] for r in rows]; Z_ = [r[4] for r in rows]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    a1.plot(Vf, Pf, color="#2a6fdb"); a1.plot(V_, P_, "o", color="#2a6fdb", ms=4)
    a1.set_xscale("log"); a1.set_xlabel("V (cm³/mol)"); a1.set_ylabel("P (bar)")
    a1.set_title("1-Propanol, 350 °C: P–V (PC-SAFT)"); a1.grid(alpha=.3, which="both")
    a2.plot(Pf, Zf, color="#d9622b"); a2.plot(P_, Z_, "o", color="#d9622b", ms=4)
    a2.axhline(1, color="gray", lw=.8, ls="--")
    a2.set_xlabel("P (bar)"); a2.set_ylabel("Z"); a2.set_title("Sıkıştırılabilirlik faktörü")
    a2.grid(alpha=.3)
    fig.tight_layout(); fig.savefig("propanol_350C.png", dpi=150)
    print("propanol_350C.png yazildi")
