# PC-SAFT yoğunluk hesabı: su + 1-propanol

`pcsaft_density.py`, su ve 1-propanol karışımının yoğunluğunu PC-SAFT durum denklemiyle hesaplar
(hard-chain + dispersiyon + 2B asosiyasyon). Yalnızca `numpy` gerekir.

```bash
python3 pcsaft_density.py
```

## Varsayılan durum
- x_su = x_propanol = 0,5 (mol bazında)
- T = 360 °C (633,15 K), P = 100 atm (10,13 MPa)
- k_ij = 0; çapraz asosiyasyon için Wolbach–Sandler kuralı

## Sonuç
| Büyüklük | Değer |
|---|---|
| Molar yoğunluk | 2,554 kmol/m³ |
| Kütlesel yoğunluk | **99,7 kg/m³** |
| Z | 0,754 |
| η | 0,048 |

PC-SAFT bu koşulda tek kök veriyor; faz yoğun gaz / süperkritik benzeri (sıvı değil).

## Parametreler (Gross & Sadowski, 2002)
| Bileşen | m | σ (Å) | ε/k (K) | ε^AB/k (K) | κ^AB |
|---|---|---|---|---|---|
| Su | 1,0656 | 3,0007 | 366,51 | 2500,7 | 0,034868 |
| 1-Propanol | 2,9997 | 3,2522 | 233,40 | 2276,8 | 0,015268 |

T, P, x veya k_ij değiştirmek için betiğin başındaki değişkenleri düzenleyin.
