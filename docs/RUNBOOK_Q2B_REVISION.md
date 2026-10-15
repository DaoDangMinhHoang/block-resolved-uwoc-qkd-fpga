# Runbook rev B — đo lại, ablation họ phân bố, dose–response Y₀, decoy-state

Ngày soạn: 15/09/2026. Bitstream: **rev B** (`verilog/output_files/top_module.sof`
build 15/09/2026). Runbook cũ `RUNBOOK_Q2_FPGA.md` vẫn đúng cho cách nối dây, công
tắc, thứ tự lệnh và các bẫy; file này chỉ ghi **cái mới** và **kế hoạch đo**.

> **Mọi con số "dự đoán" dưới đây là mô hình ROM-exact** (`python/rom_exact_model.py`,
> dùng đúng ROM và số học nguyên của RTL), **không phải phép đo**. Chúng dùng để chọn
> số block và để có tiêu chí đạt/không đạt. Bài báo chỉ được dùng số đo.

---

## 0. Bitstream rev B khác gì

| Hạng mục | rev A (campaign 09/2026) | rev B |
|---|---|---|
| Nguồn ngẫu nhiên kênh | 5 LFSR 32-bit, dịch 1 bit/clock; 4/5 đa thức **khả quy** (chuỗi h lặp sau 114 681 block, luồng lỗi lặp sau 20.1 s); hai attempt liên tiếp dùng chung 19/24 bit | 5 **xorshift64** (chu kỳ 2⁶⁴−1); 3 luồng quyết định tiến **1 bước/attempt**, 2 luồng fading tiến **1 bước/lần rút** |
| Basis PRNG (engine) | xorshift32 (lặp sau 65 536 block) | xorshift64 |
| ROM h_o | 6 hàng (L0–L5) | **8 hàng**: hàng 6 = Weibull @σ²(L3), hàng 7 = log-normal @σ²(L4) |
| Cường độ | 1 lớp | **decoy xen kẽ từng xung**: signal / decoy (= signal ≫ k) / vacuum; gain 2^g cho signal |
| Gói UART | v5, 96 byte | **v6, 144 byte** (thêm bộ đếm theo lớp), 3.13 ms/gói = 48 % block |
| Tài nguyên | 8 248 LE (44 %), slack +3.689 ns | **11 385 LE (61 %)**, 4 120 register, slack **+2.855 ns**, hold +0.445 ns |
| Testbench | 6/6 | **6/6** (thêm: gói v6 + bộ đếm lớp, lớp decoy, ablation row, từ ngẫu nhiên mới) |

Trước khi nạp: đọc `top_module.sta.summary` — setup slack phải **dương**.

### Lệnh UART mới

| Byte | Ý nghĩa | `send_security_commands.py` |
|---|---|---|
| `0x56` / `0x57` | chọn hàng ablation 6 / 7 (không cộng offset cự ly) | `--turb 6` / `--turb 7` |
| `0x18 \| k` | Y₀ = 218 × {0, ¼, ½, 1, 2, 4, 8, 16}[k] | `--y0 x0.25 … x16` |
| `0x20 \| (p<<2) \| (s−1)` | preset xác suất lớp p, decoy = signal / 2^s | `--decoy-preset p --decoy-shift s` |
| `0x14 \| g` | cường độ signal × 2^g | `--mu-gain 1/2/4/8` |
| `0x74` / `0x75` | tắt / bật decoy | `--decoy off/on` |

Preset lớp: `0` = ½,¼,¼ · `1` = 0.70, 0.20, 0.10 · `2` = 0.80, 0.15, 0.05 · `3` = ⅓ mỗi lớp.

`capture_security_blocks.py` rev E đọc được **cả v5 và v6**, tự kiểm thêm: thứ tự
`E ≤ C ≤ A ≤ N` từng lớp, tổng attempt lớp ≤ tổng, decoy tắt thì bộ đếm lớp phải 0,
vacuum không được click khi Y₀ = 0, tỉ lệ lớp gần preset (cảnh báo).

---

## 1. PHA 0 — Cổng chặn (5 phút)

```powershell
# bring-up như cũ, bắt buộc 100/100 valid
python python\send_security_commands.py --port COM5 --dist 4 --water 0 `
    --lambda-sel 0 --turb 5 --coh 11 --y0 nominal --decoy off --fast on --rate 10M --reset
python python\capture_security_blocks.py --port COM5 --out data\q2b\probe_revB.csv `
    --blocks 100 --print-every 1

# decoy bring-up: phải thấy v6, dec≈20 %, vac≈10 %, vacuum click ≈ 0.85/khối×0.1
python python\send_security_commands.py --port COM5 --dist 4 --water 0 `
    --lambda-sel 0 --turb 5 --coh 11 --y0 nominal --decoy-preset 1 --decoy-shift 3 `
    --mu-gain 4 --decoy on --fast on --rate 10M --reset
python python\capture_security_blocks.py --port COM5 --out data\q2b\probe_decoy.csv `
    --blocks 100 --print-every 1
```

Tiêu chí thêm cho gói v6: `version = 6`, `decoy_enable` đúng lệnh, **0 lỗi lớp**.
Sau probe decoy, **luôn gửi `--decoy off --mu-gain 1`** trước các pha không decoy.

---

## 2. PHA R — Đo lại campaign chính trên rev B (~89 phút)

Lý do: bài báo mô tả rev B là thiết kế chính thức. Số cũ (rev A) giữ làm kiểm tra
tái lập giữa hai bitstream.

```powershell
$levels = @{1 = 286000; 2 = 20000; 3 = 80000; 4 = 40000; 5 = 20000}
foreach ($L in 1,2,3,4,5) {
  foreach ($y0 in "nominal", "zero") {
    python python\send_security_commands.py --port COM5 --dist 4 --water 0 `
        --lambda-sel 0 --turb $L --coh 11 --y0 $y0 --decoy off --mu-gain 1 `
        --fast on --rate 10M --reset
    $tag = if ($y0 -eq "zero") { "y0zero" } else { "nom" }
    python python\capture_security_blocks.py --port COM5 `
        --out "data\q2b\${tag}_L$L.csv" --blocks $levels[$L] --print-every 5000
  }
}
```

L3 lấy **80 000** block và L4 **40 000** block (thay vì 20 000) vì chúng là cặp tham
chiếu của ablation ở pha A; thêm ~9 phút.

**Tiêu chí:** so với rev A (`data/q2/`), `Δ_cg/V` 8 bin nằm trong KTC của nhau; G-test
toàn cục của 5 control Y₀ = 0 cho p > 0.05.

---

## 3. PHA A — Ablation họ phân bố cùng σ² (~22 phút)

Thiết kế 2×2 tại cùng điểm 25 m, 450 nm, μ = 0.1, Y₀ danh định:

| | log-normal | Weibull |
|---|---|---|
| σ²_ho = 0.30 | hàng 3 = L3 (pha R, 80 k) | **hàng 6** (80 k) |
| σ²_ho = 1.00 | **hàng 7** (40 k) | hàng 4 = L4 (pha R, 40 k) |

```powershell
foreach ($row in @(@(6, 80000), @(7, 40000))) {
  foreach ($y0 in "nominal", "zero") {
    python python\send_security_commands.py --port COM5 --dist 4 --water 0 `
        --lambda-sel 0 --turb $($row[0]) --coh 11 --y0 $y0 --decoy off --mu-gain 1 `
        --fast on --rate 10M --reset
    $tag = if ($y0 -eq "zero") { "y0zero" } else { "nom" }
    $n = if ($y0 -eq "zero") { 20000 } else { $row[1] }
    python python\capture_security_blocks.py --port COM5 `
        --out "data\q2b\${tag}_row$($row[0]).csv" --blocks $n --print-every 5000
  }
}
```

**Kiểm tra cấu hình đã vào:** cột `turb` trong CSV phải là 6 / 7 (không phải 3 / 4).

**Dự đoán ROM-exact (8 bin), để so sau khi đo — không đưa vào bài như kết quả:**

| Hàng | Họ @σ² | σ² thực của ROM | Δ_cg/V | BKM/exact |
|---|---|---:|---:|---:|
| 3 | log-normal @0.30 | 0.293 | 1.81×10⁻⁴ | 0.97 |
| 6 | Weibull @0.30 | 0.299 | 2.81×10⁻⁴ | 0.68 |
| 7 | log-normal @1.00 | 0.929 | 5.39×10⁻⁴ | 0.90 |
| 4 | Weibull @1.00 | 0.985 | 1.14×10⁻³ | 0.51 |

Độ phân tách dự kiến giữa hai họ: z ≈ 4 ở σ² = 0.3 với 80 k block/case, z ≈ 8 ở
σ² = 1 với 40 k. Control Y₀ = 0 của hàng 6, 7: G-test p > 0.05.

---

## 4. PHA D — Dose–response theo Y₀ (~26 phút, không cần decoy)

```powershell
foreach ($L in 3, 5) {
  foreach ($k in "x0.25", "x0.5", "x2", "x4", "x8", "x16") {
    python python\send_security_commands.py --port COM5 --dist 4 --water 0 `
        --lambda-sel 0 --turb $L --coh 11 --y0 $k --decoy off --mu-gain 1 `
        --fast on --rate 10M --reset
    python python\capture_security_blocks.py --port COM5 `
        --out "data\q2b\dose_L${L}_$k.csv" --blocks 20000 --print-every 5000
  }
}
```

Điểm ×1 và ×0 lấy từ pha R. Kiểm tra: cột `p_noise` = 54, 109, 436, 872, 1744, 3488.

**Dự đoán ROM-exact Δ_cg/V (8 bin):**

| Y₀ code | 54 | 109 | 218 | 436 | 872 | 1744 | 3488 |
|---|---:|---:|---:|---:|---:|---:|---:|
| L3 | 1.27×10⁻⁵ | 4.94×10⁻⁵ | 1.81×10⁻⁴ | 6.17×10⁻⁴ | 1.92×10⁻³ | 5.31×10⁻³ | 1.30×10⁻² |
| L5 | 4.16×10⁻⁴ | 1.23×10⁻³ | 3.37×10⁻³ | 8.60×10⁻³ | 2.05×10⁻² | 4.57×10⁻² | 9.46×10⁻² |

---

## 5. PHA K — Decoy-state xen kẽ (~52 phút)

Cấu hình chọn theo mô hình (`decoy_security_analysis.py --model`):
μ = 0.4 (`--mu-gain 4`, `mu_level` 8), ν = μ/8 = 0.05 (`--decoy-shift 3`),
preset 1 (70 % signal, 20 % decoy, 10 % vacuum), 80 000 block/điểm.

| Điểm | Lệnh Y₀ | Lý do |
|---|---|---|
| L1, L3, L5 | `nominal` | điểm vận hành của bài |
| L5 | `x4` | giữa |
| L4, L5 | `x16` | nền mạnh — nơi phân giải theo block có ích thật |

```powershell
$points = @(@(1,"nominal"), @(3,"nominal"), @(5,"nominal"), @(5,"x4"), @(4,"x16"), @(5,"x16"))
foreach ($p in $points) {
  python python\send_security_commands.py --port COM5 --dist 4 --water 0 `
      --lambda-sel 0 --turb $($p[0]) --coh 11 --y0 $($p[1]) `
      --decoy-preset 1 --decoy-shift 3 --mu-gain 4 --decoy on --fast on --rate 10M --reset
  python python\capture_security_blocks.py --port COM5 `
      --out "data\q2b\decoy_L$($p[0])_$($p[1]).csv" --blocks 80000 --print-every 5000
}
# Kết thúc: trả board về trạng thái không decoy
python python\send_security_commands.py --port COM5 --decoy off --mu-gain 1 --y0 nominal --reset
```

**Dự đoán mô hình (rate trên attempt, bit/s ở 10 MHz = ×10⁷):**

| Điểm | pooled, asymptotic | block 8 bin, asymptotic | pooled, finite | block 8 bin, finite |
|---|---:|---:|---:|---:|
| L1, Y₀ | 4.59×10⁻⁴ | ×1.000 | 3.85×10⁻⁴ | ×0.64 |
| L5, Y₀ | 4.52×10⁻⁴ | ×1.004 | 3.78×10⁻⁴ | ×0.80 |
| L5, Y₀×4 | 4.09×10⁻⁴ | ×1.031 | 3.38×10⁻⁴ | ×0.87 |
| L4, Y₀×16 | 2.44×10⁻⁴ | ×1.19 | 1.36×10⁻⁴ | ×0.97 |
| **L5, Y₀×16** | 2.43×10⁻⁴ | **×1.35** | 1.35×10⁻⁴ | **×1.51** |

"finite" = mọi gain thay bằng cận Clopper–Pearson xấu nhất, ε tổng 10⁻⁹ (i.i.d. trong
bin, **không** phải finite-key composable). Nhãn quan sát được (click của attempt
lệch basis) cho gần như cùng kết quả với h ẩn.

---

## 6. Phân tích sau khi đo

```powershell
# Pha R + A + D: phase-2 rev E (MI/G-test, scan lồng nhau, nhãn quan sát, model ROM-exact)
python python\phase2_security_analysis.py `
    --csv data\q2b\nom_L1.csv data\q2b\nom_L2.csv data\q2b\nom_L3.csv data\q2b\nom_L4.csv `
          data\q2b\nom_L5.csv data\q2b\nom_row6.csv data\q2b\nom_row7.csv `
          data\q2b\dose_L3_x0.25.csv data\q2b\dose_L3_x0.5.csv data\q2b\dose_L3_x2.csv `
          data\q2b\dose_L3_x4.csv data\q2b\dose_L3_x8.csv data\q2b\dose_L3_x16.csv `
          data\q2b\dose_L5_x0.25.csv data\q2b\dose_L5_x0.5.csv data\q2b\dose_L5_x2.csv `
          data\q2b\dose_L5_x4.csv data\q2b\dose_L5_x8.csv data\q2b\dose_L5_x16.csv `
    --y0-csv data\q2b\y0zero_L1.csv data\q2b\y0zero_L2.csv data\q2b\y0zero_L3.csv `
             data\q2b\y0zero_L4.csv data\q2b\y0zero_L5.csv data\q2b\y0zero_row6.csv `
             data\q2b\y0zero_row7.csv `
    --bins 8 --null-reps 600 --bootstrap 600 --bin-scan --observable-label --distance 25

# Pha K: decoy-state, block vs pooled, asymptotic + finite
python python\decoy_security_analysis.py --csv data\q2b\decoy_L1_nominal.csv `
    data\q2b\decoy_L3_nominal.csv data\q2b\decoy_L5_nominal.csv data\q2b\decoy_L5_x4.csv `
    data\q2b\decoy_L4_x16.csv data\q2b\decoy_L5_x16.csv --bootstrap 200

python python\figs_q2_campaign.py
```

Nhớ: PowerShell **không** bung wildcard cho lệnh native — liệt kê từng file.

Kết quả cần điền vào `Overleaf/main.tex`: tìm `\pending{` (hiện màu đỏ trong PDF).

---

## 7. Tổng thời gian máy

| Pha | Nội dung | Acquisition |
|---|---|---:|
| 0 | cổng chặn | 5 phút |
| R | đo lại L1–L5 + control trên rev B | 89 phút |
| A | hàng 6, 7 + control | 22 phút |
| D | dose–response L3, L5 | 26 phút |
| K | decoy 6 điểm | 52 phút |
| | **tổng** | **≈ 3 h 15 phút** |

Nếu thiếu thời gian: bắt buộc **0, A, K(L5 ×1, L5 ×16)**; R có thể thay bằng
L3 80 k + L4 40 k + control L4 (≈ 20 phút) và giữ số rev A cho các mức còn lại, ghi rõ
trong bài hai bitstream.
