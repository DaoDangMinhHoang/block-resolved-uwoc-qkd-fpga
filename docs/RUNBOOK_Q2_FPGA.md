# Runbook — nâng cấp lên Q2 và đo trên DE1

Mọi con số trong file này lấy từ RTL và script thật, không phải ước lượng.

**Tiền đề:** 6/6 testbench đạt, timing đóng (+3.689 ns), fitter 8 248 LE (44 %).

> **Trạng thái 11/09/2026 — MỌI PHA ĐÃ CHẠY XONG**, kể cả 4B, 4B′ và pha 3′
> (control `Y₀ = 0` ở L2, L3). Control nay có ở cả 5/5 mức. Runbook giữ nguyên để
> tái lập. Kết quả: `README.md` §6.6, `BAO_CAO_CAMPAIGN_Q2_2026-09-08.md` (§10 cho
> pha 2, §12 cho pha 3′).

**Tổng thời gian máy chạy:**

| Pha | Nội dung | Thời gian | Trạng thái |
|---|---|---:|:--:|
| 1 | Bring-up probe, 100 block — **cổng chặn** | 1 phút | ✅ 13/13 giá trị khớp |
| 3 | Control Y₀ = 0, L4 + L5 | 5 phút | ✅ `data/q2/y0zero_L{4,5}.csv` |
| 3′ | Control Y₀ = 0, L2 + L3 (bổ sung cho đủ 5 mức) | 5 phút | ✅ 11/09 — `data/q2/y0zero_L{2,3}.csv`, KTC chứa 0, QBER = `e_pol` |
| 4A | Campaign L1–L5, 20 000 block/mức | 11 phút | ✅ `data/q2/nom_L{1..5}.csv` |
| 5 | Phân tích | 5 phút | ✅ có cả `--bin-scan` |
| | **đủ để có toàn bộ kết quả cho bài báo** | **~22 phút** | ✅ |
| 2 | Cross-check OOK vs fast (**phải dùng `fpga_collect.py`**, không phải `bb84_uwoc_measure.py`) | +15 phút | ✅ 09/09, gộp 5 mức: G lệch **1.5 %** |
| 4B | Bổ sung để phân giải L1 (~~và L2~~) | +31 phút | ✅ 09/09 — L1 phân giải, KTC hở khỏi 0 ≈**6 σ** |
| 4B′ | Control `Y₀ = 0` ghép cặp tại L1 | +31 phút | ✅ 09/09 — KTC **chứa 0**, cận trên 11.4× dưới tín hiệu |

Pha 1 là **cổng chặn** — không đạt 100/100 thì đừng chạy tiếp. Sau đó chạy pha 3
và 4A trước (đây là phần khoa học), pha 2 để sau (đây là phần phòng thủ).

> **Pha 4B rút còn một nửa — và nửa còn lại đã chạy 09/09, kết quả tốt.** Nó viết
> ra để phân giải L1 *và* L2; pha 4A đã phân giải L2 ngay ở 20 000 block
> (`Δ_cg/V = 5.60×10⁻⁵`, KTC `[3.16, 10.49]×10⁻⁵`), nên chỉ còn nhánh L1 —
> 286 000 block, 31 phút 14 giây, ghi ra file mới nên không mất gì.
>
> Lý do chạy **không** phải là "có thêm một điểm cho đẹp". Nó là: tỉ số **×4137
> đang chia cho nhiễu**. Mẫu số `8.402×10⁻⁷` nằm dưới nền null `8.24×10⁻⁶`, và
> `--bin-scan` xác nhận L1 chưa phân giải là **do thiếu block, không do cách chia
> bin**.
>
> **Kết quả:** `Δ_cg/V = 1.531×10⁻⁵`, KTC `[1.15, 2.13]×10⁻⁵`, hở khỏi 0
> **≈6 σ** (dự báo 5.1 σ), khớp mô hình 14 %. Dải động thành **×227**. Nền null
> rơi xuống `5.85×10⁻⁷` — hệ số 14.45 so với tỉ số block 14.3, đúng 1/N.
>
> **Control ghép cặp cùng 286 000 block** cho `−1.85×10⁻⁷`, KTC
> `[−7.29×10⁻⁷, +1.49×10⁻⁶]` **chứa 0**, cận trên 11.4 lần dưới tín hiệu. Chi
> tiết: `BAO_CAO_CAMPAIGN_Q2_2026-09-08.md` §11.

> ⚠️ **Bẫy ở pha 5:** PowerShell **không** bung `nom_L*.csv` cho lệnh native, và
> script cũng không gọi `glob` (`--csv` chỉ là `nargs="*"`). Viết wildcard thì
> script nhận 0 file và im lặng bỏ qua. Phải liệt kê từng file — xem §5.

---

## 0. Chuẩn bị

### 0.1 Build và nạp

```powershell
# Quartus II 13.0 GUI: mở verilog\top_module.qpf -> Compile -> Programmer -> top_module.sof
```

Trước khi nạp, **đọc lại 2 con số** trong `verilog\output_files\`:

| File | Phải thấy |
|---|---|
| `top_module.sta.summary` | Slow Model Setup slack **dương** (hiện +3.689 ns), TNS 0.000 |
| `top_module.fit.summary` | 8 248 / 18 752 LE, 8/52 multiplier |

Nếu slack âm thì dừng: mọi số đo sau đó không bảo vệ được.

### 0.2 Công tắc DE1

| Công tắc | Đặt | Lý do |
|---|:--:|---|
| `SW[9]` | **1** | PC mode. Ở fast mode `pc_input_mode = SW[9] & ~fast` tự về 0, nhưng để 1 thì lúc `--fast off` (pha 2) board mới nghe lệnh qubit |
| `SW[4]` | **1** | **Bật kênh UWOC.** Để 0 là bypass lý tưởng — sẽ ra 100 % click, 0 lỗi |
| `SW[1]` | **0** | Fixed. Bật adaptive thì controller đổi μ/λ giữa chừng → `config_changed = 1` → hỏng block |
| `SW[0]` | 0 | Automatic |
| `SW[7:5]`, `SW[3:2]` | tuỳ ý | Không dùng nữa. `[Q2]` reset về hằng số, cấu hình đến từ UART |

### 0.3 Thứ tự khởi động — quan trọng

1. **Nhấn `KEY[3]`** (reset cứng). Lúc này cấu hình về mặc định: clear ocean,
   dist 0, L1, 450 nm, **`coh_sel = 0` (bộ đếm đồng hồ cũ)**, fast off,
   `rate_sel = 6`, Y₀ nominal.
2. Chạy `send_security_commands.py`. Script gửi cấu hình trước, **`0x01` cuối cùng**.

`0x01` reset *thống kê*, không reset *cấu hình* — đó là lý do nó phải đi cuối:
nó khởi động lại `qub_index` đúng lúc `security_block_monitor` clear, và chính
sự đồng bộ đó làm `h_min == h_max` thành **bất biến** chứ không phải bộ lọc.

> **`--coh 11` là bắt buộc mỗi lần cấu hình.** Reset cứng đặt `cfg_coh = 0`.
> Quên `--coh 11` thì `h_min != h_max` trên **mọi** block và toàn bộ capture bị loại.

---

## 1. PHA 1 — Bring-up probe (cổng chặn, ~1 phút)

```powershell
# Terminal 1
python python\send_security_commands.py --port COM5 --baud 460800 `
    --dist 4 --water 0 --lambda-sel 0 --turb 5 --coh 11 `
    --y0 nominal --fast on --rate 10M --reset

# Terminal 2 (ngay sau khi lệnh trên in "done")
python python\capture_security_blocks.py --port COM5 --baud 460800 `
    --out data\q2\probe.csv --blocks 100 --print-every 1
```

`--dist 4` = **25 m** trên lưới clear ocean `(5,10,15,20,25,…)`. `--water 0`
clear ocean, `--lambda-sel 0` = 450 nm.

### 13 giá trị phải khớp

| Trường | Kỳ vọng | Sai thì nghĩa là gì |
|---|---|---|
| `attempts` | **65536 chính xác** | `ATTEMPT_LOG2` không phải 16 |
| `clicks_all + lost_all` | 65536 | Bộ đếm attempt hỏng |
| `trials_z + trials_x` | ≈ 32768, và **≫ clicks** | `trials == clicks` ⇒ **PATCH item 3 chưa áp** — `evt_bmatch` chỉ bật khi có click, mọi `G_j` sẽ ra 1.0 mà chuỗi thứ tự vẫn pass. Đây là lỗi im lặng nguy hiểm nhất |
| `clicks_z + clicks_x` | ≈ **39** ở 25 m | Lệch nhiều ⇒ sai `dist_idx` hoặc `SW[4] = 0` |
| `errors_z + errors_x` | 1–2 | |
| `h_min == h_max` | **100/100** | Sai 1 block cũng là bug. `coh_sel != 11` hoặc sai thứ tự reset. **All-or-nothing** — pass 60 % là bug, không phải 60 % dữ liệu dùng được |
| `block_synced` | 1 | Reporter chốt cờ trễ một clock |
| `attempts_per_sec` | **10 000 000 ± 1** sau khi `rate_valid = 1` | Đây là con số chứng minh luận điểm real-time |
| `rate_valid` | 1 sau giây đầu (~153 block) | |
| `engine_overruns` | **0** | Kênh không theo kịp nhịp yêu cầu → tin `attempts_per_sec`, đừng tin `rate_hz` |
| `engine_timeouts` | **0** | Kênh không trả lời |
| `dropped_blocks`, `watchdog_aborts` | **0** | UART mất quyền giữa gói |
| `protocol_errors`, `config_changed` | **0** | `qub_index` nhảy cóc / cấu hình đổi trong block |
| `block_id` | liên tục (`block_gap = 0`) | |

Dòng cuối script in `packets=100 valid=100 (100.0%)`. **Không đạt 100 % thì dừng
lại đọc cột `validation_error`.**

---

## 2. PHA 2 — Cross-check `--fast off` vs `--fast on` (~25 phút)

Đây là bằng chứng bảo vệ luận điểm *"đường tăng tốc không làm đổi thống kê kênh"*.
Đường fast **bỏ qua tầng OOK serialisation đã validate**, nên phải chứng minh bằng số.

> ### ✅ Đã chạy 09/09/2026 — ĐẠT
>
> Chân OOK `--scale 0.05` (8 block/mức) đối chiếu chân fast 2 000 block/mức.
> Phân tích: `python\xcheck_fast_vs_ook.py`.
>
> | Gộp 5 mức | OOK | fast | Tỉ số | KTC 95 % |
> |---|---:|---:|---:|---|
> | `Q_all` = click/attempt | 1.1650×10⁻³ | 1.1848×10⁻³ | 0.983 | [0.781, 1.212] |
> | `G` = sifted/attempt | 5.8365×10⁻⁴ | 5.9268×10⁻⁴ | 0.985 | [0.786, 1.210] |
> | QBER | 3.987 % | 3.691 % | +0.296 pp | z = 0.59 |
>
> Chân OOK: 2.62×10⁶ attempt, 1 530 sifted. Chân fast: 6.55×10⁸ attempt,
> 388 421 sifted, 10 000/10 000 block `capture_valid = 1`.
>
> **Phát biểu được:** hai đường không lệch nhau quá **22 %** — đó là độ phân giải
> của phép đo, không phải sai số đo được. **Không** phát biểu được "khớp trong vài
> phần trăm". Muốn siết xuống ≈ ±5 % thì đo lại chân OOK ở `--scale 1.0`
> (160 block/mức, ~4 giờ).

> ### ⛔ ĐỪNG dùng `bb84_uwoc_measure.py` cho pha này
>
> Nó chạy **~86 qubit/s**, tức 8 tiếng cho 5 mức. Không phải lỗi board.
> `FPGASource.collect()` gọi `ser.readline()` **một lần cho MỖI qubit**
> (`bb84_uwoc_measure.py:202`). Dưới nước >99.8 % qubit **không** click, nên gần
> như mọi `readline()` đều chờ hết timeout ⇒ tốc độ sập về `1/timeout`, bất kể
> FPGA chạy nhanh thế nào. Timeout còn bị chặn dưới ở 4 ms
> (`self.timeout = max(0.004, …)`) nên trần lý thuyết chỉ 250 q/s.
>
> `fpga_collect.py:130–137` ghi rõ cái bẫy này và giải bằng `_harvest()` — đọc
> **không chặn** cả buffer OS bằng `ser.in_waiting`, nên giới hạn quay về
> ~220 µs/qubit của chính FPGA (~3 300–3 900 q/s). Comment
> `fpga_collect.py:92` thậm chí ghi lại đúng con số: *"calibrated on COM10 —
> 86.6 qubits/s"*.
>
> `bb84_uwoc_measure.py` là công cụ **quét thăm dò**, không phải công cụ campaign.

```powershell
# 0) BẮT BUỘC: sao lưu data/ trước — xem cảnh báo bên dưới
Copy-Item -Recurse data data_backup_v14

# 1) Chân A — đường OOK gốc, dùng ĐÚNG công cụ campaign
python python\fpga_collect.py --port COM5 --chunk 32 --chunk-check   # đo tốc độ thật trước
python python\fpga_collect.py --port COM5 --phase fixed `
    --sections B --scale 0.05 --chunk 32 --coh 11
```

Mục **B** của `build_matrix()` chính là quét nhiễu loạn L1–L5 tại `dist_idx = 4`
(25 m), clear ocean, 450 nm — đúng điểm vận hành của pha 4. Ngân sách gốc là
`target_qubit = 160 × 65 536` = 1.05e7 attempt/mức ⇒ 47 phút/mức ở 3 700 q/s.
`--scale 0.05` rút xuống **524 288 attempt/mức ≈ 2.4 phút/mức, 12 phút cả 5 mức**.

> ### ⚠️ Nguy cơ mất dữ liệu — đọc trước khi chạy
>
> `save_clicks()` ghi `data/clicks_<tag>.csv` ở chế độ **`"w"` (ghi đè)**. Tag của
> mục B là `fixed_B_turb_L1…L5` — **trùng đúng 5 file** mà
> `defect_from_clicklog.py` đang dùng. Sao lưu `data/` trước.
>
> Ngược lại, `load_done()` đọc tag từ `data/fpga_points.csv`, nên 5 điểm này
> **đã nằm trong danh sách "done" và sẽ bị BỎ QUA**. Muốn đo lại thật thì phải xoá
> 5 dòng `fixed_B_turb_L*` khỏi `fpga_points.csv` — và lúc đó cảnh báo ghi đè ở
> trên thành hiện thực. Làm bản sao lưu trước, không có ngoại lệ.

> ### Latency Timer — **chỉ có trên FTDI**
>
> Mặc định FTDI trên Windows là **16 ms**. `fpga_collect.py:644` ghi:
> *"Cổng còn Latency Timer 16 ms (~63 qubit/s) thì dùng `--cap-scale 4`."*
> Sửa hẳn thì tốt hơn: **Device Manager → Ports (COM & LPT) → cổng của bạn →
> Properties → Port Settings → Advanced → Latency Timer = 1**.
>
> **Adapter Prolific PL2303 KHÔNG có mục này** — đừng mất công tìm. Nó cũng không
> cần: cổng Prolific đo được **3 700 qubit/s**, đúng trần ~220 µs/qubit của chính
> FPGA. Nếu `--chunk-check` ra thấp bất thường trên Prolific thì hai nút tương
> đương là **Receive Buffer → Low (1)** trong cùng hộp thoại Advanced, và tắt
> *"Allow the computer to turn off this device"* ở Power Management của USB Root
> Hub. Ô **Bits per second** trong hộp thoại là vô nghĩa — pyserial ghi đè thành
> 460 800 lúc mở cổng.

### Ba cái bẫy của riêng pha 2

**1. `Link.configure()` không bao giờ gửi `0x70`/`0x71`.** Nó chỉ gửi 5 byte
`0x30/0x40/0x50/0x60/0x01`, nên board **giữ nguyên trạng thái fast của lần chạy
trước**. Vừa chạy campaign xong mà vào thẳng chân A thì `pc_input_mode = 0`, board
nuốt sạch lệnh qubit, và 12 phút sau bạn nhận file rỗng. Gửi
`send_security_commands.py --fast off` (hoặc nhấn `KEY[3]`) trước chân A.

**2. Đừng xoá dòng checkpoint bằng `Set-Content -Encoding utf8`.** PowerShell 5.1
ghi UTF-8 **kèm BOM**; `csv.DictReader` khi đó đọc cột đầu thành `﻿tag` và
`load_done()` chết với `KeyError: 'tag'`. Dùng .NET:
> ```powershell
> $t = [System.IO.File]::ReadAllText($p, [System.Text.Encoding]::UTF8)
> [System.IO.File]::WriteAllText($p, $t, (New-Object System.Text.UTF8Encoding($false)))
> ```

**3. Đơn vị lấy lại mẫu là BLOCK, và `--scale 0.05` chỉ cho 8 block/mức.**
Sai số của phép so sánh do **số block fading** quyết định, không phải số attempt:
attempt trong cùng một block là 65 536 lần lặp của *một* lần bốc h. Dùng sai số
nhị thức thì L4 ra `|z| = 7.7` và L5 ra `5.3` cho những gain **thực ra khớp nhau**.
Nhưng bootstrap theo block trên 8 đơn vị cũng chưa đáng tin theo chiều ngược lại —
nó cho L4/QBER `TRUOT` trong khi z nhị thức chỉ 1.03. Cách đúng: **gộp cả 5 mức**
(40 đơn vị) — câu hỏi của pha 2 là về *đường dữ liệu*, mà đường dữ liệu không biết
nó đang chạy ở mức nhiễu loạn nào. `xcheck_fast_vs_ook.py` in sẵn hàng `GOP` và
cảnh báo mọi hàng dưới 20 block.

```powershell
# Chân B — đường fast, cùng 5 mức, cùng điểm vận hành
foreach ($L in 1,2,3,4,5) {
  python python\send_security_commands.py --port COM5 --dist 4 --water 0 `
      --lambda-sel 0 --turb $L --coh 11 --y0 nominal --fast on --rate 10M --reset
  python python\capture_security_blocks.py --port COM5 `
      --out "data\q2\xcheck_fast_L$L.csv" --blocks 2000 --print-every 500
}
```

2 000 block = 1.31×10⁸ attempt/mức, **13 giây/mức**. Chú ý mức chênh: đường fast
thu **250 lần** nhiều hơn đường OOK trong thời gian **ngắn hơn 11 lần** (13 s so
với 141 s) — đó chính là điểm của lớp Q2, và cũng là con số đáng đưa vào bài.

```powershell
# Khôi phục 5 log click gốc — BẮT BUỘC, chúng là đầu vào của defect_from_clicklog.py
foreach ($L in 1,2,3,4,5) {
  Move-Item "data\clicks_fixed_B_turb_L$L.csv" "data\q2\xcheck_ook_L$L.csv" -Force
  Copy-Item "data_backup_v14\clicks_fixed_B_turb_L$L.csv" "data\" -Force
}

python python\xcheck_fast_vs_ook.py
```

Bản `--scale 0.05` chỉ có 8 block/mức nên **không thay thế được** log 160 block
gốc. Kiểm tra bằng `Get-FileHash` sau khi khôi phục.

**Tiêu chí:** hàng **`GOP`** (gộp 5 mức) phải chứa 1.0 trong KTC của tỉ số cho cả
`Q_all` và `G`, và chứa 0 cho hiệu QBER. Từng mức riêng lẻ ở `--scale 0.05` có độ
phân giải ±30…68 %, **không dùng để kết luận** — script tự gắn nhãn `(it block)`
cho chúng. Ghi vào bài báo hàng `GOP` kèm độ phân giải, đừng ghi bảng 5 mức như
thể mỗi dòng là một bằng chứng độc lập.

> `bb84_uwoc_measure.py` tự gửi byte cấu hình trong `src.configure()`, nên chỉ cần
> `--fast off` trước đó. Nó **không** gửi `--coh`, nhưng ở chân A ta chỉ so pooled
> `G`/QBER nên lưới coherence không ảnh hưởng.

---

## 3. PHA 3 — Control Y₀ = 0 (~5 phút) — chạy TRƯỚC campaign

Đây là bằng chứng mạnh nhất cho phần toán, và **không tài liệu nào của thầy nghĩ ra**.

Γ₀ hạng 1 với hướng xuyên tâm nằm trong nhân ⇒ khi `Y₀ = 0` thì `T = e_pol·G`
**chính xác**, quỹ đạo `(G,T)` là một tia qua gốc, và `Δ_cg = 0` *đúng bằng 0*.
Nếu đo ra 0 thì defect ở campaign chính **bắt buộc** là hiệu ứng nền Y₀, không
thể là dư của nhiễu đếm.

```powershell
foreach ($L in 2,3,4,5) {
  python python\send_security_commands.py --port COM5 --dist 4 --water 0 `
      --lambda-sel 0 --turb $L --coh 11 --y0 zero --fast on --rate 10M --reset
  python python\capture_security_blocks.py --port COM5 `
      --out "data\q2\y0zero_L$L.csv" --blocks 20000 --print-every 500
}
```

20 000 block × 6.55 ms ≈ **2.2 phút/mức**. Bản đầu (08/09) chỉ chạy L4, L5; ngày
11/09 bổ sung L2, L3 để bảng đối chứng đủ 5 mức như bảng nominal (L1 có control riêng
286 000 block ở giai đoạn B′). **Chạy lại từ đầu thì chạy luôn 2, 3, 4, 5** như lệnh
trên — thêm 4 phút, và không ai phải hỏi "L2, L3 đâu".

**Tiêu chí:** `Δ_cg/V` của các file này phải **chứa 0** trong khoảng tin cậy.

> **Sau khi chạy giai đoạn B thì cần thêm một control nữa.** L1 phân giải ở
> 286 000 block, nên nó cũng cần control ở **cùng cỡ mẫu** — control 20 000 block
> không đối chứng được cho một điểm 286 000 block, vì nền của phương pháp giảm
> theo 1/N. Xem §4 giai đoạn B′.

---

## 4. PHA 4 — Campaign chính L1–L5 (~11 phút giai đoạn A)

Bảng planning (`phase2_results\model_expectation.csv`) nói cần:

| Mức | block cần | thời gian @10 MHz |
|---|---:|---:|
| L1 | 285 623 | 31.2 phút |
| L2 | 71 975 | 7.9 phút |
| L3 | 19 851 | 2.2 phút |
| L4 | 2 856 | 18.7 s |
| L5 | 1 089 | 7.1 s |
| **tổng** | **381 394** | **41.7 phút** |

Nhưng bảng này tự cảnh báo *"PLANNING ONLY, hệ số 20 là hiệu chỉnh kinh nghiệm,
không phải cận sample-complexity"*. Nên chia hai giai đoạn:

### Giai đoạn A — 20 000 block/mức, ~11 phút, chạy trước

```powershell
foreach ($L in 1,2,3,4,5) {
  python python\send_security_commands.py --port COM5 --dist 4 --water 0 `
      --lambda-sel 0 --turb $L --coh 11 --y0 nominal --fast on --rate 10M --reset
  python python\capture_security_blocks.py --port COM5 `
      --out "data\q2\nom_L$L.csv" --blocks 20000 --print-every 500
}
```

20 000 block đã vượt yêu cầu của L3 (19 851), L4 (2 856) và L5 (1 089) — tức là
**giai đoạn A phân giải được L3, L4, L5**, tốt hơn hẳn kết quả offline hiện tại
(chỉ phân giải L4, L5 với 160 block).

### Giai đoạn B — ~~L1 + L2~~ **chỉ còn L1** (~31 phút)

> **Nhánh L2 đã thừa.** Đoạn này viết trước khi có kết quả 4A. Ở 20 000 block, L2
> đã phân giải: `Δ_cg/V = 5.60×10⁻⁵`, KTC `[3.16, 10.49]×10⁻⁵`, `resolved = 1`.
> Chạy thêm 72 000 block cho L2 là phí 7.9 phút. **Chỉ chạy cặp lệnh L1.**

```powershell
python python\send_security_commands.py --port COM9 --dist 4 --water 0 `
    --lambda-sel 0 --turb 1 --coh 11 --y0 nominal --fast on --rate 10M --reset
python python\capture_security_blocks.py --port COM9 --out data\q2\nom_L1_full.csv `
    --blocks 286000 --print-every 5000
```

Ghi ra **file mới**, không đè `nom_L1.csv`: 154 B/block × 286 000 = **44 MB**,
286 000 × 6.5536 ms = **31 phút 14 giây**. Không cần backup, không cần `KEY[3]` —
`send_security_commands.py` ghi đè toàn bộ cấu hình rồi `0x01` cuối cùng.

**`--print-every` không được nhỏ hơn 5000.** Board phát 153 gói/giây; in dày hơn
thì console Windows chậm hơn UART, buffer serial tràn, gói vỡ.

Vì sao 286 000 là đủ (và dư): model cho `Δ_cg/V = 1.343×10⁻⁵` ở L1
(`phase2_results/model_expectation.csv`). Ở 20 000 block σ ≈ 1.0×10⁻⁵ nên KTC
chứa 0; σ giảm theo 1/√N, ở 286 000 block còn **2.64×10⁻⁶**, tức tín hiệu tách
khỏi 0 khoảng **5.1 σ**. Ngưỡng tối thiểu chỉ là ~43 000 block.

> ### ✅ Đã chạy 09/09/2026 — và tỉ số ×4137 đã đổi thành ×227
>
> `Δ_cg/V = 1.531×10⁻⁵`, KTC `[1.15, 2.13]×10⁻⁵`, `resolved = 1`, hở khỏi 0
> **≈6 σ**. Capture 286 000/286 000 block valid, `attempts_per_sec = 10 000 000`
> trên mọi block. Đo/mô hình = 1.14.
>
> **Đừng trích σ hai chữ số.** Nó đọc từ KTC bootstrap 600 lần lặp, bề rộng KTC
> mang ~15 % sai số Monte-Carlo: hai lần chạy trên **cùng** capture cho 7.1 σ và
> 6.2 σ. Trích KTC, không trích σ.
>
> ×4137 cũ chia cho điểm ước lượng L1 `8.402×10⁻⁷`, mà điểm ấy nằm **dưới** nền
> null `8.24×10⁻⁶` — nó là nhiễu, không phải phép đo. Nếu **không** chạy giai đoạn
> B thì phát biểu đúng phải dùng cận trên 95 % làm mẫu số:
> `3.4758×10⁻³ / 2.6607×10⁻⁵` ⇒ *"tăng **ít nhất 131 lần**"*. Dù đường nào,
> **×4137 không được đưa vào bài.**
>
> Phần thưởng kèm theo: nền null ở L1 rơi từ `8.24×10⁻⁶` xuống `5.85×10⁻⁷`, hệ số
> **14.45** so với tỉ số block **14.3** — đúng 1/N của thiên lệch mẫu hữu hạn, một
> kiểm chứng miễn phí cho cách trừ nền.

### Giai đoạn B′ — control `Y₀ = 0` ghép cặp tại L1 (~31 phút)

Đây là thứ đóng lỗ hổng cuối cùng của lập luận về L1: sau giai đoạn B, điểm L1 chỉ
dựa vào nền tham số cộng mức khớp mô hình, **chưa có control trực tiếp ở cùng cỡ
mẫu**. Pha 3 bản đầu chỉ chạy control ở L4, L5 (L2, L3 bổ sung
11/09 ở 20 000 block — không thay được control 286 000 block cho L1).

Đổi **đúng một byte** so với giai đoạn B — `0x72` (`--y0 zero`) thay `0x73`:

```powershell
python python\send_security_commands.py --port COM9 --dist 4 --water 0 `
    --lambda-sel 0 --turb 1 --coh 11 --y0 zero --fast on --rate 10M --reset
python python\capture_security_blocks.py --port COM9 `
    --out data\q2\y0zero_L1_full.csv --blocks 286000 --print-every 5000
```

**Dấu hiệu cấu hình đã vào:** QBER in ra console tụt về **≈3.18 %** (= `e_pol`)
chứ không phải 3.69 %. Tắt nền đếm thì lỗi chỉ còn do phân cực.

Tên file `y0zero_L1_full.csv` **không** đụng pattern `y0zero_L%d.csv` mà
`fig_control_gt` đọc (hàm đó chỉ lấy L4, L5 — và nên để vậy: ở L1 quỹ đạo (G,T)
quá ngắn, vẽ ra không nói được gì). Hình chính thì tự nhận, vì nó chọn control
theo `p_noise == 0` chứ không hardcode mức.

> ### ✅ Đã chạy 09/09/2026 — KTC chứa 0, đúng như định lý đòi
>
> | L1, cùng 286 000 block | `Δ_cg/V` | KTC 95 % | `resolved` |
> |---|---:|---|:--:|
> | `Y₀` nominal | 1.531×10⁻⁵ | [1.15, 2.13]×10⁻⁵ | ✅ hở khỏi 0 ≈6 σ |
> | **`Y₀` = 0** | **−1.85×10⁻⁷** | **[−7.29×10⁻⁷, +1.49×10⁻⁶]** | **0 — chứa 0** ✅ |
>
> Cận trên một phía của control `1.34×10⁻⁶`, tức **11.4 lần dưới** tín hiệu
> nominal — tốt hơn dự báo (≈4×10⁻⁶, 3.8 lần), vì ở L1 nhiễu loạn yếu
> (CV(h) ≈ 0.20) nên phương sai estimator nhỏ hơn nhiều so với ngoại suy từ
> control L4/L5 (CV ≈ 0.92).
>
> Điểm ước lượng **âm** không phải vi phạm: `Δ_cg ≥ 0` là định lý, nên điểm âm là
> nhiễu estimator — đúng lý do quy ước báo cáo dùng cận trên một phía.
>
> Đây là **đối chứng ghép cặp** đúng nghĩa: cùng cỡ mẫu, cùng phân bố fading, cùng
> mọi tham số, khác một byte. Bật nền thì hở khỏi 0; tắt nền thì straddle 0.

---

## 5. PHA 5 — Phân tích

```powershell
python python\phase2_security_analysis.py `
    --csv data\q2\nom_L1.csv data\q2\nom_L2.csv data\q2\nom_L3.csv `
          data\q2\nom_L4.csv data\q2\nom_L5.csv `
    --y0-csv data\q2\y0zero_L4.csv data\q2\y0zero_L5.csv `
    --bins 8 --null-reps 600 --bootstrap 600 --bin-scan --distance 25
python python\figs_q2_campaign.py
```

**Đã chạy giai đoạn B/B′ thì dùng bản này** — thay `nom_L1.csv` →
`nom_L1_full.csv` và thêm control L1 vào `--y0-csv`:

```powershell
python python\phase2_security_analysis.py `
    --csv data\q2\nom_L1_full.csv data\q2\nom_L2.csv data\q2\nom_L3.csv `
          data\q2\nom_L4.csv data\q2\nom_L5.csv `
    --y0-csv data\q2\y0zero_L1_full.csv data\q2\y0zero_L2.csv data\q2\y0zero_L3.csv `
             data\q2\y0zero_L4.csv data\q2\y0zero_L5.csv `
    --bins 8 --null-reps 600 --bootstrap 600 --bin-scan --distance 25
python python\figs_q2_campaign.py
```

**Đừng đưa cả `nom_L1.csv` và `nom_L1_full.csv` vào** — chúng là cùng một điểm vận
hành, gộp lại là đếm hai lần. Lần chạy này lâu hơn ~5 phút của giai đoạn A đáng kể:
bootstrap 600 lần trên hai file 286 000 block thay vì toàn 20 000.

**Hai script này không cần board** — chỉ đọc CSV, không `import serial`, không có
`--port`. Tháo mạch rồi vẫn chạy được.

**Liệt kê tường minh, đừng dùng `nom_L*.csv`.** PowerShell không bung wildcard cho
lệnh native và script không gọi `glob` — viết wildcard thì `--csv` nhận đúng chuỗi
`data\q2\nom_L*.csv`, không khớp file nào, và script chạy tiếp mà không báo lỗi.

Sinh ra trong `phase2_results\`:

| Sản phẩm | Nội dung |
|---|---|
| `security_summary.csv` | Bảng chính: `excess_over_V`, KTC, `resolved`, `bkm_excess_over_V`, `mean_rate` |
| `defect_vs_turbulence_*.png` | **Hình chính** — Δ_cg theo mức nhiễu loạn, kèm đường control Y₀ = 0 |
| `exact_vs_bkm_*.png` | Exact vs ½Tr(Γ₀Cov) |
| `CONTROL_GT_turb*.png` | **Quỹ đạo (G,T) của Y₀ = 0 và Y₀ nominal trên cùng trục** — hình giải thích toàn bộ kết quả |
| `GT_*.png`, `QBER_h_*.png` | Quỹ đạo và QBER theo irradiance từng mức |
| `basis_symmetry_*.png` | e_Z vs e_X — biện minh cho rút gọn 2-D bằng dữ liệu |
| `gllp_intensity_scan.csv` + hình | Quét μ: điểm tối ưu click-rate ≠ điểm tối ưu bảo mật |
| `bin_scan.csv/.png` | Độ bền theo số bin |

---

## 6. Thông số phải đo — tóm tắt

### 6.1 Nhóm validate (phải bằng 0 hoặc đúng chính xác, sai thì capture vô hiệu)

`attempts = 65536` · `clicks_all + lost_all = attempts` · `trials > clicks`
**nghiêm ngặt** · `h_min == h_max` 100 % · `block_synced = 1` ·
`config_changed = 0` · `protocol_errors = 0` · `dropped_blocks = 0` ·
`engine_timeouts = 0` · `engine_overruns = 0` · `watchdog_aborts = 0` ·
`block_gap = 0` · `attempts_per_sec = 10 000 000 ± 1`

Script tự chấm 12 cổng này và ghi lý do vào cột `validation_error`.

### 6.2 Nhóm vật lý (thứ đi vào bài báo)

| Đại lượng | Nguồn | Dùng để |
|---|---|---|
| `trials_z/x`, `clicks_z/x`, `errors_z/x`, `attempts` | mỗi block | Dựng `m_j = (G_Z, E_Z, G_X, E_X)` theo §7 của thầy Hòa |
| `h_mean` (12 bit, 256 ≙ 1.0, trần **h = 16**) | mỗi block | Chia bin irradiance — **không bão hoà** như `irrad` 8 bit của log cũ (trần h = 2.0) |
| `G = clicks/attempts`, `T = errors/attempts` | dẫn xuất | Per **attempt**, đúng trọng số Prop. 1 |
| **`Δ_cg/V` + KTC bootstrap** | 5 mức | **Hình quan trọng nhất** (§7 Technical Assessment) |
| `½Tr(Γ₀Cov)/V` | 5 mức | Exact vs BKM — chỗ tách ra là chuyển tiếp phi tuyến |
| `Δ_cg` ở Y₀ = 0 | L1–L5 | **Phải chứa 0** — chứng minh defect là vật lý. Kết quả 11/09: chứa 0 ở L1, L2, L3, L5; L4 hở biên (~1 % tín hiệu) |
| QBER ở Y₀ = 0 | L1–L5 | Phải bằng `e_pol(25 m) = 3.178 %` ở mọi mức. Kết quả 11/09: đạt, χ² = 2.64 / 5 bậc tự do |
| `e_Z` vs `e_X` | 5 mức | Emulator áp `e_pol` độc lập cơ sở ⇒ phải bằng nhau. Biện minh rút gọn 2-D |
| **`attempts_per_sec`** | mọi block | Chứng minh real-time: 10 MHz ⇒ block 6.55 ms ⇒ **vòng điều khiển 153 Hz** |
| `Q_all_pool` = `clicks_all/attempts` | mỗi mức | GLLP cần xác suất detection **trên mỗi xung phát**, không phải sifted gain |

---

## 7. Ba cái bẫy thực tế

**1. `--print-every 1` sẽ làm hỏng capture dài.** Ở 10 MHz board phát
**153 gói/giây**. In 153 dòng/giây ra console Windows chậm hơn tốc độ UART, buffer
serial tràn, gói vỡ. Dùng `--print-every 500` trở lên cho mọi lần chạy > 1 000 block.
Pha 1 (100 block) thì `--print-every 1` mới hợp lý.

**2. Không mở được hai chương trình trên cùng COM.** Windows khoá cổng độc quyền,
nên phải chạy `send_security_commands.py` xong (nó tự đóng cổng) rồi mới chạy
`capture_security_blocks.py`. Khoảng trống 1–2 giây giữa hai lệnh làm mất ~150–300
block đầu — **không sao**: `read_packet` tự đồng bộ lại theo magic `A5 5A`, và
`block_gap` tính theo gói bắt được liền trước nên không báo lỗi giả.

**3. Ghi CSV trên ổ chậm.** `capture_security_blocks.py` gọi `f.flush()` sau mỗi
dòng. Ở 153 dòng/giây trên ổ cứng cơ sẽ nghẽn. Ghi vào SSD. Dung lượng: giai đoạn A
là 100 000 block ≈ 25 MB; chạy đủ bảng planning là 381 394 block ≈ 95 MB.

---

## 8. Sau khi đo xong

`python\defect_from_clicklog.py` là đường **offline** dựng lại từ log click `[v14]`
— giữ lại làm đối chứng, đừng xoá. Khi có dữ liệu capture thật thì
`phase2_security_analysis.py` là đường chính, và nó tốt hơn ở ba điểm:

1. `h_mean` 12 bit không bão hoà (log cũ `irrad` 8 bit chạm trần 255 = h 2.0 ở
   21/160 block của L5).
2. `trials` là số đo thật, không phải giả định N/2.
3. Có control Y₀ = 0.

Hai đường cho cùng một xu hướng thì đó là một lập luận độc lập nữa cho bài báo.

---

## Phụ lục — bảng lệnh UART `[Q2]`

| Byte | Ý nghĩa |
|---|---|
| `0x30 \| dist[3:0]` | Chỉ số cự ly. Clear ocean: `(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80)` m |
| `0x40 \| {water[1:0], lam[1:0]}` | water 0/1/2 = clear/coastal/harbor · lam 0/1/2 = 450/532/650 nm |
| `0x50 \| turb[2:0]` | Mức nhiễu loạn, giữ cố định |
| `0x58 \| turb[2:0]` | Như trên, bật random walk ±1 (**đừng dùng cho campaign**) |
| `0x60 \| coh[3:0]` | `coh = 11` ⇒ 2¹⁶ = 65536 attempt/block. **Bắt buộc** |
| `0x70` / `0x71` | Fast mode off / on |
| `0x72` / `0x73` | Y₀ = 0 (control) / Y₀ nominal |
| `0x78 \| rate[2:0]` | 0 free · 1 20k · 2 100k · 3 1M · 4 2M · 5 5M · **6 = 10M (f_rep bài báo)** · 7 13.1M |
| `0x01` | Reset thống kê + lưới block. **Luôn gửi cuối cùng** |

Dải rate là **0x78–0x7F**, *không phải* 0x80–0x87 — 0x80–0x87 là lệnh qubit
`0x80 | (data<<2) | (abasis<<1) | bbasis`. Đặt nhầm thì board nuốt mọi lệnh qubit.
