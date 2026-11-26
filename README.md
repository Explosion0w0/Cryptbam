# cryptbam：NGS BAM 分區加密

把比對好的 `.bam` 依基因區間分開加密，格式相容 GA4GH Crypt4GH。同一個加密檔，每位分析人員只解得開自己被授權的區間；深度縮放讓攻擊者無法從定序深度推測拷貝數變異（CNV）。

| 資料夾 | 內容 |
|---|---|
| [`cryptbam/`](cryptbam/) | 工具本體（`python -m cryptbam keygen / encrypt / decrypt / inspect`）、模擬與評估腳本、真實資料 demo。細節見 [cryptbam/README.md](cryptbam/README.md)、[cryptbam/demo/DEMO.md](cryptbam/demo/DEMO.md) |
| [`cryptbam_demo/`](cryptbam_demo/) | 精簡 demo：一位病人的 BAM（1000 Genomes HG00684，公開資料）、參考序列、三個風險基因 BED，照 [`demo.txt`](cryptbam_demo/demo.txt) 一步步執行 |

> 私鑰（`*.sec`、`keys/`）與輸出檔不在 repo 中；`cryptbam_demo/scripts/keygen.sh` 會產生 demo 用金鑰。

---

## 快速開始

需要 Linux 或 WSL（pysam 不支援 Windows）。

```bash
conda create -y -n cryptbam -c conda-forge -c bioconda python=3.11 pysam samtools bcftools numpy cryptography scikit-learn pip
conda activate cryptbam
pip install crypt4gh                      # 官方工具，用於相容性比對（可選）
git clone https://github.com/Explosion0w0/Cryptbam.git
cd Cryptbam && export PYTHONPATH=$PWD/cryptbam
cd cryptbam_demo && bash scripts/keygen.sh    # 接著照 demo.txt 執行
```

---

## 目前進度

| 項目 | 狀態 | 說明 |
|---|---|---|
| 分區加密（多把 data key + edit list） | ✅ 完成 | 依 BED 切段，每種存取類別一把 key，再用各讀者的 public key 包進 header |
| 重疊基因處理 | ✅ 完成 | BRCA1/BRCA2 同屬乳癌與卵巢癌，資料只存一份，key 發給兩種分析人員 |
| 跨界 read 切割與 metadata 清理 | ✅ 完成 | 在邊界切開；移除 soft clip、mate 資訊、MD/NM/MC/SA 標籤、`@PG CL`、`@CO` |
| 深度縮放（K 由使用者決定） | ✅ 完成 | 每個 gene set 可設 `depth`、`quantize_af`、`softclip` |
| 真實深度另存 | ✅ 完成 | `.depth.bdg.c4gh` 用另一組 key，只發給 CNV 分析人員 |
| Crypt4GH 相容 | ✅ 完成 | 官方 `crypt4gh` 1.8.6 解出的位元組與本工具完全相同；samtools / bcftools 可直接讀 |
| 模擬資料評估 | ✅ 完成 | `scripts/evaluate.py`：已知 CNV 與變異，掃描不同 K |
| 真實資料 demo | ✅ 完成 | NA12878、HG00684（α 型海貧帶因者），`scripts/demo_real.sh` |
| 精簡 demo 資料夾 | ✅ 完成 | `cryptbam_demo/` + `demo.txt` + `view_range.sh` |
| IGV 視覺化 | ⚠️ 部分 | batch 檔與 igv.js 頁面已產生，尚未在 IGV 中實際確認畫面 |
| 檔案大小保護的量化實驗 | ❌ 未做 | 見下方指標 ③ |
| 遮蔽偵測攻擊（Mask Detection）實驗 | ❌ 未做 | 見下方指標 ④ |
| 效能量測（加解密時間） | ❌ 未做 | 只有整體 demo 時間，沒有分項計時 |

### Demo 的 gene set

| gene set | 區間（GRCh38） | 讀取權限 |
|---|---|---|
| breast_cancer | BRCA1、BRCA2、PALB2 | breast_analyst |
| ovarian_cancer | BRCA1、BRCA2、RAD51C | ovarian_analyst |
| a-thalassemia | α-globin 區 chr16:140–200 kb（HBZ–HBQ1） | thal_analyst |

三組的真實深度只發給 cnv_analyst；owner 可讀全部。原本的第三組是亨丁頓舞蹈症（HTT），因為公開資料中找不到 CAG 擴增樣本（都需申請），改成 α 型海貧：HG00684 有真實的 `--SEA` 缺失，也貼近台灣的臨床情境。

---

## 四項資安工程與做法

| 工程 | 要解決的問題 | cryptbam 的做法 | 狀態 |
|---|---|---|---|
| ① Coverage Protection | 定序深度洩漏 CNV | 授權區內重新產生每個位置剛好 K 條 reads，深度處處相同；真實深度另存、另一組 key | ✅ |
| ② Boundary Read Protection | 跨界 read 與 mate 資訊洩漏區外序列 | 在邊界切開 read；清掉 mate 位置、soft clip、MD/NM/MC/SA；header 只留 `@HD @SQ @RG` | ✅ |
| ③ File Size Protection | 從檔案或區段大小推測資料量 | 每段補齊到 64 KiB 加密區塊；縮放後區段大小只跟區間長度有關 | ⚠️ 背景區未處理 |
| ④ Zero-Coverage Protection | 0× 暴露被遮蔽的區域 | 未授權區域不是填 0，而是根本不存在於解密結果中（沒有 key 就不解密） | ⚠️ 授權區內的真實 0× 仍會變成 N |

---

## 成果報告可用的量化指標

報告建議分成三面向：**Security（洩漏降低多少）**、**Utility（仍保留多少可用資料）**、**Performance（時間與儲存成本）**。下表的「目前結果」都是已實際跑出的數字；「—」表示尚未量測。

### ① Coverage Protection（定序深度保護）

| 指標 | 定義 | 目前結果 | 重現方式 |
|---|---|---|---|
| 深度 → CNV AUROC | 以每 1 kb 視窗的 \|log2(深度 / 中位數)\| 區分 CNV 與正常視窗；0.5 = 亂猜 | 不縮放 0.959 → K=10 **0.494** | `scripts/evaluate.py`（模擬資料） |
| 等位基因比例 → CNV AUROC | 以異型合子 SNP 的 \|AF − 0.5\| 區分重複區與其他區（殘餘洩漏） | 不縮放 0.950；K=5 0.642；K=10 0.864；K=20 0.942；K=10 + quantize **0.500** | 同上 |
| 缺失區 / 兩側深度比（真實資料） | 海貧 `--SEA` 缺失區平均深度 ÷ 兩側平均深度 | 帶因者 HG00684 原始 **0.57**、正常 NA12878 1.16 → 解密後兩人都是 **1.00** | `cryptbam_demo/demo.txt` 或 `scripts/depth_bins.py` |
| BRCA1 缺失深度比（模擬） | 模擬的 BRCA1 異型合子缺失區 ÷ 兩側 | 不加密 / 一般 Crypt4GH 0.45 → cryptbam **1.00** | `scripts/compare.py` |
| Coverage Leakage Rate | 攻擊者取得的未授權位置深度數 ÷ 未授權位置總數 | 設計上為 0（深度檔每類一把 key，未授權類別解不開）；尚未以實驗統計 | 待補：對每位讀者解深度檔並計數 |
| Buffer Zone 效果 | 授權區外加緩衝區後，從授權區深度推測鄰近區 CNV 的 AUROC | — | 尚未實作 |

### ② Boundary Read Protection（跨界讀序保護）

| 指標 | 定義 | 目前結果 | 重現方式 |
|---|---|---|---|
| Sequence Leakage Rate | 洩漏的未授權比對鹼基數 ÷ 原始未授權比對鹼基總數 × 100% | 所有 K 設定下皆 **0%**；早期「不切割」基線為 0.048% | `scripts/evaluate.py` |
| 可讀的未授權鹼基數（乳癌分析人員） | 解密後落在授權區外的比對鹼基 | 不加密 / 一般 Crypt4GH 8,531,285 → cryptbam **0** | `scripts/compare.py` |
| Mate 位置洩漏 | 指向授權區外的 mate 座標數 | 56,892 → **0**（早期基線另有 248 筆） | `scripts/compare.py` |
| Metadata 洩漏 | header 中的指令路徑、病人備註 | 有 → **已移除**（`@PG CL`、`@CO`） | `scripts/compare.py` |
| Authorized Data Retention | 解密後保留的授權區比對鹼基 ÷ 原始授權區比對鹼基 | — | 待補；可用變異 recall 當替代指標（見 Utility） |

### ③ File Size Protection（檔案大小保護）

| 指標 | 定義 | 目前結果 | 重現方式 |
|---|---|---|---|
| Size Inference AUROC | 以區段大小預測該區屬於高或低資料量（或有無 CNV）；比較 padding 前後 | — | 待補：需多個樣本，可用 1000 Genomes 多人資料 |
| 區段大小與深度的相關 | 縮放 / 不縮放時，區段位元組數與真實深度的相關係數 | 定性：縮放後區段大小只跟長度有關；不縮放與 `_background` 與深度成正比 | 待補 |
| 外部人員可見資訊 | 沒有金鑰時能得到什麼 | 只有檔案總大小與 header 封包數 | `scripts/compare.py` |

### ④ Zero-Coverage Protection（零覆蓋度保護）

| 指標 | 定義 | 目前結果 | 重現方式 |
|---|---|---|---|
| 未授權區的可見性 | 解密結果中未授權區是否出現 0× 或任何 reads | 完全不存在（reads = 0，深度檔該區解不開），不會出現「人工 0×」 | `cryptbam_demo/scripts/view_range.sh` |
| 授權區真實 0× 的 N 比例 | 授權區內同型合子缺失，縮放後 reads 中 N 的比例 | UGT2B17（NA12878）缺失區 **95% 為 N**：仍可被偵測 | `scripts/depth_bins.py` |
| Mask Detection Precision / Recall / AUROC | 攻擊者從輸出判斷哪些區域被刻意遮蔽 | — | 待補；需先定義「遮蔽位置是否屬於秘密」 |

### Utility（可用性）

| 指標 | 目前結果 |
|---|---|
| 變異 recall（模擬，排除 CNV 區） | 所有 K 皆 **99.4%** |
| BRCA1 變異保留（真實 NA12878，bcftools） | 172 個中 K=10 保留 **159**（92%）；K=20 保留 165 |
| 工具相容性 | samtools quickcheck / index / view、bcftools 皆可直接使用 |
| 官方 Crypt4GH 相容 | 官方 1.8.6 解出的位元組與本工具 **完全相同** |

### Performance（成本）

| 指標 | 目前結果 |
|---|---|
| 檔案大小（模擬 30×） | 不縮放 10.30 MB；K=10 5.97 MB；K=30 6.49 MB |
| 檔案大小（真實 HG00684） | 原始 BAM 4.5 MB → `.bam.c4gh` 1.6 MB + `.depth.bdg.c4gh` 3.3 MB |
| 執行時間 | 完整真實資料 demo（9 步，含加解密、bcftools）約 35 秒；加密與解密的分項時間 — |

---

## 遇到的問題

### 已解決

| 問題 | 解法 |
|---|---|
| 官方 crypt4gh ≤ 1.8.6 在「跳過長度剛好是 64 KiB 倍數」時解密失敗 | 每個區段開頭放 1 個所有人都跳過的 byte（可回報上游） |
| htslib 要求 BGZF 的 XLEN = 6，padding 放在 gzip EXTRA 會讓索引查詢壞掉 | 改用只含空 deflate block 的 BGZF block 補齊 |
| 真實 GRCh38 header 有 3,366 條 contig，每條都切成至少 64 KiB 的段 | 只處理有 reads 或出現在 BED 的 contig |
| 計算真實深度時一次配置整條染色體的陣列（chr1 約 2 GB） | 改成每 1 Mb 分段計算 |
| BED 的 contig 不在 BAM 裡（`17` vs `chr17`）就崩潰 | 跳過並警告 |
| pysam 不支援 Windows | 在 WSL Ubuntu 的 conda 環境執行 |
| 公開資料找不到 HTT CAG 擴增樣本 | 第三組 gene set 改為 α 型海貧 |
| Windows 換行會讓 `.sh` 在 WSL 壞掉 | `.gitattributes` 強制 LF |

### 尚未解決 / 已知限制

| 限制 | 影響 | 可能方向 |
|---|---|---|
| 檔案大小洩漏 | 讀者可從 edit list 算出別人區段大小；`_background` 大小與深度成正比 | 背景區也縮放，或補齊到固定大小 |
| 等位基因比例洩漏 CNV | K=10 時 AUROC 0.864 | `quantize_af`，但不適合腫瘤低 VAF |
| 授權區真實 0× 會變成 N | 同型合子缺失仍看得出 | 尚無好解法；填參考序列會捏造結果 |
| 異型合子缺失區沒有異型合子 SNP | 可能從「一長段同型合子」推測缺失 | 尚未實測 |
| 縮放會漏掉部分變異 | K=10 漏 1 個 SNP 與數個 STR indel | 提高 K（K=20 找回該 SNP，檔案大 6%） |
| read pair 被拆開 | 結構變異 / 斷點證據消失 | 需要時該區設 `depth: 0` + `softclip: inside` |
| 竄改只在自己讀得到的段被偵測 | 別人的段被改，該分析人員不會發現 | owner 解密時可偵測 |
| header 是公開段 | `@RG` 樣本名所有人都看得到 | 可另行去識別化 |
| `owner_manifest.json` 是明文 | 含每段平均真實深度（等於明文 CNV 資訊） | 改成加密或可選輸出 |
| 深度檔未壓縮 | 3.3 MB，比 reads 檔大 | 改用 BGZF 壓縮的 bedGraph |
| `.bai` 不能明文發布 | 會透露每段資料量 | 讀者解密後自行建立 |
| samtools 直接讀 `.c4gh`（htslib-crypt4gh plugin） | 未驗證 | 目前做法是先解密再讀 |
| HG00684 的 BRCA1 變異未比對 ClinVar | 無法展示「找到致病變異」 | 比對 ClinVar；「乳癌病人」僅為情境設定 |

---

## 下一步

- [ ] 補做 Size Inference AUROC（多個 1000 Genomes 樣本）
- [ ] 補做 Mask Detection Attack，並明確定義遮蔽位置是否為秘密
- [ ] 量測 Authorized Data Retention 與 Coverage Leakage Rate
- [ ] 加密 / 解密分項計時
- [ ] 實測「異型合子 SNP 消失」能否推測缺失
- [ ] 決定 demo 預設 K（目前 10）
- [ ] `owner_manifest.json` 改為加密保存
- [ ] 在 IGV 中確認視覺化畫面
- [ ] 向 crypt4gh 上游回報 64 KiB 倍數的 bug

資料來源：1000 Genomes Project 30× 高覆蓋 GRCh38（NYGC），樣本 NA12878、HG00684；HG00684 的 `--SEA` 帶因狀態依一篇預印本的補充表，並以本專案的深度量測確認。
