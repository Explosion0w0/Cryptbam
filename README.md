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
| 量化比較（不加密 / Crypt4GH / cryptbam） | ✅ 完成 | `scripts/benchmark.py`，結果見下方 |
| 遮蔽偵測攻擊（Mask Detection）實驗 | ❌ 未做 | 見下方「尚未量測」 |
| 效能量測（加解密時間） | ✅ 完成 | 單次執行計時，見 Performance |

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

## 量化結果：不加密 vs 一般 Crypt4GH vs cryptbam

用 `python3 cryptbam/scripts/benchmark.py`（WSL、`cryptbam` 環境）重現，結果寫在 `cryptbam/out/bench/benchmark.md`。攻擊者分兩種：

* **外部人員**：拿到檔案，沒有金鑰。
* **內部分析人員**：持有自己 gene set 的合法金鑰。一般 Crypt4GH 只有一把 key，分析人員一解密就是整份原始 BAM；cryptbam 只解得開自己的區段。

資料：模擬 30× 樣本（已知 CNV 與變異）、真實 1000 Genomes HG00684（α 型海貧帶因者）與 NA12878（正常對照）。三位分析人員（乳癌、卵巢癌、海貧）結果一致，下表列出三人中最差的值。

### Security：洩漏降低多少

| 指標 | 定義 | 不加密 | 一般 Crypt4GH | cryptbam |
|---|---|---|---|---|
| 外部人員能否解密 | 用不在名單上的金鑰解密 | 不需金鑰 | 否 | 否 |
| ② Sequence Leakage Rate | 解密後落在授權區外的比對鹼基 ÷ 原始授權區外鹼基 | 100% | 100% | **0%**（模擬、真實皆是） |
| ② Mate 位置洩漏 | 指向授權區外的 mate 座標數（真實資料） | 30,296–61,711 | 同左 | **0** |
| ② Metadata 洩漏 | header 中的指令路徑 `@PG CL`、病人備註 `@CO` | 有 | 有 | **無** |
| ① Coverage Leakage Rate | 能得知深度的授權區外位置 ÷ 原始有覆蓋的授權區外位置 | 100% | 100% | **0%** |
| ① 深度 → CNV AUROC | 1 kb 視窗深度偏離中位數的程度區分 CNV（三組 gene set 全區；0.5 = 亂猜） | 0.970 | 0.970 | **0.503** |
| ① 等位比例 → CNV AUROC | 異型合子 SNP 的 \|AF − 0.5\| 區分重複區 | 0.955 | 0.955 | 0.548–0.893（殘餘洩漏，見限制） |
| ① 海貧帶因者深度比 | `--SEA` 缺失區 ÷ 兩側（HG00684 帶因 vs NA12878 正常） | 0.57 vs 1.16 | 0.57 vs 1.16 | **1.00 vs 1.00**（海貧分析人員） |
| ① 模擬 BRCA1 缺失深度比 | 缺失區 ÷ 兩側 | 0.46 | 0.46 | **1.00** |
| ③ Size Inference AUROC（外部人員） | 只看加密檔總大小，判斷有無 BRCA1 缺失（8 個模擬樣本，4 個有缺失） | — | **1.000**（10.405 vs 10.515 MB） | **0.625** |
| ③ Size Inference AUROC（內部人員） | 海貧分析人員從 edit list 得知「被跳過的位元組數」 | — | 可直接讀全部 | K=10：**0.625**；不縮放：0.938 |
| ③ 外部人員可見的 header 封包數 | 檔頭明文的封包數量 | — | 2 | 18（透露金鑰與讀者數量） |

重點：

* 一般 Crypt4GH **連外部人員都會洩漏 CNV**：缺失少了 reads，整個加密檔小了約 110 KB，AUROC 1.0。cryptbam 的基因區段經過縮放，大小只跟區間長度有關，降到 0.625（8 個樣本下接近亂猜）。
* cryptbam 若**不縮放**（K=0），內部人員仍可從 edit list 推測缺失（AUROC 0.938），所以縮放對檔案大小保護也是必要的。
* 模擬中乳癌分析人員的深度比 AUROC 雖算出 0.875，但兩組數值是 1.00002 對 1.00002–1.00003，只差 1–2 個位置，實際上沒有可用訊號；樣本數少才排得出順序。

### Utility：仍保留多少可用資料

| 指標 | 不加密 / 一般 Crypt4GH | cryptbam（K=10） |
|---|---|---|
| 授權區內有覆蓋的位置保留率 | 100% | **100%**（模擬、真實皆是） |
| 變異 recall（模擬，對照真實答案，排除 CNV 區） | 98.2–100% | **98.2–100%**（與不加密相同） |
| 變異一致性（真實 HG00684，與原始 BAM 的 bcftools 結果比對） | — | 乳癌 343/372（92%）、卵巢癌 300/322（93%）、海貧 159/173（92%） |
| samtools / bcftools 相容 | ✓ | ✓ |
| 官方 crypt4gh 1.8.6 解出相同位元組 | ✓ | ✓ |

真實資料約 8% 的差異來自重複序列的 indel 與少數 SNP；提高 K（例如 K=20）可改善，代價是等位比例洩漏增加。

### Performance：時間與儲存成本

| 指標 | 一般 Crypt4GH | cryptbam |
|---|---|---|
| 加密後大小（真實，原始 4.48 MB） | 4.48 MB | 1.58 MB reads + 3.28 MB 深度檔 |
| 加密後大小（模擬，原始 10.40 MB） | 10.41 MB | 5.97 MB reads + 4.00 MB 深度檔 |
| 加密時間（真實） | 0.08 s | 8.56 s（縮放與深度計算佔大部分） |
| 解密時間（乳癌分析人員，真實） | 0.27 s | 0.08 s（只解自己的段） |

時間為 WSL 上單次執行的結果，僅供量級參考。

### 尚未量測

* ④ Mask Detection Precision / Recall：需先定義「遮蔽位置是否屬於秘密」；目前只確認未授權區在解密結果中完全不存在，以及授權區內真實 0× 會變成 N（NA12878 UGT2B17 缺失區 95% 為 N）。
* Buffer Zone 對鄰近區 CNV 推測的影響。
* 更多樣本的 Size Inference（目前 8 個）。

## 真實規模的成本推估（尚未實測）

以下是用 demo 的實測速度往上推的**估計值**，還沒有在真實 WGS 上跑過。

### 一份 30× WGS BAM（約 40 GB）

| 項目 | 一般 Crypt4GH | cryptbam（目前程式） | cryptbam（改良後預估） |
|---|---|---|---|
| 加密時間 | 約 18 分鐘（I/O 為主） | 約 6–7 小時（單執行緒；背景區 boundary clipping 約 3.8 小時） | 約 1 小時（12 核心，依染色體平行） |
| reads 檔大小 | 約 1.0× 原始 | 約 1.0× 原始（基因區縮小，背景區不變） | 約 1.0× |
| 深度檔大小 | 無 | 約 25 GB（未壓縮 bedGraph） | 約 1–4 GB（BGZF 壓縮或分 bin） |
| 記憶體 | 很小（串流） | 每個區段的 reads 全部載入；整條染色體的背景區需 25–50 GB RAM | 串流處理，固定記憶體 |
| 解密時間（單一分析人員） | 約 15 分鐘（整份） | 目前仍逐段讀過整份檔案 | 可直接跳過不屬於自己的段，只需數秒到數分鐘 |

改良方向：區段改為串流寫出、依染色體平行、深度檔壓縮、解密時 seek 跳過別人的段。

### 用 ClinVar 的數百到數千個疾病 gene set

擔心的是「疾病越多，重疊的組合越多，需要的 key 會排列組合式暴增」。實際上**不會**：

* cryptbam 把所有 BED 疊起來切成區段（atom），只有**實際出現在某個位置上的疾病組合**才變成一個類別、一把 key。
* 每個區間只貢獻 2 個切點，所以區段數 ≤ 2 × 區間數 + 染色體數，類別數 ≤ 區段數。1000 種疾病理論上有 2^1000 種組合，但基因組上只會出現幾千種。
* 以基因為單位時，類別數約等於基因數。例如 BRCA1 出現在 50 種疾病裡，它就是「這 50 種疾病」一個類別。

真正隨規模增加的成本：

| 成本 | 跟什麼成正比 | 說明 |
|---|---|---|
| 金鑰數量 | 類別數（≈ 基因數） | 每把 32 bytes，可忽略 |
| header 大小 | （讀者，類別）配對數 | 每個配對約 108 bytes 的封包 |
| 讀 header 的時間 | 封包總數 | 每位讀者逐一嘗試，每次一次 X25519 |
| **補齊空間** | **區段數** | 每段至少 64 KiB，reads 檔與深度檔各一份 |
| 縮放時間 | 基因區間總長度（bp） | 實測約 4.8 萬 bp/秒（K=10） |
| 每段固定開銷 | 區段數 | 暫存 BAM、深度計算（估每段約 10 ms） |

推估（假設每種疾病平均 20 個基因、每個基因約 60 kb；30× WGS）：

| 情境 | 區間數 | 區段數 | 類別（key） | header | 補齊（每個檔） | 縮放時間（單執行緒） |
|---|---|---|---|---|---|---|
| 目前 demo（3 組、5 個基因） | 5 | 13 | 5 | 2 KB | 0.8 MB | 7 秒 |
| 100 種疾病、500 個基因 | 500 | 約 1,000 | 約 500 | 約 0.2 MB | 約 64 MB | 約 10 分鐘 |
| 1,000+ 種疾病、5,000 個基因（基因層級） | 5,000 | 約 10,000 | 約 5,000 | 約 2 MB | **約 650 MB** | **約 1.7 小時**（約 300 Mb） |
| 同上，改用 exon 層級 | 約 10 萬 | 約 20 萬 | 約 5,000 | 約 2 MB | **約 13 GB** | 約 5 分鐘（約 15 Mb） |

結論：

* **基因層級**：檔案多約 1.3 GB（約 3%），加密多 1–2 小時，可以接受。
* **exon 層級**：補齊空間與每段固定開銷暴增，是最該避免的設計；BAM 越小（例如 panel）越嚴重。
* 區段越多，edit list 透露的大小資訊越細（背景區可逐段推測 CNV），header 封包數也透露讀者規模。

降低成本的做法：以基因而非 exon 為單位；合併相近且同一組疾病的基因；背景區小段不另外切或一起縮放；用 HKDF 從主金鑰衍生各類 key；依「讀者角色」而非疾病分組（header 取決於讀者數 × 權限數）；依染色體平行處理。

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
| 檔案大小洩漏 | 讀者可從 edit list 算出別人區段大小；`_background` 大小與深度成正比；header 封包數透露讀者數量 | 背景區也縮放，或補齊到固定大小 |
| 等位基因比例洩漏 CNV | K=10 時 AUROC 0.864 | `quantize_af`，但不適合腫瘤低 VAF |
| 授權區真實 0× 會變成 N | 同型合子缺失仍看得出 | 尚無好解法；填參考序列會捏造結果 |
| 異型合子缺失區沒有異型合子 SNP | 可能從「一長段同型合子」推測缺失 | 尚未實測 |
| 縮放會漏掉部分變異 | K=10 漏 1 個 SNP 與數個 STR indel | 提高 K（K=20 找回該 SNP，檔案大 6%） |
| read pair 被拆開 | 結構變異 / 斷點證據消失 | 需要時該區設 `depth: 0` + `softclip: inside` |
| 竄改只在自己讀得到的段被偵測 | 別人的段被改，該分析人員不會發現 | owner 解密時可偵測 |
| header 是公開段 | `@RG` 樣本名所有人都看得到 | 可另行去識別化 |
| `owner_manifest.json` 是明文 | 含每段平均真實深度（等於明文 CNV 資訊） | 改成加密或可選輸出 |
| 深度檔未壓縮 | 3.3 MB，比 reads 檔大 | 改用 BGZF 壓縮的 bedGraph |
| 每個區段的 reads 全部載入記憶體 | 真實 WGS 的整條染色體背景區需 25–50 GB RAM | 改成串流寫出 |
| `.bai` 不能明文發布 | 會透露每段資料量 | 讀者解密後自行建立 |
| samtools 直接讀 `.c4gh`（htslib-crypt4gh plugin） | 未驗證 | 目前做法是先解密再讀 |
| HG00684 的 BRCA1 變異未比對 ClinVar | 無法展示「找到致病變異」 | 比對 ClinVar；「乳癌病人」僅為情境設定 |

---

## 下一步

- [ ] Size Inference 擴大樣本數（目前 8 個模擬樣本），並改用多個 1000 Genomes 真實樣本
- [ ] 補做 Mask Detection Attack，並明確定義遮蔽位置是否為秘密
- [ ] 實測「異型合子 SNP 消失」能否推測缺失
- [ ] 決定 demo 預設 K（目前 10）
- [ ] `owner_manifest.json` 改為加密保存
- [ ] 在 IGV 中確認視覺化畫面
- [ ] 向 crypt4gh 上游回報 64 KiB 倍數的 bug
- [ ] 加密改為串流 + 依染色體平行，深度檔壓縮，解密時跳過別人的段（WGS 規模必要）
- [ ] 下載 ClinVar 基因與疾病對照表，實際算出區段數、類別數與 header 大小

資料來源：1000 Genomes Project 30× 高覆蓋 GRCh38（NYGC），樣本 NA12878、HG00684；HG00684 的 `--SEA` 帶因狀態依一篇預印本的補充表，並以本專案的深度量測確認。
