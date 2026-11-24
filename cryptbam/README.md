# cryptbam：NGS BAM 分區加密（相容 Crypt4GH）PoC

一個 `.bam` 會產生一個 `.bam.c4gh`。每個基因區間用不同的 data key 加密，data key 再用各分析人員的 public key 包進 header。
分析人員解密後，得到的是**只含自己區間**的合法 BAM，可以直接用 samtools / bcftools。
真實深度另外存成 `.depth.bdg.c4gh`，用另一組 key 加密。官方 `crypt4gh decrypt` 解出的位元組與本工具完全相同。

## 使用

```bash
pip install pysam cryptography numpy           # 評估另需 scikit-learn、samtools、bcftools
python3 scripts/simulate.py --outdir data       # 模擬資料（含已知 SNP / indel / CNV）
python3 -m cryptbam keygen breast_analyst --outdir data/keys
python3 -m cryptbam encrypt --bam data/sample.bam --policy data/policy.json --out out/sample \
        [--depth 10] [--set-depth breast_cancer=20 --set-depth a-thalassemia=0]
python3 -m cryptbam decrypt --sk data/keys/breast_analyst.sec --in out/sample.bam.c4gh --out breast.bam
crypt4gh decrypt --sk data/keys/breast_analyst.sec < out/sample.bam.c4gh > breast.bam   # 官方工具，結果相同
python3 -m cryptbam inspect --sk data/keys/thal_analyst.sec --in out/sample.bam.c4gh      # 稽核：這把 key 能開哪些 segment
bash scripts/run_demo.sh        # 一鍵展示
python3 scripts/compare.py      # 比較：不加密 / 一般 Crypt4GH / cryptbam，攻擊者能拿到什麼
python3 scripts/evaluate.py     # 掃描縮放程度 → out/results.md
```

### 真實資料 demo（NA12878，GRCh38）

```bash
bash scripts/fetch_demo_data.sh     # breast_cancer / ovarian_cancer / a-thalassemia 區間 ±5 kb + 海貧帶因者 HG00684，約 8 MB
PAUSE=1 bash scripts/demo_real.sh   # SAM↔BAM、三種做法比較、正確/錯誤金鑰、各身分差異、海貧帶因者的深度保護
bash scripts/igv_prep.sh            # 產生 IGV batch 檔與 igv.js 頁面用的 bedGraph
```

流程與實測結果見 `demo/DEMO.md`。

## 縮放程度由使用者決定

`policy.json` 中每個 gene set 各自設定，CLI 可以臨時覆寫：

```json
"gene_sets": {
  "breast_cancer":  {"bed": "beds/breast_cancer.bed",  "depth": 10},
  "ovarian_cancer": {"bed": "beds/ovarian_cancer.bed", "depth": 10, "quantize_af": true},
  "a-thalassemia":  {"bed": "beds/a-thalassemia.bed",  "depth": 20, "softclip": "inside"}
}
```

| 參數 | 意義 |
|---|---|
| `depth: 0` | 不縮放，只做跨界切割與清理（保留真實 reads 與深度） |
| `depth: K` | 依每個位置的鹼基比例重新產生剛好 K 條 reads，深度處處等於 K；等位基因比例的解析度為 1/K |
| `quantize_af` | 再把等位基因比例四捨五入成 0 / 0.5 / 1 |
| `softclip` | `drop`（預設）全部移除；`inside` 只保留投影仍在區間內的 soft clip |
| `--depth` | policy 沒寫 `depth` 時的預設值 |
| `--set-depth GS=K` | 執行時覆寫某個 gene set |

同一段同時屬於多個 gene set、設定又不同時，採**最保守**的設定：取最小的非 0 K、任一方要求 quantize 就 quantize、soft clip 必須所有 gene set 都設 `inside` 才保留。

## 結果（模擬 30×；全部 gene set 用同一設定）

| 設定 | 檔案 MB | 區外序列洩漏 | 深度→CNV AUROC | 等位比例→CNV AUROC | 變異 recall | 官方 CLI 相同 | samtools OK |
|---|---|---|---|---|---|---|---|
| raw (0) | 10.30 | 0% | 0.959 | 0.950 | 99.4% | ✓ | ✓ |
| K=5 | 5.90 | 0% | 0.494 | 0.642 | 99.4% | ✓ | ✓ |
| K=10 | 5.97 | 0% | 0.494 | 0.864 | 99.4% | ✓ | ✓ |
| K=20 | 6.17 | 0% | 0.494 | 0.942 | 99.4% | ✓ | ✓ |
| K=30 | 6.49 | 0% | 0.494 | 0.925 | 99.4% | ✓ | ✓ |
| K=10 + quantize | 5.97 | 0% | 0.494 | 0.500 | 99.4% | ✓ | ✓ |

AUROC 0.5 代表攻擊者只能亂猜。只要開了縮放，深度就不再洩漏 CNV。K 越大，等位基因比例越精細：對腫瘤低 VAF 變異越有用，但從等位比例推 CNV 也越準。這就是要讓使用者自己取捨的地方。

## 設計重點（review 指南）

| 檔案 | 行數 | 內容 |
|---|---|---|
| `c4gh.py` | ~185 | Crypt4GH：金鑰、header 封包、多 key writer、依 edit list 解密 |
| `bgzf.py` | ~50 | BGZF 壓縮，以及「精確大小」的空 BGZF padding |
| `regions.py` | ~55 | BED → atom → 存取類別（處理重疊） |
| `transform.py` | ~220 | 跨界切割、pileup 壓縮（縮放）、真實深度 |
| `pipeline.py` | ~120 | 明文版面與加密流程、policy |
| `cli.py` | ~50 | 命令列 |

* **多 key**：Crypt4GH header 可以放多把 data key。每個讀者只打得開寄給自己的封包，也就只拿到自己類別的 key，外加一份 edit list，用來跳過沒有權限的 run。
* **重疊基因**：三個 BED 疊起來切成 atom，覆蓋集合相同的 atom 共用一把 key。例如 BRCA1 會變成類別 `breast_cancer+ovarian_cancer`，只存一份。
* **跨界 read**：在邊界切開；移除 soft clip、mate 資訊和 MD/NM/MC/SA 等 tag；header 只留 `@HD @SQ @RG`。
* **samtools 相容**：每個 run 是獨立的 BGZF block，並補齊到 64 KiB segment，任意子集合串接後仍是合法 BAM。

## 已知問題

1. 官方 crypt4gh ≤1.8.6 有 bug：要跳過的長度剛好是 64 KiB 的倍數時，會去解密被跳過的最後一段而失敗。我們在每個 run 開頭放 1 個所有人都會跳過的 byte 來繞過。
2. htslib 要求 BGZF 的 XLEN 等於 6，padding 不能放在 gzip EXTRA 欄位，否則 index 查詢會壞掉。改成用只含空 deflate block 的 BGZF block 補齊。
3. **檔案大小仍會洩漏**：讀者從 edit list 算得出別人 run 的大小。raw 和 `_background` 的大小跟深度成正比；有縮放的區間，大小只跟長度有關。
4. 縮放後，等位基因比例仍會洩漏 CNV（見上表）。`quantize_af` 可以藏住，但無法分析腫瘤的低 VAF 變異。
5. 0× 的位置會輸出 N：深度看不出來，但一整段 N 本身就透露「沒資料」。填參考序列會捏造假結果，所以不採用。
6. 切割後 read pair 變成 single-end，SV / 斷點證據會消失。需要看斷點或重複序列擴增的區間，請用 `depth: 0` 加 `softclip: inside`。
7. `.bai` 索引不可明文發佈，它會透露每段的資料量；讀者應在解密後自行建立。
8. 只支援未加密的私鑰（`crypt4gh-keygen --nocrypt` 或本工具的 `keygen`）。

模擬資料採 hg38 座標，但序列是隨機產生的。三個 gene set 為 breast_cancer、ovarian_cancer、a-thalassemia；已知變異包含 BRCA1 異型合子缺失、BRCA2 重複、α-globin 的 `--SEA` 異型合子缺失（海貧帶因）與 HBZ 5 bp 缺失。BED 中 BRCA2 exon10-11 與 400 bp 小區間是測試用的虛構區間。
