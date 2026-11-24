# cryptbam 真實資料 Demo 流程

## 資料

| 檔案 | 大小 | 來源 |
|---|---|---|
| `data/NA12878.genes.bam` (+.bai) | 5.2 MB | 1000 Genomes 30× 高覆蓋 NA12878（GRCh38，NYGC），公開、不需申請。只取 3 個 gene set 的區間 ±5 kb，約 9 萬條 reads。α-globin 區正常 |
| `data/HG00684.alpha.bam` (+.bai) | 1.1 MB | 1000 Genomes 30× HG00684（華南漢人），只取 α-globin 區 chr16:130–210 kb。α 型海洋性貧血 `--SEA` 異型合子缺失帶因者 |
| `data/NA12878.BRCA1.excerpt.sam` | 1.0 MB | NA12878 BRCA1 的前 2000 條 reads，純文字 SAM（header 只留主要染色體） |
| `data/ref.fa.gz` | 1.9 MB | GRCh38 參考序列，給 bcftools 用 |

`bash scripts/fetch_demo_data.sh` 會重新下載。它用 samtools 透過 HTTP 只抓需要的 CRAM 片段，約 10 分鐘。
原始檔：`ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239334/NA12878.final.cram`、`.../ERR324/ERR3241741/HG00684.final.cram`；參考序列：1000G `GRCh38_full_analysis_set_plus_decoy_hla.fa`。

## Gene set（`beds/`；基因座標取自 Ensembl canonical transcript）

| gene set | 區間 | reads 誰能讀 | 真實深度誰能讀 |
|---|---|---|---|
| breast_cancer | BRCA1、BRCA2、PALB2 | breast_analyst | cnv_analyst |
| ovarian_cancer | BRCA1、BRCA2、RAD51C | ovarian_analyst | cnv_analyst |
| a-thalassemia | α-globin locus chr16:140–200 kb（HBZ–HBQ1，含 HBA1/HBA2） | thal_analyst | cnv_analyst |

* BRCA1、BRCA2 同時在 breast_cancer 和 ovarian_cancer，會自動歸成 `breast_cancer+ovarian_cancer` 類別，共用一把 key。
* α-globin 用整個 locus 而不是只有 HBA1/HBA2，因為 `--SEA` 缺失超出 HBA1，延伸到約 185 kb。
* 三個 gene set 都設 `depth: 10`。thal_analyst 能讀 reads（例如找點突變），但判斷 `--SEA` 這類缺失要靠深度，深度只發給 cnv_analyst。
* owner 可讀全部，包含 BED 以外的 `_background`。`outsider` 有合法格式的金鑰，但不在 policy 裡。

## 執行

在 WSL Ubuntu 裡（已建好 conda 環境 `cryptbam`，含 samtools / bcftools 1.25、pysam、官方 crypt4gh 1.8.6）：

```bash
conda activate cryptbam
cd /mnt/d/Cryptbam/cryptbam
PAUSE=1 bash scripts/demo_real.sh      # 每一步停下來等 Enter，適合現場講解；不加 PAUSE 約 35 秒跑完
```

## 流程（每一步要講的重點）

| 步驟 | 展示 | 實際結果 |
|---|---|---|
| 0. 金鑰 | 每個身分一組 X25519 金鑰 | 6 組 |
| 1. 正常 SAM | SAM → BAM → SAM，`diff` 比對 reads | 2000 條完全相同；`samtools flagstat` 正常 |
| 2. 不加密 | 任何拿到檔案的人 | 全部可讀；α-globin 區有 16,619 條 reads；`@PG CL` 暴露分析主機路徑 |
| 3. 一般 Crypt4GH | 官方 `crypt4gh encrypt`，一把 key | outsider 被擋；但 breast_analyst 解開後 α-globin 16,619 條全都看得到（全有或全無） |
| 4. cryptbam 加密 | 依 3 個 BED 切成 13 段、4 種類別 + background | 一個 `.bam.c4gh`（1.6 MB）＋一個 `.depth.bdg.c4gh`（3.6 MB） |
| 5a. 正確金鑰 | breast_analyst 解密 | 只解自己的段、其餘跳過；官方 crypt4gh 解出的位元組相同；BRCA1 有 5,414 條、α-globin 0 條 |
| 5b. 錯誤金鑰 | outsider 解密 | 本工具與官方工具都拒絕：`no data key ... addressed to this private key` |
| 5c. 竄改 | 改最後一段的 1 個 byte | `segment 24: authentication failed` |
| 6. 不同身分 | 同一個檔案，4 個身分各自解密 | 見下表 |
| 7. 深度 | 解出的 BAM 深度恆為 10×；真實深度只有 cnv_analyst 解得開 | 原始 BRCA1 平均 34.2×，breast 視角 10.0×；breast_analyst 解深度檔被拒 |
| 8. 可用性 | BRCA1 用 bcftools call，原始 vs breast 視角 | 172 個變異中保留 159 個（92%） |
| 9. 深度保護 | 海貧帶因者 HG00684 vs 正常 NA12878，同一份 policy | 原始缺失區/兩側深度比 0.57 vs 1.16，一眼看出帶因；thal_analyst 看到的兩人都是 1.00、每格都 10×；cnv_analyst 解深度檔仍看得到 0.57；thal_analyst 解深度檔被拒 |

第 6 步的輸出（reads 數 / 平均深度）：

| | alpha-globin | BRCA1 | BRCA2 | PALB2 | RAD51C |
|---|---|---|---|---|---|
| 原始 BAM | 16,619 / 36× | 21,010 / 34× | 21,325 / 33× | 9,982 / 34× | 10,756 / 33× |
| breast_analyst | 無權限 | 5,414 / 10× | 5,660 / 10× | 2,553 / 10× | 無權限 |
| ovarian_analyst | 無權限 | 5,414 / 10× | 5,660 / 10× | 無權限 | 2,877 / 10× |
| thal_analyst | 4,008 / 10× | 無權限 | 無權限 | 無權限 | 無權限 |
| owner | 4,008 / 10× | 5,414 / 10× | 5,660 / 10× | 2,553 / 10× | 2,877 / 10× |

breast 和 ovarian 的 BRCA1 指紋相同（同一份資料、同一把 key）。

## 第 9 步的輸出（每格 5 kb 平均深度，* 為 `--SEA` 缺失區）

```
                                             *    *    *    *        del/flank   N%
HG00684 plain      32   29   35   31   31   19   15   21   18   29   33   31   0.57    0  <- 看得出帶因
NA12878 plain      35   35   35   34   31   40   35   44   39   30   36   34   1.16    0
HG00684 thal       10   10   10   10   10   10   10   10   10   10   10   10   1.00    1
NA12878 thal       10   10   10   10   10   10   10   10   10   10   10   10   1.00    0
HG00684 cnv.bdg    32   29   35   31   31   19   15   21   18   29   33   31   0.57    0  <- 授權的 CNV 分析人員
NA12878 cnv.bdg    35   35   35   34   31   40   35   44   39   30   36   34   1.16    0
```

## 縮放程度 K 與變異保留（BRCA1，breast_analyst 視角）

| K | 與原始相同的 call | 備註 |
|---|---|---|
| 10 | 159 / 172 | 漏掉 1 個 SNP（chr17:43067543 C>T），其餘差異是 STR / homopolymer indel |
| 20 | 165 / 172 | SNP 回來了；多出 12 個 call，大多是重複序列的 indel；檔案約大 6% |
| 30 | 164 / 172 | 檔案約大 17% |

可在現場用 `--set-depth breast_cancer=20 --set-depth ovarian_cancer=20` 示範。

## Demo 時可以順便講的限制

1. **竄改只在「自己讀得到的段」被偵測。** 改到別人的段時，breast_analyst 根本不會解那段，所以不會報錯；owner 會偵測到。
2. **深度檔比 BAM 還大**（3.6 MB vs 1.6 MB）。真實資料每個位置深度都在變，bedGraph 沒有壓縮。
3. **header 是公開段。** `@RG` 裡的樣本名（這裡是 NA12878）每個分析人員都看得到；`@PG CL` 和 `@CO` 已經移除。
4. **`_background` 的段大小仍和深度成正比**，edit list 會讓分析人員算出這些段的大小。
5. **0× 的同型合子缺失藏不住**：深度抹平後 reads 會是 N。先前用 NA12878 的 UGT2B17（約 100 kb 同型合子缺失）實測，缺失區 95% 的鹼基是 N。`--SEA` 的同型合子（Hb Bart's）也會這樣。
6. **異型合子缺失區裡沒有異型合子 SNP**（只剩一條染色體），攻擊者仍可能從「一長段都是同型合子」推測缺失。這點尚未實測。

## 用 IGV 視覺化

先跑完 `scripts/demo_real.sh`，再執行 `bash scripts/igv_prep.sh`（在 Git Bash 或 WSL 都可以）。它會：

* 把解密出的深度檔轉成 `out/*.bedgraph`，去掉為了對齊加密區塊而補的 `#` 空白行
* 產生 `igv/*.batch`，內含本機的絕對路徑

**IGV Desktop**：基因組選 hg38，`File > Run Batch Script...` 選下列其中一個：

| batch | 位置 | 看什麼 |
|---|---|---|
| `igv/alpha_thal.batch` | chr16:140–200 kb | HG00684 原始 BAM 在 HBA1/HBA2 附近深度掉一半；thal_analyst 版兩人都是平的 10×；cnv_analyst 的 bedGraph 仍看得到缺失 |
| `igv/access.batch` | BRCA1 | 原始 34× vs breast/ovarian 的 10×（兩者相同）；thal_analyst 沒有 reads |

**瀏覽器（igv.js，不用安裝 IGV）**：

```bash
npx http-server demo -p 8765 -c-1     # 需要支援 Range 的伺服器；python -m http.server 不支援
```

打開 `http://localhost:8765/igv/`，上方兩個按鈕對應上表兩個畫面。

提醒：

* **不要直接載入 `.c4gh` 檔。** IGV 讀不了加密檔，請載入解密後的 `.bam`（`.bai` 由 demo 在解密後自行建立）。
* **cryptbam 的合成 reads 沒有 mate。** IGV 的 paired-end 視圖、insert size 著色對它無效。
