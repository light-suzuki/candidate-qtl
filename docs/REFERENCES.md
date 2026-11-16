# 根拠文献

このワークフローは、以下の既存概念を低予算の反復ゲートとして接続したものであり、新しいQTL検出法を主張するものではない。

1. Mackay TFC. *Map-based cloning of quantitative trait loci: progress and prospects.* Genet Res. 2001. doi:10.1017/S0016672301005456
   QTL map-based cloningで、表現型・背景QTL・組換え体が中心的な制約になることの背景。

2. Salvi S, Tuberosa R. *To clone or not to clone plant QTLs: present and future challenges.* Trends Plant Sci. 2005. doi:10.1016/j.tplants.2005.04.008
   QTLのfine genetic mapping、physical mapping、候補遺伝子評価の位置づけ。

3. Keurentjes JJB, et al. *QTL isogenic recombinant analysis: a method for high-resolution mapping of QTL within a single population.* Genetics. 2007. PMID:17151265.
   flanking markerでinformative recombinantを先に選び、後続の表現型・genotypingを節約する根拠。

4. Jaganathan D, Bohra A, Thudi M, Varshney RK. *Fine mapping and gene cloning in the post-NGS era: advances and prospects.* Theor Appl Genet. 2020;133:1791-1810. doi:10.1007/s00122-020-03560-w.
   post-NGS時代でもfine mapping・marker・機能検証が不可欠であることのレビュー。

5. Kamboj D, et al. *Prioritization of candidate genes for major QTLs governing yield traits employing integrated multi-omics approach in rice.* Brief Funct Genomics. 2024. doi:10.1093/bfgp/elae035.
   QTL内遺伝子をsequence variation、発現、GO、co-expression、PPIで順位付けし、対照系統の配列・発現で確認した実例。

6. Bargsten JW, Nap JP, Sanchez-Perez GF, van Dijk ADJ. *Prioritization of candidate genes in QTL regions based on associations between traits and biological processes.* BMC Plant Biol. 2014;14:330. doi:10.1186/s12870-014-0330-3.
   形質と生物学的プロセスの対応からQTL内候補を優先順位付けする方法論。

7. Sherif A, et al. *A Pseudo-near isogenic F2 population strategy for rapid QTL cloning.* Mol Breed. 2023;43:61. doi:10.1007/s11032-023-01408-x.
   F2由来の背景を整理した集団を用いるrapid QTL cloningの実例。一般F2再利用法より強い集団設計だが、目標QTL以外を抑える重要性を示す。

8. Vesztrocy AW, Dessimoz C, Redestig H. *Prioritising candidate genes causing QTL using hierarchical orthologous groups.* Bioinformatics. 2018;34:i612-i619. doi:10.1093/bioinformatics/bty615.
   orthologyを用いて、非モデル種を含むQTL区間の候補遺伝子を順位付けする手法。マイナー作物での相同遺伝子情報の使い方に直接関係する。

9. Feng C, et al. *The genomes of Mendel's peas reveal the molecular basis of seven classical traits.* Nature. 2025. doi:10.1038/s41586-025-08891-6.
   Mendel 7形質の遺伝子同定とZW6 gene IDを、V2 benchmark truthの一次根拠として使用した。

10. Foucher F, et al. *DETERMINATE and LATE FLOWERING are two TERMINAL FLOWER1/CENTRORADIALIS homologs that control two distinct phases of flowering initiation and development in pea.* Plant Cell. 2003;15:2742-2754. doi:10.1105/tpc.015701.
    DET/PsTFL1aとLF/PsTFL1cの遺伝的・配列的根拠。公開PsTFL1a proteinから誤ったV2 aliasを修正した。

11. Hofer J, et al. *UNIFOLIATA regulates leaf and flower morphogenesis in pea.* Curr Biol. 1997. doi:10.1016/S0960-9822(06)00257-0. GenBank:AF010190.2.
    UNI/PEAFLOのcosegregation、変異、公開配列とV2 answer mappingの根拠。

12. Nomura T, et al. *Brassinosteroid deficiency due to truncated steroid 5alpha-reductase causes dwarfism in the lk mutant of pea.* Plant Physiol. 2004.
    LK/PsDET2の生化学・変異根拠。公開steroid 5-alpha-reductase proteinからV2 mappingを再確認した。

13. Liew LC, et al. *EARLY FLOWERING3 redundancy fine-tunes photoperiod sensitivity.* Plant Physiol. 2017. PMID:28202598.
    HR/ELF3aとPPD/ELF3bの根拠。PPDはV2 truth mappingを十分に確定できなかったためbenchmark分母から除外した。

## LLMの扱い

LLMは上記の文献・データベース情報に基づく仮説の整理役に限定する。因果の根拠には、元論文/DB、親間配列、実験marker、F2/F2:3、機能検証のいずれかを必須とする。
