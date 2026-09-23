# 원 데이터 출처 확정 — 2026-09-23

## 결론
로컬 `MAGs_FASTA_files/` 111개 MAG는 **Pérez-Valera E, Elhottová D. Dataset of 111 metagenome-assembled genomes from cattle manure, soil and manured soil samples. Data Brief 2025;61:111748. doi:10.1016/j.dib.2025.111748 (PMID 40599424)** 의 Zenodo 배포본과 **111/111 서열 내용 동일**(gzip 헤더 제외 서열 md5 대조, 이 폴더의 `zenodo_MAGs_data_flat.tsv`가 원 메타데이터).

- NCBI BioProject **PRJNA1231077** · Zenodo **10.5281/zenodo.15309541**(zip md5 ed4f7819…) · 관련 연구논문 Sardar et al. 2023 FEMS Microbiol Ecol 99:fiad148.
- 시료: 체코 České Budějovice, **젖소 1농장 분변(항생제 예방투여) + 유기농 2농장 토양**으로 만든 마이크로코즘, 2·14·28일차. 시료코드 GT{2,14,28}{B,S}{A,C2,EX} = 일차 × 토양 실험(B/S) × 처리.
- 🔴 **리드는 직접 메타게놈이 아니라 CHROMagar Acinetobacter 28 °C 24 h 농화배양물**의 DNA(NovaSeq 150 PE). MAG 세트는 비발효 그람음성균 편향(Pseudomonas 28·Stenotrophomonas 20·Acinetobacter 18 등 Pseudomonadota 99/111).
- 조립·비닝: metaSPAdes 개별+공동조립, MetaBAT/MaxBin/SemiBin2/COMEbin/AVAMB, dRep 95 % ANI. "ACE hybrid pipeline"은 사실이 아님.

## 🔴 MAG 접두어의 실제 의미 = 비닝 도구 (축종 아님)
| 접두어 | 원고 서술 | 실제 | n |
|---|---|---|---|
| C | cattle | COMEbin | 22 |
| M | swine | MaxBin 11 + MetaBAT 4 | 15 |
| S | sheep | SemiBin2 | 32 |
| V | poultry | AVAMB | 42 |

실제 출처 구분: manure 44 · manured soil 38 · soil 29. 여섯 후보 = S26·M1 manure(Soil S 실험) / S13·S16·S23·C22 **manured soil(Soil B 실험, GT28BC2·GT28BC2·GT2BC2·GT14BC2)**.

## 영향 범위(수정 필요)
1. 제목·Highlights·초록·키워드·§1·§2.1·§4.1·§5·윤리·Data availability: "four livestock-waste microbiomes(cattle/swine/sheep/poultry)" 전면 오류 → 단일 젖소 분변–토양 마이크로코즘 농화배양 MAG 재분석으로 재기술.
2. 축종 변수를 쓴 분석 전부 무효: §3.10 source PERMANOVA(pangenome F=1.25 P=0.11 · trait F=2.71 P=0.001 "waste source 반영") · §3.9/Fig S5 coverage-by-source · Table S2P/S2T/S3C/S3L source 열 · Fig 1 source 링·그래픽초록 source donut·Fig 5d 색 · build_v2_supS5.py PREFIX_SOURCE. → manure/manured soil/soil로 재실행하거나 삭제.
3. §3.8 "candidates came from cattle(C22), swine(M1), sheep(…)" 문장 오류.
4. MGnify 비교(cow/sheep rumen·pig/chicken gut 카탈로그 선택 근거가 4축종 전제) → 카탈로그 선택 재검토, 농화배양 편향을 한계로 명시.
5. 분류군 조성이 "livestock microbiome" 대표가 아니라 선택배지 농화물임을 §2.1·§4.4 한계에 명시.
6. Data availability: 원 accession 인용 + 파생물만 Zenodo. Ethics: 원 연구의 채취 설명으로 교체.
