# 저자 개정본(Manuscript_revised.docx) 이식 보고서 — 2026-09-23

정본은 계속 `SUBMISSION_v2/01_Manuscript.md`(→ `.docx`는 `rebuild_docx.sh`)이다. 이 날 저자 개정본의 **본문·저자·표 번호 체계**를 md 체인에 이식하고, 09-19에 확정한 참고문헌 교정을 다시 적용한 뒤 감사 2종을 통과시켰다.

## 0. 두 외부 파일의 정체

| 파일 | 정체 |
|---|---|
| `submission_v3/Manuscript.docx` (09-18 서버 복사) | 작성자 "정연 안", 09-16 생성·09-17 00:01 저장. v2 09-16 시점을 바탕으로 406문장 중 339문장을 재작성(제목 단축·초록 198단어·하이라이트 5개·저자란 한글). 그림 10종은 v2와 md5 동일. 보충 워크북은 v2의 S1/S2/S3을 **S3/S1/S2로 재번호**하고 시트명을 `S1.1`→`S1A`식으로 변경(내용은 셀 단위 동일). |
| `SUBMISSION_v2/Manuscript_revised.docx` (09-23 업로드) | 위 파일을 사용자가 09-17 15:17 저장한 판. 본문 동일(405/406문장), 차이는 **저자 5인 영문 표기 + 소속 3곳 주소**뿐. |

둘 다 09-19 참고문헌 교정 **이전** 판이라 Vancouver 63편·실재하지 않는 9편(Cheng 2017, Krawczyk 2021, Omoregie 2022, Li 2021, Zamanzadeh 2023, Suzuki 2022 등)이 그대로였다.

## 1. 이식 (`port_revised_260923.py`, 멱등)

* `Manuscript_revised.docx` → pandoc → `01_Manuscript.md` / `02_Figure_legends.md`. 수정 전 파일은 `*.pre_revised_260923.*`로 보관.
* 개정본에서 가져온 것: 제목·저자·소속·하이라이트·초록·키워드·본문 전체·범례 15종·Ethics/Competing interests/Funding/Data availability.
* 이전 정본에서 유지한 것: Table 1·2의 표 본체(제목·캡션은 개정본), 참고문헌 기계(`references/`).
* 추가한 것: 절 번호(감사·참고문헌 스크립트가 `## 1. Introduction`·`## CRediT`를 기준으로 동작) · CRediT 자리표시 · Elsevier GenAI 선언 절 · 저널 배너는 HTML 주석으로(docx에 나오지 않음, 09-11 미충족 항목 20 해소) · Methods 2.2에 괄호 하나: *"(labelled "MICP-complete" in the figures and tables)"* — 그림 10종과 Table 2가 여전히 "MICP-complete"로 표기하는데 개정 본문은 "prioritized candidates"만 써서 연결이 끊겨 있었음.
* 본문 4,320단어(v2 7,005), 초록 198단어(MicRes ≤250 충족), 절 구조 1–5 유지(Conclusions 있음).

## 2. 참고문헌 재적용 (`references/04_apply_to_master.py`, 09-23 항목 추가)

| 개정본 인용 | 처리 | 근거 |
|---|---|---|
| Cheng et al., 2017 | → Kumari et al., 2016 | 09-19 확정 |
| Krawczyk et al., 2021 | → Lapierre et al., 2020 | 09-19 확정. ⚠ 개정 문장 "MICP still relies largely on *S. pasteurii* and a narrow range of *Bacillus*/*Lysinibacillus* strains, underscoring the need for alternatives…"에서 Lapierre가 뒷받침하는 것은 *S. pasteurii*의 복합배지 의존성 부분. 앞부분 일반 진술에 리뷰(Achal and Mukherjee, 2015 등) 추가 인용을 권장 |
| Omoregie et al., 2022 | → Omoregie et al., 2019 | 실재본 |
| Dhami et al., 2014 (2곳) | → 2013 | 실재 연도 |
| Li et al., 2021 | → Chen et al., 2021 | 09-19 확정 |
| Zamanzadeh et al., 2023 | 인용 삭제; 대체본 Gupta 2016은 Discussion 4.1의 "livestock slurry is a reservoir of useful microbial functions" 문장에 부착 | Gupta는 이 문장을 뒷받침하고, 원래 자리(비교유전체 방법론 문장)와는 무관 |
| GTDB-Tk / dbCAN / yn00 / skani | Parks 2020 · Zheng 2023 · Yang and Nielsen 2000 · Shaw and Yu 2023 재부착 | 09-19 저자 결정(도구 인용) — 개정본이 누락 |
| Steinegger and Söding, 2017 | 이미 라이브러리에 있음(중복 추가 시도 후 철회) | — |

개정본이 더 이상 인용하지 않는 실재 문헌 3편(Hommel 2015, Lagesen 2007, Yuan 2015)은 목록에서 빠짐 → **최종 57편**, 전부 DOI 검증본, Harvard 알파벳순. EndNote 패키지(`references/01_Manuscript_EndNote.docx`, 임시인용 37군데·61건, `.ris` 60편) 재생성.

## 3. 보충표 번호 체계

`Supplementary_tables/`를 개정본 체계로 교체(내용은 셀 단위 동일, README 시트만 새 시트명):

| 새 이름 | 이전(v2) | 시트 |
|---|---|---|
| Table S1 reference panels and methods | Table S3 | S1A–S1O (15) |
| Table S2 per-MAG measurements | Table S1 | S2A–S2T (20) |
| Table S3 comparative statistics | Table S2 | S3A–S3M (13) |

이전 워크북은 `_superseded_260923/`. Table_S5b(DRAM)·Table_S7f(iqtree)는 그대로("provided separately"). Main Table 1·2 xlsx는 개정 제목으로 재빌드.

## 4. 감사

* `audit_consistency.py` — 정본 폴더의 워크북·`S1A` 시트 체계·`Table S2A` 콜아웃·저자 주소를 인용으로 오인하지 않도록 개정. **20항목 전부 PASS**.
* `audit_numbers.py` — 개정 문구 4건 교체 + 신규 검사 19건(AAI 93.15/93.49, RefSeq ANI 94.57/93.85/98.96/99.16, 최근연 종명, 18.9 %, RF 0.58, SH *P* < 0.001, yn00 *ureG* 0.31/0.074/7.7 × 10⁻⁸, *ureA/B/C* 비유의, gRodon *P* = 0.58, 7-유전자 단일 contig 26 MAG[정의 감사 규칙 C]). **109항목 전부 일치**.
* 실행: `UPCYCLING_MAN_DIR=/data/data/Upcycling/SUBMISSION_v2 python consolidation_260904/audit_{consistency,numbers}.py`

검증 못 한 수치(출하 표에 원자료 없음, v2에서 그대로 승계): DRAM "433,595 gene annotations across 98 modules"(98은 S5a 열 수와 일치), 형질 PERMANOVA pseudo-*F* = 2.71, 게놈 크기 1.9–6.2 Mb·GC 32–68 %.

## 5. 주장 강도 대조 (v2 09-19 정본 → 개정본)

### 5a. 적절히 완화된 것 (그대로 두면 됨)
* "statistically validated target list" → "prioritizes candidates for cultivation…"; 초록의 "engineering-friendly chassis"·4.3의 "suitable as non-pathogenic biocement chassis under BSL-1"·"defence-naive" 삭제 → "do not establish biological safety".
* antiSMASH "strongly enriched, 23-fold, *P* = 5.3 × 10⁻¹⁰" → 탐색적·미보정 패턴으로 서술.
* MGnify "≈ 30-fold enrichment over the global background" 삭제(분모가 다른 비교였음).
* 4.1 "environmental *Sphingobacterium* is typically slow-growing and outcompeted…"(무출처) 삭제; Intro "1.5 × 10⁹ t yr⁻¹"(무출처) 삭제; "to our knowledge, new for MICP bioprospecting" 삭제.
* CAZyme: "establishes … rather than an annotation artefact" → "supports a carbohydrate-utilization signal".
* 3.4 UreC: "supporting retention of catalytic activity" → "consistent with conserved urease catalytic structure"; 4.2에 "not evidence of enzymatic activity" 명시.

### 5b. 검토가 필요한 변경
1. **제목·초록 "convergent / functionally convergent lineages"** — 근거는 유전체 수준(operon 구성·활성부위·fold)이며 기능 측정은 없음. v2도 같은 표현을 썼으므로 유지 가능하나, "genomically convergent"가 더 정확. 저자 판단.
2. **2.3 신종 판정 기준 변경** — v2: GTDB 종 미배정 또는 ANI < 95 % **AND** 패널 내 AAI < 95 %. 개정: 동속 참조 전부 ANI < 95 %, AAI는 보조. 결과(3.1·3.8)와는 정합. 실제 분석 절차와 맞는지 저자 확인.
3. **2.6 선택압 데이터셋** — v2 "18 leaves" → 개정 "dataset sizes varied; S26 excluded from *ureB*". yn00 요약(hero-hero 쌍 14/8/12/12)과 부합하므로 개정본이 더 정확.
4. **투명성 항목 삭제** — v2에 있던 (i) 20-genome 참조 패널이 QC에서 전면 재구축된 사실(2.3·4.5), (ii) DRAM 0행·coverage 라벨 교정 공개(3.5·3.11·4.5), (iii) 12개 후보 규칙 검토(DEFINITION_AUDIT) 언급, (iv) geNomad MGE 부하 수치(*P* = 0.063/0.114), (v) GS–GOGAT 비유의 고갈(0.62), (vi) Nha 2.50 vs 2.19 수치. 현재 출하 표는 모두 검증판이므로 과학적 오류는 아니지만, 재구축 이력 공개 여부는 저자 결정.
5. **그림·표 라벨 "MICP-complete" vs 본문 "prioritized candidates"** — Methods 2.2 괄호로 연결해 두었음. 대안: 그림 10종·Table 2 재라벨(빌더 상수 변경 후 재빌드).
6. **Table 2의 "Tetracycline / macrolide determinants" 4.22배 enrichment** — 본문 3.5는 언급하지 않고 3.7은 "no acquired AMR genes". 키워드 스캔(product annotation) vs ResFinder(curated DB)의 차이는 Fig S1 범례에만 있음. 본문에 한 문장 권장(v2도 같은 공백).
7. **3.10 *ureC* 계통 해석** — "lineage-associated retention, gene loss, and horizontal acquisition rather than inheritance from a single recent ancestor"는 3.3의 "little evidence of recent mobilization"과 양립(오래된 HGT). 문제 없음, 기록만.

### 5c. 범례 ↔ 그림 대조
새 범례는 저자가 재작성한 것이라 그림 3종을 시각 대조: Fig 1(링 초록/흰 ✓, 후보 coral ✓), Fig 5(95 % 점선 ✓, PERMANOVA 값 ✓), Fig S4(em dash = no alignment ✓). 나머지 7종은 패널 문자·콜아웃 감사만 통과 → **제출 전 저자 시각 확인** 필요.

## 6. 저자 입력이 필요한 것 (본문에 `AUTHOR VERIFICATION REQUIRED` 표시)
1. 교신저자 이메일.
2. CRediT 역할 배정(5인).
3. GenAI 선언 문구·범위 확인(Claude: 정합성 검토·참고문헌 검증·언어 교정).
4. Acknowledgements — 개정본에서 삭제됨(v2는 자리표시). 필요하면 추가.
5. **Data availability** — 로컬 경로·`PRJNA-XXXXXXX` 잔존, 09-10에 발견한 출처 모순(공개 메타게놈 재분석 vs "우리가 기탁") 미해결. accession 확정 후 §2.1·DA 정정.
6. Highlights 5개가 143–171자로 MicRes 한도(85자)를 초과. 제안:
   * 111 livestock-waste MAGs yield six MICP candidates in two distinct lineages
   * Urease loci show little evidence of recent mobilization
   * UreC catalytic residues and predicted fold are conserved in all six candidates
   * Candidates are enriched for alkaline-stress, oxidative and CAZyme traits
   * S13 and S16 are candidate novel *Sphingobacterium* species for validation
7. Keywords 8개 → MicRes 6개 이내 권장(예: MICP; metagenome-assembled genomes; *Sphingobacterium*; urease; carbonic anhydrase; livestock waste).
8. Intro Lapierre 문장(§2 표 ⚠).
9. 5b-1·2·4·5·6 결정.
10. MicRes 커버레터 신작(현재 MB·mSystems용만 존재), PNG 500 dpi 재출력(PDF 제출 시 불필요).

## 7. 산출물
`01_Manuscript.md/.docx` · `02_Figure_legends.md/.docx` · `references/`(라이브러리·EndNote docx·`04_apply_to_master.py`) · `Supplementary_tables/`(S1–S3 새 체계) · `Main_tables/`(재빌드) · `port_revised_260923.py` · `../consolidation_260904/audit_{consistency,numbers}.py`(개정) · 보관: `*.pre_revised_260923.*`, `_superseded_260923/`.
