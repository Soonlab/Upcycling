# Microbial Genomics 규격 대조 — 저자 수정본 `MicrobialGenomics/01_Manuscript_MicrobialGenomics_workhere.docx` (2026-09-27)

대상: 저자 09-25 수정본(09-26 Vancouver 변환, 서론부터 번호 시작). 공식 사이트는 여전히 403 →
출처는 Wayback 아카이브: Prepare an article 2026-02-11 · Article types 2026-04-23(09-23 대조 때 없던 스냅샷, 단어 한도 확인됨) ·
Open Data 2026-06-18 · Ethics policies(AI 정책) 2026-04-09 · Submission and peer review 2026-03-31.
원칙: 초회 투고 format-free. 학회 서식 요청은 revision 단계.

## 통과(확인 완료)
- 단일 Word 파일, 그림 5·표 2·범례 원고 안 포함, 연속 줄번호, A4, Times New Roman, 쪽번호.
- 제목면: 제목·저자 4인·소속 3곳·교신 이메일·초록(246단어, 인용 없음, 첫 문장 주제·끝 문장 결론)·키워드 6(규정 3–6).
- Research Article 단어 한도 3,000–7,000(지시적) ↔ 본문(서론→결론) 4,941. 통과.
- Data Summary 있음: 원데이터 BioProject PRJNA1231077 + Zenodo 10.5281/zenodo.15309541, MAG별 SRA run·raw-read BioSample·MAG BioSample(S2P, 111행) → "BioProject만 말고 개별 accession" 규정 충족. 코드 GitHub Soonlab/Upcycling 공개(MIT) 확인.
- Impact Statement 143단어, 인용 없음.
- Conflicts 표준 문장 · Funding 과제번호+펀더 역할 · Ethical approval · 데이터 생산자 Acknowledgement · GenAI 공개(학회 AI 정책: 원고 또는 감사의 글에 사용 요소 명시 → 충족).
- 참고문헌 Vancouver 58편, DOI 전부, 번호 결번 없음, 첫 인용 순(1 = 서론 Fujita 2004).
- 보충자료: 단일 PDF(11쪽, Fig S1–S5+범례) + 단일 xlsx(49시트). 본문이 부른 시트 22개(S1A…S3J) 전부 존재. Fig S1–S5 콜아웃 전부 존재.
- 그림 콜아웃 순서 Fig 1→5, S1→S5 오름차순.

## Critical — 제출 자체를 막는 빈칸 (편집실이 DOI·빈칸 검사)
| # | 위치(문단) | 현재 | 규정 | 조치 |
|---|---|---|---|---|
| C1 | Data Summary ¶14 | `https://doi.org/AUTHOR VERIFICATION REQUIRED` | 파생 데이터·코드는 제출 시 DOI 존재 필수(심사 중 비공개 가능, 엠바고 불가; 편집실이 DOI 유효성 검사) | Zenodo 레코드 예약 후 DOI 기입. 커버레터 동일 |
| C2 | Author contributions ¶130 | 플레이스홀더. "five authors"라 쓰였으나 저자란은 4인, 역할 세트는 3개 | CRediT 권장(빈칸·지시문은 불가) | 4인 CRediT 역할 확정 기입. "five" 문구 제거 |
| C3 | Acknowledgements ¶132 · ref 48 ¶184 | `AUTHOR VERIFICATION REQUIRED` 2곳 | 웹사이트 참고문헌은 접속 날짜 필수 | ABRicate 접속일 기입, 추가 감사 문구 확정 또는 지시문 삭제 |

## Major
| # | 위치 | 현재 | 규정 | 조치 |
|---|---|---|---|---|
| M1 | Methods 2.x 소프트웨어 | 버전 없음: MMseqs2, ETE3, TM-align, FastTree, Biopython, gRodon2, VFDB(DB 날짜), PlasmidFinder(DB 버전) | "Please consistently cite any software used, including its version and parameters" | 저장소 확인값 Biopython 1.87. 나머지는 저자 확인 후 기입 |
| M2 | Data Summary / 참고문헌 | 원데이터는 논문 [18](Data in Brief)만 인용 | Open Data: 재사용 데이터셋 자체를 참고문헌으로 권장 — `Author. Title. Repository. Accession/DOI. (YYYY)` | 데이터셋 항목 추가 권장. 서론-첫 번호 규칙 유지하려면 Methods 2.1의 [18] 옆에 인용 → 새 19번, 이후 19–58 → 20–59 재번호(변환 스크립트로) |
| M3 | 그림 패널 문자 | Fig 1–5·S1–S5·범례·콜아웃 22곳 대문자 (A) | 소문자 괄호 (a) | revision 단계 요청 예상. 빌더에서 일괄 변경(사용자 결정) |
| M4 | 제목 번호 | "1. Introduction", "2.1 …" | 학회 템플릿은 무번호 Methods/Results/Discussion | revision 단계 |

## Minor
| # | 위치 | 현재 | 조치 |
|---|---|---|---|
| m1 | 초록 | TIM-barrel, TM-score, M0, ω 미정의 | "any abbreviations used must be defined" → (β/α)₈ triosephosphate-isomerase (TIM)-barrel, template-modeling (TM) score 등 |
| m2 | Results 3.8 ¶77 | "S. detergens", "S. multivorum" 첫 사용부터 약칭 | 첫 사용은 속명 전체 표기 |
| m3 | Funding ¶125 | 과제 수혜 저자 미명시 | "Any authors who are associated with specific funding sources should be named" → (to S.-C.K.) 등 저자 확인 |
| m4 | Table 2 | "Oxidative-stress defence" ↔ 본문 "oxidative defense" | 미국식으로 통일(Fig 4·워크북 동시, 기존 미결) |
| m5 | Results 3.1 | Pseudomonadota·Bacteroidota 비이탤릭, 속명은 이탤릭 | 학회 IJSEM 관행은 전 계급 이탤릭. 일관성 선택 |
| m6 | 보충 그림 | 범례에만 있고 본문 미인용 패널: S2C,D · S4A · S5C,D | 규정 아님. 인용 추가 또는 패널 정리 선택 |
| m7 | 보충표 인용 순 | S2 → S1 → S3 순으로 첫 인용 | 규정 없음. 선택 |
| m8 | 그림 alt-text | 없음 | 권장(≤255자). 없으면 제작사가 생성 |
| m9 | ORCID | 없음 | 투고 시스템에서 입력 권장 |
| m10 | GenAI 절 제목 | Elsevier 문구 | 학회는 위치 지정 없음(원고 또는 감사의 글). 유지 가능 |
| m11 | 교신 이메일 | gmail | 학회 규정 무관. 기관협약 무료 OA 경로만 기관 메일 필요 |

## 범위 밖 주의
- 저자 09-25 수정 수치는 44/38/29 외 미감사(09-26 README). 규격 대조와 별개로 수치 감사 필요.
- workhere가 md master보다 앞서 있음 → 빌더 재실행 금지, 이식 후 재빌드.

## 2026-09-27 Major·Minor 반영 결과 (workhere.docx 직접 편집, 변경 = 파란색 글자)
스크립트 `/tmp/.../edit_workhere.py`(세션 임시) → 수정 전 원본 `MicrobialGenomics/_backup_260927/01_Manuscript_MicrobialGenomics_workhere.pre_major_minor_260927.docx`.
검증: 본문 인용 묶음 41개 전부 원본과 동일 논문 지시(19번 삽입분 제외), 첫 인용 순 1–59 결번 없음, 참고문헌 본문 59편 텍스트 보존, 의도 밖 문단 변경 0, 본문 4,999단어.

| 항목 | 처리 |
|---|---|
| M1 버전 | MMseqs2 release 18(dram_env) · HMMER v3.4(전 env 동일) · ETE3 v3.1.3(dram_env) · TM-align 20240303(C4_esmfold install.log) · FastTree v2.2.0(설치 env 3곳 전부 동일, PATH 호출) · Biopython v1.87 · gRodon2 v2.7.2(grodon env) · **ABRicate v1.0.1 → v1.4.0 정정**(biosafety env `abricate --version`, 실행 2026-04-18) + 번들 DB 날짜 2026-04-03. CARD v3.2·dbCAN v12·ResFinder v4는 근거 미확인이라 그대로 → 🔴저자 확인 |
| M2 데이터셋 참고문헌 | Methods 2.1 `[18, 19]` + 새 19번(Zenodo 10.5281/zenodo.15309541, 레코드 제목·저자 API 확인) → 구 19–58은 20–59로 재번호 |
| M3 패널 소문자 | `_style.py` `UPCYCLING_PANEL=lower` 스위치(기본값 upper 유지) → 10장 재빌드(감사 전부 PASS) → 패키지 `Figures/`·docx 삽입 PNG·보충 PDF(범례 44개 소문자) 교체. 본문 콜아웃 25곳 `Fig. 1a` 형식(MGen 게재본 PMC13580827 관행), 범례 16곳 `(a)`. 빌더 출력 폴더는 대문자로 복원. **v4 master(`SUBMISSION_v4/Figures`, md)는 미변경** |
| M4 제목 번호 | 27개 제목에서 번호 제거 |
| m1 초록 약어 | TIM·TM·M0/ω 정의 |
| m2 종명 | Results 3.8 *Sphingobacterium detergens*·*S. multivorum* 전체 표기 |
| m3 펀딩 수혜자 | 근거 없음 → 미처리(저자) |
| m4 defence | Table 2(docx)·`Main_tables/Table_2…xlsx`·Fig 4(`build_v2_fig4.py`)·`build_main_table2.py` 전부 defense |
| m5 문 이탤릭 | Pseudomonadota(3)·Bacteroidota(1) 이탤릭 |
| m6 미인용 패널 | Fig. S2c,d(3.5)·S4a(3.9) 콜아웃 추가; **S5c,d는 3.11에 결과 문장 신설**(GC3 39.2% vs 80.5%, P=0.043; ENC 55.3 vs 37.0, P=0.0028 — S2O 시트 재계산 일치, 그림 P값과 일치) → 저자 검토 |
| m7 보충표 순서 | 규정 없음·워크북 재번호 부담 커서 미처리 |
| m8 alt-text | 5장 docPr descr(236–255자) |
| m9 ORCID·m11 이메일 | 투고 시스템/저자 |
| m10 GenAI 제목 | 유지 |

## 2026-09-28 그림 소문자 패널을 v4 master·저장소에 전파
- `new_figure/_style.py` 기본값 = 소문자 `(a)`(`UPCYCLING_PANEL=upper`로 복귀 가능) → 10장 재빌드(감사 전부 PASS) → `SUBMISSION_v4/Figures/`(pdf·png·svg) 교체.
- `02_Figure_legends.md` 38개 `**(a)**`(+ `(a–e)`·`panels c and d`), 수정 전 `02_Figure_legends.pre_lowercase_260927.md`; `02_Figure_legends.docx` 재생성. `01_Manuscript.md` 콜아웃은 원래 `Fig. 1a`라 무변경.
- `build_micgen_package.py` 대문자 변환 제거(수정 전 `_build/build_micgen_package.pre_lowercase_260927.py`), 패키지 재빌드. 감사 2종×2사본 legend 정규식 `[A-Ea-e]`.
- 이제 workhere.docx와 master의 그림·범례·콜아웃 표기가 일치. GitHub `figures/`·`figures/build/`·`scripts/submission_v4/` 동기화.

## 2026-09-28 보충표 S1↔S2 맞바꿈 (m7 처리, 사용자 지시)
- 첫 인용 순이 S2(2.1 S2P)→S1→S3였음 → 워크북 번호 교환: **S1 = per-MAG(구 S2, 시트 S1A–S1T), S2 = reference panels(구 S1, S2A–S2O)**, S3 불변. 시트 안 문자는 그대로(예: S2P→S1P).
- 적용: `Supplementary_tables/` 파일명·시트명·README 셀(원본 `_pre_swap_260928/`) · `01_Manuscript.md`(17토큰)·`02_Figure_legends.md`(블록 순서도 S1→S2→S3; `*.pre_swap_260928.md`) · docx 재생성 · 빌더 Data Summary 문구 · 감사 주석 · **workhere.docx 콜아웃 18곳(파란색, `_backup_260927/…pre_swap_260928.docx`)** · 패키지 재빌드(xlsx 시트명·Contents·보충 PDF 목차).
- 스크립트 `/tmp/.../swap_s1s2.py`(세션 임시). 감사 결과는 아래 줄 참조.
- 감사: MGen 사본 구조 ALL CHECKS PASS·수치 ALL NUMBERS AGREE, master 동일. workhere 첫 인용 순 S1→S2→S3, 인용 시트 21종 전부 존재.
