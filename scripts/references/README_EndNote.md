# EndNote 참고문헌 패키지 (2026-09-18 작성 · 09-19 대체 확정 · 09-23 저자 개정본 이식)

## 2026-09-23 저자 개정본 이식 후 상태

정본 본문이 저자 개정본(`../Manuscript_revised.docx`, `../port_revised_260923.py`로 이식)으로 바뀌었다. 09-19 판정을 개정 문장에 다시 적용했고(`04_apply_to_master.py`의 09-23 항목: Cheng→Kumari, Krawczyk→Lapierre, Omoregie 2019, Dhami 2013, Li→Chen, Zamanzadeh 삭제·Gupta는 Discussion "reservoir" 문장으로, 도구 인용 Parks 2020·Zheng 2023·Yang and Nielsen 2000·Shaw and Yu 2023 재부착), 개정본이 더 이상 인용하지 않는 Hommel 2015·Lagesen 2007·Yuan 2015는 목록에서 빠져 **본문 인용 57편**이다. 라이브러리 `.ris`는 60편 그대로(3편 미인용, EndNote가 목록에 넣지 않음). `01_Manuscript_EndNote.docx` 임시 인용 37군데·61건·미검증 0. 상세는 `../REVISED_PORT_REPORT_260923.md` §2.

`../01_Manuscript.md`가 계속 정본이다. 이 폴더의 산출물은 전부 스크립트 생성물이므로 직접 편집하지 말고,
`ref_decisions.tsv`를 고친 뒤 재빌드한다.

## 파일

| 파일 | 내용 |
|---|---|
| `Upcycling_MICP_library.ris` | **검증된 60편**(전부 본문 인용). 전 저자·권·호·쪽·DOI를 CrossRef에서 DOI로 받아 생성(수기 입력 없음). Label 필드 = 원고 번호(`MS-R03`) |
| `Upcycling_MICP_candidates.ris` | (비어 있음 — 대체 확정으로 후보 없음) |
| `01_Manuscript_EndNote.docx` | 본문 인용 전부를 EndNote 임시 인용 `{Bowers, 2017}`으로 바꾼 원고(52군데·67건, 미검증 0) |
| `ref_decisions.tsv` | 63편 각각의 판정(keep / manual / pending)·DOI·대체 후보·사유 — **원장** |
| `library_audit_original_list_260918.tsv` | (기록) 수정 전 목록의 연도·권·첫 쪽 vs CrossRef 대조 |
| `04_apply_to_master.py` | 판정을 정본 md에 반영(본문 인용 교체 + Harvard 목록 재생성), 멱등 |
| `endnote_docx_report.txt` | 변환 건수, 미인용 문헌, 연도 불일치 |
| `crossref_match_report.tsv`, `crossref_candidates.json`, `crossref_by_doi.json` | 조회 원자료·캐시 |

## PC에서 할 일 (한 번)

1. EndNote → `File ▸ New`로 새 라이브러리 생성 → `File ▸ Import ▸ File…` →
   `Upcycling_MICP_library.ris`, Import Option = **Reference Manager (RIS)**, Text Translation = **Unicode (UTF-8)**.
2. `01_Manuscript_EndNote.docx`를 Word로 열고 EndNote 탭 ▸ Style = **Microbiological Research**
   (없으면 endnote.com/downloads/styles 에서 받거나 `Elsevier Harvard (with titles)`).
3. EndNote 탭 ▸ **Update Citations and Bibliography**. `{저자, 연도}` 52군데(개별 인용 67건)가 정식 필드로 바뀌고
   참고문헌 목록이 문서 맨 끝의 `References` 제목 아래에 생긴다(그래서 이 docx에서는 표를 References 앞으로 옮겨 두었다).
   동명·동년 문헌이 없어 선택 창은 뜨지 않아야 한다.

종 이름 이탤릭(*Sporosarcina pasteurii* 등)은 RIS로 전달되지 않는다. 해당 레코드(R58 Benini 등) 제목에서 EndNote로 직접 이탤릭 지정.

## 2026-09-19 확정: 실재하지 않던 9편 처리 (정본 `01_Manuscript.md`에 반영 완료)

원 기재 9편은 CrossRef·PubMed 어디에도 없었다(이전 판정 기록은 `ref_decisions.pre_replace_260919.tsv`, 수정 전 원고는 `../01_Manuscript.pre_refs_260919.md`).

| # | 원 기재(실재 안 함) | 확정 | 근거 |
|---|---|---|---|
| 6 | Cheng 2017 *ESPR* | **Kumari et al., 2016** *Adv Appl Microbiol* 94:79–108 | 인용 문구 "immobilisation of toxic metals"를 직접 다루는 MICP 리뷰 |
| 13 | Jiménez-Martínez 2022 *WRR* | **Hommel et al., 2015** *WRR* 51:3695–3715 | 다공성 매질 MICP 수송·침전 모델 |
| 16 | Krawczyk 2021 *Front Bioeng* | **Lapierre et al., 2020** *Sci Rep* 10:22448 | *S. pasteurii* 영양요구성·복합배지 의존. ⚠ 원 문장의 "pH ≈ 7–8, freshwater"는 출처가 없어 **문장을 이 논문이 실제로 보인 내용으로 교체**(단어 수 동일): "whose auxotrophies and reliance on complex laboratory media are poorly matched to real waste streams" |
| 18 | Lee 2022 *PeerJ* | **Yuan et al., 2015** *Bioinformatics* 31:i35–i43 | 메타게놈에서 16S rRNA 유전자 조립이 어려운 이유 |
| 19 | Li 2022 *MRA* | **Chen et al., 2021** *Nat Commun* 12:1106 | 가축(돼지 장) MAG 6,339개 카탈로그, 본문 연도 2021 유지 |
| 22 | Omoregie 2019 *CBM* 225:1108 | **Omoregie et al., 2019** *CBM* 228:116828 | 실재하는 같은 저자·저널·연도 논문(모래 biocementation). 본문 연도 2022→2019 |
| 36 | Xu 2021 *AMB* | 삭제 | 본문 미인용 |
| 37 | Zamanzadeh 2023 | **Gupta et al., 2016** *Bioresour Bioprocess* 3:28 | 우분을 미생물·효소·생물정화 자원으로 정리한 리뷰 |
| 63 | Suzuki 2022 "TranslatorX update" | **Biopython (Cock et al., 2009)** 으로 대체 | 저장소 `scripts/additional/C3_dnds_codon/run_dnds_v2.py::msa_codon` 확인: PAL2NAL·TranslatorX를 돌린 적 없고 자체 Biopython 루틴으로 역번역함 |

함께 정정: Dhami 2014→2013(본문 3곳), 참고문헌 목록을 **저자-연도(Elsevier Harvard) 알파벳순 58편**으로 재생성(저자 6명 초과 시 6명 + et al., DOI 포함),
본문 6,999→6,997단어(§2.7의 중복 Biopython 인용 1개 제거로 상쇄). `audit_consistency.py`·`audit_numbers.py` 정본 기준 전항 PASS
(`UPCYCLING_MAN_DIR=/data/data/Upcycling/SUBMISSION_v2` 지정 — 기본 경로는 09-04 사본을 읽는다).

### 2026-09-19 (2) 저자 결정 반영

* **Stegen et al., 2013 삭제** — 인용 문장("축분은 유용 기능의 저장소")과 무관한 논문이라 본문·목록·라이브러리에서 제거. 해당 문장은 Gupta et al., 2016 단독 인용.
* **Methods 도구 인용 3건 추가** — GTDB r220 `(r220; Parks et al., 2020; Chaumeil et al., 2022)`, `dbCAN v12 (Zheng et al., 2023)`, `PAML yn00 (Yang and Nielsen, 2000)`.
* 결과: 참고문헌 **60편**(전부 본문 인용·DOI 검증), 본문 **7,005단어**(자체 목표 7,000은 저자 결정으로 인용에 한해 초과 허용, 감사 상한 7,050), 임시 인용 52군데·67건.

## 실재하지만 원고 기재가 틀렸던 것 (정정 완료)

* **R8 Dhami**: 2014 → **2013** (*J Microbiol Biotechnol* 23:707–714), 저자 3인은 PubMed 기준(CrossRef는 제1저자만 보유).
* **R56 TM-align**: CrossRef가 Zhang만 보유 → PubMed 기준 Zhang & Skolnick.
* **R38 dbCAN3**: 제1저자 Zhang H → **Zheng J**.
* **R32 Stegen**: 제목(2012 논문)과 서지(2013 논문)가 섞여 있던 항목 — 이후 인용 자체를 삭제.
* R3 Bowers(컨소시엄이 제1저자로 등록 → 사람 저자 우선), R7 DeJong(대문자 저자명), R59 Mitchell·R27 Schwengers(권·쪽 PubMed 보충), R9·R55(정식 제목).

## 재빌드

```bash
cd /data/data/Upcycling/SUBMISSION_v2/references
PY=/home/soon/miniconda3/bin/python
$PY 01_crossref_lookup.py      # (기록용) 수정 전 번호식 목록 → CrossRef 후보 조회
$PY 02_build_library.py        # ref_decisions.tsv → .ris 2종 + citekeys.json
$PY 04_apply_to_master.py      # 정본 md 반영(본문 인용·Harvard 목록) → 이후 ../rebuild_docx.sh
$PY 03_build_endnote_docx.py   # 01_Manuscript.md → 01_Manuscript_EndNote.docx
```

문헌을 바꾸려면 `ref_decisions.tsv`의 해당 행(`action`, `replacement_doi`)과 `04_apply_to_master.py`의 `TEXT_EDITS`/`TITLES`를 함께 고친 뒤 02 → 04 → 03 순으로 실행.
