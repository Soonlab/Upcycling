# SUBMISSION_v4 — 데이터 출처 확정에 따른 재프레이밍 보고서 (2026-09-23)

정본 = `SUBMISSION_v4/01_Manuscript.md`(v2에서 `reframe_260923.py`로 생성, 멱등). 근거 원장 = `data_provenance_260923/DATA_PROVENANCE_260923.md`.

## 1. 확정 사실
- 111 MAG = Pérez-Valera & Elhottová 2025 *Data Brief* 61:111748(PRJNA1231077 · Zenodo 10.5281/zenodo.15309541)과 111/111 서열 동일.
- 접두어 C/M/S/V = COMEbin / MaxBin+MetaBAT / SemiBin2 / AVAMB. 실제 출처 manure 44 · manured soil 38 · soil 29(체코 젖소 1농장 분변 + 유기농 토양 2종 마이크로코즘, 2·14·28일, CHROMagar Acinetobacter 농화배양).
- 후보 6: M1·S26 manure(Soil S 계열) / S13·S16·S23·C22 manured soil(Soil B 계열). S13·S16은 같은 day-28 공동조립(GT28BC2)에서 나왔으나 AAI 82.5 %로 별개 유전체.

## 2. 사용자 결정(09-23) → 실행
| 결정 | 실행 |
|---|---|
| 축종 PERMANOVA·coverage-by-source·축종 색상 삭제 | §3.10 source 문장 삭제(속 F=8.21만) · Fig 5D→5C 계통군 색 · Fig S5 B/C→B(원 저자 CoverM 매핑 coverage, 후보 중앙값 9.5× vs 17.6×, P=0.12) · S3D/S3L 삭제 · S2T Source→Origin |
| abundance proxy를 원 저자 매핑 coverage로 대체 | §2.9·§3.9·Fig S5B·Table S3J(신설)·S2P(provenance+coverage+SRA/BioSample) |
| MGnify 절 삭제, 속 내부 비교 유지 | §2.9·§3.9·§4.1·초록·결론·Fig 5B·S3G/S3H·Mitchell 2020 인용 제거; Pseudomonas_E 146 참조 비교(36.3 %)는 유지 |
| "젖소 분변–토양 마이크로코즘 농화배양 MAG 재분석"으로 재작성 | 제목·Highlights(5×≤85자)·초록(209단어)·키워드 6·서론·§2.1(출처·농화·조립·비닝 서술, "ACE hybrid" 삭제)·§3.1(Pseudomonadota 99/Bacteroidota 12·17속)·§3.5/§4.3 제목·§4.1·§4.4 한계(농화배양 편향 신설, 이후 항목 번호 이동)·§5·윤리·Data availability·Table 1 Origin 열·범례 7곳 |

## 3. 참고문헌
+ Pérez-Valera & Elhottová 2025(10.1016/j.dib.2025.111748) · + Sardar et al. 2023 FEMS Microbiol Ecol 99:fiad148(10.1093/femsec/fiad148) · − Mitchell 2020(MGnify, 미인용) → **58편**(전부 DOI 검증). EndNote 패키지 재생성(41군데·65건).

## 4. 검증
- `audit/audit_consistency.py` 20/20 PASS · `audit/audit_numbers.py` 129/129 PASS(신규 22항: 출처 수·접두어↔비닝기·후보 출처·coverage·문계·AAI·accession 인용·구 프레이밍 부재).
- 원고 본문·범례·docx에 swine/sheep/poultry/MGnify/ACE hybrid/PRJNA-XXXXXXX 0건(스크립트 assert + docx grep).
- 그림 3종 st.audit/prose_scan PASS, 시각 확인 완료(Fig 5·S5·GA).

## 5. 재검토 결과 참고(원장 판단용, 원고에는 미기재)
진짜 메타데이터로 다시 돌린 PERMANOVA: 접두어(=비닝기) trait z-score F=2.71 P=0.001은 v2 값 그대로 재현 → 구 "source 효과"는 비닝 아티팩트. 실제 출처 3군: pangenome F=2.25 P=0.001·trait F=1.89 P=0.005이나 속 조성과 교락(manure에 Acinetobacter·Stenotrophomonas 편중) → 원고에 넣지 않음.

## 6. 남은 저자 입력
교신 이메일 · CRediT · GenAI 문구 · Acknowledgements 여부 · 커버레터 빈칸(에디터·섹션·추천 리뷰어 3인) · 파생 데이터 Zenodo DOI · 변경 없는 그림 7장 ↔ 재작성 범례 시각 대조 · 공저자에게 데이터 입수 경위 확인(원 저자 인용은 완료, Zenodo 라이선스 조건 제출 전 확인).
