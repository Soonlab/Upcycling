# Microbial Genomics 규격 대조 — SUBMISSION_v4 (2026-09-23)
출처: Microbiology Society "How to prepare an article"(Wayback 2026-02-11), "Open Data"(Wayback 2026-06-18), MGen FAQ(Wayback 2025-09-11).
공식 사이트는 봇 차단(403). **Article types 페이지(단어 수 한도)는 아카이브 없음 → 미확인.**
핵심: 학회는 **format-free 초회 투고** 정책. 서식 요청은 revision 단계에서 함.

| # | 항목 | MGen 규정 | 현재 v4 | 판정 | 조치 |
|---|---|---|---|---|---|
| 1 | 제출 형태 | Word 단일 파일, 그림·표를 원고 안에 범례와 함께 | 원고·범례·그림·표 분리 | ❌ | 원고 docx에 Fig1–5·Table 1–2·범례 삽입 |
| 2 | 줄번호 | 연속 줄번호 필수 | 없음 | ❌ | docx 빌드에 줄번호 추가 |
| 3 | Impact Statement | MGen 대부분 유형에 필요 | 없음 | ❌ | 신규 작성(길이 한도 미확인) |
| 4 | Data Summary | 필수 절. 모든 데이터·코드의 DOI/accession·URL | "Data and code availability" 절, 파생 데이터는 "DOI on acceptance" | ❌ | Data Summary로 개편. **파생 데이터·코드 DOI는 초회 투고 시점에 있어야 함**(심사 중 비공개 가능, 채택 시 공개). 엠바고 불가 |
| 5 | 제3자 데이터 인용 | 생산자 인용 + DB 링크, 참고문헌에 `저자. 데이터 설명. 저장소. accession/DOI. (연도)` 형식 권장 | 논문 인용·PRJNA1231077·Zenodo 명시 | ⚠ | 데이터셋 자체를 참고문헌 항목으로 추가 |
| 6 | 개별 accession | 시료별 accession 제공, 식별자↔BioSample 연결 명확 | Table S2P에 MAG별 SRA/BioSample | ✅ | 없음 |
| 7 | 보충자료 | **단일 PDF 또는 단일 통합 Excel** | xlsx 4개 + .iqtree + Fig S1–S5 PDF 5개 | ❌ | Fig S1–S5 → 1 PDF, 표 → 1 워크북, .iqtree → Zenodo |
| 8 | 초록 | 참고문헌 금지, 약어 정의, 첫 문장 주제·마지막 문장 결론 | 209단어 | ✅ | 한도 미확인(대개 ≤250) |
| 9 | 키워드 | 3–6개 | 6 | ✅ | 없음 |
| 10 | Highlights | 규정 없음(Elsevier 항목) | 있음 | ⚠ | 삭제 |
| 11 | 그래픽 초록 | Microbiology 지만 권장, MGen 규정 없음 | 있음 | ⚠ | 제외(또는 선택) |
| 12 | 참고문헌 | 초회는 일관되면 자유. 채택 후 Vancouver로 재포맷, DOI 권장 | Harvard 58편, DOI 포함 | ✅ | 유지 |
| 13 | 패널 표기 | 소문자 괄호 (a), (b) | **본문 콜아웃 `Fig. 1a` ↔ 범례 `(A)` 불일치** | ❌ | 저널 무관 내부 오류. 범례·그림 문자를 (a)로 통일 |
| 14 | 그림 파일 | revision 때 별도 파일, ≥300 dpi, 단일 컬럼 A4 | 180 mm PDF 벡터 | ✅ | 없음 |
| 15 | 경쟁이익 | 없으면 "The author(s) declare that there are no conflicts of interest" | 다른 문장 | ⚠ | 표준 문장으로 교체, 제목 "Conflicts of interest" |
| 16 | 펀딩 | 기관·과제번호, 펀더 역할 진술 | 과제번호 있음 | ⚠ | 펀더 역할 문장 추가 |
| 17 | CRediT | 권장 | 플레이스홀더 | ⚠ | 저자 입력 |
| 18 | 윤리 | 인간·동물 연구 시 | 해당 없음 진술 있음 | ✅ | 없음 |
| 19 | GenAI 선언 | 이 페이지엔 규정 없음 | 있음 | ⚠ | 학회 AI 정책 별도 확인 후 유지 여부 결정 |
| 20 | 섹션 명칭 | Methods / Results / Discussion, 번호 없음 | "2. Materials and Methods" 등 번호 | ⚠ | revision 단계 사항, 번호 제거 권장 |
| 21 | 단어 수 | 미확인 | 본문 4,526 | ? | Article types 페이지 브라우저 확인 |
| 22 | 교신저자 | 기관 이메일(기관 협약 무료 OA 조건) | 플레이스홀더 | ⚠ | 수원대 이메일 입력, 기관 협약 여부 확인 |

## 2026-09-23 패키지 구성 후 상태 (`MicrobialGenomics/`, 빌더 `build_micgen_package.py`)
- 해소: 1 원고 단일 파일(그림·표·범례 삽입) · 2 줄번호 · 3 Impact Statement · 4 Data Summary(문안; Zenodo DOI는 저자) · 7 보충자료 PDF 1개 + Excel 1개 · 10 Highlights 삭제 · 11 GA 제외 · 13 콜아웃 대문자로 그림·범례와 통일 · 15 경쟁이익 표준 문장 · 16 펀더 역할 문장(저자 확인).
- 남음: 4 Zenodo DOI · 17 CRediT · 19 GenAI 정책 · 21 단어 수 한도 · 22 교신 이메일 · revision 단계 항목(소문자 패널, 헤딩 번호, Vancouver, "defence" 범주명).
- 감사: 구조 19/19 · 수치 129/129 PASS.
