# Microbiological Research 규격 대조 — SUBMISSION_v2 (2026-09-11)
가이드라인 원장: `/data/data/Target_Journal/MicrobiologicalResearch_Guide_for_Authors_2026-09.md`
(공식 페이지 봇차단 → 스니펫+2016 아카이브+APC 가격표로 재구성. 제출 전 브라우저 대조 1회 필요)

| # | 항목 | 규정 | 현재 (01_Manuscript.md 2026-09-09) | 판정 | 조치 |
|---|---|---|---|---|---|
| 1 | 본문 길이 | ≤ 8,000 (Intro→Conclusion) | 6,999 (README) / 7,213 (헤딩 포함 재측정) | ✅ | 없음 |
| 2 | 섹션 구조 | 1 Intro · 2 Materials and methods · 3 Results · 4 Discussion · 5 Conclusion | 동일 번호·순서 (2.1–2.10 / 3.1–3.11 / 4.1–4.5 / 5) | ✅ | 없음 |
| 3 | 초록 | ≤ 250 단어 | **09-23: 저자 개정본 198단어로 교체** | ✅ | 249단어판을 본문 초록으로 교체(Mrp 유병률 정정·S26 CA plasmid 등 09-04 이후 정정사항 반영됐는지 대조 후) |
| 4 | 하이라이트 | 3–5개 × ≤ 85자, 별도 파일 | 09-23: 5개이나 143–171자(제안문은 REVISED_PORT_REPORT §6) | ❌ | 3–5개로 압축·85자 이내 재작성, `06_Highlights.docx` 별도 파일 |
| 5 | 키워드 | ≤ 6(2016) / ≤ 7(스니펫), 영국식 철자 | 09-23: 8개 | ❌ | 6개로 축소(예: MICP · metagenome-assembled genome · *Sphingobacterium* · urease · carbonic anhydrase · alkaline tolerance) |
| 6 | 참고문헌 스타일 | 일관된 스타일(저널 스타일=Elsevier Harvard 저자–연도) | 09-19 Harvard 60편 → 09-23 개정본 이식 후 57편 | ✅ | 목록 63건을 Harvard로 재포맷(번호 제거, `Surname, I.I., Year. Title. J. Abbr. vol, pages.`), 알파벳→연도순. 메모리에 "20건은 번호로만 도달"이라 기록됐으나 본문 grep에 `[n]` 인용 0건 → 실제 미인용 여부 재검 필요 |
| 7 | 그림 형식 | 벡터 EPS/PDF(폰트 임베드) 또는 혼합 ≥ 500 dpi(전폭 3,740 px) | PDF 벡터·LiberationSans 임베드 ✅ / PNG 200 dpi 1,417 px ❌ | ⚠ | **PDF를 제출본으로**, PNG는 ≥500 dpi로 재출력(빌더 dpi 인수) 또는 제외 |
| 8 | 그림 크기 | 인쇄 규격에 맞춤 | 180 mm × ≤ 235 mm, 1장/페이지 | ✅ | 없음 |
| 9 | 캡션 | 그림과 분리, 제목+설명 | `02_Figure_legends` 별도 ✅ | ✅ | 없음 |
| 10 | 그래픽 초록 | ≥ 531×1328 px, 5×13 cm 가독, TIFF/EPS/PDF | 2,244×1,240 px(300 dpi, 190×105 mm) PDF/PNG/SVG | ✅ | 없음(비율 13:5 아니어도 무방) |
| 11 | 표 | 편집 가능 텍스트, 세로괘선 없음 | Table 1·2 markdown→docx | ✅ | docx 변환 후 괘선 확인 |
| 12 | 보충자료 | 그대로 게시, 파일별 캡션 | 워크북 3 + Fig S1–S5 + DRAM xlsx + iqtree | ✅ | 파일별 한 줄 캡션 목록 추가 |
| 13 | CRediT | 필수 | 있음(플레이스홀더 저자) | ⚠ | 실명 |
| 14 | 경쟁이익 | 필수 | 있음 | ✅ | — |
| 15 | 펀딩 | Elsevier 표준 문장 | `[Funder, Grant ID]` 플레이스홀더 | ⚠ | 과제번호 |
| 16 | GenAI 선언 | 필수(작성 과정 AI 사용 시) | 09-23 절 추가(문구 저자 확인) | ✅ | "Declaration of generative AI and AI-assisted technologies in the writing process" 절 추가 |
| 17 | Data availability | 필수 | 로컬 경로 + `PRJNA-XXXXXXX` + 09-10 발견한 출처 모순(공개 메타게놈 재분석인데 "우리가 기탁") | ❌ | 원 accession 확정 후 §2.1·DA·§4.5 정정(09-10 항목) |
| 18 | 제목 페이지 | 저자·소속·교신 이메일 | 플레이스홀더 | ⚠ | 사용자 입력 |
| 19 | 제목 | 약어 회피 | "MAG"는 풀어씀, "MICP"는 풀어씀 ✅, 길이 매우 김(≈40단어) | ⚠ | 단축 권장 |
| 20 | 본문 헤더 메모 | — | 09-23: 배너를 HTML 주석으로 → docx에 없음 | ✅ | 제출본 docx에서 제거 |
| 21 | 커버레터 | — | MB용 04·mSystems용 05만 존재 | ❌ | `06_Cover_letter_MicrobiologicalResearch` 신작(스코프 근거: Biotechnology·Environmental Microbiology 섹션) |
| 22 | 비용 | 구독형 0 / OA USD 4,790 | 게재료 상한 ₩500만 | ✅ | 구독형 선택 |
