# glmm-study

Generalized Linear Mixed Model(GLMM) 계산방법을 공부하면서 쌓은 자료 모음입니다.
"랜덤효과가 있는 우도의 적분을 어떻게 근사할 것인가"라는 하나의 문제를 중심으로,
읽은 논문의 정리, 직접 짠 시뮬레이션 코드, 관련 논문 리뷰 발표자료를 한 곳에 모았습니다.

## 구성

### [`notes/`](notes/) — 계산방법 로드맵

랜덤효과 구조(단일 grouping / 중첩 / 교차 / 고차원 / INLA)별로 어떤 근사 방법
(Laplace, adaptive Gauss-Hermite quadrature, Gaussian variational approximation, INLA)이
가능하고 왜 그런지 정리한 노트입니다. 저작권이 있는 원논문 PDF는 올리지 않고,
직접 정리한 로드맵과 구현 코드만 포함했습니다.

- [`README.md`](notes/README.md) — 구조별 계산방법 로드맵, 읽는 순서, 이론과 시뮬레이션 결과의 연결점
- [`inla-notes.md`](notes/inla-notes.md) — INLA(Integrated Nested Laplace Approximation) 정리
- [`gva_single_grouping.R`](notes/gva_single_grouping.R) — Poisson/Gamma random-intercept GLMM에 대한 GVA 추정량 구현

### [`blog/`](blog/) — 시뮬레이션 코드와 학습 기록 (Quarto 블로그)

GVA와 AGHQ를 비교하는 시뮬레이션 연구, Polya-Gamma 샘플링, 준모수 혼합모형,
교재(Jiming Jiang, McCulloch) 연습문제 풀이 등을 담은 Quarto 블로그 프로젝트입니다.
`quarto render`로 정적 사이트를 만들 수 있습니다.

주요 내용:
- `glmm_multi.R`, `gva_vs_aghq.R`, `run_K_simulation.R`, `run_dimension_study.R`,
  `run_n_scaling.R` — GVA/AGHQ 비교 시뮬레이션 (Poisson/Gamma/binomial, SQUAREM 가속)
- `epilepsy_fit.R`, `epilepsy_benchmark.R` — Ormerod & Wand의 epilepsy 예제 재현 및 벤치마크
- `discuss_prompt.md` — 시뮬레이션 중 발견한 버그를 계기로 원논문(Hall, Ormerod & Wand, 2011)의
  실제 주장을 다시 검증한 기록. 코드 검증과 논문 이해를 어떻게 교차 확인했는지 보여주는 사례
- [`posts/study/`](blog/posts/study/) — Polya-Gamma, GHQ, VGA 등 학습 노트 및 교재 연습문제 풀이

### [`paper-reviews/`](paper-reviews/) — 논문 리뷰 발표자료

- [`lmmnn/`](paper-reviews/lmmnn/) — Simchoni et al.의 LMMNN(딥러닝에 랜덤효과 결합) 논문 리뷰

## 배경

통계학과 대학원 세미나에서 GLMM의 계산방법(근사 추론)을 주제로 진행한 개인 스터디 기록입니다.
이론(원논문이 증명한 것)과 구현(직접 짠 코드의 시뮬레이션 결과)을 계속 서로 대조하면서
검증하는 방식으로 작업했습니다 — `discuss_prompt.md`가 그 과정을 보여주는 예시입니다.
