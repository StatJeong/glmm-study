# INLA (Integrated Nested Laplace Approximation)

지금까지 본 GVA/AGHQ와는 **사상 자체가 다른** 접근이다.
- GVA: ELBO라는 목적함수를 만들어 최적화한다.
- AGHQ: 각 그룹의 사후분포를 수치구적으로 직접 적분한다.
- **INLA**: 사후분포를 최적화 문제로 보지 않고, **중첩된 Laplace 근사**로 직접 근사한다.
  잠재 가우스 모형(latent Gaussian model) 전체 — GLMM은 물론 공간·시계열 랜덤효과까지 —
  를 하나의 틀로 묶는다.

## 대상 모형: 잠재 가우스 모형(LGM)

```
y_i | x_i, θ  ~  지수족(η_i = x_i)          관측모형
x | θ         ~  N(0, Q(θ)^{-1})            잠재 가우스장 (희소 정밀행렬 Q)
θ             ~  π(θ)                       초모수 (분산성분 등, 저차원)
```

GLMM의 랜덤효과 u_i, 공간 랜덤장, 시계열 성분, 스플라인 계수가 전부 이 x 안에
들어간다 — 즉 **단일/중첩/교차/공간을 전부 같은 틀에서 다룬다.**

핵심 장치는 **GMRF(Gaussian Markov Random Field)**: Q가 희소행렬이면 Cholesky 인수분해가
O(n) ~ O(n^1.5)로 되어, 우리가 계속 부딪힌 "K가 커지면 느려진다"는 문제를 아예 다른
방식으로 피해간다.

## 읽는 순서

| 순번 | 파일 | 쪽수 | 내용 |
|---|---|---|---|
| 1 | `RueMartino_2011_FastApprox_PastPresentFuture.pdf` | 8 | **여기서 시작.** 짧은 개관, 왜 만들었는지 |
| 2 | `RueRieblerSorbyeIllianSimpsonLindgren_2017_Review_AnnRevStat.pdf` | 28 | **본격 입문.** 최신 리뷰, 수식과 예제가 잘 정리됨 |
| 3 | `RueMartinoChopin_2009_INLA_original_JRSSB.pdf` | 74 | **원논문(게재본, 토론 포함).** 세 단계 중첩 Laplace 근사의 상세 유도 + 여러 저자의 토론·반박 |
| 4 | `MartinsSimpsonLindgrenRue_2013_NewFeatures_CSDA.pdf` | 16 | 소프트웨어 확장 정본(게재본, 비가우스 관측모형·사후예측 등) |
| 5 | `MartinsSimpsonLindgrenRue_2013_NewFeatures_arXivExtended.pdf` | 45 | 위 논문의 arXiv 확장판. 게재본에 없는 부록·예제 다수. 필요할 때 참조용 |
| 6 | `LindgrenRueLindstrom_2011_SPDE_JRSSB.pdf` | 76 | **SPDE 접근.** 공간 랜덤효과 ↔ GMRF를 잇는 다리. 공간모형 다룰 때 필수 |

3, 6번은 게재본(정식 페이지 매김, Royal Statistical Society 저작권 표시 있음) — 기관 접속으로
정식 구하신 것으로 보입니다.

## 세 단계 근사 (2009 논문의 핵심 구조)

INLA는 이름 그대로 "nested"다 — Laplace 근사를 세 겹으로 중첩해서 쓴다.

1. **π̃(θ|y)**: 초모수(θ, 저차원)의 사후를 Laplace 근사
2. **π̃(x_i|θ,y)**: 각 잠재변수의 조건부 사후를 Laplace 근사 (GMRF 구조로 빠르게)
3. 위 두 근사를 결합해 **π̃(x_i|y) = Σ_θ π̃(x_i|θ_k,y) π̃(θ_k|y) Δ_k**로 수치적분(중첩 격자)

이 3단계 구조가 "Integrated **Nested** Laplace Approximation"의 "Nested"다.
우리가 GLMM에서 쓴 "nested random effects"의 nested와는 다른 의미이니 헷갈리지 말 것.

## 우리 작업과의 연결점 — 읽으면서 확인할 것

- **K가 커질 때**: 우리는 GVA가 O(q) vs AGHQ의 O(q^K)로 이겼다. INLA는 아예 다른 무기
  (희소성)를 쓴다. 세 방법을 같은 축에 놓고 비교할 수 있는가?
- **교차설계**: `03_crossed`에서 본 GVA+복합가능도, Goplerud의 부분인수분해와 달리,
  INLA/GMRF는 교차 구조를 원래 잘 다루는가, 아니면 여기서도 문제가 생기는가?
- **정확도**: 논문들이 INLA를 MCMC와 얼마나 가깝다고 주장하는지, 그 근거가
  Hall-Ormerod-Wand 식의 점근이론인지 경험적 비교인지 구분해서 볼 것.
- **소프트웨어**: R-INLA 패키지로 실제로 우리 epilepsy 데이터를 적합해서
  GVA/AGHQ/glmmTMB(Laplace)와 다섯 갈래 비교가 가능하다.

## 관련 코드

R-INLA는 CRAN에 없고 별도 저장소에서 설치한다:
```r
install.packages("INLA", repos = c(getOption("repos"),
  INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```
