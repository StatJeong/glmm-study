# GLMM 계산방법 — 랜덤효과 구조별 자료 모음

모든 것이 하나의 문제에서 갈라진다:

```
ℓ(β,Σ) = Σ_i log ∫ exp{...} du_i          ← 이 적분이 닫히지 않는다
```

**폴더는 랜덤효과 구조로 나눴다.** 구조가 달라지면 적분이 분해되는 방식이 달라지고,
그에 따라 쓸 수 있는 계산법이 갈리기 때문이다.

| 구조 | 적분 분해 | AGHQ | 자료 |
|---|---|---|---|
| 단일 grouping | m개의 K차원 적분 | 가능 (q^K) | 7편 |
| 중첩 | 계층적 중첩 적분 | 원리상 가능 | 4편 |
| 교차 | **분해 안 됨** | **불가능** | 3편 |

---

## 00_foundation — 구조 무관 원전

방법론 자체의 뿌리. 어느 구조를 보든 여기서 출발한다.

| 파일 | 내용 |
|---|---|
| `BreslowClayton1993_...JASA.pdf` | **PQL의 원전.** epilepsy Model I~IV가 여기서 나왔고 이후 모든 논문의 비교 기준 |
| `Tierney_Kadane_1986_...JASA.pdf` | Laplace 근사 일반론 |
| `LiuPierce_1994_...Biometrika.pdf` | **적응 Gauss-Hermite 구적의 원전.** O&W가 로지스틱 B(μ,σ²) 계산에 인용 |
| `OrmerodWand_2010_...TAS.pdf` | **변분근사 입문.** 변분 쪽은 여기서 시작할 것 |

---

## 01_single_grouping — 단일 grouping factor

`u_i ∈ R^K`, 그룹끼리 독립. **우리가 구현·검증한 영역.**

| 파일 | 위치 |
|---|---|
| `OrmerodWand_2012_...AppendixA_JCGS.pdf` | **우리가 구현한 논문.** Appendix A(계산 상세)가 프리프린트에만 있음 |
| `PinheiroBates_1995_...JCGS.pdf` | 비선형혼합모형 우도 근사 |
| `McCulloch_1997_...JASA.pdf` | MCEM |
| `BreslowLin_1995_...Biometrika.pdf` | PQL 편향보정 (분산성분 1개) |
| `LinBreslow_1996_...JASA.pdf` | PQL 편향보정 (분산성분 여러 개) |
| `HallOrmerodWand_2011_...StatSinica.pdf` | **GVA 일관성.** 오차 O_p(m^{-1/2} + n^{-1}) |
| `HallPhamWandWang_2011_...AoS.pdf` | **점근정규성과 타당한 추론** |

**핵심**: Hall-Ormerod-Wand의 `O_p(m^{-1/2} + n^{-1})`이 우리가 시뮬레이션에서 측정한
"m은 무관, n에만 의존" 패턴의 이론적 근거다.

---

## 02_nested — 중첩 (계층)

`a_i`(학교) + `b_ij`(학급). grouping factor는 여럿이지만 **계층적으로 분해된다.**

```
L = Π_i ∫ [ Π_j ∫ (Π_k p(y_ijk | a_i, b_ij)) φ(b_ij) db_ij ] φ(a_i) da_i
```

| 파일 | 내용 |
|---|---|
| `PinheiroChao_2006_MultilevelAGQ_JCGS.pdf` | **중첩 AGQ 알고리즘. 이 폴더의 출발점** |
| `RabeHeskethSkrondalPickles_2005_...JEconometrics.pdf` | 중첩 ML 추정 (GLLAMM 계열) |
| `NolanMenictasWand_2020_...JMLR.pdf` | 상위수준 랜덤효과 streamlined 변분, 임의 층수 |
| `MaestriniBhaskaranWand_2024_...Biometrika.pdf` | 2수준 중첩 점근분산 명시형 |

읽는 순서: 2005/2006(계산) → 2020(변분) → 2024(이론)

**주의**: lme4·GLMMadaptive는 중첩에서 AGHQ를 지원하지 않는다(단일 grouping factor만).
알고리즘이 없어서가 아니라 구현이 없어서다 — Pinheiro & Chao가 2006년에 이미 냈다.

---

## 03_crossed — 교차

`a_s + b_k` (피험자 × 문항). **가능도가 분해되지 않는다.**
적분 차원 = 랜덤효과 총 개수(수백~수천). q^차원 격자는 존재 불가.

| 파일 | 접근 |
|---|---|
| `PapaspiliopoulosRobertsZanella_2020_...Biometrika.pdf` | MCMC 복잡도 이론. 평범한 Gibbs는 확장 불가, collapsed Gibbs는 가능 |
| `GhoshHastieOwen_Backfitting_...AoS.pdf` | Backfitting (빈도론) |
| `GVA_CompositeLikelihood_CrossedRE_StatSinica.pdf` | **GVA + 복합가능도. 우리 작업과 가장 가까움 — 먼저 읽을 것** |

---

## 04_general_highdim — 교차·중첩 공통

최근에는 둘을 "고차원 랜덤효과 모형" 한 틀로 묶어 다루는 흐름.

| 파일 | 내용 |
|---|---|
| `PapaspiliopoulosStumpfFetizonZanella_2023_CrossedAndNested_EJS.pdf` | 제목부터 "crossed **and** nested" |
| `GoplerudPapaspiliopoulosZanella_2025_...Biometrika.pdf` | **평균장이 고차원에서 사후분산을 심하게 과소평가**함을 보이고 처방(부분 인수분해) 제시 |

**우리 측정과의 연결**: 우리는 단일 grouping K≤3에서 GVA의 σ² 과소추정을
−0.0013 ~ −0.0026으로 실측했다. Goplerud 외는 이것이 고차원에서 심각해진다고 한다.
교차로 갈 거면 **평균장을 쓰지 않는 변분족**이 출발점이어야 한다.

---

## 이미 보유 (다른 폴더)

- McCulloch, Searle & Neuhaus, *Generalized Linear and Mixed Models* — `seminar/`
- Jiang (2007/2021), *Linear and GLMM and Their Applications* — `seminar/`
- **Jiang (2017), *Asymptotic Analysis of Mixed Effects Models: Theory, Applications,
  and Open Problems*** — `seminar/260306/` ← 부제의 "Open Problems"가 연구 방향 후보

---

## 읽는 순서 제안

1. `00_foundation/OrmerodWand_2010` → `01_single_grouping/OrmerodWand_2012`
2. `00_foundation/BreslowClayton1993` 정독 + epilepsy Model I~IV 재현
3. `01_single_grouping/HallOrmerodWand_2011` → `HallPhamWandWang_2011`
4. 방향 결정:
   - 중첩 → `02_nested/PinheiroChao_2006`
   - 교차 → `03_crossed/GVA_CompositeLikelihood` → `04_general_highdim/Goplerud2025`

---

## 관련 코드

`~/Documents/Statistic/seminar/blog/`

- `glmm_multi.R` — K차원 GVA + AGHQ (Poisson / Gamma / binomial), SQUAREM 가속
- `epilepsy_fit.R`, `epilepsy_benchmark.R` — O&W 7.2절 재현 및 벤치마크
- `run_K_simulation.R` — K=1,2,3 성능 비교

논문을 읽을 때마다 그 방법을 코드에 붙여 같은 데이터로 돌려보면 이해가 빠르다.
실제로 O&W 식 (3.1)과 코드를 한 줄씩 대조하다 offset 부호 버그를 잡았다.

## 아직 없는 것

| 논문 | DOI / 경로 |
|---|---|
| Menictas, Di Credico & Wand (2023), JCGS — 교차 LMM streamlined 변분 | https://pmc.ncbi.nlm.nih.gov/articles/PMC9983814/ (무료) |
| Gao & Owen — 교차설계 적률 계산, EJS | Project Euclid (무료, open access) |
