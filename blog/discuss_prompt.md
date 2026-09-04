# 논의 요청: Hall, Ormerod & Wand (2011)이 실제로 증명한 것

## 배경

Gamma random-intercept GLMM에서 GVA(Gaussian Variational Approximation)와
AGHQ(Adaptive Gauss-Hermite Quadrature)를 비교하는 시뮬레이션 연구를 진행 중이고,
결과를 `posts/study/Gamma_GLMM.qmd` 에 정리해 왔다.

그런데 방금 시뮬레이션 코드에서 심각한 버그를 발견해서, 지금까지의 Gamma 결과를
전부 다시 봐야 하는 상황이다. 그 과정에서 **참조 논문을 어떻게 이해하고 있었는지
자체가 흔들려서**, 논문의 실제 주장을 다시 확인하고 싶다.

## 발견한 버그

`fit_gva()` 의 M-step에서 Gamma 분기의 offset 부호가 뒤집혀 있었다.

ELBO의 $\beta$ 의존 항은

$$\sum_j\left[-\nu\,\eta_{ij}^{fix} - \nu Y_{ij}\exp\!\big(-(\eta_{ij}^{fix} + \mu_i - \lambda_i/2)\big)\right]$$

이므로 offset-GLM에 넘겨야 할 offset은 $+\mu_i - \lambda_i/2$ 인데,
코드는 $-(\mu_i - \lambda_i/2)$ 를 넘기고 있었다.

- Poisson 분기(`mu + lam/2`)는 원래 맞았다. 그래서 **Poisson 대조군이 깨끗하게
  통과했고**, 그것 때문에 코드가 검증됐다고 판단해 왔다. 버그가 대조군이 건드리지
  않는 쪽에만 있었다.
- 이 offset은 그룹 $i$ 안에서 상수이므로 절편만 밀고 기울기는 건드리지 않는다.
  그동안 관찰한 "$\beta_0,\sigma^2$는 크게 편향, $\beta_1$은 무편향" 패턴이 정확히
  이 구조였다.

## 버그 수정 후 바뀐 것

수정 후 Poisson은 소수점 12자리까지 완전히 불변임을 확인했다(회귀테스트 통과).
Gamma는 결과가 뒤집혔다:

| 조건 | GVA $\beta_0$ 편향 | AGHQ $\beta_0$ 편향 | GVA $\sigma^2$ 편향 | AGHQ $\sigma^2$ 편향 |
|---|---|---|---|---|
| $\beta_0{=}1,\sigma^2{=}1,n{=}5$ | +0.015 | +0.015 | −0.001 | +0.003 |
| $\beta_0{=}1,\sigma^2{=}1,n{=}25$ | −0.010 | −0.010 | −0.033 | −0.034 |
| $\beta_0{=}3,\sigma^2{=}4,n{=}25$ | −0.061 | −0.061 | −0.194 | −0.193 |

즉 **GVA와 AGHQ의 편향이 소수점 3자리까지 사실상 동일**하다.
기존에 보고했던 "GVA $\hat\beta_0$가 참값의 약 2배"는 버그 산물로 보인다.

추가로, ELBO를 직접 최대화한 결과와 대조해 보니 수정된 `fit_gva`가 ELBO 최적점에
도달함을 확인했다(다만 수렴에 수백~수천 회 반복이 필요해서 기존 `max_iter=150`은
한참 모자랐다. 이것도 별개의 오염 요인이었다).

시간 면에서는 GVA가 수렴까지 돌리면 AGHQ보다 **오히려 느리다**(n이 클 때 최대 5배).

## 핵심 질문 — 블로그 본문의 자체 모순

`Gamma_GLMM.qmd` 에 이렇게 써 놨는데, 두 문장이 서로 충돌한다:

> (190행) Hall, Ormerod & Wand (2011, *Statistica Sinica*)는 Poisson
> random-intercept GLMM에서 그룹당 관측치 수 $n$을 고정하고 그룹 수 $m\to\infty$로
> 보내면 GVA 추정량이 **비일관적(inconsistent)**임을 증명했다.

> (194행) Poisson은 기존 이론과 비교할 수 있는 대조군이다: $m$이 커지면 편향이
> **소멸해야** 한다.

논문이 "n 고정, m→∞에서 비일관"을 증명했다면, Poisson 대조군에서 편향이 0으로
나온 것은 이론을 **확인**하는 게 아니라 **반박**하는 결과여야 한다. 그런데 글에서는
편향 ≈ 0을 "기존 이론과 일치한다"고 서술했다. 둘 중 하나(또는 둘 다)가 틀렸다.

실제로 우리 Poisson 시뮬레이션은 $m=3000, 10000$ 에서 $\hat\beta_0$ 편향이
+0.0003, +0.0024 로 거의 0이었다.

## 알고 싶은 것

논문(Hall, P., Ormerod, J.T., Wand, M.P. (2011), *Statistica Sinica* 21, 369-389,
"Theory of Gaussian variational approximation for a Poisson mixed model")을 찾아서
다음을 확인해 주면 좋겠다:

1. **논문이 실제로 증명한 주 정리가 무엇인가?** 일관성인가 비일관성인가?
   어떤 점근 체제(asymptotic regime)에서인가 — $n$ 고정 $m\to\infty$인가,
   $m,n$ 동시 발산인가, 아니면 둘의 상대 속도에 조건이 붙는가?

2. **편향의 차수(order)가 어떻게 되는가?** $O(1/n)$ 인가 다른 차수인가.
   $m$은 편향에 영향을 주는가, 아니면 분산에만 영향을 주는가?

3. 위의 블로그 두 문장 중 어느 쪽이 논문에 부합하는가? 내가 논문 주장을
   반대로 기억하고 있었을 가능성이 높은데, 확인해 달라.

4. **Poisson 특화인가?** 논문의 논증이 Poisson의 어떤 성질에 의존하는지,
   Gamma 같은 다른 지수족으로 확장 가능한 형태인지.

5. 논문에 **계산 비용(GVA vs quadrature)** 에 대한 언급이 있는가?
   우리 결과에서는 수렴까지 돌린 GVA가 AGHQ보다 느렸는데, 논문이 GVA의 장점으로
   무엇을 내세우는지 알고 싶다.

6. 수정된 우리 결과(GVA ≈ AGHQ, 둘 다 거의 무편향)가 논문의 이론과
   **모순되는가 부합하는가?**

## 참고

- 논문을 웹에서 찾을 수 있으면 실제 내용을 확인해서 답해 달라.
  기억에만 의존한 요약이면 그렇다고 명시해 주면 좋겠다.
- 관련 후속 논문(예: Ormerod & Wand의 다른 variational 논문, 또는 GVA/GVB의
  일관성을 다룬 이후 문헌)이 있으면 같이 알려주면 도움이 된다.
- 코드는 `pilot_sim.R`, 비교 시뮬레이션은 `gva_vs_aghq.R` 에 있다.
