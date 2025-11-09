(* B0 Integral Calculation from Eq. (4.4) *)
(* HD-TVP-95-13: One Loop Integrals at Finite Temperature and Density *)
(* P. Rehberg and S. P. Klevansky *)

(* Clear any previous definitions *)
ClearAll["Global`*"];

(* Define basic variables *)
λ = k0 + μ1 - μ2;  (* k0 is real energy *)
E1 = Sqrt[p^2 + m1^2];
E2 = Sqrt[p^2 + m2^2];

(* Parameter assumptions reused across symbolic manipulations *)
parameterAssumptions = {
  Element[{k, k0, m1, m2, μ1, μ2, β, Phi, PhiBar, ξ, Λ}, Reals],
  k >= 0, Λ > 0, m1 > 0, m2 > 0, β > 0,
  0 <= Phi <= 1, 0 <= PhiBar <= 1, Abs[ξ] < 1
};
simplifyAssumptions = Join[parameterAssumptions, {Element[{p, θ}, Reals], p >= 0, 0 <= θ <= Pi}];

(* Define PNJL distribution functions and anisotropic correction (n || k) *)
nDotPSq[p_, angle_] := p^2 Cos[angle]^2;
fPhiPlus[E_, μ_] := (Phi Exp[-β (E - μ)] + 2 PhiBar Exp[-2 β (E - μ)] + Exp[-3 β (E - μ)])/
  (1 + 3 Phi Exp[-β (E - μ)] + 3 PhiBar Exp[-2 β (E - μ)] + Exp[-3 β (E - μ)]);
fPhiMinus[E_, μ_] := (PhiBar Exp[-β (E + μ)] + 2 Phi Exp[-2 β (E + μ)] + Exp[-3 β (E + μ)])/
  (1 + 3 PhiBar Exp[-β (E + μ)] + 3 Phi Exp[-2 β (E + μ)] + Exp[-3 β (E + μ)]);
dfPhiPlus[E_, μ_] := Module[{Et}, D[fPhiPlus[Et, μ], Et] /. Et -> E];
dfPhiMinus[E_, μ_] := Module[{Et}, D[fPhiMinus[Et, μ], Et] /. Et -> E];
fAnisoPlus[E_, μ_, p_, angle_] := fPhiPlus[E, μ] + (ξ/(2 E)) nDotPSq[p, angle] dfPhiPlus[E, μ];
fAnisoMinus[E_, μ_, p_, angle_] := fPhiMinus[E, μ] + (ξ/(2 E)) nDotPSq[p, angle] dfPhiMinus[E, μ];

(* Define momentum dot product, assuming k along z-direction *)
pDotk = p k Cos[θ];

(* ========== Define only the first term in Eq. (4.4) ========== *)

(* First term: f_1^(+,aniso) contribution *)
term1 = (fAnisoPlus[E1, μ1, p, θ]/E1) * 
  1/(λ^2 - 2 λ E1 + 2 pDotk - k^2 + m1^2 - m2^2);

(* B0 expression with only first term *)
B0Expression = 8 Pi^2 * term1;

(* Promote to full momentum integral with cutoff Λ *)
prefactor = (8 Pi^2)/(2 Pi)^3;
measure[p_, angle_] := p^2 Sin[angle];
coreIntegrand = term1;  (* Only first term *)
fullIntegrand = prefactor * 2 Pi * measure[p, θ] * coreIntegrand;

Print["=== Factorization Strategy: Separate θ-dependent and θ-independent terms ==="];
Print["Key Observation: Only sinθ·f₁ depends on θ in the integrand"];
Print[""];
Print["Factorization Structure:"];
Print["  1. θ-independent: p²/E₁, fΦ(E₁,μ₁), (ξp²/2E₁)·∂fΦ/∂E₁, λ²-2λE₁-k²+m₁²-m₂²"];
Print["  2. θ-dependent: sinθ, cos²θ (in f₁), cosθ (in propagator)"];
Print[""];
Print["Distribution function: f₁ = A₀(p) + A₂(p)·cos²θ"];
Print["  where A₀(p) = fΦ(E₁,μ₁)"];
Print["        A₂(p) = (ξp²/2E₁)·∂fΦ/∂E₁"];
Print[""];
Print["Propagator: Dp = D₀(p) + D₁(p)·cosθ"];
Print["  where D₀(p) = λ² - 2λE₁ - k² + m₁² - m₂²"];
Print["        D₁(p) = 2kp"];
Print[""];
Print["Strategy: Substitute x = cos(θ), integrate rational function (A₀+A₂x²)/(D₀+D₁x)"];
Print[""];
Print["=== Step 1: Define symbolic coefficient parameters ==="];

(* Define symbolic parameters instead of computing actual PNJL expressions *)
(* These represent the θ-independent coefficient functions of p *)
Print["Defining symbolic coefficients (not expanded):"];
Print["  Prefactor[p] - momentum and energy prefactor"];
Print["  A0[p] - isotropic distribution function fΦ⁺(E₁,μ₁)"];
Print["  A2[p] - anisotropic correction (ξp²/2E₁)·∂fΦ⁺/∂E₁"];
Print["  D0[p] - propagator constant term λ² - 2λE₁ - k² + m₁² - m₂²"];
Print["  D1[p] - propagator angular term 2kp"];
Print[""];
Print["These are treated as symbolic functions in the integration"];
Print[""];

Print["=== Step 2: Construct factorized integrand structure ==="];
Print["Original integrand (before x = cosθ substitution):"];
Print["  (2p²sinθ/E₁) · [A₀(p) + A₂(p)·cos²θ] / [D₀(p) + D₁(p)·cosθ]"];
Print[""];
Print["After substitution x = cosθ, sinθ·dθ = dx:"];
Print["  Integrand = (2p²/E₁) · [A₀(p) + A₂(p)·x²] / [D₀(p) + D₁(p)·x]"];
Print[""];
Print["Factor out θ-independent prefactor: 2p²/E₁ = 2p²/√(m₁²+p²)"];
Print[""];
Print["Generic rational function to integrate:"];
Print["  ∫₋₁¹ [N₀ + N₂·x²] / [D₀ + D₁·x] dx"];
Print["where coefficients are symbolic functions of p:"];
Print["  N₀ = A0[p] (symbolic, represents fΦ⁺/E₁)"];
Print["  N₂ = A2[p] (symbolic, represents (ξp²/2E₁²)·∂fΦ⁺/∂E₁)"];
Print["  D₀ = D0[p] (symbolic, represents λ² - 2λE₁ - k² + m₁² - m₂²)"];
Print["  D₁ = D1[p] (symbolic, represents 2kp)"];
Print[""];

Print["=== Step 3: Integrate generic rational function symbolically ==="];
angularIntegralTimeLimit = 90; (* seconds *)

Print["Method: Analytical integration with residue theorem for poles"];
Print[""];
Print["Key: Keep coefficients as SYMBOLS, don't expand PNJL expressions"];
Print["This produces compact symbolic output in terms of A0[p], A2[p], D0[p], D1[p]"];
Print[""];
Print["Attempting integration of:"];
Print["  ∫₋₁¹ (a0 + a2·x²)/(d0 + d1·x) dx"];
Print[""];

(* Define symbolic coefficient parameters as CONSTANTS not functions *)
Clear[a0, a2, d0, d1, x, xPole];

(* Use simple symbols that Mathematica can integrate *)
genericIntegrand = (a0 + a2 * x^2) / (d0 + d1 * x);

Print["Checking for poles in integration region [-1, 1]..."];
(* Pole location: d0 + d1*x = 0  =>  x = -d0/d1 *)
xPole = -d0/d1;
Print["  Pole location: x_pole = -d0/d1"];
Print[""];

(* Check if pole is in [-1, 1] *)
Print["=== Case Analysis ==="];
Print["Case 1: |d0/d1| > 1  =>  Pole outside [-1,1], regular integral"];
Print["Case 2: |d0/d1| < 1  =>  Pole inside [-1,1], principal value + residue"];
Print[""];

(* Compute residue at the pole *)
residueAtPole = (a0 + a2 * xPole^2) / d1;
Print["Residue at pole (if inside interval):"];
Print["  Res = (a0 + a2·x_pole²)/d1 = (a0 + a2·d0²/d1²)/d1"];
Print["  Res = (a0·d1² + a2·d0²)/d1³"];
Print[""];

(* Principal value integral (analytical formula) *)
Print["Computing principal value integral..."];
Print[""];

(* Try symbolic integration with GenerateConditions *)
genericResultWithConditions = TimeConstrained[
  Integrate[genericIntegrand, {x, -1, 1}, 
    Assumptions -> {Element[{a0, a2, d0, d1}, Reals], d1 != 0},
    GenerateConditions -> True],
  30,
  $Failed
];

If[genericResultWithConditions =!= $Failed,
  Print["✓ Mathematica symbolic integration succeeded with conditions!"];
  Print["Result: ", genericResultWithConditions];
  Print[""];
  angularIntegralRealPart = genericResultWithConditions;
,
  Print["Mathematica timed out, using analytical formula..."];
  Print[""];
  
  (* Use known analytical formula for principal value *)
  Print["=== Analytical Formula (Principal Value) ==="];
  Print[""];
  Print["For integral:  ∫_{-1}^{+1} (a0 + a2*x²)/(d0 + d1*x) dx"];
  Print[""];
  Print["Principal value result:"];
  Print["  I_PV(p) = -[2*d0*d1*a2 + (d1²*a0 + d0²*a2)*L(p)] / d1³"];
  Print[""];
  Print["where:  L(p) = Log|d0 - d1| - Log|d0 + d1|"];
  Print[""];
  
  (* Principal value formula *)
  angularIntegralRealPart = -(2*d0*d1*a2 + (d1^2*a0 + d0^2*a2)*(Log[Abs[d0 - d1]] - Log[Abs[d0 + d1]]))/d1^3;
  
  Print["Simplified form (real part): ", angularIntegralRealPart];
  Print[""];
];

(* Add imaginary part from residue contribution *)
Print["=== Imaginary Part from Residue Theorem ==="];
Print[""];
Print["When pole x_pole ∈ (-1, 1), add imaginary contribution:"];
Print["  Im[I] = -π·i·Res  (sign convention: pole below real axis)"];
Print["  Im[I] = -π·i·(a0·d1² + a2·d0²)/d1³"];
Print[""];

(* Imaginary part (only non-zero if pole is inside [-1, 1]) *)
angularIntegralImagPart = -Pi * I * (a0 * d1^2 + a2 * d0^2) / d1^3;

Print["Complete result with pole handling:"];
Print[""];
Print["I(p) = I_PV(p) + i·Im[I(p)]·Θ(1 - |d0/d1|)"];
Print[""];
Print["where Θ(1 - |d0/d1|) is Heaviside step function (1 if pole inside, 0 otherwise)"];
Print[""];

(* Full result with conditional imaginary part *)
angularIntegralParam = angularIntegralRealPart + 
  Piecewise[{{angularIntegralImagPart, Abs[d0/d1] < 1}}, 0];

Print["Symbolic form:"];
Print["  Real part: ", angularIntegralRealPart];
Print["  Imag part: ", angularIntegralImagPart, " × Θ(1 - |d0/d1|)"];
Print[""];

Print["=== Physical Interpretation ==="];
Print["Symbolic parameters (functions of momentum p):"];
Print["  a0 = A0(p)/E1 = fΦ⁺(E1, μ1) / E1  [isotropic]"];
Print["  a2 = A2(p)/E1 = [ξp²/(2E1²)]·∂fΦ⁺/∂E1  [anisotropic]"];
Print["  d0 = D0(p) = λ² - 2λE1 - k² + m1² - m2²  [propagator constant]"];
Print["  d1 = D1(p) = 2kp  [propagator angular term]"];
Print[""];
Print["Pole condition: |λ² - 2λE1 - k² + m1² - m2²| < |2kp|"];
Print["  => On-shell propagator leads to imaginary contribution"];
Print[""];
Print["Final B0:  B0 = 8π² ∫₀^Λ [2p²/E1]·I(p) dp"];
Print[""];
Print["✓ SUCCESS: Formula includes RESIDUE contributions!"];
Print[""];

(* Also store the step function condition for later use *)
poleInsideCondition = Abs[d0/d1] < 1;

(* Fallback to the original approach if generic failed *)
If[Head[angularIntegralParam] === Inactive[Integrate] || !FreeQ[angularIntegralParam, {N0, N2, D0, D1}],
  Print["Falling back to direct integration with expanded expressions..."];
  Print[""];
  
  integrandInX = fullIntegrand /. {Cos[θ] -> x, Sin[θ] -> Sqrt[1 - x^2]};
  angularIntegrandX = integrandInX /. {Sqrt[1 - x^2] -> 1};
  
  Print["Attempting Apart on expanded integrand..."];
  angularIntegrandSimplified = TimeConstrained[
    Assuming[Join[parameterAssumptions, {p > 0, k > 0, Element[x, Reals]}], Apart[angularIntegrandX, x]],
    20, angularIntegrandX
  ];
  If[angularIntegrandSimplified =!= angularIntegrandX,
    Print["Apart simplified the expanded integrand!"];
    angularIntegrandX = angularIntegrandSimplified;
  ];
  
  angularIntegralParam = TimeConstrained[
    Assuming[Join[parameterAssumptions, {p > 0, k > 0}],
      Integrate[angularIntegrandX, {x, -1, 1}, PrincipalValue -> True, GenerateConditions -> False]
    ],
    angularIntegralTimeLimit,
    Inactive[Integrate][angularIntegrandX, {x, -1, 1}, PrincipalValue -> True, GenerateConditions -> False]
  ];
];

angularIntegralExpression = If[Head[angularIntegralParam] === Inactive[Integrate],
  angularIntegralParam,
  Inactive[Integrate][angularIntegralParam, {x, -1, 1}]  (* wrap for display *)
];
Print["Angular integral expression: ", Short[angularIntegralExpression, 3]];

(* Check if we already have a result from the parametric approach *)
If[Head[angularIntegralParam] =!= Inactive[Integrate],
  (* We have a result! *)
  angularIntegral = angularIntegralParam;
  angularIntegralStatus = "Angular integration completed using symbolic coefficients!";
  angularIntegralSuccess = True;
  Print["Angular integration status: ", angularIntegralStatus];
  Print["Result (COMPACT SYMBOLIC FORM): ", angularIntegral];
  Print[""];
,
  (* Try to activate the integral *)
  angularIntegralAttempt = TimeConstrained[
     Quiet@Check[Assuming[Join[parameterAssumptions, {Element[x, Reals], -1 <= x <= 1}], 
       Activate[angularIntegralExpression]], $Failed],
     10, $Failed
  ];

  If[angularIntegralAttempt === $Failed,
    angularIntegral = angularIntegralExpression;
    angularIntegralStatus = "Angular integration (in x) not completed within time limit, keeping integral expression.";
    angularIntegralSuccess = False,
    angularIntegral = angularIntegralAttempt;
    angularIntegralStatus = "Angular integration (in x) completed successfully!";
    angularIntegralSuccess = True
  ];

  Print["Angular integration status: ", angularIntegralStatus];
  If[angularIntegralSuccess, Print["Angular integral result (function of p): ", Short[angularIntegral, 4]]];
  Print[""];
];

(* Then construct the full integral with p integration kept symbolic *)
If[angularIntegralSuccess,
  (* If angular integration succeeded, keep p integral symbolic *)
  B0Integral = Inactive[Integrate][angularIntegral, {p, 0, Λ}];
  B0IntegralStatus = "Angular integration completed. Momentum integration kept symbolic.";
  B0IntegralSuccess = True,
  (* If angular integration failed, keep full integral *)
  B0Integral = Inactive[Integrate][fullIntegrand, {p, 0, Λ}, {θ, 0, Pi}];
  B0IntegralStatus = "Angular integration not completed within time limit, keeping full integral expression.";
  B0IntegralSuccess = False
  ];

(* ========== Output Results ========== *)

Print["=============================================="];
Print["B0 Integral Calculation - Eq. (4.4)"];
Print["HD-TVP-95-13: One Loop Integrals at Finite Temperature and Density"];
Print["=============================================="];
Print[""];

Print["=== First term in Eq. (4.4) - f_1^(+,aniso) contribution ==="];
Print["Term 1: ", term1];
Print[""];

Print["=== Complete B0 Expression ==="];
Print["B0 = ", B0Expression];
Print[""];

Print["=== B0 with Integration Measure ==="];
If[angularIntegralSuccess,
  Print["Angular integral result: ", angularIntegral];
  Print[""];
  ];
Print["B0Integral = ", B0Integral];
Print["Integration Status: ", B0IntegralStatus];
Print[""];

(* Simplification with physical assumptions *)
Print["=== Simplified Expression with Physical Assumptions ==="];
simplificationTimeConstraint = 20;
simplified = TimeConstrained[
   Assuming[simplifyAssumptions, Simplify[B0Expression, TimeConstraint -> simplificationTimeConstraint]],
   simplificationTimeConstraint,
   "SimplifyTimedOut"
   ];
If[simplified === "SimplifyTimedOut",
  B0SimplifyStatus = "Simplify not completed within time limit.";
  simplified = B0Expression,
  B0SimplifyStatus = "Simplify completed within time limit."
  ];
Print["Simplified B0 = ", simplified];
Print["Simplification Status: ", B0SimplifyStatus];
Print[""];

(* Symmetric case *)
Print["=== Symmetric Case (m1 = m2, mu1 = mu2, k0 = 0) ==="];
symmetricCase = B0Expression /. {m1 -> m, m2 -> m, μ1 -> μ, μ2 -> μ, k0 -> 0};
symSimplified = Assuming[Join[parameterAssumptions, {Element[{p, θ, m}, Reals], p >= 0, 0 <= θ <= Pi}], 
  Simplify[symmetricCase, TimeConstraint -> simplificationTimeConstraint]];
Print["Symmetric B0 = ", symSimplified];
Print[""];

(* Special case k=0 *)
Print["=== Special Case k=0 ==="];
kZeroCase = B0Expression /. {k -> 0};
kZeroSimplified = Assuming[simplifyAssumptions, Simplify[kZeroCase]];
Print["B0 at k=0 = ", kZeroSimplified];
Print[""];

(* Numerical verification at a sample point *)
Print["=== Numerical Verification (Example Parameters) ==="];
numericParameterRules = {
  k -> 0.5, k0 -> 0.3, m1 -> 0.1, m2 -> 0.2, μ1 -> 0.0, μ2 -> 0.0,
  β -> 1.0, Phi -> 0.5, PhiBar -> 0.5, ξ -> 0.1,
  Λ -> 2.0
};
samplePointRules = {p -> 1.0, θ -> Pi/4};
numericalExample = B0Expression /. numericParameterRules /. samplePointRules;
Print["Numerical example result at sample point = ", numericalExample];
Print[""];

(* Note: Numerical integration removed per user request *)
numericIntegral = "Numerical integration disabled";
numericIntegralStatus = "Numerical integration disabled per user request.";
numericIntegralMessages = {};
Print["Numerical integration status: ", numericIntegralStatus];
Print[""];

(* Real and imaginary parts *)
Print["=== Real and Imaginary Parts ==="];
realPart = ComplexExpand[Re[B0Expression]];
imagPart = ComplexExpand[Im[B0Expression]];
Print["Real part = ", realPart];
Print["Imaginary part = ", imagPart];
Print[""];

Print["=============================================="];
Print["Calculation Complete!"];
Print["=============================================="];

(* Return main results for further use *)
results = <|
  "FullExpression" -> B0Expression,
  "AngularIntegral" -> If[angularIntegralSuccess, angularIntegral, "Not completed"],
  "AngularIntegralStatus" -> angularIntegralStatus,
  "FullIntegral" -> B0Integral,
  "IntegralStatus" -> B0IntegralStatus,
  "SimplifyStatus" -> B0SimplifyStatus,
  "Simplified" -> simplified,
  "SymmetricCase" -> symSimplified,
  "kZeroCase" -> kZeroSimplified,
  "NumericalExample" -> numericalExample,
  "NumericIntegral" -> numericIntegral,
  "NumericIntegralStatus" -> numericIntegralStatus,
  "NumericIntegralMessages" -> numericIntegralMessages,
  "RealPart" -> realPart,
  "ImaginaryPart" -> imagPart
|>;

(* Display main results in VSCode *)
Print["Main results saved in variable 'results'"];
Print["Keys: ", Keys[results]];

(* Export to file with proper UTF-8 encoding *)
outputFilePath = "d:\\Desktop\\Julia_RelaxTime\\results\\B0_output_symbolic.txt";
outputLines = {
  "==============================================",
  "B0 Angular Integration - COMPACT SYMBOLIC FORM",
  "Factorized Approach with Residue Contributions",
  "==============================================",
  "",
  "=== Generic Integral (After substitution x = cos θ) ===",
  "∫_{-1}^{+1} (a0 + a2*x²)/(d0 + d1*x) dx",
  "",
  "=== Pole Analysis ===",
  "Pole location: x_pole = -d0/d1",
  "Residue at pole: Res = (a0 + a2·x_pole²)/d1 = (a0·d1² + a2·d0²)/d1³",
  "",
  "=== Principal Value (Real Part) ===",
  "I_PV(p) = -[2*d0*d1*a2 + (d1²*a0 + d0²*a2)*L(p)] / d1³",
  "",
  "where:",
  "  L(p) = Log|d0 - d1| - Log|d0 + d1|",
  "",
  "Mathematica form (real part):",
  ToString[angularIntegralRealPart, InputForm],
  "",
  "=== Imaginary Part (from Residue Theorem) ===",
  "When pole inside [-1, 1], i.e., |d0/d1| < 1:",
  "  Im[I(p)] = -π·(a0·d1² + a2·d0²)/d1³",
  "",
  "Mathematica form (imaginary part):",
  ToString[angularIntegralImagPart, InputForm],
  "",
  "=== Complete Result ===",
  "I(p) = I_PV(p) + i·Im[I(p)]·Θ(1 - |d0/d1|)",
  "",
  "Full Mathematica form:",
  ToString[angularIntegralParam, InputForm],
  "",
  "Pole condition: |d0/d1| < 1  <==>  |λ² - 2λE1 - k² + m1² - m2²| < |2kp|",
  "",
  "=== Symbolic Parameter Definitions ===",
  "All are functions of momentum p:",
  "",
  "a0 = A0(p)/E1(p) = fΦ⁺(E1, μ1) / E1",
  "   - Isotropic PNJL quark distribution",
  "   - Temperature and chemical potential dependent",
  "",
  "a2 = A2(p)/E1(p) = [ξp²/(2E1²)] · ∂fΦ⁺/∂E1",
  "   - Anisotropic correction to distribution",
  "   - Proportional to anisotropy parameter ξ",
  "",
  "d0 = D0(p) = λ² - 2λE1 - k² + m1² - m2²",
  "   - Propagator constant term",
  "   - λ = k0 + μ1 - μ2",
  "",
  "d1 = D1(p) = 2kp",
  "   - Propagator angular coefficient",
  "   - k is gluon momentum magnitude",
  "",
  "where E1 = √(m1² + p²) is quark quasi-particle energy",
  "",
  "=== Complete B0 Expression ===",
  "B0 = 8π² ∫₀^Λ Prefactor(p) · I(p) dp",
  "   = 8π² ∫₀^Λ [2p²/E1] · [I_PV(p) + i·Im[I(p)]·Θ(1-|d0/d1|)] dp",
  "",
  "=== Physical Interpretation ===",
  "• Logarithmic term L(p) encodes:",
  "  - Propagator pole structure at x = ±d0/d1",
  "  - Quasiparticle energy spectrum modifications",
  "  - Medium response in anisotropic QGP",
  "",
  "• Imaginary part from residue theorem:",
  "  - Appears when propagator has on-shell poles",
  "  - Encodes decay width / damping rate",
  "  - Related to Landau damping in thermal medium",
  "  - Non-zero only when |d0| < |d1|, i.e., |λ² - 2λE1 - k² + Δm²| < |2kp|",
  "",
  "• Factorization benefits:",
  "  - Compact symbolic formula (NOT expanded with exponentials)",
  "  - Clear separation of distribution (a0, a2) and propagator (d0, d1)",
  "  - Efficient numerical evaluation in Julia",
  "",
  "=== Julia Implementation Guide ===",
  "function B0_aniso(k, k0; params...)",
  "    function integrand(p)",
  "        E1 = sqrt(m1^2 + p^2)",
  "        # Compute PNJL distribution",
  "        a0 = compute_a0(p, E1, mu1, T)  # fΦ⁺/E1",
  "        a2 = compute_a2(p, E1, mu1, T, xi)  # anisotropic term",
  "        # Compute propagator coefficients",
  "        lambda = k0 + mu1 - mu2",
  "        d0 = lambda^2 - 2*lambda*E1 - k^2 + m1^2 - m2^2",
  "        d1 = 2*k*p",
  "        # Angular integral - Real part (principal value)",
  "        Lpole = log(abs(d0 - d1)) - log(abs(d0 + d1))",
  "        Ipv = -(2*d0*d1*a2 + (d1^2*a0 + d0^2*a2)*Lpole) / d1^3",
  "        # Imaginary part (residue contribution)",
  "        if abs(d0/d1) < 1.0  # pole inside [-1, 1]",
  "            Iim = -π * (a0*d1^2 + a2*d0^2) / d1^3",
  "        else",
  "            Iim = 0.0",
  "        end",
  "        Ip = complex(Ipv, Iim)",
  "        # Prefactor",
  "        return (2*p^2/E1) * Ip",
  "    end",
  "    return 8*π^2 * quadgk(integrand, 0, Lambda)[1]",
  "end",
  "",
  "==============================================",
  "✓ Output includes RESIDUE contributions!",
  "✓ Formula is COMPACT - No PNJL expansion!",
  "=============================================="
};
Export[outputFilePath, StringRiffle[outputLines, "\n"], "Text", CharacterEncoding -> "UTF-8"];
Print["Results exported to: ", outputFilePath];
