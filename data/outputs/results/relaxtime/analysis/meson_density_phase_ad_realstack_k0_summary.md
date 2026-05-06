# Meson Phase AD Real-Stack k=0 Probe

## Goal

- Build an analysis-only generic micro-stack for a real `B0(k=0) -> polarization_with_width -> propagator` chain.
- Check float-sample consistency against the production implementation.
- Check whether `dδ/dω` through `Re/Im D` remains AD-friendly on this real-stack micro path.

## Sample

- `meson = pi`, `channel = P`, `k = 0`, `k0 = 0.2`, `gamma = 0.05`
- `m1 = 0.3`, `m2 = 0.3`, `T = 0.18`, `Phi = 0.2`, `PhiBar = 0.2`
- sample chosen below threshold so the k=0 B0 real-path integral stays non-singular in this probe

## Float Sample Consistency

- `Pi_local = 1.3860058807317754 - 0.001546477173796126im`
- `Pi_prod = 1.3860058670048467 - 0.0015464736875920697im`
- `|Pi_local - Pi_prod| = 1.416270420470015e-8`
- `D_local = -1.9908906808071385 - 0.012259840091593975im`
- `D_prod = -1.9908907899607462 - 0.012259813796671156im`
- `|D_local - D_prod| = 1.1227614636125587e-7`
- minimum sampled `|denominator|` in the local B0 quadrature: `0.040085914330948186`

## AD on Real-Stack Micro Path

- `Re D(omega0) = -1.9908906808071385`
- `Im D(omega0) = -0.012259840091593975`
- `d(Re D)/domega` via AD = `0.500315289097261`
- `d(Im D)/domega` via AD = `-0.057472340521868844`
- `d delta / domega` via component formula = `0.030414010304270874`
- `d delta / domega` via direct AD on atan(Im, Re) = `0.030414010304270878`
- `d delta / domega` via central finite difference = `0.030414010199208974`
- `|formula - direct| = 3.469446951953614e-18`
- `|formula - fd| = 1.0506190015191486e-10`

## Reading

1. A real micro-stack that already contains `B0 -> polarization_with_width -> propagator` can be made Dual-compatible without changing the production workflow contract.
2. For `d delta / domega`, `A` itself does not need to become Dual-compatible if it remains a precomputed kinematic constant with respect to `omega`.
3. The remaining production blocker is therefore narrower than 'the whole meson-density stack': it is mainly the typed signatures and helper implementations on the `B0` / polarization / propagator path.
4. This probe only covers `k=0`, isotropic, non-singular sampling; promotion to a production AD path would still need the `k>0` B0 path and interval quadrature helpers to accept generic parameter types.
