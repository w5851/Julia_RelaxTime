## Evidence boundary

- Input is the immutable v6 postprocessed package; input files are copied under `tables/v6_*`.
- Native rows are not overwritten. Derived rows carry source keys, source layer, calculation SHA and v6 manifest hash.
- The CEP bracket is retained as an interval. The estimated midpoint is a boundary estimate for the derived crossover display, not a strict CEP result.
- Maxwell rows are copied only; no Maxwell interpolation, candidate selection or geometry closure is performed.
- No ordinary extrapolation crosses a missing native support interval, a neighbouring-xi common-support boundary, or the physical `mu_q > mu_CEP` region.
- v7 is a diagnostic/derived candidate and does not promote phase-reference or unlock transport.
