# Endpoint-local geometry contract feasibility v2

verdict: `feasible_candidate`。本包只做 solver-free replay，不调用 PNJL，不修改 production/reference/transport。

- standard run: `30999377195`
- deep run: `31002704845`
- calculation SHA: `ceec2295c5c9250a3fcd45c0eceae9a6c35f4335`
- route: complete Stage-B curve + midpoint in the active left Maxwell bracket
- right-side rule: actual Stage-B outer-branch bracket; no `mu > mu_spinodal,max` requirement
- certificates observed: `endpoint_limited_first_order, endpoint_local_geometry_first_order`
- coverage: `targeted_18_plus_approved_required_three_deep`; full 24-anchor shadow remains mandatory

Deep statuses are retained as a post-route gate only; they are not used to choose support or midpoint locations. A `feasible_candidate` here authorizes creation of the focused production PR, not physical promotion.
