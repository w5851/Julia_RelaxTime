# PNJL Single-Point Benchmark
Generated at: 2025-12-02T17:41:40.552

| label | params | min (ms) | median (ms) | mean (ms) | max (ms) |
| --- | --- | ---: | ---: | ---: | ---: |
| Newton (default line search) | (method = :newton,) | 10.088 | 17.986 | 17.985 | 23.222 |
| Newton + BackTracking line search | (method = :newton, linesearch = LineSearches.BackTracking{Float64, Int64}(0.0001, 0.5, 0.1, 1000, 3, Inf, nothing)) | 18.644 | 22.194 | 23.062 | 31.235 |
| Trust-region | (method = :trust_region,) | 16.034 | 28.273 | 25.614 | 33.371 |
