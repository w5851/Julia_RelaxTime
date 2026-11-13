# Shared input contract for the standalone and full /compute servers.
function _parse_compute_request(req::HTTP.Request)
    body = JSON3.read(req.body)
    body isa AbstractDict || throw(ArgumentError("Request body must be a JSON object"))
    number(key, default=nothing) = begin
        value = get(body, key, default)
        value isa Real && !(value isa Bool) || throw(ArgumentError("$key must be a number"))
        parsed = Float64(value)
        isfinite(parsed) || throw(ArgumentError("$key must be finite"))
        parsed
    end
    p1 = [number(:p1x), number(:p1y), number(:p1z)]
    p2 = [number(:p2x), number(:p2y), number(:p2z)]
    m1, m2, m3, m4 = number(:m1), number(:m2), number(:m3), number(:m4)
    all(>(0), (m1, m2, m3, m4)) || throw(ArgumentError("masses must be positive"))
    return (p1=p1, p2=p2, m1=m1, m2=m2, m3=m3, m4=m4,
        theta_star=number(:theta_star, π / 4), phi_star=number(:phi_star, π / 6))
end
