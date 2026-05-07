function dispersion_compensation(dx::T, dt, n, λ,) where {T}
    λ0 = λ
    c = 1 / n
    λ /= n
    k = 2π / λ
    w = 2 / dt * asin(dt / dx * c * sin(dx * k / 2))
    λ = 2π / w * c
    λ *= n

    a = 0.5
    @show λ0, λ, n
    (a * λ + (1 - a) * λ0) |> T
end