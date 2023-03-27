struct HullLite{T<:Real}
    points::Vector{Point3{T}}
    simplices::Matrix{V_id}
    adj_simplices::Matrix{V_id}
    function HullLite{T}(points::Vector{Point3{T}}) where {T}
        span_points, span_ids = compute_span(points)
        if dot_cross(span_points...) < 0
            span_points[3], span_points[4] = span_points[4], span_points[3]
            span_ids[3], span_ids[4] = span_ids[4], span_ids[3]
        end
        a, b, c, d = span_ids
        simplices = Tuple{V_id,V_id,V_id}[(a, c, b), (a, b, d), (b, c, d), (c, a, d)]
        n = length(points)
        mark = zeros(V_id, n, n)
        for i in 1:n
            n_simplices = Tuple{V_id,V_id,V_id}[]
            bad_edges = Tuple{V_id,V_id}[]
            for (a, b, c) in simplices
                eps = 1e-12 # TODO: adjust?
                # dot_cross(points[[a, b, c, i]]...)
                if dot_cross(points[a], points[b], points[c], points[i]) > eps
                    # @assert a != i && b != i && c != i "a b c $(points[[a, b, c, i]]) $(dot_cross(points[[a, b, c, i]]...))"
                    append!(bad_edges, [(a, b), (b, c), (c, a)])
                else
                    push!(n_simplices, (a, b, c))
                end
            end
            for (a, b) in bad_edges
                mark[a, b] = i
            end
            for (a, b) in bad_edges
                if mark[b, a] != i
                    push!(n_simplices, (a, b, i))
                end
            end
            simplices = n_simplices
        end
        simplices_mat = reduce(hcat, collect.(simplices))
        new(points, simplices_mat, construct_adj_simplices_n2(n, simplices_mat))
    end
end

unscaled_volume_simplex(points::Vector{Point3{T}}) where {T} = det(reduce(hcat, points))
volume(p::HullLite{T}) where {T} = reduce(+, unscaled_volume_simplex(p.points[Δ]) for Δ in eachcol(p.simplices);
    init=zero(T)) / factorial(3)