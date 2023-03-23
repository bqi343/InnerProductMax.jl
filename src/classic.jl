export InnerProductMaxClassic

"""
References: Edelsbrunner  & Maurer  (1985); 
Edelsbrunner (1987, Section 9.5.3); 
Computational Geometry in C 7.10.5
"""
struct InnerProductMaxClassic{T} <: AbstractInnerProductMax{T}
    function InnerProductMaxClassic{T}(hull::Hull{T}) where T
        throw("unimplemented")
    end
end

function query(ds::InnerProductMaxClassic{T}, p::Point3{T}) where T
    throw("unimplemented")
end