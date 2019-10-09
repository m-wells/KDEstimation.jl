#TODO: currently scaling rules of thumb from normal. Should derive them directly.

@inline function rule_of_thumb1(::Type{Normal}, X::AbstractVector)
    n,σ = length(X),std(X)
    return 1.06*σ*n^(-1/5)
end

@inline function rule_of_thumb2(::Type{Normal}, X::AbstractVector)
    n,σ = length(X),std(X)
    q1,q3 = quantile(X,[0.25,0.75])
    iqr = q3 - q1
    return 1.06*min(σ,iqr/1.34)*n^(-1/5)
end

function rule_of_thumb1(::Type{D}, X::AbstractVector) where {D<:Distribution}
    h1 = rule_of_thumb1(Normal, X)
    return canonical_transform(h1, Normal, D)
end

function rule_of_thumb2(::Type{D}, X::AbstractVector) where {D<:Distribution}
    h1 = rule_of_thumb2(Normal, X)
    return canonical_transform(h1, Normal, D)
end

rule_of_thumb1(X::AbstractVector) = rule_of_thumb1(Normal, X)
rule_of_thumb2(X::AbstractVector) = rule_of_thumb2(Normal, X)
