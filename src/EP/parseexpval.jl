# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function parseexpval!(expval, siteflagave::BitArray{1}, siteflagvar::BitArray{1})

    return expave, expvar = _parseexpval!(expval, siteflagave, siteflagvar)
end

function _parseexpval!(expval::Tuple, siteflagave::BitArray{1},siteflagvar::BitArray{1})

    N = length(siteflagave)
    length(expval) == 3 || error("We expect a 3-Tuple here")

    mean, var, expsite = expval
    1 <= expsite <= N || error("expsite = $expsite not ∈ 1,...,$N")
    expave = Dict{Int,Float64}()
    expvar = Dict{Int,Float64}()

    # mean
    if mean !== nothing
        siteflagave[expsite] = false
        expave[expsite] = mean
    end

    # var
    if var !== nothing
        siteflagvar[expsite] = false
        expvar[expsite] = var
    end

    expave, expvar
end

function _parseexpval!(expval::Vector,siteflagave::BitArray{1},siteflagvar::BitArray{1})

    N = length(siteflagave)
    expave = Dict{Int,Float64}()
    expvar = Dict{Int,Float64}()
    for (mean, var, expsite) in expval

        1 <= expsite <= N || error("expsite = $expsite not ∈ 1,...,$N")

        # mean
        if mean !== nothing
            siteflagave[expsite] = false
            expave[expsite] = mean
        end

        # var
        if var !== nothing
            siteflagvar[expsite] = false
            expvar[expsite] =  var
        end
    end

    expave,expvar
end
_parseexpval!(nothing, siteflagave::BitArray{1}, siteflagvar::BitArray{1})=(Dict{Int,Float64}(),Dict{Int,Float64}())