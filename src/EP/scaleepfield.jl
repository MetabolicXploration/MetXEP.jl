# Code derived from metabolicEP (https://github.com/anna-pa-m/Metabolic-EP)

function scaleepfield!(scalefact, epfields, vals...)
    @extract epfields : μ s av va 
    rmul!(μ, scalefact)
    rmul!(s, scalefact^2)
    rmul!(av, scalefact)
    rmul!(va, scalefact^2)

    for val in vals
        rmul!(val, scalefact)
    end
    
end