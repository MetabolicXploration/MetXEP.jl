const UNCONVERGED_STATUS = :unconverged
const CONVERGED_STATUS = :converged
const UNSET_STATUS = :unset

# is_converged(out::EPOut) = (out.status == CONVERGED_STATUS)

export convergence_status
convergence_status(epm::FluxEPModelT0) = state(epm, :status, nothing)

export didconverged
didconverged(epm::FluxEPModelT0) = (convergence_status(epm) === CONVERGED_STATUS)