const UNCONVERGED_STATUS = :unconverged
const CONVERGED_STATUS = :converged
const UNSET_STATUS = :unset

# is_converged(out::EPOut) = (out.status == CONVERGED_STATUS)

convergence_status(epm::FluxEPModelT0) = state(epm, :status, nothing)

didconverged(epm::FluxEPModelT0) = (convergence_status(epm) === CONVERGED_STATUS)