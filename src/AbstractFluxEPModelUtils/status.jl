const UNCONVERGED_STATUS = :unconverged
const CONVERGED_STATUS = :converged
const UNSET_STATUS = :unset

# is_converged(out::EPOut) = (out.status == CONVERGED_STATUS)

export didconverged
didconverged(epm) = state(epm, :status, nothing) === CONVERGED_STATUS