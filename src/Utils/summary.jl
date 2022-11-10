# import MetNets._print_state_head

# _print_state_head(io::IO, out::EPOut) = printstyled(io, 
#     string(
#         " EPOut: ", 
#         "status=", status(out), ", ",
#         "curriter=", curriter(out), ", ",
#         "max_abs_av=", round(maximum(abs, av(out)); sigdigits = 4),
#         "\n"
#     ),
#     color = MetNets.INFO_COLOR
# )