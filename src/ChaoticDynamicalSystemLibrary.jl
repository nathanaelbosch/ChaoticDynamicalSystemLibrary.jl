module ChaoticDynamicalSystemLibrary

using Base: nothing_sentinel
using JSON, Markdown
using LinearAlgebra
using SciMLBase
using ComponentArrays
using SimpleUnPack

vecvec2mat(vv) = hcat(vv...)' # transpose is because of the way the data is stored in the JSON file

function format_parameters(params::Dict)
    if isempty(params)
        return nothing
    end

    intermediate_dict = Dict()
    for (k, v) in params
        val = if (v isa Vector) && (v[1] isa Vector)
            vecvec2mat(v)
        else
            v
        end
        intermediate_dict[Symbol(k)] = val
    end

    return ComponentArray(intermediate_dict)
end


function make_docstring(f)
    data = ATTRACTOR_DATA[string(f)]
    header = """
           $f()

       $(data["description"])
    """

    stats = """

       ## Stats
    """
    "correlation_dimenension" in keys(data) && (stats *= """
       - Correlation dimension: `$(data["correlation_dimension"])`
    """)
    "embedding_dimension" in keys(data) && (stats *= """
       - Embedding dimension: `$(data["embedding_dimension"])`
    """)
    "hamiltonian" in keys(data) && (stats *= """
       - Hamiltonian: `$(data["hamiltonian"])`
    """)
    "kaplan_yorke_dimension" in keys(data) && (stats *= """
       - Kaplan-Yorke dimension: `$(data["kaplan_yorke_dimension"])`
    """)
    "lyapunov_spectrum_estimated" in keys(data) && (stats *= """
       - Lyapunov spectrum (estimated): `$(data["lyapunov_spectrum_estimated"])`
    """)
    "maximum_lyapunov_estimated" in keys(data) && (stats *= """
       - Maximum Lyapunov exponent (estimated): `$(data["maximum_lyapunov_estimated"])`
    """)
    "multiscale_entropy" in keys(data) && (stats *= """
       - Multiscale Entropy: `$(data["multiscale_entropy"])`
    """)
    "nonautonomous" in keys(data) && (stats *= """
       - Non-autonomous: `$(data["nonautonomous"])`
    """)
    "period" in keys(data) && (stats *= """
       - Period: `$(data["period"])`
    """)
    "pesin_entropy" in keys(data) && (stats *= """
       - Pesin entropy: `$(data["pesin_entropy"])`
    """)

    citation = """

       ## References
       $(data["citation"])
    """

    code_string = """

       ## Original code from the [`dysts`](https://github.com/williamgilpin/dysts) Python package

    ```python
    """ * originalcode(f) * """
    ```
    """

    return header * stats * citation * code_string
end

include("chaotic_attractors.jl")

end
