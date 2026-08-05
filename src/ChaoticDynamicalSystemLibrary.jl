module ChaoticDynamicalSystemLibrary

using JSON, Markdown
using LinearAlgebra
using SciMLBase
using ComponentArrays
using SimpleUnPack

vecvec2mat(vv) = hcat(vv...)' # transpose is because of the way the data is stored in the JSON file

function format_parameters(params::JSON.Object)
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

"""
    attractor_property(name::AbstractString, property::AbstractString)

Returns the given property from the dysts dataset for the given attractor by name.
"""
function attractor_property(name::AbstractString, property::AbstractString)
    data = ATTRACTOR_DATA[name]
    if property in keys(data)
        return data[property]
    else
        throw(ArgumentError("\"$property\" not found for attractor \"$name\"."))
    end
end

# define wrapper for important propeties

attractor_property(obj, property::AbstractString) = attractor_property(string(obj), property)
correlation_dimension(obj) = attractor_property(obj, "correlation_dimension")::Float64
embedding_dimension(obj) = attractor_property(obj, "embedding_dimension")::Int
hamiltonian(obj) = attractor_property(obj, "hamiltonian")::Bool
kaplan_yorke_dimension(obj) = attractor_property(obj, "kaplan_yorke_dimension")::Float64
lyapunov_spectrum_estimated(obj) = attractor_property(obj, "lyapunov_spectrum_estimated")::Union{Nothing,Vector{Any}}
maximum_lyapunov_estimated(obj) = attractor_property(obj, "maximum_lyapunov_estimated")::Float64
multiscale_entropy(obj) = attractor_property(obj, "multiscale_entropy")::Float64
nonautonomous(obj) = attractor_property(obj, "nonautonomous")::Bool
period(obj) = attractor_property(obj, "period")::Float64
pesin_entropy(obj) = attractor_property(obj, "pesin_entropy")::Float64

function make_docstring(f)
    header = """
           $f()

       $(attractor_property(f, "description"))
    """

    stats = """

       ## Stats
    """

    for (name, property) in [
        ("Correlation dimension", "correlation_dimension"),
        ("Embedding dimension", "embedding_dimension"),
        ("Hamiltonian", "hamiltonian"),
        ("Kaplan-Yorke dimension", "kaplan_yorke_dimension"),
        ("Lyapunov spectrum (estimated)", "lyapunov_spectrum_estimated"),
        ("Maximum Lyapunov exponent (estimated)", "maximum_lyapunov_estimated"),
        ("Multiscale Entropy:", "multiscale_entropy"),
        ("Non-autonomous", "nonautonomous"),
        ("Period", "period"),
        ("Pesin entropy", "pesin_entropy")
    ]
        try
            value = attractor_property(f, property)
            stats *= """
          - `$name`: `$value`
       """
        catch
        end
    end

    citation = """

       ## References
       $(attractor_property(f, "citation"))
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

AVAILABLE_ATTRACTORS::Vector{String} = filter(x -> isdefined(@__MODULE__, Symbol(x)), collect(keys(ATTRACTOR_DATA)))

get_attractor_by_name(obj) = get_attractor_by_name(string(obj))

"""
    get_attractor_by_name(name::AbstractString)

Returns the corresponding julia object for the given attractor name.
"""
function get_attractor_by_name(name::AbstractString)
    sym = Symbol(name)
    if isdefined(@__MODULE__, Symbol(name))
        return getfield(@__MODULE__, sym)
    else
        throw(ArgumentError("No attractor \"$name\" found."))
    end
end

end