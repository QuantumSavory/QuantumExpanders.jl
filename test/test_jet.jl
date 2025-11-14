@testitem "JET checks" tags=[:jet] begin
using JET
using QuantumExpanders

using Oscar, Multigraphs

rep = report_package("QuantumExpanders";
    ignored_modules=(
        AnyFrameModule(Nemo),
        AnyFrameModule(Oscar),
        AnyFrameModule(Multigraphs),
    )
)
@show rep
@test length(JET.get_reports(rep)) <= 13
end
