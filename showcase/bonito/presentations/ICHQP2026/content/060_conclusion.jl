module ICHQP2026ConclusionDeck

using Bonito
using Markdown
using ..PageAuthoring

const ASSET_DIRECTORY = normpath(joinpath(@__DIR__, "..", "assets"))
const SHARED_IMAGE_DIRECTORY = normpath(
    joinpath(@__DIR__, "..", "..", "..", "..", "..", "assets", "img")
)
const ETCH_LOGO = joinpath(ASSET_DIRECTORY, "ETCH_LOGO_RGB_NEG.svg")
const ENERGYVILLE_LOGO = joinpath(ASSET_DIRECTORY, "ENERGYVILLE-LOGO.svg")
const KU_LEUVEN_LOGO = joinpath(ASSET_DIRECTORY, "kul_logo.svg")
const FUNDING_LOGO = joinpath(SHARED_IMAGE_DIRECTORY, "sponsorlogo_WEWIS.svg")

function conclusion_page(::Session, ::Nothing)
    remarks = webpart(
        prose(md"""
        - Uncertainty-aware models turn nominal hosting-capacity limits into
          defensible planning margins for cables and connected equipment.

        - Frequency-dependent uncertainty envelopes reveal where parameter
          variability can amplify harmonics or move critical resonances.

        - Planning around both uncertainty and sensitivity helps target
          reinforcement, protection, and monitoring for a more resilient grid.

        Thank you!
        """);
        kind = :panel,
        class = "lc-conclusion-remarks"
    )
    closing = webpart(
        logo_row(
            local_image(ETCH_LOGO; alt = "Etch logo"),
            local_image(ENERGYVILLE_LOGO; alt = "EnergyVille logo"),
            local_image(KU_LEUVEN_LOGO; alt = "KU Leuven logo"),
            local_image(FUNDING_LOGO; alt = "Funding agency logo");
            class = "lc-conclusion-logos"
        );
        kind = :panel,
        class = "lc-conclusion-closing"
    )
    body = slide(
        "Concluding remarks",
        layout_two_rows(
            remarks,
            closing;
            ratio = :tall_short,
            height = "100%",
            class = "lc-conclusion-layout"
        )
    )
    return (; body)
end

const DECK = deck_descriptor(
    id = "conclusion",
    group = "ICHQP2026",
    title = "Conclusion",
    order = 60,
    render = true,
    pages = (
        deck_page(
        "Concluding remarks";
        id = "concluding-remarks",
        class = "lc-conclusion-slide lc-fill-page lc-fluid-type",
        build = conclusion_page
    ),
    )
)

end

ICHQP2026ConclusionDeck.DECK
