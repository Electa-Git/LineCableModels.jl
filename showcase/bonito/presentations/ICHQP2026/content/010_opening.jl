module ICHQP2026OpeningDeck

using Bonito
using Markdown
using ..PageAuthoring

const ASSET_DIRECTORY = normpath(joinpath(@__DIR__, "..", "assets"))
const ETCH_LOGO = joinpath(ASSET_DIRECTORY, "ETCH_LOGO_RGB_NEG.svg")
const ENERGYVILLE_LOGO = joinpath(ASSET_DIRECTORY, "ENERGYVILLE-LOGO.svg")
const KU_LEUVEN_LOGO = joinpath(ASSET_DIRECTORY, "kul_logo.svg")

function title_page(::Session, ::Nothing)
    footer = (
        prose(md"""
        **Amauri Martins**

        KU Leuven – Etch / EnergyVille

        [amauri.martinsbritto@kuleuven.be](mailto:amauri.martinsbritto@kuleuven.be)

        August 31, 2026 - Dresden, Germany
        """),
        logo_row(
            local_image(
                ENERGYVILLE_LOGO;
                alt = "EnergyVille logo"
            ),
            local_image(
                KU_LEUVEN_LOGO;
                alt = "KU Leuven logo"
            )
        )
    )
    body = title_canvas(
        "Uncertainty quantification of frequency dependent cable parameters";
        subtitle = "Understanding the stochasticity of harmonic impedances and injections for hosting capacity calculation",
        event = "22nd International Conference on Harmonics and Quality of Power",
        logo = local_image(
            ETCH_LOGO;
            alt = "KU Leuven Electrical Energy and Computer Architectures logo"
        ),
        footer
    )
    return (; body)
end

const DECK = deck_descriptor(
    id = "opening",
    group = "ICHQP2026",
    title = "Opening",
    order = 10,
    render = true,
    expand_navigation = true,
    pages = (
        deck_page(
        "Title slide";
        id = "title",
        class = "lc-title-page",
        build = title_page
    ),
    )
)

end

ICHQP2026OpeningDeck.DECK
