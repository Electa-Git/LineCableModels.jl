const COLORBAR_WIDTH = 140
const COLORBAR_TICK_LABEL_SIZE = 12
const COLORBAR_LABEL_SIZE = 14
const COLORBAR_LABEL_GAP = 8
const COLORBAR_ROW_GAP = 8
const LEGEND_DOCK_WIDTH = 220
const LEGEND_HEIGHT_TOLERANCE = 1
const LEGEND_MARGIN = (3.0f0, 3.0f0, 3.0f0, 3.0f0)
const OUTER_WINDOW_INSET = 3.0
const GRID_ROW_GAP = 6
const GRID_COLUMN_GAP = 6
const BUTTON_SIZE = 32
const BUTTON_ICON_SIZE = 18
const BACKGROUND_INTERACTIVE = :grey90
const BACKGROUND_EXPORT = :white
const BUTTON_BACKGROUND = Makie.RGBf(0.94, 0.94, 0.94)
const ICON_COLOR = Makie.RGBAf(0.15, 0.15, 0.15, 1.0)
const MI_REFRESH = "\uE5D5"
const MI_SAVE = "\uE161"
const ICON_FONT = joinpath(
    pkgdir(PlotBuilder),
    "assets",
    "fonts",
    "material-icons",
    "MaterialIcons-Regular.ttf"
)
const EXPORT_TIMESTAMP_FORMAT = "yyyymmdd_HHMMSS"
const EXPORT_FALLBACK_DIRECTORY = "linecablemodels-exports"

function _theme(; export_mode::Bool = false, export_theme::Symbol = :default)
    validate_export_theme(export_theme)
    background = export_mode ? BACKGROUND_EXPORT : BACKGROUND_INTERACTIVE
    base = export_mode && export_theme === :publication ? Makie.theme_latexfonts() : Theme()
    custom = Theme(
        backgroundcolor = background,
        fonts = (; icons = ICON_FONT),
        Axis = (
            titlesize = 15,
            xlabelsize = 14,
            ylabelsize = 14,
            xticklabelsize = 14,
            yticklabelsize = 14,
            xminorgridvisible = false,
            yminorgridvisible = false,
            xminorticksvisible = false,
            yminorticksvisible = false
        ),
        Button = (; buttoncolor = BUTTON_BACKGROUND),
        Legend = (; fontsize = 14, labelsize = 14, margin = LEGEND_MARGIN),
        Colorbar = (; labelsize = 14, ticklabelsize = 14)
    )
    return merge(base, custom)
end
