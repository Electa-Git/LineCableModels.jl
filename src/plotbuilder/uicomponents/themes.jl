using Printf: @sprintf

# -----------------------------------------------------------------------------
# Constants
# -----------------------------------------------------------------------------

const FIG_SIZE = (800, 600)
const FIG_PADDING = (80, 60, 40, 40) # left, right, bottom, top
const CTLBAR_HEIGHT = 36
const STATUSBAR_HEIGHT = 20
const GRID_ROW_GAP = 6
const GRID_COL_GAP = 6
const LEGEND_GAP = 4
const LEGEND_WIDTH = 140
const COLORBAR_GAP = 4
const CTLBAR_GAP = 2
const BUTTON_MIN_WIDTH = 32
const BUTTON_ICON_SIZE = 18
const BUTTON_TEXT_FONT_SIZE = 15
const AXIS_TITLE_FONT_SIZE = 15
const AXIS_LABEL_FONT_SIZE = 14
const AXIS_TICK_FONT_SIZE = 14
const STATUS_FONT_SIZE = 10
const BG_COLOR_INTERACTIVE = :grey90
const BG_COLOR_EXPORT = :white
const ICON_COLOR_ACTIVE = Makie.RGBAf(0.15, 0.15, 0.15, 1.0)
const ICON_COLOR_DISABLED = Makie.RGBAf(0.55, 0.55, 0.55, 1.0)
const TICK_FMT = x -> @sprintf("%.4g", x)
const TICKFORMATTER = values -> [TICK_FMT(v) for v in values]
const EXPORT_TIMESTAMP_FORMAT = "yyyymmdd_HHMMSS"
const EXPORT_EXTENSION = "svg"

# -----------------------------------------------------------------------------
# Material UI icons
# -----------------------------------------------------------------------------
const MI_REFRESH = "\uE5D5"  # Material Icons: 'refresh'
const MI_SAVE    = "\uE161"  # Material Icons: 'save'

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

with_icon(icon::AbstractString; text::AbstractString = "",
	isize::Int = BUTTON_ICON_SIZE, tsize::Int = BUTTON_TEXT_FONT_SIZE, color = :black,
	gap::Int = 2,
	dy_icon::Float64 = -0.18, dy_text::Float64 = 0.0) =
	text == "" ?
	rich(icon; font = :icons, fontsize = isize, color = color, offset = (0, dy_icon)) :
	rich(
		rich(icon; font = :icons, fontsize = isize, color = color, offset = (0, dy_icon)),
		rich(" "^gap; font = :regular, fontsize = tsize, color = color),
		rich(text; font = :regular, fontsize = tsize, color = color, offset = (0, dy_text)),
	)

"""
	make_theme(; interactive::Bool, use_latex_fonts::Bool)

Returns the package-specific Theme delta.
"""
function make_theme(; interactive::Bool, use_latex_fonts::Bool)
	background = interactive ? BG_COLOR_INTERACTIVE : BG_COLOR_EXPORT

	# Base configuration
	config = Dict{Symbol, Any}(
		:backgroundcolor => background,
		:Axis => (
			titlesize      = AXIS_TITLE_FONT_SIZE,
			xlabelsize     = AXIS_LABEL_FONT_SIZE,
			ylabelsize     = AXIS_LABEL_FONT_SIZE,
			xticklabelsize = AXIS_TICK_FONT_SIZE,
			yticklabelsize = AXIS_TICK_FONT_SIZE,
			xtickformat    = TICKFORMATTER,
			ytickformat    = TICKFORMATTER,
		),
		:Legend => (
			fontsize = AXIS_LABEL_FONT_SIZE,
			labelsize = AXIS_LABEL_FONT_SIZE,
		),
		:Colorbar => (
			labelsize     = AXIS_LABEL_FONT_SIZE,
			ticklabelsize = AXIS_TICK_FONT_SIZE,
		),
	)

	# Conditional logic: Fonts
	# 1. Latex fonts (Export only)
	if use_latex_fonts && !interactive
		# merge! is safe on Dicts
		merge!(config, Makie.theme_latexfonts().attributes)
	end

	# 2. Icon fonts (Always try to load)
	font_path = joinpath(
		pkgdir(@__MODULE__),
		"assets",
		"fonts",
		"material-icons",
		"MaterialIcons-Regular.ttf",
	)
	if isfile(font_path)
		current_fonts = get(config, :fonts, (;))
		# Convert to NamedTuple to simple merge
		new_fonts = merge(current_fonts, (; icons = font_path))
		config[:fonts] = new_fonts
	end

	return Makie.Theme(; config...)
end
