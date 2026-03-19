using Makie

using Base: basename, mod1
using Dates: format, now
using Printf: @sprintf

import LineCableModels.PlotBuilder.BackendHandler: BackendHandler, next_fignum

using LineCableModels.PlotBuilder.PlotUIComponents:
	PlotAssembly,
	PlotBuildArtifacts,
	ControlButtonSpec,
	ControlToggleSpec,
	ControlReaction,
	_make_window,
	_run_plot_pipeline,
	with_plot_theme,
	ensure_export_background!,
	with_icon,
	MI_REFRESH,
	MI_SAVE,
	ICON_TTF,
	AXIS_LABEL_FONT_SIZE,
	clear_status!,
	TICKFORMATTER,
	EXPORT_EXTENSION,
	EXPORT_TIMESTAMP_FORMAT

using LineCableModels.PlotBuilder: AbstractPlotSpec,
	_sanitize_filename_plot,
	_default_export_path,
	_save_plot_export,
	_parse_values_expr,
	_autoscale_axis,
	_render_plot_specs,
	_build_common_plot_controls
