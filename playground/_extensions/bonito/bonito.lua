local function stringify(value)
  if value == nil then
    return ""
  end
  return pandoc.utils.stringify(value)
end

local function escape_attribute(value)
  return stringify(value)
    :gsub("&", "&amp;")
    :gsub('"', "&quot;")
    :gsub("<", "&lt;")
    :gsub(">", "&gt;")
end

local function valid_height(value)
  local number, unit = value:match("^(%d+%.?%d*)([%a%%]+)$")
  if number == nil then
    return false
  end
  return unit == "rem" or unit == "px" or unit == "vh" or
    unit == "dvh" or unit == "%"
end

return {
  bonito = function(args, kwargs)
    if not quarto.doc.is_format("html:js") then
      return pandoc.Para({
        pandoc.Emph({pandoc.Str("Interactive Bonito viewport available in HTML.")})
      })
    end

    local route = stringify(kwargs.route or args[1])
    if route == "" or route:sub(1, 1) ~= "/" then
      return quarto.shortcode.error_output(
        "bonito",
        "route must be an absolute same-origin path such as /widgets/slider",
        "block"
      )
    end

    local height = stringify(kwargs.height)
    if height == "" then
      height = "20rem"
    elseif not valid_height(height) then
      return quarto.shortcode.error_output(
        "bonito",
        "height must use rem, px, vh, dvh, or %",
        "block"
      )
    end

    local title = stringify(kwargs.title)
    if title == "" then
      title = "Live Bonito widget"
    end

    local html = string.format(
      '<iframe class="lc-widget-frame" src="%s" title="%s" loading="lazy" allowfullscreen style="--lc-widget-height: %s;"></iframe>',
      escape_attribute(route),
      escape_attribute(title),
      escape_attribute(height)
    )
    return pandoc.RawBlock("html", html)
  end
}
