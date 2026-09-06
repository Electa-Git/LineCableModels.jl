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

local function valid_public_url(value)
  if value:sub(1, 1) == "/" and value:sub(1, 2) ~= "//" then
    return true
  end
  return value:match("^https://") ~= nil or value:match("^http://") ~= nil
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

    local public_url = stringify(kwargs["public-url"])
    if public_url ~= "" and not valid_public_url(public_url) then
      return quarto.shortcode.error_output(
        "bonito",
        "public-url must be an absolute same-origin path or an http(s) URL",
        "block"
      )
    end

    if public_url ~= "" then
      local html = string.format(
        '<figure class="lcm-live-viewport" data-lcm-live="true" data-lcm-title="%s">' ..
        '<iframe class="lc-widget-frame" data-lcm-src="%s" title="%s" allowfullscreen ' ..
        'style="--lc-widget-height: %s;"></iframe>' ..
        '<figcaption class="lcm-live-placeholder">' ..
        '<strong>%s</strong>' ..
        '<span>Interactive view omitted from this static presentation surface.</span>' ..
        '<a href="%s">Open in the LineCableModels playground</a>' ..
        '</figcaption></figure>',
        escape_attribute(title),
        escape_attribute(route),
        escape_attribute(title),
        escape_attribute(height),
        escape_attribute(title),
        escape_attribute(public_url)
      )
      return pandoc.RawBlock("html", html)
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
