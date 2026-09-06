local layouts = {
  ["full-canvas"] = 1,
  ["balanced"] = 2,
  ["feature-sidebar"] = 2,
  ["top-split"] = 3,
  ["dashboard-grid"] = 4,
  ["media-story"] = 2,
}

local function has_class(element, name)
  for _, class in ipairs(element.classes) do
    if class == name then
      return true
    end
  end
  return false
end

local function slot_count(div)
  local count = 0
  for _, block in ipairs(div.content) do
    if block.t == "Div" and has_class(block, "lcm-slot") then
      count = count + 1
    end
  end
  return count
end

local function layout_name(div)
  for _, class in ipairs(div.classes) do
    local name = class:match("^lcm%-layout%-(.+)$")
    if name ~= nil then
      return name
    end
  end
  return nil
end

local function fail(message)
  error("lcm-deck: " .. message, 0)
end

function Div(div)
  local name = layout_name(div)
  if name == nil then
    return nil
  end

  local expected = layouts[name]
  if expected == nil then
    fail("unknown layout '" .. name .. "'")
  end

  local actual = slot_count(div)
  if actual ~= expected then
    fail("layout '" .. name .. "' requires " .. expected ..
      " direct .lcm-slot children; found " .. actual)
  end

  if not has_class(div, "lcm-layout") then
    div.classes:insert("lcm-layout")
  end
  div.attributes["data-lcm-layout"] = name
  return div
end

function Meta(meta)
  quarto.doc.add_html_dependency({
    name = "lcm-deck-controller",
    version = "0.1.0",
    scripts = { "deck-init.js", "deck-controller.js" },
    stylesheets = { "deck-controller.css" },
  })
  meta["lcm-deck-contract"] = "1"
  return meta
end

return {
  { Div = Div },
  { Meta = Meta },
}
