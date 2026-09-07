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

-- Read literal \cssId anchors, not the mathematics inside their second argument.
-- TeX control sequences and comments must be skipped so examples/comments cannot
-- accidentally claim a target. Macro-generated IDs are deliberately unsupported.
local function math_ids(text)
  local ids, offset = {}, 1
  while offset <= #text do
    local char = text:sub(offset, offset)
    if char == "%" then
      offset = (text:find("\n", offset, true) or #text) + 1
    elseif char == "\\" then
      local command = text:match("^([%a]+)", offset + 1)
      if command then
        offset = offset + #command + 1
        if command == "cssId" then
          local tail = text:sub(offset)
          local id, finish = tail:match("^%s*{([^{}]*)}%s*{()")
          if not id or not id:match("^[A-Za-z][A-Za-z0-9_-]*$") then
            fail("\\cssId requires a literal ID (letter, then letters/digits/_/-) and a braced term")
          end
          ids[#ids + 1] = id
          offset = offset + finish - 1
        end
      else
        offset = offset + 2
      end
    else
      offset = offset + 1
    end
  end
  return ids
end

local function annotations(document)
  local anchors, notes, identifiers, slide = {}, {}, {}, 0
  document:walk({
    Inline = function(element)
      if element.identifier and element.identifier ~= "" then identifiers[element.identifier] = true end
    end,
    Block = function(element)
      if element.identifier and element.identifier ~= "" then identifiers[element.identifier] = true end
    end,
  })
  local blocks = document.blocks
  for index, block in ipairs(blocks) do
    if block.t == "Header" and block.level == 2 then slide = slide + 1 end
    blocks[index] = pandoc.Pandoc({ block }):walk({
      Math = function(math)
        for _, id in ipairs(math_ids(math.text)) do
          if anchors[id] or identifiers[id] then fail("duplicate equation anchor '" .. id .. "'") end
          anchors[id] = slide
        end
      end,
      Div = function(div)
        if not has_class(div, "lcm-math-note") then return nil end
        local target = div.attributes.target
        if not target or not target:match("^[A-Za-z][A-Za-z0-9_-]*$") then
          fail(".lcm-math-note requires a literal target ID")
        end
        if notes[target] then fail("duplicate explanation for '" .. target .. "'") end
        if div.attributes.trigger and div.attributes.trigger ~= "click" then
          fail("math notes use click-toggle; trigger must be omitted or 'click'")
        end
        if pandoc.utils.stringify(div.content):match("^%s*$") then fail("empty explanation for '" .. target .. "'") end
        div:walk({
          Header = function() fail("use paragraphs or bold text, not headings, inside math notes") end,
          RawBlock = function() fail("raw markup is not allowed inside math notes") end,
          RawInline = function() fail("raw markup is not allowed inside math notes") end,
          Image = function() fail("math notes accept text, lists, code and inline math, not images") end,
          Math = function(math)
            if #math_ids(math.text) > 0 then fail("equation anchors cannot be inside math notes") end
          end,
          Div = function(child)
            if child ~= div and has_class(child, "lcm-math-note") then fail("math notes cannot be nested") end
          end,
        })
        notes[target] = slide
        div.attributes.target = nil
        div.attributes.trigger = nil
        div.attributes["data-lcm-math-target"] = target
        div.attributes.hidden = ""
        return div
      end,
    }).blocks[1]
  end
  for target, owner in pairs(notes) do
    if not anchors[target] then fail("missing \\cssId anchor for '" .. target .. "'") end
    if anchors[target] ~= owner then fail("math note and anchor '" .. target .. "' must be on the same slide") end
  end
  document.blocks = blocks
  return document
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
    scripts = { "deck-init.js", "math-notes.js", "deck-controller.js" },
    stylesheets = { "deck-controller.css", "math-notes.css" },
  })
  meta["lcm-deck-contract"] = "1"
  return meta
end

return {
  { Pandoc = annotations },
  { Div = Div },
  { Meta = Meta },
}
