-- Exercise the production AST filter through Pandoc, without a Quarto build.
quarto = { doc = { add_html_dependency = function() end } }
return dofile(pandoc.path.join({ pandoc.path.directory(PANDOC_SCRIPT_FILE),
  "..", "..", "_extensions", "lcm-deck", "deck.lua" }))
