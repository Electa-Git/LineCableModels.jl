# Input fixtures

`mv_cable_design.json` is UTF-8 JSON using fixture schema version 1. Its root object is
an encoded `LineCableModels.DataModel.CablesLibrary`. Nested records use fully
qualified `__julia_type__` discriminators and SI-valued numeric fields accepted by the
LineCableModels 0.2 serialisation schema. JSON object ordering does not affect
loading.

The file describes one 18/30 kV three-component cable design. Tests read it
through a fixture factory so every caller receives fresh mutable model objects.
A schema change requires an explicit fixture migration and matching
round-trip and malformed-input tests.
