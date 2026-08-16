# Input fixtures

`mv_cable_design.json` is UTF-8 JSON using fixture schema version 1. Its root object is
a serialized `LineCableModels.DataModel.CablesLibrary`; nested records use fully
qualified `__julia_type__` discriminators and SI-valued numeric fields accepted by the
LineCableModels 0.2 serialization contract. JSON object ordering is not significant.

The file represents one 18/30 kV three-component cable design and is immutable test
input. Tests must load it through a fixture factory so every caller receives fresh
mutable model objects. Schema changes require an explicit fixture migration and matching
round-trip and malformed-input contract updates.
