# Parameter-set JSON

ChargeFW parameter sets are JSON documents consumed by the native
`parameters::load_parameter_set_json*()` functions. Bundled sets use this format. The CLI and Python
package expose only the installed bundled catalog. Custom parameter-set loading is available through the
native C++ API only.

A parameter set declares values for one method. Method requirements determine which common, atom, and
bond parameter names must be present. During assessment, ChargeFW verifies those names, then classifies
every relevant atom and bond against the keys below.

## Document structure

```json
{
  "metadata": {
    "id": "example",
    "method": "sqeqp",
    "name": "Example",
    "publication": "10.0000/example",
    "notes": "Optional provenance",
    "priority": 0
  },
  "common": {
    "names": ["scale"],
    "values": [1.0]
  },
  "atom": {
    "names": ["electronegativity", "hardness", "width", "q0"],
    "data": [
      {
        "key": {"element": "C", "classifier": "hbo", "type": "1"},
        "values": [5.0, 10.0, 1.0, 0.0]
      }
    ]
  },
  "bond": {
    "names": ["kappa"],
    "data": [
      {
        "key": {
          "atoms": [
            {"element": "C", "classifier": "hbo", "type": "1"},
            {"element": "O", "classifier": "hbo", "type": "1"}
          ],
          "bond": {"classifier": "bo", "type": "1"}
        },
        "values": [0.5]
      }
    ]
  }
}
```

`metadata.id` and `metadata.method` are required. The non-empty ID is the parameter set's stable identity
and must be unique within a loaded catalog; it does not depend on the source filename. `name`,
`publication`, `notes`, and `priority` are optional. Higher priorities are selected before lower priorities;
equal priorities use the stable ID as the tie-breaker. `priority` must be an integer from `0` through
`65535`.

`common`, `atom`, and `bond` are optional. Each section has a `names` array. Its names define the order of
that section's numeric `values`; every values array must contain exactly one finite number per name.

## Atom parameters

Each `atom.data` entry has an atom key and a values array:

```json
{
  "key": {"element": "C", "classifier": "bonded", "type": "CNO"},
  "values": [2.6846, 1.571, 0.5771, 0.8167]
}
```

The key fields are all required:

- `element` is a supported element symbol, or `"*"` to match every element.
- `classifier` selects how ChargeFW derives the atom's `type`.
- `type` is the expected string produced by that classifier.

## Bond parameters

Each `bond.data` entry uses two atom keys and one bond key:

```json
{
  "key": {
    "atoms": [
      {"element": "C", "classifier": "bonded", "type": "CNO"},
      {"element": "N", "classifier": "bonded", "type": "CCH"}
    ],
    "bond": {"classifier": "bo", "type": "1"}
  },
  "values": [0.5]
}
```

`atoms` must contain exactly two atom-key objects. Bond matching is unordered: either molecular bond
endpoint may match either object. The `bond` object requires `classifier` and `type` fields and determines
the bond property to match.

## Classifiers

Atom keys support these classifiers:

| Classifier | Derived type | Example |
| --- | --- | --- |
| `plain` | Always `"*"` | `{"element": "O", "classifier": "plain", "type": "*"}` |
| `hbo` | Highest bond order among bonds incident to the atom, as a decimal string | An atom incident to a double bond has type `"2"` |
| `bonded` | Concatenation of neighboring element symbols, sorted lexically | Carbon bonded to C, N, and O has type `"CNO"` |

Bond keys support these classifiers:

| Classifier | Derived type | Example |
| --- | --- | --- |
| `plain` | Always `"*"` | `{"classifier": "plain", "type": "*"}` |
| `bo` | Bond order as `"1"`, `"2"`, or `"3"` | A double bond has type `"2"` |

`plain` keys therefore require `"*"` as their type to match. `bonded` ignores bond order and formal
charge; it uses only the sorted immediate-neighbor element symbols.

## Entry selection

ChargeFW evaluates entries in their `data` order and selects the first matching entry. Parameter files can
therefore express precedence by placing a more specific key before a broader key. Duplicate or overlapping
keys are accepted; later matching entries are not used.

## Strict and permissive matching

Strict classification is the default. With permissive classification (`--permissive-types` in the CLI or
`parameter_matching="permissive"` in Python), `hbo` and `bo` lower a derived order above one by one before
matching. For example, a double-bond order/type may match `"1"` and a triple-bond order/type may match
`"2"`. `plain` and `bonded` matching do not change. ChargeFW tries strict matching first and tries
permissive matching only when no strict entry matches.
