# Data Processing Patterns for BiOmics-BRICK

When analyzing `brick qr` outputs, especially for biomedical QA (e.g., BioASQ), use these patterns in an `execute_code` block.

## Efficient Evidence Search

If the answer is not in a structured KG relationship (for example marker evidence queried as `produce`), search the `original_relation` column which contains source text snippets.

```python
import pandas as pd
df = pd.read_csv('results.csv')

# Search for multiple keywords across all columns (fuzzy matching)
keywords = ['resistance', 'sensitive', 'sensitiz']
pattern = '|'.join(keywords)
hits = df[df.astype(str).apply(lambda x: x.str.contains(pattern, case=False)).any(axis=1)]

# Display relevant snippets
if not hits.empty:
    for i, row in hits.iterrows():
        print(f"Target: {row['path.2.name']}")
        print(f"Evidence: {row.get('path.1.original_relation', 'N/A')}")
```

## Common Pitfalls

### Python Syntax in Sandbox
Avoid literal newlines inside `print()` strings or f-strings that span multiple lines without proper delimiters, as they trigger `SyntaxError: unterminated string literal`.

**Bad**:
```python
print("
Target: ...")
```

**Good**:
```python
print(f"
Target: {row['path.2.name']}")
# OR
print("
" + "Target: " + str(row['path.2.name']))
```

### Column Naming
`brick qr` uses a path-based naming convention:
- `path.0.name`: Source entity
- `path.1.relation`: Relation type
- `path.1.original_relation`: Source text snippet (the most valuable for "Why?" questions)
- `path.1.condition`: Experimental conditions or additional context
- `path.2.name`: Target entity

## Multi-Target Path Validation
To validate a mechanism (e.g., "LKB1 loss -> Leptin sensitivity"), query the source (LKB1) and multiple targets (Leptin, UCP2, HIF1A) simultaneously.

```bash
brick qr "LKB1" -t "LEP,UCP2,HIF1A" --st Gene --tt Gene -n 50 -o validation.csv
```
