"""
Read table files (CSV, TSV, etc.) and convert to 2D list.
"""
import csv
from pathlib import Path
from langchain_core.tools import tool
from pydantic import BaseModel


class ReadTableInput(BaseModel):
    file_path: str = ""
    delimiter: str = None


def read_table_to_list(file_path: str, delimiter: str = None) -> list[list]:
    """
    Read a table file and return as a 2D list.
    
    Args:
        file_path: Path to the table file (CSV, TSV, etc.)
        delimiter: Column delimiter. If None, auto-detect based on file extension
                 (.csv -> ',', .tsv -> '\t', others -> '\t')
    
    Returns:
        2D list where each row is a list of values
        
    Example:
        >>> data = read_table_to_list('data.csv')
        >>> print(data[0])  # Header row
        ['Name', 'Age', 'City']
        >>> print(data[1])  # First data row
        ['Alice', '25', 'NYC']
    """
    file_path = Path(file_path)
    
    if delimiter is None:
        if file_path.suffix.lower() == '.csv':
            delimiter = ','
        elif file_path.suffix.lower() == '.tsv':
            delimiter = '\t'
        else:
            delimiter = '\t'
    
    result = []
    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter=delimiter)
        for row in reader:
            result.append(row)
    
    return result


@tool(args_schema=ReadTableInput)
def read_table_tool(file_path: str, delimiter: str = None) -> dict:
    """
    Read a CSV or TSV table file and return its content as a 2D list.
    Use this tool when the user provides a .csv or .tsv file.
    
    Args:
        file_path: Path to the table file (CSV, TSV, etc.)
        delimiter: Column delimiter. If None, auto-detect based on file extension
    
    Returns:
        A dictionary containing the table data, dimensions, and preview
    """
    try:
        data = read_table_to_list(file_path, delimiter)
        num_rows = len(data)
        num_cols = len(data[0]) if data else 0
        
        # Get header and first few rows as preview
        header = data[0] if data else []
        preview_rows = data[1:6] if len(data) > 1 else []  # First 5 data rows
        
        # Analyze column types and missing values
        column_info = []
        for col_idx, col_name in enumerate(header):
            col_values = [row[col_idx] if col_idx < len(row) else "" for row in data[1:]]
            non_empty = [v for v in col_values if v.strip() != ""]
            missing_count = len(col_values) - len(non_empty)
            
            # Infer column type
            col_type = "string"
            if non_empty:
                try:
                    [float(v) for v in non_empty[:100]]  # Check first 100 values
                    col_type = "numeric"
                except ValueError:
                    pass
            
            column_info.append({
                "name": col_name,
                "type": col_type,
                "non_empty": len(non_empty),
                "missing": missing_count
            })
        
        # Format preview as markdown table
        preview_str = "| " + " | ".join(str(cell)[:20] for cell in header) + " |\n"
        preview_str += "| " + " | ".join("---" for _ in header) + " |\n"
        for row in preview_rows:
            row_cells = [str(row[i])[:20] if i < len(row) else "" for i in range(len(header))]
            preview_str += "| " + " | ".join(row_cells) + " |\n"
        
        # Build description in markdown format
        description = (
            f"**Table:** {Path(file_path).name}\n\n"
            f"**Dimensions:** {num_rows} rows × {num_cols} columns (including header)\n\n"
            f"**Columns:** {', '.join(header)}\n\n"
            f"### Column Details\n\n"
        )
        for info in column_info:
            description += f"- **{info['name']}**: {info['type']}, {info['non_empty']} values, {info['missing']} missing\n"
        description += f"\n### Data Preview (first 5 rows)\n\n{preview_str}"
        
        return {
            "success": True,
            "file_path": file_path,
            "num_rows": num_rows,
            "num_cols": num_cols,
            "header": header,
            "column_info": column_info,
            "preview": preview_rows,
            "description": description
        }
    except Exception as e:
        return {
            "success": False,
            "file_path": file_path,
            "error": str(e),
            "description": f"Failed to read table: {str(e)}"
        }
