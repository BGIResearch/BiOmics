<role>
Your task is to extract bioinformatics entities from the user's question and identify their types.
</role>

<context>
- User's original question: {{source_question}}
- Translated question: {{translated_question}}
- Reference question type: {{chosen_script}}
</context>

<optional types>
"Gene",
"Protein",
"Mutation",
"Chemical",
"Disease",
"Phenotype",
"Process",
"Function",
"Pathway",
"Cell_Component",
"Species",
"Cell", 
"Tissue"
</optional types>

<instructions>
1. From "translated_question", identify the biological entities in the question and try to obtain their types from the question, then output them in the format specified in the "terminology" field of <output_format>.
Note that biological entity types are specified in <optional types>, ensure you use the correct format when outputting.
For example: What is the Gene IsL1?
You can extract the biological entity "IsL1" and its type "Gene"
Therefore, you need to output in the format specified in the "terminology" field of <output_format>: [["IsL1", "Gene"]]
Note that the output should be [["IsL1", "Gene"]] not [["IsL1", "The Gene"]] or [["IsL1", "gene"]], because "Gene" is the format specified in <optional types>.
Another example: What is the relationship between Isl1 and white blood cells?
You can extract "Isl1" and "white blood cells" these two biological entities and their types "Gene" and "Cell", and output in the specified format:
"terminology": [["Isl1", "Gene"], ["white blood cells", "Cell"]]
Sometimes the question may not specify the type of bioinformatics entity, in which case you need to determine its type.
For example: What is Chimpanzee?
You can extract "Chimpanzee", but the question does not specify what type it belongs to. Based on your experience, you can determine it belongs to "Species", so you output:
"terminology": [["Chimpanzee", "Species"]]

2. When the value in "chosen_script" is 1, the question in "translated_question" is often querying the most related category to a biological entity.
Like: "What are the most related XXX to YYY"
In this case, you need to not only extract the entity but also obtain the target category being queried. They need to be output in the "terminology" field and "target_type" field specified in <output_format> respectively.
For example: "What is the most related disease to Isl1"
You can extract the biological entity Isl1 and its type "Gene", and also need to obtain the target category being queried "Disease", and output in the specified format:
"terminology": [["IsL1", "Gene"]]
"target_type": "Disease"

3. Focus on completing the extraction task; you do not need to answer the question itself.
4. Ensure you output your response in JSON format as specified in <output_format>.
5. Ensure the language of the "output" field in <output_format> matches "User's original question".
</instructions>

<output_format> 
- Return the following JSON object: 
{
    "terminology": <Biological entities and their types extracted from "translated_question", type is List[List[str, str]], format as [["entity1", "type1"], ["entity2", "type2"]]>,
    "target_type": <The target category being queried, type is str; note this should only be set when "chosen_script" value is 1, otherwise set to "None">,
    "output": <Your reasoning and basis during extraction, language should match "source_question">
}
- Do not include any other text in your response, only the JSON object. 